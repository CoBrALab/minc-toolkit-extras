#!/bin/bash
#set -e

rounds=6
output=$1
tmpdir=$(mktemp -d)
shift 1

for file in "$@"
do
    mkdir $tmpdir/$(basename $file)
    ImageMath 3 $tmpdir/$(basename $file)/trunc.mnc TruncateImageIntensity $file 0 0.995 2048
    ImageMath 3 $tmpdir/$(basename $file)/hotspots.mnc TruncateImageIntensity $file 0.95 1 2048
    ThresholdImage 3 $tmpdir/$(basename $file)/hotspots.mnc $tmpdir/$(basename $file)/hotmask_inv.mnc Otsu 1
    ImageMath 3 $tmpdir/$(basename $file)/hotmask.mnc Neg $tmpdir/$(basename $file)/hotmask_inv.mnc

    bestlinreg_g -clobber -nmi -lsq12 -target_mask $mnimask $tmpdir/$(basename $file)/trunc.mnc $mnimodel $tmpdir/$(basename $file)/to_mni.xfm
    itk_resample --clobber --invert_transform --labels --like $tmpdir/$(basename $file)/trunc.mnc --transform $tmpdir/$(basename $file)/to_mni.xfm $mnimask $tmpdir/$(basename $file)/brainmask.mnc
    mincmorph -successive EE $tmpdir/$(basename $file)/brainmask.mnc $tmpdir/$(basename $file)/brainmaskE.mnc

    #cp $tmpdir/$(basename $file)/trunc.mnc $tmpdir/$(basename $file)/input.mnc

    ThresholdImage 3 $tmpdir/$(basename $file)/trunc.mnc $tmpdir/$(basename $file)/otsu.mnc Otsu 1
    mincdefrag $tmpdir/$(basename $file)/otsu.mnc $tmpdir/$(basename $file)/otsu_defrag.mnc 1 6
    autocrop -isoexpand 50mm $tmpdir/$(basename $file)/otsu_defrag.mnc $tmpdir/$(basename $file)/otsu_expanded.mnc
    itk_morph --exp 'D[25] E[25]' $tmpdir/$(basename $file)/otsu_expanded.mnc $tmpdir/$(basename $file)/otsu_expanded_closed.mnc
    itk_resample --clobber --labels --like $tmpdir/$(basename $file)/otsu.mnc $tmpdir/$(basename $file)/otsu_expanded_closed.mnc $tmpdir/$(basename $file)/headmask.mnc

    ImageMath 3 $tmpdir/$(basename $file)/brainmaskE_novent.mnc m $tmpdir/$(basename $file)/brainmaskE.mnc $tmpdir/$(basename $file)/otsu.mnc
    ImageMath 3 $tmpdir/$(basename $file)/brainmaskE_novent_nohot.mnc m $tmpdir/$(basename $file)/brainmaskE_novent.mnc $tmpdir/$(basename $file)/hotmask.mnc

    N4BiasFieldCorrection -d 3 -i $tmpdir/$(basename $file)/trunc.mnc -w $tmpdir/$(basename $file)/brainmaskE_novent_nohot.mnc -b [400] -c [200x200x200x200x200,0.0] -v -r 1 -x $tmpdir/$(basename $file)/headmask.mnc -o $tmpdir/$(basename $file)/input.mnc

    # for i in {1..3}
    # do
    #     itk_minc_nonlocal_filter --clobber --beta 0.5 --anlm $tmpdir/$(basename $file)/input.mnc $tmpdir/$(basename $file)/denoise.mnc
    #     ThresholdImage 3 $tmpdir/$(basename $file)/denoise.mnc $tmpdir/$(basename $file)/kmeans.mnc Kmeans 2 $tmpdir/$(basename $file)/brainmaskE.mnc
    #     minclookup -clobber -discrete -lut_string "3 1" $tmpdir/$(basename $file)/kmeans.mnc $tmpdir/$(basename $file)/wm.mnc
    #     N4BiasFieldCorrection -d 3 -i $tmpdir/$(basename $file)/trunc.mnc -w $tmpdir/$(basename $file)/wm.mnc -b [400] -c [200x200x200x200x200,0.0] -v -r 1 -x $tmpdir/$(basename $file)/headmask.mnc -o $tmpdir/$(basename $file)/input.mnc
    # done
    ImageMath 3 $tmpdir/$(basename $file).N4.mnc m $tmpdir/$(basename $file)/input.mnc $tmpdir/$(basename $file)/headmask.mnc
    #cp $tmpdir/$(basename $file)/input.mnc $tmpdir/$(basename $file).N4.mnc
done

mincbigaverage --robust --sdfile $tmpdir/sd0.mnc $tmpdir/*.N4.mnc $tmpdir/average0.mnc
target=$tmpdir/average0.mnc

iters="2000"
shrink=$(expr $rounds - 1)
sigma=$(expr $rounds - 2)

for ((n=1;n<$rounds;n++))
do
  for file in $tmpdir/*.N4.mnc
  do
    antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 --minc --verbose --interpolation BSpline \
                --output [$tmpdir/$(basename $file)-$(basename $target),$tmpdir/$(basename $file)-$(basename $target)] \
                --transform Rigid[0.1] --metric GC[$target,$file,1] \
                --convergence [$iters,1e-6,10] --shrink-factors $shrink --smoothing-sigmas $sigma
  done
  mincaverage -normalize -sdfile $tmpdir/sd$n.mnc $tmpdir/*-$(basename $target) $tmpdir/average$n.mnc
  target=$tmpdir/average$n.mnc
  iters+="x2000"
  shrink+="x$(expr $rounds - $n - 1)"
  sigma+="x$(expr $rounds - $n - 2)"
done

#cp $target $output
mincmath -clamp -const2 0 $(mincstats -quiet -max $target) $target $output
rm -r $tmpdir
