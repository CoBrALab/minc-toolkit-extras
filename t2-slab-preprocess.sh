#!/bin/bash
#Preprocessing script which registers slab to wholebrain, back transforms mask and does bias field correction, and does brain extraction
set -euo pipefail
set -x

tmpdir=$(mktemp -d)

slab=$1
wholebrain=$2
wholebrain_mask=$3
output=$4

minccalc -byte -unsigned -expression '1' ${slab} ${tmpdir}/initmask.mnc
ImageMath 3 ${tmpdir}/slab_mask.mnc ThresholdAtMean ${slab} 0.5

N4BiasFieldCorrection --verbose -d 3 -s 4 -w ${tmpdir}/slab_mask.mnc -x ${tmpdir}/initmask.mnc \
  -b [200] -c [300x300x300x300,1e-5] \
  -i ${slab} \
  -o ${tmpdir}/slab_N4.mnc -r 0

antsRegistration_affine.sh --fixed-mask ${wholebrain_mask} ${tmpdir}/slab_N4.mnc ${wholebrain} ${tmpdir}/towholebrain.xfm


#antsApplyTransforms -d 3 -i ${wholebrain_mask} -r ${slab} -t [${tmpdir}/towholebrain.xfm,1] -n GenericLabel -o ${tmpdir}/slabmask.mnc --verbose
mincresample -unsigned -byte -labels -keep -near -like ${slab} -transform ${tmpdir}/towholebrain.xfm ${wholebrain_mask} ${tmpdir}/slabmask.mnc


#Calculate bins for N4 with Freedman-Diaconis’s Rule
n4weight=${tmpdir}/slabmask.mnc
n4input=${slab}
n4initmask=${tmpdir}/initmask.mnc
n4corrected=${tmpdir}/N4.mnc
n4bias=${tmpdir}/bias.mnc
min=$(mincstats -quiet -min -mask ${n4weight} -mask_range 1e-9,inf ${n4input})
max=$(mincstats -quiet -max -mask ${n4weight} -mask_range 1e-9,inf ${n4input})
npoints=$(mincstats -quiet -count -mask ${n4weight} -mask_range 1e-9,inf ${n4input})
pct25=$(mincstats -quiet -pctT 25 -hist_bins 10000 -mask ${n4weight} -mask_range 1e-9,inf ${n4input})
pct75=$(mincstats -quiet -pctT 75 -hist_bins 10000 -mask ${n4weight} -mask_range 1e-9,inf ${n4input})
histbins=$(python -c "print( int((float(${max})-float(${min}))/(2.0 * (float(${pct75})-float(${pct25})) * float(${npoints})**(-1.0/3.0)) ))" )

N4BiasFieldCorrection --verbose -d 3 -s 4 -w ${n4weight} -x ${n4initmask} \
  -b [200] -c [300x300x300x300,1e-5] --histogram-sharpening [0.05,0.01,${histbins}] \
  -i ${n4input} \
  -o [${n4corrected},${n4bias}] -r 0

#Demean bias field estimate and recorrect file
ImageMath 3 ${n4bias} / ${n4bias} $(mincstats -quiet -mean ${n4bias})
ImageMath 3 ${n4corrected} / ${n4input} ${n4bias}

mincresample -keep -near -like ${tmpdir}/N4.mnc ${tmpdir}/slabmask.mnc ${tmpdir}/slabmask2.mnc
minccalc -clobber -expression 'A[0]*A[1]' ${tmpdir}/N4.mnc ${tmpdir}/slabmask2.mnc ${tmpdir}/extracted.mnc

ImageMath 3 ${tmpdir}/extracted.mnc PadImage ${tmpdir}/extracted.mnc 50
mincresample -like ${tmpdir}/extracted.mnc -labels -keep -near ${tmpdir}/slabmask2.mnc ${tmpdir}/slabmask3.mnc

ExtractRegionFromImageByMask 3 ${tmpdir}/extracted.mnc ${output} ${tmpdir}/slabmask3.mnc 1 10



rm -rf ${tmpdir}
