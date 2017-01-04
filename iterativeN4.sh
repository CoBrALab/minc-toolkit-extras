#!/bin/bash
set -euo pipefail

BEASTLIBRARY_DIR="${QUARANTINE_PATH}/resources/BEaST_libraries/combined"
REGISTRATIONMODEL="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c.mnc"
REGISTRATIONBRAINMASK="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc"
WMPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_wm_tal_nlin_sym_09c.mnc"
GMPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_gm_tal_nlin_sym_09c.mnc"
CSFPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_csf_tal_nlin_sym_09c.mnc"

tmpdir=$(mktemp -d)

#Calculate shrink values
#Target is 4mm steps for first round, 2mm for second
dx=$(mincinfo -attvalue xspace:step $1)
dy=$(mincinfo -attvalue yspace:step $1)
dz=$(mincinfo -attvalue zspace:step $1)

shrinkround1=$(python -c "import math; print max(4,int(math.ceil(4 / ( ( abs($dx) + abs($dy) + abs($dz) ) / 3.0))))")
shrinkround2=$(python -c "import math; print max(2,int(math.ceil(2 / ( ( abs($dx) + abs($dy) + abs($dz) ) / 3.0))))")

maxval=$(mincstats -max -quiet $1)

#Generate a whole-image mask
minccalc -unsigned -byte -expression 'A[0]?1:1' $1 $tmpdir/initmask.mnc


#Very lightly truncate image intensity and re-normalize
ImageMath 3 $tmpdir/trunc.mnc TruncateImageIntensity $1 0.025 0.995 256
ImageMath 3 $tmpdir/truncnorm0.mnc Normalize $tmpdir/trunc.mnc
ImageMath 3 $tmpdir/trunc.mnc m $tmpdir/truncnorm0.mnc ${maxval}

rm $tmpdir/truncnorm0.mnc

#Correct entire image domain
N4BiasFieldCorrection -d 3 -s $shrinkround1 --verbose -x $tmpdir/initmask.mnc \
  -b [200] -c [300x300x300x300,0.0] --histogram-sharpening [0.05,0.01,200] -i $tmpdir/trunc.mnc -o $tmpdir/precorrect.mnc

#Normalize and rescale intensity
ImageMath 3 $tmpdir/norm0.mnc Normalize $tmpdir/precorrect.mnc
ImageMath 3 $tmpdir/precorrect.mnc m $tmpdir/norm0.mnc ${maxval}

rm $tmpdir/norm0.mnc

#First iteration
n=0

#Denoise for registration and beast
minc_anlm  --rician $tmpdir/precorrect.mnc $tmpdir/denoise${n}.mnc

#Register to MNI space
bestlinreg_g -nmi -lsq12 -target_mask ${REGISTRATIONBRAINMASK} $tmpdir/denoise${n}.mnc ${REGISTRATIONMODEL} $tmpdir/0_${n}.xfm $tmpdir/mni${n}.mnc

#Intensity normalize
volume_pol --order 1 --min 0 --max 100 --noclamp $tmpdir/mni${n}.mnc ${REGISTRATIONMODEL} --source_mask ${REGISTRATIONBRAINMASK} --target_mask ${REGISTRATIONBRAINMASK}  $tmpdir/mni${n}.norm.mnc

rm -f $tmpdir/mni${n}.mnc

#Run a quick beast to get a brain mask
mincbeast -clobber -verbose -fill -median -same_res -flip -conf ${BEASTLIBRARY_DIR}/default.4mm.conf ${BEASTLIBRARY_DIR} $tmpdir/mni${n}.norm.mnc $tmpdir/beastmask${n}.mnc

#Resample beast mask and MNI mask to native space
itk_resample  --labels --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm ${REGISTRATIONBRAINMASK} $tmpdir/mnimask${n}.mnc
itk_resample  --labels --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm $tmpdir/beastmask${n}.mnc $tmpdir/bmask${n}.mnc
#Combine the masks because sometimes beast misses badly biased cerebellum
mincmath  -or $tmpdir/mnimask${n}.mnc $tmpdir/bmask${n}.mnc $tmpdir/mask${n}.mnc

#Expand the mask a bit
mincmorph -3D26 -successive DDD $tmpdir/mask${n}.mnc $tmpdir/mask${n}_D.mnc

#Do first round of masked bias field correction, use brain mask as weight
N4BiasFieldCorrection -d 3 -s $shrinkround1 --verbose -w $tmpdir/mask${n}_D.mnc -x $tmpdir/initmask.mnc \
-b [200] -c [300x300x300x300,0.0] --histogram-sharpening [0.05,0.01,200] -i $tmpdir/trunc.mnc -o $tmpdir/round${n}.mnc

#Normalize and rescale intensity
ImageMath 3 $tmpdir/norm${n}.mnc Normalize $tmpdir/round${n}.mnc
ImageMath 3 $tmpdir/round${n}.mnc m $tmpdir/norm${n}.mnc ${maxval}

rm -f $tmpdir/norm${n}.mnc

n=1

#Denoise new output
minc_anlm  --rician $tmpdir/round$((${n} - 1)).mnc $tmpdir/denoise${n}.mnc

#Register to MNI space
bestlinreg_g -nmi -lsq12 -target_mask ${REGISTRATIONBRAINMASK} -source_mask $tmpdir/mask$((${n} -1 ))_D.mnc $tmpdir/denoise${n}.mnc ${REGISTRATIONMODEL} $tmpdir/0_${n}.xfm $tmpdir/mni${n}.mnc

xfminvert $tmpdir/0_${n}.xfm $tmpdir/0_${n}_inv.xfm

antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 -a 1 --verbose --minc \
  --output $tmpdir/nonlin${n} \
  --initial-moving-transform $tmpdir/0_${n}_inv.xfm \
  --transform SyN[0.25,3,0] --metric MI[${REGISTRATIONMODEL},$tmpdir/denoise${n}.mnc,1,32] --convergence [200x100x70x50,1e-6,10] --shrink-factors 10x8x4x2 --smoothing-sigmas 4x3x2x1mm --masks [${REGISTRATIONBRAINMASK},NULL]

#Intensity normalize
volume_pol --order 1 --min 0 --max 100 --noclamp $tmpdir/mni${n}.mnc ${REGISTRATIONMODEL} --source_mask $tmpdir/beastmask$((${n} - 1)).mnc --target_mask ${REGISTRATIONBRAINMASK}  $tmpdir/mni${n}.norm.mnc

rm -f $tmpdir/mni${n}.mnc

#Quick Beast
mincbeast  -verbose -fill -median -same_res -flip -conf ${BEASTLIBRARY_DIR}/default.4mm.conf ${BEASTLIBRARY_DIR} $tmpdir/mni${n}.norm.mnc $tmpdir/beastmask${n}.mnc

#Resample and join brain masks
#itk_resample  --labels --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm ${REGISTRATIONBRAINMASK} $tmpdir/mnimask${n}.mnc
itk_resample  --labels --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm $tmpdir/beastmask${n}.mnc $tmpdir/bmask${n}.mnc

antsApplyTransforms -d 3 -i ${REGISTRATIONBRAINMASK} -t $tmpdir/nonlin${n}_inverse.xfm -n NearestNeighbor -r $tmpdir/denoise${n}.mnc -o $tmpdir/mnimask${n}.mnc

mincmath  -or $tmpdir/mnimask${n}.mnc $tmpdir/bmask${n}.mnc $tmpdir/mask${n}.mnc

#Expand the mask a bit
mincmorph -3D26 -successive DD $tmpdir/mask${n}.mnc $tmpdir/mask${n}_D.mnc

#Resample MNI Priors to Native space for classification
antsApplyTransforms -i ${WMPRIOR} -t $tmpdir/nonlin${n}_inverse.xfm -r $tmpdir/denoise${n}.mnc -o $tmpdir/round${n}SegmentationPrior3.mnc -d 3 -n Linear
antsApplyTransforms -i ${GMPRIOR} -t $tmpdir/nonlin${n}_inverse.xfm -r $tmpdir/denoise${n}.mnc -o $tmpdir/round${n}SegmentationPrior2.mnc -d 3 -n Linear
antsApplyTransforms -i ${CSFPRIOR} -t $tmpdir/nonlin${n}_inverse.xfm -r $tmpdir/denoise${n}.mnc -o $tmpdir/round${n}SegmentationPrior1.mnc -d 3 -n Linear

#Do an initial classification using the MNI priors, remove outliers
Atropos -d 3 -x $tmpdir/mask${n}_D.mnc -c [5,0.0] -a $tmpdir/round$((${n} - 1)).mnc --verbose 1 -i  PriorProbabilityImages[3,$tmpdir/round${n}SegmentationPrior%d.mnc,0.25] -k Gaussian -m [0.1,1x1x1] \
  -o [$tmpdir/classify${n}.mnc,$tmpdir/round${n}SegmentationPosteriors%d.mnc] -r 1 -p Socrates[0] -w GrubbsRosner

rm -f $tmpdir/denoise${n}.mnc
rm -f $tmpdir/round${n}SegmentationPrior3.mnc $tmpdir/round${n}SegmentationPrior2.mnc $tmpdir/round${n}SegmentationPrior1.mnc
rm -f $tmpdir/nonlin*.xfm $tmpdir/nonlin*.mnc

#Combine GM and WM probably images into a N4 mask,
ImageMath 3 $tmpdir/weight${n}.mnc PureTissueN4WeightMask $tmpdir/round${n}SegmentationPosteriors2.mnc $tmpdir/round${n}SegmentationPosteriors3.mnc

#Perform bias field correction with weight mask
N4BiasFieldCorrection -d 3 -s $shrinkround2 --verbose -w $tmpdir/weight${n}.mnc -x $tmpdir/initmask.mnc \
  -b [200] -c [300x300x300x300,0.0] --histogram-sharpening [0.05,0.01,200] -i $tmpdir/trunc.mnc -o $tmpdir/round${n}.mnc

rm -f $tmpdir/weight${n}.mnc

#Normalize intensity
ImageMath 3 $tmpdir/norm${n}.mnc Normalize $tmpdir/round${n}.mnc
ImageMath 3 $tmpdir/round${n}.mnc m $tmpdir/norm${n}.mnc ${maxval}

rm -f $tmpdir/norm${n}.mnc

#Repeat until happy
for (( n = 2; n <= 3; n++ ))
do

  minc_anlm  --rician $tmpdir/round$((${n} - 1)).mnc $tmpdir/denoise${n}.mnc

  bestlinreg_g -nmi -lsq12 -target_mask ${REGISTRATIONBRAINMASK} -source_mask $tmpdir/mask$((${n} -1 ))_D.mnc $tmpdir/denoise${n}.mnc ${REGISTRATIONMODEL} $tmpdir/0_${n}.xfm $tmpdir/mni${n}.mnc

  volume_pol --order 1 --min 0 --max 100 --noclamp $tmpdir/mni${n}.mnc ${REGISTRATIONMODEL} --source_mask $tmpdir/beastmask$((${n} - 1)).mnc --target_mask ${REGISTRATIONBRAINMASK}  $tmpdir/mni${n}.norm.mnc

  rm -f $tmpdir/mni${n}.mnc

  mincbeast  -verbose -fill -median -same_res -flip -conf ${BEASTLIBRARY_DIR}/default.4mm.conf ${BEASTLIBRARY_DIR} $tmpdir/mni${n}.norm.mnc $tmpdir/beastmask${n}.mnc

  itk_resample  --labels --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm ${REGISTRATIONBRAINMASK} $tmpdir/mnimask${n}.mnc
  itk_resample  --labels --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm $tmpdir/beastmask${n}.mnc $tmpdir/bmask${n}.mnc
  mincmath -or $tmpdir/mnimask${n}.mnc $tmpdir/bmask${n}.mnc $tmpdir/mask${n}.mnc

  mincmorph -3D26 -successive DD $tmpdir/mask${n}.mnc $tmpdir/mask${n}_D.mnc

  Atropos -d 3 -x $tmpdir/mask${n}_D.mnc -c [5,0.0] -a $tmpdir/round$((${n} - 1)).mnc --verbose 1 -i PriorProbabilityImages[3,$tmpdir/round$((${n} -1 ))SegmentationPosteriors%d.mnc,0.25] \
    -k Gaussian -m [0.1,1x1x1] -o [$tmpdir/classify${n}.mnc,$tmpdir/round${n}SegmentationPosteriors%d.mnc] -r 1 -p Socrates[1]

  rm -f $tmpdir/denoise${n}.mnc

  rm -f $tmpdir/round$((${n} -1 ))SegmentationPosteriors?.mnc

  ImageMath 3 $tmpdir/weight${n}.mnc PureTissueN4WeightMask $tmpdir/round${n}SegmentationPosteriors2.mnc $tmpdir/round${n}SegmentationPosteriors3.mnc

  N4BiasFieldCorrection -d 3 -s $shrinkround2 --verbose -w $tmpdir/weight${n}.mnc -x $tmpdir/initmask.mnc \
    -b [200] -c [300x300x300x300,0.0] --histogram-sharpening [0.05,0.01,200] -i $tmpdir/trunc.mnc -o $tmpdir/round${n}.mnc

  rm -f $tmpdir/round$((${n} -1 )).mnc
  rm -f $tmpdir/weight${n}.mnc

  ImageMath 3 $tmpdir/norm${n}.mnc Normalize $tmpdir/round${n}.mnc
  ImageMath 3 $tmpdir/round${n}.mnc m $tmpdir/norm${n}.mnc ${maxval}

  rm -f $tmpdir/norm${n}.mnc

done

cp -f $tmpdir/round3.mnc $2
rm -rf $tmpdir
