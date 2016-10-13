#!/bin/bash
set -e
BEASTMODEL_DIR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2"
BEASTMODEL_NAME="mni_icbm152_t1_tal_nlin_sym_09c"
BEASTLIBRARY_DIR="${QUARANTINE_PATH}/resources/BEaST_libraries/combined"
REGISTRATIONMODEL="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c.mnc"
REGISTRATIONBRAINMASK="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc"
REGISTRATIONHEADMASK="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c_headmask.mnc"
REGISTRATIONWMMASK="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_wm_tal_nlin_sym_09c.mnc"
REGISTRATIONGMMASK="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_gm_tal_nlin_sym_09c.mnc"
WMPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_wm_tal_nlin_sym_09c.mnc"
GMPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_gm_tal_nlin_sym_09c.mnc"
CSFPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_csf_tal_nlin_sym_09c.mnc"

tmpdir=$(mktemp -d)

#Calculate shrink values
#Target is 4mm steps for first round, 2mm for second
dx=$(mincinfo -attvalue xspace:step $1)
dy=$(mincinfo -attvalue yspace:step $1)
dz=$(mincinfo -attvalue zspace:step $1)

shrinkround1=$(python -c "import math; print int(math.floor(4 / ( ( abs($dx) + abs($dy) + abs($dz) ) / 3.0)))")
shrinkround2=$(python -c "import math; print int(math.floor(2 / ( ( abs($dx) + abs($dy) + abs($dz) ) / 3.0)))")


#Generate a whole-image mask
minccalc -clobber -unsigned -byte -expression 'A[0]?1:1' $1 $tmpdir/initmask.mnc

#Bootstrap mask for first intensity normalization
cp ${REGISTRATIONBRAINMASK} $tmpdir/beastmask0.mnc

#Very lightly truncate image intensity and re-normalize
ImageMath 3 $tmpdir/trunc.mnc TruncateImageIntensity $1 0.025 0.995 256
ImageMath 3 $tmpdir/truncnorm0.mnc Normalize $tmpdir/trunc.mnc
ImageMath 3 $tmpdir/trunc.mnc m $tmpdir/truncnorm0.mnc 1000

#Correct entire image domain
N4BiasFieldCorrection -d 3 -s $shrinkround1 --verbose -x $tmpdir/initmask.mnc \
-b [200] -c [300x300x300x300,0.0] --histogram-sharpening [0.05,0.01,200] -i $tmpdir/trunc.mnc -o $tmpdir/precorrect.mnc

#Normalize and rescale intensity
ImageMath 3 $tmpdir/norm0.mnc Normalize $tmpdir/precorrect.mnc
ImageMath 3 $tmpdir/precorrect.mnc m $tmpdir/norm0.mnc 1000

#First iteration
n=0

#Denoise for registration and beast
minc_anlm --clobber --rician $tmpdir/precorrect.mnc $tmpdir/denoise${n}.mnc

#Register to MNI space
bestlinreg_g -nmi -lsq12 -target_mask ${REGISTRATIONBRAINMASK} $tmpdir/denoise${n}.mnc ${REGISTRATIONMODEL} $tmpdir/0_${n}.xfm $tmpdir/mni${n}.mnc -clobber

#Intensity normalize
volume_pol --order 1 --min 0 --max 100 --noclamp $tmpdir/mni${n}.mnc ${REGISTRATIONMODEL} --source_mask $tmpdir/beastmask0.mnc --target_mask ${REGISTRATIONBRAINMASK} --clobber $tmpdir/mni${n}.norm.mnc

#Run a quick beast to get a brain mask
mincbeast -clobber -verbose -fill -median -same_res -flip -conf /opt/quarantine/resources/BEaST_libraries/combined/default.4mm.conf /opt/quarantine/resources/BEaST_libraries/combined $tmpdir/mni${n}.norm.mnc $tmpdir/beastmask${n}.mnc

#Resample beast mask and MNI mask to native space
itk_resample --clobber --labels --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm ${REGISTRATIONBRAINMASK} $tmpdir/mnimask${n}.mnc
itk_resample --clobber --labels --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm $tmpdir/beastmask${n}.mnc $tmpdir/bmask${n}.mnc
#Combine the masks because sometimes beast misses badly biased cerebellum
mincmath -clobber -or $tmpdir/mnimask${n}.mnc $tmpdir/bmask${n}.mnc $tmpdir/mask${n}.mnc

#Expand the mask a bit
mincmorph -3D26 -successive DD $tmpdir/mask${n}.mnc $tmpdir/mask${n}_D.mnc -clobber

#Do first round of masked bias field correction, use brain mask as weight
N4BiasFieldCorrection -d 3 -s $shrinkround1 --verbose -w $tmpdir/mask${n}_D.mnc -x $tmpdir/initmask.mnc \
-b [200] -c [300x300x300x300,0.0] --histogram-sharpening [0.05,0.01,200] -i $tmpdir/trunc.mnc -o $tmpdir/round${n}.mnc

#Normalize and rescale intensity
ImageMath 3 $tmpdir/norm${n}.mnc Normalize $tmpdir/round${n}.mnc
ImageMath 3 $tmpdir/round${n}.mnc m $tmpdir/norm${n}.mnc 1000

n=1

#Denoise new output
minc_anlm --clobber --rician $tmpdir/round$((${n} - 1)).mnc $tmpdir/denoise${n}.mnc

#Register to MNI space
bestlinreg_g -nmi -lsq12 -target_mask ${REGISTRATIONBRAINMASK} $tmpdir/denoise${n}.mnc ${REGISTRATIONMODEL} $tmpdir/0_${n}.xfm $tmpdir/mni${n}.mnc -clobber

antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 -a 1 --verbose --minc \
  --output $tmpdir/nonlin${n} \
  --initial-moving-transform $tmpdir/0_${n}.xfm \
  --transform SyN[0.1,3,0] --metric MI[$tmpdir/denoise${n}.mnc,${REGISTRATIONMODEL},1,32] --convergence [100x70x50,1e-6,10] --shrink-factors 8x4x2 --smoothing-sigmas 3x2x1mm

#Intensity normalize
volume_pol --order 1 --min 0 --max 100 --noclamp $tmpdir/mni${n}.mnc ${REGISTRATIONMODEL} --source_mask $tmpdir/beastmask$((${n} - 1)).mnc --target_mask ${REGISTRATIONBRAINMASK} --clobber $tmpdir/mni${n}.norm.mnc

#Quick Beast
mincbeast -clobber -verbose -fill -median -same_res -flip -conf /opt/quarantine/resources/BEaST_libraries/combined/default.4mm.conf /opt/quarantine/resources/BEaST_libraries/combined $tmpdir/mni${n}.norm.mnc $tmpdir/beastmask${n}.mnc

#Resample and join brain masks
itk_resample --clobber --labels --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm ${REGISTRATIONBRAINMASK} $tmpdir/mnimask${n}.mnc
itk_resample --clobber --labels --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm $tmpdir/beastmask${n}.mnc $tmpdir/bmask${n}.mnc
mincmath -clobber -or $tmpdir/mnimask${n}.mnc $tmpdir/bmask${n}.mnc $tmpdir/mask${n}.mnc

#Expand the mask a bit
mincmorph -3D26 -successive DD $tmpdir/mask${n}.mnc $tmpdir/mask${n}_D.mnc -clobber

#Resample MNI Priors to Native space for classification
antsApplyTransforms -i ${WMPRIOR} -t $tmpdir/nonlin${n}.xfm -r $tmpdir/denoise${n}.mnc -o $tmpdir/round${n}SegmentationPrior3.mnc -d 3 -n Linear
antsApplyTransforms -i ${GMPRIOR} -t $tmpdir/nonlin${n}.xfm -r $tmpdir/denoise${n}.mnc -o $tmpdir/round${n}SegmentationPrior2.mnc -d 3 -n Linear
antsApplyTransforms -i ${CSFPRIOR} -t $tmpdir/nonlin${n}.xfm -r $tmpdir/denoise${n}.mnc -o $tmpdir/round${n}SegmentationPrior1.mnc -d 3 -n Linear

#itk_resample --clobber --order 1 --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm ${WMPRIOR} $tmpdir/round${n}SegmentationPrior3.mnc
#itk_resample --clobber --order 1 --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm ${GMPRIOR} $tmpdir/round${n}SegmentationPrior2.mnc
#itk_resample --clobber --order 1 --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm ${CSFPRIOR} $tmpdir/round${n}SegmentationPrior1.mnc

#Do an initial classification using the MNI priors, remove outliers
Atropos -d 3 -x $tmpdir/mask${n}_D.mnc -c [5,0.0] -a $tmpdir/round$((${n} - 1)).mnc --verbose 1 -i  PriorProbabilityImages[3,$tmpdir/round${n}SegmentationPrior%d.mnc,0.25] -k Gaussian -m [0.1,1x1x1] \
-o [$tmpdir/classify${n}.mnc,$tmpdir/round${n}SegmentationPosteriors%d.mnc] -r 1 -p Socrates[0] -w GrubbsRosner

#Combine GM and WM probably images into a N4 mask,
ImageMath 3 $tmpdir/weight${n}.mnc PureTissueN4WeightMask $tmpdir/round${n}SegmentationPosteriors2.mnc $tmpdir/round${n}SegmentationPosteriors3.mnc

#Perform bias field correction with weight mask
N4BiasFieldCorrection -d 3 -s $shrinkround2 --verbose -w $tmpdir/weight${n}.mnc -x $tmpdir/initmask.mnc \
-b [200] -c [300x300x300x300,0.0] --histogram-sharpening [0.05,0.01,200] -i $tmpdir/trunc.mnc -o $tmpdir/round${n}.mnc

#Normalize intensity
ImageMath 3 $tmpdir/norm${n}.mnc Normalize $tmpdir/round${n}.mnc
ImageMath 3 $tmpdir/round${n}.mnc m $tmpdir/norm${n}.mnc 1000

#Repeat until happy
for (( n = 2; n <= 3; n++ ))
do

  minc_anlm --clobber --rician $tmpdir/round$((${n} - 1)).mnc $tmpdir/denoise${n}.mnc

  bestlinreg_g -nmi -lsq12 -target_mask ${REGISTRATIONBRAINMASK} $tmpdir/denoise${n}.mnc ${REGISTRATIONMODEL} $tmpdir/0_${n}.xfm $tmpdir/mni${n}.mnc -clobber

  volume_pol --order 1 --min 0 --max 100 --noclamp $tmpdir/mni${n}.mnc ${REGISTRATIONMODEL} --source_mask $tmpdir/beastmask$((${n} - 1)).mnc --target_mask ${REGISTRATIONBRAINMASK} --clobber $tmpdir/mni${n}.norm.mnc

  mincbeast -clobber -verbose -fill -median -same_res -flip -conf /opt/quarantine/resources/BEaST_libraries/combined/default.4mm.conf /opt/quarantine/resources/BEaST_libraries/combined $tmpdir/mni${n}.norm.mnc $tmpdir/beastmask${n}.mnc

  itk_resample --clobber --labels --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm ${REGISTRATIONBRAINMASK} $tmpdir/mnimask${n}.mnc
  itk_resample --clobber --labels --like $tmpdir/denoise${n}.mnc --invert_transform --transform $tmpdir/0_${n}.xfm $tmpdir/beastmask${n}.mnc $tmpdir/bmask${n}.mnc
  mincmath -or $tmpdir/mnimask${n}.mnc $tmpdir/bmask${n}.mnc $tmpdir/mask${n}.mnc -clobber

  mincmorph -3D26 -successive DD $tmpdir/mask${n}.mnc $tmpdir/mask${n}_D.mnc -clobber

  Atropos -d 3 -x $tmpdir/mask${n}_D.mnc -c [5,0.0] -a $tmpdir/round$((${n} - 1)).mnc --verbose 1 -i PriorProbabilityImages[3,$tmpdir/round$((${n} -1 ))SegmentationPosteriors%d.mnc,0.25] \
  -k Gaussian -m [0.1,1x1x1] -o [$tmpdir/classify${n}.mnc,$tmpdir/round${n}SegmentationPosteriors%d.mnc] -r 1 -p Socrates[1]

  ImageMath 3 $tmpdir/weight${n}.mnc PureTissueN4WeightMask $tmpdir/round${n}SegmentationPosteriors2.mnc $tmpdir/round${n}SegmentationPosteriors3.mnc

  N4BiasFieldCorrection -d 3 -s $shrinkround2 --verbose -w $tmpdir/weight${n}.mnc -x $tmpdir/initmask.mnc \
  -b [200] -c [300x300x300x300,0.0] --histogram-sharpening [0.05,0.01,200] -i $tmpdir/trunc.mnc -o $tmpdir/round${n}.mnc

  ImageMath 3 $tmpdir/norm${n}.mnc Normalize $tmpdir/round${n}.mnc
  ImageMath 3 $tmpdir/round${n}.mnc m $tmpdir/norm${n}.mnc 1000

done

#mincstats -quiet -histogram $(dirname $2)/$(basename $2 .mnc).hist -mask $tmpdir/mask3_D.mnc -mask_binvalue 1 $tmpdir/round3.mnc

cp -f $tmpdir/round3.mnc $2
rm -rf $tmpdir
