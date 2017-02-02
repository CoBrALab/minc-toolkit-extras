#!/bin/bash
set -euo pipefail

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${THREADS_PER_COMMAND:-$(nproc)}

BEASTLIBRARY_DIR="${QUARANTINE_PATH}/resources/BEaST_libraries/combined"
REGISTRATIONMODEL="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c.mnc"
REGISTRATIONBRAINMASK="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc"
WMPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_wm_tal_nlin_sym_09c.mnc"
GMPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_gm_tal_nlin_sym_09c.mnc"
CSFPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_csf_tal_nlin_sym_09c.mnc"

input=$1
output=$2

tmpdir=$(mktemp -d)

#Calculate shrink values
#Target is 4mm steps for first round, 2mm for second
dx=$(mincinfo -attvalue xspace:step ${input})
dy=$(mincinfo -attvalue yspace:step ${input})
dz=$(mincinfo -attvalue zspace:step ${input})

shrinkround1=$(python -c "import math; print max(4,int(math.ceil(4 / ( ( abs($dx) + abs($dy) + abs($dz) ) / 3.0))))")
shrinkround2=$(python -c "import math; print max(2,int(math.ceil(2 / ( ( abs($dx) + abs($dy) + abs($dz) ) / 3.0))))")

maxval=$(mincstats -max -quiet ${input})

set -x

#Generate a whole-image mask
minccalc -quiet -unsigned -byte -expression 'A[0]?1:1' ${input} ${tmpdir}/initmask.mnc

################################################################################################################################################################################################################################
#Round 0, N4 with Otsu Mask
################################################################################################################################################################################################################################
n=0
mkdir -p ${tmpdir}/${n}

#Very lightly truncate image intensity and re-normalize
#minccalc -quiet -expression "if(A[0]<${bottom}){out=${bottom};}else{if(A[0]>${top}){out=${top};}else{out=A[0];}}" ${input} ${tmpdir}/trunc.mnc
ImageMath 3 ${tmpdir}/trunc.mnc TruncateImageIntensity ${input} 0.025 0.9995 256
ImageMath 3 ${tmpdir}/trunc.mnc Normalize ${tmpdir}/trunc.mnc
ImageMath 3 ${tmpdir}/trunc.mnc m ${tmpdir}/trunc.mnc ${maxval}

ThresholdImage 3 ${tmpdir}/trunc.mnc ${tmpdir}/${n}/headmask.mnc Otsu 1
ImageMath 3 ${tmpdir}/${n}/weight.mnc GetLargestComponent ${tmpdir}/${n}/headmask.mnc


#Correct entire image domain
N4BiasFieldCorrection -d 3 -s $shrinkround1 -x ${tmpdir}/initmask.mnc -w ${tmpdir}/${n}/weight.mnc \
-b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${input} -o [${tmpdir}/${n}/corrected.mnc,${tmpdir}/bias${n}.mnc] -r 0

ImageMath 3 ${tmpdir}/${n}/corrected.mnc TruncateImageIntensity ${tmpdir}/${n}/corrected.mnc 0.025 0.9995 256

#Normalize and rescale intensity
ImageMath 3 ${tmpdir}/${n}/corrected.mnc Normalize ${tmpdir}/${n}/corrected.mnc
ImageMath 3 ${tmpdir}/${n}/corrected.mnc m ${tmpdir}/${n}/corrected.mnc ${maxval}

################################################################################################################################################################################################################################
#Round 1, N4 with brain mask intersected with Otsu mask
################################################################################################################################################################################################################################
((++n))
mkdir -p ${tmpdir}/${n}

#Denoise for registration and beast
minc_anlm --mt $ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS --rician ${tmpdir}/$((${n} - 1))/corrected.mnc ${tmpdir}/${n}/denoise.mnc

rm -rf ${tmpdir}/$((${n} - 1))

#Register to MNI space
bestlinreg_g -noverbose -nmi -lsq12 -target_mask ${REGISTRATIONBRAINMASK} ${tmpdir}/${n}/denoise.mnc ${REGISTRATIONMODEL} ${tmpdir}/${n}/linear.xfm ${tmpdir}/${n}/mni.mnc

#Intensity normalize
volume_pol --order 1 --min 0 --max 100 --noclamp ${tmpdir}/${n}/mni.mnc ${REGISTRATIONMODEL} --source_mask ${REGISTRATIONBRAINMASK} --target_mask ${REGISTRATIONBRAINMASK} ${tmpdir}/${n}/mni.norm.mnc

#Run a quick beast to get a brain mask
mincbeast -clobber -fill -median -same_res -flip -conf ${BEASTLIBRARY_DIR}/default.4mm.conf ${BEASTLIBRARY_DIR} ${tmpdir}/${n}/mni.norm.mnc ${tmpdir}/${n}/beastmask.mnc

#Resample beast mask and MNI mask to native space
itk_resample  --labels --like ${tmpdir}/${n}/denoise.mnc --invert_transform --transform ${tmpdir}/${n}/linear.xfm ${REGISTRATIONBRAINMASK} ${tmpdir}/${n}/mnimask.mnc
itk_resample  --labels --like ${tmpdir}/${n}/denoise.mnc --invert_transform --transform ${tmpdir}/${n}/linear.xfm ${tmpdir}/${n}/beastmask.mnc ${tmpdir}/${n}/bmask.mnc

#Combine the masks because sometimes beast misses badly biased cerebellum
mincmath -byte -unsigned -or ${tmpdir}/${n}/mnimask.mnc ${tmpdir}/${n}/bmask.mnc ${tmpdir}/${n}/mask.mnc

#Expand the mask a bit
mincmorph -successive D ${tmpdir}/${n}/mask.mnc ${tmpdir}/${n}/mask_D.mnc

ThresholdImage 3 ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/weight.mnc Otsu 1 ${tmpdir}/${n}/mask_D.mnc


#Do first round of masked bias field correction, use brain mask as weight
N4BiasFieldCorrection -d 3 -s $shrinkround1 -w ${tmpdir}/${n}/weight.mnc -x ${tmpdir}/initmask.mnc \
-b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${input} -o [${tmpdir}/${n}/corrected.mnc,${tmpdir}/bias${n}.mnc] -r 0

ImageMath 3 ${tmpdir}/${n}/corrected.mnc TruncateImageIntensity ${tmpdir}/${n}/corrected.mnc 0.025 0.9995 256 ${tmpdir}/${n}/mask.mnc

#Normalize and rescale intensity
ImageMath 3 ${tmpdir}/${n}/corrected.mnc Normalize ${tmpdir}/${n}/corrected.mnc
ImageMath 3 ${tmpdir}/${n}/corrected.mnc m ${tmpdir}/${n}/corrected.mnc ${maxval}

itk_similarity -msq ${tmpdir}/bias$((${n} - 1)).mnc ${tmpdir}/bias${n}.mnc > ${tmpdir}/convergence.txt

################################################################################################################################################################################################################################
#Round 2, N4 with nonlinearly MNI-bootstrapped WM/GM segmentation proabilities
################################################################################################################################################################################################################################
((++n))
mkdir -p ${tmpdir}/${n}

#Denoise for registration and beast
minc_anlm --mt $ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS --rician ${tmpdir}/$((${n} - 1))/corrected.mnc ${tmpdir}/${n}/denoise.mnc

#Register to MNI space
bestlinreg_g -init_xfm ${tmpdir}/$((${n} - 1))/linear.xfm -close -noverbose -nmi -lsq12 -target_mask ${REGISTRATIONBRAINMASK} -source_mask ${tmpdir}/$((${n} - 1))/mask.mnc ${tmpdir}/${n}/denoise.mnc ${REGISTRATIONMODEL} ${tmpdir}/${n}/linear.xfm ${tmpdir}/${n}/mni.mnc

#Intensity normalize
volume_pol --order 1 --min 0 --max 100 --noclamp ${tmpdir}/${n}/mni.mnc ${REGISTRATIONMODEL} --source_mask ${tmpdir}/$((${n} - 1))/beastmask.mnc --target_mask ${REGISTRATIONBRAINMASK} ${tmpdir}/${n}/mni.norm.mnc

rm -rf ${tmpdir}/$((${n} - 1))

#Run a quick beast to get a brain mask
mincbeast -clobber -fill -median -same_res -flip -conf ${BEASTLIBRARY_DIR}/default.4mm.conf ${BEASTLIBRARY_DIR} ${tmpdir}/${n}/mni.norm.mnc ${tmpdir}/${n}/beastmask.mnc

mincmath -byte -unsigned -or ${tmpdir}/${n}/beastmask.mnc ${REGISTRATIONBRAINMASK} ${tmpdir}/${n}/mnicombinedmask.mnc

minc_nuyl ${tmpdir}/${n}/mni.mnc ${REGISTRATIONMODEL} --source-mask ${tmpdir}/${n}/mnicombinedmask.mnc --target-mask ${REGISTRATIONBRAINMASK} ${tmpdir}/${n}/mni.nuyl.mnc --clobber

antsRegistration --dimensionality 3 --float 1 --collapse-output-transforms 1 -a 0 --minc -u 0 \
--output [${tmpdir}/${n}/nonlin,${tmpdir}/${n}/mninonlin.mnc] \
--transform SyN[0.25,3,0] --metric Mattes[${REGISTRATIONMODEL},${tmpdir}/${n}/mni.nuyl.mnc,1,64] \
--convergence [200x100x70x50,1e-6,10] --shrink-factors 8x4x4x2 --smoothing-sigmas 3.3973x1.6986x0.8493x0.4247mm --masks [${REGISTRATIONBRAINMASK},NULL]

#Resample MNI Priors to Native space for classification
antsApplyTransforms -i ${WMPRIOR} -t ${tmpdir}/${n}/linear.xfm -t ${tmpdir}/${n}/nonlin0_inverse_NL.xfm -r ${tmpdir}/${n}/denoise.mnc -o ${tmpdir}/${n}/SegmentationPrior3.mnc -d 3 -n Linear
antsApplyTransforms -i ${GMPRIOR} -t ${tmpdir}/${n}/linear.xfm -t ${tmpdir}/${n}/nonlin0_inverse_NL.xfm -r ${tmpdir}/${n}/denoise.mnc -o ${tmpdir}/${n}/SegmentationPrior2.mnc -d 3 -n Linear
antsApplyTransforms -i ${CSFPRIOR} -t ${tmpdir}/${n}/linear.xfm -t ${tmpdir}/${n}/nonlin0_inverse_NL.xfm -r ${tmpdir}/${n}/denoise.mnc -o ${tmpdir}/${n}/SegmentationPrior1.mnc -d 3 -n Linear
#Masks
antsApplyTransforms -i ${REGISTRATIONBRAINMASK} -t ${tmpdir}/${n}/linear.xfm -r ${tmpdir}/${n}/denoise.mnc -o ${tmpdir}/${n}/mnimask.mnc -d 3 -n NearestNeighbor
antsApplyTransforms -i ${tmpdir}/${n}/beastmask.mnc -t ${tmpdir}/${n}/linear.xfm -r ${tmpdir}/${n}/denoise.mnc -o ${tmpdir}/${n}/bmask.mnc -d 3 -n NearestNeighbor

#Combine the masks because sometimes beast misses badly biased cerebellum
mincmath -byte -unsigned -or ${tmpdir}/${n}/mnimask.mnc ${tmpdir}/${n}/bmask.mnc ${tmpdir}/${n}/mask.mnc

#Expand the mask a bit
mincmorph -successive D ${tmpdir}/${n}/mask.mnc ${tmpdir}/${n}/mask_D.mnc

#Do an initial classification using the MNI priors, remove outliers
Atropos -d 3 -x ${tmpdir}/${n}/mask_D.mnc -c [25,0.0] -a ${tmpdir}/${n}/denoise.mnc -i PriorProbabilityImages[3,${tmpdir}/${n}/SegmentationPrior%d.mnc,0] -k Gaussian -m [0.1,1x1x1] \
-o [${tmpdir}/${n}/classify.mnc,${tmpdir}/${n}/SegmentationPosteriors%d.mnc] -r 1 -p Socrates[0,1,0.1,0.1] -w GrubbsRosner

#Combine GM and WM probably images into a N4 mask,
ImageMath 3 ${tmpdir}/${n}/weight.mnc PureTissueN4WeightMask ${tmpdir}/${n}/SegmentationPosteriors2.mnc ${tmpdir}/${n}/SegmentationPosteriors3.mnc

#Perform bias field correction with weight mask
N4BiasFieldCorrection -d 3 -s $shrinkround1 -w ${tmpdir}/${n}/weight.mnc -x ${tmpdir}/initmask.mnc \
-b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${input} -o [${tmpdir}/${n}/corrected.mnc,${tmpdir}/bias${n}.mnc] -r 0

ImageMath 3 ${tmpdir}/${n}/corrected.mnc TruncateImageIntensity ${tmpdir}/${n}/corrected.mnc 0.025 0.9995 256 ${tmpdir}/${n}/mask.mnc

#Normalize and rescale intensity
ImageMath 3 ${tmpdir}/${n}/corrected.mnc Normalize ${tmpdir}/${n}/corrected.mnc
ImageMath 3 ${tmpdir}/${n}/corrected.mnc m ${tmpdir}/${n}/corrected.mnc ${maxval}

itk_similarity -msq ${tmpdir}/bias$((${n} - 1)).mnc ${tmpdir}/bias${n}.mnc >> ${tmpdir}/convergence.txt

################################################################################################################################################################################################################################
#Remaining rounds, N4 with segmentations bootstrapped from prior run
#Stop when max iterations or when normalized sum of squares of difference of bias
#Fields less than 1e-4
################################################################################################################################################################################################################################
while true
do
  ((++n))
  mkdir -p ${tmpdir}/${n}

  #Denoise for registration and beast
  minc_anlm --mt $ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS --rician ${tmpdir}/$((${n} - 1))/corrected.mnc ${tmpdir}/${n}/denoise.mnc

  #Register to MNI space
  bestlinreg_g -init_xfm ${tmpdir}/$((${n} - 1))/linear.xfm -close -noverbose -nmi -lsq12 -target_mask ${REGISTRATIONBRAINMASK} -source_mask ${tmpdir}/$((${n} - 1))/mask.mnc ${tmpdir}/${n}/denoise.mnc ${REGISTRATIONMODEL} ${tmpdir}/${n}/linear.xfm ${tmpdir}/${n}/mni.mnc

  #Intensity normalize
  volume_pol --order 1 --min 0 --max 100 --noclamp ${tmpdir}/${n}/mni.mnc ${REGISTRATIONMODEL} --source_mask ${tmpdir}/$((${n} - 1))/beastmask.mnc --target_mask ${REGISTRATIONBRAINMASK} ${tmpdir}/${n}/mni.norm.mnc

  #Run a quick beast to get a brain mask
  mincbeast -clobber -fill -median -same_res -flip -conf ${BEASTLIBRARY_DIR}/default.2mm.conf ${BEASTLIBRARY_DIR} ${tmpdir}/${n}/mni.norm.mnc ${tmpdir}/${n}/beastmask.mnc

  #Resample beast mask and MNI mask to native space
  itk_resample  --labels --like ${tmpdir}/${n}/denoise.mnc --invert_transform --transform ${tmpdir}/${n}/linear.xfm ${REGISTRATIONBRAINMASK} ${tmpdir}/${n}/mnimask.mnc
  itk_resample  --labels --like ${tmpdir}/${n}/denoise.mnc --invert_transform --transform ${tmpdir}/${n}/linear.xfm ${tmpdir}/${n}/beastmask.mnc ${tmpdir}/${n}/bmask.mnc

  #Combine the masks because sometimes beast misses badly biased cerebellum
  mincmath -unsigned -byte -or ${tmpdir}/${n}/mnimask.mnc ${tmpdir}/${n}/bmask.mnc ${tmpdir}/${n}/mask.mnc

  #Do an initial classification using the MNI priors, remove outliers
  Atropos -d 3 -x ${tmpdir}/${n}/mask.mnc -c [25,0.0] -a ${tmpdir}/${n}/denoise.mnc -i PriorProbabilityImages[3,${tmpdir}/$((${n} - 1))/SegmentationPosteriors%d.mnc,0.25] -k Gaussian -m [0.1,1x1x1] \
  -o [${tmpdir}/${n}/classify.mnc,${tmpdir}/${n}/SegmentationPosteriors%d.mnc] -r 1 -p Socrates[1] -w GrubbsRosner

  rm -rf ${tmpdir}/$((${n} - 1))

  #Combine GM and WM probably images into a N4 mask,
  ImageMath 3 ${tmpdir}/${n}/weight.mnc PureTissueN4WeightMask ${tmpdir}/${n}/SegmentationPosteriors2.mnc ${tmpdir}/${n}/SegmentationPosteriors3.mnc

  #Perform bias field correction with weight mask
  N4BiasFieldCorrection -d 3 -s $shrinkround2 -w ${tmpdir}/${n}/weight.mnc -x ${tmpdir}/initmask.mnc \
  -b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${input} -o [${tmpdir}/${n}/corrected.mnc,${tmpdir}/bias${n}.mnc] -r 0

  ImageMath 3 ${tmpdir}/${n}/corrected.mnc TruncateImageIntensity ${tmpdir}/${n}/corrected.mnc 0.025 0.9995 256 ${tmpdir}/${n}/mask.mnc

  #Normalize and rescale intensity
  ImageMath 3 ${tmpdir}/${n}/corrected.mnc Normalize ${tmpdir}/${n}/corrected.mnc
  ImageMath 3 ${tmpdir}/${n}/corrected.mnc m ${tmpdir}/${n}/corrected.mnc ${maxval}

  itk_similarity -msq ${tmpdir}/bias$((${n} - 1)).mnc ${tmpdir}/bias${n}.mnc >> ${tmpdir}/convergence.txt

  [[ ( ${n} -lt 4 ) || ( (${n} -lt 10) && ( $(python -c "print $(itk_similarity -msq ${tmpdir}/bias$((${n} - 1)).mnc ${tmpdir}/bias${n}.mnc) > 0.0001") == "True" ) ) ]] || break

done

cp -f ${tmpdir}/${n}/corrected.mnc ${output}
#cp -f ${tmpdir}/${n}/mask.mnc $(dirname $output)/$(basename $output .mnc).mask.mnc
#cp -f ${tmpdir}/${n}/classify.mnc $(dirname $output)/$(basename $output .mnc).classify.mnc

rm -rf ${tmpdir}
