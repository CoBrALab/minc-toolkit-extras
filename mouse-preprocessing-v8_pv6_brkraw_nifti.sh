#!/bin/bash

set -euo pipefail
set -x

calc() { awk "BEGIN{ print $* }"; }

function pad_recrop_image() {
  local input=$1
  local padding=$2
  local output=$3

  local localtmpdir=$(mktemp -d)

  local res=$(calc "(sqrt(($(mincinfo -attvalue xspace:step ${input}))**2) + sqrt(($(mincinfo -attvalue yspace:step ${input}))**2) + sqrt(($(mincinfo -attvalue zspace:step ${input}))**2))/3")
  volpad -noauto -distance $(calc "int(5/${res})") ${input} ${localtmpdir}/padded.mnc

  ThresholdImage 3 ${localtmpdir}/padded.mnc ${localtmpdir}/bgmask.mnc $(mincstats -quiet -pctT 0.1 -hist_floor 1e-12 ${input}) Inf 1 0
  ThresholdImage 3 ${localtmpdir}/padded.mnc ${localtmpdir}/otsu.mnc Otsu 4 ${localtmpdir}/bgmask.mnc
  ThresholdImage 3 ${localtmpdir}/otsu.mnc ${localtmpdir}/otsu.mnc 2 Inf 1 0

  ExtractRegionFromImageByMask 3 ${localtmpdir}/padded.mnc ${output} ${localtmpdir}/otsu.mnc 1 $(calc "int(${padding}/${res})")

  rm -rf ${localtmpdir}
}

function isotropize_downsample() {
  # Need smoothing for downsampling to avoid aliasing
  # Ideas stolen from https://discourse.itk.org/t/resampling-to-isotropic-signal-processing-theory/1403
  local input=$1
  local outputres=$2
  local output=$3
  local inputres=$(python -c "print('\n'.join([str(abs(x)) for x in [float(x) for x in \"$(PrintHeader ${input} 1)\".split(\"x\")]]))")
  local blurs=""

  local tmpdir=$(mktemp -d)

  for dim in ${inputres}; do
    if [[ $(python -c "print(${dim}>(${outputres}-1e-6))") == True ]]; then
      #Special casing for zero/negative blurs
      blurs+=1e-12x
    else
      blurs+=$(python -c "import math; print(math.sqrt((${outputres}**2.0 - ${dim}**2.0)/(2.0*math.sqrt(2.0*math.log(2.0)))**2.0))")x
    fi
  done

  SmoothImage 3 ${input} "${blurs%?}" ${tmpdir}/smoothed.mnc 1 0
  ResampleImage 3 ${tmpdir}/smoothed.mnc ${tmpdir}/isotropized.mnc ${outputres}x${outputres}x${outputres} 0 4

  mincmath -clobber -quiet -unsigned -short -clamp -const2 0 65535 ${tmpdir}/isotropized.mnc ${output}

  rm -rf ${tmpdir}
}

function weighted_n4() {

  local input=$1
  local weight=$2
  local output=$3
  local biasout=$4

  local tmpdir=$(mktemp -d)

  minc_anlm --short --mt $(nproc) ${input} ${tmpdir}/denoise.mnc

  # Downsample the image
  isotropize_downsample ${tmpdir}/denoise.mnc ${_downsample_res} ${tmpdir}/downsample.mnc

  # Produce a foreground mask for use in N4
  minccalc -unsigned -byte -expression '1' ${tmpdir}/downsample.mnc ${tmpdir}/fullmask.mnc

  # Resample the weight image to the downsampled image
  antsApplyTransforms -d 3 -i ${weight} -o ${tmpdir}/weight_downsample.mnc \
    -n GenericLabel -r ${tmpdir}/downsample.mnc --verbose

  N4BiasFieldCorrection -d 3 -i ${tmpdir}/downsample.mnc -b [ 40 ] \
    -c [ 50x50x50x50x50,0.00000001 ] -r 1 -s 1 \
    -x ${tmpdir}/fullmask.mnc -w ${tmpdir}/weight_downsample.mnc \
    -o [ ${tmpdir}/corrected.mnc,${tmpdir}/bias.mnc ] --verbose

  # Upsample the bias field image
  antsApplyTransforms -d 3 -i ${tmpdir}/bias.mnc \
    -r ${input} \
    -n BSpline[5] \
    --default-value 1 \
    -o ${tmpdir}/bias_upsample.mnc

  local biasmean=$(mincstats -quiet -mean -mask ${weight} -mask_binvalue 1 ${tmpdir}/bias_upsample.mnc)

  minccalc -unsigned -short -expression "clamp(A[0]/(A[1]/${biasmean}),0,65535)" ${input} \
    ${tmpdir}/bias_upsample.mnc ${tmpdir}/corrected.mnc -clobber

  minccalc -clobber -expression "A[0]/${biasmean}" ${tmpdir}/bias_upsample.mnc ${biasout}

  # Denoise the image
  minc_anlm --short --mt $(nproc) ${tmpdir}/corrected.mnc ${tmpdir}/corrected_denoise.mnc

  # Renormalize the image and output
  mincnorm -clobber -clamp -cutoff 0 -short -out_floor 0 -out_ceil 65535 ${tmpdir}/corrected_denoise.mnc ${output}

  rm -rf ${tmpdir}
}

function threshold_image_to_mask() {

  local input=$1
  local output=$2
  local tmpdir=$(mktemp -d)

  # Construct a background mask (just above zero)
  ThresholdImage 3 ${input} ${tmpdir}/bgmask.mnc $(mincstats -quiet -pctT 0.1 -hist_floor 1e-12 ${input}) Inf 1 0
  # Do Otsu threshold in the mask
  ThresholdImage 3 ${input} ${tmpdir}/otsu.mnc Otsu 4 ${tmpdir}/bgmask.mnc
  # Convert the mask into a foreground mask
  ThresholdImage 3 ${tmpdir}/otsu.mnc ${tmpdir}/weight.mnc 2 Inf 1 0
  ImageMath 3 ${output} GetLargestComponent ${tmpdir}/weight.mnc

  rm -rf ${tmpdir}
}

tmpdir=$(mktemp -d)

mkdir -p ${tmpdir}

input=$1
output=$2

# Fixed downsample resolution
_downsample_res=0.2

#REGMODEL=${QUARANTINE_PATH}/resources/in-vivo-MEMRI_90micron/recrop/Dorr_2008_on_MEMRI_C57BL6_P43_average.mnc
#REGMASK=${QUARANTINE_PATH}/resources/in-vivo-MEMRI_90micron/recrop/Dorr_2008_on_MEMRI_C57BL6_P43_mask.mnc

REGMODEL=${QUARANTINE_PATH}/resources/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/ex-vivo/N4/DSURQE_40micron.mnc
REGMASK=${QUARANTINE_PATH}/resources/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/ex-vivo/N4/DSURQE_40micron_mask.mnc

minres=$(python -c "print(max([str(abs(x)) for x in [float(x) for x in \"$(PrintHeader ${input} 1)\".split(\"x\")]]))")
isotropize_downsample ${REGMODEL} ${minres} ${tmpdir}/model.mnc
REGMODEL=${tmpdir}/model.mnc
minccalc -expression "A[0]>$(mincstats -quiet -pctT 1 -floor 1e-12 ${REGMODEL})?1:0" -unsigned -byte ${tmpdir}/model.mnc ${tmpdir}/model_fov.mnc

mincresample -keep -near -label -like ${REGMODEL} ${REGMASK} ${tmpdir}/model_mask.mnc
REGMASK=${tmpdir}/model_mask.mnc

# Fix sometimes malformed niftis
ConvertImage 3 ${input} ${tmpdir}/input.nii.gz
nii2mnc ${tmpdir}/input.nii.gz ${tmpdir}/input.mnc

input=${tmpdir}/input.mnc

# Change the slicing direction to be the same as what ITK outputs by default
mincreshape -dimorder zspace,yspace,xspace +direction ${input} ${tmpdir}/reshape1.mnc
mincreshape -dimorder zspace,yspace,xspace +direction ${tmpdir}/reshape1.mnc ${tmpdir}/reshape2.mnc
rm -f ${tmpdir}/reshape1.mnc

# Store the direction cosine transform because we're going to zero out the file
dircos_to_xfm ${tmpdir}/reshape2.mnc ${tmpdir}/transform_to_input.xfm

#Square up direction cosines
cp -f ${tmpdir}/reshape2.mnc ${tmpdir}/standardized_cosines.mnc
minc_modify_header -dinsert xspace:direction_cosines=1,0,0 ${tmpdir}/standardized_cosines.mnc
minc_modify_header -dinsert yspace:direction_cosines=0,1,0 ${tmpdir}/standardized_cosines.mnc
minc_modify_header -dinsert zspace:direction_cosines=0,0,1 ${tmpdir}/standardized_cosines.mnc

input=${tmpdir}/standardized_cosines.mnc

# Use a histogram cutoff to trim values, but don't clamp
mincnorm -double -cutoff 0.001 ${input} ${tmpdir}/trim.mnc
# Clamp below zero after
mincmath -clamp -const2 0 inf ${tmpdir}/trim.mnc ${tmpdir}/trim_clamp.mnc
# This time renormalize 0-65535 with no histogram cutoff
mincnorm -clamp -cutoff 0 -short -out_floor 0 -out_ceil 65535 ${tmpdir}/trim_clamp.mnc ${tmpdir}/norm.mnc

pad_recrop_image ${tmpdir}/norm.mnc 2 ${tmpdir}/norm.mnc

threshold_image_to_mask ${tmpdir}/norm.mnc ${tmpdir}/weight1.mnc

weighted_n4 ${tmpdir}/norm.mnc ${tmpdir}/weight1.mnc ${tmpdir}/corrected_denoise_renorm.mnc ${tmpdir}/bias1.mnc

antsRegistration_affine_SyN.sh --float \
  --skip-nonlinear \
  --fixed-mask ${REGMASK} \
  ${tmpdir}/corrected_denoise_renorm.mnc \
  ${REGMODEL} \
  ${tmpdir}/to_model_

antsApplyTransforms -d 3 --verbose -n GenericLabel \
  -i ${REGMASK}  \
  -t [ ${tmpdir}/to_model_0_GenericAffine.xfm,1 ] \
  -r ${tmpdir}/corrected_denoise_renorm.mnc \
  -o ${tmpdir}/affine_mask.mnc -n GenericLabel

ImageMath 3 ${tmpdir}/precorrect.mnc / ${tmpdir}/norm.mnc ${tmpdir}/bias1.mnc

weighted_n4 ${tmpdir}/precorrect.mnc ${tmpdir}/affine_mask.mnc ${tmpdir}/corrected_denoise_renorm.mnc ${tmpdir}/bias2.mnc

antsRegistration_affine_SyN.sh --float --clobber --close \
  --initial-transform ${tmpdir}/to_model_0_GenericAffine.xfm \
  --skip-nonlinear \
  --fixed-mask ${REGMASK} \
  --moving-mask ${tmpdir}/affine_mask.mnc \
  ${tmpdir}/corrected_denoise_renorm.mnc \
  ${REGMODEL} \
  ${tmpdir}/to_model_

antsApplyTransforms -d 3 --verbose -n GenericLabel \
  -i ${tmpdir}/model_fov.mnc \
  -r ${tmpdir}/corrected_denoise_renorm.mnc \
  -t [ ${tmpdir}/to_model_0_GenericAffine.xfm,1 ] \
  -o ${tmpdir}/subject_fov.mnc

antsApplyTransforms -d 3 --verbose -n GenericLabel \
  -i ${REGMASK}  \
  -t [ ${tmpdir}/to_model_0_GenericAffine.xfm,1 ] \
  -r ${tmpdir}/corrected_denoise_renorm.mnc \
  -o ${tmpdir}/affine_mask.mnc -n GenericLabel

ImageMath 3 ${tmpdir}/corrected_denoise_renorm.mnc m ${tmpdir}/corrected_denoise_renorm.mnc ${tmpdir}/subject_fov.mnc

antsRegistration_affine_SyN.sh --float --clobber \
  --skip-linear \
  --fast \
  --fixed-mask ${REGMASK} \
  --moving-mask ${tmpdir}/affine_mask.mnc \
  --initial-transform ${tmpdir}/to_model_0_GenericAffine.xfm \
  ${tmpdir}/corrected_denoise_renorm.mnc \
  ${REGMODEL} \
  ${tmpdir}/to_model_

antsApplyTransforms -d 3 --verbose -n GenericLabel \
  -i ${REGMASK} \
  -t [ ${tmpdir}/to_model_0_GenericAffine.xfm,1 ] \
  -t ${tmpdir}/to_model_1_inverse_NL.xfm \
  -r ${tmpdir}/corrected_denoise_renorm.mnc \
  -o ${tmpdir}/mask.mnc

weighted_n4 ${tmpdir}/precorrect.mnc ${tmpdir}/mask.mnc ${tmpdir}/corrected_denoise_renorm.mnc ${tmpdir}/bias3.mnc

ImageMath 3 ${tmpdir}/corrected_denoise_renorm.mnc m ${tmpdir}/corrected_denoise_renorm.mnc ${tmpdir}/subject_fov.mnc

pad_recrop_image ${tmpdir}/corrected_denoise_renorm.mnc 2 ${tmpdir}/corrected_denoise_renorm.mnc

mincresample -clobber -tfm_input_sampling -transform ${tmpdir}/transform_to_input.xfm ${tmpdir}/corrected_denoise_renorm.mnc ${output}
mincresample -keep -near -label -transform ${tmpdir}/transform_to_input.xfm -like ${output} ${tmpdir}/mask.mnc $(dirname ${output})/$(basename ${output} .mnc)_mask.mnc

rm -rf ${tmpdir}
