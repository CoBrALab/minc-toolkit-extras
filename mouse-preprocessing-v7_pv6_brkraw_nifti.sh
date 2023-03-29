#!/bin/bash
#Preprocessing script for mouse brains, using Dorr-steadman template and MINC files converted from DICOM from Paravision 5.x

#If using conversion from PV6, might need to remove the volflip -y command


set -euo pipefail
set -x

calc(){ awk "BEGIN { print $* }" ;}

tmpdir=$(mktemp -d)

input=$1
output=$2

origdistance=20
distance=${origdistance}
levels=4
cycles=3
iters=25
lambda=2e-6
shrink=4.0
fwhm=0.1
stop=1e-5

do_correct() {
  distance=${origdistance}
  j=0
  while (( j < levels )); do
    i=0
    while (( i < cycles )); do
      nu_correct -clobber -normalize_field \
        -stop ${stop} -distance ${distance} -iterations ${iters} -fwhm ${fwhm} -shrink ${shrink} -lambda ${lambda} \
        -mask ${tmpdir}/weight.mnc ${n3input} ${tmpdir}/corrected_${distance}_${i}.mnc
      n3input=${tmpdir}/corrected_${distance}_${i}.mnc
      ((++i))
    done
    distance=$(calc  ${distance}/2)
    ((++j))
  done

}

modelfile=${QUARANTINE_PATH}/resources/in-vivo-MEMRI_90micron/recrop/Dorr_2008_on_MEMRI_C57BL6_P43_average.mnc
modelmask=${QUARANTINE_PATH}/resources/in-vivo-MEMRI_90micron/recrop/Dorr_2008_on_MEMRI_C57BL6_P43_mask.mnc

ConvertImage 3 ${input} ${tmpdir}/input.nii.gz
nii2mnc ${tmpdir}/input.nii.gz ${tmpdir}/input.mnc

input=${tmpdir}/input.mnc

min_res=$(python -c "print(min([abs(x) for x in [float(x) for x in \"$(PrintHeader ${input} 1)\".split(\"x\")]]))")

# Clamp range to avoid negative numbers, rescale to 0-65535
mincnorm -noclamp -short -out_floor 0 -out_ceil 65535 ${input} ${tmpdir}/input2.mnc
mincmath -clamp -const2 0 inf ${tmpdir}/input2.mnc ${tmpdir}/norm.mnc
mv -f ${tmpdir}/norm.mnc ${tmpdir}/input.mnc


# Pad image for processing
volpad -noauto -distance $(calc "int(5/${min_res})") ${tmpdir}/input.mnc ${tmpdir}/pad.mnc
mv -f ${tmpdir}/pad.mnc ${tmpdir}/input.mnc

cp -f ${tmpdir}/input.mnc ${tmpdir}/originput.mnc


#Construct Otsu Mask of entire image
input=${tmpdir}/originput.mnc

ThresholdImage 3 ${input} ${tmpdir}/bgmask.mnc 1e-12 Inf 1 0
ThresholdImage 3 ${input} ${tmpdir}/weight.mnc Otsu 4 ${tmpdir}/bgmask.mnc
ThresholdImage 3 ${tmpdir}/weight.mnc ${tmpdir}/weight.mnc 2 Inf 1 0
iMath 3 ${tmpdir}/weight.mnc ME ${tmpdir}/weight.mnc 1 1 ball 1
ImageMath 3 ${tmpdir}/weight.mnc GetLargestComponent ${tmpdir}/weight.mnc
iMath 3 ${tmpdir}/weight.mnc MD ${tmpdir}/weight.mnc 1 1 ball 1
ExtractRegionFromImageByMask 3 ${input} ${tmpdir}/repad.mnc ${tmpdir}/weight.mnc 1 $(calc "int(2/${min_res})")
ExtractRegionFromImageByMask 3 ${tmpdir}/originput.mnc ${tmpdir}/repad_orig.mnc ${tmpdir}/weight.mnc 1 $(calc "int(2/${min_res})")
mv -f ${tmpdir}/repad_orig.mnc ${tmpdir}/originput.mnc
mv -f ${tmpdir}/repad.mnc ${input}
ThresholdImage 3 ${input} ${tmpdir}/bgmask.mnc 1e-12 Inf 1 0
mincresample -like ${input} -keep -near -labels ${tmpdir}/weight.mnc ${tmpdir}/weight2.mnc
mv -f ${tmpdir}/weight2.mnc ${tmpdir}/weight.mnc

minc_anlm --clobber --mt $(nproc) ${input} ${tmpdir}/denoise.mnc

n3input=${tmpdir}/denoise.mnc
do_correct

ThresholdImage 3 ${n3input} ${tmpdir}/weight.mnc Otsu 4 ${tmpdir}/bgmask.mnc
ThresholdImage 3 ${tmpdir}/weight.mnc ${tmpdir}/weight.mnc 2 Inf 1 0
iMath 3 ${tmpdir}/weight.mnc ME ${tmpdir}/weight.mnc 1 1 ball 1
ImageMath 3 ${tmpdir}/weight.mnc GetLargestComponent ${tmpdir}/weight.mnc
iMath 3 ${tmpdir}/weight.mnc MD ${tmpdir}/weight.mnc 1 1 ball 1
mincresample -like ${input} -keep -near -labels ${tmpdir}/weight.mnc ${tmpdir}/weight2.mnc
mv -f ${tmpdir}/weight2.mnc ${tmpdir}/weight.mnc

n3input=${tmpdir}/denoise.mnc
do_correct

for file in ${tmpdir}/*imp; do
  echo nu_evaluate -clobber -mapping ${file} -mask ${tmpdir}/weight.mnc -field ${tmpdir}/$(basename $file .imp)_field.mnc ${input} ${tmpdir}/$(basename $file .imp).mnc
done | parallel

mincmath -clobber -mult ${tmpdir}/*field.mnc ${tmpdir}/precorrect_field_combined.mnc
mincmath -clobber -copy_header -zero -div ${tmpdir}/originput.mnc ${tmpdir}/precorrect_field_combined.mnc ${tmpdir}/precorrect.mnc

minc_anlm --clobber --mt $(nproc) ${tmpdir}/precorrect.mnc ${tmpdir}/denoise.mnc
n3input=${tmpdir}/denoise.mnc

antsRegistration_affine_SyN.sh --verbose --fast --fixed-mask ${modelmask} \
  ${n3input} ${modelfile} ${tmpdir}/tomodel

antsApplyTransforms -d 3 -i ${modelmask} -t [ ${tmpdir}/tomodel0_GenericAffine.xfm,1 ] -t ${tmpdir}/tomodel1_inverse_NL.xfm \
  -o ${tmpdir}/newmask.mnc -n GenericLabel -r ${n3input} --verbose

iMath 3 ${tmpdir}/newmask.mnc MD ${tmpdir}/newmask.mnc 1 1 ball 1
ThresholdImage 3 ${n3input} ${tmpdir}/weight.mnc Otsu 4 ${tmpdir}/newmask.mnc
ThresholdImage 3 ${tmpdir}/weight.mnc ${tmpdir}/weight.mnc 2 Inf 1 0
ImageMath 3 ${tmpdir}/weight.mnc m ${tmpdir}/newmask.mnc ${tmpdir}/weight.mnc
iMath 3 ${tmpdir}/weight.mnc ME ${tmpdir}/weight.mnc 1 1 ball 1
ImageMath 3 ${tmpdir}/weight.mnc GetLargestComponent ${tmpdir}/weight.mnc
iMath 3 ${tmpdir}/weight.mnc MD ${tmpdir}/weight.mnc 1 1 ball 1
mincresample -like ${input} -keep -near -labels ${tmpdir}/weight.mnc ${tmpdir}/weight2.mnc
mv -f ${tmpdir}/weight2.mnc ${tmpdir}/weight.mnc

do_correct

for file in ${tmpdir}/*imp; do
  echo nu_evaluate -clobber -mapping ${file} -mask ${tmpdir}/weight.mnc -field ${tmpdir}/$(basename $file .imp)_field.mnc ${input} ${tmpdir}/$(basename $file .imp).mnc
done | parallel

mincmath -clobber -mult ${tmpdir}/*field.mnc ${tmpdir}/precorrect_field_combined.mnc ${tmpdir}/field_final.mnc
mincmath -clobber -copy_header -zero -div ${tmpdir}/originput.mnc ${tmpdir}/field_final.mnc ${tmpdir}/correct.mnc

minc_anlm --mt $(nproc) ${tmpdir}/correct.mnc ${tmpdir}/denoise_correct.mnc

ExtractRegionFromImageByMask 3 ${tmpdir}/denoise_correct.mnc ${tmpdir}/recrop.mnc ${tmpdir}/weight.mnc 1 $(calc "int(2/${min_res})")
mv -f ${tmpdir}/recrop.mnc ${tmpdir}/denoise_correct.mnc
mincresample -keep -near -like ${tmpdir}/denoise_correct.mnc ${tmpdir}/weight.mnc ${tmpdir}/weight2.mnc
mv -f ${tmpdir}/weight2.mnc ${tmpdir}/weight.mnc

cp ${tmpdir}/weight.mnc $(dirname ${output})/$(basename ${output} .mnc)_mask.mnc
cp ${tmpdir}/denoise_correct.mnc ${output}

xfminvert ${tmpdir}/tomodel0_GenericAffine.xfm ${tmpdir}/tomodel0_GenericAffine_invert.xfm
param2xfm $(xfm2param ${tmpdir}/tomodel0_GenericAffine_invert.xfm | grep -E 'scale|shear') ${tmpdir}/scaleshear.xfm
xfminvert ${tmpdir}/scaleshear.xfm ${tmpdir}/unscaleshear.xfm
xfmconcat ${tmpdir}/tomodel0_GenericAffine_invert.xfm ${tmpdir}/unscaleshear.xfm ${tmpdir}/lsq6.xfm

mincresample -tfm_input_sampling -transform ${tmpdir}/lsq6.xfm ${tmpdir}/denoise_correct.mnc ${tmpdir}/lsq6.mnc

mincmath -clamp -const2 0 $(mincstats -quiet -max ${tmpdir}/lsq6.mnc) ${tmpdir}/lsq6.mnc $(dirname ${output})/$(basename ${output} .mnc)_lsq6.mnc

mincresample -transform ${tmpdir}/lsq6.xfm  -like $(dirname ${output})/$(basename ${output} .mnc)_lsq6.mnc -keep -near -labels ${tmpdir}/weight.mnc $(dirname ${output})/$(basename ${output} .mnc)_lsq6_mask.mnc

rm -rf ${tmpdir}
