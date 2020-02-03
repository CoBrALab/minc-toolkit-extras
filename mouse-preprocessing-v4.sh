#!/bin/bash
#Preprocessing script for mouse brains, using Dorr-steadman template and MINC files converted from DICOM from Paravision 5.x

#If using conversion from PV6, might need to remove the volflip -y command


set -euo pipefail
set -x

calc(){ awk "BEGIN { print "$*" }"; }

tmpdir=$(mktemp -d)

input=$1
output=$2

cp -f ${input} ${tmpdir}/input.mnc
input=${tmpdir}/input.mnc

minimum_resolution=$(python -c "print(min([abs(x) for x in [float(x) for x in \"$(PrintHeader ${input} 1)\".split(\"x\")]]))")

minc_modify_header $input -sinsert :history=''

#Clamp negatives
mincmath -clamp -const2 0 $(mincstats -max -quiet ${input}) ${input} ${tmpdir}/clamp.mnc

#Reorient to MINC RAS
volmash -swap zy ${tmpdir}/clamp.mnc ${tmpdir}/mash.mnc
volflip -x ${tmpdir}/mash.mnc ${tmpdir}/flip.mnc

#Center FOV
volcentre -zero_dircos ${tmpdir}/flip.mnc ${tmpdir}/centre.mnc

cp ${tmpdir}/centre.mnc ${tmpdir}/centre_orig.mnc
minc_anlm --mt $(nproc) ${tmpdir}/centre.mnc ${tmpdir}/denoise.mnc
mv -f ${tmpdir}/denoise.mnc ${tmpdir}/centre.mnc

ImageMath 3 ${tmpdir}/weight.mnc ThresholdAtMean ${tmpdir}/centre.mnc 0.5
iMath 3 ${tmpdir}/weight.mnc ME ${tmpdir}/weight.mnc 2 1 ball 1
ImageMath 3 ${tmpdir}/weight.mnc GetLargestComponent ${tmpdir}/weight.mnc
iMath 3 ${tmpdir}/weight.mnc MD ${tmpdir}/weight.mnc 2 1 ball 1

minccalc -clobber -unsigned -byte -expression '1' ${tmpdir}/centre.mnc ${tmpdir}/initmask.mnc

N4BiasFieldCorrection -d 3 -i ${tmpdir}/centre.mnc -o ${tmpdir}/N4.mnc -s $(calc "int(4*0.1/${minimum_resolution}+0.5)") -b [ 20 ] -c [ 1000x1000x1000,1e-6 ] -w ${tmpdir}/weight.mnc -x ${tmpdir}/initmask.mnc --verbose

minccalc -clobber -unsigned -byte \
  -expression "A[0]>0.5*$(mincstats -quiet -biModalT -floor $(mincstats -quiet -pctT 0.1 ${tmpdir}/N4.mnc) -ceil $(mincstats -quiet -pctT 99.9 ${tmpdir}/N4.mnc) ${tmpdir}/N4.mnc)?1:0" \
  ${tmpdir}/N4.mnc ${tmpdir}/weight.mnc
iMath 3 ${tmpdir}/weight.mnc ME ${tmpdir}/weight.mnc 2 1 ball 1
ImageMath 3 ${tmpdir}/weight.mnc GetLargestComponent ${tmpdir}/weight.mnc
iMath 3 ${tmpdir}/weight.mnc MD ${tmpdir}/weight.mnc 2 1 ball 1
mincresample -clobber -like ${tmpdir}/N4.mnc ${tmpdir}/weight.mnc ${tmpdir}/weight2.mnc
mv -f ${tmpdir}/weight2.mnc ${tmpdir}/weight.mnc
minccalc -clobber -unsigned -byte -expression "A[0]&&A[1]<$(mincstats -quiet -mask ${tmpdir}/weight.mnc -mask_binvalue 1 -pctT 99.9 ${tmpdir}/N4.mnc)" ${tmpdir}/weight.mnc ${tmpdir}/N4.mnc ${tmpdir}/weight2.mnc
mv -f ${tmpdir}/weight2.mnc ${tmpdir}/weight.mnc

minccalc -clobber -unsigned -byte -expression '1' ${tmpdir}/centre.mnc ${tmpdir}/initmask.mnc

N4BiasFieldCorrection -d 3 -i ${tmpdir}/centre.mnc -o ${tmpdir}/N4.mnc -s $(calc "int(4*0.1/${minimum_resolution}+0.5)") -b [ 20 ] -c [ 1000x1000x1000,1e-6 ] -w ${tmpdir}/weight.mnc -x ${tmpdir}/initmask.mnc --verbose

minccalc -clobber -unsigned -byte \
  -expression "A[0]>0.5*$(mincstats -quiet -biModalT -floor $(mincstats -quiet -pctT 0.1 ${tmpdir}/N4.mnc) -ceil $(mincstats -quiet -pctT 99.9 ${tmpdir}/N4.mnc) ${tmpdir}/N4.mnc)?1:0" \
  ${tmpdir}/N4.mnc ${tmpdir}/weight.mnc
iMath 3 ${tmpdir}/weight.mnc ME ${tmpdir}/weight.mnc 2 1 ball 1
ImageMath 3 ${tmpdir}/weight.mnc GetLargestComponent ${tmpdir}/weight.mnc
iMath 3 ${tmpdir}/weight.mnc MD ${tmpdir}/weight.mnc 2 1 ball 1
minccalc -clobber -unsigned -byte -expression "A[0]&&A[1]<$(mincstats -quiet -mask ${tmpdir}/weight.mnc -mask_binvalue 1 -pctT 99.9 ${tmpdir}/N4.mnc)" ${tmpdir}/weight.mnc ${tmpdir}/N4.mnc ${tmpdir}/weight2.mnc
mv -f ${tmpdir}/weight2.mnc ${tmpdir}/weight.mnc

#Clip image for registration
minccalc -expression "clamp(A[0],0,($(mincstats -quiet -median ${tmpdir}/N4.mnc) + 1.5*($(mincstats -quiet -pctT 75 ${tmpdir}/N4.mnc) - $(mincstats -quiet -pctT 25 ${tmpdir}/N4.mnc))))" \
  ${tmpdir}/N4.mnc ${tmpdir}/clampreg.mnc


#Model for brainmask
fixedfile=${QUARANTINE_PATH}/resources/in-vivo-MEMRI_90micron/Dorr_2008_on_MEMRI_C57BL6_P43_average_recenter_upsample_N4.mnc
movingfile=${tmpdir}/clampreg.mnc
fixedmask=${QUARANTINE_PATH}/resources/in-vivo-MEMRI_90micron/Dorr_2008_on_MEMRI_C57BL6_P43_mask_recenter_upsample.mnc
movingmask=NOMASK

#Optimized multi-stage affine registration
antsRegistration --dimensionality 3 --verbose --minc \
--use-histogram-matching 0 \
--output [ ${tmpdir}/reg ] \
--transform Rigid[ 0.1 ] \
  --metric Mattes[ ${fixedfile},${movingfile},1,32,None ] \
  --convergence [ 2025x2025x2025x2025x2025x2025x2025x2025x2025x2025x2025x2025x2025x2025x2025x2025x675x225x75,1e-6,10 ] \
  --shrink-factors 17x17x17x17x16x15x14x13x12x11x10x9x8x7x6x5x4x3x2 \
  --smoothing-sigmas 0.67945744023x0.645484568219x0.611511696207x0.577538824196x0.543565952184x0.509593080173x0.475620208161x0.44164733615x0.407674464138x0.373701592127x0.339728720115x0.305755848104x0.271782976092x0.237810104081x0.203837232069x0.169864360058x0.135891488046x0.101918616035x0.067945744023mm \
  --masks [ NOMASK,NOMASK ] \
--transform Similarity[ 0.1 ] \
  --metric Mattes[ ${fixedfile},${movingfile},1,32,None ] \
  --convergence [ 2025x2025x2025x2025x2025x2025x675x225x75,1e-6,10 ] \
  --shrink-factors 10x9x8x7x6x5x4x3x2 \
  --smoothing-sigmas 0.339728720115x0.305755848104x0.271782976092x0.237810104081x0.203837232069x0.169864360058x0.135891488046x0.101918616035x0.067945744023mm \
  --masks [ NOMASK,NOMASK ] \
--transform Similarity[ 0.1 ] \
  --metric Mattes[ ${fixedfile},${movingfile},1,32,None ] \
  --convergence [ 2025x2025x2025x2025x2025x2025x675x225x75,1e-6,10 ] \
  --shrink-factors 10x9x8x7x6x5x4x3x2 \
  --smoothing-sigmas 0.339728720115x0.305755848104x0.271782976092x0.237810104081x0.203837232069x0.169864360058x0.135891488046x0.101918616035x0.067945744023mm \
  --masks [ ${fixedmask},${movingmask} ] \
--transform Affine[ 0.1 ] \
  --metric Mattes[ ${fixedfile},${movingfile},1,51,None ] \
  --convergence [ 2025x675x225x75x25x25,1e-6,10 ] \
  --shrink-factors 5x4x3x2x1x1 \
  --smoothing-sigmas 0.169864360058x0.135891488046x0.101918616035x0.067945744023x0.0339728720115x0mm \
  --masks [ ${fixedmask},${movingmask} ]
  
#Resample mask
antsApplyTransforms -d 3 -i ${fixedmask} -r ${movingfile} -t [${tmpdir}/reg0_GenericAffine.xfm,1] \
  -n GenericLabel --verbose -o ${tmpdir}/mask.mnc

ImageMath 3 ${tmpdir}/weight.mnc m ${tmpdir}/weight.mnc ${tmpdir}/mask.mnc
ImageMath 3 ${tmpdir}/weight.mnc GetLargestComponent ${tmpdir}/weight.mnc
iMath 3 ${tmpdir}/weight.mnc ME ${tmpdir}/weight.mnc 2 1 ball 1
ImageMath 3 ${tmpdir}/weight.mnc GetLargestComponent ${tmpdir}/weight.mnc
iMath 3 ${tmpdir}/weight.mnc MD ${tmpdir}/weight.mnc 2 1 ball 1

#Redo bias field correction
N4BiasFieldCorrection -d 3 -i ${tmpdir}/centre.mnc -b [ 20 ] -c [ 300x300x300,1e-6 ] -r 0 -w ${tmpdir}/weight.mnc -x ${tmpdir}/initmask.mnc \
  -o [ ${tmpdir}/N4.mnc, ${tmpdir}/bias.mnc ] -s $(calc "int(2*0.1/${minimum_resolution}+0.5)") --verbose --histogram-sharpening [ 0.05,0.01,256 ]

ImageMath 3 ${tmpdir}/bias.mnc / ${tmpdir}/bias.mnc $(mincstats -quiet -mean -mask ${tmpdir}/mask.mnc -mask_binvalue 1 ${tmpdir}/bias.mnc)

ImageMath 3 ${tmpdir}/N4.mnc / ${tmpdir}/centre_orig.mnc ${tmpdir}/bias.mnc

volcentre -zero_dircos -centre $(mincstats -com -quiet -world_only ${tmpdir}/weight.mnc) ${tmpdir}/N4.mnc $(dirname ${output})/$(basename ${output} .mnc)_uncropped.mnc

ImageMath 3 ${tmpdir}/N4.mnc PadImage ${tmpdir}/N4.mnc 50
antsApplyTransforms -d 3 -i ${tmpdir}/weight.mnc -o ${tmpdir}/cropmask.mnc -r ${tmpdir}/N4.mnc --verbose

ExtractRegionFromImageByMask 3 ${tmpdir}/N4.mnc ${tmpdir}/N4.crop.mnc ${tmpdir}/cropmask.mnc 1 10

iMath 3 ${tmpdir}/weight.mnc MD ${tmpdir}/weight.mnc 1 1 ball 1
ImageMath 3 ${tmpdir}/weight.mnc FillHoles ${tmpdir}/weight.mnc 2
iMath 3 ${tmpdir}/weight.mnc ME ${tmpdir}/weight.mnc 1 1 ball 1

antsApplyTransforms -d 3 -i ${tmpdir}/weight.mnc -r ${tmpdir}/N4.crop.mnc -o ${tmpdir}/weight.mnc

volcentre -zero_dircos -centre $(mincstats -com -quiet -world_only ${tmpdir}/weight.mnc) ${tmpdir}/N4.crop.mnc ${output}
volcentre -zero_dircos -com ${tmpdir}/weight.mnc $(dirname ${output})/$(basename ${output} .mnc)_mask.mnc

param2xfm $(xfm2param ${tmpdir}/reg0_GenericAffine.xfm | grep -E 'scale|shear') ${tmpdir}/scaleshear.xfm
xfminvert ${tmpdir}/scaleshear.xfm ${tmpdir}/unscaleshear.xfm
xfmconcat ${tmpdir}/reg0_GenericAffine.xfm ${tmpdir}/unscaleshear.xfm ${tmpdir}/lsq6.xfm
xfminvert ${tmpdir}/lsq6.xfm ${tmpdir}/lsq6_invert.xfm

mincresample -tfm_input_sampling -transform ${tmpdir}/lsq6_invert.xfm ${tmpdir}/N4.crop.mnc ${tmpdir}/lsq6.mnc

mincmath -clamp -const2 0 $(mincstats -quiet -max ${tmpdir}/lsq6.mnc) ${tmpdir}/lsq6.mnc $(dirname ${output})/$(basename ${output} .mnc)_lsq6.mnc

mincresample -transform ${tmpdir}/lsq6_invert.xfm  -like $(dirname ${output})/$(basename ${output} .mnc)_lsq6.mnc -keep -near -labels ${tmpdir}/weight.mnc $(dirname ${output})/$(basename ${output} .mnc)_lsq6_mask.mnc

#antsApplyTransforms -d 3 -i ${tmpdir}/N4.mnc -r ${fixedfile} -t ${tmpdir}/reg0_GenericAffine.xfm \
#  -n BSpline[5] --verbose -o ${tmpdir}/DSUQE.mnc

#mincmath -clamp -const2 0 $(mincstats -max -quiet ${tmpdir}/DSUQE.mnc) ${tmpdir}/DSUQE.mnc $(dirname ${output})/$(basename ${output} .mnc)_DSUQE.mnc



rm -rf ${tmpdir}
