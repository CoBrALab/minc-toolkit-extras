#!/bin/bash
#Preprocessing script for rat brains, using Fischer344 template and MINC files converted from DICOM from Paravision 5.x

#If using conversion from PV6, might need to remove the volflip -y command


set -euo pipefail
set -x

tmpdir=$(mktemp -d)

input=$1
output=$2

#Clamp negatives
mincmath -clamp -const2 0 $(mincstats -max -quiet ${input}) ${input} ${tmpdir}/clamp.mnc

#Reorient to MINC RAS
volmash -swap zy ${tmpdir}/clamp.mnc ${tmpdir}/mash.mnc
volflip -x ${tmpdir}/mash.mnc ${tmpdir}/flip.mnc

#Center FOV
volcentre -zero_dircos -com ${tmpdir}/flip.mnc ${tmpdir}/centre.mnc

#Denoise
minc_anlm --rician --mt $(nproc) ${tmpdir}/centre.mnc ${tmpdir}/denoise.mnc

#Foreground/background threshold
ImageMath 3 ${tmpdir}/weight.mnc ThresholdAtMean ${tmpdir}/denoise.mnc 0.5
ImageMath 3 ${tmpdir}/weight.mnc GetLargestComponent ${tmpdir}/weight.mnc

#Make a crop mask
ImageMath 3 ${tmpdir}/cropmask.mnc FillHoles ${tmpdir}/weight.mnc 2
ImageMath 3 ${tmpdir}/denoise.mnc PadImage ${tmpdir}/denoise.mnc 50
antsApplyTransforms -d 3 -i ${tmpdir}/cropmask.mnc -r ${tmpdir}/denoise.mnc -n GenericLabel --verbose -o ${tmpdir}/cropmask.mnc

#Extract data with cropmask
ImageMath 3 ${tmpdir}/denoise.mnc m ${tmpdir}/denoise.mnc ${tmpdir}/cropmask.mnc
ExtractRegionFromImageByMask 3 ${tmpdir}/denoise.mnc ${tmpdir}/denoise.crop.mnc ${tmpdir}/cropmask.mnc 1 5
mv -f ${tmpdir}/denoise.crop.mnc ${tmpdir}/denoise.mnc

minccalc -unsigned -byte -expression '1' ${tmpdir}/denoise.mnc ${tmpdir}/initmask.mnc
antsApplyTransforms -d 3 -i ${tmpdir}/weight.mnc -r ${tmpdir}/denoise.mnc -n GenericLabel --verbose -o ${tmpdir}/weight.mnc

#Round 1 bias field correction
N4BiasFieldCorrection -d 3 -i ${tmpdir}/denoise.mnc -b [30] -c [300x300x300,1e-5] -r 0 -w ${tmpdir}/weight.mnc -x ${tmpdir}/initmask.mnc \
-o ${tmpdir}/N4.mnc -s 4 --verbose

#Calculate an otsu mask for multiplying brainmask
ThresholdImage 3 ${tmpdir}/N4.mnc ${tmpdir}/otsu.mnc Otsu 1

#Model for brainmask
fixedfile=${QUARANTINE_PATH}/resources/Fischer344/Fischer344_template.mnc
movingfile=${tmpdir}/N4.mnc
fixedmask=${QUARANTINE_PATH}/resources/Fischer344/Fischer344_mask.mnc
movingmask=NOMASK

#Optimized multi-stage affine registration
antsRegistration --dimensionality 3 --verbose --minc \
  --output [ ${tmpdir}/reg ] \
  --use-histogram-matching 0 \
  --initial-moving-transform [ ${fixedfile},${movingfile},1 ] \
--transform Translation[ 0.5 ] \
  --metric Mattes[ ${fixedfile},${movingfile},1,32,Regular,0.5 ] \
  --convergence [ 2025x2025x2025,1e-6,10 ] \
  --shrink-factors 8x8x4 \
  --smoothing-sigmas 0.849321800288x0.424660900144x0.212330450072mm \
  --masks [ NOMASK,NOMASK ] \
--transform Rigid[ 0.5 ] \
  --metric Mattes[ ${fixedfile},${movingfile},1,64,Regular,0.5 ] \
  --convergence [ 2025x2025x986,1e-6,10 ] \
  --shrink-factors 8x4x2 \
  --smoothing-sigmas 0.424660900144x0.212330450072x0.106165225036mm \
  --masks [ NOMASK,NOMASK ] \
--transform Similarity[ 0.25 ] \
  --metric Mattes[ ${fixedfile},${movingfile},1,128,Regular,0.75 ] \
  --convergence [ 2025x986x314,1e-6,10 ] \
  --shrink-factors 4x2x1 \
  --smoothing-sigmas 0.212330450072x0.106165225036x0.053082612518mm \
  --masks [ NOMASK,NOMASK ] \
--transform Affine[ 0.1 ] \
  --metric Mattes[ ${fixedfile},${movingfile},1,256,None ] \
  --convergence [ 986x314x177x25,1e-6,10 ] \
  --shrink-factors 2x1x1x1 \
  --smoothing-sigmas 0.106165225036x0.053082612518x0.026541306259x0mm \
  --masks [ ${fixedmask},NOMASK ]

#Resample mask
antsApplyTransforms -d 3 -i ${fixedmask} -r ${movingfile} -t [${tmpdir}/reg0_GenericAffine.xfm,1] \
-n GenericLabel --verbose -o ${tmpdir}/mask.mnc

#Zero out non-brain tissue in mask using otsu
ImageMath 3 ${tmpdir}/weight.mnc m ${tmpdir}/otsu.mnc ${tmpdir}/mask.mnc

#Redo bias field correction
N4BiasFieldCorrection -d 3 -i ${tmpdir}/denoise.mnc -b [30] -c [300x300x300,1e-5] -r 1 -w ${tmpdir}/weight.mnc -x ${tmpdir}/initmask.mnc \
-o ${tmpdir}/N4.mnc -s 4 --verbose

volcentre -zero_dircos -com ${tmpdir}/N4.mnc ${output}

rm -rf ${tmpdir}
