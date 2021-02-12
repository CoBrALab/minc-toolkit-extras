#!/bin/bash
#Script to do optimal preprocessing on in-vivo/ex-vivo structural scans
#Taken using the CIC Bruker 7T
#usage:
#mouse-preprocessing-v3.sh input.mnc output.mnc

#Operations
# registers to DSURQE atlas
# gets approximate brain mask from atlas
# Bias field correction with N4

set -euo pipefail

REGTARGET=${QUARANTINE_PATH}/resources/P7_ECox_2015/N4/P7_ECox_average.mnc
REGMASK=${QUARANTINE_PATH}/resources/P7_ECox_2015/N4/P7_ECox_mask.mnc

tmpdir=$(mktemp -d)

input=$1
output=$2

# Added by Lani Cupo
# Reorient scans from CIC to MINC RAS
volmash -swap zy $input ${tmpdir}/mash.mnc
volflip -x ${tmpdir}/mash.mnc ${tmpdir}/flip.mnc

volcentre ${tmpdir}/flip.mnc ${tmpdir}/centered.mnc

mincmath -clamp -const2 0 $(mincstats -max -quiet ${tmpdir}/centered.mnc) ${tmpdir}/centered.mnc ${tmpdir}/centered.clamp.mnc
mv -f ${tmpdir}/centered.clamp.mnc  ${tmpdir}/centered.mnc

ThresholdImage 3 ${tmpdir}/centered.mnc ${tmpdir}/thresholdmask.mnc Otsu 1
minccalc -unsigned -byte -expression '1' ${tmpdir}/centered.mnc ${tmpdir}/fullmask.mnc

N4BiasFieldCorrection -d 3 -s 8 -i ${tmpdir}/centered.mnc -b [20] -c [300x300x300,1e-6] \
  -w ${tmpdir}/thresholdmask.mnc -x ${tmpdir}/fullmask.mnc -o ${tmpdir}/N4.mnc --verbose

volcentre -clobber ${tmpdir}/N4.mnc ${tmpdir}/centered.mnc

ThresholdImage 3 ${tmpdir}/centered.mnc ${tmpdir}/thresholdmask.mnc Otsu 1

autocrop -bbox ${tmpdir}/thresholdmask.mnc -isoexpand 20v ${tmpdir}/centered.mnc ${output}

#ExtractRegionFromImageByMask 3 ${tmpdir}/centered.mnc ${output} ${tmpdir}/thresholdmask.mnc 1 10

rm -rf ${tmpdir}

exit

antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 --minc --verbose \
  --output ${tmpdir}/trans \
  --use-histogram-matching 0 \
  --initial-moving-transform [${REGTARGET},${tmpdir}/N4.mnc,1] \
  --transform Rigid[0.1] \
    --metric Mattes[${REGTARGET},${tmpdir}/N4.mnc,1,64,None] \
    --convergence [2000x2000x2000,1e-6,10,1] --shrink-factors 8x6x4 --smoothing-sigmas 4x3x2 --masks [NULL,NULL] \
  --transform Similarity[0.1] \
    --metric Mattes[${REGTARGET},${tmpdir}/N4.mnc,1,64,None] \
    --convergence [2000x2000x2000,1e-6,10,1] --shrink-factors 6x4x2 --smoothing-sigmas 3x2x1 --masks [NULL,NULL] \
  --transform Affine[0.1] \
    --metric Mattes[${REGTARGET},${tmpdir}/N4.mnc,1,64,None] \
    --convergence [2000x2000x2000,1e-6,10,1] --shrink-factors 6x4x2 --smoothing-sigmas 3x2x1 --masks [NULL,NULL] \
  --transform Affine[0.1] \
    --metric Mattes[${REGTARGET},${tmpdir}/N4.mnc,1,64,None] \
    --convergence [2000x2000x50,1e-6,10,1] --shrink-factors 4x2x1 --smoothing-sigmas 2x1x0 --masks [${REGMASK},NULL]

antsApplyTransforms -d 3 -i ${REGMASK} -o ${tmpdir}/mask.mnc -t [${tmpdir}/trans0_GenericAffine.xfm,1] -r ${tmpdir}/centered.mnc -n GenericLabel --verbose

N4BiasFieldCorrection -d 3 -s 8 -i ${tmpdir}/centered.mnc -b [20] -c [200x200x200x200,0.0] -w ${tmpdir}/mask.mnc -r 1 -x ${tmpdir}/fullmask.mnc -o ${output} --verbose

cp ${tmpdir}/mask.mnc $(dirname ${output})/$(basename ${output} .mnc)_mask.mnc

rm -rf ${tmpdir}
