#!/bin/bash
#Script to do optimal preprocessing on in-vivo/ex-vivo structural scans
#Taken using the CIC Bruker 7T
#usage:
#mouse-preprocessing-v2.sh input.mnc output.mnc

#Operations
# centers brain in space
# denoises
# registers to Ratlas
# gets approximate brain mask from atlas
# expands mask
# Bias field correction with N4

set -euo pipefail

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${THREADS_PER_COMMAND:-$(nproc)}

REGTARGET=${QUARANTINE_PATH}/resources/Fischer344/Fischer344_template.mnc
REGMASK=${QUARANTINE_PATH}/resources/Fischer344/Fischer344_mask.mnc

tmpdir=$(mktemp -d)

input=$1
output=$2

set -x

maxval=$(mincstats -max -quiet ${input})

volflip -y $input $tmpdir/flip.mnc
volmash -swap zy $tmpdir/flip.mnc $tmpdir/mash.mnc

minc_anlm --mt ${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS} $tmpdir/mash.mnc $tmpdir/denoise.mnc

ImageMath 3 ${tmpdir}/weight.mnc ThresholdAtMean $tmpdir/denoise.mnc 1

ImageMath 3 ${tmpdir}/weight.mnc GetLargestComponent ${tmpdir}/weight.mnc

minccalc -byte -unsigned -expression '1' $tmpdir/denoise.mnc $tmpdir/fullmask.mnc

N4BiasFieldCorrection -d 3 -i $tmpdir/denoise.mnc -b [30] -c [200x200x200x200,0.0] -r 0 -x $tmpdir/fullmask.mnc -w ${tmpdir}/weight.mnc -o $tmpdir/corrected.mnc -s 4 --verbose

ThresholdImage 3 $tmpdir/corrected.mnc $tmpdir/otsu.mnc  Otsu 1

ImageMath 3 $tmpdir/otsu.mnc GetLargestComponent $tmpdir/otsu.mnc

N4BiasFieldCorrection -d 3 -i $tmpdir/denoise.mnc -b [30] -c [200x200x200x200,0.0] -r 0 -w $tmpdir/otsu.mnc -x $tmpdir/fullmask.mnc -o $tmpdir/corrected.mnc -s 2 --verbose

antsRegistration --dimensionality 3 --float 1 --collapse-output-transforms 1 --minc \
--output $tmpdir/trans \
--use-histogram-matching 0 \
--initial-moving-transform [${REGTARGET},$tmpdir/corrected.mnc,1] \
--transform Rigid[0.1] --metric Mattes[${REGTARGET},$tmpdir/corrected.mnc,1,32,None] --convergence [2000x2000x2000x2000x2000x0,1e-6,10,1] --shrink-factors 10x8x6x4x2x1 --smoothing-sigmas 5x4x3x2x1x0 --masks [NULL,NULL] \
--transform Similarity[0.1] --metric Mattes[${REGTARGET},$tmpdir/corrected.mnc,1,32,None] --convergence [2000x2000x2000,1e-6,10,1] --shrink-factors 6x4x2 --smoothing-sigmas 3x2x1 --masks [NULL,NULL] \
--transform Affine[0.1]     --metric Mattes[${REGTARGET},$tmpdir/corrected.mnc,1,32,None] --convergence [2000x2000x0,1e-6,10,1] --shrink-factors 4x2x1 --smoothing-sigmas 2x1x0 --masks [${REGMASK},NULL] --verbose

antsApplyTransforms -d 3 -i ${REGMASK} -o $tmpdir/mask.mnc -t [$tmpdir/trans0_GenericAffine.xfm,1] -r $tmpdir/denoise.mnc -n NearestNeighbor

ThresholdImage 3 $tmpdir/denoise.mnc ${tmpdir}/weight.mnc Otsu 1 $tmpdir/mask.mnc

ThresholdImage 3 ${tmpdir}/weight.mnc ${tmpdir}/weight.mnc 2 2 1 0

ImageMath 3 ${tmpdir}/weight.mnc FillHoles ${tmpdir}/weight.mnc 2

N4BiasFieldCorrection -d 3 -i $tmpdir/denoise.mnc -b [30] -c [200x200x200,0.0] -w $tmpdir/weight.mnc -r 0 -x $tmpdir/fullmask.mnc -o $tmpdir/corrected.mnc -s 2 --verbose

ImageMath 3 $tmpdir/corrected.mnc TruncateImageIntensity $tmpdir/corrected.mnc 0.0005 0.9995 1024 $tmpdir/mask.mnc

ImageMath 3 ${tmpdir}/corrected.mnc RescaleImage $tmpdir/corrected.mnc 0 ${maxval}

autocrop -bbox $tmpdir/mask.mnc -isoexpand 10v ${tmpdir}/corrected.mnc ${tmpdir}/cropped.mnc
mincresample -like ${tmpdir}/cropped.mnc -keep -near -byte -unsigned $tmpdir/mask.mnc ${tmpdir}/mask2.mnc

affine_xfm_to_rigid.sh $tmpdir/trans0_GenericAffine.xfm ${tmpdir}/rigid.xfm

mincresample -tfm_input_sampling -transform ${tmpdir}/rigid.xfm ${tmpdir}/cropped.mnc $(dirname $output)/$(basename $output .mnc)_lsq6.mnc
antsApplyTransforms -d 3 -i ${tmpdir}/mask2.mnc -r $(dirname $output)/$(basename $output .mnc)_lsq6.mnc -t [${tmpdir}/rigid.xfm,1] -o $(dirname $output)/$(basename $output .mnc)_lsq6_mask.mnc -n GenericLabel

volcentre ${tmpdir}/cropped.mnc ${output}
volcentre ${tmpdir}/mask2.mnc $(dirname $output)/$(basename $output .mnc)_mask.mnc

rm -rf $tmpdir
