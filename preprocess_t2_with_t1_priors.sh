#!/bin/bash

set -euo pipefail

tmpdir=$(mktemp -d)

if [[ $# -lt 5 ]] ; then
    echo 'Usage: preprocess.sh T1.mnc T1_brainvolume.mnc T1_brainmask.mnc T2_N4_denoise_extracted.mnc'
    exit 1
fi

t1=$1
t1brainvolume=$2
t1brainmask=$3
t2=$4
output=$5

antsRegistration_affine_SyN.sh --close --skip-nonlinear --linear-type rigid $t2 $t1 ${tmpdir}/t2_to_t1

antsApplyTransforms -d 3 -i ${t1brainvolume} \
-t [ ${tmpdir}/t2_to_t10_GenericAffine.xfm,1 ] \
-o ${tmpdir}/weight.mnc -n GenericLabel -r ${t2} --verbose

antsApplyTransforms -d 3 -i ${t1brainmask} \
-t [ ${tmpdir}/t2_to_t10_GenericAffine.xfm,1 ] \
-o ${tmpdir}/mask.mnc -n GenericLabel -r ${t2} --verbose


N4BiasFieldCorrection -d 3 -s 4 -b [ 200 ] -c [ 50x50x50,0.0 ] -w  ${tmpdir}/weight.mnc -i ${t2} -o ${tmpdir}/corrected.mnc --verbose

minc_anlm --mt $(nproc) ${tmpdir}/corrected.mnc ${tmpdir}/corrected.denoise.mnc

ImageMath 3 ${output} m ${tmpdir}/corrected.denoise.mnc ${tmpdir}/mask.mnc

rm -rf ${tmpdir}
