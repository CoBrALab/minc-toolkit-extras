#!/bin/bash
#Preprocessing script which does intrasubject registration and uses brain mask to bias field correct and extract scan
set -euo pipefail
set -x

tmpdir=$(mktemp -d)

slab=$1
wholebrain=$2
wholebrain_mask=$3
output=$4

bestlinreg_g -nmi -lsq12 ${slab} ${wholebrain} -target_mask ${wholebrain_mask} ${tmpdir}/towholebrain.xfm
mincresample -like ${slab} -keep -near -byte -unsigned -invert_transform -transform ${tmpdir}/towholebrain.xfm ${wholebrain_mask} ${tmpdir}/slabmask.mnc

minccalc -byte -unsigned -expression '1' ${slab} ${tmpdir}/initmask.mnc

N4BiasFieldCorrection -d 3 -i ${slab} -w ${tmpdir}/slabmask.mnc -x ${tmpdir}/initmask.mnc -r 1 -b [200] -c [200x200x200x200,0.0] -o ${tmpdir}/N4.mnc --verbose

mincresample -keep -near -like ${tmpdir}/N4.mnc ${tmpdir}/slabmask.mnc ${tmpdir}/slabmask2.mnc
minccalc -clobber -expression 'A[0]*A[1]' ${tmpdir}/N4.mnc ${tmpdir}/slabmask2.mnc ${output}

rm -rf ${tmpdir}
