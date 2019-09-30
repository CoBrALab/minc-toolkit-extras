#!/bin/bash
#Script to do optimal preprocessing on in-vivo/ex-vivo structural scans
#Taken using the CIC Bruker 7T
#usage:
#mouse-preprocessing-v2.sh input.mnc output.mnc

#Operations
# swaps zy to re-orient mouse
# flips x to fix left-right mixup
# centers brain in space
# denoises
# registers to DSUR_2016 atlas
# gets approximate brain mask from atlas
# expands mask
# Bias field correction with N4

set -euo pipefail
set -v
REGTARGET=${QUARANTINE_PATH}/resources/DSUR_2016/DSUR_40micron_average.mnc
REGMASK=${QUARANTINE_PATH}/resources/DSUR_2016/DSUR_40micron_mask_version2.mnc

tmpdir=$(mktemp -d)

cp $1 $tmpdir/input.mnc

input=$tmpdir/input.mnc
output=$2

#Fix "too long history" from some dcm2mnc conversion
minc_modify_header $input -sinsert :history=‘’

#Swap orientations so that scans start oriented the same way as MiCE atlases
volmash -swap zy $input $tmpdir/mash.mnc
volflip $tmpdir/mash.mnc $tmpdir/flip.mnc

#Center in space and cleanup slicing
clean_and_center_minc.pl $tmpdir/flip.mnc $tmpdir/centered.mnc

#Denoise
minc_anlm $tmpdir/centered.mnc $tmpdir/denoise.mnc

#Make a full-scan mask to short-circuit some N4 sillyness
minccalc -byte -unsigned -expression 'A[0]?1:1' $tmpdir/denoise.mnc $tmpdir/fullmask.mnc

#Register to atlas
antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 --verbose --minc \
--output $tmpdir/trans \
--use-histogram-matching 0 \
--transform Rigid[0.1] --metric Mattes[$tmpdir/denoise.mnc,${REGTARGET},1] --convergence [2000x2000,1e-6,10] --shrink-factors 6x4 --smoothing-sigmas 3x2 --masks [NULL,NULL] \
--transform Similarity[0.1] --metric Mattes[$tmpdir/denoise.mnc,${REGTARGET},1] --convergence [2000x2000,1e-6,10] --shrink-factors 4x2 --smoothing-sigmas 2x1 --masks [NULL,NULL] \
--transform Affine[0.1]     --metric Mattes[$tmpdir/denoise.mnc,${REGTARGET},1] --convergence [2000x2000,1e-6,10] --shrink-factors 2x1 --smoothing-sigmas 2x1 --masks [NULL,NULL]  \
--transform Affine[0.1]     --metric Mattes[$tmpdir/denoise.mnc,${REGTARGET},1] --convergence [2000x2000,1e-6,10] --shrink-factors 2x1 --smoothing-sigmas 1x0.5 --masks [NULL,${REGMASK}]

#Transform brain mask back
antsApplyTransforms -d 3 -i ${REGMASK} -o $tmpdir/mask.mnc -t $tmpdir/trans0_GenericAffine.xfm -r $tmpdir/denoise.mnc -n NearestNeighbor

#Bias field first round
N4BiasFieldCorrection -d 3 -i $tmpdir/denoise.mnc -b [20] -c [200x200x200,0.0] -w $tmpdir/mask.mnc -r 1 -x $tmpdir/fullmask.mnc -o $tmpdir/denoise.N4.mnc

#Re-regiser to atlas
antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 --verbose --minc \
--output $tmpdir/trans \
--use-histogram-matching 0 \
--transform Rigid[0.1] --metric Mattes[$tmpdir/denoise.N4.mnc,${REGTARGET},1] --convergence [2000x2000,1e-6,10] --shrink-factors 6x4 --smoothing-sigmas 3x2 --masks [NULL,NULL] \
--transform Similarity[0.1] --metric Mattes[$tmpdir/denoise.N4.mnc,${REGTARGET},1] --convergence [2000x2000,1e-6,10] --shrink-factors 4x2 --smoothing-sigmas 2x1 --masks [NULL,NULL] \
--transform Affine[0.1]     --metric Mattes[$tmpdir/denoise.N4.mnc,${REGTARGET},1] --convergence [2000x2000,1e-6,10] --shrink-factors 2x1 --smoothing-sigmas 2x1 --masks [NULL,NULL]  \
--transform Affine[0.1]     --metric Mattes[$tmpdir/denoise.N4.mnc,${REGTARGET},1] --convergence [2000x2000,1e-6,10] --shrink-factors 2x1 --smoothing-sigmas 1x0.5 --masks [NULL,${REGMASK}]

#Transform new mask
antsApplyTransforms -d 3 -i ${REGMASK} -o $tmpdir/mask.mnc -t $tmpdir/trans0_GenericAffine.xfm -r $tmpdir/denoise.N4.mnc -n NearestNeighbor

#Correct again
N4BiasFieldCorrection -d 3 -i $tmpdir/denoise.mnc -b [20] -c [200x200x200,0.0] -w $tmpdir/mask.mnc -r 1 -x $tmpdir/fullmask.mnc -o $output

#Remove scale/shear parts to make an lsq6 transform
xfminvert $tmpdir/trans0_GenericAffine.xfm $tmpdir/trans0_GenericAffine_inverse.xfm
param2xfm $(xfm2param $tmpdir/trans0_GenericAffine_inverse.xfm | grep -E 'scale|shear') $tmpdir/scale.xfm
xfminvert $tmpdir/scale.xfm $tmpdir/unscale.xfm
xfmconcat $tmpdir/trans0_GenericAffine_inverse.xfm $tmpdir/unscale.xfm $tmpdir/lsq6.xfm

#Output masks
cp $tmpdir/mask.mnc $(dirname $output)/$(basename $output .mnc)_mask.mnc

#Provide lsq6 aligned scans in addition to native-space scans
mincresample -tfm_input_sampling -invert_transform -transform $tmpdir/lsq6.xfm $output $(dirname $output)/$(basename $output .mnc)_lsq6.mnc
mincresample -tfm_input_sampling -invert_transform -transform $tmpdir/lsq6.xfm $tmpdir/mask.mnc $(dirname $output)/$(basename $output .mnc)_lsq6_mask.mnc

rm -rf $tmpdir
