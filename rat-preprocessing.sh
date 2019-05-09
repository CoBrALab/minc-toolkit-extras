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

REGTARGET=/opt/quarantine/resources/Ratlas/cropped/rescale_centered.mnc
REGMASK=/opt/quarantine/resources/Ratlas/cropped/rat_mask_centered.mnc

tmpdir=$(mktemp -d)

input=$1
output=$2

set -x

maxval=$(mincstats -max -quiet ${input})

volflip -y $input $tmpdir/flip.mnc
volmash -swap zy $tmpdir/flip.mnc $tmpdir/mash.mnc


clean_and_center_minc.pl $tmpdir/mash.mnc $tmpdir/centered.mnc

minccalc -byte -unsigned -expression 'A[0]?1:1' $tmpdir/centered.mnc $tmpdir/fullmask.mnc

N4BiasFieldCorrection -d 3 -i $tmpdir/centered.mnc -b [30] -c [200x200x200,0.0] -r 0 -x $tmpdir/fullmask.mnc -o $tmpdir/corrected.mnc -s 4 --verbose

ThresholdImage 3 $tmpdir/corrected.mnc $tmpdir/otsu.mnc  Otsu 1

N4BiasFieldCorrection -d 3 -i $tmpdir/centered.mnc -b [30] -c [200x200x200,0.0] -r 0 -w $tmpdir/otsu.mnc -x $tmpdir/fullmask.mnc -o $tmpdir/corrected.mnc -s 2 --verbose

ImageMath 3 $output RescaleImage $tmpdir/corrected.mnc 0 ${maxval}


rm -rf $tmpdir
