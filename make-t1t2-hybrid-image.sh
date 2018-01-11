#!/bin/bash

set -euo pipefail

t1=$1
t2=$2
classify=$3
mask=$4
output=$5

t1median=$(mincstats -quiet -median -mask $classify -mask_binvalue 2 $t1)
t2median=$(mincstats -quiet -median -mask $classify -mask_binvalue 2 $t2)

minccalc -expression "A[2]>0?(A[0] - ${t1median}/${t2median}*A[1])/(A[0] + ${t1median}/${t2median}*A[1]):0" $t1 $t2 $mask $output
