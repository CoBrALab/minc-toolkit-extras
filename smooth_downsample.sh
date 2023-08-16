#!/bin/bash
#Isotropize, and normalize intensity range, this is the file that will be processed in the pipeline
#Need smoothing for downsampling to avoid aliasing
#Ideas stolen from https://discourse.itk.org/t/resampling-to-isotropic-signal-processing-theory/1403

set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "You must enter exactly 3 command line arguments: isostep input output"
fi

tmpdir=$(mktemp -d)

input=$1
isostep=$2
output=$3

inputres=$(python -c "print('\n'.join([str(abs(x)) for x in [float(x) for x in \"$(PrintHeader ${input} 1)\".split(\"x\")]]))")

blurs=""
for dim in ${inputres}; do
  if [[ $(python -c "print(${dim}>(${isostep}-1e-6))") == True ]]; then
    blurs+=1e-12x
  else
    blurs+=$(python -c "import math; print(math.sqrt((${isostep}**2.0 - ${dim}**2.0)/(2.0*math.sqrt(2.0*math.log(2.0)))**2.0))")x
  fi
done

SmoothImage 3 ${input} "${blurs%?}" ${tmpdir}/smoothed.h5 1 0
ResampleImage 3 ${tmpdir}/smoothed.h5 ${tmpdir}/downsampled.h5 ${isostep}x${isostep}x${isostep} 0 4

ThresholdImage 3 ${tmpdir}/downsampled.h5 ${tmpdir}/mask.h5 0 Inf 1 0
ImageMath 3 ${output} m ${tmpdir}/downsampled.h5 ${tmpdir}/mask.h5

rm -rf ${tmpdir}
