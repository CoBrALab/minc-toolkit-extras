#!/bin/bash

set -euo pipefail

input=$1
output=$2

tmpdir=$(mktemp -d)

param2xfm $(xfm2param ${input} | grep -E 'scale|shear') ${tmpdir}/scaleshear.xfm
xfminvert ${tmpdir}/scaleshear.xfm ${tmpdir}/unscaleshear.xfm
xfmconcat ${input} ${tmpdir}/unscaleshear.xfm ${output}

rm -rf ${tmpdir}
