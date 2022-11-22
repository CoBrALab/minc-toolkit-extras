#!/bin/bash
# Script to collect volumes in parallel, and use names from a csv if available, see https://github.com/CobraLab/atlases for csv examples
# CSV Format
# 1,name
# 2,name
# etc

set -euo pipefail
IFS=$'\n\t'

tmpdir=$(mktemp -d)

firstarg=${1:-}
if ! [[  ( $firstarg = *csv ) ||  ( $firstarg = *mnc )  ]]; then
    echo "usage: $0 labels.mnc jacobian1.nii(.gz) [ jacobian2.nii(.gz) ... jacobianN.nii(.gz) ]"
    exit 1
fi


AWKCMD='BEGIN{FS=","}

{for(i=1;i<=NF;i++)a[NR,i]=$i}
END{
  for(i=1;i<=NF;i++){
    for(j=1;j<=NR;j++){
      line=line sep a[j,i]
      sep=FS
    }
    print line
    line=sep=""
  }
}
'

labels=$1
shift

{ echo -n "Subject,"; nii2mnc -quiet ${1} ${tmpdir}/temp.mnc && label_volumes_from_jacobians ${labels} ${tmpdir}/temp.mnc; } | cut -f1 -d"," | awk -f<(echo "$AWKCMD")
for file in "$@"
do
    sem -j+0 "nii2mnc -quiet ${file} ${tmpdir}/$(basename ${file} .nii.gz).mnc && label_volumes_from_jacobians ${labels} ${tmpdir}/$(basename ${file} .nii.gz).mnc | tail -n +2 | cut -f2 -d"," | { echo -n "$file,"; awk -f<(echo '$AWKCMD') ; } | tr -d '[:space:]' && rm ${tmpdir}/$(basename ${file} .nii.gz).mnc"
done
sem --wait

rm -rf ${tmpdir}
