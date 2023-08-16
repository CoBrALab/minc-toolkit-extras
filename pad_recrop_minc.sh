#!/bin/bash

calc() { awk "BEGIN{ print $* }"; }

 function pad_recrop_image() {
   local input=$1
   local padding=$2
   local output=$3

   local localtmpdir=$(mktemp -d)
   local res=$(calc "(sqrt(($(mincinfo -attvalue xspace:step ${input}))**2) + sqrt(($(mincinfo -attvalue yspace:step ${input}))**2) + sqrt(($(mincinfo -attvalue zspace:step ${input}))**2))/3")
   volpad -noauto -distance $(calc "int(50/${res})") ${input} ${localtmpdir}/padded.mnc

   ThresholdImage 3 ${localtmpdir}/padded.mnc ${localtmpdir}/bgmask.mnc 1e-12 Inf 1 0
   ThresholdImage 3 ${localtmpdir}/padded.mnc ${localtmpdir}/otsu.mnc Otsu 4 ${localtmpdir}/bgmask.mnc
   ThresholdImage 3 ${localtmpdir}/otsu.mnc ${localtmpdir}/otsu.mnc 2 Inf 1 0

   ExtractRegionFromImageByMask 3 ${localtmpdir}/padded.mnc ${output} ${localtmpdir}/otsu.mnc 1 $(calc "int(${padding}/${res})")

   rm -rf ${localtmpdir}
 }

function _main() {
  pad_recrop_image $1 $2 $3
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    _main "$@"
fi
