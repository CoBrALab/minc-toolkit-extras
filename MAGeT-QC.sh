#!/bin/bash
#Generate auto-scaled MaGeT QC image
#TODO
#-option for label on-off slices
#-option for number of slices
#-option for image size
#-clobber check

set -euo pipefail

image=$1
label=$2
output=$3
tmpdir=$(mktemp -d)

if [[ ${label} == "none" ]]; then
    ThresholdImage 3 ${image} ${tmpdir}/bgmask.mnc 1e-12 Inf 1 0
    ThresholdImage 3 ${image} ${tmpdir}/bgmask.mnc Otsu 4 ${tmpdir}/bgmask.mnc
    ThresholdImage 3 ${tmpdir}/bgmask.mnc ${tmpdir}/bgmask.mnc 2 Inf 1 0
    mincresample -keep -near -like ${image} ${tmpdir}/bgmask.mnc ${tmpdir}/resample.mnc
    mincresample $(mincbbox -mincresample ${tmpdir}/resample.mnc) ${tmpdir}/resample.mnc ${tmpdir}/label-crop.mnc
    minccalc -expression "1" ${tmpdir}/label-crop.mnc ${tmpdir}/bounding.mnc
else
    mincresample $(mincbbox -mincresample ${label}) ${label} ${tmpdir}/label-crop.mnc
    minccalc -expression "1" ${tmpdir}/label-crop.mnc ${tmpdir}/bounding.mnc
    #Labelled layers
    create_verify_image -align_com -range_floor 0 ${tmpdir}/$(basename ${image} .mnc)_t.rgb \
        -width 1920 -autocols 10 -autocol_planes t \
        -bounding_volume ${tmpdir}/bounding.mnc \
        -row ${image} color:gray \
        volume_overlay:${label}:0.3

    create_verify_image -align_com -range_floor 0 ${tmpdir}/$(basename ${image} .mnc)_s.rgb \
        -width 1920 -autocols 10 -autocol_planes s \
        -bounding_volume ${tmpdir}/bounding.mnc \
        -row ${image} color:gray \
        volume_overlay:${label}:0.3

    create_verify_image -align_com -range_floor 0 ${tmpdir}/$(basename ${image} .mnc)_c.rgb \
        -width 1920 -autocols 10 -autocol_planes c \
        -bounding_volume ${tmpdir}/bounding.mnc \
        -row ${image} color:gray \
        volume_overlay:${label}:0.3
fi





#Trasverse
create_verify_image -align_com -range_floor 0 ${tmpdir}/$(basename ${image} .mnc)_t2.rgb \
    -width 1920 -autocols 10 -autocol_planes t \
    -bounding_volume ${tmpdir}/bounding.mnc \
    -row ${image} color:gray

#Saggital
create_verify_image -align_com -range_floor 0 ${tmpdir}/$(basename ${image} .mnc)_s2.rgb \
    -width 1920 -autocols 10 -autocol_planes s \
    -bounding_volume ${tmpdir}/bounding.mnc \
    -row ${image} color:gray

#Coronal
create_verify_image -align_com -range_floor 0 ${tmpdir}/$(basename ${image} .mnc)_c2.rgb \
    -width 1920 -autocols 10 -autocol_planes c \
    -bounding_volume ${tmpdir}/bounding.mnc \
    -row ${image} color:gray

convert -background black -strip -interlace Plane -sampling-factor 4:2:0 -quality "85%"  -append -trim ${tmpdir}/*.rgb $output

rm -rf ${tmpdir}
