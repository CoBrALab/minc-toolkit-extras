#!/bin/bash
#Generate auto-scaled MaGeT QC image
#TODO
#-option for label on-off slices
#-option for number of slices
#-option for image size
#-clobber check

set -e

image=$1
label=$2
output=$3
tmpdir=$(mktemp -d)

mincresample $(mincbbox -mincresample $label) $label $tmpdir/label-crop.mnc
minccalc -expression "A[0]?1:1" $tmpdir/label-crop.mnc $tmpdir/bounding.mnc

#Trasverse
create_verify_image $tmpdir/$(basename $image .mnc)_t.rgb \
-width 3840 -autocols 20 -autocol_planes t \
-bounding_volume $tmpdir/bounding.mnc \
-row $image color:gray \
volume_overlay:$label:0.4

create_verify_image $tmpdir/$(basename $image .mnc)_t2.rgb \
-width 3840 -autocols 20 -autocol_planes t \
-bounding_volume $tmpdir/bounding.mnc \
-row $image color:gray

#Saggital
create_verify_image $tmpdir/$(basename $image .mnc)_s.rgb \
-width 3840 -autocols 20 -autocol_planes s \
-bounding_volume $tmpdir/bounding.mnc \
-row $image color:gray \
volume_overlay:$label:0.4

create_verify_image $tmpdir/$(basename $image .mnc)_s2.rgb \
-width 3840 -autocols 20 -autocol_planes s \
-bounding_volume $tmpdir/bounding.mnc \
-row $image color:gray

#Coronal
create_verify_image $tmpdir/$(basename $image .mnc)_c.rgb \
-width 3840 -autocols 20 -autocol_planes c \
-bounding_volume $tmpdir/bounding.mnc \
-row $image color:gray \
volume_overlay:$label:0.4

create_verify_image $tmpdir/$(basename $image .mnc)_c2.rgb \
-width 3840 -autocols 20 -autocol_planes c \
-bounding_volume $tmpdir/bounding.mnc \
-row $image color:gray

convert -append -trim $tmpdir/*.rgb $output

rm -rf $tmpdir
