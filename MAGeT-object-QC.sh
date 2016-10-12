#!/bin/bash
#Generate auto-scaled MaGeT QC image
#TODO
#-option for label on-off slices
#-option for number of slices
#-option for image size
#-clobber check

set -e

image=$1
shift
output="${@: -1}"
tmpdir=$(mktemp -d)


#Join objects with objconcat
objconcat ${@:1:$#-1} $(for((i=1;i<=$# - 1;i++));do printf "%s" "none ";done) $tmpdir/combined.obj none
scan_object_to_volume  $image  $tmpdir/combined.obj  $tmpdir/label.mnc
mincresample $(mincbbox -mincresample $tmpdir/label.mnc) $tmpdir/label.mnc $tmpdir/label-crop.mnc
minccalc -expression "A[0]?1:1" $tmpdir/label-crop.mnc $tmpdir/bounding.mnc

colours=(red green blue magenta orange)

overlay_list=""
i=0
for object in ${@:1:$#-1}
do
  overlay_list+="surface_overlay:$object:${colours[$i]} "
((i++)) || true
done

#Trasverse
echo create_verify_image $tmpdir/$(basename $image .mnc)_t.rgb \
-width 3840 -autocols 20 -autocol_planes t \
-bounding_volume $tmpdir/bounding.mnc \
-row $image color:gray \
$overlay_list

create_verify_image $tmpdir/$(basename $image .mnc)_t2.rgb \
-width 3840 -autocols 20 -autocol_planes t \
-bounding_volume $tmpdir/bounding.mnc \
-row $image color:gray

#Saggital
create_verify_image $tmpdir/$(basename $image .mnc)_s.rgb \
-width 3840 -autocols 20 -autocol_planes s \
-bounding_volume $tmpdir/bounding.mnc \
-row $image color:gray \
$overlay_list

create_verify_image $tmpdir/$(basename $image .mnc)_s2.rgb \
-width 3840 -autocols 20 -autocol_planes s \
-bounding_volume $tmpdir/bounding.mnc \
-row $image color:gray

#Coronal
create_verify_image $tmpdir/$(basename $image .mnc)_c.rgb \
-width 3840 -autocols 20 -autocol_planes c \
-bounding_volume $tmpdir/bounding.mnc \
-row $image color:gray \
$overlay_list

create_verify_image $tmpdir/$(basename $image .mnc)_c2.rgb \
-width 3840 -autocols 20 -autocol_planes c \
-bounding_volume $tmpdir/bounding.mnc \
-row $image color:gray

convert -append -trim $tmpdir/*.rgb $output

rm -rf $tmpdir
