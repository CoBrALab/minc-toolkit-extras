#!/bin/bash
#USAGE
#create_ciet_image.sh left_average.obj left_statmap left_column left_FDR right_average.obj right_statmap right_column right_FDR output.png

#Automatically generates a standardized view of a brain, using stats files for left and right hemispheres (from RMINC)
#Colourizes objects using stats files, ray traces standard views and them merges them together

set -euo pipefail
IFS=$'\n\t'

if [[ $# -ne 9 ]]
then
echo Usage: create_ciet_image.sh left_average.obj left_statmap left_column left_FDR right_average.obj right_statmap right_column right_FDR output.png
exit
fi


TMPDIR=$(mktemp -d)
leftbrain=$1
leftstatmap=$2
leftcolumn=$3
leftFDR=$4
rightbrain=$5
rightstatmap=$6
rightcolumn=$7
rightFDR=$8
output=$9

#Pick out column, multiply by mask and store in new tempfile
cut -d " " -f $leftcolumn $leftstatmap | paste -d\* /opt/quarantine/resources/CIVET_2.0_mask_left_short.txt - | sed 's/[eE]+\{0,1\}/*10^/g' | bc > $TMPDIR/$(basename $leftstatmap)
cut -d " " -f $rightcolumn $rightstatmap | paste -d\* /opt/quarantine/resources/CIVET_2.0_mask_left_short.txt - | sed 's/[eE]+\{0,1\}/*10^/g' | bc > $TMPDIR/$(basename $rightstatmap)

#Redefine statmap to masked tempfile
leftstatmap=$TMPDIR/$(basename $leftstatmap)
rightstatmap=$TMPDIR/$(basename $rightstatmap)

#Colourize left
if (( $(echo "$leftFDR <= $(sort -g $leftstatmap | tail -1)" | bc) ))
then
    colour_object $leftbrain $leftstatmap $TMPDIR/left_pos.obj red_metal_inv $leftFDR $(sort -g $leftstatmap | tail -1) transparent black composite
else
    cp $leftbrain $TMPDIR/left_pos.obj
fi

if (( $(echo "-$leftFDR >= $(sort -g $leftstatmap  | head -1)" | bc) ))
then
    colour_object $TMPDIR/left_pos.obj $leftstatmap $TMPDIR/left_comb.obj cold_metal $(sort -g $leftstatmap  | head -1) -$leftFDR black transparent composite
else
    cp  $TMPDIR/left_pos.obj $TMPDIR/left_comb.obj
fi

#Colourize right
if (( $(echo "$rightFDR <= $(sort -g $rightstatmap | tail -1)" | bc) ))
then
    colour_object $rightbrain $rightstatmap $TMPDIR/right_pos.obj red_metal_inv $rightFDR $(sort -g $rightstatmap | tail -1) transparent black composite
else
    cp $rightbrain $TMPDIR/right_pos.obj
fi

if (( $(echo "-$rightFDR >= $(sort -g $rightstatmap  | head -1)" | bc) ))
then
    colour_object $TMPDIR/right_pos.obj $rightstatmap $TMPDIR/right_comb.obj cold_metal $(sort -g $rightstatmap  | head -1) -$rightFDR black transparent composite
else
    cp $TMPDIR/right_pos.obj $TMPDIR/right_comb.obj
fi

#Generate views
echo """ray_trace -output $TMPDIR/bottom.rgb $TMPDIR/left_comb.obj $TMPDIR/right_comb.obj -size 1000 1000 -crop -sup 3 -bg black -bottom
ray_trace -output $TMPDIR/top.rgb $TMPDIR/left_comb.obj $TMPDIR/right_comb.obj -size 1000 1000 -crop -sup 3 -bg black -top
ray_trace -output $TMPDIR/front.rgb $TMPDIR/left_comb.obj $TMPDIR/right_comb.obj -size 1000 1000 -crop -sup 3 -bg black -front
ray_trace -output $TMPDIR/back.rgb $TMPDIR/left_comb.obj $TMPDIR/right_comb.obj -size 1000 1000 -crop -sup 3 -bg black -back
ray_trace -output $TMPDIR/left_medial.rgb $TMPDIR/left_comb.obj -size 1000 1000 -crop -sup 3 -bg black -right
ray_trace -output $TMPDIR/left_lateral.rgb $TMPDIR/left_comb.obj -size 1000 1000 -crop -sup 3 -bg black -left
ray_trace -output $TMPDIR/right_medial.rgb $TMPDIR/right_comb.obj -size 1000 1000 -crop -sup 3 -bg black -left
ray_trace -output $TMPDIR/right_lateral.rgb $TMPDIR/right_comb.obj -size 1000 1000 -crop -sup 3 -bg black  -right""" | parallel -j8

#Create pngs
for file in $TMPDIR/*.rgb
do
    echo convert $file $TMPDIR/$(basename $file rgb)png
done | parallel -j8

#Create colourbars
convert xc:black xc:rgb\(0,0,128\)  xc:rgb\(0,128,255\) xc:rgb\(128,255,255\) xc:white  -append -filter Cubic -resize 30x600\! -bordercolor black -border 2x2 -bordercolor white -border 2x2 $TMPDIR/coldmetalinv.png
convert xc:black xc:rgb\(128,0,0\)  xc:rgb\(255,0,128\) xc:rgb\(255,128,255\) xc:white  -append -filter Cubic -resize 30x600\! -bordercolor black -border 2x2 -bordercolor white -border 2x2 $TMPDIR/redmetalinv.png

#Create combined figure
convert -gravity Center -background black +append $TMPDIR/left_medial.png $TMPDIR/right_medial.png $TMPDIR/1.png
convert -gravity Center -background black +append $TMPDIR/left_lateral.png $TMPDIR/right_lateral.png $TMPDIR/2.png
convert -gravity Center -background black +append $TMPDIR/top.png $TMPDIR/bottom.png $TMPDIR/3.png
convert -gravity Center -background black +append $TMPDIR/front.png $TMPDIR/back.png $TMPDIR/4.png
convert $TMPDIR/1.png $TMPDIR/2.png $TMPDIR/3.png $TMPDIR/4.png -gravity Center -background black -append $TMPDIR/unlabelled.png

convert $TMPDIR/unlabelled.png \
-page +1364+1685 $TMPDIR/coldmetalinv.png \
-page +199+1685 $TMPDIR/redmetalinv.png \
-stroke  none   -fill white  -pointsize 45 -annotate +1410+2285 'FDR 5%' \
-stroke  none   -fill white  -pointsize 45 -annotate +1410+1725 'FDR<1%' \
-stroke  none   -fill white  -pointsize 45 -annotate +5+2285 'FDR 5%' \
-stroke  none   -fill white  -pointsize 45 -annotate +5+1725 'FDR<1%' \
-flatten $output

rm -r $TMPDIR
