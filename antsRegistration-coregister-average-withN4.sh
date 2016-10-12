#!/bin/bash
set -e

rounds=4
output=$1
tmpdir=$(mktemp -d)
shift 1

for file in "$@"
do

iterativeN4.sh $file $tmpdir/$(basename $file).N4.mnc

done

set +e

mincbigaverage --robust --sdfile $tmpdir/sd0.mnc $tmpdir/*.N4.mnc $tmpdir/average0.mnc

target=$tmpdir/average0.mnc

iters="2000"
shrink=$(expr $rounds - 1)
sigma=$(expr $rounds - 2)

for ((n=1;n<$rounds;n++))
do
  for file in $tmpdir/*.N4.mnc
  do
    antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 --minc --verbose --interpolation Linear \
                --output [$tmpdir/$(basename $file)-$(basename $target),$tmpdir/$(basename $file)-$(basename $target)] \
                --transform Rigid[0.1] --metric GC[$target,$file,1] \
                --convergence [$iters,1e-7,10] --shrink-factors $shrink --smoothing-sigmas $sigma
  done
  mincaverage -normalize -sdfile $tmpdir/sd$n.mnc $tmpdir/*-$(basename $target) $tmpdir/average$n.mnc

  target=$tmpdir/average$n.mnc
  iters+="x2000"
  shrink+="x$(expr $rounds - $n - 1)"
  sigma+="x$(expr $rounds - $n - 2)"
done

#cp $target $output
mincmath -clamp -const2 0 $(mincstats -quiet -max $tmpdir/average0.mnc) $target $output
rm -r $tmpdir
