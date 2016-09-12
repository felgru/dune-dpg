#!/bin/bash

N="10 20 40 60 80 100 120 140 160 200 250 300"
for n in $N
do
  echo -n "$n "
  src/profile_testspacecoefficientmatrix_ls2_ks3 $n | tr '\n' ' ' | \
    sed 's/^[^0-9]*\([0-9]*\)us[^0-9]*\([0-9]*\)us[^0-9]*\([0-9]*\)us. /\1 \2 \3\n/'
done
