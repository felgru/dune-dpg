#!/bin/bash

N="$@"
for n in $N
do
  echo -n "$n "
  src/profile_testspacecoefficientmatrix_ls2_ks3 $n | tr '\n' ' ' | \
    sed 's/^[^0-9]*\([0-9]*\)us[^0-9]*\([0-9]*\)us[^0-9]*\([0-9]*\)us. /\1 \2 \3\n/'
done
