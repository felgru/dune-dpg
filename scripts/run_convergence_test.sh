#!/bin/bash

N="10 20 40 60 80 100 120 140 160 200 250 300"
for n in $N
do
  for CONVERGENCE_TEST in src/convergence_test_ls*
  do
    echo ${CONVERGENCE_TEST} $n
    eval ${CONVERGENCE_TEST} $n
  done
done
