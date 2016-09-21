#!/bin/bash

N="$@"
for n in $N
do
  for CONVERGENCE_TEST in src/convergence_test_ls*
  do
    echo ${CONVERGENCE_TEST} $n
    eval ${CONVERGENCE_TEST} $n
  done
done
