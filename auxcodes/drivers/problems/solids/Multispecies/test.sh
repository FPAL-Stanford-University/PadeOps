#!/bin/sh
declare -i nxval

for pow in {1..4}
do
  echo Now working $pow
  rm -f input.in
  nxval=50*${pow}
  nxval=${nxval}+1
  echo $nxval
  sed "s/NXVAR/${nxval}/g" input_template.in > input.in
  ./Multispecies input.in
  mv exact_solution.dat exact_solution_${pow}.dat
  mv input.in input_${pow}.in
done
