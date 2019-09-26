#!/bin/sh
declare -i nxval

for pow in {0..4}
do
  echo Now working $pow
  rm -f input.in
  nxval=2**${pow}
  nxval=16*${nxval}
  
  for iprob in {0..1}
  do

      #echo ${nxval} XB_1[$icase] $icase #${XB_$iprob}[$icase] #{XA_${iprob}} ${XA_${iprob}[$icase]} $iprob
      #echo "${pxa[$icase]} ${pxb[$icase]}"
      sed "s/NXVAR/${nxval}/g" input_template.in > input_interm1.in
      sed "s/IPROB/$iprob/g"   input_interm1.in  > input_interm2.in

      sed "s/XCEN/.true./g"    input_interm2.in  > input_${iprob}_T.in
      cp input_${iprob}_T.in input.in
      ./test_filters_NP
      
      sed "s/XCEN/.false./g"   input_interm2.in  > input_${iprob}_F.in
      cp input_${iprob}_F.in input.in
      ./test_filters_NP
  done
done

rm input_interm1.in
rm input_interm2.in

#TTXX
#TFXX
#FTXX
#FFXX
#
#TTXA
#TFXA
#FTXA
#FFXA
#
#TTXS
#TFXS
#FTXS
#FFXS
#
#TTAX
#TFAX
#FTAX
#FFAX
#
#TTAA
#TFAA
#FTAA
#FFAA
#
#TTAS
#TFAS
#FTAS
#FFAS
#
#TTSX
#TFSX
#FTSX
#FFSX
#
#TTSA
#TFSA
#FTSA
#FFSA
#
#TTSS
#TFSS
#FTSS
#FFSS
