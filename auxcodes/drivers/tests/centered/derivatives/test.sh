#!/bin/sh
declare -i nxval

XA_1=(0 1   0  1)
XB_1=(0 0  -1 -1)

XA_2=(0 -1  0 -1)
XB_2=(0  0  1  1)

XA_3=(0  1  0  1)
XB_3=(0  0  1  1)

XA_4=(0 -1  0 -1)
XB_4=(0  0 -1 -1)


for pow in {0..4}
do
  echo Now working $pow
  rm -f input.in
  nxval=2**${pow}
  nxval=16*${nxval}
  
  for iprob in {1..4}
  do

    for icase in {0..3}
    do

      #echo ${nxval} XB_1[$icase] $icase #${XB_$iprob}[$icase] #{XA_${iprob}} ${XA_${iprob}[$icase]} $iprob
      #echo "${pxa[$icase]} ${pxb[$icase]}"
      sed "s/NXVAR/${nxval}/g"             input_template.in > input_interm1.in
      sed "s/IPROB/$iprob/g"             input_interm1.in  > input_interm2.in
      if [ $iprob -eq 1 ]
      then
        sed "s/XA/${XA_1[$icase]}/g"  input_interm2.in  > input_interm1.in
        sed "s/XB/${XB_1[$icase]}/g"  input_interm1.in  > input_interm2.in
      elif [ $iprob -eq 2 ]
      then
        sed "s/XA/${XA_2[$icase]}/g"  input_interm2.in  > input_interm1.in
        sed "s/XB/${XB_2[$icase]}/g"  input_interm1.in  > input_interm2.in
      elif [ $iprob -eq 3 ]
      then
        sed "s/XA/${XA_3[$icase]}/g"  input_interm2.in  > input_interm1.in
        sed "s/XB/${XB_3[$icase]}/g"  input_interm1.in  > input_interm2.in
      elif [ $iprob -eq 4 ]
      then
        sed "s/XA/${XA_4[$icase]}/g"  input_interm2.in  > input_interm1.in
        sed "s/XB/${XB_4[$icase]}/g"  input_interm1.in  > input_interm2.in
      fi

      sed "s/XCEN/.true./g"                input_interm2.in  > input_${iprob}_${icase}_T.in
      cp input_${iprob}_${icase}_T.in input.in
      ./testcd10np
      
      sed "s/XCEN/.false./g"                input_interm2.in  > input_${iprob}_${icase}_F.in
      cp input_${iprob}_${icase}_F.in input.in
      ./testcd10np
    done
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
