#!/bin/bash

# 
# This script extracts lists of input parameters from all models and creates an overview that can be copied to our wiki.
# This way, comments on the input parameters in the code are directly reflected in the wiki.
# (c) Matthias Hoelzl, 2019
# 

startdir=`pwd`

outfile="parameter-overview.txt"

models=""
for i in `ls -1d models/model*/initialise_parameters.f90`; do
  model=`echo $i | sed -e 's|models/model||' -e 's|/.*||'`
  models="$models $model"
  echo $model
  grep /in1/ models/model$model/initialise_parameters.f90 -A 999 | grep "&" -A 1 | tr -d '\n' | sed -e 's/[&,]/ /g' -e 's/\t/ /g' -e 's/  */ /g' -e 's/,//g' -e 's|^.*/ ||' | tr ' ' '\n' | sort | uniq > tmp_${model}_$$
done

cat tmp_*_$$ | sort | uniq > tmp_$$

function insert_header() {
  echo "^  parameter  ^  default  ^  description  ^  models  ^" >> $outfile
}

rm -f $outfile
insert_header
k=0
for param in `cat tmp_$$`; do
  k=$[k+1]
  m=$(( $k % 20 ))
  if [ $m -eq 0 ]; then
    insert_header
  fi
  description=`egrep " $param[( ]" models/phys_module.f90 | sed -e 's/^.*![< ]*//' -e 's/\\\f//g' | tr '\n' ';' | sed -e 's/;$//'`
  default=`grep " $param[ (=]" models/preset_parameters.f90 | sed -e 's/^.*= *//' -e 's/!.*$//'`
  echo -n "| $param | $default | $description |  " >> $outfile
  for model in $models; do
    matches=`egrep -x $param tmp_${model}_$$ | wc -l | sed -e 's/0//' -e "s/1/$model/"`
    echo -n "$matches " >> $outfile
  done
  echo "|" >> $outfile
done

rm -f tmp*_$$
