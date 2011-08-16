#!/bin/bash

#
# Purpose: Set the model in the makefile and/or certain parameters in the respective
#   mod_parameters file.
#
# Data: 2011-04-07
# Author: Matthias Hoelzl, IPP Garching
#


function usage() {
  echo ""
  echo "Usage: `basename $0` <key>=<value> [...]"
  echo ""
  echo "Purpose: Manipulate model in config.in and/or Makefile.inc and"
  echo "  the parameters n_tor, n_period, and n_plane in the corresponding"
  echo "  mod_parameters file."
  echo ""
  echo "If called witout parameters, the current configuration is printed."
  echo ""
  echo "Example: `basename $0` model=302 n_tor=3 n_period=8 n_plane=4"
  echo ""
  echo "Remark: Call this script from the JOREK trunk."
  echo ""
}

tmp="/tmp/tmp_sc_$$"

SCRIPTDIR=`dirname $0`

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  usage
  exit
fi

function key() {
  echo "`echo $1 | sed -e 's/=.*$//'`"
}

function val() {
  echo "`echo $1 | sed -e 's/^.*=//'`"
}

function setmodel() {
  model=$1

  # --- Basic checks for the specified model
  if [ ${#model} -eq 3 ]; then
    model="model$model"
  fi
  if [ ! ${#model} -eq 8 ] || [[ ! ${model:5:3} =~ ^[0-9]+$ ]]; then
    echo "ERROR: Illegal model specified: '$model'."
  fi

  # --- Set MODEL = modelXXX in the makefile configuration files
  for file in config.in Makefile.inc; do
    if [ -f $file ]; then
      cp $file $tmp
      cat $tmp | sed -e "s/\(^ *MODEL *= *\)[^ ]*\(.*$\)/\1$model\2/" > $file
    fi
  done
}

function getmodel() {
  model=""
  for file in config.in Makefile.inc; do
    if [ -f $file ]; then
      model2=`egrep "MODEL *= *model[0-9]*" $file | sed -e "s/^ *MODEL *= *\(model[0-9]*\).*$/\1/"`
      if [ "$model" != "" ]; then
        if [ "$model" != "$model2" ]; then
          echo "ERROR: Models in makefile configuration files do not agree." >&2
          echo "Call setconfig.sh model=modelXXX to fix that." >&2
          model="AMBIGUOUS"
        fi
      else
        model=$model2
      fi
    fi
  done
  echo $model
}

function setparam() {
  key=$1
  val=$2
  file="models/$model/mod_parameters.f90"
  if [ ! -f $file ]; then
    echo "ERROR: File '$file' does not exist."
    exit 1
  fi
  cp $file $tmp
  cat $tmp | sed -e "s/\(^.*:: *$key *= *\)[^! ]*\(.*$\)/\1$val\2/" > $file
}

function print_info() {
  config_file="$1"
  echo ""
  echo "======================="
  echo "$config_file"
  echo "-----------------------"
  
  if [ -f "$config_file" ]; then
    
    model=`egrep "MODEL *= *model[0-9]*" "$config_file" | sed -e 's/^ *MODEL *= *//' -e 's/ *#.*$//'`
    echo "  $model"
    modeldir="models/$model"
    param="$modeldir/mod_parameters.f90"
    if [ -d "$modeldir" ]; then
      egrep "integer.*n_tor"    $param | sed -e 's/^ *//' -e 's/ *!.*$//' -e 's/integer *, *parameter *:: *//' -e 's/  */ /g' -e 's/^/  /'
      egrep "integer.*n_period" $param | sed -e 's/^ *//' -e 's/ *!.*$//' -e 's/integer *, *parameter *:: *//' -e 's/  */ /g' -e 's/^/  /'
      egrep "integer.*n_plane"  $param | sed -e 's/^ *//' -e 's/ *!.*$//' -e 's/integer *, *parameter *:: *//' -e 's/  */ /g' -e 's/^/  /'
    else
      echo "WARNING: Directory $modeldir does not exist."
    fi
  else
    echo "(file not found)"
  fi
  echo "======================="
  echo ""
}

# --- First set the model (if there is a respective command line argument)
for arg in $@; do
  if [ `key $arg` == "model" ]; then
    setmodel `val $arg`
  fi
done

# --- Determine the model
model=`getmodel`

if [ "$model" = "AMBIGUOUS" ]; then
  exit 1
fi

# --- Then set the parameters
for arg in $@; do
  if [ `key $arg` != "model" ]; then
    setparam `key $arg` `val $arg`
  fi
done

# --- Print the configuration
for file in config.in Makefile.inc; do
  print_info $file
done

rm -f $tmp
