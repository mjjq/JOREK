#!/bin/bash

#
# Purpose: Show which model is configured in the makefiles and what parameters are set in
#   the respective mod_parameters file.
#
# Data: 2011-04-07
# Author: Matthias Hoelzl, IPP Garching
#

function usage() {
  echo ""
  echo "Usage: `basename $0`"
  echo ""
  echo "Purpose: Print model from config.in and/or Makefile.inc and"
  echo "  the parameters n_tor, n_period, and n_plane in the corresponding"
  echo "  mod_parameters file."
  echo ""
  echo "Remark: Call this script from the JOREK trunk."
  echo ""
}

# --- Evaluate command line parameters
while [ $# -gt 0 ]; do
  if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    usage
    exit
  else
    echo ""
    echo "ERROR: Unkown option '$1'."
    usage
    exit 1
  fi
done

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

for file in config.in Makefile.inc; do
  print_info $file
done
