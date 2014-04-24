#!/bin/bash

#
# Purpose: Set the model in the makefile and/or certain parameters in the respective
#   mod_parameters file.
#
# Date: 2011-04-07
# Author: Matthias Hoelzl, IPP Garching
#

function usage() {
  echo ""
  echo "Purpose: Modify or print physics model in config.in and/or Makefile.inc and"
  echo "  parameters like n_tor in the corresponding mod_parameters file."
  echo ""
  echo "Usage: `basename $0` [<key1>=<value1> [...]]   Modify model and/or parameters"
  echo "       `basename $0` -p <key>                  Print the value for <key> and exit"
  echo ""
  echo "Examples:"
  echo "  `basename $0` model=302 n_tor=3 n_period=8 n_plane=4"
  echo "  `basename $0` -p model"
  echo ""
}

SCRIPTDIR=`dirname $0`; SCRIPTDIR=`readlink -f $SCRIPTDIR`
make_config_files=`ls config.in Makefile.inc 2>/dev/null`
params="n_tor n_period n_plane n_vertex_max n_nodes_max n_elements_max n_boundary_max n_pieces_max"

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  usage && exit
elif [ -z "$make_config_files" ]; then
  echo "Could not find a make configuration file. Are you in the JOREK trunk?" >&2
  exit 1
fi

function key() {
  echo $1 | sed -e 's/=.*$//' 
}

function val() {
  echo $1 | sed -e 's/^.*=//' 
}

function setmodel() {
  model=$1
  # --- Some checks
  if [ ${#model} -eq 3 ]; then
    model="model$model"
  elif [ ! ${#model} -eq 8 ] || [[ ! ${model:5:3} =~ ^[0-9]+$ ]]; then
    echo "ERROR: Illegal model specified: '$model'." >&2
  fi
  # --- Set model in makefile configuration files
  for file in $make_config_files; do
    sed -i -e "s/\(^ *MODEL *= *\)[^ ]*\(.*$\)/\1$model\2/" $file
  done
}

function getmodel() {
  model=""
  for file in $make_config_files; do
    model2=`egrep "MODEL *= *model[0-9]*" $file | sed -e "s/^ *MODEL *= *\(model[0-9]*\).*$/\1/"`
    if [ "$model" != "" ]; then
      if [ "$model" != "$model2" ]; then
        echo "ERROR: Models in makefile configuration files do not agree. Call 'config.sh model=modelXXX' to fix that." >&2
        exit 1
      fi
    else
      model=$model2
    fi
  done
  echo $model
}

function setparam() {
  key=$1
  val=$2
  matches=`grep -c ":: *$key" $param_file`
  if [ "$matches" -ne 1 ]; then
    echo "ERROR: Could not set parameter $key in $param_file." >&2
    exit 1
  else
    sed -i -e "s/\(^.*:: *$key *= *\)[^ !\t]*\(.*$\)/\1$val\2/" $param_file
  fi
}

function getparam() {
  key=$1
  grep ":: *$key[ =]" $param_file | sed -e "s/^.*:: *$key *= *\([^ !\t]*\).*$/\1/"
}

function print_info() {
  file="$1"
  echo ""
  echo "======================="
  echo "$file:"
  
  if [ -f "$file" ]; then
    echo "-----------------------"
    echo "  `getmodel`"
    for param in $params; do
      echo "  $param = `getparam $param`"
    done
  else
    echo "  (file not found)"
  fi
  echo "======================="
  echo ""
}

function check_param_file() {
  param_file="models/$model/mod_parameters.f90"
  if [ ! -f $param_file ]; then
    echo "ERROR: File '$param_file' does not exist." >&2
    exit 1
  fi
}

# --- Determine the model
model=`getmodel` && check_param_file

# --- If argument -p is given, just print the requested parameter value and exit
if [ "$1" == "-p" ]; then
  if [ "$2" == "model" ]; then
    echo `getmodel | sed -e 's/model//'`
  else
    value=`getparam $2`
    if [ -z "$value" ]; then
      echo "ERROR: Could not find parameter '$2'." >&2
      exit 1
    fi
    echo "$value"
  fi
  exit
fi

# --- First set the model (if it is specified as a command line argument)
for arg in $@; do
  if [ `key $arg` == "model" ]; then
    setmodel `val $arg`
  fi
done
model=`getmodel` && check_param_file

# --- Set the parameters
for arg in $@; do
  if [ `key $arg` != "model" ]; then
    setparam `key $arg` `val $arg`
  fi
done

# --- Print the configuration
for file in $make_config_files; do
  print_info $file
done
echo "('`basename $0` -h' for help)"
