#!/bin/bash

#
# Date: 2011-04-07
# Author: Matthias Hoelzl, IPP Garching
#

function usage() {
  echo ""
  echo "Usage: `basename $0` <inputfile> <key1> [...]"
  echo "       `basename $0` -p <inputfile> <key>"
  echo ""
  echo "Purpose: Print parameter values from a JOREK namelist input file."
  echo ""
  echo "Example: `basename $0` input eta nstep_n"
  echo ""
}

if [ $# -lt 2 ] || [ "$1" == "-h" ]; then
  usage
  exit
fi

option_p="0"
if [ "$1" == "-p" ]; then
  option_p="1"
  shift
fi

inputfile="$1"
shift
if [ ! -f $inputfile ]; then
  echo ""
  echo "ERROR: Inputfile '$inputfile' does not exist."
  usage
  exit 1
fi

# --- Show input parameters
if [ $option_p -eq 0 ]; then # (human-readable)
  echo ""
  echo "==============================================================="
  echo "$inputfile"
  echo "---------------------------------------------------------------"
  while [ $# -gt 0 ]; do
    param="$1"
    shift
    NBMATCHES=$(sed -n "/^[ \t]*$param[ \t]*=.*$/p" $inputfile | wc -l)
    if [ $NBMATCHES -eq 0 ]; then
      echo "$param = [NOT FOUND]"
    else
      sed -n "s/^[ \t]*\($param\)[ \t]*=[ \t]*\([^!]*\).*/\1 = \2/p" $inputfile
    fi
  done
  echo "==============================================================="
  echo ""
else # (machine-readable with option -p)
  param="$1"
  NBMATCHES=$(sed -n "/^[ \t]*$param[ \t]*=.*$/p" $inputfile | wc -l)
  if [ $NBMATCHES -eq 0 ]; then
    echo "ERROR: Input parameter $param not found." >&2
    exit 1
  fi
  sed -n "s/^[ \t]*$param[ \t]*=[ \t]*\([^!]*\).*/\1/p" $inputfile
fi
