#!/bin/bash

#
# Data: 2011-04-07
# Author: Matthias Hoelzl, IPP Garching
#

function usage() {
  echo ""
  echo "Usage: `basename $0` <inputfile> <key> [...]"
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

inputfile="$1"
shift
if [ ! -f $inputfile ]; then
  echo ""
  echo "ERROR: Inputfile '$inputfile' does not exist."
  usage
  exit 1
fi

# --- Show input parameters
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
