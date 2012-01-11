#!/bin/bash

#
# Purpose: Set parameters in a namelist input file.
#
# Date: 2011-04-07
# Author: Matthias Hoelzl, IPP Garching
#

function usage() {
  echo ""
  echo "Usage: `basename $0` <inputfile> <key>=<value> [...]"
  echo ""
  echo "Purpose: Manipulate parameter values in a JOREK namelist input file."
  echo ""
  echo "Example: `basename $0` input eta=2.d-6 'nstep_n=20, 40'"
  echo ""
}

function setparam() {
  key=$1
  val="${2%,}" # Strip comma at the end
  NBMATCHES=$(sed -n "/^[ \t]*$1[ \t]*=.*$/p" $inputfile | wc -l)
  if [ $NBMATCHES -eq 1 ]; then
    cp $inputfile ${inputfile}.bck
    cat ${inputfile}.bck | sed -e "s/\(^[ \t]*$key[ \t]*=[ \t]*\)[^!]*/\1$val/" > $inputfile
  else
    echo ""
    echo "ERROR: Found $NBMATCHES entries for '$key' in the input file." >&2
    usage
    exit 1
  fi
}

SCRIPTDIR=`dirname $0`; SCRIPTDIR=`readlink -f $SCRIPTDIR`

if [ $# -lt 2 ] || [ "$1" == "-h" ]; then
  usage
  exit
fi

inputfile="$1"
shift
if [ ! -f $inputfile ]; then
  echo ""
  echo "ERROR: Inputfile '$inputfile' does not exist." >&2
  usage
  exit 1
fi

# --- Set input parameters
allkeys=""
while [ $# -gt 0 ]; do
  thekey="${1%%=*}"
  theval="${1#*=}"
  if [ "$thekey" == "$1" ] || [ -z "$theval" ]; then
    echo ""
    echo "ERROR: Parameter '$1' is not in the form <value>=<key>." >&2
    usage
    exit 1
  fi
  setparam "$thekey" "$theval"
  allkeys="$allkeys $thekey"
  shift
done

$SCRIPTDIR/showinput.sh $inputfile $allkeys
