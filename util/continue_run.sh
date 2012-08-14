#!/bin/bash

function usage () {
  echo ""
  echo "Restart a JOREK run."
  echo ""
  echo "Usage: `basename $0` jobscript infile [key=value ...]"
  echo ""
}

SCRIPTDIR=`dirname $0`; SCRIPTDIR=`readlink -f $SCRIPTDIR`

if [ $# -lt 2 ] || [ "$1" == '-h' ] || [ "$1" == '--help' ]; then
  usage
  exit
fi

jobscript="$1"
shift
infile="$1"
shift

setparams="$@ restart=.true."
$SCRIPTDIR/setinput.sh $infile $setparams

newest_restart=`ls -1rt jorek[0-9]*.rst | tail -n 1`
if [ "$newest_restart" == "" ]; then
  echo
  echo "ERROR: Could not find a restart file."
  echo
  exit 1
fi

cp $newest_restart jorek_restart.rst

sbatch $jobscript
