#!/bin/bash

#
# Purpose: Analyze a JOREK logfile to identify typical problems
#
# Date: 2019
# Author: Matthias Hoelzl, IPP Garching
#

function usage() {
  echo ""
  echo "Usage: `basename $0` [-h] logfile"
  echo ""
  echo "  -h         Print this usage information"
  echo "  logfile    Name of the logfile to be analyzed"
  echo ""
}

function cleanup_blanks() {
  sed -e 's/^[ \t*]*//' | sed -e 's/[ \t*]*$//' | sed -e 's/  */ /g'
}

function indent() {
  sed -e 's/^/    /'
}

SCRIPTDIR=`dirname $0`; SCRIPTDIR=`readlink -f $SCRIPTDIR`

if [ $# -ne 1 ] || [ "$1" == "-h" ] || [ ! -e "$1" ]; then
  usage
  exit
fi

logfile_orig="$1"

nruns=`grep "JOREK_2.0" $logfile_orig | wc -l`
echo "The logfile contains the output of $nruns JOREK simulations"

rm -f xx??
csplit -s $logfile_orig "/JOREK_2.0/" {*}
mv xx01 xx01.bck; mv xx00 xx01; cat xx01.bck >> xx01; rm xx01.bck

i=0
for logfile in `ls xx??`; do
  i=$[i+1]
  
  echo ""
  echo "************************"
  echo "Run $i"
  echo "************************"
  
  echo "Preprocessor switches:"
  grep "USE" $logfile | cleanup_blanks | grep -v "defines" | indent
  
  grep "jorek_model  " $logfile | cleanup_blanks
  
  grep "freeboundary " $logfile | cleanup_blanks
  
  egrep "n_tor |version of element_matrix" $logfile | tr -d '\n' | sed -e 's/$/\n/' | \
    cleanup_blanks
  
  grep "n_period  " $logfile | cleanup_blanks
  
  grep "n_plane  " $logfile | cleanup_blanks
  
  grep "MPI processes" $logfile | cleanup_blanks
  
  grep "#MPI id" $logfile | cleanup_blanks | indent | sort -n
  
  grep "OpenMP threads" $logfile | cleanup_blanks
  
  grep "restart time  " $logfile | cleanup_blanks
  
  echo "Last time steps:"
  grep "  time step : " $logfile | tail -n 5 | cleanup_blanks | indent
  
  echo "Last number of iterations:"
  grep "info(2)" $logfile | tail -n 5 | cleanup_blanks | indent
  
  grep "iter_precon" $logfile | cleanup_blanks | indent
  
  grep "gmres_max_iter" $logfile | cleanup_blanks | indent
  
  echo "Elapsed time in last time steps:"
  grep "Elapsed" $logfile | grep "  0 " | sed -e 's/^.*time[,]* *//' | tail -n 30 | \
    cleanup_blanks | sed -e 's/^\([^I]\)/    \1/' | indent
    
  echo "Lines with NaNs:"
  grep NaN $logfile | cleanup_blanks | indent
  
  echo "Warning messages:"
  grep -v "are displayed" $logfile | grep -i -A 2 warning | cleanup_blanks | indent
  
  echo "Error messages:"
  grep -v "are displayed" $logfile | sed -e 's/^ *error *[0123456789.e-]* *$//' | \
    egrep -i -A 2 "^ *error|^ *fatal" | cleanup_blanks | indent
done
