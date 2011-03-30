#!/bin/bash

#
# Purpose: Plot energies and growth rates during or after a JOREK code run
#
# Data: 2011-03-30
# Author: Matthias Hoelzl, IPP Garching
#

function usage() {
  echo ""
  echo "Usage: $0 [-h -q <qtty> -ps -(no)log]"
  echo ""
  echo "  -h         Print this usage information"
  echo "  -q <qtty>  Plot the given quantity (default: 'energies')"
  echo "               Currently possible: energies growth_rates"
  echo "  -ps        Plot to .ps files (default: plot to screen)"
  echo "  -(no)log   (Non-)Logarithmic y-axis (default: logarithmic)"
  echo ""
}

SCRIPTDIR=`dirname $0`

# --- Evaluate command line parameters
ps=0
qtty="energies"
logy=1
while [ $# -gt 0 ]; do
  if [ "$1" == "-ps" ]; then
    ps=1
    shift
  elif [ "$1" == "-q" ]; then
    qtty="$2"
    shift 2
  elif [ "$1" == "-log" ]; then
    logy=1
    shift
  elif [ "$1" == "-nolog" ]; then
    logy=0
    shift
  elif [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    usage
    exit
  else
    echo ""
    echo "ERROR: Unknown option: '$1'."
    usage
    exit 1 
  fi
done

# --- Check that required scripts are available
if [ ! -f "$SCRIPTDIR/extract_live_data.sh" ]; then
  echo "ERROR: The script extract_live_data.sh must be available at the same position"
  echo "as the plot_live_data.sh script."
  exit 1
elif [ ! -f "$SCRIPTDIR/plot_live_data.gnuplot" ]; then
  echo "ERROR: The file plot_live_data.gnuplot must be available at the same position"
  echo "as the plot_live_data.sh script."
  exit 1
fi

# --- Extract necessary data from macroscopic_vars.dat first
extract_live_data="$SCRIPTDIR/extract_live_data.sh"
n_tor=`$extract_live_data n_tor -`
ncols=`$extract_live_data n_${qtty} -`
$extract_live_data ${qtty} ${qtty}.dat

# --- Plot energies and growth rates
cat $SCRIPTDIR/plot_live_data.gnuplot                                                     \
  | sed -e "s/<ncols>/$ncols/" -e "s/<ps>/$ps/" -e "s/<qtty>/$qtty/" -e "s/<logy>/$logy/" \
  > local_plot.gnuplot
gnuplot local_plot.gnuplot
