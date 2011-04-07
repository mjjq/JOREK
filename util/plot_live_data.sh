#!/bin/bash

#
# Purpose: Plot energies and growth rates during or after a JOREK code run
#
# Data: 2011-03-30
# Author: Matthias Hoelzl, IPP Garching
#

function usage() {
  echo ""
  echo "Usage: `basename $0` [-h -q <qtty> -ps -(no)log]"
  echo ""
  echo "  -h         Print this usage information"
  echo "  -f <file>  Take data from <file> instead of 'macroscopic_vars.dat'"
  echo "  -l         List plottable quantities"
  echo "  -q <qtty>  Plot the given quantity (default: '-q energies')"
  echo "  -ps        Plot to .ps files (default: plot to screen)"
  echo "  -(no)log   (Non-)Logarithmic y-axis (default: '-log')"
  echo ""
  echo "Remarks:"
  echo "  * The beginning of a quantity name is sufficient, e.g., the command"
  echo "    '`basename $0` -q gr' will plot growth_rates."
  echo ""
}

SCRIPTDIR=`dirname $0`

# --- Some sanity checks
#   --- Required scripts available?
extract_live_data="$SCRIPTDIR/extract_live_data.sh"
if [ ! -f "$extract_live_data" ]; then
  echo "ERROR: The script extract_live_data.sh must be available in the same folder"
  echo "as the plot_live_data.sh script."
  exit 1
elif [ ! -f "$SCRIPTDIR/plot_live_data.gnuplot" ]; then
  echo "ERROR: The file plot_live_data.gnuplot must be available in the same folder"
  echo "as the plot_live_data.sh script."
  exit 1
fi
#   --- Gnuplot available?
tmp="/tmp/tmp_pld_$$"
gnuplot --version > $tmp
if [ $? -ne 0 ]; then
  echo "ERROR: You need to have gnuplot to use this script. If you have it installed,"
  echo "you may need to add the installation folder to your PATH environment variable."
  exit 1
fi
echo ""
echo "Gnuplot found: `cat $tmp`"
echo "(`which gnuplot`)"
rm $tmp

# --- Evaluate command line parameters
ps=0
qtty="energies"
logy=1
file="macroscopic_vars.dat"
while [ $# -gt 0 ]; do
  if [ "$1" == "-ps" ]; then
    ps=1
    shift
  elif [ "$1" == "-f" ]; then
    file="$2"
    shift 2
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
  elif [ "$1" == "-l" ]; then
    plottable=`$extract_live_data plottable -f $file`
    echo ""
    echo "Plottable quantities are: $plottable"
    echo ""
    exit
  else
    echo ""
    echo "ERROR: Unknown option: '$1'."
    usage
    exit 1 
  fi
done

# --- Check that the selected quantity is marked plottable
plottable=`$extract_live_data plottable -f $file`
is_plottable=0
for s in $plottable; do
  if [[ $s =~ ^${qtty} ]]; then
    is_plottable=1
    qtty=$s
  fi
done
if [ $is_plottable -eq 0 ]; then
  echo ""
  echo "ERROR: Quantity '$qtty' is not plottable."
  echo "Plottable quantities are: $plottable"
  usage
  exit 1
fi

# --- Extract necessary data from macroscopic_vars.dat
n_tor=`$extract_live_data n_tor -f $file`
ncols=`$extract_live_data n_${qtty} -f $file`
$extract_live_data ${qtty} ${qtty}.dat -f $file

# --- Plot the quantity
cat $SCRIPTDIR/plot_live_data.gnuplot                                                     \
  | sed -e "s/<ncols>/$ncols/" -e "s/<ps>/$ps/" -e "s/<qtty>/$qtty/" -e "s/<logy>/$logy/" \
  > local_plot.gnuplot
gnuplot local_plot.gnuplot
