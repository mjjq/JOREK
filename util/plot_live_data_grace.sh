#!/bin/bash

#
# Purpose: Plot energies and growth rates during or after a JOREK code run using GRACE
#
# Date: 2011-04-08
# Author: Matthias Hoelzl, IPP Garching
#

function usage() {
  echo ""
  echo "Usage: `basename $0` [-h -f <file> -l -q <qtty> -(no)log]"
  echo ""
  echo "  -h         Print this usage information"
  echo "  -f <file>  Take data from <file> instead of 'macroscopic_vars.dat'"
  echo "  -l         List plottable quantities"
  echo "  -q <qtty>  Plot the given quantity (default: '-q energies')"
  echo "  -(no)log   (Non-)Logarithmic y-axis (default: '-log')"
  echo ""
  echo "Remarks:"
  echo "  * The beginning of a quantity name is sufficient, e.g., the command"
  echo "    '`basename $0` -q gr' will plot growth_rates."
  echo ""
}

SCRIPTDIR=`dirname $0`; SCRIPTDIR=`readlink -f $SCRIPTDIR`

# --- Some sanity checks
#   --- Required scripts available?
extract_live_data="$SCRIPTDIR/extract_live_data.sh"
if [ ! -f "$extract_live_data" ]; then
  echo "ERROR: The script extract_live_data.sh must be available in the same folder"
  echo "as the plot_live_data.sh script."
  exit 1
fi
#   --- Grace available?
if [ ! `which xmgrace` ]; then
  echo ""
  echo "ERROR: You need to have xmgrace to use this script. If you have it installed,"
  echo "you may need to add the installation folder to your PATH environment variable."
  echo ""
  echo "NOTE: If you have a recent gnuplot version installed, you may use the script"
  echo "plot_live_data.sh instead."
  echo ""
  exit 1
fi

# --- Evaluate command line parameters
qtty="energies"
logy=1
file="macroscopic_vars.dat"
while [ $# -gt 0 ]; do
  if [ "$1" == "-f" ]; then
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
xlabel=`$extract_live_data ${qtty}_xlabel -f $file`
ylabel=`$extract_live_data ${qtty}_ylabel -f $file`
$extract_live_data ${qtty} -f $file | sed -e 's/\(^ *[^0-9 ]\)/#\1/' > ${qtty}_grace.dat

options="-nxy"
if [ $logy -eq 1 ]; then
  options="-log y $options"
  echo "NOTE: You may have to zoom out manually in GRACE as the program cannot autozoom"
  echo "  with a logarithmic axis if there are zero values in the data."
fi

# --- Plot with grace
xmgrace $options ${qtty}_grace.dat
