#!/bin/bash

#
# Purpose: Plot energies and growth rates during or after a JOREK code run using GNUPLOT
#
# Date: 2011-03-30
# Author: Matthias Hoelzl, IPP Garching
#

function usage() {
  echo ""
  echo "Usage: `basename $0` [-h -f <file> -l -q <qtty> -ps -(no)log -title <title>]"
  echo ""
  echo "  -h                 Print this usage information"
  echo "  -f <file>          Take data from <file> instead of 'macroscopic_vars.dat'"
  echo "  -l                 List plottable quantities"
  echo "  -q <qtty>          Plot the given quantity (default: '-q energies')"
  echo "  -ps                Plot to .ps files (default: plot to screen)"
  echo "  -si                Plot with si units"
  echo "  -no0               Plot without the n=0 mode"
  echo "  -(no)log           (Non-)Logarithmic y-axis (default: '-log')"
  echo "  -title <title>     Add this to title (after given quantity)"
  echo "  -xrange <x0>:<xf>  X range for plotting."
  echo "  -yrange <y0>:<yf>  Y range for plotting."
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
elif [ ! -f "$SCRIPTDIR/plot_live_data.gnuplot" ]; then
  echo "ERROR: The file plot_live_data.gnuplot must be available in the same folder"
  echo "as the plot_live_data.sh script."
  exit 1
fi
#   --- Gnuplot available?
$SCRIPTDIR/check_gnuplot.sh 4.4
if [ $? -ne 0 ]; then
  echo ""
  echo "NOTE: If you have xmgrace installed, you may use the script"
  echo "plot_live_data_grace.sh instead."
  echo ""
  exit 1
fi

# --- Evaluate command line parameters
ps=0
qtty="energies"
si=0
no0=0
logy2=""
file="macroscopic_vars.dat"
title=""
xrange="[:]";  yrange="[:]"
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
  elif [ "$1" == "-title" ] || [ "$1" == '-t' ]; then
    title="$2"
    shift 2
  elif [ "$1" == "-xrange" ]; then
    if [[ $2 = *":"* ]]; then xrange="[$2]"
    else echo "Wrong format for range. Should be #:#"; exit 1; fi
    shift 2
  elif [ "$1" == "-yrange" ]; then
    if [[ $2 = *":"* ]]; then yrange="[$2]"
    else echo "Wrong format for range. Should be #:#"; exit 1; fi
    shift 2
  elif [ "$1" == "-si" ]; then
    si=1
    shift
  elif [ "$1" == "-nosi" ]; then
    si=0
    shift
  elif [ "$1" == "-no0" ]; then
    no0=1
    shift
  elif [ "$1" == "-log" ]; then
    logy2=1
    shift
  elif [ "$1" == "-nolog" ]; then
    logy2=0
    shift
  elif [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    usage
    exit
  elif [ "$1" == "-l" ]; then
    echo ""
    echo "Plottable quantities are:"
    $extract_live_data plottable -f $file | sed -e 's/\s\+/\n/g' | sort | sed -e 's/^/\* /'
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
elif [ $no0 -eq 1 ] && (! [[ $qtty =~ "energies"  ]] && ! [[ $qtty =~ "growth_rates" ]]); then
  #elif [ $no0 -eq 1 ] && ([ "$qtty" != *"energies"* ] && [ "$qtty" != *"growth_rates"* ]); then
  echo ""
  echo "ERROR: -no0 flag can only be used with energies and growth_rates (kinetic and magnetic)."
  exit 1
fi

# --- Extract necessary data from macroscopic_vars.dat
n_tor=`$extract_live_data n_tor -f $file`
ncols=`$extract_live_data n_${qtty} -f $file`
if [ $si -eq 1 ]; then
  xlabel=`$extract_live_data ${qtty}_xlabel_si -f $file`
  ylabel=`$extract_live_data ${qtty}_ylabel_si -f $file`
  qtty_x2si=`$extract_live_data ${qtty}_x2si -f $file`
  qtty_y2si=`$extract_live_data ${qtty}_y2si -f $file`
else
  xlabel=`$extract_live_data ${qtty}_xlabel -f $file`
  ylabel=`$extract_live_data ${qtty}_ylabel -f $file`
  qtty_x2si=1
  qtty_y2si=1
fi
logy=`$extract_live_data ${qtty}_logy -f $file`
if [ -n "$logy2" ]; then
  logy=$logy2
fi
if [ -z "$logy" ]; then
  logy="1"
fi
$extract_live_data ${qtty} ${qtty}.dat -f $file

# --- Plot the quantity
cat $SCRIPTDIR/plot_live_data.gnuplot                                                     \
  | sed -e "s|<ncols>|$ncols|"         -e "s|<ps>|$ps|"         -e "s|<qtty>|$qtty|"      \
        -e "s|<title>|$title|"         -e "s|<no0>|$no0|"       -e "s|<logy>|$logy|"      \
        -e "s|<qtty_x2si>|$qtty_x2si|" -e "s|<ylabel>|$ylabel|" -e "s|<xlabel>|$xlabel|"  \
        -e "s|<qtty_y2si>|$qtty_y2si|" -e "s|<xrange>|$xrange|" -e "s|<yrange>|$yrange|"  \
                                                                       > local_plot.gnuplot

gnuplot local_plot.gnuplot &
gnuplot_pid=$!

# --- Refresh live data as long as gnuplot is running -- pressing E in gnuplot will update
if [ $ps -eq 0 ]; then
  while true; do
    sleep 10
    ps | grep -q $gnuplot_pid # See if Gnuplot is still running
    if [ $? -ne 0 ]; then
      break # If Gnuplot is not running, exit the loop
    fi
    # Refresh the data 
    $extract_live_data ${qtty} ${qtty}.dat.tmp -f ${file}
    mv ${qtty}.dat.tmp ${qtty}.dat
  done &
fi
