#!/bin/bash

#
# Purpose: Plot JOREK grids
#
# Date: 2011-04-12
# Author: Matthias Hoelzl, IPP Garching
#

function usage() {
  echo ""
  echo "Usage: `basename $0` [options]"
  echo ""
  echo "  -h         Print this usage information"
  echo "  -o <grid>  Plot only the specified grid"
  echo "  -ps        Plot to a .ps file (default: plot to screen)"
  echo "  -png       Plot to a .png file (default: plot to screen)"
  echo "  -r <resol> Resolution of a png plot (default: '800x800')"
  echo ""
  echo "Example: `basename $0` -o xpoint -png -r 1600x1600"
  echo ""
}

SCRIPTDIR=`dirname $0`; SCRIPTDIR=`readlink -f $SCRIPTDIR`

# --- Some sanity checks
#   --- Gnuplot available?
$SCRIPTDIR/check_gnuplot.sh 4.2
if [ $? -ne 0 ]; then
  exit 1
fi

# --- Evaluate command line parameters
plotto="screen"
plotonly=""
resol="800,800"
while [ $# -gt 0 ]; do
  if [ "$1" == "-ps" ]; then
    plotto="ps"
    shift
  elif [ "$1" == "-png" ]; then
    plotto="png"
    shift
  elif [ "$1" == "-r" ]; then
    resol="$2"
    resol=`echo $resol | sed -e 's/x/,/'`
    shift 2
  elif [ "$1" == "-o" ]; then
    plotonly="$2"
    shift 2
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

tmp="./local_plot.gnuplot"
rm -f $tmp
if [ "$plotto" == "screen" ]; then
  echo "# plot to screen"                      >> $tmp
elif [ "$plotto" == "ps" ]; then
  echo "set term postscript enhanced color"    >> $tmp
  echo "set output 'grids.ps'"                 >> $tmp
elif [ "$plotto" == "png" ]; then
  echo "set term png size $resol crop"         >> $tmp
  echo "set output 'grids.png'"                >> $tmp
  echo "set border lw 2"                       >> $tmp
else
  echo "ERROR: Unexpected error: Illegal value plotto='$plotto'."
  exit 1
fi
echo "set size ratio -1"                       >> $tmp
echo "set title 'JOREK GRIDS'"                 >> $tmp
echo "set xlabel 'R'"                          >> $tmp
echo "set ylabel 'Z'"                          >> $tmp
echo "set style line 1 lt 1 lc rgb '#000000'"  >> $tmp
echo "set style line 2 lt 1 lc rgb '#0000dd'"  >> $tmp
echo "set style line 3 lt 1 lc rgb '#00cc00'"  >> $tmp
echo "plot \\"                                 >> $tmp

if [ "$plotonly" == "" ]; then
  grids=`ls -rt1 grid_*.dat` 2> /dev/null
else
  grids=`ls -rt1 grid_*$plotonly*.dat` 2> /dev/null
fi
if [ -z "$grids" ]; then
  echo ""
  echo "ERROR: No grid data found."
  usage
  exit 1
fi
i=0
for grid in $grids; do
  i=$[i+1]
  if [ $i -gt 1 ]; then
    echo ",  \\"                               >> $tmp
  fi
  echo -n " '$grid' w l ls $i"                 >> $tmp
done
echo ""                                        >> $tmp

if [ "$plotto" == "screen" ]; then
  echo "pause mouse button1"                   >> $tmp
fi

gnuplot $tmp
