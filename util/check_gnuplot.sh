#!/bin/bash

function usage() {
  echo ""
  echo "Usage: `basename $0` <required_version>"
  echo ""
}

if [ $# -lt 1 ] || [ "$1" == -h ]; then
  usage
  exit
fi

required_version="$1"
required_version2=`echo $required_version | sed -e 's/\.//'`

which_gnuplot=`which gnuplot`
if [ ! `which gnuplot` ]; then
  echo ""
  echo "ERROR: You need gnuplot to use this script. If you get this message although"
  echo "you have it installed, you may need to add the installation folder to your"
  echo "PATH environment variable."
  echo ""
  exit 1
fi
gnuplot_version_string=`gnuplot --version`
gnuplot_version=`echo $gnuplot_version_string | cut -d' ' -f 2`
gnuplot_version2=`echo $gnuplot_version | sed -e 's/\.//'`
if [ $gnuplot_version2 -lt $required_version2 ]; then
  echo ""
  echo "ERROR: Gnuplot $required_version or newer is required."
  echo "Found $gnuplot_version_string instead."
  echo ""
  exit 1
fi
