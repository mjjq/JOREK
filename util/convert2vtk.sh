#!/bin/bash

#
# Purpose: Produces .vtk files that can be plotted, e.g., by Paraview from JOREK restart files.
#
# Date: 2011-08-22
# Author: Matthias Hoelzl, IPP Garching
#

trap cleanup 1 2 3 6

function cleanup () {
  if [ "$1" != "0" ]; then
    echo ""
    echo "ABORTING. Waiting for unfinished threads..."
    echo ""
  fi
  wait
  for i in `seq $nthreads`; do
    if [ ! -z "${tmpdir[$i]}" ]; then
      rm -rf ${tmpdir[$i]}
    fi
  done
  if [ "$1" == "0" ]; then
    echo "done."
  fi
  exit
}

function usage () {
  echo ""
  echo "Usage: `basename $0` [options] jorek2vtk infile [extra-files]"
  echo ""
  echo "Converts JOREK restart files into VTK files that can be"
  echo "visualized using, e.g., VISIT or PARAVIEW."
  echo ""
  echo "Options:"
  echo "  -j <nthreads>   Convert parallel using <nthreads> threads [default: serial]."
  echo "  -min <minstep>  Minimum step number to convert"
  echo "  -max <maxstep>  Maximum step number to convert"
  echo ""
  echo "  jorek2vtk       jorek2vtk(3d) executable"
  echo "  infile          Input file of the corresponding JOREK run"
  echo "  extra-files     Additional files that are required for running jorek2vtk"
  echo ""
}

function mark_running () {
  ithread="$1"
  touch ${tmpdir[$ithread]}/ISRUNNING
}

function unmark_running () {
  ithread="$1"
  rm -f ${tmpdir[$ithread]}/ISRUNNING
}

function is_running () {
  ithread="$1"
  if [ -f ${tmpdir[$ithread]}/ISRUNNING ]; then
    echo 'yes'
  else
    echo 'no'
  fi
}

function get_available_thread () {
  while true; do
    for i in `seq $nthreads`; do
      if [ `is_running $i` == 'no' ]; then
        echo "$i"
	return
      fi
    done
    sleep 1
  done
}

function do_convert () {
  file="$1"
  ithread="$2"
  
  cd ${tmpdir[$ithread]}
  
  targetFile=${file##*/} # Remove directory from filename
  stepnum=${targetFile:5:5}
  targetFile="jorek.${targetFile:5:5}.vtk" # Target filename with same number as source
  targetFile="$targetDir/$targetFile" # Target filename with full path
  
  # Convert only new restart files in the range between minstep and maxstep
  if ( [ ! -e $targetFile ] || [ "$file" -nt "$targetFile" ] ) \
    && [ $stepnum -ge $minstep ] && [ $stepnum -le $maxstep ]; then
    
    echo "CONVERTING '${file##*/}'"
    
    cp $file jorek_restart.rst
    for copyfile in $copyfiles; do
      cp $startDir/$copyfiles .
    done
    $jorek2vtk < $infile > ./log
    if [ $? -ne 0 ]; then
      echo "AN ERROR OCCURED!"
      cat ./log
      return 1
    else
      grep -i "warning" ./log
    fi
    mv jorek_tmp.vtk $targetFile
  fi
  
  unmark_running $ithread
}



# --- Process command line parameters
nthreads="1"
minstep="0"
maxstep="99999"
while [ $# -gt 1 ]; do
  if [ "$1" == "-j" ]; then
    nthreads="$2"
    shift 2
  elif [ "$1" == "-min" ]; then
    minstep="$2"
    shift 2
  elif [ "$1" == "-max" ]; then
    maxstep="$2"
    shift 2
  elif [ "$1" == "-h" ] || [ "$1" = "--help" ]; then
    usage
    exit 0
  elif [ "${1:0:1}" == "-" ]; then
    echo "ERROR: Unknown option "$1"."
    usage
    exit 1
  else
    break
  fi
done

if [ $# -lt 2 ]; then
  echo "ERROR: Not enough parameters."
  usage
  exit 1
fi

jorek2vtk=`readlink -f $1`
shift
infile=`readlink -f $1`
shift
sourceDir=`readlink -f .`
copyfiles="$@"



# --- Some basic checks
if [ ! -f $jorek2vtk ]; then
  echo "ERROR: $jorek2vtk does not exist."
  usage
  exit 1
elif [ ! -f $infile ]; then
  echo "ERROR: $infile does not exist."
  usage
  exit 1
elif [ ! -d $sourceDir ]; then
  echo "ERROR: $sourceDir does not exist."
  usage
  exit 1
elif [ ! -d $targetDir ]; then
  echo "ERROR: $targetDir does not exist."
  usage
  exit 1
fi



# --- Create vtk subfolder
startDir=`pwd`
mkdir -p ./vtk
targetDir=`readlink -f ./vtk`



# --- Prepare thread temporary folders
for i in `seq $nthreads`; do
  tmpdir[$i]="/tmp/tmp_c2v_$$_$i"
  mkdir ${tmpdir[$i]}
done



# --- Convert files (in parallel)
files=`ls $sourceDir/jorek?????.rst 2> /dev/null`

for file in $files; do
  ithread=`get_available_thread`
  mark_running $ithread
  do_convert $file $ithread &
done



cleanup 0
