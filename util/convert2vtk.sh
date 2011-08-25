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
    echo "ABORTING."
    echo ""
  fi
  echo ""
  echo "Waiting for threads to finish..."
  wait
  if [ ! -z "$local_tmp_dir" ]; then
    rm -rf $local_tmp_dir
  fi
  echo "...done."
  exit
}

function usage () {
  echo ""
  echo "Usage: `basename $0` [options] binary infile [extra-files]"
  echo ""
  echo "Converts JOREK restart files into VTK files that can be"
  echo "visualized using, e.g., VISIT or PARAVIEW."
  echo ""
  echo "Options:"
  echo "  -j <nthreads>   Convert parallel using <nthreads> threads [default: serial]."
  echo "  -min <minstep>  Minimum step number to convert"
  echo "  -max <maxstep>  Maximum step number to convert"
  echo "  -dir <dir>      Write vtk files to the specified directory [default: ./vtk]."
  echo ""
  echo "Options passed to jorek2vtk via namelist input (see code for details):"
  echo "  -nsub <nsub>    Number of finite element subdivisions [default: 5]."
  echo "  -i_tor <i_tor>  Select a single toroidal mode (-1 means all) [default: -1]."
  echo "  -i_plane <i_plane> Select the toroidal plane [default: 1]."
  echo ""
  echo "  binary          jorek2vtk(3d) executable"
  echo "  infile          Input file of the corresponding JOREK run"
  echo "  extra-files     Additional files that are required for running"
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
    
    rm -f jorek_restart.rst
    ln -s $file jorek_restart.rst
    for copyfile in $copyfiles; do
      cp $startDir/$copyfile .
    done
    $binary < $infile > ./log 2>&1
    if [ $? -ne 0 ]; then
      if [ ! -f $ERROR_STOP_FILE ]; then
        touch $ERROR_STOP_FILE
        cat ./log
        echo ""
        echo "ithread = $ithread"
        echo "file    = `basename $file`"
        echo ""
        echo "ERRORS OCCURED EXECUTING THE BINARY. SEE ABOVE."
      fi
      unmark_running $ithread
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
dir="./vtk"
nsub=""
i_tor=""
i_plane=""
writenml="no"
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
  elif [ "$1" == "-dir" ]; then
    dir="$2"
    shift 2
  elif [ "$1" == "-nsub" ]; then
    nsub="$2"
    shift 2
    writenml="yes"
  elif [ "$1" == "-i_tor" ]; then
    i_tor="$2"
    shift 2
    writenml="yes"
  elif [ "$1" == "-i_plane" ]; then
    i_plane="$2"
    shift 2
    writenml="yes"
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

binary=`readlink -f $1`
shift
infile=`readlink -f $1`
shift
sourceDir=`readlink -f .`
copyfiles="$@"



# --- Some basic checks
if [ ! -f $binary ]; then
  echo "ERROR: $binary does not exist."
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



# --- Create directory for vtk files
startDir=`pwd`
mkdir -p $dir || exit 1
targetDir=`readlink -f $dir`



# --- Create local temporary directory
local_tmp_dir0="tmp_vtk_$$"
mkdir -p $local_tmp_dir0
local_tmp_dir=`readlink -f $local_tmp_dir0` # absolute path
ERROR_STOP_FILE="$local_tmp_dir/ERROR_STOP"



# --- Create the namelist file with jorek2vtk parameters
if [ "$writenml" == "yes" ]; then
  vtk_nml="$local_tmp_dir0/vtk.nml"
  echo "&vtk_params"               > $vtk_nml
  if [ ! -z "$nsub" ]; then
    echo "  nsub    = $nsub"      >> $vtk_nml
  fi
  if [ ! -z "$i_tor" ]; then
    echo "  i_tor   = $i_tor"     >> $vtk_nml
  fi
  if [ ! -z "$i_plane" ]; then
    echo "  i_plane = $i_plane"   >> $vtk_nml
  fi
  echo "/"                        >> $vtk_nml
  copyfiles="$copyfiles $vtk_nml"
elif [ -f "vtk.nml" ]; then
  # If parameters -nsub, -i_tor, -i_plane were not provided, but
  # a vtk.nml file exists, include it automatically
  copyfiles="$copyfiles vtk.nml"
fi



# --- Prepare thread temporary folders
for i in `seq $nthreads`; do
  tmpdir[$i]="$local_tmp_dir/thread_$i"
  mkdir -p ${tmpdir[$i]}
done



# --- Parallel file conversion
echo ""
files=`ls $sourceDir/jorek?????.rst 2> /dev/null`
for file in $files; do
  if [ -f "$ERROR_STOP_FILE" ]; then cleanup; fi
  ithread=`get_available_thread`
  if [ ! -f "$ERROR_STOP_FILE" ]; then
    mark_running $ithread
    do_convert $file $ithread &
  fi
done



cleanup 0
