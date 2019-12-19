#!/bin/bash

#
# Purpose: Prints out a list of restart files for given step numbers or times
#          This script is required by convert2*.sh. It can also be used to specify
#          restart files, the user wants to copy or move e.g.
#
# Date: 2019-11-03
# Author: Fabian Wieschollek, IPP Garching
#

function usage() {
  echo ""
  echo "Choose JOREK restart files for given step numbers or times"
  echo ""
  echo "Usage: `basename $0` [options]"
  echo ""
  echo "  -h                          Print this help text"
  echo "  -only <step>,<step>         Convert only listed time steps"
  echo "  -only <step>-<step>         Convert only time steps in the given range"
  echo "  -only <step>-<dstep>-<step> Convert only time steps in the given range with given interval"
  echo "  -donly <dstep>              Equivalent to -only 0-<dstep>-99999"
  echo "  -time <time>,<time>         Selects time step roughly at <time> (JOREK-units)"
  echo "  -time <time>-<dtime>-<time> Selects time step within given time range with given interval"
  echo "  -dtime <dtime>              Equivalent to -time 0-<dtime>-infinity"
  echo "  -ms                         -time is given in milli seconds instead of in JOREK-units"
  echo "  -s                          Prints only the selected stepnumbers, otherwise the full paths (default)"
  echo "  -l                          Creates a file containing all selected timesteps and times,"
  echo "                              if parameter -(d)time is used (default:off)"
  echo ""
  echo "Remarks:"
  echo "  * Option -only will only select those restart files,"
  echo "    whose stepnumber are explictly given by the user."
  echo "  * Option -time will try to identifiy for each time step a restart file,"
  echo "    that matches the time most likely."
  echo ""
}

# --- Process command line parameters
selected_steps="0-99999"
selected_times=""
full_path=1
ms=0
list=0
while [ $# -gt 0 ]; do
  if [ "$1" == "-only" ]; then
    selected_steps="$2"
    shift 2
  elif [ "$1" == "-donly" ]; then
    selected_steps="0-$2-99999"
    shift 2
  elif [ "$1" == "-dtime" ]  ; then
    selected_times="0-$2-infinity"
    shift 2
  elif [ "$1" == "-time" ]; then
    selected_times="$2"
    shift 2
  elif [ "$1" == "-ms" ]; then
    ms=1
    shift 1
  elif [ "$1" == "-s" ]; then
    full_path=0
    shift 1 
  elif [ "$1" == "-l" ]; then
    list=1
    shift 1 
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



# --- Some parameter checks
if [ $# -gt 0 ]; then
  echo "ERROR: Unkown paramaters given."
  usage
  exit 1
fi
if [ -z "$selected_times" ] && ( [ $ms -eq 1 ] || [ $list -eq 1 ] ); then
  echo "WARNING: -l and -ms parameters will be ignored, if -(d)time is not set."
fi
regexp_steps="^[0-9]{1,5}(-[0-9]{1,5}){0,2}(,[0-9]{1,5}(-[0-9]{1,5}){0,2})*$"
if [[ ! "$selected_steps" =~ $regexp_steps   ]]; then
  echo "ERROR: -(d)only-parameter given in wrong format."
  usage
  exit 1
fi
if [ ! -z $selected_times ]; then
  rfloat="[0-9]*(\.)?[0-9]*" #regexp for an integer or real
  regexp_times="^"$rfloat"(-"$rfloat"-("$rfloat"|infinity))?(,"$rfloat"(-"$rfloat"-("$rfloat"|infinity))?)*$"
  if [[ ! "$selected_times" =~ $regexp_times ]]; then
    echo "ERROR: -(d)time-parameter given in wrong format."
    usage
    exit 1
  fi
fi



SCRIPTDIR=`dirname $0`; SCRIPTDIR=`readlink -f $SCRIPTDIR`
export sourceDir=`readlink -f .`




if [ -z $RST_TYPE ] ; then
  . ${SCRIPTDIR}/detect_rst_type.sh > /dev/null 2>&1
  if [ "$RST_TYPE" != "h5" ] && [ "$RST_TYPE" != "rst" ]; then
    echo "ERROR: RST_TYPE not detected properly: $RST_TYPE"
    usage
    exit 1
  fi
fi



# --- Creates list of restart files
avail_steps=`find . -regextype sed -regex  "\./jorek[0-9]*\.h5*" | \
sort -n |  sed  -e "s/\.\/jorek\([0-9]*\)\.h5/\1/g" `



# --- Checks, if restart files are present, 
#     and if macroscopic_vars are present, if time points are selected
if  [ -z "$avail_steps" ] || ( [ ! -f "./macroscopic_vars.dat" ] && [ ! -z $selected_times ] ) ; then
  echo "ERROR: No restart files or macroscopic_vars.dat could be detected."
  usage
  exit 1
fi



# --- Gets normalization value
unit=`${SCRIPTDIR}/extract_live_data.sh sqrt_mu0_rho0`


# --- If time points are given, select restart files via the Python-script
if [ ! -z $selected_times ]; then
  name_times="times_tmp_$$"
  name_steps="steps_tmp_$$"
  echo "0." > $name_times #Set t=0. for the 0th Step
  ${SCRIPTDIR}/extract_live_data.sh times | tail -n +2 | sed  -e "s/[0-9 ]* //g"  >> $name_times
  echo $avail_steps > $name_steps
  python ${SCRIPTDIR}/select_restart_files.py $name_times $name_steps $full_path $selected_times $list $unit $ms
  rm $name_times
  rm $name_steps
	
# --- If step numbers are given, select restart files as below
else
  for avail_step in $avail_steps; do
    step_number=`echo $avail_step | sed -e 's/^[0 ]*\(.*.\)$/\1/'` # remove leading zeros
    step_ranges=`echo $selected_steps | tr ',' ' '` # split selected_steps, e.g., 1-3,5-7 -> 1-3 5-7
    for step_range in $step_ranges; do
      step_numbers=(`echo $step_range | tr '-' ' '`) # split step_range, e.g., 1-3 -> 1 3
      if [[ ${#step_numbers[*]} -eq 1 && ${step_numbers[0]} -eq $step_number ]] || \
         [[ ${#step_numbers[*]} -eq 2 && ${step_numbers[0]} -le $step_number && ${step_numbers[1]} -ge $step_number ]] || \
         [[ ${#step_numbers[*]} -eq 3 && ${step_numbers[0]} -le $step_number && $(($step_number % ${step_numbers[1]})) -eq 0 && ${step_numbers[2]} -ge $step_number ]] ; then

        # --- Prints list of absolute paths of the selected restart files or only their step numbers
        if [ $full_path -eq 1 ] ; then
          echo "$sourceDir/jorek$avail_step.${RST_TYPE}"
        else
          echo $avail_step
        fi
      fi
    done
  done
fi
