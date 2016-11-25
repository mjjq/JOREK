#!/bin/bash
# 
# This script prepares a complete JOREK simulation:
#   * Compile the binaries
#   * Prepare the input files
#   * Prepare the job script files
# 
# You need to create a configuration file which contains all information about the simulation.
# An example is found in 'trunk/util/jobscript_templates/EXAMPLE/'. To test the example, execute
# '../../prepare_run.sh config.run' in the EXAMPLE folder.
# 
# The jobscript template is taken from trunk/util/jobscript_templates/. The input namelist
# template 'input.template' must be provided by the user.
# 
# Date: 2012-01-09
# Author: Matthias Hoelzl, IPP Garching
# 

function usage() {
  echo ""
  echo "Usage:"
  echo "  `basename $0` <configuration file>"
  echo ""
}

start_dir=`pwd`
SCRIPTDIR=`dirname $0`; SCRIPTDIR=`readlink -f $SCRIPTDIR`
tmp="/tmp/tmp_pr_$$"
if [ $# -ne 1 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  usage && exit
fi

function add_param() {
  the_key="${1%%=*}"
  the_val="${1#*=}"
  eval "$prefix$the_key=$the_val"
  params="$params $the_key"
}

function fillin() {
  sed -i -e "s|$2|$3|g" $1
}

# --- Load the configuration file
echo ""
echo "Loading configuration file."
PREPARE_RUN_ALREADY_STARTED="yes"
config_file="$1"
source $config_file || exit 1
template_dir=`readlink -f $template_dir`
echo "  machine=$machine"
echo "  nsteps=$nsteps"
echo "done."
echo ""

# --- Copy include files
echo "Copying include files."
cd $template_dir
echo "  $include_files"
if [ -n "$include_files" ]; then
  cp $include_files $start_dir || exit 1
fi
echo "done."
echo ""

# --- Create link to util folder
echo "Creating link to util folder."
cd $start_dir
rm -f util
ln -s $SCRIPTDIR util || exit 1
echo "done."
echo ""
  
# --- Prepare the jobsteps
cd $start_dir
for job_step in `seq $nsteps`; do
  echo "Preparing job step $job_step."
  
  # --- Get parameters for this job step
  source $config_file
  template_dir=`readlink -f $template_dir`
  code_dir=`readlink -f $code_dir`
  
  # --- Set parameters in Makefile.inc and mod_parameters.f90; compile the binary
  echo "  Compiling binary."
  cd $code_dir
  $SCRIPTDIR/config.sh $code_params >> /dev/null || exit 1
  jorek_model=`$SCRIPTDIR/config.sh -p model`
  binary="jorek_model${jorek_model}_step$job_step"
  rm -f $start_dir/$binary; echo "    $binary"
  if [ "$binaries_from_template_folder" == "1" ]; then
    echo "    COPYING BINARIES FROM TEMPLATE FOLDER!"
    cp $template_dir/$binary $start_dir
  else
    if [ "$code_params_old" == "$code_params" ]; then
      cd $start_dir
      binary_previous="jorek_model${jorek_model}_step$[job_step-1]"
      echo "    Using $binary_previous (same parameters; soft link)."
      ln -s $binary_previous $binary || exit 1
    else
      make cleanall > /dev/null && make -j $compile_threads > /dev/null || exit 1
      cp jorek_model${jorek_model} $start_dir/$binary || exit 1
    fi
    code_params_old="$code_params"
  fi
  echo "  done."
  echo ""
  
  # --- Prepare the input file
  echo "  Preparing input."
  cd $start_dir
  input="input$job_step"
  rm -f $input; echo "    $input"
  cp $template_dir/input.template $input || exit 1
  if [ $job_step -gt 1 ]; then restart=".true."; else restart=".false."; fi
  $SCRIPTDIR/setinput.sh $input restart=$restart $input_params > /dev/null || exit 1; rm -f ${input}.bck
  fillin $input "<<input_comment>>" "Created by prepare_run.sh `date '+%F %T'` from the template $template_dir/input."
  echo "  done."
  echo ""
  
  # --- Prepare the jobscript
  echo "  Preparing jobscript."
  cd $start_dir
  jobscript="jobscript$job_step"
  rm -f $jobscript; echo "    $jobscript"
  cp $SCRIPTDIR/jobscript_templates/$machine/job.template $jobscript
  fillin $jobscript "<<jobscript_comment>>" "Created by prepare_run.sh `date '+%F %T'` for the machine $machine."
  prefix="__JP__"; params=""
  for param in $job_params "main_dir=$start_dir" "job_step=$job_step" "model=$jorek_model"; do
    add_param "$param"
  done
  source $SCRIPTDIR/jobscript_templates/$machine/machine.sh || exit 1
  for param in $params; do
    value=`eval echo "\\\${$prefix$param}"`
    fillin $jobscript "<<$param>>" "$value" || exit 1
  done
  echo "  done."
  echo ""
  
  echo "done."
  echo ""
  
done

# --- Compile diagnostic tools
echo "Compiling diagnostic tools."
echo "  $diagnostics"
cd $code_dir
if [ "$binaries_from_template_folder" == "1" ]; then
  echo "    COPYING BINARIES FROM TEMPLATE FOLDER!"
  cd $template_dir
  cp $diagnostics $start_dir
else
  make -j $compile_threads $diagnostics > /dev/null || exit 1
  cp $diagnostics $start_dir
fi
echo "done."
echo ""
