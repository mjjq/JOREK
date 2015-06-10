#!/bin/bash


# This script carries out a non regression test.
# * run_test.sh -h prints usage information.
# * It compiles the code.
# * It runs the test case.
# * It extracts data using jorek_extract_data.
# * It compares the data to the reference data of the test case.
# * The return code will be zero for a passed test, otherwise non-zero.


startdir=`pwd`
codedir=`readlink -f ..` # Assumption about source code location

# --- Define default configuration for launching MPI runs
if [ -z "$PRERUN" ]; then
    export PRERUN='echo "MPI configuration"'
fi
if [ -z "$MPIRUN" ]; then
    export MPIRUN="mpirun -n"
fi

# --- Usage printing function
function printusage() {
  echo ""
  echo "Usage:"
  echo "  $0 [options] testcase"
  echo ""
  echo "Options:"
  echo "  -h            Print this help information"
  echo "  -k            Keep temporary run directory"
  echo "  -l            List available test cases."
  echo "  -n            Do not compile (assume executables already exist)"
  echo "  -t tempdir    Specify a temp directory used for the test run"
  echo "                (default: current directory)"
  echo ""
}

# --- Process command line options
testcase="NONE" # (preset) 
compile="yes"   # (preset)
keep="no"       # (preset)
tmpdir="$startdir/$$"
while [ $# -gt 0 ]; do
  option="$1"
  if [ "$option" == "-h" ]; then
    printusage
    exit 1
  elif [ "$option" == "-k" ]; then
    keep="yes"
    shift
  elif [ "$option" == "-l" ]; then
    echo ""
    echo "Available test cases:"
    cases=`ls -1 testcases`
    for case in $cases; do
      echo ""
      echo "*** $case ***"
      source testcases/$case/settings.sh
      echo "$description"
    done
    echo ""
    exit 1
  elif [ "$option" == "-n" ]; then
    compile="no"
    shift
  elif [ "$option" == "-t" ]; then
    tmpdir="$2/$$"
    shift 2
  elif [ "${option:0:1}" != "-" ]; then
    if [ "${testcase}" != "NONE" ]; then
      echo ""
      echo "ERROR: Two test case names not supported at present."
      printusage
      exit 1
    else
      testcase="$1"
    fi
    shift
  fi
done

# --- Check if the testcase really exists
if [ ! -d  "testcases/$testcase" ]; then
  echo ""
  echo "ERROR: Test case '$testcase' does not exist. Valid test cases are:"
  cases="`ls -1 testcases | tr '\n' ' ' | sed -e 's/  */ /g'`"
  echo "  $cases"
  echo ""
  printusage
  exit 1
fi
testcasedir=`readlink -f testcases/$testcase`


# --- Read test case information
relativeaccuracy="0.0" # (preset)
absoluteaccuracy="0.0" # (preset)
source $testcasedir/settings.sh


# --- Set hard-coded parameters and compile
#Remark(GL): We need several executables sometime for one single run (ntor=1, ntor=XX), 
#Remark(GL): we need a loop here with several jorek hard coded parameters: the solution 
#Remark(GL): proposed is to fill the compile_jorek function into the setting.sh script
if [ "$compile" == "yes" ]; then
  cd $codedir
  compile_jorek
fi


# --- Create run directory and copy files there
mkdir -p $tmpdir
cd $tmpdir || exit 1
#Remark(GL): TODO, we need several executables sometime, copy all of them with the regexp
#Remark(GL): "jorek_model${jorekmodel}_*" 
cp $codedir/jorek_model${jorekmodel}* $codedir/jorek_extract_data $requiredfiles . || exit 1


# --- Run the test case
cd $tmpdir || exit 1
#Remark(GL): Sequence of mpirun calls is now handled in the run_jorek function.
#Remark(GL): The run_jorek function is defined in setting.sh
#Remark(GL): The command to launch job is not always mpirun, should be flexible enough, 
#Remark(GL): so the env. variable $MPIRUN is used instead of mpirun.
run_jorek

#Remark(GL): The (extraction+comparison) process is encapuslated into a bash function written in setting.sh
#Remark(GL): I think about using h5diff also to compare the results which is a good way to go in many cases.
# --- Extract data
cd $tmpdir || exit 1
compare_jorek_res
returncode=$?

# --- Remove the temporary directory
if [ ! "$keep" == "yes" ]; then
  rm -rf $tmpdir
fi

# --- Output the return code
if [ $returncode -eq 0 ]; then
  echo "Test '$testcase' passed."
else
  echo "Test '$testcase' failed."
fi
echo ""
exit $returncode
