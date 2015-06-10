#!/bin/bash

# This script carries out a non regression test.
# * run_test.sh -h prints usage information.
# * It compiles the code.
# * It runs the test case.
# * It compares the HDF5 results data to the reference data of the test case.
# * see Wiki pages at:
#       http://jorek.eu/wiki/doku.php?id=nrt
# * The return code will be zero for a successful test, otherwise non-zero.

# --- define some Colors
NO_COL="\x1b[0m"         # black
ERROR_COL="\x1b[31;01m"  # bold red
ERROR_COL="\033[0;31m"   # red
OK_COL="\x1b[32;02m"     # green#

startdir=`readlink -f $(dirname $0)`
codedir=`readlink -f ${startdir}/..` # Assumption about source code location

if [ -z "$PRERUN" ]; then
    export PRERUN=""
fi
if [ -z "$MPIRUN" ]; then
    export MPIRUN="mpirun -n"
fi

# --- Test directory if 'non_regression_tests' exists
if [ ! -d "non_regression_tests" ]; then
    printf "$ERROR_COL Run the script from the trunk, non_regression_tests directory should exist !! \n $NO_COL"
    exit 1
fi

# --- Verify that $MPIRUN can be executed
MPIRUN_cmd=`echo $MPIRUN | cut -d' ' -f1`
stringarray=($MPIRUN)
MPIRUN_cmd=${stringarray[0]}
which $MPIRUN_cmd 1>/dev/null || exit 1

# --- Verify that 'h5diff' can be executed
which h5diff 1>/dev/null || exit 1

# --- Verify that 'Makefile.inc' exist
if [ ! -f "Makefile.inc" ]; then
    printf "$ERROR_COL Please provide a Makefile.inc file !! \n $NO_COL"
    exit 1
fi

# --- Usage printing function
function printusage() {
    echo ""
    echo " Usage:"
    echo "   $0 [options] testcase"
    echo ""
    echo " Options:"
    echo "   -h            Print this help information"
    echo "   -k            Keep temporary run directory"
    echo "   -i            Launch the inital, full-length run (not a short run starting from a restart file)"
    echo "   -l            List available test cases."
    echo "   -n            Do not compile (assume executables already exist)"
    echo "   -p            Prepare the case but do not run it"
    echo "   -t tempdir    Specify a temp directory used for the test run"
    echo "                 (default: name chosen randomly)"
    echo ""
}

# --- Process command line options
testcase="NONE" # (preset) 
compile="yes"   # (preset)
keep="no"       # (preset)
runit="yes"     # (preset)
initialrun="no" # (preset)
tmpdir="$startdir/tmp$$"

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
	cases=`ls -1 -d ${startdir}/testcases/*/`
	for i in $cases; do
	    case=$(basename $i)
	    echo ""
	    echo "*** $case ***"
	    source ${startdir}/testcases/$case/settings.sh
	    echo "$description"
	done
	echo ""
	exit 1
    elif [ "$option" == "-i" ]; then
	initialrun="yes"
	shift
    elif [ "$option" == "-p" ]; then
	runit="no"
	shift
    elif [ "$option" == "-n" ]; then
	compile="no"
	shift
    elif [ "$option" == "-t" ]; then
	tmpdir="$2"
	echo " tmpdir = " $tmpdir
	shift 2
    elif [ "${option:0:1}" != "-" ]; then
	if [ "${testcase}" != "NONE" ]; then
	    echo ""
	    printf "$ERROR_COL ERROR: Two test case names not supported at present. \n $NO_COL"
	    printusage
	    exit 1
	else
	    testcase="$1"
	fi
	shift
    else
	printf "$ERROR_COL ERROR: option $option NOT valid ! \n $NO_COL"
	printusage
	exit 1
    fi
done
echo " tmpdir = " $tmpdir

# --- Check if the testcase really exists
if [ ! -d  "${startdir}/testcases/$testcase" ]; then
  echo ""
  printf "$ERROR_COL ERROR: Test case '$testcase' does not exist.\n $NO_COL"
  printf " Valid test cases are:\n"
  cases="`ls -1 ${startdir}/testcases | tr '\n' ' ' | sed -e 's/  */ /g'`"
  echo "  $cases"
  echo ""
  printusage
  exit 1
fi
testcasedir=`readlink -f ${startdir}/testcases/$testcase`


# --- Read test case information
source $testcasedir/settings.sh


# --- Set hard-coded parameters and compile
if [ "$compile" == "yes" ]; then
  cd $codedir
  compile_jorek || exit 1
fi


# --- Create run directory and copy files there
returncode=0
if [ "$runit" == "yes" ]; then
  # --- Copy files
  mkdir -p $tmpdir
  cd $tmpdir || exit 1
  echo " Copied files " $requiredfiles
  cp $requiredfiles . || exit 1
    
  # --- Run the test case
  cd $tmpdir || exit 1
  if [ "$initialrun" == "no" ]; then
    restart_run || exit 1
  else
    initial_run || exit 1
  fi

  
  # --- Extract data
  cd $tmpdir || exit 1
  compare_results
  returncode=$?
  if [ $returncode -eq 0 ]; then
    echo "Test '$testcase' passed."
  else
    echo "Test '$testcase' failed."
  fi

  # --- Remove the temporary directory
  if [ ! "$keep" == "yes" ]; then
      rm -rf $tmpdir
  fi
fi


exit $returncode
