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
cd $codedir || exit 1

# --- Usage printing function
function printusage() {
    echo ""
    echo " Usage:"
    echo "   $0 [options] testcase"
    echo ""
    echo " Options:"
    echo "   -i            Launch the inital, full-length run (not only the test itself)"
    echo "                 NOTE: IF '-i' IS PRESENT, IT MUST BE THE FIRST OPTION"
    echo "   -h            Print this help information"
    echo "   -k            Keep temporary run directory"
    echo "   -j nthreads   Set the number of compile threads (default 1)"
    echo "   -d            Compilation with debugging options (DEBUG=1)"
    echo "   -l            List available test cases using long format."
    echo "   -L            List available test cases without any description (short format)"
    echo "   -n            Do not compile (assume executables already exist)"
    echo "   --diff        Print difference between results and reference (using Python + numpy)"
    echo "   -p            Prepare the case but do not run it"
    echo "   -t tempdir    Specify a temp directory used for the test run"
    echo "                 (default: random name in current directory)"
    echo ""
}

# --- Generic function for comparing test results (values of end.h5 versus jorek_restart.h5).
#     This function may be used for all testcases, but it is also possible to define specific
#     functions in the testcase settings.sh files.
function compare_results_generic() {
  threshold=$1
  ./rst_bin2hdf5 < ./input                                                           || exit 1
  ln -s ${testcasedir}/end.h5 end.h5                                                 || exit 1
  if [ "$printdiff" == "yes" ]; then
    echo "Difference of 'values' between result and reference: `python $startdir/tools/maximum-difference.py` (threshold: $threshold)"
  fi
  h5diff -d $threshold jorek_restart.h5 end.h5 values                                || exit 1
}

if [ -z "$PRERUN" ]; then
    export PRERUN=""
fi
if [ -z "$MPIRUN" ]; then
    export MPIRUN="mpirun -n"
fi

# --- Test if directory 'non_regression_tests' exists
if [ ! -d "${codedir}/non_regression_tests" ]; then
    printf "\n$ERROR_COL ERROR: Run the script from the trunk. \n $NO_COL"
    printusage
    exit 1
fi

# --- Verify that $MPIRUN can be executed
MPIRUN_cmd=`echo $MPIRUN | cut -d' ' -f1`
stringarray=($MPIRUN)
MPIRUN_cmd=${stringarray[0]}
which $MPIRUN_cmd >/dev/null 2>&1
if [ $? -ne 0 ]; then
  printf "\nERROR: $MPIRUN_cmd not found\n"
  exit 1
fi

# --- Verify that 'h5diff' can be executed
which h5diff >/dev/null 2>&1
if [ $? -ne 0 ]; then
  printf "\nERROR: h5diff not found\n"
  exit 1
fi

# --- Verify that 'Makefile.inc' exist
if [ ! -f "Makefile.inc" ]; then
    printf "\n$ERROR_COL Please provide a Makefile.inc file.\n $NO_COL"
    exit 1
fi

# --- Process command line options
testcase="NONE"         # (preset) 
compile="yes"           # (preset)
keep="no"               # (preset)
runit="yes"             # (preset)
initialrun="no"         # (preset)
printdiff="no"          # (preset)
debugoptions=""         # (preset)
if [ -z "$compilethreads" ]; then
    compilethreads="8"  # (preset)
fi
tmpdir="$startdir/tmp$$"

firstoption="yes"
while [ $# -gt 0 ]; do
    option="$1"
    if [ "$option" == "-h" ]; then
	printusage
	exit 1
    elif [ "$option" == "-j" ]; then
	compilethreads="$2"
	shift 2
    elif [ "$option" == "-d" ]; then
	debugoptions="DEBUG=1"
	shift
    elif [ "$option" == "-k" ]; then
	keep="yes"
	shift
    elif [ "$option" == "--diff" ]; then
        printdiff="yes"
        shift
    elif [ "$option" == "-l" ]; then
	echo ""
	echo "Available test cases:"
        echo ""
	cases=`ls -1 -d ${startdir}/testcases/*/ `
	for i in $cases; do
	    if [ -e ${i}/settings.sh ]; then
  	      case=$(basename $i)
	      source ${startdir}/testcases/$case/settings.sh
	      printf "$OK_COL %-45s $NO_COL%s\n" "$case" "$description"
              echo ""
            fi
	done
	echo ""
	exit 1
    elif [ "$option" == "-L" ]; then
	cases=`ls -1 -d ${startdir}/testcases/*/ `
	for i in $cases; do
	    if [ -e ${i}/settings.sh ]; then
              case=$(basename $i)
              echo $case
            fi
	done
	exit 0
    elif [ "$option" == "-i" ]; then
        if [ "$firstoption" == "no" ]; then
          printf "$ERROR_COL ERROR: When providing the option '-i', it needs to be the first option. \n $NO_COL"
          printusage
          exit -1
        fi
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
    firstoption="no"
done
echo " tmpdir = " $tmpdir

# --- Check if the testcase exists
if [ ! -d  "${startdir}/testcases/$testcase" ]; then
  printf "\n$ERROR_COL ERROR: Testcase '$testcase' does not exist. Use command line option -l to list available test cases.$NO_COL\n"
  printusage
  exit 1
fi
testcasedir=`readlink -f ${startdir}/testcases/$testcase`


# --- Read test case information
source $testcasedir/settings.sh


# --- Set hard-coded parameters and compile
if [ "$compile" == "yes" ]; then
  cd $codedir
  compilopt="-j $compilethreads"
  compile_jorek
  make cleanall
  if [ $? -ne 0 ]; then
    printf "\n$ERROR_COL ERROR: Compilation failed.$NO_COL\n"
    exit 1
  fi
  if [ "$initialrun" == "yes" ] && [ "$binaries_initial" != "" ]; then
    mv $binaries_initial $testcasedir/ || exit 1
  fi
  mv $binaries $testcasedir/ || exit 1
fi


# --- Create run directory and copy files there
returncode=0
if [ "$runit" == "yes" ]; then
  mkdir -p $tmpdir
  # --- Copy files
  cd $testcasedir
  echo " requiredfiles=" $requiredfiles
  cp $requiredfiles $tmpdir || exit 1
  cp $binaries $tmpdir || exit 1
  if [ "$initialrun" == "yes" ] && [ "$binaries_initial" != "" ]; then
    cp $binaries_initial $tmpdir || exit 1
  fi
  cd $tmpdir
    
  # --- Some preparations
  if [ -n "$PRERUN" ]; then
    eval $PRERUN                                          || exit 1
  fi
  if [ -n "$ompthreads" ]; then
    export OMP_NUM_THREADS=$ompthreads
  fi
  cd $tmpdir                                              || exit 1

  # --- Run the test case
  if [ "$initialrun" == "no" ]; then
    cp ${testcasedir}/begin.h5 jorek_restart.h5           || exit 1
    ./rst_hdf52bin < ./input                              || exit 1
    restart_run                                           || exit 1
    
    cd $tmpdir                                              || exit 1
    compare_results
    returncode=$?
    if [ $returncode -eq 0 ]; then
      echo "Test '$testcase' passed."
    else
      echo "Test '$testcase' failed."
    fi
  else
    initial_run                                           || exit 1
    ./rst_bin2hdf5 < ./input                              || exit 1
    cp jorek_restart.h5 ${testcasedir}/begin.h5           || exit 1
    
    sleep 3s # to avoid strange "tee: write error" problems
    
    restart_run                                           || exit 1
    ./rst_bin2hdf5 < ./input                              || exit 1
    cp jorek_restart.h5 ${testcasedir}/end.h5             || exit 1
  fi

  # --- Remove the temporary directory
  if [ ! "$keep" == "yes" ]; then
      rm -rf $tmpdir
  fi
fi

exit $returncode
