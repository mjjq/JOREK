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
testcase="NONE"
compile="yes"
keep="no"
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
relativeaccuracy="0.0"; absoluteaccuracy="0.0" # (presets)
source $testcasedir/settings.sh


# --- Set hard-coded parameters and compile
#Remark(GL): we need several executables sometime for one single run (ntor=1, ntor=XX), 
#Remark(GL): that's why util/compile_test.sh generete several ones
if [ "$compile" == "yes" ]; then
  cd $codedir
  ./util/config.sh model=$jorekmodel $jorekparameters
  make clean
  make -j 3 || exit 1
  make -j 3 jorek_extract_data || exit 1
fi


# --- Create run directory and copy files there
mkdir -p $tmpdir
cd $tmpdir || exit 1
cp $codedir/jorek_model$jorekmodel $codedir/jorek_extract_data $requiredfiles . || exit 1


# --- Run test case
cd $tmpdir || exit 1
export OMP_NUM_THREADS=$ompthreads
mpirun -n $mpitasks ./jorek_model$jorekmodel < input > logfile || exit 1


# --- Extract data
cd $tmpdir || exit 1

echo "&extract" > extract_data.nml
comparedata2=`echo "$comparedata" | sed -e 's/^ *//' -e 's/ *$//' -e 's/  */ /g' -e "s/^/'/" -e "s/$/'/" -e "s/ /', '/g"`
echo "  extract_data = $comparedata2" >> extract_data.nml
echo "/" >> extract_data.nml
./jorek_extract_data < input

cp extracted_data.dat $testcasedir/extracted_data_`date "+%Y-%m-%d_%H-%M-%S"`.dat ### TODO: Remove later on


# --- Compare the result to the reference data
# ### TODO: This has to be implemented cleaner later on...
echo "print \\" > comparison.dat
paste $testcasedir/reference.dat extracted_data.dat | sed -e 's/\t/ /' -e                                                     \
  "s|^ *\([^ ]*\) *\([^ ]*\)|( abs(\1-\2)/abs(\1+1e-99) <= $relativeaccuracy or abs(\1-\2) <= $absoluteaccuracy ) and \\\\|"  \
  >> comparison.dat || exit 1
echo "True" >> comparison.dat
okay=`python comparison.dat`
returncode=1
if [ "$okay" == "True" ]; then
  returncode=0
fi


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
