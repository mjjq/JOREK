#!/bin/bash

description="Test case for time evolution of model 303: Ballooning mode in X-point plasma at n_tor=3."

# --- Model used
jorekmodel="303"

# --- Files required to run the code (executables copied automatically)
requiredfiles="$testcasedir/input $testcasedir/jorek_restart.rst"

# --- How many MPI tasks and OpenMP threads are required?
mpitasks=2
ompthreads=4

# --- Compare which kind of data? (separate keywords with blanks)
comparedata="energies"

# --- Tolerance for the comparison. A comparison fails if both thresholds are exceeded.
relativeaccuracy="1.0e-10"
absoluteaccuracy="0.0e+00" # i.e., check only for relative differences

function compile_jorek () {
  ./util/config.sh model=$jorekmodel "n_tor=3 n_plane=4 n_period=6"
  make clean
  make -j 2 jorek_model${jorekmodel} && make -j 2 jorek_extract_data 
}

function run_jorek () {
  $PRERUN
  export OMP_NUM_THREADS=$ompthreads
  $MPIRUN $mpitasks ./jorek_model$jorekmodel < input > logfile || exit 1
}

function compare_jorek_res () {
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
  return $returncode
}

