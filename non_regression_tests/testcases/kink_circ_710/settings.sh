# --- General settings
jorekmodel="710"
description="Internal kink, circular plasma, model$jorekmodel, n_tor=3."
mpitasks=2
binaries="jorek_model${jorekmodel}_3 rst_bin2hdf5 rst_hdf52bin"
requiredfiles="$binaries input"


# --- Compile the code for the test case
function compile_jorek () {
  ./util/config.sh model=$jorekmodel n_tor=3 n_plane=4 n_period=1                    || exit 1
  make cleanall                                                                      || exit 1
  make $compilopt jorek_model${jorekmodel} rst_bin2hdf5 rst_hdf52bin                 || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_3                             || exit 1
}


# --- Re-run the whole case from scratch into the non-linear phase
function initial_run () {
  ${codedir}/util/setinput.sh input nstep_n=21 tstep_n=300.                          || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile                  || exit 1
}


# --- Carry out the test case, i.e., run a single time step in the non-linear phase
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=300. nout=1        || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile                  || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  ./rst_bin2hdf5                                                                     || exit 1
  h5diff -d 1e-07 jorek_restart.h5 ${testcasedir}/end.h5 values                      || exit 1
}
