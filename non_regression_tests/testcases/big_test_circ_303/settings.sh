# --- General settings
jorekmodel="303"
description="Tearing mode, circular plasma, model$jorekmodel, n_tor=7."
mpitasks=10
binaries="jorek_model${jorekmodel}_7 rst_bin2hdf5 rst_hdf52bin"
binaries_initial="jorek_model${jorekmodel}_1 jorek_model${jorekmodel}_3"
requiredfiles="input"
extra_remote_files=""


# --- Compile the code for the test case
function compile_jorek () {
  if [ "$initialrun" == "yes" ]; then
    ./util/config.sh model=$jorekmodel n_tor=1 n_plane=1 n_period=1                    || exit 1
    make cleanall                                                                      || exit 1
    make $compilopt jorek_model${jorekmodel}                                           || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                             || exit 1

    ./util/config.sh model=$jorekmodel n_tor=3 n_plane=32 n_period=1                    || exit 1
    make cleanall                                                                      || exit 1
    make $compilopt jorek_model${jorekmodel} rst_bin2hdf5 rst_hdf52bin                 || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_3                             || exit 1
  fi
  ./util/config.sh model=$jorekmodel n_tor=7 n_plane=32 n_period=1                    || exit 1
  make cleanall                                                                      || exit 1
  make $compilopt jorek_model${jorekmodel} rst_bin2hdf5 rst_hdf52bin                 || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_7                             || exit 1
}


# --- Re-run the whole case from scratch into the non-linear phase
function initial_run () {
  ${codedir}/util/setinput.sh input nstep_n=0  tstep_n=1.          || exit 1
  ${MPIRUN} 1 ./jorek_model${jorekmodel}_1 < input | tee -a logfile          || exit 1
  ${codedir}/util/setinput.sh input nstep_n=10,10,10 restart=.t. tstep_n=1.,100.,3000.          || exit 1
  ${MPIRUN} 1 ./jorek_model${jorekmodel}_1 < input | tee -a logfile          || exit 1
  ${codedir}/util/setinput.sh input nstep_n=30 tstep_n=3000. restart=.t.             || exit 1
  ${MPIRUN} $mpitasks ./jorek_model${jorekmodel}_3 < input | tee -a logfile          || exit 1
  ${codedir}/util/setinput.sh input nstep_n=1 tstep_n=300. restart=.t.             || exit 1
  ${MPIRUN} $mpitasks ./jorek_model${jorekmodel}_7 < input | tee -a logfile          || exit 1
}


# --- Carry out the test case, i.e., run a single time step in the non-linear phase
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=300. nout=1       || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_7 < input | tee -a logfile                  || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  ./rst_bin2hdf5 < ./input                                                           || exit 1
  h5diff -d 1e-12 jorek_restart.h5 ${testcasedir}/end.h5 values                      || exit 1
}
