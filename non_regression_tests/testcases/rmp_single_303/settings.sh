# --- General settings
jorekmodel="303"
description="after eq with flows single rmp n_tor=3."
mpitasks=2
binaries="jorek_model${jorekmodel}_1 jorek_model${jorekmodel}_3 rst_bin2hdf5 rst_hdf52bin"
requiredfiles="$binaries JET_EQ_FLOWS_303_NTOR1_NPER1_NPLA12 JET_RMP3_303_NTOR3_NPER3_NPLA32 RMP_psi_cos_JET_N3.txt RMP_psi_sin_JET_N3.txt RMP_start_time.dat "


# --- Compile the code for the test case
function compile_jorek () {
  ./util/config.sh model=$jorekmodel n_tor=1 n_plane=1 n_period=12                    || exit 1
  make cleanall                                                                      || exit 1
  make $compilopt jorek_model${jorekmodel}                                           || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                             || exit 1
  ./util/config.sh model=$jorekmodel n_tor=3 n_plane=32 n_period=3                    || exit 1
  make cleanall                                                                      || exit 1
  make $compilopt jorek_model${jorekmodel} rst_bin2hdf5 rst_hdf52bin                 || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_3                             || exit 1
}


# --- Re-run the whole case from scratch into the non-linear phase
function initial_run () {
  #${codedir}/util/setinput.sh input nstep_n=10,10,10, tstep_n=1.,100.,3000.          || exit 1
  ./jorek_model${jorekmodel}_1 < JET_EQ_FLOWS_303_NTOR1_NPER1_NPLA12 | tee -a logfile                                    || exit 1
 # ${codedir}/util/setinput.sh input nstep_n=30 tstep_n=3000. restart=.t.             || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < JET_RMP3_303_NTOR3_NPER3_NPLA32  | tee -a logfile                  || exit 1
}


# --- Carry out the test case, i.e., run a single time step in the non-linear phase
function restart_run () {
  ${codedir}/util/setinput.sh JET_RMP3_303_NTOR3_NPER3_NPLA32 restart=.t. nstep_n=1 nout=1       || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < JET_RMP3_303_NTOR3_NPER3_NPLA32 | tee -a logfile                  || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  ./rst_bin2hdf5 < ./JET_RMP3_303_NTOR3_NPER3_NPLA32                                                           || exit 1
  h5diff -d 1e-12 jorek_restart.h5 ${testcasedir}/end.h5 values                      || exit 1
}
