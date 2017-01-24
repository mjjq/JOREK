# --- General settings
jorekmodel="303"
description="after eq with flows single rmp n_tor=3."
mpitasks=2
binaries="jorek_model${jorekmodel}_1 jorek_model${jorekmodel}_3 rst_bin2hdf5 rst_hdf52bin"
requiredfiles="$binaries input JET_RMP3_303_NTOR3_NPER3_NPLA32 RMP_psi_cos_JET_N3.txt RMP_psi_sin_JET_N3.txt RMP_start_time.dat "


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
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input | tee -a logfile                                    || exit 1
  echo "Equil done"
  ${codedir}/util/setinput.sh input restart=.t. tstep_n=5. nstep_n=1 iter_precon=0 gmres_4=1.d4 RMP_on=.t. RMP_psi_cos_file=RMP_psi_cos_JET_N3.txt RMP_psi_sin_file=RMP_psi_sin_JET_N3.txt nout=1       || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input  | tee -a logfile                  || exit 1
}


# --- Carry out the test case, i.e., run a single time step in the non-linear phase
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.t. tstep_n=5. nstep_n=1 iter_precon=0 gmres_4=1.d4 RMP_on=.t. RMP_psi_cos_file=RMP_psi_cos_JET_N3.txt RMP_psi_sin_file=RMP_psi_sin_JET_N3.txt nout=1       || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input | tee -a logfile                  || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  ./rst_bin2hdf5 < ./input                                                           || exit 1
  h5diff -d 1e-12 jorek_restart.h5 ${testcasedir}/end.h5 values                      || exit 1
}
