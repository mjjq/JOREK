description="Test case for time evolution of model 199: Tearing mode in circular plasma at n_tor=3."

# --- Hard coded parameters
jorekmodel="199"
jorekparameters="n_tor=3 n_plane=4 n_period=1"

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
