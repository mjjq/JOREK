description="Test case for X-point equilibrium."

# --- Hard coded parameters
jorekmodel="199"
jorekparameters="n_tor=1 n_plane=1 n_period=1"

# --- Files required to run the code (executables copied automatically)
requiredfiles="$testcasedir/input"

# --- How many MPI tasks and OpenMP threads are required?
mpitasks=1
ompthreads=1

# --- Compare which kind of data? (separate keywords with blanks)
comparedata="special_points"

# --- Tolerance for the comparison. A comparison fails if both thresholds are exceeded.
relativeaccuracy="0.0e+00" # i.e., check only for absolute differences
absoluteaccuracy="1.0e-08"
