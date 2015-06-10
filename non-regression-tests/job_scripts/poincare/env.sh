module purge
module load intel/12.1.0 mvapich2/1.8.1 mkl/10.3  hwloc/1.6.2_intel
export PRERUN='export KMP_AFFINTY=verbose'
export MPIRUN="mpiexec -launcher-exec /opt/ibmll/LoadL/scheduler/full/bin/llspawn.stdio -f $LOADL_HOSTFILE  -n "
