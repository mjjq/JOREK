
module purge
module load intel/15.0.0.090 bullxmpi/1.2.8.3 hdf5/1.8.14
export LANG=C
export HDF5HOME=/opt/software/libraries/hdf5/1.8.14
export JOREK_HOST=occigen
export BATCHCOMMAND="sbatch"
export KMP_STACK_SIZE=16M
export KMP_AFFINITY="granularity=fine,compact,1,0,verbose"
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=0
export MPIRUN="srun --mpi=pmi2 -K1 -m block:block -c 24 --resv-ports -n "
