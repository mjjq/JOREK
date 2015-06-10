module purge
module load intel/14.0.4.211
module load bullxmpi/1.2.8.2
module load hdf5/1.8.14

export JOREK_HOST=helios
export PRERUN='export OMP_NUM_THREADS=16'
export MPIRUN='mpirun -np '


