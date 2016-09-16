title: Example
This page contains instructions on how to run a simple particle case.

### You will need
* A [compiled `JOREK2_particles` binary](installation.html)
* One or more JOREK restart files
* A JOREK input file named `in_jorek`, used to generate the above restart files
* A JOREK2_particles input file.
* ADAS ADF11-coefficients. Get them [here](http://open.adas.ac.uk/adf11?element=w&metastable_unresolved=1&acd=1&scd=1&plt=1&prb=1&year=&searching=2#searchbutton)

### Preparing the simulation
Your directory will look something like this:
```bash
acd50_w.dat  in_jorek  in_particles  job_mpi  jorek_restart  plt50_w.dat  prb50_w.dat  scd50_w.dat
```

Your input file can be something like
```ini
&in2
  species(1) = 74
  atomic_mass(1) = 183.84
  N_particles(1) = 10000
  adas_suffix(1) = '50_w'

!  location_accept_function = 'joined_gaussian'
!  location_accept_parameters(1,1) = 1 ! Use psi
!  location_accept_parameters(2,1) = -1.00 ! Center of first
!  location_accept_parameters(3,1) = 0.001 ! stdev of first
!  location_accept_parameters(4,1) = -0.47 ! center of second
!  location_accept_parameters(5,1) = 0.001 ! stdev of second

  t_particles_begin = -1
  t_step_particles  = 0.0015
  n_step_particles  = 5000
  nout_particles    = 1000
&end
```

### Run the simulation
To run, either submit a job using your queueing system or run locally:
```bash
jorek2_particles < in_particles
```
which will use all of your cores in OpenMP mode by default.

#### Example jobscript
```bash
#!/bin/sh
#PBS -N part_test
#PBS -V
#PBS -j oe
#PBS -l walltime="02:59:00"
#PBS -l nodes=4:ppn=8

cd $PBS_O_WORKDIR

rm -f ${PBS_JOBNAME}.{out,err}
export I_MPI_EXTRA_FILESYSTEM=enable
export I_MPI_EXTRA_FILESYSTEM_LIST=lustre

export OMP_PROC_BIND=true
mpirun -ppn 1 ~/jorek/jorek2_particles < in_particles 2>${PBS_JOBNAME}.err >${PBS_JOBNAME}.out
```
