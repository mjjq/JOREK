!> The camera_example.f90 contains an executable programme which
!> tests and shows how to run the synthetic camera diagnostic
!> for generating images of plasma radiation from a JORE-particle
!> population. The test proposed below concerns the imaging of the
!> sychrotron radiation of a runaway electron population.
!> More specifically, we consider the JET tokamak configuration
!> having an IR and fast visible camera looking at the RE beam
!> from octant 5. Note that the RE toroidal motion
!> is counter-clockeise. Only one snapshot of the RE beam is 
!> considered hereater hence, only one particle population and 
!> JOREK MHD fields are used.
program camera_RE_example
use mod_mpi_tools,   only: init_mpi_threads,finalize_mpi_threads
use particle_tracer
use mod_particle_io, only: read_simulation_hdf5,get_simulation_hdf5_time
implicit none
!> Variables
type(event)     :: field_reader
integer         :: n_groups,my_id,n_cpus,ierr
real*8          :: sim_time
character(len=15) :: particle_filename
character(len=17) :: fields_filename

!> Variable definitions -------------------------------------------------------------------
particle_filename = 'part_restart.h5'
fields_filename   = 'jorek_equilibrium' 
n_groups = 1 !< number of particle groups

!> Initialisation  ------------------------------------------------------------------------
!> Initialise MPI communicator and simulation
call init_mpi_threads(my_id,n_cpus,ierr)
call sim%initialize(n_groups,.true.,my_id,n_cpus)

!> read particle, time and MHD fields data
write(*,*) 'Reading particle data ...'
call read_simulation_hdf5(sim,particle_filename)
write(*,*) 'Reading particle data completed!'
write(*,*) 'Read particle simulation time ...'
sim_time = get_simulation_hdf5_time(particle_filename)
write(*,*) 'Reading particle simulation time completed! Simulation time is: ',sim_time
field_reader = event(read_jorek_fields_interp_linear(basename=trim(fields_filename),i=-1))
call with(sim,field_reader)

!> Finalisation ---------------------------------------------------------------------------
call finalize_mpi_threads(ierr)

end program camera_RE_example
