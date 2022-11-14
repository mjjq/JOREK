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
use mod_mpi_tools,                  only: init_mpi_threads,finalize_mpi_threads
use particle_tracer
use mod_particle_io,                only: read_simulation_hdf5,get_simulation_hdf5_time
use mod_spectra_deterministic,      only: spectrum_integrator_2nd
use mod_pinhole_lens,               only: pinhole_lens
use mod_camera_perspective_static,  only: camera_perspective_static
use mod_synchrotron_light_vertices, only: synchrotron_light_vertices

implicit none

!> Variables -------------------------------------------------------------------------------
type(event)                       :: field_reader
type(camera_perspective_static)   :: camera
type(spectrum_integrator_2nd)     :: spectra
type(pinhole_lens)                :: lens
type(synchrotron_light_vertices)  :: synch_sources
type(particle_sim),dimension(:),allocatable :: sims
integer                           :: ii
integer                           :: n_groups,my_id,n_cpus,n_x,ierr
integer                           :: n_wavelenghts,n_spectra
integer                           :: n_int_camera_param,n_real_camera_param
integer                           :: n_times
integer,dimension(:),allocatable  :: int_camera_param
real*8,dimension(:),allocatable   :: min_spectra,max_spectra,pinhole_positions
real*8,dimension(:),allocatable   :: real_camera_param,sim_times
character(len=15)                 :: particle_filename
character(len=17)                 :: fields_filename

!> Variable definitions -------------------------------------------------------------------
particle_filename = 'part_restart.h5'
fields_filename   = 'jorek_equilibrium' 
n_x = 3 !< number of spatial coordinates
n_groups = 1 !< number of particle groups
n_spectra     = 1
n_wavelenghts = 100
n_int_camera_param  = 3
n_real_camera_param = 8
n_times = 1
allocate(min_spectra(n_spectra)); min_spectra = [0d0];
allocate(max_spectra(n_spectra)); max_spectra = [1d0];
allocate(pinhole_positions(n_x)); pinhole_positions = [0d0,0d0,0d0];
allocate(int_camera_param(n_int_camera_param)); int_camera_param = [0,0,0];
allocate(real_camera_param(n_real_camera_param)); 
real_camera_param = [0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0];
allocate(sim_times(n_times)); allocate(sims(n_times));

!> Initialisation  ------------------------------------------------------------------------
!> Initialise MPI communicator
call init_mpi_threads(my_id,n_cpus,ierr)

!> read particle, time and MHD fields data
write(*,*) 'Reading particle data ...'
field_reader = event(read_jorek_fields_interp_linear(basename=trim(fields_filename),i=-1))
do ii=1,n_times
  call read_simulation_hdf5(sims(ii),particle_filename)
  sim_times(ii) = get_simulation_hdf5_time(particle_filename)
  call with(sims(ii),field_reader)
enddo
write(*,*) 'Reading particle data: completed!'

!> Initialise synthetic diagnostics
write(*,*) 'Initialise synthetic camera and light sources'
spectra = spectrum_integrator_2nd(n_wavelenghts,n_spectra,min_spectra,max_spectra)
call lens%init_pinhole(n_x,pinhole_positions)
call camera%init_camera(lens,spectra,n_int_camera_param,n_real_camera_param,&
int_camera_param,real_camera_param)
call synch_sources%init_lights_from_particles(n_times,sims)
write(*,*) 'Initialise synthetic camera and light sources: completed'

!> Finalisation ---------------------------------------------------------------------------
if(allocated(min_spectra))       deallocate(min_spectra)
if(allocated(max_spectra))       deallocate(max_spectra)
if(allocated(pinhole_positions)) deallocate(pinhole_positions)
if(allocated(int_camera_param))  deallocate(int_camera_param)
if(allocated(real_camera_param)) deallocate(real_camera_param)
if(allocated(sim_times))         deallocate(sim_times)
if(allocated(sims))              deallocate(sims)
call finalize_mpi_threads(ierr)

end program camera_RE_example
