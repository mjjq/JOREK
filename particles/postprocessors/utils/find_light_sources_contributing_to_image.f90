!> find_light_sources_contributing_to_image identifies and 
!> writes in a file all the point light sources contributing
!> to a given image. The linght intensity emitted towards the 
!> camera is also stored for generating colormaps
program find_light_sources_contributing_to_image
use constants,                                       only: PI
use mod_mpi_tools,                                   only: init_mpi_threads,finalize_mpi_threads
use particle_tracer
use mod_particle_io,                                 only: read_simulation_hdf5,get_simulation_hdf5_time
use mod_spectra_deterministic,                       only: spectrum_integrator_2nd
use mod_filter_unity,                                only: filter_unity
use mod_pinhole_lens,                                only: pinhole_lens
use mod_filter_unity,                                only: filter_unity
use mod_camera_perspective_static,                   only: camera_perspective_static
use mod_gyroaverage_synchrotron_light_dist_vertices, only: gyroaverage_synchrotron_light_dist

implicit none

!> Variables -----------------------------------------------------------------------------------------------
type(event)                                 :: field_reader
type(pinhole_lens)                          :: lens
type(spectrum_integrator_2nd)               :: spectra
type(camera_perspective_static)             :: camera
type(gyroaverage_synchrotron_light_dist)    :: synch_sources
type(filter_unity),dimension(:),allocatable :: filter_spectra
type(particle_sim),dimension(:),allocatable :: sims
integer                                     :: ii,t0,t1,n_1d
integer                                     :: n_groups,my_id,n_cpus,n_x,ierr
integer                                     :: n_times,n_wavelengths,n_spectra
integer                                     :: n_int_camera_param,n_real_camera_param
integer,dimension(:),allocatable            :: int_camera_param
real*8,dimension(:),allocatable             :: min_spectra,max_spectra,pinhole_positions
real*8,dimension(:),allocatable             :: real_camera_param,sim_times
real*8,dimension(:,:),allocatable           :: active_light_source_positions
real*8,dimension(:,:),allocatable           :: active_light_source_intensities
character(len=15)                           :: particle_filename
character(len=17)                           :: fields_filename
character(len=26)                           :: light_intensity_filename
!> Variables definitions -----------------------------------------------------------------------------------
particle_filename        = 'part_restart.h5'
fields_filename          = 'jorek_equilibrium'
light_intensity_filename = 'active_light_intesities.h5'
n_1d                = 1 !< size one dimension
n_x                 = 3 !< size of the spatial coordinates array
n_groups            = 1 !< number of particle groups
n_spectra           = 1 !< number of spectra
n_wavelengths       = 40 !< number of wavelengths per spectrum
n_int_camera_param  = 5
n_real_camera_param = 9
n_times             = 1
!> JET KLDT-E5WC wavlength: 3d-6 , 3.5d-6 [m]
allocate(min_spectra(n_spectra)); min_spectra = [3d-6];
allocate(max_spectra(n_spectra)); max_spectra = [3.5d-6];
allocate(filter_spectra(n_spectra));
allocate(pinhole_positions(n_x)); pinhole_positions = [-8.86d-1,-4.002,-3.32d-1];
!> one pinhole => n_lens_samples=1, JET KLDT-E5WC pixels nx=120,ny=176
allocate(int_camera_param(n_int_camera_param)); int_camera_param = [1,120,176,0,1]
!> JET KLDT-E5WC camera inputs:
!> 1:3 -> half width, half height and orientation of the visual plane angle
!> 4:6 -> camera focal direction: focal distance, latitude, azimuth
!> 7:9 -> camera position in the tokamak reference system
allocate(real_camera_param(n_real_camera_param));
real_camera_param = [5.23d-1,5.23d-1,5d-1*PI,9.998025d-1,1.5807985,2.09801,-8.86d-1,-4.002,-3.32d-1]
allocate(sim_times(n_times)); allocate(sims(n_times));
!> Initialisation ------------------------------------------------------------------------------------------
!> initialise the MPI communicator
call init_mpi_threads(my_id,n_cpus,ierr)

!> read particle, time and MHD fields
write(*,*) "Reading particle and MHD data ..."
field_reader = event(read_jorek_fields_interp_linear(basename=trim(fields_filename),i=-1))
do ii=1,n_times
  call sims(ii)%initialize(n_groups,.true.,my_id,n_cpus)
  call read_simulation_hdf5(sims(ii),particle_filename)
  sim_times(ii) = get_simulation_hdf5_time(particle_filename)
  call with(sims(ii),field_reader)
enddo
write(*,*) "Reading particle and MHD data: completed!"

!> Initialise synthetic camera and light sources
write(*,*) 'Initialise synthetic camera and light sources ...'
call system_clock(t0)
spectra = spectrum_integrator_2nd(n_wavelengths,n_spectra,min_spectra,max_spectra)
call spectra%generate_spectrum()
do ii=1,n_spectra
  call filter_spectra(ii)%init_filter(n_1d)
enddo
call lens%init_pinhole(n_x,pinhole_positions)
call camera%init_camera(lens,spectra,n_int_camera_param,n_real_camera_param,&
int_camera_param,real_camera_param)
call synch_sources%init_lights_from_particles(n_times,sims)
call system_clock(t1)
write(*,*) 'Initialise synthetic camera and light sources: completed!'
write(*,*) my_id, 'System time fast camera initialisation (s): ',real(t1-t0,kind=8)/1d3

!> Find active light sources -------------------------------------------------------------------------------
!> Write active light sources ------------------------------------------------------------------------------
!> Finalisation --------------------------------------------------------------------------------------------
if(allocated(min_spectra))       deallocate(min_spectra);
if(allocated(max_spectra))       deallocate(max_spectra);
if(allocated(filter_spectra))    deallocate(filter_spectra);
if(allocated(pinhole_positions)) deallocate(pinhole_positions);
if(allocated(int_camera_param))  deallocate(int_camera_param);
if(allocated(real_camera_param)) deallocate(real_camera_param);
if(allocated(sim_times))         deallocate(sim_times);
if(allocated(sims))              deallocate(sims);
call finalize_mpi_threads(ierr)
end program find_light_sources_contributing_to_image
