!> find_light_sources_contributing_to_image identifies and 
!> writes in a file all the point light sources contributing
!> to a given image. The linght intensity emitted towards the 
!> camera is also stored for generating colormaps
program find_light_sources_contributing_to_image
use constants,                                        only: PI
use mod_mpi_tools,                                    only: init_mpi_threads,finalize_mpi_threads
use particle_tracer
use mod_particle_io,                                  only: read_simulation_hdf5,get_simulation_hdf5_time
use mod_spectra_deterministic,                        only: spectrum_integrator_2nd
use mod_pinhole_lens,                                 only: pinhole_lens
use mod_filter_unity,                                 only: filter_unity
use mod_camera_perspective_static,                    only: camera_perspective_static
use mod_gyroaverage_synchrotron_lights_dist_vertices, only: gyroaverage_synchrotron_light_dist

implicit none

!> Variables -----------------------------------------------------------------------------------------------
type(event)                                 :: field_reader
type(pinhole_lens)                          :: lens
type(spectrum_integrator_2nd)               :: spectra
type(camera_perspective_static)             :: camera
type(gyroaverage_synchrotron_light_dist)    :: synch_sources
type(particle_sim),dimension(:),allocatable :: sims,sims_gc
integer                                     :: ii,t0,t1
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
n_x                 = 3 !< size of the spatial coordinates array
n_groups            = 1 !< number of particle groups
n_spectra           = 1 !< number of spectra
n_wavelengths       = 40 !< number of wavelengths per spectrum
n_int_camera_param  = 5
n_real_camera_param = 9
n_times             = 1

!> Initialisation ------------------------------------------------------------------------------------------
!> Find active light sources -------------------------------------------------------------------------------
!> Write active light sources ------------------------------------------------------------------------------
!> Finalisation --------------------------------------------------------------------------------------------

end program find_light_sources_contributing_to_image
