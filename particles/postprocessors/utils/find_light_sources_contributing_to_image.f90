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
use mod_pinhole_lens,                                only: pinhole_lens
use mod_camera_perspective_static,                   only: camera_perspective_static
use mod_gyroaverage_synchrotron_light_dist_vertices, only: gyroaverage_synchrotron_light_dist

implicit none

!> Variables -----------------------------------------------------------------------------------------------
type(event)                                 :: field_reader
type(pinhole_lens)                          :: lens
type(spectrum_integrator_2nd)               :: spectra
type(camera_perspective_static)             :: camera
type(gyroaverage_synchrotron_light_dist)    :: synch_sources
type(particle_sim),dimension(:),allocatable :: sims
integer                                     :: ii,t0,t1
integer                                     :: n_groups,my_id,n_cpus,n_x,ierr
integer                                     :: n_times,n_wavelengths,n_spectra
integer                                     :: n_int_camera_param,n_real_camera_param
integer,dimension(:),allocatable            :: int_camera_param
real*8,dimension(:),allocatable             :: min_spectra,max_spectra,pinhole_positions
real*8,dimension(:),allocatable             :: real_camera_param,sim_times
real*8,dimension(:,:,:),allocatable         :: active_light_source_positions
real*8,dimension(:,:,:),allocatable         :: active_light_source_intensities
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
!> JET KLDT-E5WC wavlength: 3d-6 , 3.5d-6 [m]
allocate(min_spectra(n_spectra)); min_spectra = [3d-6];
allocate(max_spectra(n_spectra)); max_spectra = [3.5d-6];
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
call lens%init_pinhole(n_x,pinhole_positions)
call camera%init_camera(lens,spectra,n_int_camera_param,n_real_camera_param,&
int_camera_param,real_camera_param)
call synch_sources%init_lights_from_particles(n_times,sims)
call system_clock(t1)
allocate(active_light_source_positions(n_x,sum(synch_sources%n_active_vertices),camera%n_vertices)) 
active_light_source_positions = 0d0;
allocate(active_light_source_intensities(n_spectra,sum(synch_sources%n_active_vertices),camera%n_vertices)); 
active_light_source_intensities = 0d0;
write(*,*) 'Initialise synthetic camera and light sources: completed!'
write(*,*) my_id, 'System time fast camera initialisation (s): ',real(t1-t0,kind=8)/1d3

!> Find active light sources -------------------------------------------------------------------------------
write(*,*) "Findig light sources contributing to an image and computing spectral intensities: ..."
call system_clock(t0)
call compute_contribution_light_source_to_image_static(camera,synch_sources,spectra,&
active_light_source_positions,active_light_source_intensities)
call system_clock(t1)
write(*,*) "Findig light sources contributing to an image and computing spectral intensities: completed!"
write(*,*) my_id,'System time finding contributing light sources and computing intensities (s): ',real(t1-t0,kind=8)

!> Write active light sources ------------------------------------------------------------------------------

!> Finalisation --------------------------------------------------------------------------------------------
if(allocated(min_spectra))       deallocate(min_spectra);
if(allocated(max_spectra))       deallocate(max_spectra);
if(allocated(pinhole_positions)) deallocate(pinhole_positions);
if(allocated(int_camera_param))  deallocate(int_camera_param);
if(allocated(real_camera_param)) deallocate(real_camera_param);
if(allocated(sim_times))         deallocate(sim_times);
if(allocated(sims))              deallocate(sims);
call finalize_mpi_threads(ierr)

contains

!> Tools ---------------------------------------------------------------------------------------------------
!> methods used for computing the contribution of a light source
!> to an image (for each point on a lens and for all time)
!> inputs:
!>   camera_inout:  (camera_static)  class defining a static camera
!>   lights_inout:  (light_vertices) class containing all active lights
!>   spectra_inout: (spectral_base)  camera spectrum class
!> outputs:
!>   active_light_source_positions:   (real8)(n_x,n_active_lights,n_lens_points) position
!>                                    of the light sources contributing to an image
!>   active_light_source_intensities: (real8)(n_spectra,n_contributing_lights,n_lens_points)
!>                                    integrated spectral intensities of the light sources
!>                                    contributing to an image
subroutine compute_contribution_light_source_to_image_static(camera_inout,lights_inout,&
spectra_inout,active_light_source_positions,active_light_source_intensities)
  use mod_spectra,        only: spectrum_base
  use mod_light_vertices, only: light_vertices
  use mod_camera_static,  only: camera_static
  implicit none
  !> inputs-outputs:
  class(camera_static),intent(inout)  :: camera_inout
  class(light_vertices),intent(inout) :: lights_inout
  class(spectrum_base),intent(inout)  :: spectra_inout
  !> outputs:
  real*8,dimension(:,:,:),allocatable,intent(out) :: active_light_source_positions
  real*8,dimension(:,:,:),allocatable,intent(out) :: active_light_source_intensities
  !> variables:3
  logical :: intersect
  integer :: ii,jj,kk,counter
  real*8  :: material_value,visibility_geometry
  real*8,dimension(3)                                              :: plane_line_coords
  real*8,dimension(spectra_inout%n_spectra)                        :: integrated_irradiance
  real*8,dimension(spectra_inout%n_points,spectra_inout%n_spectra) :: spectral_irradiance
  !> initialisation
  counter = 1
  do jj=1,lights_inout%n_times !< loop on the times
    !$omp parallel do default(private) firstprivate(jj) &
    !$omp shared(lights_inout,camera_inout,spectra_inout,&
    !$omp active_light_source_positions,active_light_source_intensities) &
    !$omp reduction(+:counter) collapse(2)
    !> compute the contribution of each light for each time
    do ii=1,camera_inout%n_vertices                !< loop on the camera lens points
      do kk=1,lights_inout%n_active_vertices(jj) !< loop on the lights
        !> compute the material function skip is <= 0
        call camera_inout%physical_material_funct(lights_inout%x(:,kk,jj),ii,material_value)
        if(material_value.le.0.d0) cycle
        !> compute the intersection camera - light source, skip if not
        call camera_inout%find_ray_image_plane_intersection(lights_inout%x(:,kk,jj),ii,&
        intersect,plane_line_coords); if(.not.intersect) cycle;
        !> compute the geometry and visibility functions
        call camera_inout%visibility_geometry_funct(lights_inout,1,jj,ii,kk,visibility_geometry)
        !> compute the spectral irradiance
        call lights_inout%spectral_irradiance(spectra_inout,jj,kk,camera_inout%x(:,ii,1),&
        spectral_irradiance)
        call spectra_inout%integrate_data(spectral_irradiance,integrated_irradiance)
        !> store the light position and light contribution
        active_light_source_positions(:,counter,ii) = lights_inout%x(:,kk,jj)
        active_light_source_intensities(:,counter,ii) = integrated_irradiance*&
        material_value*visibility_geometry
        !> update counter
        counter = counter + 1
      enddo
    enddo
    !$omp end parallel do
  enddo
end subroutine compute_contribution_light_source_to_image_static
!> ---------------------------------------------------------------------------------------------------------
end program find_light_sources_contributing_to_image
