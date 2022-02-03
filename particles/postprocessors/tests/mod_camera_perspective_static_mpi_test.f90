!> mod_camera_perspective_static_mpi_test cantains all variables
!> and procedures used for testing mpi enabled procedures of the
!> camera_perspective_static class.
module mod_camera_perspective_static_mpi_test
use fruit
use fruit_mpi
use constants,                           only: PI,TWOPI
use mod_spectra_deterministic,           only: spectrum_integrator_2nd
use mod_filter_unity,                    only: filter_unity
use mod_pinhole_lens,                    only: pinhole_lens
use mod_camera_perspective_static,       only: camera_perspective_static
use mod_omnidirectional_gaussian_lights, only: omnidirectional_gaussian_lights
implicit none
private
public :: run_fruit_camera_perspective_static_mpi

!> Variable and data types ---------------------------------------------------
integer,parameter                           :: n_x_sol=3
integer,parameter                           :: n_spectra_sol=2
integer,parameter                           :: n_lines_per_spectrum_sol=53
integer,parameter                           :: n_times=3
integer,parameter                           :: n_pixels_x=64
integer,parameter                           :: n_pixels_y=32
integer,parameter                           :: n_particles_time=1345
real*8                                      :: plane_distance_sol=4.23d0
real*8,dimension(2),parameter               :: half_angle_lowbnd=(/PI/1.d1,PI/4.d0/)
real*8,dimension(2),parameter               :: half_angle_uppbnd=(/PI/1.9d0,PI/3.d0/)
real*8,dimension(2),parameter               :: costheta_interval=(/-1.d0,1.d0/)
real*8,dimension(2),parameter               :: phi_interval=(/0.d0,TWOPI/)
real*8,dimension(n_spectra_sol),parameter   :: min_wlen=(/3.d-6,2.5d-7/)
real*8,dimension(n_spectra_sol),parameter   :: max_wlen=(/3.5d-6,4.2d-7/)
real*8,dimension(n_x_sol),parameter         :: center_pos_lowbnd=(/-9.d-1,6.d1,-3.2d0/)
real*8,dimension(n_x_sol),parameter         :: center_pos_uppbnd=(/3.5d0,9.5d1,1.2d0/)
type(filter_unity)                          :: filter_pixel_sol
type(filter_unity)                          :: filter_time_sol
type(spectrum_integrator_2nd)               :: spectra_sol
type(pinhole_lens)                          :: pinhole_sol
type(camera_perspective_static)             :: camera_sol
type(omnidirectional_gaussian_lights)       :: lights_sol
type(filter_unity),dimension(n_spectra_sol) :: filter_spectra_sol

!> Interfaces ----------------------------------------------------------------
contains
!> Fruit basket --------------------------------------------------------------
!> fruit basket executes all set-up, test and tear-down procedures
subroutine run_fruit_camera_perspective_static_mpi(rank,n_tasks,ifail)
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  !> inputs:
  integer,intent(in) :: rank,n_tasks
  if(rank.eq.0) write(*,*) "  ... setting-up: camera perspective static mpi tests"
  call setup_spectra_filters(rank,n_tasks,ifail)
  call setup_camera(rank,n_tasks,ifail)
  if(rank.eq.0) write(*,*) "  ... running: camera perspective static mpi tests"
  if(rank.eq.0) write(*,*) "  ... tearing-down: camera perspective static mpi tests"
  call teardown(rank,n_tasks,ifail)
end subroutine run_fruit_camera_perspective_static_mpi

!> Set-up and tear-down ------------------------------------------------------
!> set-up spectra and filter types
subroutine setup_spectra_filters(rank,n_tasks,ifail)
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  !> inputs:
  integer,intent(in)    :: rank,n_tasks
  !> variables
  integer :: ii,n_1d,n_2d
  !> initialisations
  n_1d = 1; n_2d = 2;
  spectra_sol = spectrum_integrator_2nd(n_lines_per_spectrum_sol,&
  n_spectra_sol,min_wlen,max_wlen)
  call filter_time_sol%init_filter(n_1d)
  call filter_pixel_sol%init_filter(n_2d)
  do ii=1,n_spectra_sol
    call filter_spectra_sol(ii)%init_filter(n_1d)
  enddo
end subroutine setup_spectra_filters

!> set-up the pinhole and camera features
subroutine setup_camera(rank,n_tasks,ifail)
  use mpi
  use mod_gnu_rng,  only: gnu_rng_interval
  use mod_sampling, only: sample_uniform_sphere
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  !> inputs:
  integer,intent(in)    :: rank,n_tasks
  !> variables:
  integer :: n_int_param,n_real_param
  integer,dimension(:),allocatable :: int_param
  real*8,dimension(n_x_sol)        :: center
  real*8,dimension(:),allocatable  :: real_param
  !> initialisation
  n_int_param = 3; n_real_param = 8;
  if(.not.allocated(int_param))  allocate(int_param(n_int_param))
  if(.not.allocated(real_param)) allocate(real_param(n_real_param))
  !> initialise pinhole object
  int_param = (/n_lines_per_spectrum_sol,n_pixels_x,n_pixels_y/)
  if(rank.eq.0) call gnu_rng_interval(n_x_sol,center_pos_lowbnd,center_pos_uppbnd,center)
  call MPI_Bcast(center,n_x_sol,MPI_DOUBLE,0,MPI_COMM_WORLD,ifail)
  call pinhole_sol%init_pinhole(n_x_sol,center)
  !> initialise camera object
  int_param = (/n_lines_per_spectrum_sol,n_pixels_x,n_pixels_y/)
  if(rank.eq.0) then
    call gnu_rng_interval(2,half_angle_lowbnd,half_angle_uppbnd,real_param(1:2))
    call random_number(real_param(3:5))
    real_param(3:5) = sample_uniform_sphere(plane_distance_sol,&
    costheta_interval,phi_interval,real_param(3:5))
    call gnu_rng_interval(n_x_sol,center_pos_lowbnd,center_pos_uppbnd,real_param(6:8))
  endif
  call MPI_Bcast(real_param,n_real_param,MPI_DOUBLE,0,MPI_COMM_WORLD,ifail)
  call camera_sol%init_camera(pinhole_sol,spectra_sol,n_int_param,&
  n_real_param,int_param,real_param)
  !> cleanup
  if(allocated(int_param))  deallocate(int_param)
  if(allocated(real_param)) deallocate(real_param)
end subroutine setup_camera
!> set-up the particles features

!> tear-down all test features
subroutine teardown(rank,n_tasks,ifail)
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  !> inputs:
  integer,intent(in)    :: rank,n_tasks
  !> variables:
  integer :: ii
  call spectra_sol%deallocate_spectrum !< clean spectra
  !> clean filter
  call filter_time_sol%deallocate_filter
  call filter_pixel_sol%deallocate_filter
  do ii=1,n_spectra_sol
    call filter_spectra_sol(ii)%deallocate_filter
  enddo
  call pinhole_sol%deallocate_lens
  call camera_sol%deallocate_camera_perspective_static
end subroutine teardown
!> Tests ---------------------------------------------------------------------
!> Tools ---------------------------------------------------------------------
!>----------------------------------------------------------------------------
end module mod_camera_perspective_static_mpi_test

