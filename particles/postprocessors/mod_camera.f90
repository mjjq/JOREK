!> mod_camera extends mod_vertices containing the basic
!> variables, datatypes and procedures common to the
!> different camera models. Each camera is supposed to have
!> one pupil and one visual plane hence, camera systems having
!> multiple pupils with multiple visual planes must be
!> simulated as different synthetic cameras
module mod_camera
use mod_vertices, only: n_x, vertices
implicit none

private
public :: camera

!> Variable and type definitions -----------------------------------

type,abstract,extends(vertices) :: camera
   !> the x variables contains the global position of the 
   !> visual points on the lens / film for each frame
   !> number of spectra and pixels
   integer,dimension(3) :: n_pixels_spectra
   real*8  :: exposure_time !< exposure time for each camera frame
   !> array containing the pixel intensity per each time
   real*8,dimension(:,:,:),allocatable :: pixel_intensities
  contains
   procedure,pass(camera_inout) :: initialise_camera
   procedure,pass(camera_inout) :: allocate_camera 
end type camera

!> Interfaces ------------------------------------------------------
contains
!> Procedures ------------------------------------------------------
!> allocate camera arrays and initialise them to zero
!> inputs:
!>   camera_inout: (camera) camera to be allocated
!>   n_times:      (integer) number of times
!>   n_vertices:   (integer) number of vertices
!>   n_spectra:    (integer) number of spectrum intervals
!>   n_pixels_x:   (integer) number of pixels in the x-direction
!>   n_pixels_y:   (integer) number of pixels in the y-direction
!> outputs:
!>   camera_inout: (camera) allocated camera
subroutine allocate_camera(camera_inout,n_times,n_vertices,&
n_spectra,n_pixels_x,n_pixels_y)
  implicit none
  !> inputs-outputs:
  class(camera),intent(inout) :: camera_inout
  !> inputs:
  integer,intent(in) :: n_times,n_vertices
  integer,intent(in) :: n_spectra,n_pixels_x,n_pixels_y
  
  !> allocate vertices
  call camera_inout%allocate_vertices(n_times,n_vertices)
  !> allocate camera
  if(allocated(camera_inout%pixel_intensities)) &
  deallocate(camera_inout%pixel_intensities)
  allocate(camera_inout%pixel_intensities(n_spectra,&
  n_pixels_x*n_pixels_y,n_times))
  camera_inout%pixel_intensities = 0.d0;
  camera_inout%n_pixels_camera = (/n_spectra,n_pixels_x,n_pixels_y/)
end subroutine allocate_camera

!> deallocate camera arrays and reset counters to 0
!> inputs:
!>   camera_inout: (camera) camera to deallocated
!> outputs:
!>   camera_inout: (camera) deallocated camera
subroutine deallocate_camera(camera_inout)
  implicit none
  !> inputs-outputs
  class(camera),intent(inout) :: camera_inout
  !> deallocare everything and reset counters
  if(allocated(camera_inout%pixel_intensities)) &
  deallocate(camera_inout%pixel_intensities)
  camera%n_pixels_spectra = 0; camera%exposure_time = 0.d0;
end subroutine deallocate_camera
!>------------------------------------------------------------------
end module mod_camera
