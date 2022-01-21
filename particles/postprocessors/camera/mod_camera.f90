!> mod_camera extends mod_vertices containing the basic
!> variables, datatypes and procedures common to the
!> different camera models. Each camera is supposed to have
!> one pupil and one visual plane hence, camera systems having
!> multiple pupils with multiple visual planes must be
!> simulated as different synthetic cameras
module mod_camera
use mod_vertices, only: vertices
implicit none

private
public :: camera

!> Variable and type definitions -----------------------------------
type,abstract,extends(vertices) :: camera
   !> the x variables contains the global position of the 
   !> visual points on the lens / film for each frame
   !> number of spectra and pixels
   integer              :: size_point_on_lens_pdf !< size a point on a lens
   integer              :: n_points_on_lens_pdf !< size and number of points on a lens
   integer,dimension(3) :: n_pixels_spectra
   real*8  :: exposure_time !< exposure time for each camera frame
   !> array containing the pixel intensity per each time
   real*8,dimension(:,:,:),allocatable :: pixel_intensities
  contains
   procedure(int_init_camera),pass(camera_inout),deferred     :: init_camera
   procedure(int_gen_points_lens),pass(camera_inout),deferred :: generate_points_on_lens
   procedure,pass(camera_inout)                               :: allocate_camera
   procedure,pass(camera_inout)                               :: deallocate_camera 
end type camera

!> Interfaces ------------------------------------------------------
interface
  !> interface of the deferred procedure init_camera which initialises a 
  !> camera type.
  !> inputs:
  !>   camera_inout: (camera) camera to be initialised
  !>   lens_inout:   (lens)   the lens to be used for the camera
  !>   n_int_param:  (integer) number of integer parameters to be passed
  !>   n_real_param: (integer) number of real parameters to be passed
  !>   int_param:    (integer)(n_int_param) integer parameter array
  !>   real_param:   (real8)(n_real_param) real parameter array
  !> outputs:
  !>   camera_inout: (camera) initialised camera
  subroutine int_init_camera(camera_inout,lens_inout,n_int_param,&
  n_real_param,int_param,real_param)
    use mod_lens, only: lens
    implicit none
    !> inputs-outputs
    class(camera),intent(inout) :: camera_inout
    class(lens),intent(inout)   :: lens_inout
    !> inputs
    integer,intent(in)                          :: n_int_param,n_real_param
    integer,dimension(:),allocatable,intent(in) :: int_param
    real*8,dimension(:),allocatable,intent(in)  :: real_param
  end subroutine int_init_camera

  !> generate a set of points on a lens
  !> inputs:
  !>   camera_inout: (camera) camera model used for sampling
  !>   lens_inout:   (lens) lens from which points are sampled
  !> outputs:
  !>   camera_inout: (camera) camera with sampled points on lens
  !>   lens_inout:   (lens) lens from which points are sampled
  subroutine int_gen_points_lens,pass(camera_inout,lens_inout)
    use mod_lens, only: lens
    implicit none
    !> inputs-outputs:
    class(camera),intent(inout) :: camera_inout
    class(lens),intent(inout)   :: lens_inout
  end subroutine int_gen_points_lens
end interface

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
  camera_inout%n_pixels_spectra = (/n_spectra,n_pixels_x,n_pixels_y/)
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
  camera_inout%n_pixels_spectra = 0; camera_inout%exposure_time = 0.d0;
end subroutine deallocate_camera

!>------------------------------------------------------------------
end module mod_camera
