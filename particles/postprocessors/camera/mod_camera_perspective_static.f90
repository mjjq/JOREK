!> mod_camera_perspective_static contains all variables and
!> procedures defining a perspective camera which does not
!> vary with time (static)
module mod_camera_perspective_static
use mod_camera, only: camera
implicit none
  
private
public :: camera_perspective_static

!> Variables and type definitions ----------------------------
type,extends(camera) :: camera_perspective_static
  integer :: n_plane_points !< number of points of a plane
  real*8,dimension(2) :: pixel_size !< width and height of a pixel
  real*8,dimension(:,:),allocatable :: image_plane !< vertices of the image plane
  contains
  procedure,pass(camera_inout) :: init_camera => init_camera_perspective_static
  procedure,pass(camera_inout) :: generate_points_on_lens_pdf => &
  generate_points_on_lens_static_perspective
  procedure,pass(camera_inout) :: allocate_camera_perspective_static
  procedure,pass(camera_inout) :: deallocate_camera_perspective_static
  procedure,pass(camera_inout) :: define_image_plane_pixel_size
  procedure,pass(camera_inout) :: plane_to_pixel_local_coord
end type camera_perspective_static

!> Interfaces ------------------------------------------------
contains
!> Procedures ------------------------------------------------
!> inputs:
!>   camera_inout:   (camera) unallocated camera
!>   lens_inout:     (lens) camera lens object
!>   spectrum_inout: (spectrum_inout) camera spectrum object
!>   n_int_param:    (integer) number of integer parameter: 3
!>   n_real_param:   (integer) number of real parameter: 8
!>   int_praram:     (integer)(n_int_param) integer parameters:
!>                            1) number of lens samples
!>                            2) number of pixels in the x-direction
!>                            3) number of pixels in the y-direction
!>  real_param:      (real8)(n_int_param) real parameters:
!>                            1:2) image plane half widht and half
!>                                 height angles in the focal reference
!>                            3:5) distance betweem the image plane
!>                                 and the camera focal point
!>                            6:8) camera focal point position
!> outputs:
!>   camera_inout:   (camera) allocated camera
!>   lens_inout:     (lens) camera lens object
!>   spectrum_inout: (spectrum_inout) camera spectrum object
subroutine init_camera_perspective_static(camera_inout,lens_inout,&
spectrum_inout,n_int_param,n_real_param,int_param,real_param)
  use mod_spectra,  only: spectrum_base
  use mod_geometry, only: define_plane_from_half_angles
  use mod_lens,     only: lens
  implicit none
  !> inputs-outputs:
  class(camera_perspective_static),intent(inout) :: camera_inout
  class(lens),intent(inout)                      :: lens_inout
  class(spectrum_base),intent(inout)             :: spectrum_inout
  !> inputs:
  integer,intent(in)                          :: n_int_param,n_real_param
  integer,dimension(n_int_param),intent(in)  :: int_param
  real*8,dimension(n_real_param),intent(in)  :: real_param

  !> set variables
  camera_inout%n_property_vertex=1; camera_inout%n_plane_points=3;
  !> initialise attributes
  !> lens samples are stored as x positions of the
  !> vertices while their pdfs as property
  call camera_inout%allocate_camera_perspective_static(spectrum_inout,&
  n_int_param,int_param)
  !> define the image plane characteristics
  call camera_inout%define_image_plane_pixel_size(n_real_param,real_param)
  !> sample the lens
  call camera_inout%generate_points_on_lens_pdf(lens_inout,int_param(1))
  !> initialise the image plane vertices
  call define_plane_from_half_angles(real_param(1:2),real_param(3:5),&
  real_param(6:8),camera_inout%image_plane)
end subroutine init_camera_perspective_static

!> procedure used for allocating all attributes of camera perspective static
!> inputs:
!>   camera_inout:   (camera) unallocated camera
!>   spectrum_inout: (spectrum_inout) camera spectrum object
!>   n_int_param:    (integer) number of integer parameter: 3
!>   int_praram:     (integer)(n_int_param) integer parameters:
!>                            1) number of lens samples
!>                            2) number of pixels in the x-direction
!>                            3) number of pixels in the y-direction
!> outputs:
!>   camera_inout:   (camera) allocated camera
!>   spectrum_inout: (spectrum_inout) camera spectrum object
subroutine allocate_camera_perspective_static(camera_inout,&
spectrum_inout,n_int_param,int_param)
  use mod_spectra, only: spectrum_base
  implicit none
  !> inputs-outputs
  class(camera_perspective_static),intent(inout) :: camera_inout
  class(spectrum_base),intent(inout)             :: spectrum_inout
  !> inputs
  integer,intent(in)                        :: n_int_param
  integer,dimension(n_int_param),intent(in) :: int_param
  !> allocate camera base type
  call camera_inout%allocate_camera(1,int_param(1),&
  int_param(2),int_param(3),spectrum_inout%n_spectra)
  !> allocate attributes specific to camera_perspective_static
  if(allocated(camera_inout%image_plane)) then
    if(size(camera_inout%image_plane,2).ne.camera_inout%n_plane_points) &
    deallocate(camera_inout%image_plane)
  endif
  allocate(camera_inout%image_plane(camera_inout%n_x,&
  camera_inout%n_plane_points))
end subroutine allocate_camera_perspective_static

!> procedure used for deallocating all attributes of camera perspective static
!> inputs:
!>   camera_inout: (camera) allocated camera
!> outputs:
!>   camera_inout: (camera) deallocated camera
subroutine deallocate_camera_perspective_static(camera_inout)
  implicit none
  !> inputs-outputs
  class(camera_perspective_static),intent(inout) :: camera_inout
  !> deallocate variables
  call camera_inout%deallocate_camera
  if(allocated(camera_inout%image_plane)) deallocate(camera_inout%image_plane)
end subroutine deallocate_camera_perspective_static

!> generate points on lens and retrive their pdf
!> inputs:
!>   camera_inout: (camera) camera with initialised points on lens
!>   lens_inout:   (lens) lens model for generating points
!>   n_points_in:  (integer)(optional) number of points to sample
!>                                     default: 1000, pinhole: 1
!> outputs:
!>   camera_out:   (camera) camera with initialised points on lens
subroutine generate_points_on_lens_static_perspective(camera_inout,&
lens_inout,n_points_in)
  use mod_lens,         only: lens
  use mod_pinhole_lens, only: pinhole_lens
  implicit none
  !> inputs-outputs:
  class(camera_perspective_static),intent(inout) :: camera_inout
  class(lens),intent(inout)                      :: lens_inout
  !> inputs:
  integer,intent(in),optional                    :: n_points_in
  !> initialisations
  camera_inout%n_vertices = 1000
  if(present(n_points_in)) camera_inout%n_vertices = n_points_in
  !> check if the lens is a pinhole, overwrite n_points_lens
  select type(ln=>lens_inout)
    type is(pinhole_lens)
    camera_inout%n_vertices = 1
  end select
  camera_inout%n_active_vertices = camera_inout%n_vertices
  !> sample the lens
  if(allocated(camera_inout%x)) then
    if(camera_inout%n_vertices.ne.size(camera_inout%x,dim=2)) &
    deallocate(camera_inout%x)
  endif
  if(allocated(camera_inout%properties)) then
    if(camera_inout%n_vertices.ne.size(camera_inout%properties,dim=2)) &
    deallocate(camera_inout%properties)
  endif
  if(.not.allocated(camera_inout%x)) &
  allocate(camera_inout%x(camera_inout%n_x,camera_inout%n_vertices,1))
  if(.not.allocated(camera_inout%properties)) &
  allocate(camera_inout%properties(camera_inout%n_property_vertex,camera_inout%n_vertices,1))
  call lens_inout%sampling(camera_inout%n_vertices,camera_inout%x(:,:,1))
  call lens_inout%pdf(camera_inout%n_vertices,camera_inout%x(:,:,1),&
  camera_inout%properties(1,:,1))
end subroutine generate_points_on_lens_static_perspective

!> generate the image plane points and compute the pixel width and height
!> inputs:
!>  camera_inout: (camera_perspective_static) camera with unallocated image plane
!>  n_real_param: (integer) number of real parameters
!>  real_param:   (real8)(n_real_param) real parameters, order:
!>                1:2) plane width and height half angles
!>                3:5) plane position w.r.t. the pupil in spherical coordinates
!>                6:8) position of the pupil in cartesian coordinates
!> outputs:
!>  camera_inout: (camera_perspective_static) camera with defined image plane
subroutine define_image_plane_pixel_size(camera_inout,n_real_param,real_param)
  use mod_geometry, only: define_plane_from_half_angles
  implicit none
  !> inputs-outputs:
  class(camera_perspective_static),intent(inout) :: camera_inout
  !> inputs:
  integer,intent(in)                             :: n_real_param
  real*8,dimension(n_real_param),intent(in)      :: real_param
  !> define the plane from the width/height half angles, the distance
  !> from the pupil and the pupil position
  call define_plane_from_half_angles(real_param(1:2),real_param(3:5),&
  real_param(6:8),camera_inout%image_plane)
  !> compute the pixel width and heigh in the pixel reference system
  camera_inout%pixel_size = (/1.d0,1.d0/)/real(camera_inout%n_pixels_spectra(2:3),kind=8)
end subroutine define_image_plane_pixel_size

!> compute the pixel number and the position in the pixel local coordinates
!> of a point on the image plane (in the plane local coordinates)
!> inputs:
!>  camera_inout: (camera_perspective_static) camera with unallocated image plane
!>  st_plane:     (real8)(2) position in the plane local coordinates
!> outputs:
!>  camera_inout: (camera_perspective_static) camera with defined image plane
!>  i_pixel:      (integer)(2) pixel indices of the point on the plane (s,t)
!>  st_pixel:     (real8)(2) position in the pixel local coordinates
subroutine plane_to_pixel_local_coord(camera_inout,st_plane,i_pixel,st_pixel)
  implicit none
  !> inputs-outputs
  class(camera_perspective_static),intent(inout) :: camera_inout
  !> inputs
  real*8,dimension(2),intent(in)   :: st_plane
  !> outputs
  integer,dimension(2),intent(out) :: i_pixel
  real*8,dimension(2),intent(out)  :: st_pixel
  !> find the local pixel coordinates and find the position in the pixel local coordinates
  st_pixel = st_plane/camera_inout%pixel_size
  i_pixel = floor(st_pixel)
  st_pixel = st_pixel - real(i_pixel,kind=8)
  i_pixel = i_pixel + 1
end subroutine plane_to_pixel_local_coord

!>------------------------------------------------------------
end module mod_camera_perspective_static
