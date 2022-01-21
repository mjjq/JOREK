!> mod_camera_perspective_static contains all variables and
!> procedures defining a perspective camera which does not
!> vary with time (static)
module mod_camera_perspective_static
use mod_camera_static, only: camera_static
implicit none
  
private
public :: camera_perspective_static

!> Variables and type definitions ----------------------------
type,extends(camera_static) :: camera_perspective_static
  real*8,dimension(:,:),allocatable :: points_on_lens
  real*8,dimension(:),allocatable   :: points_on_lens_pdf 
  contains
  procedure,pass(camera_inout) :: init_camera => init_camera_perspective_static
  procedure,pass(camera_inout) :: generate_points_on_lens_pdf => &
  generate_points_on_lens_static_perspective
end type camera_perspective_static

!> Interfaces ------------------------------------------------
contains
!> Procedures ------------------------------------------------
subroutine init_camera_perspective_static(camera_inout,lens_inout,&
n_int_param,n_real_param,int_param,real_param)
  use mod_camera,   only: camera
  use mod_geometry, only: define_plane_from_half_angles
  use mod_lens,     only: lens
  implicit none
  !> inputs-outputs:
  class(camera),intent(inout) :: camera_inout
  class(lens),intent(inout)   :: lens_inout
  !> inputs:
  integer,intent(in)                          :: n_int_param,n_real_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,dimension(:),allocatable,intent(in)  :: real_param

  !> initialise variables
end subroutine init_camera_perspective_static

!> generate points on lens and retrive their pdf
!> inputs:
!>   camera_inout: (camera) camera with initialised points on lens
!>   lens_inout:   (lens) lens model for generating points
!>   n_points_in:  (integer),optional
!> outputs:
!>   camera_out:   (camera) camera with initialised points on lens
subroutine generate_points_on_lens_static_perspective(camera_inout,&
lens_inout,n_points_in)
  use mod_lens,         only: lens
  use mod_pinhole_lens, only: pinhole_lens
  implicit none
  !> inputs-outputs:
  class(camera),intent(inout) :: camera_inout
  class(lens),intent(inout)   :: lens_inout
  !> sample the lens
  select type (cam=>camera_inout)
    type is(camera_perspective_static)
    call lens%sampling(cam%n_points_on_lens,cam%points_on_lens)
    call lens%pdf(cam%n_points_on_lens,cam%pdf_points_on_lens)
  end select
end subroutine generate_points_on_lens_static_perspective

!>------------------------------------------------------------
end module mod_camera_perspective_static
