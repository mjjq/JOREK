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
  contains
  procedure,pass(camera_inout) :: init_camera => init_camera_perspective_static
  procedure,pass(camera_inout) :: generate_points_on_lens_pdf => &
  generate_points_on_lens_static_perspective
end type camera_perspective_static

!> Interfaces ------------------------------------------------
contains
!> Procedures ------------------------------------------------
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
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,dimension(:),allocatable,intent(in)  :: real_param

  !> initialise variables
  !> lens samples are stored as x positions of the
  !> vertices while their pdfs as property
  camera_inout%n_property_vertex=1;
  call camera_inout%allocate_camera(1,int_param(1),&
  int_param(2),int_param(3),spectrum_inout%n_spectra)
  !> sample the lens
  call camera_inout%generate_points_on_lens_pdf(lens_inout,int_param(1))
end subroutine init_camera_perspective_static

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
  allocate(camera_inout%properties(camera_inout%n_property_vertex,camera_inout%n_points_on_lens,1))
  call lens_inout%sampling(camera_inout%n_vertices,camera_inout%x(:,:,1))
  call lens_inout%pdf(camera_inout%n_vertices,camera_inout%x(:,:,1),&
  camera_inout%properties(1,:,1))
end subroutine generate_points_on_lens_static_perspective

!>------------------------------------------------------------
end module mod_camera_perspective_static
