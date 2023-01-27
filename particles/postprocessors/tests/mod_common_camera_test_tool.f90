!> mod_common_camera_test_tool contains variables and
!> procedures used by multiple camera unit test modules
module mod_common_camera_test_tool
implicit none

private
public :: generate_image_plane_variables

!> Variable and data types -------------------------
!> Interfaces --------------------------------------
contains
!> Procedures --------------------------------------
!> sample the image plane variables
!> inputs:
!>   n_x:                (integer) size of the position vector
!>   n_planes:           (integer) number of image plane to test
!>   mirror_xy_interval: (integer)(2) interval (0,1) for selecting if a plane
!>                       should be mirrored or not
!>   plane_distance:     (real8) distance of the image plane from the pupil
!>   costheta_interval:  (real8)(2) selection interval of the plane center colatitude
!>   phi_interval:       (real8)(2) selection interval of the plane center azimuth
!>   half_angle_lowbnd:  (real8)(3) image plane dimension lowerbounds
!>   half_angle_uppbnd:  (real8)(3) image plane dimension upperbounds
!>   center_pos_lowbnd:  (real8)(3) pupil position lowerbounds
!>   center_pos_uppbnd:  (real8)(3) pupil position upperbounds
!> outputs:
!>   mirror_xy:         (integer)(2,n_planes) if a plane should be mirrored or not
!>   half_angle:        (real8)(3,n_panes) image plane angular coordinates
!>   image_plane_coord: (real8)(n_x,n_planes) positions of the image plane
!>                       centers in spherical coordinates w.r.t. the pupil
!>   pupil_positions:   (real8)(n_x,n_planes) position in xyz coordinates
!>                      of the camera pupil for each plane
subroutine generate_image_plane_variables(n_x,n_planes,mirror_xy_interval,&
plane_distance,costheta_interval,phi_interval,half_angle_lowbnd,&
half_angle_uppbnd,center_pos_lowbnd,center_pos_uppbnd,mirror_xy,&
half_angle,image_plane_coords,pupil_positions)
  use mod_gnu_rng,  only: gnu_rng_interval
  use mod_sampling, only: sample_uniform_sphere
  implicit none
  !> inputs:
  integer,intent(in)               :: n_x,n_planes
  integer,dimension(2),intent(in)  :: mirror_xy_interval
  real*8,intent(in)                :: plane_distance
  real*8,dimension(2),intent(in)   :: costheta_interval,phi_interval
  real*8,dimension(3),intent(in)   :: half_angle_lowbnd,half_angle_uppbnd 
  real*8,dimension(n_x),intent(in) :: center_pos_lowbnd,center_pos_uppbnd
  !> outputs:
  integer,dimension(2,n_planes),intent(out)  :: mirror_xy
  real*8,dimension(3,n_planes),intent(out)   :: half_angle
  real*8,dimension(n_x,n_planes),intent(out) :: image_plane_coords
  real*8,dimension(n_x,n_planes),intent(out) :: pupil_positions
  !> variables:
  integer             :: ii
  real*8,dimension(3) :: rand
  !> generation routine
  call gnu_rng_interval(2,n_planes,mirror_xy_interval,mirror_xy)
  do ii=1,n_planes
    call random_number(rand)
    image_plane_coords(:,ii) = sample_uniform_sphere(&
    plane_distance,costheta_interval,phi_interval,rand)
    call gnu_rng_interval(3,half_angle_lowbnd,&
    half_angle_uppbnd,half_angle(:,ii))
    call gnu_rng_interval(n_x,center_pos_lowbnd,&
    center_pos_uppbnd,pupil_positions(:,ii))
  enddo
end subroutine generate_image_plane_variables
!> -------------------------------------------------
end module mod_common_camera_test_tool
