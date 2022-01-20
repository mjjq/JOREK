!> mod_geometry contains variables and procedures for solving
!> basic / common geometrical problems such as intersections
module mod_geometry
implicit none

private
public :: compute_global_cart_coord_plane_points
public :: compute_plane_line_intersection_cart
public :: define_plane_from_half_angles

!> Interfaces -------------------------------------------------
interface compute_global_cart_coord_plane_points
  module procedure compute_global_cart_coord_plane_points_r8
end interface compute_global_cart_coord_plane_points

interface compute_plane_line_intersection_cart
  module procedure compute_plane_line_intersect_cart_points_r8
end interface compute_plane_line_intersection_cart

interface define_plane_from_half_angles
  module procedure define_direction_origin_plane_from_half_angles
  module procedure define_direction_plane_from_half_angles
  module procedure define_standard_plane_from_half_angles
end interface define_plane_from_half_angles

contains
!> Procedures -------------------------------------------------
!> compute the global cartesian coordinates of a point on 
!> a plane given the nodes defining the plane and the point
!> local coordinates (double precision)
!> P1 -> s -> P2
!> |           |
!> v           v
!> t           t
!> |           |
!> v           v
!> P3 -> s -> P4
!> inputs:
!>   pp: (real8)(3,3) points defining a plane in cartesian 
!>                    coordinates: pp(:,1) -> origin
!>                                 pp(:,2) -> s=1,t=0 node
!>                                 pp(:,3) -> s=0,t=1 node
!>   st: (real8)(2) position of a point on a plane in the plane
!>                  local coordinates
!> outpus:
!>   x: (real8)(3) position of a point on a plane in the global
!>                 cartesian coordinates
pure function compute_global_cart_coord_plane_points_r8(pp,st) result(x)
  implicit none
  real*8,dimension(2),intent(in)   :: st
  real*8,dimension(3,3),intent(in) :: pp
  real*8,dimension(3)              :: x
  x = pp(:,1)*(1.d0-st(1)-st(2))+pp(:,2)*st(1)+pp(:,3)*st(2)
end function compute_global_cart_coord_plane_points_r8

!> compute_plane_line_intersect_cart_points_r8 computes the 
!> intersection between a line and a plane in double precision
!> the plane must be defined by three points while the line
!> must be defined by two points. The function returns the
!> intersection coordinates in the local plane (s,t) and lines
!> (q) coordinate system. The intersection is found if the
!> local coordinates of the intersection points are in [0,1].
!> The north-east node is the origin of the plane
!> P1 -> s -> P2
!> |           |
!> v           v
!> t           t
!> |           |
!> v           v
!> P3 -> s -> P4 
!> inputs:
!>   pp: (real8)(3,3) points defining a plane in cartesian 
!>                    coordinates: pp(:,1) -> origin
!>                                 pp(:,2) -> s=1,t=0 node
!>                                 pp(:,3) -> s=0,t=1 node
!>   pl: (real8)(3,2) points defining a line in cartesian
!>                    coordinates: pl(:,1) -> origin
!>                                 pl(:,2) -> q=1 node
!> outputs:
!>   intersect: (bool) if true intersection found
!>   stq:       (real8)(3) plane(s,t) and line (1) local
!>              coordinates of the intersection
subroutine compute_plane_line_intersect_cart_points_r8(pp,pl,&
intersect,stq)
  use mod_math_operators, only: solve_3x3_linear_problem
  implicit none
  real*8,dimension(3,3),intent(in) :: pp
  real*8,dimension(3,2),intent(in) :: pl
  logical,intent(out)              :: intersect
  real*8,dimension(3),intent(out)  :: stq
  real*8,dimension(3,3) :: matrix
  
  matrix(:,1) = pp(:,2)-pp(:,1); matrix(:,2) = pp(:,3)-pp(:,1);
  matrix(:,3) = pl(:,1)-pl(:,2);
  call  solve_3x3_linear_problem(matrix,pl(:,1)-pp(:,1),stq)
  intersect = (all(stq.ge.0.d0)).and.(all(stq.le.1.d0))
end subroutine compute_plane_line_intersect_cart_points_r8

!> define a plane vertices given the half width, half height,
!> the spherical coordinates of the plane mid point in
!> the origin coordinate system and the origin coordinates
!> inputs:
!>   half_angles: (real8)(2) 1:half_width 2:half_height
!>   rthetaphi:   (real8)(3) spherical coordinates of the
!>                plane midpoin: rthetaphi(1): distance
!>                               rthetaphi(2): colatitude
!>                               rthetaphi(3): azimuth
!>   origin:      (real8)(3) position of the sphere origin
!> outputs
!>   pp: (real8)(3,3) points defining a plane in cartesian 
!>                    coordinates: pp(:,1) -> P1
!>                                 pp(:,2) -> P2
!>                                 pp(:,3) -> P3
subroutine define_direction_origin_plane_from_half_angles(half_angles,&
rthetaphi,origin,pp)
  implicit none
  real*8,dimension(2),intent(in) :: half_angles
  real*8,dimension(3),intent(in) :: rthetaphi,origin
  !> outputs:
  real*8,dimension(3,3),intent(out) :: pp
  !> variables
  integer :: ii
  !> compute directional plane vertices
  call define_direction_plane_from_half_angles(half_angles,rthetaphi,pp)
  !> translate vertices
  do ii=1,3
    pp(:,ii) = origin + pp(:,ii)
  enddo
end subroutine define_direction_origin_plane_from_half_angles


!> define a plane vertices given the half width, half height
!> and the spherical coordinates of the plane mid point in
!> the origin coordinate system
!> inputs:
!>   half_angles: (real8)(2) 1:half_width 2:half_height
!>   rthetaphi:   (real8)(3) spherical coordinates of the
!>                plane midpoin: rthetaphi(1): distance
!>                               rthetaphi(2): colatitude
!>                               rthetaphi(3): azimuth
!> outputs
!>   pp: (real8)(3,3) points defining a plane in cartesian 
!>                    coordinates: pp(:,1) -> P1
!>                                 pp(:,2) -> P2
!>                                 pp(:,3) -> P3
subroutine define_direction_plane_from_half_angles(half_angles,&
rthetaphi,pp)
  implicit none
  !> inputs:
  real*8,dimension(2),intent(in) :: half_angles
  real*8,dimension(3),intent(in) :: rthetaphi
  !> outputs:
  real*8,dimension(3,3),intent(out) :: pp
  !> variables
  integer :: ii
  real*8,dimension(2) :: cos_thetaphi,sin_thetaphi
  real*8,dimension(3,3) :: rot
  !> compute rotation transform
  cos_thetaphi = cos(rthetaphi(2:3));
  sin_thetaphi = sin(rthetaphi(2:3));
  rot(:,1) = (/-sin_thetaphi(2),cos_thetaphi(2),0.d0/)
  rot(:,2) = (/cos_thetaphi(1)*cos_thetaphi(2),&
             cos_thetaphi(1)*sin_thetaphi(2),-sin_thetaphi(1)/)
  rot(:,3) = (/sin_thetaphi(1)*cos_thetaphi(2),&
             sin_thetaphi(1)*sin_thetaphi(2),cos_thetaphi(1)/)
  !> compute plane vertices
  call define_standard_plane_from_half_angles(half_angles,pp)
  !> transform the vertices in the new positions
  do ii=1,3
    pp(:,ii) = rthetaphi(1)*matmul(rot,pp(:,ii))
  enddo
end subroutine define_direction_plane_from_half_angles

!> define a plane vertices given the half width and half height
!> angles in the standard reference systme: origin = (/0,0,0/)
!> and plane normal z=(/0,0,1/)
!> P1 -> width -> P2
!>   |           |
!>   v           v
!> height      height
!>   |           |
!>   v           v
!> P3 -> width -> P4 
!> inputs:
!>   half_angles: (real8)(2) 1:half_width 2:half_height
!> outputs
!>   pp: (real8)(3,3) points defining a plane in cartesian 
!>                    coordinates: pp(:,1) -> P1
!>                                 pp(:,2) -> P2
!>                                 pp(:,3) -> P3
subroutine define_standard_plane_from_half_angles(half_angles,pp)
  implicit none
  !> inputs:
  real*8,dimension(2),intent(in) :: half_angles
  !> outputs:
  real*8,dimension(3,3),intent(out) :: pp
  !> variables
  real*8,dimension(2) :: tan_angles
  !> compute plane vertex coordinates
  tan_angles = tan(half_angles)
  pp(:,1) = (/-tan_angles(1),tan_angles(2),1.d0/)
  pp(:,2) = (/tan_angles(1),tan_angles(2),1.d0/)
  pp(:,3) = (/-tan_angles(1),-tan_angles(2),1.d0/)
end subroutine define_standard_plane_from_half_angles

!>-------------------------------------------------------------
end module mod_geometry
