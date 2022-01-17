!> mod_geometry contains variables and procedures for solving
!> basic / common geometrical problems such as intersections
module mod_geometry
implicit none

private
public :: compute_global_cart_coord_plane_points
public :: compute_plane_line_intersection_cart

!> Interfaces -------------------------------------------------
interface compute_global_cart_coord_plane_points
  module procedure compute_global_cart_coord_plane_points_r8
end interface compute_global_cart_coord_plane_points

interface compute_plane_line_intersection_cart
  module procedure compute_plane_line_intersect_cart_points_r8
end interface compute_plane_line_intersection_cart

contains
!> Procedures -------------------------------------------------
!> compute the global cartesian coordinates of a point on 
!> a plane given the nodes defining the plane and the point
!> local coordinates (double precision)
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
!>-------------------------------------------------------------
end module mod_geometry
