!> mod_geometry_test contains variables and procedures
!> used for testing the routines contained in mod_geometry
module mod_geometry_test
use fruit
implicit none

private
public :: run_fruit_geometry

!> Variables ----------------------------------------------
integer,parameter :: n_planes=11
integer,parameter :: n_nodes_plane=3
integer,parameter :: n_lines_per_test=23
integer,parameter :: n_tests=5
integer,parameter :: n_lines=n_tests*n_lines_per_test
real*8,parameter  :: tol_real8=7.5d-10
real*8,dimension(2),parameter  :: len_interval=(/2.d0,5.d2/)
real*8,dimension(3),parameter  :: ppxyz_lowbound=(/-1.d0,2.d0,-4.d1/)
real*8,dimension(3),parameter  :: ppxyz_uppbound=(/5.d0,2.5d1,-3.2d1/)
real*8,dimension(3),parameter  :: plxyz_lowbound=(/-3.d0,1.d1,4.d1/)
real*8,dimension(3),parameter  :: plxyz_uppbound=(/7.d0,5.d0,3.2d1/)
real*8,dimension(3),parameter  :: pl0_lowbound=(/-2.1d1,4.5d0,-9.3d1/)
real*8,dimension(3),parameter  :: pl0_uppbound=(/8.d1,1.1d1,5.2d1/)
real*8,dimension(2,n_tests),parameter :: st_lowbound=&
       reshape((/0.d0,0.d0,-3.d1,-4.d0,-5.1d1,1.1d0,&
       1.2d0,-8.d1,1.1d0,1.05d0/),shape(st_lowbound))
real*8,dimension(2,n_tests),parameter :: st_uppbound=&
       reshape((/1.d0,1.d0,-5.d-2,-4.d-3,-1.d-1,5.d1,&
       9.2d1,-1.d-3,1.1d1,3.05d0/),shape(st_uppbound))
logical,dimension(n_lines,n_planes)        :: intersect_sol
real*8,dimension(3,n_nodes_plane,n_planes) :: pp_sol
real*8,dimension(3,n_lines,n_planes)       :: stq_sol,pos_intersect_sol
real*8,dimension(3,2,n_lines,n_planes)     :: pl_sol
!> Interfaces ---------------------------------------------

contains
!> Fruit basket -------------------------------------------
!> basket having all set-up, tests and tearing-down routines
subroutine run_fruit_geometry()
  implicit none
  write(*,'(/A)') "  ... setting-up: geometry tests"
  call setup
  write(*,'(/A)') "  ... running: geometry tests"
  call test_compute_test_line_intersect_cart_points
  write(*,'(/A)') "  ... tearing-down: geometry tests"
end subroutine run_fruit_geometry

!> Set-up and tear-down -----------------------------------
subroutine setup()
  use mod_geometry, only: compute_global_cart_coord_plane_points
  use mod_gnu_rng,  only: gnu_rng_interval
  implicit none
  !> variables
  integer :: ii,jj,kk,id
  real*8              :: length
  real*8,dimension(2) :: st
  real*8,dimension(3) :: pos

  do ii=1,n_planes
    !> generate plane nodes
    do jj=1,n_nodes_plane
      call gnu_rng_interval(3,ppxyz_lowbound,&
      ppxyz_uppbound,pp_sol(:,jj,ii))
    enddo
    !> generate lines
    do jj=1,n_tests
      do kk=1,n_lines_per_test
        !> generate first line node
        id = (jj-1)*n_lines_per_test+kk
        call gnu_rng_interval(3,plxyz_lowbound,&
        plxyz_uppbound,pl_sol(:,1,id,ii))
        !> generate second line node
        call gnu_rng_interval(2,st_lowbound(:,jj),&
        st_uppbound(:,jj),st)
        call gnu_rng_interval(len_interval,length)
        pos = compute_global_cart_coord_plane_points(pp_sol(:,:,ii),st)
        pl_sol(:,2,id,ii) = pl_sol(:,1,id,ii) + &
        (pos-pl_sol(:,1,id,ii))*length
        !> check if intersection
        intersect_sol(id,ii)=(all(st.ge.0.d0).and.all(st.le.1.d0))
        !> store intersection position values
        pos_intersect_sol(:,id,ii) = pos
        stq_sol(:,id,ii) = (/st(1),st(2),norm2(pos-pl_sol(:,1,id,ii))/&
        norm2(pl_sol(:,2,id,ii)-pl_sol(:,1,id,ii))/)
      enddo
    enddo
  enddo
end subroutine setup

!> Tests --------------------------------------------------
!> test the procedure for finding the intersection between
!> a line and a plane
subroutine test_compute_test_line_intersect_cart_points()
  use mod_assert_equals_tools, only: assert_equals_rel_error
  use mod_geometry, only: compute_global_cart_coord_plane_points
  use mod_geometry, only: compute_plane_line_intersection_cart
  implicit none
  integer :: ii,jj
  logical :: intersect
  real*8,dimension(3) :: stq,pos

  do ii=1,n_planes
    do jj=1,n_lines
      !> compute intersection point in local and global coordinates
      call compute_plane_line_intersection_cart(&
      pp_sol(:,:,ii),pl_sol(:,:,jj,ii),intersect,stq)
      pos = compute_global_cart_coord_plane_points(&
      pp_sol(:,:,ii),stq(1:2))
      !> check solutions
      call assert_equals_rel_error(3,stq,stq_sol(:,jj,ii),tol_real8,&
      "Error compute plane line intersection cart: stq mismatch!")
      call assert_equals_rel_error(3,pos,pos_intersect_sol(:,jj,ii),tol_real8,&
      "Error compute plane line intersection cart: pos intersect mismatch!")
      call assert_equals(intersect,intersect_sol(jj,ii),&
      "Error compute plane line intersection cart: intersect mismatch!")
    enddo
  enddo
end subroutine test_compute_test_line_intersect_cart_points

!> Tools --------------------------------------------------
!>---------------------------------------------------------
end module mod_geometry_test

