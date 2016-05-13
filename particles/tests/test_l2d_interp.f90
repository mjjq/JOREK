!> This program runs some testcases for the linear 2d interpolation
!> TODO: port to FRUIT or another framework
program test_l2d_interp
use openadas
implicit none

integer, parameter :: nx = 2, ny = 3
integer :: ix, iy
real*8, parameter :: ax = 2.5, ay = 7.5 !< Slopes in x and y
real*8 :: tx(nx)
real*8 :: ty(ny)
real*8 :: f(nx,ny)

! Prepare function and nodes
do ix=1,nx
  tx(ix) = real(ix,8)
  do iy=1,ny
    ty(iy) = real(iy,8)
    f(ix,iy) = tx(ix)*ax + ty(iy)*ay
  enddo
enddo

! Test a few cases
! Fully within the domain
call test_interp(1.5d0,1.5d0)

! Outside the domain
call test_interp(0.d0,0.d0)

! On one of the nodes
call test_interp(1.5d0,2.d0)

! On both of the nodes
call test_interp(2.d0,2.d0)

contains
subroutine test_interp(x,y)
implicit none

real*8, parameter :: tol = 1.d-12
real*8, intent(in) :: x, y
real*8 :: err

err = L2Dinterp(tx,ty,f,x,y) - (x*ax + y*ay)
if (err < tol) then
  !write(*,*) "Test succes: x=", x, " y=", y, " error=", err
else
  write(*,*) "Test failed: x=", x, " y=", y, " error=", err
endif
end subroutine test_interp
end program test_l2d_interp
