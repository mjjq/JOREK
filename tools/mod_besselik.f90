! The module besselIK implements routines for computing the fractional
! modified bessel function of the first (I) and second (K) type.
! As first test the method reported in:
! W.H. Press et al., Numerical Recipes in FORTRAN 77, vol.1, 2nd Ed., 
! Cambridge University Press, 1992
! Eventually, performance might be improved resolving the integral 
! formulation of the modified bessel functions, as reported in:
! M. Abramowitz and I. A. Stegun, Handbook of Mathematical
! tables, National Bureau of Standards, 10th Ed, 1972 and in further
! publications integrating using the Gauss-Laguerre integration method.
! At the moment, only methods for variable x>0 and order nu>0 are implemented.
module mod_besselik
use constants,only :: PI
implicit none

private
public besselik

! module variables
real,constant :: fp_min=1.d-30  !< value near the minimum floating point value
real,constant :: x_val_min=2.d0 !< x value smaller than which the small x 
                                !< approximation of the bessel function is used

contains

! besselk_posxnu computes the modified bessel function of the
! second kind with fracional order for positive x and order nu
! inputs
!		Nx:			   (int) number of x positions
!		x:			   (real8)(Nx) position at which Inu and Knu are evaluated
!		nu:			   (real8) bessel function order
!		tol_in:		 (real8)(optional) stop tolerance: default 1.e-16
!		max_it_in: (int)(optional) maximum number of iterations (default: 100000)
! outputs:
!		Knu:	(real8)(Nx) modified bessel function of the 2nd kind
!		dKnu:	(real8)(Nx) derivative of the modified bessel function 2nd kind
!   ierr: (integer) error code, 0: success, 1: wrong inputs
!                   2: Lentz's algorithm did not converge
subroutine besselk_posxnu(Nx,x,nu,Knu,dKnu,ierr,tol_in,max_it_in)
  use constants,only :: PI
  implicit none

  ! inputs
  integer,intent(in)              :: Nx
	real*8,dimension(Nx),intent(in) :: x
	real*8,intent(in)               :: nu
	real*8,intent(in),optional      :: tol_in
  integer*8,intent(in),optional   :: max_it_in
  ! outputs
	real*8,dimension(Nx),intent(out) :: Knu,dKnu
  integer                          :: ierr
  ! internal variables
  integer                         :: ii,it
  integer                         :: N_ge_minx,N_lt_minx !< N of x larger and smaller than x_val_min
  real*8                          :: err
  real*8                          :: mu !< order of the recurrence
  real*8,dimension(Nx)            :: x_ge_minx,x_lt_minx !< values larger and smaller than x_val_min

  ! check inputs
  max_it = 100000; tol_in = 1.d-16;
  if(present(max_it_in)) max_it = max_it_in
  if(present(tol_in)) tol = tol_in
  if(any(x.le.0.d0)) then
    ierr = 1
    return
  endif
  if(any(nu.lt.0.d0)) then
    ierr = 1
    return
  endif 

  ! separate variables in larger and smaller than x_val_min
  mu=nu-int(nu+5.d-1)
  call  split_array_value(Nx,x,x_val_min,N_lt_minx,N_ge_minx,x_lt_minx,x_ge_minx)
  ! compute the modified bessel function 2nd kind for x<x_val_min
  ! compute the modified bessel function 2nd kind for x>=x_val_min

end subroutine besselik_posxnu

pure subroutine split_array_value(Nx,x,x_val,N_lt_minx,&
N_ge_minx,x_lt_minx,x_ge_minx)
  implicit none
  
  ! inputs
  integer,intent(in)              :: Nx
  real*8,dimension(Nx),intent(in) :: x
  real*8,intent(in)               :: x_val
  ! outputs
  integer,intent(out) :: N_lt_minx,N_ge_minx
  real*8,dimension(Nx),intent(out) :: x_lt_minx,x_ge_minx

  N_lt_minx = count(x.lt.x_val_min)
  N_ge_minx = Nx-N_lt_minx
  where(x.lt.x_val_min)
    x_lt_minx(1:N_lt_minx) = x
  elsewhere
    x_ge_minx(1:N_ge_minx) = x
  endwhere

end subroutine split_array_value

end module mod_besselik_posxnu
