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
use constants,only: PI
implicit none

private
public besselk_posxnu
public split_array_value

! module variables
real*8,parameter :: fp_min=1.d-30  !< value near the minimum floating point value
real*8,parameter :: x_val_min=2.d0 !< x value smaller than which the small x 
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
  integer               :: max_it
  integer               :: N_ge_minx,N_lt_minx     !< N of x larger and smaller than x_val_min
  integer,dimension(Nx) :: ids_ge_minx,ids_lt_minx !< x-ids for array reconstruction
  real*8                :: err,tol
  real*8                :: mu                      !< order of the recurrence
  real*8,dimension(Nx)  :: x_ge_minx,x_lt_minx     !< values larger and smaller than x_val_min

  ! check inputs
  max_it = 100000; tol = 1.d-16;
  if(present(max_it_in)) max_it = max_it_in
  if(present(tol_in)) tol = tol_in
  if(any(x.le.(0.d0))) then
    ierr = 1
    return
  endif
  if(nu.lt.0.d0) then
    ierr = 1
    return
  endif 

  ! separate variables in larger and smaller than x_val_min
  mu=nu-int(nu+5.d-1)
  call  split_array_value(Nx,x,x_val_min,N_lt_minx,N_ge_minx,&
	x_lt_minx,x_ge_minx,ids_lt_minx,ids_ge_minx)
  ! compute the modified bessel function 2nd kind for x<x_val_min
  ! compute the modified bessel function 2nd kind for x>=x_val_min

end subroutine besselk_posxnu

! This procedure split an array into two different arrays
! containing x values below and above a threshold value.
! The element ids in the original array are returned as well.
! There is probably a much better way for doing it!!!!
! inputs:
! outputs:
!   N_lt_minx:   (integer) number of x below threshold
!   N_ge_minx:   (integer) number of x above threshold
!   x_lt_minx:   (double)(Nx) array containing x below threshold
!   x_ge_minx:   (double)(Nx) array containing x above threshold
!   ids_lt_minx: (integer)(Nx) x-ids of x below threshold
!   ids_ge_minx: (integer)(Nx) x-ids of x above threshold
pure subroutine split_array_value(Nx,x,x_val,N_lt_minx,&
N_ge_minx,x_lt_minx,x_ge_minx,ids_lt_minx,ids_ge_minx)
  implicit none
  
  ! inputs
  integer,intent(in)              :: Nx
  real*8,dimension(Nx),intent(in) :: x
  real*8,intent(in)               :: x_val
  ! outputs
  integer,intent(out) :: N_lt_minx,N_ge_minx
  real*8,dimension(Nx),intent(out) :: x_lt_minx,x_ge_minx
  integer,dimension(Nx),intent(out) :: ids_lt_minx,ids_ge_minx
	! variables
	integer :: ii,it1,it2
	logical,dimension(Nx) :: mask_lt_minx

  ! extract in two different arrays x values abobe
	! and below the threshold, store their ids as well
  mask_lt_minx = (x.lt.x_val)
  N_lt_minx = 0; N_ge_minx = 0;
	do ii=1,Nx
	  if(mask_lt_minx(ii)) then
		  N_lt_minx = N_lt_minx + 1
			ids_lt_minx(N_lt_minx) = ii
		  x_lt_minx(N_lt_minx) = x(ii)
		else
		  N_ge_minx = N_ge_minx + 1
			ids_ge_minx(N_ge_minx) = ii
			x_ge_minx(N_ge_minx) = x(ii)
		endif
	enddo
  
end subroutine split_array_value

end module mod_besselik
