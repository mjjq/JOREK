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
use constants, only: PI
implicit none

private
public besselik
public split_array_value

! module variables
real*8,parameter :: fp_min=1.d-30  !< value near the minimum floating point value
real*8,parameter :: x_val_min=2.d0 !< x value smaller than which the small x 
                                   !< approximation of the bessel function is used

! interface for functios computing the modified bessel functions
! of fractional order
interface besselik
  module procedure besselk
end interface besselik

contains

! besselk computes the modified bessel function of the
! second kind with fracional order to be extended for treating 
! x<0 and xnu < 0 if required
! inputs
!		Nx:			   (int) number of x positions
!		x:			   (real8)(Nx) position at which Inu and Knu are evaluated
!		nu:			   (real8) bessel function order
!		tol_in:		 (real8) stop tolerance: default 1.e-16
!		max_it_in: (int) maximum number of iterations (default: 100000)
! outputs:
!		Knu:	(real8)(Nx) modified bessel function of the 2nd kind
!		dKnu:	(real8)(Nx) derivative of the modified bessel function 2nd kind
!   ierr: (integer) error code, 0: success, 1: wrong inputs
!                   2: Lentz's algorithm did not converge
subroutine besselk(Nx,x,nu,Knu,dKnu,ierr,tol_in,max_it_in)
  implicit none

  ! inputs
  integer,intent(in)              :: Nx
	real*8,dimension(Nx),intent(in) :: x
	real*8,intent(in)               :: nu
	real*8,intent(in),optional      :: tol_in
  integer,intent(in),optional     :: max_it_in
  ! outputs
	real*8,dimension(Nx),intent(out) :: Knu,dKnu
  integer                          :: ierr
	! variables
	integer :: max_it
	real*8  :: tol

  ! check inputs
  max_it = 1000000; tol = 1.d-16;
  if(present(max_it_in)) max_it = max_it_in
  if(present(tol_in)) tol = tol_in

  ! compute bessel function of the 2nd kind for x>0 and nu>0
	if((any(x.gt.0.d0)).and.(nu.gt.0)) then
	  call besselk_posxnu(Nx,x,nu,Knu,dKnu,ierr,tol,max_it)
	else
	  ierr = 1
    write(*,*) "Error besselk for x<=0 or nu<=0 not implemented yet"
	endif	
  
end subroutine besselk

! besselk_posxnu computes the modified bessel function of the
! second kind with fracional order for positive x and order nu
! inputs
!		Nx:			   (int) number of x positions
!		x:			   (real8)(Nx) position at which Inu and Knu are evaluated
!		nu:			   (real8) bessel function order
!		tol_in:		 (real8) stop tolerance: default 1.e-16
!		max_it_in: (int) maximum number of iterations (default: 100000)
! outputs:
!		Knu:	(real8)(Nx) modified bessel function of the 2nd kind
!		dKnu:	(real8)(Nx) derivative of the modified bessel function 2nd kind
!   ierr: (integer) error code, 0: success, 1: wrong inputs
!                   2: Lentz's algorithm did not converge
subroutine besselk_posxnu(Nx,x,nu,Knu,dKnu,ierr,tol,max_it)
  implicit none

  ! inputs
  integer,intent(in)              :: Nx
	real*8,dimension(Nx),intent(in) :: x
	real*8,intent(in)               :: nu
	real*8,intent(in)               :: tol
  integer,intent(in)              :: max_it
  ! outputs
	real*8,dimension(Nx),intent(out) :: Knu,dKnu
  integer                          :: ierr
  ! internal variables
	integer               :: N_rec  !< number of downward recursions
  integer               :: Nx_ge_minx,Nx_lt_minx   !< N of x larger and smaller than x_val_min
  integer,dimension(Nx) :: ids_ge_minx,ids_lt_minx !< x-ids for array reconstruction
  real*8                :: mu                      !< order of the recurrence
  real*8,dimension(Nx)  :: xi,x_ge_minx,x_lt_minx  !< values larger and smaller than x_val_min
  real*8,dimension(Nx)  :: kmu,k1,kmu_ge_minx
	real*8,dimension(Nx)  :: k1_ge_minx,kmu_lt_minx
	real*8,dimension(Nx)  :: k1_lt_minx

  ! separate variables in larger and smaller than x_val_min
  N_rec=int(nu+5.d-1); mu=nu-N_rec; xi = 1.d0/x;
  call  split_array_value(Nx,x,x_val_min,Nx_lt_minx,Nx_ge_minx,&
	x_lt_minx,x_ge_minx,ids_lt_minx,ids_ge_minx)
  ! compute the modified bessel function 2nd kind for x<x_val_min
	if(Nx_lt_minx.gt.0) then
	  !TODO
	  !kmu(ids_lt_minx(1:Nx_lt_minx)) = kmu_lt_minx(1:Nx_lt_minx)
		!k1(ids_lt_minx(1:Nx_lt_minx)) = k1_lt_minx(1:Nx_lt_minx)
	endif
  ! compute the modified bessel function 2nd kind for x>=x_val_min
	if(Nx_ge_minx.gt.0) then
	  call besselk_posxnu_ge_xval(Nx_ge_minx,x_ge_minx,mu,&
	  kmu_ge_minx(1:Nx_ge_minx),k1_ge_minx(1:Nx_ge_minx),&
	  tol,max_it,ierr)
	  kmu(ids_ge_minx(1:Nx_ge_minx)) = kmu_ge_minx(1:Nx_ge_minx)
		k1(ids_ge_minx(1:Nx_ge_minx)) = k1_ge_minx(1:Nx_ge_minx)
	endif
	! compute the modified bessel function of the 2nd kind
	call comp_besselk_posxnu(N_rec,Nx,mu,nu,xi,kmu,k1,knu,dknu)

end subroutine besselk_posxnu

! comp_besselk_posxnu computes the modified bessel function of the
! 2nd kind after having resolved the continuous fractional equation
subroutine comp_besselk_posxnu(N_rec,Nx,mu,nu,xi,kmu,k1,knu,dknu)
  implicit none

	! input
	integer,intent(in) :: N_rec,Nx
	real*8,intent(in)  :: mu,nu
	real*8,dimension(Nx),intent(in) :: xi
	! output
	real*8,dimension(Nx),intent(out) :: knu,dknu
	! input-output
	real*8,dimension(Nx),intent(inout) :: kmu,k1
	! variables
	integer :: ii
	real*8,dimension(Nx) :: kmup,ktemp

	! compute the modified bessel function 2nd kind
  kmup = mu*xi*kmu-k1
	do ii=1,N_rec
	  ktemp = 2.d0*(mu+ii)*xi*k1 + kmu
		kmu = k1 
		k1 = ktemp
	enddo
	knu = kmu
	dknu = nu*xi*kmu-k1

end subroutine comp_besselk_posxnu

! bessekl_posxnu_ge_xval computes the continuous fractional equation 
! of the modified bessel function for the branch x>=xval_min
subroutine besselk_posxnu_ge_xval(Nx_ge,x_ge,mu,kmu,k1,tol,max_it,ierr)
  implicit none

	! inputs
	integer,intent(in) :: Nx_ge,max_it
	real*8,dimension(Nx_ge),intent(in) :: x_ge
	real*8,intent(in) :: mu,tol
	! outputs
	integer,intent(out) :: ierr
	real*8,dimension(Nx_ge),intent(out) :: kmu,k1
	! variables
	integer :: ii
	real*8,dimension(Nx_ge) :: a,b,c,d,h
	real*8,dimension(Nx_ge) :: q,q1,q2,qnew,s
	real*8,dimension(Nx_ge) :: delh,dels

  ierr = 0
	! init the continuous fractional solver
	b = 2.d0*(1.d0+x_ge); d = 1.d0/b; delh = d;
	h = delh; q1 = 0.d0; q2 = 1.d0; c = 2.5d-1 - mu*mu;
	q = c; a = -c; s = 1.d0 + q*delh; ii = 2; dels = 1.d10;
	! compute the fractional equation (as in numerical recipes)
	do while((maxval(abs(dels/s)).ge.tol).and.(ii.le.max_it))
	  a = a - 2.d0*(ii-1)
		c = -a*c/ii
		qnew = (q1-q2)/a
		q1 = q2; q2 = qnew;
		q = q + c*qnew; b = b + 2.d0
		d = 1.d0/(b + a*d)
		delh = (b*d-1.d0)*delh; h = h + delh;
    dels = q*delh; s = s + dels; 
		ii = ii + 1;
	enddo
	! check for convergence
	if(ii.gt.max_it) then
	  ierr = 1
		write(*,*) "Error fractional continuous solver for x>=x_val: not converged! "
	endif
	! compute the contribution to the modified bessel function 2nd kind
	! omit the exp(-x) to scale all returned values by exp(x)
	kmu = sqrt(PI/(2.d0*x_ge))*exp(-x_ge)/s
	k1 = kmu*(mu+x_ge+5.d-1-(2.5d-1 - mu*mu)*h)/x_ge
end subroutine besselk_posxnu_ge_xval

! split_array_value an array into two different arrays
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
