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
! Other effective methods using recursive fractional equations are reported in:
! B.R. Fabijonas et al. J. Comp. App. Math., vol.161, p. 179, 2003
! J. Rappoport, Mathematical Software for Modified Bessel Functions, in
! Mathematical Software - ICMS 2014, p. 352, Springer, 2014
! I.J. Thompson & A.R. Barnett, Comp. Phys. Comm., vol.47, p. 245, 1987
! N.M. Temme, J. Comp. Phys. vol.19, p.324, 1975
! W.J. Lentz, App. Optics, vol.15, p.668, 1976
! T. Takekawa, ArXiv: 2108.11560v2,  2021
! D.E. Amos, ACM. Trans. Math. Soft., vol.12, p.265, 1986
! D.E. Amos, ACM. Trans. Math. Soft., vol.21, p.388, 1995
! Other possibilities which my improve performance is to use external libraries:
! BOOST C++: www.boost.org <- requires isobinding between C++ and F90
! ALGORITHM 644: www.netlib.org/amos <- implement D.E. Amos in F77
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
!   ierr: (integer) error code, 0: success, 1: wrong inputs
!                   2: Lentz's algorithm did not converge
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
  integer,intent(inout)            :: ierr
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
!   ierr: (integer) error code, 0: success, 2: Lentz's algorithm did not converge
! outputs:
!		Knu:	(real8)(Nx) modified bessel function of the 2nd kind
!		dKnu:	(real8)(Nx) derivative of the modified bessel function 2nd kind
!   ierr: (integer) error code, 0: success, 2: Lentz's algorithm did not converge
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
  integer,intent(inout)            :: ierr
  ! internal variables
	integer               :: N_rec  !< number of downward recursions
  integer               :: Nx_ge_minx,Nx_lt_minx   !< N of x larger and smaller than x_val_min
  integer,dimension(Nx) :: ids_ge_minx,ids_lt_minx !< x-ids for array reconstruction
  real*8                :: mu                      !< order of the recurrence
  real*8,dimension(Nx)  :: x_ge_minx,x_lt_minx  !< values larger and smaller than x_val_min
  real*8,dimension(Nx)  :: kmu,k1,kmu_ge_minx
	real*8,dimension(Nx)  :: k1_ge_minx,kmu_lt_minx
	real*8,dimension(Nx)  :: k1_lt_minx

  ! separate variables in larger and smaller than x_val_min
  N_rec=floor(nu+5.d-1); mu=nu-N_rec;
  call  split_array_value(Nx,x,x_val_min,Nx_lt_minx,Nx_ge_minx,&
	x_lt_minx,x_ge_minx,ids_lt_minx,ids_ge_minx)
  ! compute the modified bessel function 2nd kind for x<x_val_min
	if(Nx_lt_minx.gt.0) then
	  call besselk_posxnu_lt_xval(Nx_lt_minx,x_lt_minx,mu,&
		kmu_lt_minx(1:Nx_lt_minx),k1_lt_minx(1:Nx_lt_minx),&
		tol,max_it,ierr)
	  kmu(ids_lt_minx(1:Nx_lt_minx)) = kmu_lt_minx(1:Nx_lt_minx)
		k1(ids_lt_minx(1:Nx_lt_minx)) = k1_lt_minx(1:Nx_lt_minx)
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
	call comp_besselk_posxnu(N_rec,Nx,mu,nu,1.d0/x,kmu,k1,knu,dknu)

end subroutine besselk_posxnu

! comp_besselk_posxnu computes the modified bessel function of the
! 2nd kind after having resolved the continuous fractional equation
! inputs:
!   N_rec: (integer) number of recursion
!   Nx:    (integer) number of x-values
!   mu:    (real8) fractional part of the bessel function order
!   nu:    (real8) bessel function order
!   xi:    (real8)(Nx) inverse of the vector x: 1/x
!   kmu:   (real8)(Nx) modfied bessel fct. 2nd kind at order mu
!   k1:    (real8)(Nx) derivative modified bessel fct 2nd kind order mu
! outputs:
!   knu:    (real8)(Nx) modified bessel fct 2nd kind order nu
!   dknu:   (real8)(Nx) derivative modified bessel fct 2nd kind order nu
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
	do ii=1,N_rec
	  ktemp = 2.d0*(mu+ii)*xi*k1 + kmu
		kmu = k1 
		k1 = ktemp
	enddo
	knu = kmu
	dknu = nu*xi*kmu-k1

end subroutine comp_besselk_posxnu

! bessekl_posxnu_lt_val computes the continuous fractional equation
! of the modified bessel function for the branch x<xval_min
! inputs:
!   Nx_lt:  (integer) number of x elements < x_val
!   x_lt:   (real8)(Nx) x values < x_val
!   mu:     (real8) fractional part of the bessel function order
!   tol:    (real8) tolerance for convergence
!   max_it: (integer) maximum number of iterations
!   ierr: (integer) error code, 0: success, 2: Lentz's algorithm did not converge
! outputs:
!   kmu:   (real8)(Nx) modfied bessel fct. 2nd kind at order mu
!   k1:    (real8)(Nx) derivative modified bessel fct 2nd kind order mu
!   ierr:  (integer) error code, 0: success, 2: Lentz's algorithm did not converge! 
subroutine besselk_posxnu_lt_xval(Nx_lt,x_lt,mu,kmu,k1,tol,max_it,ierr)
  implicit none
	
	! inputs
	integer,intent(in) :: Nx_lt,max_it
	real*8,dimension(Nx_lt),intent(in) :: x_lt
	real*8,intent(in) :: mu,tol
	! outputs
	integer,intent(inout) :: ierr
	real*8,dimension(Nx_lt),intent(out) :: kmu,k1
	! variables
	integer :: ii
  real*8 :: fact,gam1,gam2,gampl,gammi
	real*8,dimension(Nx_lt) :: c,d,e,p,q,ff,fact2,summ,summ1
  real*8,dimension(Nx_lt) :: del,del1
	
  ! evaluation of Gamm1 and Gamm2
	gammi = 1.d0/gamma(1.d0-mu)
	gampl = 1.d0/gamma(1.d0+mu)
	gam1 = 5.d-1*(gammi-gampl)/mu
	gam2 = 5.d-1*(gammi+gampl)

	! initialize the variables for the recurrent fraction solver
	fact = 1.d0
	if((PI*mu).ge.tol) fact = PI*mu/sin(PI*mu)
  d=-log(5.d-1*x_lt); e = mu*d;
	fact2 =1.d0
	where(abs(e).ge.tol) fact2 = sinh(e)/e
  ff = fact*(gam1*cosh(e)+gam2*fact2*d); 
	summ = ff; e = exp(e); d = 2.5d-1*x_lt*x_lt;
	p = 5.d-1*e/gampl; q = 5.d-1/(e*gammi);
	c = 1.d0; summ1 = p; ii=1; del=1.d10;

	! compute the recurrent fraction
  do while((maxval(abs(del)).ge.(maxval(abs(summ))*tol)).and.(ii.le.max_it))
	  ff = (ii*ff+q+p)/(ii*ii-mu*mu)
		c = c*d/real(ii); p=p/(ii-mu); q=q/(ii+mu);
		del = c*ff; summ = summ + del;
		del1 = c*(p-ii*ff); summ1 = summ1+del1;
		ii = ii+1;
	enddo
  ! check for convergence
	if(ii.gt.max_it) then
	  ierr = 2
		write(*,*) "Error fractional continuous solver for x<x_val: not converged! "
	endif

	! comput kmu and its derivatives
	kmu = summ;
	k1 = 2.d0*summ1/x_lt;

end subroutine besselk_posxnu_lt_xval

! bessekl_posxnu_ge_xval computes the continuous fractional equation 
! of the modified bessel function for the branch x>=xval_min
! inputs:
!   Nx_ge:  (integer) number of x elements >= x_val
!   x_ge:   (real8)(Nx) x values >= x_val
!   mu:     (real8) fractional part of the bessel function order
!   tol:    (real8) tolerance for convergence
!   max_it: (integer) maximum number of iterations
!   ierr: (integer) error code, 0: success, 2: Lentz's algorithm did not converge
! outputs:
!   kmu:   (real8)(Nx) modfied bessel fct. 2nd kind at order mu
!   k1:    (real8)(Nx) derivative modified bessel fct 2nd kind order mu
!   ierr:  (integer) error code, 0: success, 2: Lentz's algorithm did not converge!
subroutine besselk_posxnu_ge_xval(Nx_ge,x_ge,mu,kmu,k1,tol,max_it,ierr)
  implicit none

	! inputs
	integer,intent(in) :: Nx_ge,max_it
	real*8,dimension(Nx_ge),intent(in) :: x_ge
	real*8,intent(in) :: mu,tol
	! outputs
	integer,intent(inout) :: ierr
	real*8,dimension(Nx_ge),intent(out) :: kmu,k1
	! variables
	integer :: ii
	real*8,dimension(Nx_ge) :: a,b,c,d,h
	real*8,dimension(Nx_ge) :: q,q1,q2,qnew,s
	real*8,dimension(Nx_ge) :: delh,dels

	! init the continuous fractional solver
	b = 2.d0*(1.d0+x_ge); d = 1.d0/b; delh = d;
	h = delh; q1 = 0.d0; q2 = 1.d0; c = 2.5d-1 - mu*mu;
	q = c; a = -c; s = 1.d0 + q*delh; ii = 2; dels = 1.d10;
	! compute the fractional equation (as in numerical recipes)
	do while((maxval(abs(dels/s)).ge.tol).and.(ii.le.max_it))
	  a = a - 2.d0*(ii-1)
		c = -a*c/ii
		qnew = (q1-b*q2)/a
		q1 = q2; q2 = qnew;
		q = q + c*qnew; b = b + 2.d0
		d = 1.d0/(b + a*d)
		delh = (b*d-1.d0)*delh; h = h + delh;
    dels = q*delh; s = s + dels; 
		ii = ii + 1;
	enddo
	! check for convergence
	if(ii.gt.max_it) then
	  ierr = 2
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
	ids_lt_minx = -1; ids_ge_minx = -1;
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
