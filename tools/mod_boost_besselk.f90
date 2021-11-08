! The module boost_besselk contains all interfaces
! and procedures for importing and bining the 
! boost function cyl_bessel_k computing
! the modified bessel function of the second kind
! with fractional order.
! WARNING: cyl_bessel_k(nu,x) accept only scalal
!          inputs
module mod_boost_besselk

implicit none

private
public ::  boost_besselk

! binding interface to the boost cyl_bessel_k
interface 
  function besselk_cpp(nu,x) bind(C,name="boost_besselk_cpp")
    use iso_c_binding, only: c_double
    implicit none
    real(c_double) :: besselk_cpp
    real(c_double) :: nu,x
  end function besselk_cpp
end interface

! interface for the different procedure in the module
interface boost_besselk
    module procedure besselk, besselk_x_array, &
    besselk_nu_array, besselk_x_nu_array
end interface boost_besselk

contains

! besselk is a specialization of the cyl_bessel_k
! function of boost to double datatypes
! inputs:
!   nu: (real8) bessel function fractional order
!   x:  (real8) value at which the bessel function
!       is computed
! outputs:
!   bknu: (real8) value of the bessel functions
subroutine besselk(nu,x,bknu)
  implicit none 
  real*8, intent(in) :: nu,x
  ! outputs
  real*8,intent(out) :: bknu
  bknu = besselk_cpp(nu,x)
end subroutine besselk

! besselk_x_array computes the modified bessel
! function of the second kind for an array of x
! inputs:
!   Nx: (integer) number of x values
!   nu: (real8) bessel function fractional order
!   x:  (real8)(Nx) array of x values
! outputs:
!   bknu: (real8)(Nx) array of bessel functions
subroutine besselk_x_array(Nx,nu,x,bknu)
  implicit none
  ! inputs
  integer,intent(in)              :: Nx
  real*8,intent(in)               :: nu
  real*8,dimension(Nx),intent(in) :: x
  ! outputs
  real*8,dimension(Nx),intent(out) :: bknu
  ! variables
  integer :: ii

  do ii=1,Nx
    bknu(ii) = besselk_cpp(nu,x(ii))
  enddo
end subroutine besselk_x_array

! besselk_nu_array computes the modified bessel
! function of the second kind for an array of nu
! inputs:
!   Nnu: (integer) number of fractional orders
!   nu:  (real)(Nu) array of fractional orders
!   x:   (real8) value at which the bessel function
!        is computed
! outputs:
!   bknu: (real8)(Nu) array of bessel functions
subroutine besselk_nu_array(Nnu,nu,x,bknu)
  implicit none
  ! inputs
  integer,intent(in)               :: Nnu
  real*8,dimension(Nnu),intent(in) :: nu
  real*8,intent(in)                :: x
  !outputs
  real*8,dimension(Nnu),intent(out) :: bknu
  ! variables
  integer :: ii

  do ii=1,Nnu
    bknu(ii) = besselk_cpp(nu(ii),x)
  enddo
end subroutine besselk_nu_array

! besselk_x_nu_array computes the modified bessel
! function of the second kind for an array of x
! and nu being the x the first index
! inputs:
!   Nx:  (integer) number of x values
!   Nnu: (integer) number of fractional orders
!   nu:  (real)(Nu) array of fractional orders
!   x:   (real8)(Nx) array of x values
! outputs:
!   bknu: (real8)(Nx,Nu) array of bessel functions
subroutine besselk_x_nu_array(Nx,Nnu,nu,x,bknu)
  implicit none
  ! inputs
  integer,intent(in)               :: Nx,Nnu
  real*8,dimension(Nnu),intent(in) :: nu
  real*8,dimension(Nx),intent(in)  :: x
  ! outputs
  real*8,dimension(Nx,Nnu),intent(out) :: bknu
  ! variables
  integer :: ii,jj

  do jj=1,Nnu
    do ii=1,Nx
      bknu(ii,jj) = besselk_cpp(nu(jj),x(ii))
    enddo
  enddo

end subroutine besselk_x_nu_array

end module mod_boost_besselk
