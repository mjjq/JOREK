!> Module to perform different kinds of sampling tricks
!> including sampling from gaussian distributions and 
!> sampling in cylindrical coordinates.
!>
module mod_sampling
  implicit none
  private
  public :: transform_uniform_cylindrical
  public :: boxmueller_transform
  public :: sample_chi_squared_3
contains
  !> Transform three uniform random numbers in [0,1] to 
  !> Uniform random numbers in cylindrical coordinates (R,Z,Phi)
  pure subroutine transform_uniform_cylindrical(ran3, Rbox, Zbox, Phibox, R, Z, Phi)
    implicit none
    real*8, dimension(3), intent(in) :: ran3
    real*8, intent(in), dimension(2) :: Rbox, Zbox, Phibox
    real*8, intent(out) :: R, Z, Phi
    ! Use inversion sampling to correct for cylindrical coordinates
    ! r = sqrt(rand() (B^2 - A^2) + A^2) for min and max radius A and B
    R   = sqrt(ran3(1) * (Rbox(2)**2-Rbox(1)**2) + Rbox(1)**2)
    Z   = (Zbox(2)-Zbox(1))*ran3(2) + Zbox(1)
    phi = (Phibox(2)-Phibox(1))*ran3(3) + Phibox(1)
  end subroutine transform_uniform_cylindrical

  !> Transform 2N uniform random numbers in [0,1] to
  !> gaussian-distributed random numbers with mean 0 and sigma 1
  !> Using the box-muller method (very slow!)
  pure function boxmueller_transform(ran) result(out)
    implicit none
    real*8, dimension(:), intent(in) :: ran
    real*8, parameter :: TWOPI = 6.2831853071795864769d0 
    real*8, dimension(size(ran,1)) :: out
    integer :: i

    do i=1,size(ran,1),2
      out(i:i+1) = sqrt(-2*log(ran(i))) * &
          (/cos(TWOPI*ran(i+1)), sin(TWOPI*ran(i+1))/)
    end do
  end function boxmueller_transform

  !> Transform a uniformly distributed number u on [0,1] into a chi^2(3)-distributed
  !> number by inverse transform sampling.
  !> This function might be a bit expensive. There are cheaper ways of making
  !> chi-squared distributions from squares of normally distributed variables.
  !> If you are using quasi-random sequences those are not suitable however.
  pure function sample_chi_squared_3(u) result(x)
    use mod_rootfinding, only: newtons_method
    real*8, intent(in)   :: u !< Uniformly distributed number in [0,1]
    real*8               :: x !< Chi-squared(3) distributed number
    real*8               :: x0 !< Initial guess
    integer              :: ierr
    real*8, parameter    :: guess_a = 0.24851051d0, guess_b = 1.11289237d0

    ! Generate a guess by inverting a nearby distribution:
    ! u = tanh(a*x)**b
    x0 = atanh(u**(1.d0/guess_b))/guess_a
    
    call newtons_method(f=CDF_chi_squared_3, &
                       df=PDF_chi_squared_3, &
                       y0=u, x0=x0, x=x, ierr=ierr)
    ! Dangerous: ignore ierr for now
  end function sample_chi_squared_3

  pure function PDF_chi_squared_3(x) result(P)
    real*8, intent(in) :: x
    real*8             :: P
    real*8, parameter  :: PI = atan2(0.d0, -1.d0)
    ! Continue for negative values
    P = sign(sqrt(abs(x))*exp(-abs(x)*0.5d0)/sqrt(2.d0*PI), x)
  end function PDF_chi_squared_3

  pure function CDF_chi_squared_3(x) result(C)
    real*8, intent(in) :: x
    real*8             :: C
    real*8, parameter  :: PI = atan2(0.d0, -1.d0)
    ! Continue for negative values, mirror about y=-x
    C = sign(erf(sqrt(abs(x))/sqrt(2.d0)) - sqrt(2.d0/PI)*exp(-abs(x)*0.5d0)*sqrt(abs(x)), x)
  end function CDF_chi_squared_3

  !> Second derivative of the CDF for a chi_squared(3) distribution.
  !> Use this if you would like to use Halley's method to sample.
  !> In my tests it is faster to use Newton-Raphson iteration though.
  pure function PDF_prime_chi_squared_3(x) result(P)
    real*8, intent(in) :: x
    real*8             :: P
    real*8, parameter  :: PI = atan2(0.d0, -1.d0)
    ! Workaround if we go out of domain
    P = sign((1-abs(x))*exp(-abs(x)*0.5d0)/(2.d0*sqrt(2.d0*PI*abs(x))), x)
  end function PDF_prime_chi_squared_3
end module mod_sampling
