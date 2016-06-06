!> Module to perform different kinds of sampling tricks
!! including sampling from gaussian distributions and 
!! sampling in cylindrical coordinates.
module mod_sampling
  implicit none
  public :: transform_uniform_cylindrical
  public :: boxmueller_transform
  private ! all other methods are private
contains
  !> Transform three uniform random numbers in [0,1] to 
  !! Uniform random numbers in cylindrical coordinates (R,Z,Phi)
  subroutine transform_uniform_cylindrical(ran3, Rbox, Zbox, Phibox, R, Z, Phi)
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
  !! gaussian-distributed random numbers with mean 0 and sigma 1
  !! Using the box-muller method (very slow!)
  function boxmueller_transform(ran) result(out)
    implicit none
    real*8, dimension(:), intent(in) :: ran
    real*8, parameter :: TWOPI = 6.2831853071795864769d0 
    real*8, dimension(size(ran,1)) :: out
    integer :: i

    do i=1,size(ran,1)/2
      out(i:i+1) = sqrt(-2*log(ran(2*i))) * &
          (/cos(TWOPI*ran(2*i+1)), sin(TWOPI*ran(2*i+1))/)
    end do
  end function boxmueller_transform
end module mod_sampling
