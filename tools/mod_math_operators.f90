!> This module contains multi-purposes mathematical operators
module mod_math_operators

  implicit none
  private
  public :: cross_product

interface cross_product
  module procedure cross_product_r4
  module procedure cross_product_r8
end interface cross_product
  
contains

!> The vector cross product single precision
pure function cross_product_r4(a, b)
  real*4, dimension(3) :: cross_product_r4
  real*4, dimension(3), intent(in) :: a, b

  cross_product_r4(1) = a(2) * b(3) - a(3) * b(2)
  cross_product_r4(2) = a(3) * b(1) - a(1) * b(3)
  cross_product_r4(3) = a(1) * b(2) - a(2) * b(1)

end function cross_product_r4
  
!> The vector cross product double precision
pure function cross_product_r8(a, b)
  real*8, dimension(3) :: cross_product_r8
  real*8, dimension(3), intent(in) :: a, b

  cross_product_r8(1) = a(2) * b(3) - a(3) * b(2)
  cross_product_r8(2) = a(3) * b(1) - a(1) * b(3)
  cross_product_r8(3) = a(1) * b(2) - a(2) * b(1)

end function cross_product_r8
  
end module mod_math_operators
