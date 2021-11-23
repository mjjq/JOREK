!> mod_test_functions contains some functions and their
!> analytical integrals used for testing
module mod_test_functions
implicit none

private
public :: sin2x, int_sin2x

!> Interfaces ----------------------------------------------------------------
interface sin2x
  module procedure sin2x_serial,sin2x_vector
end interface

interface int_sin2x
  module procedure int_sin2x_serial,int_sin2x_vector
end interface

contains

!> sin(x)^2 ------------------------------------------------------------------

!> sin^2(x) function
function sin2x_serial(x)
  implicit none
  real*8,intent(in) :: x
  real*8 :: sin2x_serial
  sin2x_serial = sin(x)*sin(x)
end function sin2x_serial

!> integral of the sin^2(x) function
function int_sin2x_serial(x)
  implicit none
  real*8,intent(in) :: x
  real*8 :: int_sin2x_serial
  int_sin2x_serial = 5.d-1*(x-sin(x)*cos(x))
end function int_sin2x_serial

!> sin^2(x) function
function sin2x_vector(N,x)
  implicit none
  integer,intent(in) :: N
  real*8,dimension(N),intent(in) :: x
  real*8,dimension(N) :: sin2x_vector
  sin2x_vector = sin(x)*sin(x)
end function sin2x_vector

!> integral of the sin^2(x) function
function int_sin2x_vector(N,x)
  implicit none
  integer,intent(in) :: N
  real*8,dimension(N),intent(in) :: x
  real*8,dimension(N) :: int_sin2x_vector
  int_sin2x_vector = 5.d-1*(x-sin(x)*cos(x))
end function int_sin2x_vector

!> ---------------------------------------------------------------------------

end module mod_test_functions
