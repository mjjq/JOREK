!> Module testing accuracy of calculating sines by applying De Moivre's formula
module test_moivre_series
use fruit
implicit none
contains

subroutine test_moivre_sincos
  integer, parameter :: N_max = 40, N_phi = 4
  real*8, parameter :: tol = 1d-14
  real*8, parameter :: phi0(N_phi) = [0.1d0, 0.5d0, 1.d0, 12.d0]
  character(len=3) :: i_s
  integer :: i, j
  real*8 :: e1r, e1i, enr, eni, enrtmp

  do i=1,N_phi
    e1r = cos(phi0(i))
    e1i = sin(phi0(i))
    enr = 1.d0
    eni = 0.d0
    do j=1,N_max
      enrtmp = enr*e1r - eni*e1i
      eni = eni*e1r + enr*e1i
      enr = enrtmp
      write(i_s,'(i3)') j
      call assert_equals(cos(real(j)*phi0(i)), enr, tol, 'Cosine match '//i_s)
      call assert_equals(sin(real(j)*phi0(i)), eni, tol, 'Sine match '//i_s)
    end do
  end do
end subroutine test_moivre_sincos
end module test_moivre_series
