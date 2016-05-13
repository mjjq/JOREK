!> Print interpolated values of the GRC coefficients at many points, for plotting and comparison
program test_linear_spline
use openadas
implicit none

integer, parameter :: n_points = 10000
integer, parameter :: z_test = 46
character(len=6), parameter :: suffix = "50_w"
real*8, parameter  :: density = 20.d0 ! log10 density in m^-3
real*8 :: temperature

type(type_ADF11_all) :: adf11
integer :: i

adf11   = read_adf11(suffix)           ! read openadas data for ionisation, recombination and radiation rates

open(unit=10,status="replace",file="test_interp.txt")
do i=-n_points/10,n_points+n_points/10
  temperature = real(i,8)/real(n_points,8)* &
      (adf11%ACD%temperature(ubound(adf11%ACD%temperature,1)) - adf11%ACD%temperature(1)) &
      + adf11%ACD%temperature(1)
  write(10,"(2g16.8)") temperature, GRC(adf11%ACD, z_test, density, temperature)
end do
close(unit=10)
write(*,*) "Output written to test_interp.txt"

open(unit=10,status="replace",file="test_interp_real_values.txt")
do i=lbound(adf11%ACD%temperature,1),ubound(adf11%ACD%temperature,1)
  write(10,"(2g16.8)") adf11%ACD%temperature(i), 10.d0**adf11%ACD%GRC(1,i,z_test)
end do
close(unit=10)
write(*,*) "Output written to test_interp_real_values.txt"
end program test_linear_spline
