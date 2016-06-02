!> This program calculates the coronal equilibrium and writes output suitable for gnuplot
program test_coronal_eq
use openadas
use mod_coronal
implicit none

character(len=6), parameter :: suffix = "50_w"

type(type_ADF11_all) :: adf11
type(type_coronal)   :: cor

integer :: id, it
real*8 :: Te(4), Z

adf11   = read_adf11(suffix)                                    ! read openadas data for ionisation, recombination and radiation rates
cor     = coronal_equilibrium(adf11)                            ! calculate the coronal equilibria from the adas data

! Write output in format:
! ROW1x COL1y VAL
! (blank)
! ROW2x COL1y VAL

open(unit=10,status="replace",file="coronal_eq.txt")
do id=1,size(cor%density,1)
  do it=1,size(cor%temperature,1)
    write(10,"(3g16.8)") cor%density(id), cor%temperature(it), cor%Z(id, it)
  enddo
  write(10,*)
enddo
close(unit=10)

Te = (/5.42635, 6.88417, 6.89714, 7.48787/)
do id=1,4
  call interpolate_coronal(cor, 19.d0, Te(id), Z)
  write(*,*) Z
end do
end program test_coronal_eq
