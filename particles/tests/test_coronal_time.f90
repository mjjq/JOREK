!> This program calculates the coronal equilibrium and writes output suitable for gnuplot
program test_coronal_time
use openadas
use mod_coronal
implicit none

character(len=6), parameter :: suffix = "50_w"
real*8, parameter  :: tstep = 1.d-12
integer, parameter :: nstep = 10000
real*8, parameter  :: density = 20.d0 ! log10 density in m^-3
real*8, parameter  :: temperature = 7.677 ! log10 temperature in K: 7.677 = 4.1 keV

type(type_ADF11_all) :: adf11
real*8, dimension(:), allocatable :: p ! population in one equilibrium
integer :: i, it

adf11   = read_adf11(suffix)                                    ! read openadas data for ionisation, recombination and radiation rates

allocate(p(0:adf11%n_Z))

open(unit=10,status="replace",file="coronal_time.txt")
p = 1.d0

! Write output in format:
! ROW1x COL1y VAL
! (blank)
! ROW2x COL1y VAL
do i=1,nstep
  call coronal_timestep(adf11, p, tstep, density, temperature)
  write(*,"(i4,g16.8)") i, maxloc(p)
  do it=lbound(p,1),ubound(p,1)
    write(10,"(3g16.8)") i, it, p(it)/sum(p)
  enddo
  write(10,*)
enddo
close(unit=10)
write(*,*) "Output written to coronal_time.txt"
end program test_coronal_time
