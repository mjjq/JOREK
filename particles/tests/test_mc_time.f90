!> This program calculates the charge state equilibrium with a MC method
program test_mc_time
use openadas
use mod_ionisation_recombination
implicit none

character(len=6), parameter :: suffix = "50_w"
real*8, parameter  :: tstep = 1.d-12
integer, parameter :: nstep = 10000
integer, parameter :: n_particles = 1000
real*8, parameter  :: density = 20.d0 ! log10 density in m^-3
real*8, parameter  :: temperature = 7.677 ! log10 temperature in K: 7.677 = 4.1 keV

type(type_ADF11_all) :: adf11
integer, dimension(:), allocatable :: z ! Charge states of all particles
integer, dimension(:), allocatable :: p ! population at one time
integer :: i, it

adf11   = read_adf11(suffix)            ! read openadas data for ionisation, recombination and radiation rates
call random_seed()

allocate(p(0:adf11%n_Z))
allocate(z(n_particles))

open(unit=10,status="replace",file="mc_time.txt")
do i=1,n_particles
  z(i) = int(adf11%n_Z * real(i,8)/real(n_particles,8))
enddo

! Write output in format:
! ROW1x COL1y VAL
! (blank)
! ROW2x COL1y VAL
do i=1,nstep
  p = 0
  do it=1,n_particles
    z(it) = new_charge(z(it), adf11, density, temperature, tstep)
    p(z(it)) = p(z(it)) + 1
  enddo
  do it=lbound(p,1),ubound(p,1)
    write(10,"(3g16.8)") i, it, real(p(it),8)/real(sum(p),8)
  enddo
  if (mod(i,1000) .eq. 0) write(*,*) i, maxloc(p), maxval(p)
enddo
close(unit=10)
write(*,*) "Output written to mc_time.txt"
end program test_mc_time
