#define DSET PLT
#define DSET_S 'ne_PLT_'
program print_GRC
use mpi_mod
use mod_openadas
use constants, only: K_BOLTZ, EL_CHG
implicit none

integer :: ierr, provided, i, j, u
type(adf11_all) :: ad
character(len=2) :: s

integer, parameter :: n_T = 10000
real*8, parameter :: T0 = 1d-1 ! eV
real*8, parameter :: T1 = 1d4 ! eV

real*8 :: T, out

call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)

ad = read_adf11(0, '96_ne')

do i=0,ad%n_Z
  write(s,'(i2)') i
  open(newunit=u, file=DSET_S//trim(s)//'.dat')
  do j=1,size(ad%DSET%temperature,1)
    write(u,'(2g16.8)') ad%DSET%temperature(j), ad%DSET%GRC(size(ad%DSET%density,1),j,i)
  end do
  close(unit=u)

  ! Spline fit
  open(newunit=u, file=DSET_S//trim(s)//'_spline.dat')
  do j=1,n_T
    T = T0 * (T1/T0)**(real(j-1,8)/real(n_T-1,8)) * (EL_CHG / K_BOLTZ)
    call ad%DSET%interp(i, ad%DSET%density(size(ad%DSET%density,1)), log10(T), out)
    write(u,'(2g16.8)') T, out
  end do
  close(unit=u)
end do
end program print_GRC
