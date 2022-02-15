program print_Prad
use mpi_mod
use mod_openadas
use mod_coronal
implicit none

integer :: ierr, provided, i, j, u
type(adf11_all) :: ad
type(coronal) :: cor
character(len=2) :: s

call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)

ad = read_adf11(0, '96_ne')
cor = coronal(ad)

open(newunit=u,file='ne_Prad.dat')
do j=1,size(cor%temperature,1)
  write(u,'(2g16.8)') cor%temperature(j), cor%Prad(size(cor%density,1),j)
end do
close(unit=u)
end program print_Prad
