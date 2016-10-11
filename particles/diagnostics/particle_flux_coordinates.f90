!> Calculate a histogram of flux (p(1)) at all particle positions
!> Run as particle_flux_coordinates jorek_file.rst particle_file.h5 < in_jorek
!> Output is written to flux_filenum.txt
program particle_flux_coordinates
use phys_module
use data_structure
use basis_at_gaussian
use mod_particles
use mod_particle_io
use mod_particle_diagnostics
use mpi
implicit none

real*8, parameter :: binstart = -0.95, &
                     binend   = -0.15
integer, parameter :: n_bins = 100

type (type_particle_list) :: particle_list
integer    :: ierr, provided
character*25 :: particle_file, restart_file, output_file
integer    :: i_tor, i
type (type_node_list)   , pointer :: node_list
type (type_element_list), pointer :: element_list
real*8, allocatable, dimension(:) :: fluxcoords
logical, allocatable, dimension(:) :: mask
integer, allocatable, dimension(:) :: tmp
integer, allocatable, dimension(:) :: histogram

call MPI_Init_thread(MPI_THREAD_SINGLE, provided, ierr)
write(*,*) '***************************************'
write(*,*) '* JOREK particle flux coordinates     *'
write(*,*) '***************************************'

call initialise_parameters(0, "__NO_FILENAME__")

! Get the filename as the first cli argument
if (command_argument_count() < 2) then
  write(*,*) "Expected two filename arguments: restart_file.rst and particle_file.rst"
  call exit(1)
endif
call get_command_argument(1, restart_file)
call get_command_argument(2, particle_file)

allocate(node_list)
allocate(element_list)
do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
enddo
call import_binary_restart(node_list,element_list, restart_file, rst_format, ierr)
if (ierr .ne. 0) call exit(1)
call initialise_basis                              ! define the basis functions at the Gaussian points



call import_particles(particle_file, particle_list)
output_file = 'flux'//particle_file(5:index(particle_file,'.h5',.true.))//'txt' !  .true. searches backwards
allocate(fluxcoords(particle_list%n_particles),mask(particle_list%n_particles))

! Count number of particles in each element
call get_particle_flux_coordinates(node_list,element_list,particle_list,fluxcoords, mask)

allocate(tmp(count(mask,1)))

! Find bins
tmp = nint((pack(fluxcoords,mask) - binstart)/(binend - binstart)*real(n_bins,8))

! Calculate histogram
allocate(histogram(n_bins))
histogram = 0
!$omp parallel do default(none) &
!$omp shared(tmp) &
!$omp private(i) &
!$omp reduction(+:histogram)
do i=1,size(tmp,1)
  if (tmp(i) .ge. 1 .and. tmp(i) .le. n_bins) histogram(tmp(i)) = histogram(tmp(i)) + 1
enddo
!$omp end parallel do

open(file=output_file,status="replace",unit=21,access="stream",form='formatted')
write(21,'(A)') "# Flux, Count"
do i=1,n_bins
  write(21,'(g16.8,i6)') (binend-binstart)*((real(i,8)-0.5d0)/real(n_bins,8))+binstart, histogram(i)
enddo
close(21)
write(*,*) "Done counting, wrote output to ", output_file
call MPI_Finalize(ierr)
end program particle_flux_coordinates
