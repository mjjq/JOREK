!> Calculate change in flux coordinate bins between two particle files
!! Run as particle_flux_coordinate_diffusion jorek_file.rst particles_reference.rst part$filenum.h5 < in_jorek
!! Output is written to diff_filenum.txt
program particle_flux_coordinate_diffusion
use phys_module
use data_structure
use basis_at_gaussian
use mod_particles
use mod_particle_io
use mod_particle_diagnostics
use mpi_mod
implicit none

integer, parameter :: n_bins = 100
real*8, parameter  :: binstart = -0.95d0, binend = -0.15d0

type (type_particle_list) :: particle_list_reference
type (type_particle_list) :: particle_list
integer    :: ierr, provided
character*25 :: particle_file, particle_reference_file, restart_file, output_file
integer    :: i_tor, i, bin
type (type_node_list)   , pointer :: node_list
type (type_element_list), pointer :: element_list
real*8, dimension(:), allocatable :: fluxcoord, fluxcoord_reference
logical, allocatable, dimension(:) :: mask, mask_reference
integer, allocatable, dimension(:) :: tmp

real*8, dimension(n_bins) :: mean, stddev

call MPI_Init_thread(MPI_THREAD_SINGLE, provided, ierr)
write(*,*) '***************************************'
write(*,*) '* JOREK particle flux coordinate diff *'
write(*,*) '***************************************'

! Get the filename as the first cli argument
if (command_argument_count() < 3) then
  write(*,*) "Expected three filename arguments: restart_file.rst particles_reference.rst part$filenum.h5"
  call exit(1)
endif
call get_command_argument(1, restart_file)
call get_command_argument(2, particle_reference_file)
call get_command_argument(3, particle_file)
call initialise_parameters(0, "__NO_FILENAME__")


allocate(node_list)
allocate(element_list)
do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
enddo
call import_binary_restart(node_list,element_list, restart_file, rst_format, ierr)
if (ierr .ne. 0) call exit(1)
call initialise_basis                              ! define the basis functions at the Gaussian points


call import_particles(particle_file, particle_list)
call import_particles(particle_reference_file, particle_list_reference)
output_file = 'diff'//particle_file(5:index(particle_file,'.h5',.true.))//'txt' !  .true. searches backwards
allocate(fluxcoord(particle_list%n_particles),mask(particle_list%n_particles))
allocate(fluxcoord_reference(particle_list_reference%n_particles),mask_reference(particle_list_reference%n_particles))

call get_particle_flux_coordinates(node_list,element_list,particle_list, fluxcoord, mask)
call get_particle_flux_coordinates(node_list,element_list,particle_list_reference, fluxcoord_reference, mask_reference)

if (.not. any(mask_reference)) then
  write(*,*) "No valid particles found in reference file, exiting"
  call MPI_Finalize(ierr)
  call exit(-1)
endif

allocate(tmp(particle_list_reference%n_particles))
! Calculate in which bin particles are
tmp = nint((fluxcoord_reference - binstart)/(binend - binstart)*real(n_bins,8))

mean = 0.d0
! Calculate mean per bin
!$omp parallel do default(none) &
!$    shared(particle_list_reference, fluxcoord, fluxcoord_reference, mask_reference, tmp) &
!$    private(i) &
!$    reduction(+:mean)
do i=1,particle_list_reference%n_particles
  if (mask_reference(i) .and. tmp(i) .ge. 1 .and. tmp(i) .le. n_bins) then
    mean(tmp(i)) = mean(tmp(i)) + (fluxcoord(i) - fluxcoord_reference(i))
  endif
enddo
!$end omp parallel do
mean = mean/count(mask_reference)

stddev = 0.d0
! Calculate stddev per bin
!$omp parallel do default(none) &
!$    shared(particle_list_reference, fluxcoord, fluxcoord_reference, mask_reference, mean, tmp) &
!$    private(i) &
!$    reduction(+:stddev)
do i=1,particle_list_reference%n_particles
  if (mask_reference(i) .and. tmp(i) .ge. 1 .and. tmp(i) .le. n_bins) then
    stddev(tmp(i)) = stddev(tmp(i)) + ((fluxcoord(i) - fluxcoord_reference(i))-mean(tmp(i)))**2
  endif
enddo
!$omp end parallel do
stddev = sqrt(stddev/count(mask_reference))

open(file=output_file,status="replace",unit=21,access="stream",form='formatted')
write(21,'(A)') "# bin_center, mean, stddev"
do i=1,n_bins
  write(21,'(3g16.8)') (binend-binstart)*((real(i,8)-0.5d0)/real(n_bins,8))+binstart, mean(i), stddev(i)
enddo
close(21)


write(*,*) "Done, wrote output to ", output_file
call MPI_Finalize(ierr)
end program particle_flux_coordinate_diffusion
