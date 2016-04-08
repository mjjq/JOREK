!> Calculate flux (p(1)) at all particle positions
!! Run as particle_flux_coordinates jorek_file.rst particle_file.rst < in_jorek
!! Output is written to flux_filenum.dat
program particle_flux_coordinates
use phys_module
use data_structure
use basis_at_gaussian
use mod_particles
use mod_import_export_particles
use mpi_mod
implicit none

type (type_particle_list) :: particle_list
integer    :: ierr, provided
character*25 :: particle_file, restart_file, output_file
integer    :: i_tor
type (type_node_list)   , pointer :: node_list
type (type_element_list), pointer :: element_list

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
call import_binary_restart(node_list,element_list, restart_file, rst_format, ierr)
if (ierr .ne. 0) call exit(1)
call initialise_basis                              ! define the basis functions at the Gaussian points

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
enddo


call import_particles(particle_file, particle_list)
output_file = 'flux'//particle_file(5:index(particle_file,'.rst',.true.))//'txt' !  .true. searches backwards
call write_particle_flux_coordinates(node_list,element_list,particle_list,output_file)
write(*,*) "Done counting, wrote output to ", output_file

call MPI_Finalize(ierr)
contains
subroutine write_particle_flux_coordinates(node_list,element_list,particle_list,filename)
use phys_module
use data_structure
use basis_at_gaussian ! for HZ (initialise_basis must be called before use)
use mod_particles
use constants
implicit none

!> Input parameters
type(type_node_list), intent(in)      :: node_list
type(type_element_list), intent(in)   :: element_list
type (type_particle_list), intent(in) :: particle_list
character*(*), intent(in)             :: filename

integer :: i, j
real*8 :: R, R_s, R_t, Z, Z_s, Z_t
real*8, dimension(1) :: P, P_s, P_t, P_phi
integer :: i_elm
real*8, allocatable :: flux(:)

real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)

real*8  :: total_volume, total_area, volume, area, xjac, wst, phif, weight
integer :: ms, mt

type(type_node) :: node
type(type_element) :: element


allocate(flux(particle_list%n_particles))
! Count number of particles in each element
do i=1,particle_list%n_particles
  if (particle_list%particle(i)%lost) cycle ! skip this iteration
  i_elm  = particle_list%particle(i)%i_elm
  if (i_elm .lt. 1) cycle
  call interp_PRZ(node_list, element_list, i_elm, (/1/), 1, particle_list%particle(i)%st(1),particle_list%particle(i)%st(2), particle_list%particle(i)%x(3), P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
  write(*,"(7g16.8)") P
enddo

! Write output
open(file=filename,status="replace",unit=21,access="stream",form='formatted')
do i=1,particle_list%n_particles
  if (particle_list%particle(i)%lost) cycle ! skip this iteration
  i_elm  = particle_list%particle(i)%i_elm
  if (i_elm .lt. 1) cycle
  write(21,'(g16.8)') flux(i)
enddo
close(21)

end subroutine write_particle_flux_coordinates
end program particle_flux_coordinates
