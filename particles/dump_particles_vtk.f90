!> Dump all particles into a vtk file
!! Run as dump_particles_vtk part[0-9]*.rst < jorek_in
!! It uses the elements of file jorek_restart.rst
program dump_particles_vtk
use phys_module
use nodes_elements
use mod_particles
use mod_import_export_particles
use mod_vtk
use mpi_mod
implicit none

type (type_particle_list) :: particle_list
integer    :: ierr, provided
character*17 :: particle_file, restart_file
integer    :: i_t, i



call MPI_Init_thread(MPI_THREAD_SINGLE, provided, ierr)

write(*,*) '***************************************'
write(*,*) '* JOREK dump particles '
write(*,*) '***************************************'

call initialise_parameters(0, "__NO_FILENAME__")

do i_t=1, n_tor
  mode(i_t) = + int(i_t / 2) * n_period
  write(*,*) ' toroidal mode numbers : ',i_t,mode(i_t)
enddo

restart_file = 'jorek_restart.rst'
call import_binary_restart(node_list,element_list, restart_file, rst_format, ierr)
if (ierr .ne. 0) call exit(1)

! Get the filename as the first cli argument
if (command_argument_count() < 1) then
  write(*,*) "Expected a filename argument"
  call exit(1)
endif

do i=1,command_argument_count()
  ! Get particle filename from commandline
  call get_command_argument(i, particle_file)
  call import_particles(particle_file, particle_list)

  ! Write the first value of node_list to a vtk file
  particle_file = 'pos'//particle_file(5:index(particle_file,'.rst',.true.))//'vtk' !  .true. searches backwards
  call particles_vtk(particle_list,particle_file)
enddo

call MPI_FINALIZE(IERR)
end program dump_particles_vtk
