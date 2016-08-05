!> Dump all particle positions into a vtk file
!> Run as `dump_particles_vtk part[0-9]*.rst`
program dump_particles_vtk
use mod_particles
use mod_particle_io
use mod_vtk
use_mpi
implicit none

type (type_particle_list) :: particle_list
integer    :: ierr, provided
character*17 :: particle_file
integer    :: i


call MPI_Init_thread(MPI_THREAD_SINGLE, provided, ierr)

write(*,*) '***************************************'
write(*,*) '* JOREK dump particles '
write(*,*) '***************************************'

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
