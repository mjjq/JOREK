module mod_import_export_particles
  use mod_particles
contains

!> Export all particles using MPI File IO
!! Writes a file with first the number of particles in an integer
!! And then this many copies of the particle datatype from
!! get_particle_derived_type
subroutine export_particles(particle_list,particle_file)
use mpi_mod
implicit none

type (type_particle_list), intent(in) :: particle_list
character*(*)            , intent(in) :: particle_file

integer              :: my_id, n_cpu, ierr, ierr2
integer, allocatable, dimension(:) :: particles_per_proc

! For MPI writing
integer                       :: fh, dtype, status(MPI_STATUS_SIZE)
character*(*), parameter      :: datarep = 'native'

! Debugging output
integer :: resultlen
character(len=MPI_MAX_ERROR_STRING) :: string

call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs

allocate(particles_per_proc(0:n_cpu-1))

if (my_id .eq. 0) then
  write(*,*) '***********************************'
  write(*,*) '*       export particles          *'
endif

! Find the number of particles on each node
call MPI_AllGather(particle_list%n_particles,1,MPI_INTEGER,particles_per_proc,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

call MPI_File_open(MPI_COMM_WORLD,particle_file,MPI_MODE_CREATE + MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
if (ierr .gt. 0) then
  call MPI_ERROR_STRING(ierr, string, resultlen, ierr2)
  write(*,*) "file open failed:", ierr, trim(string)
endif
call MPI_File_set_errhandler(fh, MPI_ERRORS_ARE_FATAL, ierr)

! Write header containing only the number of particles
if (my_id .eq. 0) then
  call MPI_File_write(fh, sum(particles_per_proc), 1, MPI_INTEGER, status, ierr)
endif

dtype = get_particle_derived_type()

! Calculate the displacement after the header = one unit of size MPI_INTEGER
call MPI_File_set_view(fh, 1*MPI_INTEGER, dtype, dtype, datarep, MPI_INFO_NULL, ierr)

! write_at sets the displacement to the number of particles already written in units of dtype
call MPI_File_write_at(fh, int(sum(particles_per_proc(0:my_id-1)),MPI_OFFSET_KIND), particle_list%particle, particle_list%n_particles, dtype, status, ierr)

call MPI_File_close(fh,ierr)

write(*,*) '*       particles exported        *'
write(*,*) '***********************************'
end subroutine export_particles


!> Export all particles using MPI File IO
!! Reads the number of particles from a file, determines a 
!! particle distribution over all processors and read this many
!! particles per processor.
subroutine import_particles(particle_file, particle_list)
use mpi_mod
implicit none

character*(*)            , intent(in)  :: particle_file
type (type_particle_list), intent(out) :: particle_list

integer              :: my_id, n_cpu, ierr, ierr2, n_particles
integer, allocatable, dimension(:) :: particles_per_proc

! For MPI writing
integer                       :: fh, dtype, status(MPI_STATUS_SIZE)
character*(*), parameter      :: datarep = 'native'

! Debugging output
integer :: resultlen
character(len=MPI_MAX_ERROR_STRING) :: string

call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs

allocate(particles_per_proc(0:n_cpu-1))

if (my_id .eq. 0) then
  write(*,*) '***********************************'
  write(*,*) '*       import particles          *'
endif

call MPI_File_open(MPI_COMM_WORLD,particle_file,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
if (ierr .gt. 0) then
  call MPI_ERROR_STRING(ierr, string, resultlen, ierr2)
  write(*,*) "file open failed:", ierr, trim(string)
endif
call MPI_File_set_errhandler(fh, MPI_ERRORS_ARE_FATAL, ierr)

! Read number of particles from header
if (my_id .eq. 0) then
  call MPI_File_read(fh, n_particles, 1, MPI_INTEGER, status, ierr)

  ! Divide particles over processors
  particles_per_proc(0)         = n_particles/n_cpu + modulo(n_particles, n_cpu)
  particles_per_proc(1:n_cpu-1) = n_particles/n_cpu
endif
! Broadcast this array to all processors
call MPI_Bcast(particles_per_proc, n_cpu, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

! Set n_particles for all cpus
particle_list%n_particles = particles_per_proc(my_id)
! And allocate the required space
allocate(particle_list%particle(particle_list%n_particles), stat=ierr)
if (ierr .gt. 0) write(*,"(i3,a,i8,a)") my_id, "unable to allocate particle_list%particle(", particle_list%n_particles, ")"

dtype = get_particle_derived_type()

! Calculate the displacement after the header = one unit of size MPI_INTEGER
call MPI_File_set_view(fh, 1*MPI_INTEGER, dtype, dtype, datarep, MPI_INFO_NULL, ierr)

! write_at sets the displacement to the number of particles already written in units of dtype
call MPI_File_read_at(fh, int(sum(particles_per_proc(0:my_id-1)),MPI_OFFSET_KIND), particle_list%particle, particle_list%n_particles, dtype, status, ierr)

call MPI_File_close(fh,ierr)

write(*,*) '*       particles imported        *'
write(*,*) '***********************************'
endsubroutine import_particles

end module mod_import_export_particles
