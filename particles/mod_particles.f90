!> Mod_particles contains the particle(_list) type,
!> logic for reading particle input namelists and broadcasting parameters
!> as well as a function to create an MPI derived particle type.
module mod_particles
  use parameters

  type type_particle
    real*8    :: st(n_dim)        !< particle position in the finite element (i_elm)
    real*8    :: x(3)             !< particle position in real space (R,Z,phi) [m] at t
    real*8    :: v(3)             !< particle velocity in (R,Z,phi) [m/s * sqrt(mu0 rho0)] at t-dt/2
    real*4    :: mass             !< mass [atomic mass units]
    real*4    :: weight           !< weight (i.e. number of particles)
    integer   :: i_elm            !< the index of the element containing the particle in the element_list
    integer*1 :: q                !< charge [e]
    integer*1 :: label            !< Particle type number (i in species(i))
    logical*1 :: lost             !< particle is active or lost
  end type type_particle

  type type_particle_list
    integer                                         :: n_particles  !< the number of particles in the list
    type (type_particle), dimension(:), allocatable :: particle     !< an allocatable list of particles
  end type type_particle_list

  ! Particle input file parameters
  integer, parameter :: N_species = 9   !< Maximum number of particle species
  integer :: species(N_species) = 0     !< Atomic numbers of each species (-1) for electrons
  real*4  :: atomic_mass(N_species)     !< Atomic mass of each species, in a.m.u.
  integer :: N_particles(N_species) = 0 !< Number of particles to initialize
  logical :: particle_GC(N_species) = .false. !< Is this a guiding center species?
  logical :: particle_ion_rec(N_species) = .true. !< Calculate ionisation/recombination events
  character(len=6)  :: adas_suffix(N_species) = '' !< Suffix for adas files to read in (ex: scd50_w.dat => 50_w)
  character(len=80) :: location_accept_function(N_species) = 'location_accept_any' !< Which function to use for particle position rejection sampling
  real*4  :: location_accept_parameters(1:9,1:N_species) = 0 !< Extra arguments for this function
  integer :: particle_seed(N_species) = 0 !< Seed for PCG random sequence used for particle init
  character(len=6) :: particle_initializer(N_species) = 'pcg32' !< Method to use for seeding particles (options: pcg32, sobol)

  real*8  :: t_step_particles !< the time step for the advance of the particles
  integer :: n_step_particles !< the number of time steps for the particles
  integer :: nout_particles   !< number of particle timestep between each output file
  logical :: write_energies   !< Output energies of all the particles every nout
  logical :: write_momenta    !< Output generalized toroidal momentum of all the particles every nout
  integer :: t_particles_begin = -1 !< Number of first JOREK restart file (if -1 use only jorek_restart)
  integer :: t_particles_end  !< Number of last JOREK restart file
  character(len=80) :: particle_restart_file = '' !< Particle restart file to read from
contains


!> Append a single particle to the list and grow it if needed
subroutine append_particle_to_list(particle_list, particle)
  implicit none

  type(type_particle_list), intent(inout) :: particle_list
  type(type_particle), intent(in) :: particle

  if (size(particle_list%particle,1) .lt. particle_list%n_particles+1) then
    call grow_particle_list(particle_list)
  endif
  particle_list%n_particles = particle_list%n_particles + 1
  particle_list%particle(particle_list%n_particles) = particle
end subroutine append_particle_to_list


!> Append a list of particles to the list and grow it if necessary
subroutine append_particles_to_list(particle_list, particles)
  implicit none

  type(type_particle_list), intent(inout) :: particle_list
  type(type_particle), dimension(:), intent(in) :: particles

  ! Grow it until it fits
  do while (size(particle_list%particle,1) .lt. particle_list%n_particles+size(particles,1))
    call grow_particle_list(particle_list)
  enddo
  particle_list%particle(particle_list%n_particles+1:particle_list%n_particles+size(particles,1)) = particles
  particle_list%n_particles = particle_list%n_particles + size(particles,1)
end subroutine append_particles_to_list


!> Grow the particle list by a factor of two
subroutine grow_particle_list(particle_list)
  implicit none
  type(type_particle_list), intent(inout) :: particle_list
  type(type_particle), dimension(:), allocatable :: temp

  allocate(temp(lbound(particle_list%particle,1):ubound(particle_list%particle,1)+size(particle_list%particle,1)))
  temp(lbound(particle_list%particle,1):ubound(particle_list%particle,1)) = particle_list%particle
  call move_alloc(from=temp,to=particle_list%particle) ! deallocates temp as well
end subroutine grow_particle_list


!> Read particle parameters from filename or stdin if my_id == 0
subroutine initialise_particle_parameters(my_id, filename)
implicit none
character(len=*), intent(in) :: filename
integer, intent(in) :: my_id
integer :: ierr

namelist /in2/ species, atomic_mass, N_particles, particle_GC, &
    adas_suffix, particle_seed, &
    location_accept_function, location_accept_parameters, &
    particle_initializer, &
    n_step_particles, t_step_particles, nout_particles, &
    write_energies, write_momenta, t_particles_begin, t_particles_end, &
    particle_restart_file, particle_ion_rec

if (my_id .eq. 0) then
  ! --- Read input parameters from namelist file or stdin.
  if (trim(filename) .ne. "__NO_FILENAME__" ) then
     open(42, file=filename, status='old', action='read', iostat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'ERROR: COULD NOT OPEN NAMELIST FILE "', trim(filename), '".'
      stop
    end if
    read(42,in2)
    close(42)
  else
    read(5,in2)
  endif
endif
end subroutine initialise_particle_parameters

subroutine broadcast_particle_parameters(my_id)
use mpi_mod
implicit none

! --- Routine parameters
integer, intent(in) :: my_id

! --- internal variables
integer                :: ierr, INT_EXT, IDBL_EXT, ISGL_EXT, ILOG_EXT, CHAR_EXT, position, bufsize
character, allocatable :: buffer(:)

!----------------------------------- one line would be enough if only MPI_TYPE_STRUCT would work on IXIA
!call MPI_BCAST(phys_list,1,MPI_phys,0,MPI_COMM_WORLD,ierr)
call MPI_PACK_SIZE(1,MPI_REAL8,MPI_COMM_WORLD,IDBL_EXT,ierr)
call MPI_PACK_SIZE(1,MPI_REAL4,MPI_COMM_WORLD,ISGL_EXT,ierr)
call MPI_PACK_SIZE(1,MPI_INTEGER,MPI_COMM_WORLD,INT_EXT,ierr)
call MPI_PACK_SIZE(1,MPI_LOGICAL,MPI_COMM_WORLD,ILOG_EXT,ierr)
call MPI_PACK_SIZE(1,MPI_CHARACTER,MPI_COMM_WORLD,CHAR_EXT,ierr)

bufsize = ( (N_species*3+4)*INT_EXT + &
            (N_species*(9+1))*ISGL_EXT + &
            1*IDBL_EXT + &
            (N_species*2+2)*ILOG_EXT + &
            (80*(N_species+1)+2*6*N_species)*CHAR_EXT)
allocate(buffer(bufsize))

if (my_id .eq. 0) then
  position = 0
  call MPI_PACK(species,                     N_species,MPI_INT      ,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(atomic_mass,                 N_species,MPI_REAL4    ,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(N_particles,                 N_species,MPI_INT      ,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(particle_GC,                 N_species,MPI_LOGICAL  ,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(particle_ion_rec,            N_species,MPI_LOGICAL  ,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(particle_seed,               N_species,MPI_INT      ,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(adas_suffix,               N_species*6,MPI_CHARACTER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(location_accept_function, N_species*80,MPI_CHARACTER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(location_accept_parameters,  N_species*9,  MPI_REAL4,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(particle_initializer,      N_species*6,MPI_CHARACTER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(t_step_particles,                    1,MPI_REAL8    ,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(n_step_particles,                    1,MPI_INT      ,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(nout_particles,                      1,MPI_INT      ,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(write_energies,                      1,MPI_LOGICAL  ,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(write_momenta,                       1,MPI_LOGICAL  ,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(t_particles_begin,                   1,MPI_INT      ,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(t_particles_end,                     1,MPI_INT      ,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(particle_restart_file,              80,MPI_CHARACTER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
endif
call MPI_BCAST(buffer,bufsize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
if (my_id .ne. 0) then
  position = 0
  call MPI_UNPACK(buffer,bufsize,position,species,                      N_species,MPI_INT      ,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,atomic_mass,                  N_species,MPI_REAL4    ,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,N_particles,                  N_species,MPI_INT      ,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,particle_GC,                  N_species,MPI_LOGICAL  ,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,particle_ion_rec,             N_species,MPI_LOGICAL  ,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,particle_seed,                N_species,MPI_INT      ,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,adas_suffix,                N_species*6,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,location_accept_function,  N_species*80,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,location_accept_parameters, N_species*9,MPI_REAL4    ,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,particle_initializer,       N_species*6,MPI_CHARACTER,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,t_step_particles,             1,MPI_REAL8    ,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,n_step_particles,             1,MPI_INT      ,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,nout_particles,               1,MPI_INT      ,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,write_energies,               1,MPI_LOGICAL  ,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,write_energies,               1,MPI_LOGICAL  ,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,t_particles_begin,            1,MPI_INT      ,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,t_particles_end,              1,MPI_INT      ,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,particle_restart_file,       80,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
endif

deallocate(buffer)
end subroutine broadcast_particle_parameters


!#
! Cross product in the JOREK coordinate system (which is left-handed)
! input  :  (\(a_R,a_Z,a_\phi), (b_R,b_Z,b_\phi)\)
! output :  (\(a \times b)_{(R,Z,\phi)}\)
pure function cross_product(a,b)
implicit none
real*8 :: cross_product(3)
real*8, intent(in) :: a(3), b(3)

cross_product(1) = a(2)*b(3) - a(3)*b(2)
cross_product(2) = a(3)*b(1) - a(1)*b(3)
cross_product(3) = a(1)*b(2) - a(2)*b(1)
end function cross_product

!> This function creates a derived MPI type for the particle and returns it
!! If it already exists the old handle is returned
function get_particle_derived_type() result(dtype_out)
  use mpi_mod
  use parameters

  implicit none

  integer               :: ierr, dtype_out
  integer, save         :: dtype
  logical, save         :: dtype_set = .false.

  integer :: len(9) = (/n_dim,3,3,1,1,1,1,1,1/), t(9) = (/ &
    MPI_REAL8,MPI_REAL8,MPI_REAL8,MPI_REAL4,MPI_REAL4, &
    MPI_INTEGER,MPI_INTEGER1,MPI_INTEGER1,MPI_INTEGER1/) ! MPI_INTEGER1 == MPI_LOGICAL1

  integer(kind=MPI_ADDRESS_KIND) :: base, disp(9)
  type(type_particle) :: particle

  dtype_out = dtype
  if (dtype_set) return

  ! Get memory addresses in the type
  call MPI_Get_address(particle,        base,    ierr)
  call MPI_Get_address(particle%st,     disp(1), ierr)
  call MPI_Get_address(particle%x,      disp(2), ierr)
  call MPI_Get_address(particle%v,      disp(3), ierr)
  call MPI_Get_address(particle%mass,   disp(4), ierr)
  call MPI_Get_address(particle%weight, disp(5), ierr)
  call MPI_Get_address(particle%i_elm,  disp(6), ierr)
  call MPI_Get_address(particle%q,      disp(7), ierr)
  call MPI_Get_address(particle%label,  disp(8), ierr)
  call MPI_Get_address(particle%lost,   disp(9), ierr)

  ! Rebase to particle memory beginning
  disp = disp - base

  ! Commit the structured type
  call MPI_Type_create_struct(9, len, disp, t, dtype, ierr)
  if (ierr .ne. 0) write(*,*) "Error creating particle datatype: ", ierr
  call MPI_Type_commit(dtype, ierr)
  if (ierr .ne. 0) write(*,*) "Error committing particle datatype: ", ierr

  ! Set the save bit
  dtype_set = .true.
  dtype_out = dtype
  return
end function get_particle_derived_type

end module mod_particles
