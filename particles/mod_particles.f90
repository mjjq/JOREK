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
    integer*1 :: species          !< Atomic number
    logical*1 :: lost             !< particle is active or lost
  end type type_particle

  type type_particle_list
    integer                                         :: n_particles  !< the number of particles in the list
    type (type_particle), dimension(:), allocatable :: particle     !< an allocatable list of particles
  end type type_particle_list

  integer, parameter :: n_species      = 9         !< Maximum number of particle species
  !> @name particle initialization parameters
  integer :: species(N_species) = 0
  real*4  :: atomic_mass(N_species)
  integer :: N_particles(N_species) = 0
  logical :: particle_GC(N_species) = .false.
  character(len=80) :: location_accept_function(N_species) = 'location_accept_any'
  real*4  :: location_accept_parameters(1:9,1:N_species) = 0
  !> @name particle timestepping input parameters
  real*8  :: t_step_particles !< the time step for the advance of the particles
  integer :: n_step_particles !< the number of time steps for the particles
  integer :: nout_particles   !< number of particle timestep between each output file
  logical :: write_energies !< Output energies of all the particles every nout
  logical :: write_momenta  !< Output generalized toroidal momentum of all the particles every nout
  integer :: t_particles_begin = -1 !< Start JOREK restart file postprocessing (with special value to use only jorek_restart here)
  integer :: t_particles_end !< end JOREK restart file postprocessing
  character(len=80) :: particle_restart_file !< Particle restart file to read at the beginning of simulation
contains

subroutine initialise_particle_parameters(filename)
implicit none
character(len=*), intent(in) :: filename
integer :: ierr

namelist /in2/ species, atomic_mass, N_particles, particle_GC, &
    location_accept_function, location_accept_parameters, &
    n_step_particles, t_step_particles, nout_particles, &
    write_energies, write_momenta, t_particles_begin, t_particles_end

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

end subroutine initialise_particle_parameters


function cross_product(a,b)
!----------------------------------------
! input  :  (a_R,a_Z,A_phi), (b_R,b_Z,b_phi)
! output :  (a x b)_(R,Z,phi)
!----------------------------------------
implicit none
real*8 :: cross_product(3)
real*8 :: a(3), b(3)

cross_product(1) = a(2)*b(3) - a(3)*b(2)
cross_product(2) = a(3)*b(1) - a(1)*b(3)
cross_product(3) = a(1)*b(2) - a(2)*b(1)

return
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
    MPI_INTEGER,MPI_INTEGER1,MPI_INTEGER1,MPI_LOGICAL/)

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
  call MPI_Get_address(particle%species,disp(8), ierr)
  call MPI_Get_address(particle%lost,   disp(9), ierr)

  ! Rebase to particle memory beginning
  disp = disp - base

  ! Commit the structured type
  call MPI_Type_create_struct(9, len, disp, t, dtype, ierr)
  call MPI_Type_commit(dtype, ierr)

  ! Set the save bit
  dtype_set = .true.
  dtype_out = dtype
  return
end function get_particle_derived_type

end module mod_particles
