!> Module containing a datatype for simulation parameters
module mod_particle_sim
use mod_particle_types
use mod_fields
use mod_openadas
implicit none
private
public particle_group, particle_sim

!> A group of particles, implemented as an allocatable array.
!> It must contain particles of the same species (charge number).
type :: particle_group
  integer :: Z !< Atomic number of al particles in the group (-1 for electrons, 0 for fieldline-following)
  real*8  :: mass !< Mass of all the particles in the group
  type(ADF11_all) :: ad !< OPEN-ADAS datafiles for this species
  class(particle_base), dimension(:), allocatable :: particles
end type particle_group

!> Particle simulation type, containing all variables pertaining to a simulation.
type :: particle_sim
  real*8                                          :: time = 0.d0 !< time of the simulation. Only accurate when in events with sync or at
  !< the start of the simulation
  class(fields_base), allocatable                 :: fields
  logical                                         :: stop_now = .false.
  type(particle_group), dimension(:), allocatable :: groups
  !< MPI settings
  integer :: my_id
  integer :: n_cpu
  real*8 :: wtime_start !< Clock time at the start of the program
contains
  procedure :: finalize
  procedure :: initialize
end type particle_sim

contains
!> Actions to perform when setting up a simulation
subroutine initialize(sim, num_groups)
  use mpi
  use mod_parameters, only: n_tor, n_period
  use phys_module, only: mode
  class(particle_sim), intent(inout) :: sim !< why is this class() and not type()?
  integer, intent(in) :: num_groups
  integer :: provided, ierr, i_tor
  character(len=MPI_MAX_PROCESSOR_NAME) :: name
  integer :: resultlength

  call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  if (ierr .ne. 0) write(*,*) "Error ", ierr, " in MPI_Init_thread"
  call MPI_COMM_RANK(MPI_COMM_WORLD, sim%my_id, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, sim%n_cpu, ierr)
  ! Synchronize clocks
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  sim%wtime_start = MPI_Wtime() ! accurate up to the network latency (fine for times measured in seconds)

  if (provided .ne. MPI_THREAD_MULTIPLE .and. sim%my_id .eq. 0) write(*,*) "WARNING: provided(", provided, ") != MPI_THREAD_MULTIPLE"
  allocate(sim%groups(num_groups))
  call MPI_GET_PROCESSOR_NAME(name,resultlength,ierr)
  write(*,'(A,I5,2A)') '#MPI id, ProcessorName ', sim%my_id, ': ', name

  if (sim%my_id .eq. 0) call print_version

  ! Initialise mode numbers
  do i_tor=1, n_tor
    mode(i_tor) = + int(i_tor / 2) * n_period
    if (sim%my_id .eq. 0) write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
  enddo

  ! Initialise parameters
  call initialise_and_broadcast_parameters(sim%my_id, "__NO_FILENAME__")

  ! Initialise the gaussian points at basis functions
  call initialise_basis
end subroutine

!> Actions to perform when stopping the simulation.
subroutine finalize(sim)
  class(particle_sim), intent(in) :: sim !< why is this class() and not type()?
  integer :: ierr
  if (sim%stop_now) then
    write(*,"(A,g14.6,A)") "INFO: Stop requested at ", sim%time, " , exiting"
  else
    write(*,"(A,g14.6,A)") "INFO: End of events at ", sim%time, " , exiting"
  end if
  call MPI_Finalize(ierr)
end subroutine
end module mod_particle_sim
