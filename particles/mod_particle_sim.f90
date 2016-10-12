!> Module containing a datatype for simulation parameters
module mod_particle_sim
use mod_particle_types
implicit none
!private
public particle_group, particle_sim

!> A group of particles, implemented as an allocatable array
type :: particle_group
  class(particle_base), dimension(:), allocatable :: particles
end type particle_group

!> Particle simulation type, containing all variables pertaining to a simulation.
type :: particle_sim
  real*8                                          :: time = 0.d0 !< time of the simulation. Only accurate when in events with sync or at
  !< the start of the simulation
  logical                                         :: stop_now = .false.
  type(particle_group), dimension(:), allocatable :: groups
contains
  procedure :: finalize
  procedure :: initialize
end type particle_sim

contains
!> Actions to perform when setting up a simulation
subroutine initialize(sim, num_groups)
#include "version.h"
  use mpi
  class(particle_sim), intent(inout) :: sim !< why is this class() and not type()?
  integer, intent(in) :: num_groups
  integer :: provided, ierr, my_id
  call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  if (ierr .ne. 0) write(*,*) "Error ", ierr, " in MPI_Init_thread"
  if (provided .ne. MPI_THREAD_MULTIPLE) write(*,*) "WARNING: provided(", provided, ") != MPI_THREAD_MULTIPLE"
  allocate(sim%groups(num_groups))
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

  111 format(2x,a,': ',a)
  if (my_id .eq. 0) then
    write(*,*) ' ', trim(adjustl(RCS_VERSION))
    write(*,111) 'compile_time        ', trim(adjustl(compile_time))
    write(*,111) 'compile_user        ', trim(adjustl(compile_user))
    write(*,111) 'compile_machine     ', trim(adjustl(compile_machine))
    write(*,111) 'compile_dir         ', trim(adjustl(compile_dir))
    write(*,111) 'compile_command     ', trim(adjustl(compile_command))
    write(*,111) 'compile_flags       ', trim(adjustl(compile_flags))
    write(*,111) 'compile_includes    ', trim(adjustl(compile_includes))
    write(*,111) 'compile_defines     ', trim(adjustl(compile_defines))
    write(*,111) 'compile_libs        ', trim(adjustl(compile_libs))
    write(*,111) 'compile_modules     ', trim(adjustl(compile_modules))
  end if
end subroutine

!> Actions to perform when stopping the simulation.
subroutine finalize(sim)
  class(particle_sim), intent(in) :: sim !< why is this class() and not type()?
  integer :: ierr
  if (sim%stop_now) then
    write(*,"(A,g12.6,A)") "INFO: Stop requested at ", sim%time, " , exiting"
  else
    write(*,"(A,g12.6,A)") "INFO: End of events at ", sim%time, " , exiting"
  end if
  call MPI_Finalize(ierr)
end subroutine
end module mod_particle_sim
