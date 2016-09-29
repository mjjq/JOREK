!> Module containing a datatype for simulation parameters
module mod_particle_sim
use mod_fields
use mod_particle_types
implicit none
private
public particle_group, particle_sim

type :: particle_group
  integer :: pusher !< which pusher is responsible for this group
  class(particle_base), dimension(:), allocatable :: particles
  real*8 :: time !< time at which this group is currently
end type particle_group

!> Particle simulation type, containing all variables pertaining to a simulation.
type :: particle_sim
  real*8                                          :: time = 0.d0 !< time of the simulation. Only accurate when in events with sync or at
  !< the start of the simulation
  logical                                         :: stop_now = .false.
  type(particle_group), dimension(:), allocatable :: groups
  class(fields_base), allocatable                 :: fields
end type particle_sim
end module mod_particle_sim
