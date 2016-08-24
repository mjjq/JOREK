!> Module containing a datatype for simulation parameters
module mod_particle_sim
use mod_fields
use mod_particle_group
use mod_particle_io
implicit none

!> Particle simulation type, containing all variables pertaining to a simulation.
type particle_sim
  real*8                                          :: time = 0.d0 !< time of the simulation. Only accurate when in events with sync or at
  !< the start of the simulation
  type(particle_hdf5_io)                          :: io
  type(particle_group), dimension(:), allocatable :: groups
  class(fields_base), allocatable                 :: fields
end type particle_sim
end module mod_particle_sim
