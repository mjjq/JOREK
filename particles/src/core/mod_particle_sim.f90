!> Module containing a datatype for simulation parameters and routines for
!> reading this from an input namelist
module mod_particle_sim
use mod_fields
use mod_particle_group
use mod_particle_io
implicit none

!> Particle simulation type, containing all variables pertaining to a simulation.
type particle_sim
  type(particle_hdf5_io)                          :: io
  type(particle_group), dimension(:), allocatable :: group
  class(fields_base), allocatable                 :: fields
end type particle_sim
end module mod_particle_sim
