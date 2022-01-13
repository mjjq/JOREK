!> the mod_particles_parameters contains paramters common to multiple
!> particle modules
module mod_particle_parameters
  implicit none
  !> Parameters -----------------------------------------------------
  integer,parameter :: n_particles_per_tile=25 !< number of particles per OpenMP memory tile
end module mod_particle_parameters

