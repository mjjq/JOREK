!> Particle pusher module with the Boris scheme
module mod_particle_boris
  use mod_particle_base

  type, extends(particle_base) :: particle_boris
    real*8, dimension(3) :: v !< Velocity in real space
  end type particle_boris
end module mod_particle_boris
