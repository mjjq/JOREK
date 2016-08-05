!> Particle pusher module with the Boris scheme
module mod_boris
  use mod_particle_boris
  use mod_particle_io_boris
  use mod_pusher_boris

  type, extends(pusher_params_base) :: pusher_params_boris
  end type pusher_params_boris
end module mod_boris
