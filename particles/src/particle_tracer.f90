module particle_tracer
  ! Base section
  use mod_constants
  use mod_particle_sim
  use mod_particle_types
  use mod_event
  use mod_action
  use mod_main_loop
  use mod_pusher
  use mod_hook

  ! IO
  use mod_io_actions

  ! Fields section
  use mod_prescribed_fields

  ! Pushers
  use mod_boris
  implicit none
end module particle_tracer
