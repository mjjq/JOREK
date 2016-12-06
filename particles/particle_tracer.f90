module particle_tracer
! Base section
use constants
use mod_particle_sim
use mod_particle_types
use mod_event
use mod_action
use mod_initialise_particles

! IO
use mod_io_actions

! Pushers
use mod_boris

! Fields
use mod_fields
use mod_fields_linear

! ADAS
use mod_openadas
use mod_coronal

! RNGs
use mod_pcg32_rng
use mod_sobseq_rng

! Diagnostics
use mod_diag_print_kinetic_energy
use mod_project_particles

! JOREK
use data_structure

! Default variables
implicit none
type(particle_sim) :: sim
type(event), dimension(:), allocatable :: events
end module particle_tracer
