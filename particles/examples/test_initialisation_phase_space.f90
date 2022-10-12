program test_initialisation_phase_space
!> program used for testing the initialisation method: 
!>   initialise_particles_in_phase_space 
!> contained in mod_initialise_particles
!> the idea is to provide a reject_uniform function
!> which reject a uniformly distributed random sample
!> between [0 1] as a function of a user defined 
!> 6D distribution. The distribution function as to be
!> in R,Z,phi coordinates for the physical space, 
!> momentum, pitch and gyro angles coordinates for the
!> momentum space and the distribution of electric charges
!> be aware that the electric charge should be a double 
!use mod_particle_sim
!use mod_particle_types,       only: particle_base,particle_kinetic_relativistic 
!use mod_fields_linear,        only: jorek_field_interp_linear,read_jorek_fields_interp_linear
!use mod_initialise_particles, only: initialise_particles_in_phase_space
use constants, only: TWOPI,PI
use mod_random_seed
use particle_tracer
implicit none

!> Variable declarations --------------------------------------------------------------------
integer,parameter            :: n_groups=1
type(sobseq_rng)             :: sob_rngs
type(event)                  :: field_reader
integer                      :: n_particles
real*8                       :: start_time,mass
real*8,dimension(2)          :: Rbound,Zbound,Phibound,Ekinbound,Pitchbound,Chibound,Chargebound
character(len=:),allocatable :: jorek_filename

!> MPI and groups initialisation ------------------------------------------------------------
call sim%initialize(num_groups=n_groups)

!> Define inputs ----------------------------------------------------------------------------
n_particles = 1000
start_time  = 0.d0
mass        = 5.48579909065d-4 !< electron mass in AMU
Rbound      = [0.d0,9.99d2]
Zbound      = [-9.99d2,9.99d2]
Phibound    = 5.d-1*[PI,3.d0*PI]
Ekinbound   = 2d7
Pitchbound  = [PI-2.95d-1,PI]
Chibound    = [0.d0,TWOPI]
Chargebound = -1.d0
allocate(character(len=25)::jorek_filename)
jorek_filename = 'jorek_equilibrium' 


!> Clean-up
call sim%finalize()

contains

!>  Phase space distributions for initialisation
function reject_uniform(n_x,x,st,i_elm,rng,fields)
  use mod_fields, only: fields_base
  implicit none
  !> constants
  real*8,parameter                 :: threshold=5.d-1
  !> inputs:
  integer,intent(in)               :: n_x,i_elm
  real*8,intent(in)                :: rng
  real*8,dimension(n_x),intent(in) :: x
  real*8,dimension(2),intent(in)   :: st
  class(fields_base),intent(in)    :: fields
  !> outputs:
  logical                          :: reject_uniform

  !> check if the particle is valid
  reject_uniform = .true.
  if(i_elm.le.0) return
  if((st(1).lt.0.d0).or.(st(1).gt.1.d0)) return
  if((st(2).lt.0.d0).or.(st(2).gt.1.d0)) return
  !> rject of accept solution
  if(rng.le.threshold) reject_uniform = .false.

end function reject_uniform

end program test_initialisation_phase_space
