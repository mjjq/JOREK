program static_field_example

use mod_particle_sim, only: particle_sim
use mod_event, only: event
use mod_action, only: stop_action
use mod_prescribed_fields, only: prescribed_fields
implicit none

type(particle_sim) :: sim
type(event), dimension(:), allocatable :: events

integer :: i

allocate(prescribed_fields::sim%fields)
sim%fields%electric_field => E
sim%fields%magnetic_field => B

events = [ &
  event(stop_action(), start=1.d0) & ! Stop the sim after 1 second
]

write(*,*) "not implemented yet"

contains
pure function E(x, t)
  real*8, dimension(3), intent(in) :: x
  real*8, intent(in) :: t
  real*8, dimension(3) :: E
  E = (/0.d0, 0.d0, 0.d0/)
end function E
pure function B(x, t)
  real*8, dimension(3), intent(in) :: x
  real*8, intent(in) :: t
  real*8, dimension(3) :: B
  B = (/0.d0, 0.d0, 1.d0/)
end function B
end program static_field_example
