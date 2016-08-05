!> Module to formalize performing an action every now and then
!> Particle events extend normal events
!> istep is now in units of the particle timestep
module mod_particle_events
use mod_events
implicit none

type, extends(type_event) :: type_particle_event
  procedure(particle_action), pointer, public :: action => NULL()
contains
  procedure, public  :: check_and_run => check_and_run_particle_action
end type type_particle_event

interface
  subroutine particle_action(this, node_list, element_list, particle_list)
    use data_structure
    use mod_particles
    implicit none

    class(type_particle_event), intent(inout) :: this
    type(type_node_list), intent(inout)     :: node_list
    type(type_element_list), intent(inout)  :: element_list
    type(type_particle_list), intent(inout) :: particle_list
  end subroutine particle_action
end interface

contains
  subroutine check_and_run_particle_action(this, index_now, t_now, node_list, element_list, particle_list)
    use data_structure
    use mod_particles
    implicit none

    class(type_particle_event), intent(inout) :: this
    integer, intent(in) :: index_now
    real*8 , intent(in) :: t_now
    type(type_node_list), intent(inout)       :: node_list
    type(type_element_list), intent(inout)    :: element_list
    type(type_particle_list), intent(inout)   :: particle_list

    real*8 :: t_start, t_end

    if (can_start(this, index_now, t_now)) then
      call cpu_time(t_start)
      call this%action(node_list, element_list, particle_list)
      call cpu_time(t_end)
      this%walltime = this%walltime + (t_end - t_start)
    end if
  end subroutine check_and_run_particle_action
end module mod_particle_events
