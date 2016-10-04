!> Calculate the particle trajectories with the Boris method in
!> static JOREK fields, file jorek_restart.h5 or jorek_restart.rst
program ex1_jorek
use particle_tracer
!use jorek_integration
use data_structure
implicit none

real*8             :: timestep = [1d-6] !< timesteps might be slightly changed by fix_timesteps
integer, parameter :: num_events = 2

type(particle_sim) :: sim
type(type_node_list), target :: node_list
type(type_element_list), target :: element_list
type(event), dimension(num_events) :: events
logical, dimension(num_events) :: to_run

integer :: i, j, k, n_steps, i_elm_old, ifail
real*8 ::  next_event_time, t
real*8, dimension(3) :: E, B
real*8, dimension(2) :: rz_old, st_old

call sim%initialize(num_groups=1)

call with(sim, read_jorek_fields(node_list, element_list, basename='jorek_restart'))
call with(sim, init_particles(particle_kinetic_leapfrog(), num_groups=1, uniform_density=1.d0, num_particles=100000))

events = [event(write_action(basename='test', decimal_digits=0), step=1d-4), &
          event(stop_action(), start=1d-3)]
call check_and_fix_timesteps(timesteps, events)
call with(sim, events, mask=abs(events%start - 0.d0) .le. TICK)

do while (.not. sim%stop_now)
  call next_event_at(events, sim%time, to_run, next_event_time)

  do i=1,1
    call sim%groups(i)%before_push
    n_steps = nint((next_event_time - sim%time)/timestep)

    select type (particles => sim%groups(i)%particles)
    type is (particle_kinetic_leapfrog)
      !$omp parallel do default(none) &
      !$omp shared(n_steps, timesteps, node_list, element_list) &
      !$omp private(j, k, t, E, B, rz_old, st_old, i_elm_old)
      do j=1,size(particles)
        if (particles(j)%lost) cycle
        do k=1,n_steps
          t = sim%time + k*timestep
          !call sim%fields%at(particles(j), E, B, t)
          rz_old = particles(j)%x(1:2)
          st_old = particles(j)%st
          i_elm_old = particles(j)%i_elm
          call boris_push_cylindrical(particles(j), E, B, timestep))
          call find_RZ_nearby(node_list, element_list, particles(j)%x(1:2), rz_old, &
              st_old, particles(j)%st, i_elm_old, particles(j)%i_elm, ifail)
          if (ifail .eq. -1) then
            particles(j)%lost = .true.
            exit ! the time-loop
          end if

          ! Check new charge state every 3 iterations
          if (mod(k, 3) .eq. 0) then
            particles(j)%q = new_charge(particles(j)%q, ad, electron_density, electron_temperature, timestep)
          end if
          ! Check collisions every 100 iterations
          if (mod(k, 100) .eq. 0) then
          end if


        end do
      end do
      !$omp end parallel do
    end select
    call sim%groups(i)%after_push(i, next_event_time)
  end do

  sim%time = next_event_time
  call with(sim, events, to_run)
end do

call sim%finalize
end program ex1_jorek
