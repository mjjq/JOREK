!> Calculate the particle trajectories with the Boris method in
!> static JOREK fields, file jorek_restart.h5 or jorek_restart.rst
program ex1_jorek
use particle_tracer
!use jorek_integration
use data_structure
implicit none

real*8             :: timestep = [1d-6] !< timesteps might be slightly changed by fix_timesteps

type(particle_sim) :: sim
type(type_node_list), target :: node_list
type(type_element_list), target :: element_list
type(event), dimension(:), allocatable :: events

integer :: i, j, k, n_steps, i_elm_old, ifail
real*8 ::  target_time, t
real*8, dimension(3) :: E, B
real*8, dimension(2) :: rz_old, st_old

call sim%initialize(num_groups=1)

call with(sim, read_jorek_fields(node_list, element_list, basename='jorek_restart'))
call with(sim, init_particles(particle_kinetic_leapfrog(), num_groups=1, uniform_density=1.d0, num_particles=100000))

events = [event(write_action(basename='test', decimal_digits=0), step=1d-4), &
          event(stop_action(), start=1d-3)]
call check_and_fix_timesteps(timesteps, events)
call with(sim, events, at=0.d0)

do while (.not. sim%stop_now)
  call next_event_at(events, sim%time, target_time)

  do i=1,1 ! loop over groups
    n_steps = nint((next_event_time - sim%time)/timestep)

    select type (particles => sim%groups(i)%particles)
    type is (particle_kinetic_leapfrog)
      !$omp parallel do default(private) &
      !$omp shared(n_steps, timesteps, i, node_list, element_list)
      do j=1,size(particles)
        if (particles(j)%lost) cycle
        do k=1,n_steps
          t = sim%time + k*timestep
          !call sim%fields%at(particles(j), E, B, t)
          rz_old = particles(j)%x(1:2)
          st_old = particles(j)%st
          i_elm_old = particles(j)%i_elm
          call boris_push_cylindrical(particles(j), E, B, timestep)
          call find_RZ_nearby(node_list, element_list, particles(j)%x(1:2), rz_old, &
              st_old, particles(j)%st, i_elm_old, particles(j)%i_elm, ifail)
          if (ifail .eq. -1) then
            particles(j)%lost = .true.
            exit ! the time-loop
          end if

          ! Check new charge state every 3 iterations
          if (mod(k, 3) .eq. 0) then
            particles(j)%q = new_charge(particles(j)%q, ad, electron_density, electron_temperature, timesteps(i))
          end if
          ! Check collisions every 100 iterations
          if (mod(k, 100) .eq. 0) then
          end if
        end do ! steps
      end do ! particles
      !$omp end parallel do
    end select
  end do ! groups
  sim%time = target_time
  call with(sim, events, at=sim%time)
end do
call sim%finalize
end program ex1_jorek
