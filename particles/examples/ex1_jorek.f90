!> Calculate the particle trajectories with the Boris method in
!> static JOREK fields, file jorek_restart.h5 or jorek_restart.rst
program ex1_jorek
use particle_tracer
implicit none

real*8 :: timesteps(1) = [1d-7] !< seconds

type(adf11_all) :: adas
type(coronal)   :: cor
integer :: i, j, k, n_steps, i_elm_old, ifail
real*8 :: target_time, t
real*8 :: E(3), B(3), rz_old(2), st_old(2), psi, U

! Start up MPI, jorek
call sim%initialize(num_groups=1)

! Set up the field reader
fields = read_jorek_fields_interp_linear(basename='jorek_restart', i=-1)
call with(sim, fields)

! Set up particles
allocate(particle_kinetic_leapfrog::sim%groups(1)%particles(100000))
sim%groups(1)%Z    = 74
sim%groups(1)%mass = 183.84 !< atomic mass units

! Prepare the coronal equilibrium
adas = read_adf11('50_w')
cor  = coronal(adas)

! Distribute particles uniformly throughout the domain
call seed_positions(sim%groups(1)%particles, &
    fields%node_list, fields%element_list, pcg32_rng())
call adjust_particle_weights(sim%groups(1)%particles, num_atoms_total=1d23)
call set_velocity_from_T(sim%groups(1)%particles, sim%groups(1)%mass, &
    fields%node_list, fields%element_list, pcg32_rng(), cor, v_par=.true.)

events = [event(write_action(basename='test'),   step=1d-4), &
          !event(diag_print_kinetic_energy(),     step=1d-6), &
          event(project_to_vtk(fields%node_list, fields%element_list, smoothing=1d-3, basename='proj'), step=1d-5), &
          event(stop_action(), start=1d-3)]
call check_and_fix_timesteps(timesteps, events)
call with(sim, events, at=0.d0)

do while (.not. sim%stop_now)
  target_time = next_event_at(sim, events)

  do i=1,1 ! loop over groups
    n_steps = nint((target_time - sim%time)/timesteps(i))

    select type (particles => sim%groups(i)%particles)
    type is (particle_kinetic_leapfrog)
      !$omp parallel do default(private) &
      !$omp shared(sim, n_steps, timesteps, i, fields)
      do j=1,size(particles,1)
        do k=1,n_steps
          if (particles(j)%i_elm .eq. 0) exit
          t = sim%time + k*timesteps(i)
          call EM_fields_interp_linear(fields, particles(j)%i_elm, &
              particles(j)%st, particles(j)%x(3), E, B, psi, U)
          rz_old    = particles(j)%x(1:2)
          st_old    = particles(j)%st
          i_elm_old = particles(j)%i_elm

          call boris_push_cylindrical(particles(j), sim%groups(i)%mass, E, B, timesteps(i))
          call find_RZ_nearby(fields%node_list, fields%element_list, particles(j)%x(1:2), rz_old, &
              st_old, particles(j)%st, i_elm_old, particles(j)%i_elm, ifail)
          if (particles(j)%i_elm .eq. 0) write(*,*) "particle ", j, " lost at step ", k
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
