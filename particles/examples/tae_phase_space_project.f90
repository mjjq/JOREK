!> Testing the coupling of the projections of particles to JOREK
program tae_loop

  use particle_tracer
  use mod_particle_diagnostics
  use mpi
  use mod_atomic_elements
  use mod_particle_io
  use mod_event
  use mod_project_particles
  use mod_particle_loop
  use mod_jorek_timestepping
  use mod_random_seed
  use mod_interp, only: mode_moivre, interp_RZ, interp_0
  use mod_basisfunctions
  use nodes_elements
  use phys_module, only: tstep, restart, t_start, restart_particles
  use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, xcase, xpoint,F0
  use phys_module, only: n_particles, nstep_particles, nsubstep_particles, tstep_particles
  use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0
  use phys_module, only: n_mode_families

  use constants,   only: MU_ZERO, MASS_PROTON, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG
  use mod_math_operators, only: cross_product
  use mod_particle_sputtering, only: particle_sputter, sample_fluid_particle_energy
  use mod_projection_functions, only: proj_f_combined_density, &
      proj_f_combined_energy, proj_f_combined_par_momentum
  use mod_edge_domain
  use mod_edge_elements
  use data_structure, only: type_bnd_element_list, type_bnd_node_list
  use equil_info
  use mod_boundary, only: boundary_from_grid
  use mod_phase_space_project

!$ use omp_lib

  implicit none

  type(event)                                       :: fieldreader, partreader, partwriter
  !type(adf11_all)                                   :: adas
  type(pcg32_rng), dimension(:), allocatable        :: rng
  type(count_action)                                :: counter
  type(projection), target                          :: jorek_feedback, project_density, project_current
  type(jorek_timestep_action), target               :: jorek_stepper
  type(type_edge_domain), allocatable, dimension(:) :: edge_domains
  type(edge_elements)                               :: D_edge
  type(write_particle_diagnostics)                  :: diag
  type(phase_space_projection)                      :: test_phase
  real*8, parameter  :: binding_energy = 2.18d-18 ! ionization energy of a hydrogen atom [J] (= 13.6 eV)
  real*8    :: target_time
  real*8    :: physical_particles, weight
  real*8    :: oldtime, step_rest_time, particle_step_time, particle_start_time, diag_time
  real*8    :: rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, tstep_si, timesteps
  real*8    :: v_kin_temp, E(3), B(3), psi, U, B_norm(3)
  real*8    :: rescale_coef, T_axis(1), E_axis, E_hot, rho_part, v2, tstart_jorek
!$ real*8 :: w0, w1, mmm(3)

  real*8 :: test_x(3), test_E(3), test_B(3), test_psi, test_U
  !real*8 :: b_norm_r, b_norm_z, b_norm_phi

  integer   :: n_particles_local,n_reflect,ifail, ino
  integer   :: i, j, k, l, m, n_steps, i_elm_old, i_diagno
  integer   :: seed, i_rng, n_stream,i_tor,ierr

  ! Start up MPI, jorek
  call sim%initialize(num_groups=1)

  rho_part    = 1.195d19 !(corrected value to obtain density=1.441e17 (as in benchmark, for original profile with toroidal flux)
  n_particles_local = int(n_particles/sim%n_cpu)
  timesteps         = tstep_particles
  write(*,*) tstep
  ! Set up the field reader
  fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1,mode_divisor=100))
  call with(sim, fieldreader)

  write(*,*) 'main : t_start = ',t_start


  if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)

  call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)

  call update_equil_state(sim%my_id,sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase)

  n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
  rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
  t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
  write(*,*) 'Alfven =', F0/10.d0/t_norm
  tstep_si  = tstep * t_norm
  n_steps   = floor(tstep_si / timesteps)
  timesteps = tstep_si / n_steps
  n_steps   = tstep_si / timesteps

  if (sim%my_id .eq.0) then
    write(*,*) ' adapt time step to be multiple of jorek time step'
    write(*,*) "tstep = ", tstep_si, n_steps, timesteps
    write(*,*) "check :", n_steps, tstep_si - n_steps*timesteps

    i_diagno =  sim%fields%node_list%n_nodes / 3
    write(*,'(A,6f8.4)') ' probe at : ',sim%fields%node_list%node(i_diagno)%x(1,1,1:2)
    open(111,file='diagno.txt')
  endif
  write(*,*) "until start phase" , sim%my_id
  if (.not. restart_particles) then
    ! Set up particles
    sim%groups(1)%Z    = 1
    sim%groups(1)%mass = atomic_weights(-2) !< atomic mass units
    !sim%groups(1)%ad   = adas

    allocate(particle_kinetic_leapfrog::sim%groups(1)%particles(n_particles_local))

    !< If projecting phase space, it is vital to use a by construction phi independent initial distribtuion, i.e. phi planes.
    call initialise_particles_H_mu_psi_phiplanes(sim%groups(1)%particles, sim%fields, pcg32_rng(),sim%groups(1)%mass, &
         uniform_space=.true., uniform_space_rej_f=f_toroidal_flux, &
         uniform_space_rej_vars=[1], charge = 1, T_maxwell = 4d5,n_phi_planes_in=12)


    call adjust_particle_weights(sim%groups(1)%particles, rho_part)
    if (sim%my_id .eq. 0) write(*,*) "Particle density was adjusted to:", rho_part, sim%groups(1)%particles(1)%weight
    diag = write_particle_diagnostics(filename='diag.h5', append=.false.)
    call with(sim,diag)
    select type (p => sim%groups(1)%particles)
      type is (particle_kinetic_leapfrog)

        call boris_all_initial_half_step_backwards_RZPhi(p, sim%groups(1)%mass, sim%fields, sim%time, timesteps)

    end select


  else  ! restarting particles

    if (sim%my_id .eq. 0) write(*,*) 'restarting particles: reading part_restart.h5'

    deallocate(sim%groups)
    allocate(sim%groups(0))

    partreader = event(read_action(filename='part_restart.h5'))
    call with(sim, partreader)

  endif


  jorek_feedback = new_projection(sim%fields%node_list, sim%fields%element_list, &
      filter    = filter_perp, filter_hyper    = filter_hyper, filter_parallel    = filter_par, &
      filter_n0 = filter_perp, filter_hyper_n0 = filter_hyper, filter_parallel_n0 = filter_par_n0, &
      calc_integrals=.false., to_vtk=.true., to_h5 = .false., basename='projections')

  ! A few examples of phase space projections.

  ! Simple 2D density histogram
  !test_phase = new_phase_space_projection(sim,ndim=2,res=[200,210],start=[9.d0,-1.d0],end=[11.d0,1.d0],f_proj=proj_f(proj_one, group = 1),&
  !                                        f_grids= [proj_f(proj_R,group =1 ),proj_f(proj_Z,group=1)],bandwidths=[0.2,0.2])
  ! Power versus minor radius (f_proj is not relevant, as we fill the value arrays manually by averaging over particle orbits.)
  ! If you want to use this, you'll have to use a restart file with some mode structure (i.e. a linear-phase TAE mode from previous simulations)
  test_phase = new_phase_space_projection(sim,ndim=2,res=[150,150],start=[-5.d-19,-0.1d0],end=[5.d-19,1500.d0],f_proj=proj_f(proj_one, group = 1),&
                                         f_grids= [proj_f(proj_R,group = 1 ),proj_f(proj_Z,group = 1 )],bandwidths=[1d-19,100d0])

 
  ! Initial density of mu
  !test_phase = new_phase_space_projection(sim,ndim=1,res=[200],start=[0.d0],end=[1.d6],f_proj=proj_f(proj_one, group = 1),&
  !                                         f_grids= [proj_f(proj_mu,group = 1 )],bandwidths=[0.05d6])

  ! Example calling whole sim projection and outputting
  !call with(sim,test_phase)
  !call output_phase_project(test_phase)
  ! This only works w/ nearest neighbour projection (not shaped kernels), so be careful! Manual filling (as in particle loop)
  ! works for both.
  
  ! Full tensor + density (for density flattening was the idea)
  allocate(jorek_feedback%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 7))

  jorek_feedback%rhs = 0.d0

  aux_node_list => jorek_feedback%node_list

  ! For proper timestepping, the projections need to be defined before the jorek timestepper
  jorek_stepper = new_jorek_timestep_action(jorek_feedback%node_list)

  if (restart) then
    tstart_jorek = sim%time + tstep_si
  else
    tstart_jorek = sim%time
  endif

  diag_time = timesteps
  events = [new_event_ptr(jorek_feedback,  start = tstart_jorek),  &
      new_event_ptr(jorek_stepper,   start = tstart_jorek), &
      event(stop_action(), start=1d12)  ]

  jorek_stepper%extra_event => events(1)

  ! Set up random numbers for ionisation probability
  seed = random_seed()
  n_stream = 1
!$ n_stream = omp_get_max_threads()
  allocate(rng(n_stream))
  do i=1,n_stream
    call rng(i)%initialize(1, seed, n_stream, i)
  end do

  ! Call events at sim%time once to help event scheduler, before entering particle loop
  step_rest_time = 0.d0
  ino = 0
  if (.not. restart_particles) call with(sim, events, at=sim%time)

  do while (.not. sim%stop_now)

    target_time = next_event_at(sim, events)
    particle_start_time = (sim%time - step_rest_time)
    particle_step_time  = target_time - particle_start_time
    n_steps             = particle_step_time/timesteps
    step_rest_time      = particle_step_time - real(n_steps,8) * timesteps

    if (sim%my_id .eq. 0) then
      if (n_steps < 10) write(*,*) 'low n_steps,', n_steps
      write(*,*) 'Time difference between particles and jorek: ', step_rest_time
      write(*,*) "PARTICLE : target time         : ",target_time
      write(*,*) "PARTICLE : timesteps           : ",timesteps
      write(*,*) "PARTICLE : sim%time            : ",sim%time
      write(*,*) "PARTICLE : particle_start_time : ",particle_start_time
      write(*,*) "PARTICLE : particle_step_time  : ",particle_step_time
      write(*,*) "PARTICLE : n_steps             : ",n_steps
      write(*,*) "PARTICLE : step_rest_time      : ",step_rest_time
    endif
    test_phase%values = 0.d0
    call loop_particle_kinetic_local(sim, jorek_feedback, rng, timesteps, n_steps , particle_start_time,test_phase)
    call output_phase_project(test_phase,ino)
    ino = ino +1
    ! Output 2D pressure projection as well for completeness (& verify the initialization worked, to_vtk = .true.)
    sim%time = target_time
    call with(sim,events, at=sim%time)

  end do



  partwriter = event(write_action(filename='part_restart.h5'))
  call with(sim, partwriter)

  call sim%finalize

  if (sim%my_id == 0) close(111)

contains


subroutine loop_particle_kinetic_local(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time,test_phase)
  use mod_project_particles
  use mod_random_seed
  use mod_interp, only: mode_moivre
  use mod_basisfunctions
  use mod_particle_types, only: copy_particle_kinetic_leapfrog

  implicit none

  class(particle_sim), target, intent(inout)                :: sim
  type(projection), target, intent(inout)                   :: jorek_feedback
  type(phase_space_projection), target, intent(inout)       :: test_phase
  type(count_action)                                        :: counter
  type(pcg32_rng), dimension(:), allocatable, intent(inout) :: rng
  type(particle_kinetic_leapfrog)                           :: particle_tmp

  real*8, intent(in)     :: timesteps, particle_start_time
  real*8    :: n_norm, rho_norm, t_norm, v_norm, E_norm, M_norm
  real*8    :: t, E(3), B(3), psi, U, n_e, T_e, rz_old(2), st_old(2),rzp_old(3),vcart_old(3),vcart_new(3)
  real*8    :: v_temp(3), T_eV, K_eV, v_kin_temp, B_norm(3), v,E_diff
  real*8    :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac, HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4), E_tot,E_tot_red, E_after, E_after_red
  real*8    :: b_norm_r, b_norm_z, b_norm_phi, vr_tilde,vz_tilde,v_par, p_par, p_perp, p_atrop,val_tmp(test_phase%totsupport),val_tmp2
  real*8    :: Zephi
  real*8    :: Zephi_red
!$ real*8 :: w0, w1, mmm(3)


  integer, intent(in)   :: n_steps
  integer   :: i, j, k, l, m, i_elm_old, i_elm,i_phase
  integer   :: seed, i_rng, n_stream, ierr, nthreads
  integer   :: i_tor, index_lm, i_elm_temp
  integer   :: n_particles, ifail, index_phase_tmp(test_phase%totsupport),index_phase_tmp2
  real*8,allocatable :: feedback_rhs(:,:,:,:,:)
  real*8, allocatable :: phase_proj(:)

!$ w0 = omp_get_wtime()

  n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
  rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
  t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
  v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
  E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
  M_norm   = rho_norm * v_norm                                    ! momentum normalisation

  jorek_feedback%rhs_gather_time = jorek_feedback%rhs_gather_time + n_steps * timesteps
  allocate(feedback_rhs,source=jorek_feedback%rhs)
  allocate(phase_proj,source=test_phase%values)

  jorek_feedback%rhs = 0.d0
  feedback_rhs       = 0.d0

  call with(sim, counter)
  E_tot = 0.d0
  E_tot_red=0.d0
  E_after =0.d0
  E_after_red = 0.d0
  Zephi=0.d0
  select type (particles => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)
       do j=1,size(particles,1)
          if (particles(j)%i_elm > 0) then
             call sim%fields%calc_EBpsiU(t, particles(j)%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
             Zephi = Zephi + particles(j)%q*particles(j)%weight*el_chg*U*F0
          endif
        E_tot = E_tot + 0.5d0*particles(j)%weight*sim%groups(1)%mass*mass_proton*dot_product(particles(j)%v, particles(j)%v)
      enddo
  end select
  call MPI_REDUCE(E_tot,E_tot_red,1,MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(Zephi, Zephi_red,1,MPI_REAL8, MPI_SUM,0,MPI_COMM_WORLD,ierr)
  write(*,*) "On process", sim%my_id, "tot energy: ", E_tot
  write(*,*) "On process", sim%my_id, "tot ezphi: ", Zephi
  if(sim%my_id .eq. 0 ) write(*,*) "Total energy:", E_tot_red
  if(sim%my_id .eq. 0) write(*,*) "Total ezphi:", Zephi_red

  select type (particles => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)
#ifdef __GFORTRAN__
      !$omp parallel do default(shared) &
#else
      !$omp parallel do default(none) &
      !$omp shared(sim, particles, n_steps, timesteps, rng, particle_start_time,        &
      !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm,                           &
      !$omp jorek_feedback, CENTRAL_DENSITY, CENTRAL_MASS,test_phase)                              &
#endif
      !$omp private(particle_tmp, i_rng, i,j,k,l,m, t, E, B, psi, U, rz_old, st_old, index_phase_tmp,  val_tmp, i_phase, index_phase_tmp2,val_tmp2, &
      !$omp i_elm_old, i_elm, n_e, T_e, b_norm_r, b_norm_z,b_norm_phi, vr_tilde, vz_tilde,v_par, p_par, p_perp, p_atrop,&
      !$omp R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, HH, HH_s, HH_t, HZ, index_lm, ifail, v,rzp_old,vcart_old,vcart_new,E_diff) &
      !$omp schedule(dynamic,10) &
      !$omp reduction(+:feedback_rhs)&
      !$omp reduction(+:phase_proj)
      do j=1,size(particles,1)

        call copy_particle_kinetic_leapfrog(particles(j),particle_tmp)

        !      i_rng = 1
!$      i_rng = omp_get_thread_num()+1

        do k=1,n_steps

          if (particle_tmp%i_elm .le. 0) exit

          t = particle_start_time + (k-1)*timesteps

          call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)

          rz_old    = particle_tmp%x(1:2)
          rzp_old   = particle_tmp%x
          st_old    = particle_tmp%st
          i_elm_old = particle_tmp%i_elm

          if (particle_tmp%i_elm .gt. 0) then
            ! Do phase space projection before pushing
            ! For completeness, this is the nearest neighbour implementation.
            ! call calc_index_val_phaseproj(test_phase,particle_tmp,index_phase_tmp2,sim)
            ! if(index_phase_tmp2 > 0) then
            !  phase_proj(index_phase_tmp2)=phase_proj(index_phase_tmp2)+1.d0/n_steps*particle_tmp%weight*dot_product(particle_tmp%v,E)!
            ! endif

            ! Push the particle and determine its new location.
            E_diff = particle_tmp%weight*0.5d0*sim%groups(1)%mass*MASS_PROTON*dot_product(particle_tmp%v,particle_tmp%v)
            call boris_push_cylindrical(particle_tmp, sim%groups(1)%mass, E, B, timesteps)

            call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
                     particle_tmp%x(1), particle_tmp%x(2), particle_tmp%st(1), particle_tmp%st(2), particle_tmp%i_elm, ifail)
            E_diff =-E_diff+ particle_tmp%weight*0.5d0*sim%groups(1)%mass*MASS_PROTON*dot_product(particle_tmp%v,particle_tmp%v)
            !Particle perpendicular & parallel pressure averaging
            !Normalized b vector
            b_norm_r= B(1)/sqrt(B(1)**2+B(2)**2+B(3)**2)
            b_norm_z= B(2)/sqrt(B(1)**2+B(2)**2+B(3)**2)
            b_norm_phi= B(3)/sqrt(B(1)**2+B(2)**2+B(3)**2)
            !Orthonormal (compared to magnetic field and eachother)
            vr_tilde=(-(b_norm_r)*particle_tmp%v(3)+b_norm_phi*particle_tmp%v(1))/sqrt(b_norm_phi**2+b_norm_r**2)
            vz_tilde = (particle_tmp%v(2)-b_norm_z*(b_norm_phi*particle_tmp%v(3)+b_norm_r*particle_tmp%v(1)+b_norm_z*particle_tmp%v(2)))/sqrt(b_norm_phi**2+b_norm_r**2)
            v_par= b_norm_phi*particle_tmp%v(3)+b_norm_r*particle_tmp%v(1)+b_norm_z*particle_tmp%v(2)
            !write(*,*) v_par, dot_product(particle_tmp%v, B)/norm2(B)

            !Parallel & perpendicular pressures
            p_perp = 1.d0/2.d0*(vr_tilde**2+vz_tilde**2)
            p_par = v_par**2
            p_atrop = p_par-p_perp

            ! Calculating indices of the support in the values array of the test_phase, the array is the [p_phi, E] value of the particle.
            ! This is a bit trial and error to get good 'looking' bandwidths. (call loop w/ 1 step & E_diff ->pctls weight for dist function
            call calc_index_shaped_part_x(test_phase,particle_tmp,index_phase_tmp, val_tmp,sim,[particle_tmp%q*el_chg*psi+sim%groups(1)%mass*MASS_PROTON*particle_tmp%v(3)*particle_tmp%x(1), dot_product(particle_tmp%v,particle_tmp%v)*0.5d0*sim%groups(1)%mass*MASS_PROTON/EL_CHG/1d3])

            ! Adding to main test_phase array of all the particle contributions.
            do i_phase=1, test_phase%totsupport
              if(index_phase_tmp(i_phase) > 0)then
                phase_proj(index_phase_tmp(i_phase))=phase_proj(index_phase_tmp(i_phase))+1.d0*val_tmp(i_phase)*E_diff
              endif
            enddo
            if(particle_tmp%i_elm .gt.0) then





              call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
              call mode_moivre(particle_tmp%x(3), HZ)

              i_elm=particle_tmp%i_elm

              do l=1,n_vertex_max
                do m=1,n_order+1

                  index_lm = (l-1)*(n_order+1) + m

                  v = HH(l,m) * sim%fields%element_list%element(i_elm)%size(l,m)

                  do i_tor=1,n_tor


                    feedback_rhs(m,l,i_elm,i_tor,1) = feedback_rhs(m,l,i_elm,i_tor,1) &

                                                           + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &

                                                            * (p_perp+b_norm_r**2*p_atrop) * mu_zero !PI_RR
                    feedback_rhs(m,l,i_elm,i_tor,2) = feedback_rhs(m,l,i_elm,i_tor,2) &

                                                            + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &

                                                            * ( p_perp+b_norm_z**2*p_atrop ) * mu_zero                       !PI_ZZ
                    feedback_rhs(m,l,i_elm,i_tor,3) = feedback_rhs(m,l,i_elm,i_tor,3) &

                                                            + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &

                                                            * (p_perp+b_norm_phi**2*p_atrop) * mu_zero !PI_PHIPHI
                    feedback_rhs(m,l,i_elm,i_tor,4) = feedback_rhs(m,l,i_elm,i_tor,4) &

                                                            + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &

                                                            * ( b_norm_r*b_norm_z*p_atrop ) * mu_zero                            !PI_RZ
                    feedback_rhs(m,l,i_elm,i_tor,5) = feedback_rhs(m,l,i_elm,i_tor,5) &

                                                           + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &

                                                            * (b_norm_r*b_norm_phi*p_atrop ) * mu_zero !PI_RPHI
                    feedback_rhs(m,l,i_elm,i_tor,6) = feedback_rhs(m,l,i_elm,i_tor,6) &

                                                           + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &

                                                            *(b_norm_z*b_norm_phi*p_atrop ) * mu_zero !PI_ZPHI
                    feedback_rhs(m,l,i_elm,i_tor,7) = feedback_rhs(m,l,i_elm,i_tor,7) &

                                                           + HZ(i_tor) * v * particle_tmp%weight  !Density


                  enddo! <tor harmonic

                enddo   !< order
              enddo     !< vertex

            end if !<particle_temp%i_elm gt 0 after pushing!
          end if !<particle_temp%i_elm gt 0

        end do ! steps

        call copy_particle_kinetic_leapfrog(particle_tmp, particles(j))





      end do   ! particles
      !$omp end parallel do

      if (sim%my_id .eq. 0) write(*,*) "End of the particle loop"


  end select
  select type(particles => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)
      do j=1,size(particles,1)
        E_after = E_after + 0.5d0*particles(j)%weight*sim%groups(1)%mass*mass_proton*dot_product(particles(j)%v, particles(j)%v)
      enddo
  end select
  call MPI_REDUCE(E_after,E_after_red,1,MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  write(*,*) "On process", sim%my_id, "tot energy after: ", E_tot
  if(sim%my_id .eq. 0 ) write(*,*) "Total energy after:", E_after_red, " Energy diff = ", E_after_red - E_tot_red
  jorek_feedback%rhs = feedback_rhs/n_steps
  test_phase%values=test_phase%values+phase_proj
  deallocate(feedback_rhs)
  deallocate(phase_proj)

end subroutine


pure function f_adapted(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 ::s, coeff(0:4)
  real*4 :: f

  coeff(0)=0.53
  coeff(1)=0.3
  coeff(2)=0.2
  coeff(3)=0.52
  coeff(4)=0.26

  s = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

  f = (f - coeff(4)) / (1.d0 - coeff(4))

end function f_adapted

pure function f_original(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 ::s, coeff(0:3)
  real*4 :: f

  ! central densiy should be 1.44131x10^17

  coeff(0)=0.49123
  coeff(1)=0.298228
  coeff(2)=0.198739
  coeff(3)=0.521298

  s = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

end function f_original

pure function f_toroidal_flux(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 :: s, psi_norm, coeff(0:3)
  real*4 :: f

  ! central densiy should be 1.44131x10^17

  coeff(0)=0.49123
  coeff(1)=0.298228
  coeff(2)=0.198739
  coeff(3)=0.521298

  psi_norm = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  s = 0.957 * psi_norm + 0.043 * psi_norm**2

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

end function f_toroidal_flux

end program tae_loop

