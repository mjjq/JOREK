!> particles/examples/RE_full_kinetic_example.f90
!> Requires: JOREK input file
!> Compile: make RE_example, run with: RE_example < jorek_inputfile
!> Specify the jorek restart to be read in the variable jorek_filename
program RE_kinetic_example
use mod_model_settings, only: n_var
use mod_random_seed
use particle_tracer
use mod_kinetic_relativistic
use mod_particle_diagnostics
!$ use omp_lib
implicit none

!> Variable declarations -----------------------------------------------
type(sobseq_rng),dimension(:),allocatable  :: sob_rngs !< sobolev seq. rng
type(event)                                :: field_reader
type(write_particle_diagnostics)           :: diag
integer(kind=1) :: q_e
integer         :: ii,jj,kk
integer         :: n_steps,n_write_steps,n_max_threads
integer         :: n_groups,n_particles,n_mhd_fields
integer         :: seed,thread_id,ifail
integer,dimension(:),allocatable :: mhd_field_ids
real*8                           :: t_step,stop_time,time,mass_e,t_target
real*8,dimension(4)              :: rands
real*8,dimension(:),allocatable  :: t_steps
character(len=:),allocatable     :: jorek_filename,diag_filename

!> MPI initialization --------------------------------------------------
n_groups = 1              !< number of particle groups
call sim%initialize(num_groups=n_groups) !< open the MPI communicator

!> Define the inputs ----------------------------------------------------
ifail = 0
n_particles = 1000000     !< number of particles per group
n_mhd_fields = 1          !< number of required mhd fields for particle init
n_write_steps = 100       !< number of time steps between writien actions
t_step = 1.d-13           !< time step
stop_time = 1.d-9;        !< time at which the simulation is stop
mass_e = 5.48579909065d-4 !< electron mass in AMU
q_e = -1                  !< electron charge
!> allocate time steps
allocate(t_steps(n_groups)); t_steps = t_step;
!> allocate and set the mhd field ids
allocate(mhd_field_ids(n_mhd_fields)); mhd_field_ids = (/1/);
allocate(character(len=25)::jorek_filename); 
jorek_filename = 'jorek_equilibrium';
allocate(character(len=10)::diag_filename)
diag_filename = 'RE_diag.h5'

!> Initialise the simulation --------------------------------------------
!> initialise sobolev sequence rng
n_max_threads = 1; thread_id = 0;
!$ n_max_threads = omp_get_max_threads()
allocate(sob_rngs(0:n_max_threads-1))
do ii=0,n_max_threads-1
  call sob_rngs(ii)%initialize(size(rands),random_seed(),&
  sim%n_cpu*n_max_threads,sim%my_id*n_max_threads+ii+1,ifail)
  if(ifail.ne.0) call MPI_ABORT(MPI_COMM_WORLD,-1,ifail)
enddo
!> read jorek restart field
field_reader = event(read_jorek_fields_interp_linear(&
basename=trim(jorek_filename),i=-1),start=0.d0) !< read the jorek fields
call with(sim,field_reader)
!> write diagnostics
diag = write_particle_diagnostics(filename=trim(diag_filename))
!> set the events
events = [field_reader,event(write_action(),step=t_step*real(n_write_steps,kind=8)),&
         event(diag,step=t_step*real(n_write_steps,kind=8)),&
         event(stop_action(),start=stop_time)]

!> Initialise particle population --------------------------------------
call check_and_fix_timesteps(t_steps,events) !< check the time step
call with(sim,events,at=0.d0)
!> Initialise particle groups
do ii=1,n_groups
  !> initialising runaway positions and mass
  sim%groups(ii)%mass = mass_e
  allocate(particle_kinetic_relativistic::sim%groups(ii)%particles(n_particles))
  call initialise_particles_H_mu_psi(sim%groups(ii)%particles,&
  sim%fields,sobseq_rng(),sim%groups(ii)%mass,uniform_space=.true.,&
  uniform_space_rej_f=f_accept_mhd_fields_value,&
  uniform_space_rej_vars=mhd_field_ids,charge=int(q_e),&
  normalise_uniform_space_rej_vars_in=.true.)
  !> initialise runaway electron energy
  !$omp parallel default(private) shared(sim,sob_rngs) &
  !$omp firstprivate(ii,thread_id,q_e)
  !$ thread_id = omp_get_thread_num()
  select type (particles=>sim%groups(ii)%particles)
    type is (particle_kinetic_relativistic)
    !$omp do
    do jj=1,size(particles)
      if(particles(ii)%i_elm.le.0) cycle
      call sob_rngs(thread_id)%next(rands)
    enddo
    !$omp end do
  end select
  !$omp end parallel 
enddo

!> first dummy initialisation for checking if the integration loop makes sens
do ii=1,size(sim%groups)

enddo

!> Integrate particle trajectory ---------------------------------------
!> loop until the simulation should stop
do while(.not.sim%stop_now)
  t_target = next_event_at(sim,events)
  !> loop on the groups
  do ii=1,size(sim%groups)
    n_steps = nint((time-sim%time)/t_steps(ii));
    !$omp parallel default(private) firstprivate(n_steps,t_steps,ii,ifail) &
    !$omp shared(sim)
    select type(particles=>sim%groups(ii)%particles)
      type is (particle_kinetic_relativistic)
      !$omp do
      do jj=1,size(sim%groups(ii)%particles)
        do kk=1,n_steps
          if(particles(ii)%i_elm.le.0) exit
          time = sim%time + kk*t_steps(ii)
          call volume_preserving_push_jorek(particles(jj),&
          sim%fields,sim%groups(ii)%mass,time,t_steps(ii),ifail)          
        enddo
      enddo
      !$omp end do
    end select
    !$omp end parallel
  enddo
  sim%time = t_target
  call with(sim,events,at=sim%time)
enddo

!> Tear down the simulation ---------------------------------------------
deallocate(mhd_field_ids);  deallocate(jorek_filename);
deallocate(diag_filename); deallocate(sob_rngs);
call sim%finalize()
contains
!> ----------------------------------------------------------------------
!> acceptance rejection function for particle initialisation
!> be aware: a quick a dirty solution has been implemented 
!> inputs:
!>   n: (integer) size of the MHD field vector
!>   P: (real8)(3) normalised MHD field vector
!>   gradP: (real8)(3,n) not used
!> outputs:
!>   f_accept_mhd_fields_value: (real4) value of the distribution
pure function f_accept_mhd_fields_value(n,P,gradP) result(f_accept)
  implicit none
  integer,intent(in)               :: n
  real*8,dimension(n),intent(in)   :: P
  real*8,dimension(3,n),intent(in) :: gradP
  real*4 :: f_accept
  f_accept = max(minval(P),0.d0);
end function f_accept_mhd_fields_value
!> ----------------------------------------------------------------------
end program RE_kinetic_example
