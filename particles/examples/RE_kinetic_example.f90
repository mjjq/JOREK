!> particles/examples/RE_full_kinetic_example.f90
!> Requires: JOREK input file
!> Compile: make RE_example, run with: RE_example < jorek_inputfile
!> Specify the jorek restart to be read in the variable jorek_filename
program RE_kinetic_example
use constants, only: PI,TWOPI,SPEED_OF_LIGHT
use mod_model_settings, only: var_rho,n_var
use mod_random_seed
use particle_tracer
use mod_kinetic_relativistic
use mod_particle_diagnostics
use mod_pusher_tools, only: get_orthonormals
!$ use omp_lib
implicit none

!> Variable declarations -----------------------------------------------
type(sobseq_rng),dimension(:),allocatable  :: sob_rngs !< sobolev seq. rng
type(event)                                :: field_reader
type(write_particle_diagnostics)           :: diag
logical         :: diag_append
integer(kind=1) :: q_e
integer         :: ii,jj,kk
integer         :: n_steps,n_write_steps,n_max_threads
integer         :: n_groups,n_particles,n_mhd_fields
integer         :: seed,thread_id,ifail,writename_decimal_digits
integer         :: writename_fractional_digits
integer,dimension(:),allocatable :: mhd_field_ids
real*8                           :: t_step,stop_time,time,mass_e,t_target
real*8                           :: rest_energy,B_norm,psi,U,start_time
real*8,dimension(2)              :: energy_kin,pitch_angles,gyro_angles
real*8,dimension(2)              :: p_int,cos_pitch_int
real*8,dimension(3)              :: rands,b,E,e1,e2
real*8,dimension(:),allocatable  :: t_steps
character(len=:),allocatable     :: jorek_filename,diag_filename

!> MPI initialization --------------------------------------------------
n_groups = 1              !< number of particle groups
call sim%initialize(num_groups=n_groups) !< open the MPI communicator

!> Define the inputs ----------------------------------------------------
ifail = 0
diag_append = .true.             !< append diagnostic values
writename_decimal_digits    = 5  !< number of decimal digits of the restart file filename
writename_fractional_digits = 13 !< number of fractional digits of the restart file filename
n_particles = 1000               !< number of particles per group
n_mhd_fields = 1                 !< number of required mhd fields for particle init
n_write_steps = 100              !< number of time steps between writien actions
t_step = 1.d-13                  !< time step
start_time = 0.d0                !< initial simulation time
stop_time = 1.d-9;               !< time at which the simulation is stop
mass_e = 5.48579909065d-4        !< electron mass in AMU
rest_energy = 0.51099895         !< electron rest energy in MeV/c^2
energy_kin = (/2.d1,2.d1/)       !< kinetic energy in MeV
pitch_angles  = (/PI-0.289,PI/)  !< pitch angle in radians
gyro_angles = (/0.d0,TWOPI/)     !< gyro angles in radians
q_e = -1                         !< electron charge
!> allocate time steps
allocate(t_steps(n_groups)); t_steps = t_step;
!> allocate and set the mhd field ids
allocate(mhd_field_ids(n_mhd_fields)); mhd_field_ids = (/var_rho/);
allocate(character(len=25)::jorek_filename); 
jorek_filename = 'jorek_equilibrium';
allocate(character(len=10)::diag_filename)
diag_filename = 'RE_diag.h5'

!> Initialise the simulation --------------------------------------------
!> compute initialisation values
p_int = mass_e*SPEED_OF_LIGHT*sqrt((((energy_kin/rest_energy)+1)**2)-1.d0)
cos_pitch_int = cos(pitch_angles)
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
basename=trim(jorek_filename),i=-1)) !< read the jorek fields
call with(sim,field_reader)
!> write diagnostics
diag = write_particle_diagnostics(filename=trim(diag_filename),append=diag_append)
!> set the events
sim%time = start_time
events = [event(write_action(decimal_digits=writename_decimal_digits,&
         fractional_digits=writename_fractional_digits),&
         start=sim%time,step=t_step*real(n_write_steps,kind=8)),&
         event(diag,start=sim%time,step=t_step*real(n_write_steps,kind=8)),&
         event(stop_action(),start=stop_time)]

!> Initialise particle population --------------------------------------
write(*,*) t_step*real(n_write_steps,kind=8)
call check_and_fix_timesteps(t_steps,events) !< check the time step
write(*,*) "----------------------------------"
write(*,*) "Running particle initialisation"
write(*,*) "----------------------------------"
!> Initialise particle groups
do ii=1,n_groups
  !> initialising runaway positions and mass
  sim%groups(ii)%mass = mass_e
  allocate(particle_kinetic_relativistic::sim%groups(ii)%particles(n_particles))
  call initialise_particles(sim%groups(ii)%particles,sim%fields%node_list,&
  sim%fields%element_list,sobseq_rng(),variables=mhd_field_ids,&
  transform=transform_accept_mhd_fields_value,normalise_uniform_space_rej_vars_in=.true.)

  !> initialise runaway electron energy
  !$omp parallel default(private) shared(sim,sob_rngs) &
  !$omp firstprivate(ii,thread_id,q_e,p_int,cos_pitch_int,gyro_angles)
  !$ thread_id = omp_get_thread_num()
  select type (particles=>sim%groups(ii)%particles)
    type is (particle_kinetic_relativistic)
    !$omp do
    do jj=1,size(particles)
      if(particles(ii)%i_elm.le.0) cycle
      call sob_rngs(thread_id)%next(rands)
      call sim%fields%calc_EBpsiU(0.d0,particles(jj)%i_elm,particles(jj)%st,&
      particles(jj)%x(3),E,b,psi,U)
      B_norm = norm2(b); b = b/B_norm;
      call get_orthonormals(b,e1,e2)
      particles(jj)%p = uniform_init_kinetic_relativistic(p_int,cos_pitch_int,&
      gyro_angles,rands,b,e1,e2)
      particles(jj)%q = q_e
    enddo
    !$omp end do
  end select
  !$omp end parallel 
enddo
call with(sim,events,at=sim%time) !< store the initialised particles
write(*,*) "----------------------------------"
write(*,*) "Terminated particle initialisation"
write(*,*) "----------------------------------"

!> Integrate particle trajectory ---------------------------------------
write(*,*) "----------------------------------"
write(*,*) "Running particle simulation"
write(*,*) "----------------------------------"
!> loop until the simulation should stop
do while(.not.sim%stop_now)
  t_target = next_event_at(sim,events)
  write(*,*) "t_target: ",t_target
  !> loop on the groups
  do ii=1,size(sim%groups)
    n_steps = nint((t_target-sim%time)/t_steps(ii));
    !$omp parallel default(private) firstprivate(n_steps,t_steps,ii,ifail) &
    !$omp shared(sim)
    select type(particles=>sim%groups(ii)%particles)
      type is (particle_kinetic_relativistic)
      !$omp do
      do jj=1,size(sim%groups(ii)%particles)
        do kk=1,n_steps
          if(particles(jj)%i_elm.le.0) exit
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
write(*,*) "----------------------------------"
write(*,*) "Particle simulation terminated"
write(*,*) "----------------------------------"

!> Tear down the simulation ---------------------------------------------
deallocate(mhd_field_ids);  deallocate(jorek_filename);
deallocate(diag_filename); deallocate(sob_rngs);
call sim%finalize()
contains
!> ----------------------------------------------------------------------
!> initialise the energy of a relativistic particle: uniform initialisation
!> inputs:
!> outputs:
function uniform_init_kinetic_relativistic(p_int,cos_theta_int,phi_int,&
u,b,e1,e2) result(pxpypz)
  use mod_sampling, only: sample_uniform_sphere_corona_rthetaphi
  implicit none
  real*8,dimension(2),intent(in) :: p_int,cos_theta_int,phi_int
  real*8,dimension(3),intent(in) :: u,b,e1,e2
  real*8,dimension(3)            :: pxpypz
  real*8,dimension(3)            :: pthetaphi

  !> compute unifrom sample spherical coordinates  
  pthetaphi = sample_uniform_sphere_corona_rthetaphi(p_int,cos_theta_int,phi_int,u)
  !> compute momentum in cartesian coordinates
  pxpypz = pthetaphi(1)*(b*cos(pthetaphi(2)) + &
           sin(pthetaphi(2))*(e1*cos(pthetaphi(3)) + &
           e2*sin(pthetaphi(3))))
end function uniform_init_kinetic_relativistic

!> acceptance rejection function for particle initialisation
!> be aware: a quick a dirty solution has been implemented 
!> inputs:
!>   P: (real8)(3) normalised MHD field vector
!> outputs:
!>   f_accept_mhd_fields_value: (real4) value of the distribution
pure function transform_accept_mhd_fields_value(n_mhd,P) result(f_accept)
  implicit none
  integer,intent(in) :: n_mhd
  real*8,dimension(n_mhd),intent(in)   :: P
  real*8 :: f_accept
  f_accept = max(minval(P),0.d0);
end function transform_accept_mhd_fields_value

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
