!> mod_gyroaverage_synchrotron_light_dist_test contains all variables and
!> procedures for testing the gyroaveraged synchrotron light vertices
module mod_gyroaverage_synchrotron_light_dist_vertices_test
use fruit
use mod_particle_types,                              only: particle_gc_relativistic_id
use mod_particle_types,                              only: particle_kinetic_relativistic_id
use mod_particle_types,                              only: particle_kinetic_id
use mod_particle_types,                              only: particle_gc_vpar_id
use mod_particle_sim,                                only: particle_sim
use mod_spectra_deterministic,                       only: spectrum_integrator_2nd
use mod_gyroaverage_synchrotron_light_dist_vertices, only: gyroaverage_synchrotron_light_dist
implicit none

private
public :: run_fruit_gyroaverage_synchrotron_light_dist_vertices

!> Variables ---------------------------------------------------------
!> general parameters
real*8,parameter :: tol_real8=1d-16
real*8,parameter :: tol2_real8=1d-16
real*8,parameter :: mass_RE=5.48579909065d-4
!> parameters for generating synchrotron lights
integer,parameter :: n_mhd_sol=16
integer,parameter :: n_properties=10
integer,parameter :: n_x=3
integer,parameter :: fill_type_base=1 !< use cylindrical initialisation
integer,parameter :: n_times_sol=2
integer,parameter :: n_particle_types_check_sol=1
integer,dimension(n_times_sol),parameter              :: n_groups_per_sim=(/4,3/)
integer,parameter                                     :: n_groups_max=maxval(n_groups_per_sim)
integer,dimension(n_groups_max,n_times_sol),parameter :: n_particles_per_group=&
           reshape((/132,324,42,10,237,143,23,34/),shape(n_particles_per_group))
integer,parameter :: n_particles_max=maxval(n_particles_per_group)
integer,dimension(n_groups_max,n_times_sol),parameter :: particle_types_sol=&
           reshape((/particle_kinetic_relativistic_id,particle_gc_relativistic_id,&
           particle_kinetic_id,particle_gc_vpar_id,particle_gc_relativistic_id,&
           particle_gc_vpar_id,particle_kinetic_id,0/),shape(particle_types_sol))
real*8,parameter                                      :: survival_threshold=0.45
!> parameters for generating spectra
integer,parameter                                     :: n_spectra=2
integer,parameter                                     :: n_lines_per_spectrum=16
real*8,dimension(n_spectra),parameter                 :: min_wlen=(/3d-6,2.5d-7/)
real*8,dimension(n_spectra),parameter                 :: max_wlen=(/3.5d-6,4.2d-7/)
!> parameters for generating shadowed points
integer,parameter                                     :: n_shaded_points=57
integer,parameter                                     :: n_shaded_points_per_particle=7
real*8,dimension(2),parameter                         :: length_shadowed=(/2d-2,7d1/) 
!> variables for generating synchrotron lights
type(gyroaverage_synchrotron_light_dist)              :: vertex_sol
type(particle_sim),dimension(n_times_sol)             :: sims_particles
integer                                               :: n_gc_RE_max
integer,dimension(n_particle_types_check_sol)         :: particle_types_check_sol
integer,dimension(n_times_sol)                        :: n_active_vertices_sol
integer,dimension(n_groups_max,n_times_sol)           :: n_active_particles_sol
integer,dimension(n_particles_max,n_groups_max,n_times_sol) :: active_particle_ids_sol
real*8,dimension(n_times_sol)                         :: time_vector_sol
real*8,dimension(n_x,n_particles_max*n_groups_max,n_times_sol)          :: x_cart_sol
real*8,dimension(n_properties,n_particles_max*n_groups_max,n_times_sol) :: properties_sol
!> variables for generating spectra
type(spectrum_integrator_2nd)                         :: spectrum
!> variables for generating shadowed points
real*8,dimension(:,:,:,:),allocatable                 :: x_shadowed

!> Interfaces --------------------------------------------------------
contains
!> Fruit test basket -------------------------------------------------
!> fruit basket having all set-up, tests and tearing-down procedures
subroutine run_fruit_gyroaverage_synchrotron_light_dist_vertices()
  implicit none
  write(*,'(/A)') "  ... setting-up: gyroaverage synchrotron light vertices tests"
  call setup
  write(*,'(/A)') "  ... running: gyroaverage synchrotron light vertices tests"
  write(*,'(/A)') "  ... tearing-down: gyroaverage synchrotron light vertices tests"
  call teardown
end subroutine run_fruit_gyroaverage_synchrotron_light_dist_vertices

!> Set-up and tear-down procedures------------------------------------
!> allocate and initialise the unit test features
subroutine setup()
  use mod_rng,                        only: type_rng
  use mod_pcg32_rng,                  only: pcg32_rng
  use mod_gnu_rng,                    only: gnu_rng_interval
  use mod_common_test_tools,          only: omp_initialize_rngs
  use mod_particle_common_test_tools, only: sim_time_interval
  use mod_particle_common_test_tools, only: allocate_one_particle_list_type
  use mod_particle_common_test_tools, only: fill_groups,fill_mass_RE 
  use mod_particle_common_test_tools, only: fill_particles_tokamak
  use mod_particle_common_test_tools, only: invalidate_particles
  use mod_particle_common_test_tools, only: obtain_active_particle_ids
  use mod_fields_linear,              only: jorek_fields_interp_linear
  !$ use omp_lib
  implicit none
  !> variables
  integer :: ii,jj,ifail,n_gc_RE_max_loc,n_threads
  class(type_rng),dimension(:),allocatable :: rngs
  !> initialisation
  particle_types_check_sol = (/particle_gc_relativistic_id/)
  vertex_sol%n_property_vertex = n_properties; ifail = 0;
  n_gc_RE_max = 0; n_active_particles_sol = 0; n_threads = 1;
  !$ n_threads = omp_get_max_threads()
  call gnu_rng_interval(n_times_sol,sim_time_interval,time_vector_sol)
  call vertex_sol%allocate_vertices(n_times_sol,n_particles_max*n_groups_max)
  !> allocate and initialise particle lists
  do ii=1,n_times_sol
    sims_particles(ii)%time = time_vector_sol(ii)
    allocate(jorek_fields_interp_linear::sims_particles(ii)%fields)
    allocate(sims_particles(ii)%groups(n_groups_per_sim(ii)))
    call allocate_one_particle_list_type(n_groups_per_sim(ii),&
    n_particles_per_group(1:n_groups_per_sim(ii),ii),&
    particle_types_sol(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups,ifail)
    call fill_groups(n_groups_per_sim(ii),sims_particles(ii)%groups)
    call fill_mass_RE(n_groups_per_sim(ii),sims_particles(ii)%groups)
    call fill_particles_tokamak(n_groups_per_sim(ii),sims_particles(ii)%groups,fill_type_base)
    call invalidate_particles(n_groups_per_sim(ii),n_particles_max,survival_threshold,&
    n_active_particles_sol(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    call obtain_active_particle_ids(n_groups_per_sim(ii),n_particles_max,&
    active_particle_ids_sol(:,1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    n_active_vertices_sol(ii) = sum(n_active_particles_sol(:,ii))
    where(particle_types_sol(:,ii).ne.particle_gc_relativistic_id)
      n_active_particles_sol(:,ii) = 0
    endwhere
    n_gc_RE_max_loc = 0
    do jj=1,n_groups_per_sim(ii)
      if(particle_types_sol(jj,ii).ne.particle_gc_relativistic_id) cycle
      n_gc_RE_max_loc = n_gc_RE_max_loc + n_particles_per_group(jj,ii)
    enddo
      n_gc_RE_max = max(n_gc_RE_max,n_gc_RE_max_loc)
  enddo
  !> initialise positions and properties tables
  call compute_gyroavg_synch_x_properties_ana()
  !> initialise deterministic spectra
  spectrum = spectrum_integrator_2nd(n_lines_per_spectrum,n_spectra,min_wlen,max_wlen)
  call spectrum%generate_spectrum()
  !> generate shadowed point positions
  allocate(x_shadowed(n_x,n_shaded_points_per_particle,n_gc_RE_max,n_times_sol))
  call compute_x_shadowed_gc
end subroutine setup 

!> destroy all test features
subroutine teardown()
  implicit none
  call vertex_sol%deallocate_vertices; call spectrum%deallocate_spectrum;
  if(allocated(x_shadowed)) deallocate(x_shadowed)
end subroutine teardown

!> Tests -------------------------------------------------------------

!> Tools -------------------------------------------------------------
!> generate shadowed points for each light. The shadowed point position
!> is taken within the emission cone of the synchrotron radiation.
!> The emission cone half angle is approximated with 
!> cos(theta) = 1/(2*rel_fact). The relativistic factor is defined as:
!> rel_fact = sqrt(1+(p/(mass*c))**2) with p the total gc momentum and
!> c the speed of light
subroutine compute_x_shadowed_gc()
  use mod_sampling, only: sample_uniform_cone
  implicit none
  !> variables
  integer :: ii,jj,kk,n_gc_time
  real*8  :: cos_half_angle
  real*8,dimension(n_x) :: v_gc_dir,x_gc,rng
  !> generate shadowed points
  do kk=1,n_times_sol
    n_gc_time = sum(n_active_particles_sol(:,jj))
    do jj=1,n_gc_time
      x_gc = x_cart_sol(:,jj,kk)
      v_gc_dir = properties_sol(1:3,jj,kk)
      cos_half_angle = cos(5d-1/properties_sol(4,jj,kk))
      do ii=1,n_shaded_points_per_particle
        call random_number(rng)
        x_shadowed(:,ii,jj,kk) = sample_uniform_cone(cos_half_angle,&
        rng,v_gc_dir,x_gc,length_shadowed)
      enddo
    enddo
  enddo
end subroutine compute_x_shadowed_gc

!> compute & fill the gyroaverage synchrotron light positions and properties
subroutine compute_gyroavg_synch_x_properties_ana()
  use mod_coordinate_transforms, only: cylindrical_to_cartesian
  use mod_particle_types,        only: particle_gc_relativistic
  implicit none
  !> variables
  integer :: ii,jj,kk,counter
  !> initialise positions and properties arrays
  x_cart_sol = 0d0; properties_sol = 0d0;
  !> fill property table
  do kk=1,n_times_sol
    counter = 0
    do jj=1,n_groups_per_sim(kk)
      select type(p_list=>sims_particles(kk)%groups(jj)%particles)
        type is (particle_gc_relativistic)
        do ii=1,n_particles_per_group(jj,kk)
          if(p_list(ii)%i_elm.eq.0) cycle
          counter = counter + 1
          x_cart_sol(:,counter,kk) = cylindrical_to_cartesian(p_list(ii)%x)
          call compute_gyroavg_synch_properties_ana_1p(p_list(ii),&
          sims_particles(kk)%groups(jj)%mass,properties_sol(:,counter,kk))
        enddo
      end select
    enddo
  enddo
end subroutine compute_gyroavg_synch_x_properties_ana

!> compute gyroaverage synchrotron light properties using the analytical
!> tokamak-like MHD fields for one guiding center
subroutine compute_gyroavg_synch_properties_ana_1p(gc_in,mass,properties)
  use constants,                      only: PI,EL_CHG,EPS_ZERO
  use constants,                      only: ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
  use mod_coordinate_transforms,      only: vector_cylindrical_to_cartesian
  use mod_coordinate_transforms,      only: cylindrical_to_cartesian_velocity
  use mod_particle_types,             only: particle_gc_relativistic
  use mod_gc_relativistic,            only: compute_relativistic_gc_rhs
  use mod_particle_common_test_tools, only: compute_test_E_B_normB_gradB_curlb_Dbdt_fields
  implicit none
  !> inputs:
  type(particle_gc_relativistic),intent(in)  :: gc_in
  real*8,intent(in)                          :: mass 
  !> outputs:
  real*8,dimension(n_properties),intent(out) :: properties
  !> variables:
  real*8                :: normB,rel_fact,rel_fact_parallel
  real*8                :: thetap,charge,p_perp
  real*8,dimension(n_x) :: E_field,b_field,gradB,curlb,dbdt
  real*8,dimension(4)   :: x_gc_velocity
  
  !> compute the analytical MHD fields for computing the GC velocity
  call compute_test_E_B_normB_gradB_curlb_Dbdt_fields(gc_in%x,E_field,&
       b_field,normB,gradB,curlb,dbdt)
  E_field = vector_cylindrical_to_cartesian(gc_in%x(3),E_field)
  b_field = vector_cylindrical_to_cartesian(gc_in%x(3),b_field)
  gradB   = vector_cylindrical_to_cartesian(gc_in%x(3),gradB)
  curlb   = vector_cylindrical_to_cartesian(gc_in%x(3),curlb)
  dbdt    = vector_cylindrical_to_cartesian(gc_in%x(3),dbdt)
  !> compute the guiding center velocity
  x_gc_velocity = compute_relativistic_gc_rhs(int(gc_in%q,kind=4),mass,gc_in%x(2),&
  gc_in%x(1),gc_in%p(1),normB,E_field,b_field,gradB,curlb,dbdt)
  x_gc_velocity(1:3) = cylindrical_to_cartesian_velocity(&
  gc_in%x(1),gc_in%x(3),x_gc_velocity(1:3))
  properties(1:3) = x_gc_velocity(1:3)/norm2(x_gc_velocity(1:3))
  !> compute the relativistic factor
  rel_fact = 1d0 + ((gc_in%p(1)/(mass*SPEED_OF_LIGHT))**2) 
  rel_fact_parallel = sqrt(rel_fact)
  rel_fact = sqrt(rel_fact + ((2d0*gc_in%p(2)*normB)/(mass*(SPEED_OF_LIGHT**2))))
  properties(4) = rel_fact
  !> compute the beta
  properties(5) = sqrt(1d0 - (1d0/(rel_fact**2)));
  !> compute the pitch angle
  p_perp = sqrt(2d0*mass*gc_in%p(2)*normB); thetap = atan2(p_perp,gc_in%p(1));
  properties(6:7) = (/cos(thetap),sin(thetap)/)
  !> critical wavelength
  charge = real(gc_in%q,kind=8)*EL_CHG
  properties(8) = (4d0*PI*mass*ATOMIC_MASS_UNIT*SPEED_OF_LIGHT*rel_fact_parallel)/&
                  (3d0*charge*normB*(rel_fact**2))
  !> compute the directionality function intensity
  properties(9) = (27d0*charge*normB*(rel_fact**7))/&
                  (128d0*(PI**2)*mass*ATOMIC_MASS_UNIT*SPEED_OF_LIGHT*&
                  (sin(thetap)**2)*(rel_fact_parallel**4))
  !> compute the synchrotron power normalisation
  p_perp = p_perp/(mass*rel_fact) !< perpendicular velocity
  properties(10) = ((charge**4)*((normB*rel_fact*rel_fact_parallel*p_perp)**2))/&
                   (6d0*PI*EPS_ZERO*((mass*ATOMIC_MASS_UNIT)**2)*(SPEED_OF_LIGHT**3))
end subroutine compute_gyroavg_synch_properties_ana_1p

!>--------------------------------------------------------------------
end module mod_gyroaverage_synchrotron_light_dist_vertices_test
 
