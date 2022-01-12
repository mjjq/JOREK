!> mod_synchrotron_light_vertices_test contains all variables and
!> procedures for testing the synchrotron light vertices
module mod_synchrotron_light_vertices_test
use fruit
use mod_particle_sim,               only: particle_sim
use mod_vertices,                   only: n_x
use mod_synchrotron_light_vertices, only: n_properties
use mod_synchrotron_light_vertices, only: synchrotron_light_vertices
implicit none

private
public :: run_fruit_synchrotron_light_vertices

!> Variables ---------------------------------------------------------
integer,parameter :: fill_type_base=1 !< use cylindrical initialisation
integer,parameter :: n_times_sol=2
integer,dimension(n_times_sol),parameter :: n_groups_per_sim=(/3,2/)
integer,parameter                        :: n_groups_max=maxval(n_groups_per_sim)
integer,dimension(n_groups_max,n_times_sol),parameter :: n_particles_per_group=&
           reshape((/135,247,512,367,413,0/),shape(n_particles_per_group))
integer,parameter                        :: n_particles_max=maxval(n_particles_per_group)
real*8,parameter                         :: survival_threshold=0.33
real*8,parameter                         :: tol_real8=2.5d-11
real*8,parameter                         :: mass_RE=5.48579909065d-4
type(synchrotron_light_vertices)            :: vertex_sol
type(particle_sim),dimension(n_times_sol)   :: sims_particles
integer,dimension(n_times_sol)              :: n_active_vertices_sol
integer,dimension(n_groups_max,n_times_sol) :: n_active_particles_sol
integer,dimension(n_particles_max,n_groups_max,n_times_sol) :: active_particle_ids_sol
real*8,dimension(n_times_sol)               :: time_vector_sol
real*8,dimension(n_x,n_particles_max*n_groups_max,n_times_sol) :: x_cart_sol
real*8,dimension(n_properties,n_particles_max*n_groups_max,n_times_sol) :: properties_sol

!> Interfaces --------------------------------------------------------
contains
!> Fruit test basket -------------------------------------------------
!> fruit basket having all set-up, tests and tearing-down procedures
subroutine run_fruit_synchrotron_light_vertices()
  implicit none
  write(*,'(/A)') "  ... setting-up: synchrotron light vertices tests"
  call setup
  write(*,'(/A)') "  ... running: synchrotron light vertices tests"
  call test_compute_synchrotron_light_properties
  call test_fill_synchrotron_lights_from_particles
  write(*,'(/A)') "  ... tearing-down: synchrotron light vertices tests"
end subroutine run_fruit_synchrotron_light_vertices

!> Set-up and tear-down procedures------------------------------------
!> allocate and initialise the unit test features
subroutine setup()
  use mod_particle_types,                   only: particle_gc_vpar_id
  use mod_particle_types,                   only: particle_kinetic_id
  use mod_particle_types,                   only: particle_kinetic_relativistic_id
  use mod_gnu_rng,                          only: gnu_rng_interval
  use mod_particle_common_test_tools,       only: sim_time_interval
  use mod_particle_common_test_tools,       only: allocate_one_particle_list_type 
  use mod_particle_common_test_tools,       only: fill_groups,fill_mass_RE
  use mod_particle_common_test_tools,       only: fill_particles_tokamak
  use mod_particle_common_test_tools,       only: invalidate_particles
  use mod_particle_common_test_tools,       only: obtain_active_particle_ids
  implicit none
  !> variables
  integer,dimension(n_groups_max,n_times_sol),parameter :: particle_types=&
  reshape((/particle_kinetic_relativistic_id,particle_gc_vpar_id,&
  particle_kinetic_relativistic_id,particle_kinetic_id,&
  particle_kinetic_relativistic_id,-1/),shape(particle_types))
  integer :: ii,ifail
  !> initialisation
  vertex_sol%n_property_vertex = n_properties; ifail = 0; n_active_particles_sol = 0;
  call gnu_rng_interval(n_times_sol,sim_time_interval,time_vector_sol)
  call vertex_sol%allocate_vertices(n_times_sol,n_particles_max*n_groups_max)

  !> allocate and initialise particle list
  do ii=1,n_times_sol
    sims_particles(ii)%time = time_vector_sol(ii)
    allocate(sims_particles(ii)%groups(n_groups_per_sim(ii)))
    call allocate_one_particle_list_type(n_groups_per_sim(ii),&
    n_particles_per_group(1:n_groups_per_sim(ii),ii),&
    particle_types(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups,ifail)
    call fill_groups(n_groups_per_sim(ii),sims_particles(ii)%groups)
    call fill_mass_RE(n_groups_per_sim(ii),sims_particles(ii)%groups)
    call fill_particles_tokamak(n_groups_per_sim(ii),sims_particles(ii)%groups,fill_type_base)
    call invalidate_particles(n_groups_per_sim(ii),n_particles_max,survival_threshold,&
    n_active_particles_sol(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    call obtain_active_particle_ids(n_groups_per_sim(ii),n_particles_max,&
    active_particle_ids_sol(:,1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    n_active_vertices_sol(ii) = sum(n_active_particles_sol(:,ii))
  enddo
  !> initialise positions and properties tables
  call compute_synch_x_properties_ana()
end subroutine setup

!> Tests -------------------------------------------------------------
!> test fill synchrotron lights from particles
subroutine test_fill_synchrotron_lights_from_particles()
  use mod_assert_equals_tools,        only: assert_equals_rel_error
  use mod_synchrotron_light_vertices, only: fill_synchrotron_lights_from_particles
  implicit none
  !> variables
  !> fill synchrotron lights from particles
  call fill_synchrotron_lights_from_particles(vertex_sol,sims_particles,&
  n_groups_max,n_particles_max,n_groups_per_sim,n_active_particles_sol,&
  active_particle_ids_sol)
  call assert_equals_rel_error(n_x,n_particles_max*n_groups_max,n_times_sol,&
  vertex_sol%x,x_cart_sol,tol_real8,&
  "Error fill synchrotron lights from particles: properties errors too large!")
  call assert_equals_rel_error(n_properties,n_particles_max*n_groups_max,&
  n_times_sol,vertex_sol%properties,properties_sol,tol_real8,&
  "Error fill synchrotron lights from particles: properties errors too large!")
end subroutine test_fill_synchrotron_lights_from_particles

!> test the property function of synchrotron light properties
subroutine test_compute_synchrotron_light_properties()
  use mod_coordinate_transforms,      only: vector_cylindrical_to_cartesian
  use mod_particle_types,             only: particle_kinetic_relativistic
  use mod_synchrotron_light_vertices, only: compute_synchrotron_light_properties
  use mod_particle_common_test_tools, only: compute_test_E_B_fields 
  implicit none
  !> variables
  integer :: ii,jj,kk,counter
  real*8,dimension(3) :: E_field,B_field
  real*8,dimension(n_properties) :: properties
  real*8,dimension(n_properties,n_particles_max*n_groups_max) :: error,zeros
  !> loop for computing the properties and testing
  zeros = 0.d0
  do kk=1,n_times_sol
    counter = 0; error = 0.d0;
    do jj=1,n_groups_per_sim(kk)
      select type(p_list=>sims_particles(kk)%groups(jj)%particles)
        type is(particle_kinetic_relativistic)
        do ii=1,n_particles_per_group(jj,kk)
          if(p_list(ii)%i_elm.le.0) cycle
          counter = counter + 1
          call compute_test_E_B_fields(p_list(ii)%x,E_field,B_field)
          B_field = vector_cylindrical_to_cartesian(p_list(ii)%x(3),B_field)
          E_field = vector_cylindrical_to_cartesian(p_list(ii)%x(3),E_field)
          call compute_synchrotron_light_properties(3,n_properties,&
          p_list(ii),sims_particles(kk)%groups(jj)%mass,E_field,B_field,properties)
          error(:,counter) = abs((properties - properties_sol(:,counter,kk))/&
          properties_sol(:,counter,kk))
        enddo
      end select
    enddo
    !> check if the properties arrays are equal
    call assert_equals(error,zeros,n_properties,n_particles_max*n_groups_max,&
    tol_real8,"Error synchrotron light compute properties: too large errors!")
  enddo
end subroutine test_compute_synchrotron_light_properties

!> Tools -------------------------------------------------------------
!> compute and fill particle positions and properties for RE
subroutine compute_synch_x_properties_ana()
  use mod_coordinate_transforms, only: cylindrical_to_cartesian
  use mod_particle_types,        only: particle_kinetic_relativistic
  implicit none
  !> variables
  integer :: ii,jj,kk,counter
  !> initialise positions and properties arrays
  x_cart_sol = 0.d0; properties_sol = 0.d0;
  !> fill property table
  do kk=1,n_times_sol
    counter = 0
    do jj=1,n_groups_per_sim(kk)
      select type(p_list=>sims_particles(kk)%groups(jj)%particles)
      type is(particle_kinetic_relativistic)
        do ii=1,n_particles_per_group(jj,kk)
          if(p_list(ii)%i_elm.le.0) cycle
          counter = counter + 1
          x_cart_sol(:,counter,kk) = cylindrical_to_cartesian(p_list(ii)%x)
          call compute_synch_properties_ana_1p(p_list(ii)%x(3),&
          sims_particles(kk)%groups(jj)%mass,p_list(ii),&
          properties_sol(:,counter,kk))
        enddo
      end select
    enddo
  enddo
end subroutine compute_synch_x_properties_ana

!> compute synchrotron electron properties using the analytical
!> tokamak like electric and magnetic fields for one particle
subroutine compute_synch_properties_ana_1p(phi,mass,particle,property)
  use constants,                      only: PI,EL_CHG,ATOMIC_MASS_UNIT
  use constants,                      only: SPEED_OF_LIGHT,EPS_ZERO
  use mod_math_operators,             only: cross_product
  use mod_coordinate_transforms,      only: vector_cylindrical_to_cartesian
  use mod_particle_types,             only: particle_kinetic_relativistic
  use mod_particle_common_test_tools, only: compute_test_E_B_fields 
  implicit none
  !> inputs-outputs
  type(particle_kinetic_relativistic),intent(inout) :: particle
  !> inputs
  real*8,intent(in) :: mass,phi
  !> outputs
  real*8,dimension(n_properties),intent(out) :: property
  !> variables
  real*8 :: velocity,beta,rel_fact,kappa,P_rad
  real*8,dimension(3) :: vel_vec,E_field,B_field,T_vec,N_vec,B_vec
  real*8,dimension(3) :: vec_real_size3
  !> compute relativistic factor, velocity and beta
  rel_fact = sqrt(1.d0+(dot_product(particle%p,particle%p)/&
             (SPEED_OF_LIGHT*SPEED_OF_LIGHT*mass*mass)))
  vel_vec = particle%p/(mass*rel_fact)
  velocity = norm2(vel_vec)
  beta     =  velocity/SPEED_OF_LIGHT
  T_vec   = vel_vec/velocity
  !> compute electric and magnetic field
  call compute_test_E_B_fields(particle%x,E_field,B_field)
  B_field = vector_cylindrical_to_cartesian(phi,B_field)
  E_field = vector_cylindrical_to_cartesian(phi,E_field)
  !> compute normal and binormal vectors
  N_vec = E_field + cross_product(vel_vec,B_field) - dot_product(T_vec,E_field)*T_vec
  N_vec = N_vec/norm2(N_vec)
  B_vec = cross_product(T_vec,N_vec); B_vec = B_vec/norm2(B_vec);
  !> compute orbit curvature (L. Carbajal, PPCF, 2017)
  vec_real_size3 = E_field + cross_product(vel_vec,B_field)
  vec_real_size3 = cross_product(vel_vec,vec_real_size3)
  kappa = (norm2(vec_real_size3)*EL_CHG*abs(real(particle%q,kind=8)))/&
  (rel_fact*mass*ATOMIC_MASS_UNIT*velocity**3.d0)
  !> compute total radiated power (L. Carbajal, PPCF, 2017)
  P_rad = (((rel_fact*velocity)**4.d0)*((kappa*EL_CHG*real(particle%q,kind=8))**2.d0))/&
  (6.d0*PI*EPS_ZERO*(SPEED_OF_LIGHT**3.d0))
  !> Store all values in the array
  property(1:3) = T_vec; property(4:6) = N_vec; property(7:9) = B_vec;
  property(10) = beta; property(11) = rel_fact; property(12) = kappa;
  property(13) = P_rad;
end subroutine compute_synch_properties_ana_1p

!>--------------------------------------------------------------------
end module mod_synchrotron_light_vertices_test
