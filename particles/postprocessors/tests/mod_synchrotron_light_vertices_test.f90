!> mod_synchrotron_light_vertices_test contains all variables and
!> procedures for testing the synchrotron light vertices
module mod_synchrotron_light_vertices_test
use fruit
use mod_particle_sim,               only: particle_sim
use mod_spectra_monte_carlo,        only: spectrum_rng_uniform
use mod_synchrotron_light_vertices, only: synchrotron_light_vertices
implicit none

private
public :: run_fruit_synchrotron_light_vertices

!> Variables ---------------------------------------------------------
!> general parameters
real*8,parameter                            :: tol_real8=2.5d-11
real*8,parameter                            :: tol2_real8=3.5d-8
real*8,parameter                            :: mass_RE=5.48579909065d-4
real*8,parameter                            :: exclusion_values=1.d-100
!> parameters for generating synhrotron lights
integer,parameter :: n_x=3
integer,parameter :: n_properties=13
integer,parameter :: fill_type_base=1 !< use cylindrical initialisation
integer,parameter :: n_times_sol=2
integer,dimension(n_times_sol),parameter    :: n_groups_per_sim=(/3,2/)
integer,parameter                           :: n_groups_max=maxval(n_groups_per_sim)
integer,dimension(n_groups_max,n_times_sol),parameter :: n_particles_per_group=&
           reshape((/135,247,512,367,413,0/),shape(n_particles_per_group))
integer,parameter                           :: n_particles_max=maxval(n_particles_per_group)
real*8,parameter                            :: survival_threshold=0.33
!> parameters for generating spectra
integer,parameter                           :: n_spectra=2
integer,parameter                           :: n_lines_per_spectrum=13
real*8,dimension(n_spectra),parameter       :: min_wlen=(/3.0d-6,2.5d-7/)
real*8,dimension(n_spectra),parameter       :: max_wlen=(/3.5d-6,4.2d-7/)
!> parameters for generating shadowed points
integer,parameter                           :: n_shadowed_per_particle=7
real*8,dimension(2),parameter               :: length_shadowed=(/2.d-1,7.d0/)
!> variables for generating synchrotron lights
type(synchrotron_light_vertices)            :: vertex_sol
type(particle_sim),dimension(n_times_sol)   :: sims_particles
integer                                     :: n_particles_RE_max
integer,dimension(n_times_sol)              :: n_active_vertices_sol
integer,dimension(n_groups_max,n_times_sol) :: n_active_particles_sol
integer,dimension(n_particles_max,n_groups_max,n_times_sol) :: active_particle_ids_sol
real*8,dimension(n_times_sol)               :: time_vector_sol
real*8,dimension(n_x,n_particles_max*n_groups_max,n_times_sol) :: x_cart_sol
real*8,dimension(n_properties,n_particles_max*n_groups_max,n_times_sol) :: properties_sol
!> variables for generating spectra
type(spectrum_rng_uniform)                  :: spectrum
!> variables for generating shadowed points
real*8,dimension(:,:,:,:),allocatable       :: x_shadowed

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
  call test_init_synchrotron_lights_from_particles
  call test_synchrotron_irradiance_directional_func
  write(*,'(/A)') "  ... tearing-down: synchrotron light vertices tests"
  call teardown
end subroutine run_fruit_synchrotron_light_vertices

!> Set-up and tear-down procedures------------------------------------
!> allocate and initialise the unit test features
subroutine setup()
  use mod_rng,                              only: type_rng
  use mod_pcg32_rng,                        only: pcg32_rng
  use mod_particle_types,                   only: particle_gc_vpar_id
  use mod_particle_types,                   only: particle_kinetic_id
  use mod_particle_types,                   only: particle_kinetic_relativistic_id
  use mod_gnu_rng,                          only: gnu_rng_interval
  use mod_common_test_tools,       only: omp_initialize_rngs
  use mod_particle_common_test_tools,       only: sim_time_interval
  use mod_particle_common_test_tools,       only: allocate_one_particle_list_type 
  use mod_particle_common_test_tools,       only: fill_groups,fill_mass_RE
  use mod_particle_common_test_tools,       only: fill_particles_tokamak
  use mod_particle_common_test_tools,       only: invalidate_particles
  use mod_particle_common_test_tools,       only: obtain_active_particle_ids
  !$ use omp_lib
  implicit none
  !> variables
  integer,dimension(n_groups_max,n_times_sol),parameter :: particle_types=&
  reshape((/particle_kinetic_relativistic_id,particle_gc_vpar_id,&
  particle_kinetic_relativistic_id,particle_kinetic_id,&
  particle_kinetic_relativistic_id,0/),shape(particle_types))
  integer :: ii,jj,ifail,n_particles_RE_max_loc,n_threads
  class(type_rng),dimension(:),allocatable :: rngs
  !> initialisation
  vertex_sol%n_property_vertex = n_properties; ifail = 0; 
  n_particles_RE_max = 0; n_active_particles_sol = 0;
  n_threads = 1
  !$ n_threads = omp_get_max_threads()
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
    !> set to zero if particles are not kinetic relativistic
    where(particle_types(:,ii).ne.particle_kinetic_relativistic_id) 
      n_active_particles_sol(:,ii)=0
    endwhere
    n_particles_RE_max_loc = 0
    do jj=1,n_groups_per_sim(ii)
      if(particle_types(jj,ii).ne.particle_kinetic_relativistic_id) cycle
      n_particles_RE_max_loc = n_particles_RE_max_loc + n_particles_per_group(jj,ii)
    enddo
    n_particles_RE_max = max(n_particles_RE_max,n_particles_RE_max_loc)
  enddo
  !> initialise positions and properties tables
  call compute_synch_x_properties_ana()

  !> initialise monte-carlo spectra
  allocate(pcg32_rng::rngs(n_threads))
  call omp_initialize_rngs(n_lines_per_spectrum,n_threads,rngs)
  spectrum = spectrum_rng_uniform(n_lines_per_spectrum,n_spectra,min_wlen,max_wlen)
  call spectrum%generate_spectrum(rngs)
  deallocate(rngs)

  !> generate shadowed points positions
  allocate(x_shadowed(n_x,n_shadowed_per_particle,n_particles_RE_max,n_times_sol))
  call compute_x_shadowed_particles() 
end subroutine setup

!> destroy all test features
subroutine teardown()
  implicit none
  call vertex_sol%deallocate_vertices; call spectrum%deallocate_spectrum;
  deallocate(x_shadowed);
end subroutine teardown

!> Tests -------------------------------------------------------------
!> test the synchrotron light irradiance and directional functions
subroutine test_synchrotron_irradiance_directional_func()
  use mod_assert_equals_tools, only: assert_equals_rel_error
  implicit none
  !> variables
  integer :: ii,jj,kk
  integer,dimension(n_times_sol) :: n_particles_time
  real*8,dimension(n_x,n_particles_RE_max,n_times_sol)          :: x_cart_loc
  real*8,dimension(n_properties,n_particles_RE_max,n_times_sol) :: properties_loc
  real*8,dimension(spectrum%n_points,spectrum%n_spectra)        :: dir_fun,irradiance
  real*8,dimension(spectrum%n_points,spectrum%n_spectra)        :: dir_fun_sol,irradiance_sol
  !> initialise the synchrotron lights (for safety)
   x_cart_loc = 0.d0; properties_loc = 0.d0; n_particles_time = 0;
  do ii=1,n_times_sol
    n_particles_time(ii) = sum(n_active_particles_sol(:,ii))
    x_cart_loc(:,1:n_particles_time(ii),ii) = x_cart_sol(:,1:n_particles_time(ii),ii)
    properties_loc(:,1:n_particles_time(ii),ii) = properties_sol(:,1:n_particles_time(ii),ii)   
  enddo 
  call vertex_sol%init_lights_from_particles(n_times_sol,sims_particles)

  !> compute the irradiance and directionality function
  do kk=1,n_times_sol
    do jj=1,n_particles_time(kk)
      do ii=1,n_shadowed_per_particle
        call compute_synch_directionality_irradiance(x_shadowed(:,ii,jj,kk),&
        x_cart_loc(:,jj,kk),properties_loc(:,jj,kk),dir_fun_sol,irradiance_sol)
        call vertex_sol%directionality_funct(spectrum,kk,jj,x_shadowed(:,ii,jj,kk),dir_fun)
        call vertex_sol%spectral_irradiance(spectrum,kk,jj,x_shadowed(:,ii,jj,kk),irradiance)
        !> set to 1 values which are too small 
        where(abs(irradiance).lt.exclusion_values)     
          irradiance = 1.d0; irradiance_sol = 1.d0;
        endwhere
        where(abs(dir_fun).lt.exclusion_values)        
          dir_fun = 1.d0; dir_fun_sol = 1.d0;
        endwhere
        !> check the solution via relative error
        call assert_equals_rel_error(spectrum%n_points,spectrum%n_spectra,&
        dir_fun,dir_fun_sol,tol2_real8,&
        "Error synchrotron directionality function: directionality function mismatch!")
        call assert_equals_rel_error(spectrum%n_points,spectrum%n_spectra,&
        irradiance,irradiance_sol,tol2_real8,&
        "Error synchrotron irradiance: irradiance mismatch!")
      enddo
    enddo
  enddo
end subroutine test_synchrotron_irradiance_directional_func

!> test the initialisation of synchrotron lights from particles
subroutine test_init_synchrotron_lights_from_particles()
  use mod_assert_equals_tools, only: assert_equals_rel_error
  implicit none
  !> variables
  integer,parameter                                             :: n_sync_fail=1
  integer                                                       :: ii,n_particles_time
  real*8,dimension(n_x,n_particles_RE_max,n_times_sol)          :: x_cart_loc
  real*8,dimension(n_properties,n_particles_RE_max,n_times_sol) :: properties_loc
  !> test light structure initialisation with given size
  call vertex_sol%init_lights_from_particles(n_times_sol,&
  sims_particles,n_particles_max*n_groups_max)
  call assert_equals_rel_error(n_x,n_particles_max*n_groups_max,n_times_sol,&
  vertex_sol%x,x_cart_sol,tol_real8,&
  "Error init synchrotron lights from particles set n lights large: positions errors too large!")
  call assert_equals_rel_error(n_properties,n_particles_max*n_groups_max,&
  n_times_sol,vertex_sol%properties,properties_sol,tol_real8,&
  "Error init synchrotron lights from particles set n lights large: properties errors too large!")

  !> copy valued of x and properties
  x_cart_loc = 0.d0; properties_loc = 0.d0;
  do ii=1,n_times_sol
    n_particles_time = sum(n_active_particles_sol(:,ii))
    x_cart_loc(:,1:n_particles_time,ii) = x_cart_sol(:,1:n_particles_time,ii)
    properties_loc(:,1:n_particles_time,ii) = properties_sol(:,1:n_particles_time,ii)   
  enddo
  !> test vertices initialisation no inputs
  call vertex_sol%init_lights_from_particles(n_times_sol,sims_particles)
  call assert_equals(shape(vertex_sol%x),shape(x_cart_loc),3,&
  "Error init synchrotron lights from particles: positions shape mismatch!")
  call assert_equals(shape(vertex_sol%properties),shape(properties_loc),3,&
  "Error init synchrotron lights from particles: properties shape mismatch!")
  call assert_equals_rel_error(n_x,n_particles_RE_max,n_times_sol,&
  vertex_sol%x,x_cart_loc,tol_real8,&
  "Error init synchrotron lights from particles: positions errors too large!")
  call assert_equals_rel_error(n_properties,n_particles_RE_max,&
  n_times_sol,vertex_sol%properties,properties_loc,tol_real8,&
  "Error init synchrotron lights from particles: properties errors too large!")

  !> test vertices initilisation with too small number of lights
  call vertex_sol%init_lights_from_particles(n_times_sol,sims_particles,n_sync_fail)
  call assert_equals(shape(vertex_sol%x),shape(x_cart_loc),3,&
  "Error init synchrotron lights from particles set n lights small: positions shape mismatch!")
  call assert_equals(shape(vertex_sol%properties),shape(properties_loc),3,&
  "Error init synchrotron lights from particles set n lights small: properties shape mismatch!")
  call assert_equals_rel_error(n_x,n_particles_RE_max,n_times_sol,&
  vertex_sol%x,x_cart_loc,tol_real8,&
  "Error init synchrotron lights from particles set n lights small: positions errors too large!")
  call assert_equals_rel_error(n_properties,n_particles_RE_max,&
  n_times_sol,vertex_sol%properties,properties_loc,tol_real8,&
  "Error init synchrotron lights from particles set n lights small: properties errors too large!")
end subroutine test_init_synchrotron_lights_from_particles 

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
  "Error fill synchrotron lights from particles: positions errors too large!")
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
!> compute the positions of points shadowed by particle lights.
!> The shadowed positions are taken within the synchrotron emission
!> cone of each relativistic particle. It is assumed that 
!> the half cone of radiation emission is 1/(2.0*gamma) 
!> gamma = sqrt(1+(p/m*c)**2) = relativistic factor
subroutine compute_x_shadowed_particles()
  use mod_gnu_rng,  only: gnu_rng_interval
  use mod_sampling, only: sample_uniform_cone
  implicit none
  !> variables 
  integer             :: ii,jj,kk,n_particles_time
  real*8              :: cos_half_angle
  real*8,dimension(n_x) :: p_dir,x_part,rng
  !> generate shadowed particle positions
  do kk=1,n_times_sol
    n_particles_time = sum(n_active_particles_sol(:,kk))
    do jj=1,n_particles_time
      x_part = x_cart_sol(:,jj,kk)
      p_dir = properties_sol(1:3,jj,kk)
      cos_half_angle = cos(1.d0/(2.d0*properties_sol(11,jj,kk)))
      do ii=1,n_shadowed_per_particle
        call random_number(rng)
        x_shadowed(:,ii,jj,kk) = sample_uniform_cone(cos_half_angle,&
        rng,p_dir,x_part,length_shadowed)
      enddo
    enddo
  enddo 
end subroutine compute_x_shadowed_particles

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

!> compute the synchrotron lights directionaly function and irradiance
subroutine compute_synch_directionality_irradiance(x_illum,&
x_light,property,dir_func,irradiance)
  use constants,                 only: TWOPI,PI,EL_CHG,EPS_ZERO,SPEED_OF_LIGHT
  use mod_boost_besselk,         only: besselk
  use mod_coordinate_transforms, only: cartesian_to_spherical_latitude
  implicit none
  !> inputs
  real*8,dimension(n_x),intent(in)          :: x_illum,x_light
  real*8,dimension(n_properties),intent(in) :: property
  !> outputs
  real*8,dimension(spectrum%n_points,spectrum%n_spectra),intent(out) :: dir_func,irradiance
  !> variables
  integer               :: ii,jj
  real*8                :: onethird=1.d0/3.d0
  real*8                :: twothird=2.d0/3.d0
  real*8                :: z,zeta,besselk13,besselk23
  real*8                :: one_z2,z_z3,factor,factor_2
  real*8,dimension(n_x) :: rpsichi

  !> compute the spherical coordinate variables
  rpsichi = cartesian_to_spherical_latitude(x_illum,x_light,&
  property(1:3),property(4:6),property(7:9))
  !> compute the particle dependent variables
  z = (property(11)*rpsichi(3))/sqrt(1.d0+(property(11)*rpsichi(2))**2.d0)
  one_z2 = 5.d-1*(1.d0+z*z); z_z3 = 1.5d0*(z+((z**3.d0)/3.d0))
  factor_2 = ((property(11)*rpsichi(2))**2.d0)/(1.d0+(property(11)*rpsichi(2))**2.d0)
  factor = (1.d0+(property(11)*rpsichi(2))**2.d0)**2.d0
  factor = (factor*SPEED_OF_LIGHT*(EL_CHG**2.d0))/&
           (sqrt(3.d0)*EPS_ZERO*property(12)*(property(11)**4.d0))
  !> compute the irradiance
  !> compute the directionaly function
  do ii=1,spectrum%n_spectra
    do jj=1,spectrum%n_points
      zeta = TWOPI*(((1.d0/property(11)**2.d0)+rpsichi(2)**2.d0)**1.5d0)/&
      (3.d0*spectrum%points(jj,ii)*property(12))
      call besselk(onethird,zeta,besselk13)
      call besselk(twothird,zeta,besselk23)
      irradiance(jj,ii) = factor*((factor_2-one_z2)*besselk13*cos(zeta*z_z3)+&
      besselk23*z*sin(zeta*z_z3))/(spectrum%points(jj,ii)**4.d0)
      dir_func(jj,ii) = irradiance(jj,ii)/property(13)
    enddo
  enddo
end subroutine compute_synch_directionality_irradiance

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
