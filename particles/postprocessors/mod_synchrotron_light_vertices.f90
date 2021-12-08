!> the mod_synchrotron_light implements variable and
!> procedures defining a synchrotron light source
module mod_synchrotron_light_vertices
use mod_light_vertices, only: light_vertices
implicit none

private
public :: synchrotron_light_vertices

!> Variables ---------------------------------------
real*8,parameter :: onethird=1.d0/3.d0
real*8,parameter :: twothirds=2.d0/3.d0
real*8,parameter :: sqrt3=sqrt(3.d0)
type,extends(light_vertices) :: synchrotron_light_vertices
  contains
  procedure,pass(sync_light_vertices) :: init_lights_from_particles =>&
                                           init_synchrotron_light_from_particles
  procedure,pass(sync_light_vertices) :: directionality_funct => &
                                           synchrotron_directionality_funct
  procedure,pass(sync_light_vertices) :: spectral_irradiance => &
                                           synchrotron_spectral_irradiance
end type synchrotron_light_vertices
!> Interfaces --------------------------------------

contains

!> Procedures --------------------------------------
!> init_synchrotron_light_from_particles initilises compute the properties
!> of the synchrotron light for each particle and stores them 
!> in the proprerties array
!> inputs:
!>   sync_lights:      (synchrotron_light_vertices) empty synchrotron lights
!>   n_times:          (integer) number of simulation times
!>   sims_particles:   (particle_sim)(n_times) array of particle simulations
!>   n_sync_lights_in: (integer)(optional) number of requested synchrotron lights
!> outputs:
!>   sync_lights: (synchrotron_light_vertices) initialised synchrotron lights
subroutine init_synchrotron_lights_from_particles(sync_lights,n_times,&
sims_particles,n_sync_lights_in)
  use mod_particle_types,        only: particle_kinetic_relativistic_id
  use mod_particle_sim,          only: particle_sim
  implcit none
  !> inputs-outputs
  type(synchrotron_light_vertices),intent(inout) :: sync_lights
  !> inputs
  type(partice_sim),dimension(n_times),intent(in) :: sims_particles
  integer,intent(in)                              :: n_times
  integer,intent(in),optional                     :: n_sync_lights_in
  !> variables
  integer :: n_sync_lights,n_groups_max
  integer,dimension(n_times) :: n_groups,n_particle_relativistics
  integer,dimension(:,:)     :: n_particles,particle_types

  sync_lights%n_property_vertex = 13 !< set number of synchrotron vertex properties
  !> initialise time vector
  call sync_lights%allocate_time_vector(n_times)
  call sync_lights%fill_time_vector(sims_particles)
  !> allocate extract number of particles and particles type
  call sync_lights%extract_n_groups_all_particle_sims(sims_particles,n_groups)
  n_groups_max = maxval(n_groups)
  allocate(n_particles(n_groups_max,sync_lights%n_times)) 
  allocate(particle_types(n_groups_max,sync_lights%n_times))
  call sync_lights%extract_n_particles_all_particle_sims(sims_particles,&
  n_groups_max,n_particles)
  call sync_lights%extract_particle_types_all_particle_sims(sims_particles,&
  n_groups_max,particle_types)
  !> compute the number of relativistic particles per each time
  do ii=1,sync_lights%n_times
    n_particle_relativistics(ii) = sum(n_particles(:,ii),&
    mask=particle_types(:,ii)==particle_kinetic_relativistic_id)
  enddo
  n_sync_lights = maxval(n_particle_relativistics)
  if(present(n_sync_lights_in)) then
    if(n_sync_lights.lt.n_sync_lights_in) then
      write(*,*) "Error initialise synchrotron lights from particles"
      write(*,*) "Requested number of lights < number of particles,use: ",n_sync_lights
    endif
    n_sync_lights = max(n_sync_lights,n_sync_lights_in)
  endif
  !> allocate vertices
  call sync_lights%allocate_x_properties(n_sync_lights)

  !> if no omp fill the synchrotron lights sequentially
  call fill_synchrotron_lights_from_particles_serial(sync_lights,sims_particles,&
  n_groups_max,n_groups,n_particles)

  !> cleanup 
  deallocate(n_particles); deallocate(particle_types);
end subroutine init_synchrotron_lights_from_particles

!> synchrotron_directionality_funct computes the directionaliy function
!> for synchrotron lights which is the full angular-spectral distribution
!> divided by the total synchrotron radiation (L. Carbajal, PPCF, 2017)
!> inputs:
!>   sync_lights: (synchrotron_light vertices) synchrotron light sources
!>   spectra:     (spectrum_base) spectral intervals and integrators
!>   time_id:     (integer) the time index
!>   light_id:    (integer) the light index
!>   x_shaded:    (real8)(3) shaded point position in cartesian coord
!> outputs: 
!>   sync_lights: (synchrotron_light vertices) synchrotron light sources
!>   light_dstb:  (real*8)(n_points,n_intervals) synchrotron full spectral
!>                angular distribution per unit of total power towards
!>                the shaded point x_shaded
subroutine synchrotron_directionality_funct(sync_lights,spectra,time_id,&
light_id,x_shaded,light_dstb)
  use constants,               only: PI,SPEED_OF_LIGHT
  use mod_boost_bessek,        only: besselk
  use mod_coordinate_transform,only: cartesian_to_spherical_latitude
  use mod_spectra,             only: spectra_base
  implicit none
  !> inputs:
  type(synchrotron_light_vertices),intent(inout) :: sync_lights
  class(spectrum_base),intent(in)                :: spectra
  integer,intent(in)                             :: time_id,light_id
  real*8,dimension(synch_lights%n_x),intent(in)  :: x_shaded
  !> outputs:
  real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_dstb
  !> variables
  real*8,dimension(3) :: rpsichi !< spherical coordinates
  integer :: ii,jj
  real*8,dimension(sync_lights%n_property_vertex) :: light_properties
  real*8  :: zeta,one_over_gamma,z_value,factor_1,factor_2,z_cos

  !> compute the spherical coordinates of the light-point ray
  light_properties = sync_lights%properties(:,light_id,time_id)
  rpsichi = cartesian_to_sphetical_latitude(x_shaded,sync_lights%x(:,light_id,time_id),&
  light_properties(1:3,light_properties(4:6),light_properties(7:9))
  !> compute the factors and the value of z
  one_over_gamma = 1.d0/light_properties(11) !< 1/gamma
  factor_2 = light_properties(11)*light_properties(11)*rpsichi(2)*rpsichi(2) !< gamma**2 * psi**2
  factor_1 = 1.d0+factor_2 !< 1 + gamma**2 * psi**2
  factor_2 = factor_2/factor_1 !< (gamma**2 * psi**2) / (1 + gamma**2 * psi**2)
  !> z = gamma*chi / sqrt(1 + gamma**2 * psi**2)
  z_value = (light_properties(11)*rpsichi(3))/sqrt(factor_1) 
  z_cos = 1.5d0*z_value*(1.d0+z_value*z_value/3.d0) !< z_cos = (3/2)*z*(1 + (z**2)/3)
  z_value = 5.d-1**(1.d0+z_value*z_value) !< z = 0.5*(1+z**2)
  !> I = Power*( 1 + gamma**2 * psi**2)**2 / Power_tot = 
  !> (6*PI / (sqrt(3)) * beta**4 * gamma**8 * kappa**3 )*( 1 + gamma**2 * psi**2)**2
  !> beta = v/c; kappa = (|q|/(gamma*mass*v**3))||v X (E + v X B)||
  !> Power_tot = (q**2/(6*PI*eps0*c**3))*gamma**4 * v**4 * kappa**2
  factor_1 = factor_1*factor_1*((6.d0*PI)/(sqrt3*(properties(10)**4.d0)*&
  (properties(11)**8.d0)*(properties(12)**3.d0)))
  !$omp parallel do default(private) shared(spectra,light_dstb,&
  !$omp rpsichi,one_over_gamma,light_properties,factor_1,factor_2,&
  !$omp z_value,z_cos) collapse(2)
  do ii=1,spectra%n_spectra
    do jj=1,spectra%n_points
      !> zeta = (2*PI*(1/gamma**2 + psi**2)**(3/2))/(3*kappa*lambda)
      zeta = 2.d0*PI*((one_over_gamma*one_over_gamma+rpsichi(2)*rpsichi(2))**1.5d0)/&
      (3.d0*spectra%points(jj,ii)*light_properties(12))
      !> funct = K_1/3(zeta)*cos(zeta*z_cos)*(((gamma**2 * psi**2)/(1 + gamma**2 * psi**2)) -
      !>  0.5*(1+z**2)) + K_(2/3)(zeta)*sin(zeta*z_cos)
      light_dstb(jj,ii) = besselk(onethird,zeta)*cos(zeta*z_cos)*(factor_2-z_value)+&
      besselk(twothirds,zeta)*sin(zeta*z_cos)
      !> SAD/P_tot = I*funct/lambda**4
      light_dstb(jj,ii) = fact_1*light_dstb(jj,ii)/&
      (spectra%points(jj,ii)*spectra%points(jj,ii)*spectra%points(jj,ii)*spectra%points(jj,ii))
    enddo
  enddo
  !$end omp parallel do
end subroutine synchrotron_directionality_funct

!> synchrotron_spectral_irradiance computes the full power spectral angular
!> angular distribution for synchrotron lights which (L. Carbajal, PPCF, 2017)
!> emitted towards the shaded point x_shaded
!> inputs:
!>   sync_lights: (synchrotron_light vertices) synchrotron light sources
!>   spectra:     (spectrum_base) spectral intervals and integrators
!>   time_id:     (integer) the time index
!>   light_id:    (integer) the light index
!>   x_shaded:    (real8)(3) shaded point position in cartesian coord
!> outputs: 
!>   sync_lights:            (synchrotron_light vertices) synchrotron light sources
!>   light_spec_irradiance:  (real*8)(n_points,n_intervals) synchrotron full spectral
!>                           angular distribution per unit of total power at the
!>                           at the shaded point x_shaded
subroutine synchrotron_spectral_irradiance(sync_lights,spectra,time_id,&
light_id,x_shaded,light_spec_irradiance)
  use mod_spectra, only: spectra_base
  implicit none
  !> inputs:
  type(synchrotron_light_vertices),intent(inout) :: sync_lights
  class(spectrum_base),intent(in)                :: spectra
  integer,intent(in)                             :: time_id,light_id
  real*8,dimension(synch_lights%n_x),intent(in)  :: x_shaded
  !> outputs:
  real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_spec_irradiance

  !> compute the directionality function
  call sync_lights%directionality_funct(sync_lights,spectra,time_id,&
  light_id,x_shaded,light_spec_irradiance)
  !> multiply the directionality function by the total synchrotron power
  light_spec_irradiance = light_spec_irradiance*sync_lights%properties(13,light_id,time_id)
end subroutine synchrotron_spectral_irradiance


!> Tools ------------------------------------------
!> fill_synchrotron_lights_from_particles_serial fill the
!> x and properties array of synchrotron light from
!> particle list sequentially
!> inputs:
!>   sync_lights:      (synchrotron_light_vertices) empty synchrotron lights
!>   sims_particles:   (particle_sim)(n_times) array of particle simulations
!>   n_groups_max:     (integer) maximum size of groups
!>   n_groups:         (integer)(n_times) size of each group 
!>   n_particles:      (integer)(n_group_max,n_times) size of the particle list
!>                     for each group and for each time
!> outputs:
!>   sync_lights: (synchrotron_light_vertices) initialised synchrotron lights
subroutine fill_synchrotron_lights_from_particles_serial(&
sync_lights,sims_particles,n_groups_max,n_groups,n_particles)
  use mod_particle_sim,   only: particle_sim
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  implicit none
  !> inputs-outputs
  class(synchrotron_light_vertices),intent(inout) :: sync_lights
  !> inputs:
  type(particle_sim),dimension(sync_lights%n_times),intent(in)   :: sims_particles
  integer,intent(in)                                             :: n_groups_max
  integer,dimension(sync_lights%n_times),intent(in)              :: n_groups
  integer,dimension(n_groups_max,sync_lights%n_times),intent(in) :: n_particles
  !> variables
  integer :: ii,jj,kk,pp
  real*8  :: psi,U
  real*8,dimension(sync_lights%n_x) :: E_field,B_field

  !> compute synchrotron light properties from particle simulations
  do ii=1,sync_lights%n_times
    pp = 0
    do jj=1,n_groups(ii)
      select type (p_list=>sims_particles(ii)%groups(jj)%particles)
        type is (particle_kinetic_relativistic)
        do kk=1,n_particles(jj,ii)
          if(p_list(kk)%i_elm.lt.1) cycle !< skip invalid particle
          pp = pp + 1 !< increase the light index
          call sync_lights%store_light_x_from_particle_id(pp,ii,p_list(kk)) !< store position
          !> compute E,B fields
          call sims_particles(ii)%fields%calc_EBpsiU(sync_lights%times(ii),&
          p_list(kk)%i_elm,p_list(kk)%st,p_list(kk)%x(3),E_field,B_field,psi,U)
          !> compute synchrotron light properties
          call compute_synchrotron_light_properties(&
          sync_lights%n_x,sync_lights%n_property_vertex,p_list(kk),&
          sims_particles(ii)%groups(jj)%mass,&
          vector_cylindrical_to_cartesian(p_list(kk)%x(3),E_field),&
          vector_cylindrical_to_cartesian(p_list(kk)%x(3),B_field),&
          synch_lights%properties(:,pp,ii))
        enddo
      end select
    enddo
  enddo
end subroutine fill_synchrotron_lights_from_particles_serial

!> fill_synchrotron_lights_from_particles_openmp fill the
!> x and properties array of synchrotron light from
!> particle list parallelised by openmp
!> inputs:
!>   sync_lights:      (synchrotron_light_vertices) empty synchrotron lights
!>   sims_particles:   (particle_sim)(n_times) array of particle simulations
!>   n_groups_max:     (integer) maximum size of groups
!>   n_groups:         (integer)(n_times) size of each group 
!>   n_particles:      (integer)(n_group_max,n_times) size of the particle list
!>                     for each group and for each time
!> outputs:
!>   sync_lights: (synchrotron_light_vertices) initialised synchrotron lights
subroutine fill_synchrotron_lights_from_particles_openmp(&
sync_lights,sims_particles,n_groups_max,n_groups,n_particles)
  use mod_particle_sim,   only: particle_sim
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  implicit none
  !> inputs-outputs
  class(synchrotron_light_vertices),intent(inout) :: sync_lights
  !> inputs:
  type(particle_sim),dimension(sync_lights%n_times),intent(in)   :: sims_particles
  integer,intent(in)                                             :: n_groups_max
  integer,dimension(sync_lights%n_times),intent(in)              :: n_groups
  integer,dimension(n_groups_max,sync_lights%n_times),intent(in) :: n_particles
  !TODO parallel fill synchrotron lights sources
end subroutine fill_synchrotron_lights_from_particles_openmp

!> compute_synchrotron_light_properties computes the
!> synchrotron radiation properties from a
!> kinetic relativistic particle. Variables
!> inputs:
!>   field_size:    (integer) size of the field vectors
!>   property_size: (integer) size of the property vector
!>   particle_in:   (particle_kinetic_relativistic) jorek particle
!>   mass:          (real8) mass of the particle
!>   E_field_cart:  (real8)(field_size) electric field at particle position
!>                  in the cartesian reference system
!>   B_field_cart:  (real8)(field_size) magnetic field at particle position
!>                  in the cartesian reference system
!> outputs:
!>   sync_properties: (real8)(property_size) synchrotron radiation properties
!>                    1:3 -> component of the velocity direction (cartesian)
!>                        -> T = v/||v||
!>                    4:6 -> second orthonormal basis cartesian coordinates
!>                        -> N = E + v X B - v*E
!>                    7:9 -> components of the third orthonormal basis (cartesian)
!>                        -> B = T X N
!>                    10  -> beta -> velocity/speed of light = v/c
!>                    11  -> relativistic factor gamma = sqrt(1+(p/(mass*c))**2)
!>                    12  -> orbit curvature (L. Carbakal, PPCF, 2017)
!>                        -> kappa = (|q|/(gamma*mass*v**3))||v X (E + v X B)||
!>                    13  -> total radiation power (L. Carbajal, PPCF, 2017)i
!>                        -> P_tot = (q**2/(6*PI*eps0*c**3))*gamma**4 * v**4 * kappa**2
subroutine compute_synchrotron_light_properties(field_size,property_size,particle_in,&
,mass,E_field_cart,B_field_cart,snchrt_properties)
  use constants,                 only: PI,EPS_ZERO,EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
  use mod_math_operators,        only: cross_product
  use mod_coordinate_transforms, only: vectors_to_orthonormal_basis
  use mod_particle_types,        only: particle_kinetic_relativistic
  implicit none
  !> inputs
  type(particle_kinetic_relativistic),intent(in) :: particle_in
  integer                                        :: field_size,property_size
  real*8,intent(in)                              :: mass
  real*8,dimension(field_size),intent(in)        :: B_field,E_field
  !> outputs
  real*8,dimension(property_size),intent(out) :: sync_properties
  !> variables
  real*8 :: velocity
  real*8,dimension(3) :: vector_1d_3,vector_1d_3_2,vector_1d_3_3

  !> compute velocity, velocity direction and relativistic factor
  velocity =  sqrt(particle_in%x(1)*particle_in%x(1)+&
              particle_in%x(2)*particle_in%x(2)+&
              particle_in%x(3)*particle_in%x(3)) 
  sync_properties(1:3) = particle_in%p/velocity
  sync_properties(10)   = velocity/SPEEF_OF_LIGHT
  sync_properties(11)   = sqrt(1.d0 + (sync_properties(10)*sync_properties(10))/(mass*mass))
  sync_properties(10)   = sync_property(10)/(mass*sync_properties(11))
  !> compute orbit curvature
  synch_properties(4:6) = E_field+cross_product(particle_in%p/(mass*sync_properties(8)),B_field)
  vector_1d_3(4:6) = cross_product(sync_properties(1:3),synch_properties(4:6))
  sync_properties(9) = (abs(real(particle_in%q,kind=8))*EL_CHG*&
                       sqrt(vector_1d_3(1)*vector_1d_3(1)+&
                       vector_1d_3(2)*vector_1d_3(2)+&
                       vector_1d_3(3)*vector_1d_3(3)))/&
                       (sync_properties(11)*mass*ATOMI_MASS_UNIT*&
                       sync_properties(10)*sync_properties(7)*&
                       SPEED_OF_LIGHT*SPEED_OF_LIGHT)
  !> compute total synchrotron power
  sync_properties(13) = (real(particle_in%q*particle_in%q,kind=8)*&
                       EL_CHG*EL_CHG*SPEED_OF_LIGHT*sync_properties(10)*&
                       sync_properties(10)*sync_properties(10)*&
                       sync_properties(10)*sync_properties(11)*&
                       sync_properties(11)*sync_properties(11)*&
                       sync_properties(11)*sync_properties(12)*&
                       sync_properties(12))/(6.d0*PI*EPS_ZERO)
  !> compute orthonormal basis
end subroutine compute_synchrotron_light_properties

!> construct and store the orthonormal basis
call vectors_to_orthonormal_basis(sync_properties(1:3),sync_properties(4:6),&
vector_1d_3,vector_1d_3_2,vector_1d_3_3)
sync_properties(1:3) = vector_1d_3; sync_properties(4:6) = vector_1d_3_2
sync_properties(7:9) = vector_1d_3_3

!>-------------------------------------------------
end module mod_synchrotron_light_vertices
