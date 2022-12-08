!> the mod_synchrotron_light implements variables and
!> procedures defining a synchrotron light source
module mod_synchrotron_light_vertices
use mod_light_vertices, only: light_vertices
implicit none

private
public :: synchrotron_light_vertices

!> Variables ---------------------------------------
real*8,parameter  :: onethird=1.d0/3.d0
real*8,parameter  :: twothirds=2.d0/3.d0
real*8,parameter  :: sqrt3=sqrt(3.d0)
type,extends(light_vertices) :: synchrotron_light_vertices
  contains
  procedure,pass(light_vert) :: directionality_funct => &
                                synchrotron_directionality_funct
  procedure,pass(light_vert) :: spectral_irradiance => &
                                synchrotron_spectral_irradiance
  procedure,pass(light_vert) :: compute_mhd_fields => &
                                compute_synchrotron_mhd_fields
  procedure,pass(light_vert) :: compute_light_properties => &
                                compute_synchrotron_light_properties
  procedure,pass(light_vert) :: setup_light_class => &
                                setup_synchrotron_light_class
  procedure,nopass           :: check_x_shaded_in_emission_zone => &
                                check_shaded_x_in_synchrotron_cone
end type synchrotron_light_vertices
!> Interfaces --------------------------------------

contains

!> Procedures --------------------------------------

!> synchrotron_directionality_funct computes the directionaliy function
!> for synchrotron lights which is the full angular-spectral distribution
!> divided by the total synchrotron radiation (L. Carbajal, PPCF, 2017)
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> WARNING: the current implementation is prone to roundoff
!> error due to large differences in the terms exponents 
!> TODO: find a more performant normalisation for the 
!> TODO: full spectral-angular power distribution
!> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> inputs:
!>   light_vert: (synchrotron_light vertices) synchrotron light sources
!>   spectra:     (spectrum_base) spectral intervals and integrators
!>   time_id:     (integer) the time index
!>   light_id:    (integer) the light index
!>   x_shaded:    (real8)(3) shaded point position in cartesian coord
!> outputs: 
!>   light_vert: (synchrotron_light vertices) synchrotron light sources
!>   spectra:     (spectrum_base) spectral intervals and integrators
!>   light_dstb:  (real*8)(n_points,n_intervals) synchrotron full spectral
!>                angular distribution per unit of total power towards
!>                the shaded point x_shaded
subroutine synchrotron_directionality_funct(light_vert,spectra,time_id,&
light_id,x_shaded,light_dstb)
  use constants,                only: TWOPI,SPEED_OF_LIGHT
  use mod_besselk,        only: f_besselk
  use mod_coordinate_transforms,only: cartesian_to_spherical_latitude
  use mod_spectra,              only: spectrum_base
  !$ use omp_lib
  implicit none
  !> inputs-outputs:
  class(synchrotron_light_vertices),intent(inout) :: light_vert
  class(spectrum_base),intent(inout)              :: spectra
  !> inputs:
  integer,intent(in)                           :: time_id,light_id
  real*8,dimension(light_vert%n_x),intent(in)  :: x_shaded
  !> outputs:
  real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_dstb
  !> variables
  logical :: in_parallel
  integer :: ii,jj 
  integer,dimension(0) :: int_param
  real*8  :: zeta,one_over_gamma,z_value,z2_value,factor_1,factor_2,z_cos
  real*8,dimension(light_vert%n_x) :: rpsichi !< spherical coordinates
  real*8,dimension(light_vert%n_property_vertex) :: light_properties

  !> initialisations
  in_parallel = .false.; light_dstb = 0d0;
  !$ in_parallel = omp_in_parallel()
  !> check if the shaded point is in the synchrotron emission cone
  light_properties = light_vert%properties(:,light_id,time_id)
  if(.not.light_vert%check_x_shaded_in_emission_zone(light_vert%n_x,x_shaded,&
  light_vert%x(:,light_id,time_id),0,4,int_param,[light_properties(1),&
  light_properties(2),light_properties(3),light_properties(11)])) return
  !> compute the spherical coordinates of the light-point ray
  rpsichi = cartesian_to_spherical_latitude(x_shaded,light_vert%x(:,light_id,time_id),&
  light_properties(1:3),light_properties(4:6),light_properties(7:9))
  !> compute the factors and the value of z
  one_over_gamma = 1.d0/(light_properties(11)**2) !< 1/gamma
  factor_2 = (light_properties(11)*rpsichi(2))**2 !< gamma**2 * psi**2
  factor_1 = 1.d0+factor_2 !< 1 + gamma**2 * psi**2
  factor_2 = factor_2/factor_1 !< (gamma**2 * psi**2) / (1 + gamma**2 * psi**2)
  !> z = gamma*chi / sqrt(1 + gamma**2 * psi**2)
  z_value = (light_properties(11)*rpsichi(3))/sqrt(factor_1) !< chi*gamma/sqrt(1+gamma**2 * psi**2)
  z_cos = 1.5d0*z_value*(1.d0+((z_value**2)/3.d0)) !< z_cos = (3/2)*z*(1 + (z**2)/3)
  z2_value = 5.d-1*(1.d0+z_value**2) !< z = 0.5*(1+z**2)
  !> I = Power*( 1 + gamma**2 * psi**2)**2 / Power_tot = 
  !> (6*PI / (sqrt(3)) * beta**4 * gamma**8 * kappa**3 )*( 1 + gamma**2 * psi**2)**2
  !> beta = v/c; kappa = (|q|/(gamma*mass*v**3))||v X (E + v X B)||
  !> Power_tot = (q**2/(6*PI*eps0i*c**3))*gamma**4 * v**4 * kappa**2
  factor_1 = ((3.d0*TWOPI*(factor_1**2))/(sqrt3*(light_properties(10)**4)*&
  (light_properties(11)**8)*(light_properties(12)**3)))
  if(in_parallel) then
    !$omp taskloop default(shared) private(ii,jj) &
    !$omp firstprivate(one_over_gamma,rpsichi,z_cos,&
    !$omp factor_2,z_value,z2_value,factor_1) collapse(2)
    do ii=1,spectra%n_spectra
      do jj=1,spectra%n_points
        call compute_synchrotron_directionality_funct(light_vert%n_x,&
        spectra%points(jj,ii),light_properties(12), one_over_gamma,rpsichi,&
        z_cos,factor_2,z_value,z2_value,factor_1,light_dstb(jj,ii))
      enddo
    enddo
    !$omp end taskloop
  else
    !$omp parallel do default(shared) private(ii,jj) &
    !$omp firstprivate(one_over_gamma,rpsichi,z_cos,&
    !$omp factor_2,z_value,z2_value,factor_1) collapse(2)
    do ii=1,spectra%n_spectra
      do jj=1,spectra%n_points
        call compute_synchrotron_directionality_funct(light_vert%n_x,&
        spectra%points(jj,ii),light_properties(12),one_over_gamma,rpsichi,&
        z_cos,factor_2,z_value,z2_value,factor_1,light_dstb(jj,ii))
      enddo
    enddo
    !$omp end parallel do
  endif
end subroutine synchrotron_directionality_funct

!> synchrotron_spectral_irradiance computes the full spectral angular
!> power distribution for synchrotron lights which (L. Carbajal, PPCF, 2017)
!> emitted towards the shaded point x_shaded
!> inputs:
!>   light_vert: (synchrotron_light vertices) synchrotron light sources
!>   spectra:    (spectrum_base) spectral intervals and integrators
!>   time_id:    (integer) the time index
!>   light_id:   (integer) the light index
!>   x_shaded:   (real8)(3) shaded point position in cartesian coord
!> outputs: 
!>   light_vert: (synchrotron_light vertices) synchrotron light sources
!>   spectra:    (spectrum_base) spectral intervals and integrators
!>   light_spec_irradiance:  (real8)(n_points,n_spectra) synchrotron full spectral
!>                           angular distribution per unit of total power at the
!>                           at the shaded point x_shaded
subroutine synchrotron_spectral_irradiance(light_vert,spectra,time_id,&
light_id,x_shaded,light_spec_irradiance)
  use mod_spectra,  only: spectrum_base
  implicit none
  !> inputs-outputs:
  class(synchrotron_light_vertices),intent(inout) :: light_vert
  class(spectrum_base),intent(inout)              :: spectra
  !> inputs:
  integer,intent(in)                              :: time_id,light_id
  real*8,dimension(light_vert%n_x),intent(in)     :: x_shaded
  !> outputs:
  real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_spec_irradiance

  !> compute the directionality function
  call light_vert%directionality_funct(spectra,time_id,light_id,x_shaded,light_spec_irradiance)
  !> multiply the directionality function by the total synchrotron power
  light_spec_irradiance = light_spec_irradiance*light_vert%properties(13,light_id,time_id)
end subroutine synchrotron_spectral_irradiance

!> interpolate the JOREK MHD fields required for computing the
!> synchrotron radiation properties
!> inputs:
!>   light_vert: (synchrotron_light_vertices) empty synchrotron lights
!>   fields:     (fields_base) JOREK MHD fields
!>   particle:   (particle_base) JOREK particle base structure
!>   mass:       (real8) particle mass
!> outputs:
!>   mhd_fields: (real8)(n_mhd) JOREK MHD fields in cartesian coordinates
!>               1-3: x,y,z electric field componenets
!>               4-6: x,y,z magnetic field componenets
subroutine compute_synchrotron_mhd_fields(light_vert,fields,&
particle,mass,mhd_fields)
  use mod_fields,                only: fields_base
  use mod_particle_types,        only: particle_base
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  !> used only for unit testing but required for compilation
#ifdef UNIT_TESTS_AFIELDS
  use mod_particle_common_test_tools, only: compute_test_E_B_fields
#endif
  implicit none
  !> Inputs:
  class(synchrotron_light_vertices),intent(in) :: light_vert
  class(fields_nase),intent(in)                :: fields
  class(particle_base),intent(in)              :: particle_base
  real*8,intent(in)                            :: mass
  !> Outputs:
  real*8,dimension(light_vert%n_mhd),intent(out) :: mhd_fields
  !> Variables:
  real*8 :: psi,U
  !> compute the MHD fields
#ifndef UNIT_TESTS_AFIELDS
  !> compute JOREK electric and magnetic JOREK fields in cartesian coordinates
  call fields%calc_EBpsiU(light_vert%times(ii),particle%i_elm,&
  particle%st,particle%x(3),mhd_fields(1:3),mhd_fields(4:6),psi,U)
  mhd_fields(1:3) = vector_cylindrical_to_cartesian(particle%x(3),mhd_fields(1:3))
  mhd_fields(4:6) = vector_cylindrical_to_cartesian(particle%x(3),mhd_fields(4:6))
#else
  !> analytical fields only for unit testing
  call compute_test_E_B_fields(particle%x,mhd_fields(1:3),mhd_fields(4:6))
#endif
end subroutine compute_synchrotron_mhd_fields

!> compute_synchrotron_light_properties computes the
!> synchrotron radiation properties from a
!> kinetic relativistic particle. Variables
!> inputs:
!>   light_vert:  (synchrotron_light_vertices) empty synchrotron lights
!>   property_id: (integer) index of the property to be initialised
!>   time_id:     (integer) time index
!>   particle_in: (particle_kinetic_relativistic) jorek particle
!>   mass:        (real8) mass of the particle
!>   mhd_fields:  (real8)(n_mhd) JOREK MHD fields in cartesian coordinates
!>                1-3: x,y,z electric field componenets
!>                4-6: x,y,z magnetic field componenets
!> outputs:
!>   light_vert: (synchrotron_light_vertices) synchrotron lights with
!>               initialised properties. First dimension of the properties are:
!>                 1:3 -> component of the velocity direction (cartesian)
!>                     -> T = v/||v||
!>                 4:6 -> second orthonormal basis cartesian coordinates
!>                     -> N = E + v X B - v*E
!>                 7:9 -> components of the third orthonormal basis (cartesian)
!>                     -> B = T X N
!>                 10  -> beta -> velocity/speed of light = v/c
!>                 11  -> relativistic factor gamma = sqrt(1+(p/(mass*c))**2)
!>                 12  -> orbit curvature (L. Carbakal, PPCF, 2017)
!>                     -> kappa = (|q|/(gamma*mass*v**3))||v X (E + v X B)||
!>                 13  -> total radiation power (L. Carbajal, PPCF, 2017)i
!>                     -> P_tot = (q**2/(6*PI*eps0*c**3))*gamma**4 * v**4 * kappa**2
subroutine compute_synchrotron_light_properties(light_vert,&
property_id,time_id,particle_in,mass,mhd_fields)
  use constants,                 only: PI,EPS_ZERO,EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
  use mod_math_operators,        only: cross_product
  use mod_coordinate_transforms, only: vectors_to_orthonormal_basis
  use mod_particle_types,        only: particle_base,particle_kinetic_relativistic
  implicit none
  !> inputs-outputs
  class(synchrotron_light_vertices),intent(inout) :: light_vert
  !> inputs
  class(particle_base),intent(in)                 :: particle_in
  integer,intent(in)                              :: property_id,time_id
  real*8,intent(in)                               :: mass
  real*8,dimension(light_vert%n_mhd),intent(in)   :: mhd_fields
  !> variables
  real*8 :: velocity
  real*8,dimension(field_size) :: vector_1d_3,vector_1d_3_2,vector_1d_3_3

  select type(p_in=>particle_in)
    type is (particle_kinetic_relativistic)
    !> compute velocity, velocity direction and relativistic factor
    velocity = sqrt(p_in%p(1)**2+p_in%p(2)**2+p_in%p(3)**2) 
    light_vert%properties(1:3,property_id,time_id) = p_in%p/velocity
    light_vert%properties(10,property_id,time_id)  = velocity/SPEED_OF_LIGHT
    light_vert%properties(11,property_id,time_id)  = sqrt(1.d0 + &
                           (light_vert%properties(10,property_id,time_id)**2)/(mass**2))
    light_vert%properties(10,property_id,time_id)  = light_vert%properties(10,property_id,time_id)/&
                           (mass*light_vert%properties(11,property_id,time_id))
  !> compute orbit curvature
    light_vert%properties(4:6,property_id,time_id) = mhd_fields(1:3)+cross_product(p_in%p/&
                           (mass*light_vert%properties(11,property_id,time_id)),mhd_fields(4:6))
    vector_1d_3 = cross_product(light_vert%properties(1:3,property_id,time_id),&
                           light_vert%properties(4:6,property_id,time_id))
    light_vert%properties(12,property_id,time_id)  = (abs(real(p_in%q,kind=8))*EL_CHG*&
                           sqrt(vector_1d_3(1)**2+vector_1d_3(2)**2+vector_1d_3(3)**2))/&
                           (light_vert%properties(11,property_id,time_id)*mass*ATOMIC_MASS_UNIT*&
                           ((light_vert%properties(10,property_id,time_id)*SPEED_OF_LIGHT)**2))
    !> compute total synchrotron power
    light_vert%properties(13,property_id,time_id)  = (((EL_CHG*real(p_in%q,kind=8))**2)*&
                           SPEED_OF_LIGHT*(light_vert%properties(10,property_id,time_id)**4)*&
                           (light_vert%properties(11,property_id,time_id)**4)*&
                           (light_vert%properties(12,property_id,time_id)**2))/(6.d0*PI*EPS_ZERO)
  end select
  !> construct and store the orthonormal basis
  call vectors_to_orthonormal_basis(light_vert%properties(1:3,property_id,time_id),&
  light_vert%properties(4:6,property_id,time_id),vector_1d_3,vector_1d_3_2,vector_1d_3_3)
  light_vert%properties(1:3,property_id,time_id) = vector_1d_3; 
  light_vert%properties(4:6,property_id,time_id) = vector_1d_3_2;
  light_vert%properties(7:9,property_id,time_id) = vector_1d_3_3
end subroutine compute_synchrotron_light_properties

!> initialise and allocate synchrotron light variables
!> inputs:
!>   light_vert: (synchrotron_light_vertices) synchrotron lights class
!> outputs:
!>   light_vert: (synchrotron_light_vertices) synchrotron lights class
subroutine setup_synchrotron_light_class(light_vert)
  use mod_particle_types, only: particle_kinetic_relativistic_id
  implicit none
  !> inputs-outputs
  class(synchrotron_light_vertices),intent(inout) :: light_vert
  !> set-up the synchrotron light variables 
  light_vert%n_property_vertex = 13; light_vert%n_mhd = 6;
  light_vert%n_particle_types = 1;
  if(allocate(light_vert%particle_types)) deallocate(light_vert%particle_types)
  allocate(light_vert%particle_types(light_vert%n_particle_types))
  light_vert%particle_types = [particle_kinetic_relativistic_id]
end subroutine setup_synchrotron_light_class

!> Tools ------------------------------------------
!> compute the synchrotron radiation directionality function
subroutine compute_synchrotron_directionality_funct(n_x,&
wavelength,orbit_curvature,one_over_gamma,rpsichi,z_cos,&
factor_2,z_value,z2_value,factor_1,dir_funct)
  use constants,         only: TWOPI
  use mod_besselk, only: f_besselk
  implicit none
  !> inputs:
  integer,intent(in)               :: n_x
  real*8,intent(in)                :: one_over_gamma,z_value,z2_value
  real*8,intent(in)                :: wavelength,orbit_curvature,factor_1
  real*8,intent(in)                :: factor_2,z_cos
  real*8,dimension(n_x),intent(in) :: rpsichi
  !> outputs:
  real*8,intent(out) :: dir_funct
  !> variables:
  real*8 :: zeta,besselk_1,besselk_2
  !> compute the directionality function
  dir_funct = 0.d0;  
  !> zeta = (2*PI*(1/gamma**2 + psi**2)**(3/2))/(3*kappa*lambda)
  zeta = TWOPI*((one_over_gamma+(rpsichi(2)**2))**1.5d0)/&
  (3.d0*wavelength*orbit_curvature)
  besselk_1 = f_besselk(onethird,zeta); besselk_2 = f_besselk(twothirds,zeta);
  if(isnan(besselk_1)) return; if(isnan(besselk_2)) return;
  !> funct = I*(K_1/3(zeta)*cos(zeta*z_cos)*(((gamma**2 * psi**2)/(1 + gamma**2 * psi**2)) -
  !>  0.5*(1+z**2)) + K_(2/3)(zeta)*sin(zeta*z_cos))/lambda**4
  dir_funct = factor_1*(besselk_1*cos(zeta*z_cos)*(factor_2-z2_value)+&
  besselk_2*z_value*sin(zeta*z_cos))/(wavelength**4)
end subroutine compute_synchrotron_directionality_funct

!> check if the shaded point of the synchrotron light is inside
!> the synchrotron radiation cone of half width sin(theta) ≃ 1/gamma
!> where gamma is the relativistic factor. 
!> inputs:
!>   n_x:          (integer) size of the coordinate system
!>   x_shaded:     (real8)(n_x) position of the shaded point
!>   x_light:      (real8)(n_x) position of the point light
!>   n_int_param:  (integer) number of integer parameters: 0
!>   n_real_param: (real8) number of real parameters: 4
!>   int_param:    (integer)(n_int_param) integer parameters
!>   real_param:   (real8)(n_real_param) real_parameters:
!>                 1- x component of the emission cone (momentum) direction
!>                 2- y component of the emission cone (momentum) direction
!>                 3- z component of the emission cone (momentum) direction
!>                 4- relativistic factor
!> outouts:
!>   in_code: (logical) if true the gather point is in the synchrotron cone
function check_shaded_x_in_synchrotron_cone(n_x,x_shaded,x_light,&
n_int_param,n_real_param,int_param,real_param) result(in_range)
  !> Inputs:
  integer,intent(in)                        :: n_x,n_int_param,n_real_param
  integer,dimension(n_int_param),intent(in) :: int_param
  real*8,dimension(n_x),intent(in)          :: x_shaded,x_light
  real*8,dimension(n_real_param),intent(in) :: real_param
  !> Outputs:
  logical :: in_range
  !> Variables:
  real*8 :: costheta
  !> initialisation
  in_range = .false.
  !> check if the shaded point is in the synchrotron conede
  costheta = dot_product(x_shaded-x_light,real_param(1:3))/norm2(x_shaded-x_light)
  if((costheta.ge.0).and.((sqrt(1d0-costheta**2)*real_param(4)).le.1d0)) in_range=.true.
end function check_shaded_x_in_synchrotron_cone

!>-------------------------------------------------
end module mod_synchrotron_light_vertices
