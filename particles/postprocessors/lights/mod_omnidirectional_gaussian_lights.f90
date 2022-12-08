!> the mod_omnidirectional_gaussian_lights implements variables
!> and procedure defining an omnidirectional light with a
!> gaussian spectrum
module mod_omnidirectional_gaussian_lights
use mod_light_vertices, only: light_vertices
implicit none

private
public :: omnidirectional_gaussian_lights

!> Variables ----------------------------------------------------
type,extends(light_vertices) :: omnidirectional_gaussian_lights
  real*8 :: light_intensity=1.d0
  contains
  procedure,pass(light_vert) :: init_lights_from_particles => &
                                init_omnidir_gaussian_lights_from_particles
  procedure,pass(light_vert) :: directionality_funct => &
                                omnidir_gaussian_directionality_funct
  procedure,pass(light_vert) :: spectral_irradiance => &
                                omnidir_gaussian_spectral_irradiance
  procedure,pass(light_vert) :: compute_mhd_fields => &
                                compute_omnidirectional_mhd_fields
  procedure,pass(light_vert) :: compute_light_properties => & 
                                compute_omnidirectional_light_properties
  procedure,pass(light_vert) :: setup_light_class => &
                                setup_omnidirectional_light_class
  procedure,nopass           :: check_x_shaded_in_emission_zone => &
                                check_shaded_x_omnidirectional_light
end type omnidirectional_gaussian_lights
!> Interfaces ---------------------------------------------------

contains

!> Procedures ---------------------------------------------------
!> interpolate the JOREK MHD fields required for computing the
!> omnidirectional radiation properties
!> inputs:
!>   light_vert: (synchrotron_light_vertices) empty synchrotron lights
!>   fields:     (fields_base) JOREK MHD fields
!>   particle:   (particle_base) JOREK particle base structure
!>   mass:       (real8) particle mass
!> outputs:
!>   mhd_fields: (real8)(n_mhd) JOREK MHD fields in cartesian coordinates
!>               1 -> intensity of the magnetic field
subroutine compute_omnidirectional_mhd_fields(light_vert,fields,&
particle_in,mass,mhd_fields)
  use mod_fields,                only: fields_base
  use mod_particle_types,        only: particle_base
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  !> used only for unit testing but required for compilation
#ifdef UNIT_TESTS_AFIELDS
  use mod_particle_common_test_tools, only: compute_test_E_B_fields
#endif
  implicit none
  !> Inputs:
  class(omnidirectional_gaussian_lights),intent(in) :: light_vert
  class(fields_nase),intent(in)                     :: fields
  class(particle_base),intent(in)                   :: particle_in
  real*8,intent(in)                                 :: mass
  !> Outputs:
  real*8,dimension(light_vert%n_mhd),intent(out)    :: mhd_fields
  !> Variables:
  real*8 :: psi,U
  real*8,dimension(3) :: E,B
  !> compute the MHD fields
#ifndef UNIT_TESTS_AFIELDS
  !> compute JOREK electric and magnetic JOREK fields in cartesian coordinates
  call fields%calc_EBpsiU(light_vert%times(ii),particle_in%i_elm,&
  particle_in%st,particle_in%x(3),E,B,psi,U)
#else
  !> analytical fields only for unit testing
  call compute_test_E_B_fields(particle_in%x,E,B)
#endif
  mhd_fields(1) = norm2(B)
end subroutine compute_omnidirectional_mhd_fields

!> compute omnidirectional light properties from kinetic relativistic
!> and guiding center particles.
!> inputs:
!>   light_vert:  (omnidirectional_light_vertices) empty synchrotron lights
!>   property_id: (integer) index of the property to be initialised
!>   time_id:     (integer) time index
!>   particle_in: (particle_kinetic_relativistic) jorek particle
!>   mass:        (real8) mass of the particle
!>   mhd_fields:  (real8)(n_mhd) JOREK MHD fields in cartesian coordinates
!>                  1 -> norm of the magnetic field
!> outputs:
!>   light_vert: (omnidirectional_light_vertices) synchrotron lights with
!>               initialised properties. First dimension of the properties are:
!>                 1 -> relativistic factor gamma = sqrt(1+(p/(mass*c))**2)
!>                 2 -> total power 1/sqrt(2*pi*relativistic_factor)
subroutine compute_omnidirectional_light_properties(light_vert,property_id,&
time_id,particle_in,mass,mhd_fields)
  use constants,                 only: TWOPI,SPEED_OF_LIGHT
  use mod_math_operators,        only: cross_product
  use mod_coordinate_transforms, only: vectors_to_orthonormal_basis
  use mod_particle_types,        only: particle_base,particle_kinetic_relativistic
  use mod_particle_types,        only, particle_gc_relativistic
  implicit none
  !> inputs-outputs
  class(omnidirectional_gaussian_lights),intent(inout) :: light_vert
  !> inputs
  class(particle_base),intent(in)                 :: particle_in
  integer,intent(in)                              :: property_id,time_id
  real*8,intent(in)                               :: mass
  real*8,dimension(light_vert%n_mhd),intent(in)   :: mhd_fields
  select type (p_in => particle_in)
    type is (particle_kinetic_relativistic)
    light_vert%properties(1,property_id,time_id) = sqrt(1d0 + &
    (dot_product(p_in%p,p_in%p)/((mass*SPEED_OF_LIGHT)**2)))
    type is (particle_gc_relativistic)
    light_vert%properties(1,property_id,time_id) = sqrt(1d0 + &
    ((p_in%p(q)*p_in%p(1))/((mass*SPEED_OF_LIGHT)**2)) + &
    ((2d0*p_in%p(2)*mhd_fields(1))/(mass*(SPEED_OF_LIGHT**2))))
  end select
  light_vert%properties(2,property_id,time_id) = 1d0/sqrt(TWOPI*&
  light_vert%properties(1,property_id,time_id))
end subroutine compute_omnidirectional_light_properties 

!> check if the shaded point is within the emission range.
!> Given that it is an omnidirectional light, the function returns true
!> inputs:
!>   n_x:          (integer) size of the coordinate system
!>   x_shaded:     (real8)(n_x) position of the shaded point
!>   x_light:      (real8)(n_x) position of the point light
!>   n_int_param:  (integer) number of integer parameters: 0
!>   n_real_param: (real8) number of real parameters: 0
!>   int_param:    (integer)(n_int_param) integer parameters
!>   real_param:   (real8)(n_real_param) real_parameters
!> outouts:
!>   in_range: (logical) always true
function check_shaded_x_omnidirectional_light(n_x,x_shaded,x_light,&
n_int_param,n_real_param,int_param,real_param) result(in_range)
  implicit none
  !> Inputs:
  integer,intent(in)                        :: n_x,n_int_param,n_real_param
  integer,dimension(n_int_param),intent(in) :: int_param
  real*8,dimension(n_x),intent(in)          :: x_shaded,x_light
  real*8,dimension(n_real_param),intent(in) :: real_param
  !> Outputs:
  logical :: in_range
  in_range = .true.
end function check_shaded_x_omnidirectional_light

!> omnidir_gaussian_spectral_irradiance computes the full spectral anguler
!> power distribution for omnidirectional gaussian lights
!> inputs:
!>   light_vert: (omnidirectional_gaussian_lights) omnidirectional gaussian lights
!>   spectra:    (spectrum_base) spectral intervals and integrators
!>   time_id:    (integer) the time index
!>   light_id:   (integer) the light index
!>   x_shaded:   (real8)(3) shaded point position in cartesian coord
!> outputs:
!>   light_vert:            (omnidirectional_gaussian_lights) omnidirectional gaussian lights
!>   spectra:               (spectrum_base) spectral intervals and integrators  
!>   light_spec_irradiance: (real8)(n_points,n_spectra) omnidirectional gaussian
!>                          full spectral angular distribution
subroutine omnidir_gaussian_spectral_irradiance(light_vert,spectra,time_id,light_id,&
x_shaded,light_spec_irradiance)
  use mod_spectra, only: spectrum_base
  !$ use omp_lib
  implicit none
  !> inputs-outputs:
  class(omnidirectional_gaussian_lights),intent(inout) :: light_vert
  class(spectrum_base),intent(inout)                   :: spectra
  !> inputs:
  integer,intent(in)                                   :: time_id,light_id
  real*8,dimension(light_vert%n_x),intent(in)          :: x_shaded
  !> outputs:
  real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_spec_irradiance
  !> variables
  logical :: in_parallel
  integer :: ii,jj
  real*8,dimension(spectra%n_spectra) :: spectra_midpoint,spectra_bin
  !> initialisation
  in_parallel = .false.
  !$ in_parallel = omp_in_parallel()
  spectra_bin = (spectra%points(spectra%n_points,:)-spectra%points(1,:))
  spectra_midpoint = spectra%points(1,:)+5.d-1*spectra_bin
  !> compute irradiance
  if(in_parallel) then
    !$omp taskloop default(shared) private(ii,jj) collapse(2)
    do ii=1,spectra%n_spectra
      do jj=1,spectra%n_points
        light_spec_irradiance(jj,ii) = light_vert%light_intensity*&
        exp(-((spectra%points(jj,ii)-spectra_midpoint(ii))*&
        (spectra%points(jj,ii)-spectra_midpoint(ii)))/(2.d0*&
        spectra_bin(ii)*light_vert%properties(1,light_id,time_id)))
      enddo
    enddo
    !$omp end taskloop
  else
    !$omp parallel do default(shared) private(ii,jj) collapse(2)
    do ii=1,spectra%n_spectra
      do jj=1,spectra%n_points
        light_spec_irradiance(jj,ii) = light_vert%light_intensity*&
        exp(-((spectra%points(jj,ii)-spectra_midpoint(ii))*&
        (spectra%points(jj,ii)-spectra_midpoint(ii)))/(2.d0*&
        spectra_bin(ii)*light_vert%properties(1,light_id,time_id)))
      enddo
    enddo
    !$omp end parallel do
  endif
end subroutine omnidir_gaussian_spectral_irradiance

!> omnidir_gaussian_directionality_funct computes the directionality function
!> for omnidirectional gaussian lights
!> inputs:
!>   light_vert: (omnidirectional_gaussian_lights) omnidirectional gaussian lights
!>   spectra:    (spectrum_base) spectral intervals and integrators
!>   time_id:    (integer) the time index
!>   light_id:   (integer) the light index
!>   x_shaded:   (real8)(3) shaded point position in cartesian coord
!> outputs:
!>   light_vert: (omnidirectional_gaussian_lights) omnidirectional gaussian lights
!>   spectra:    (spectrum_base) spectral intervals and integrators  
!>   light_dstb: (real8)(n_points,n_spectra) omnidirectional gaussian light
!>               directionality function
subroutine omnidir_gaussian_directionality_funct(light_vert,spectra,time_id,light_id,&
x_shaded,light_dstb)
  use mod_spectra, only: spectrum_base
  implicit none
  !> inputs-outputs:
  class(omnidirectional_gaussian_lights),intent(inout) :: light_vert
  class(spectrum_base),intent(inout)                   :: spectra
  !> inputs:
  integer,intent(in)                                   :: time_id,light_id
  real*8,dimension(light_vert%n_x),intent(in)          :: x_shaded
  !> outputs:
  real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_dstb
  !> compute the directionality function
  call light_vert%spectral_irradiance(spectra,time_id,light_id,x_shaded,light_dstb)
  light_dstb = light_dstb/(light_vert%properties(2,light_id,time_id)*light_vert%light_intensity)
end subroutine omnidir_gaussian_directionality_funct

!> initialise and allocate synchrotron light variables
!> inputs:
!>   light_vert: (omnidirectional_gaussian_lights) omnidirectional lights class
!> outputs:
!>   light_vert: (omnidirectional_gaussian_lights) omnidirectional lights class
subroutine setup_omnidirectional_light_class(light_vert)
  use mod_particle_types, only: particle_kinetic_relativistic_id
  use mod_particle_types, only: particle_gc_relativistic_id
  implicit none
  !> inputs-outputs
  class(omnidirectional_gaussian_lights),intent(inout) :: light_vert
  !> set-up the omnidirectional light variables 
  light_vert%n_property_vertex = 2; light_vert%n_mhd = 1;
  light_vert%n_particle_types = 2;
  if(allocated(light_vert%particle_types)) deallocate(light_vert%particle_types)
  allocate(light_vert%particle_types(light_vert%n_particle_types))
  light_vert%particle_types = [particle_kinetic_relativistic_id,&
                               particle_gc_relativistic_id]
end subroutine setup_omnidirectional_light_class

!> Tools --------------------------------------------------------

!>---------------------------------------------------------------
end module mod_omnidirectional_gaussian_lights

