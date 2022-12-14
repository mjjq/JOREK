!> the mod_gyroaverage_synchrotron_dist_light implements
!> variables and procedures for defining the gyroaverage
!> synchrotron light distribution of guiding center light
!> sources. The model used is:
!> M. Hoppe et al., Nucl. Fusion, vol.58, p.026032, 2018
module mod_gyroaverage_synchrotron_dist_light_vertices
use mod_synchrotron_light_vertices, only: synchrotron_light
implicit none

private
public :: gyroaverage_synchrotron_light_dist

!> Variables ---------------------------------------------
real*8,parameter :: onethirds=1d0/3d0
real*8,parameter :: twothirds=2d0/3d0
type,extends(synchrotron_light) :: gyroaverage_synchrotron_light_dist
  contains
  procedure,pass(light_vert) :: directionality_funct => &
                                gyroaverage_synchrotron_directionality_funct
  procedure,pass(light_vert) :: spectral_irradiance => &
                                gyroaverage_synchrotron_spectral_irradiance
  procedure,pass(light_vert) :: compute_mhd_fields => &
                                compute_gyroaverage_synchrotron_mhd_fields
  procedure,pass(light_vert) :: compute_light_properties => & 
                                compute_gyroaverage_synchrotron_light_properties
  procedure,pass(light_vert) :: setup_light_class => &
                                setup_gyroaverage_synchrotron_light_class
end type gyroaverage_synchrotron_light_dist
!> Interfaces --------------------------------------------

contains

!> Procedures --------------------------------------------

!> directionality function for gyroaverage synchrotron radiation emission
!> as explained in M. Hoppe, Nucl. Fuson, vol.58, p.026032, 2018
!> inputs:
!>   light_vert:  (gyroaverage_synchrotron_light_dist) 
!>                gyroaverage synchrotron light class
!>   spectra:     (spectrum_base) spectral intervals and integrators
!>   time_id:     (integer) time index
!>   light_id:    (integer) light index
!>   x_shaded:    (real8)(n_x) shaded point position in cartesian coords
!> outputs:
!>   light_vert:  (gyroaverage_synchrotron_light_dist) 
!>                gyroaverage synchrotron light class
!>   spectra:     (spectrum_base) spectral intervals and integrators
!>   light_dstb:  (real8)(n_points,n_intervals) gyroaverage synchrotron
!>                radiation spectrum per unit of total power emitted
!>                towards the shaded point x_shaded
subroutine gyroaverage_synchrotron_directionality_funct(light_vert,spectra,&
time_id,light_id,x_shaded,light_dstb)
  use mod_spectra, only: spectrum_base
  !$ use omp_lib
  implicit none
  !> Inputs-Outputs:
  class(gyroaverage_synchrotron_light_dist),intent(in) :: light_vert
  class(spectrum_base),intent(inout)                   :: spectra
  !> Inputs:
  integer,intent(in)                          :: time_id,light_id
  real*8,dimension(light_vert%n_x),intent(in) :: x_shaded
  !> Outputs:
  real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_dstb
  !> Variables:
  logical :: in_parallel
  integer :: ii,jj
  integer,dimension(0) :: int_param
  real*8 :: cosmu,cospsi,fact_1,fact_2,fact_3
  real*8,dimension(light_vert%n_property_vertex) :: light_properties

  !> initialisations
  in_parallel = .false.; light_dstb = 0d0;
  !$ in_parallel = omp_in_parallel()
  light_properties = light_vert%properties(:,light_id,time_id)
  !> check if the shaded point is in the synchrotron emission cone
  if(.not.light_vert%check_x_shaded_in_emission_zone(light_vert%n_x,&
  x_shaded,light_vert%x(:,light_id,time_id),0,4,int_param,[light_properties(1),&
  light_properties(2),light_properties(3),light_properties(4)])) return
  !> the mu and psi angle cosinues, note that psi = mu-thetap hence
  !> cos(mu-thetap) = cos(mu)*cos(thetap) + sin(mu)*sin(thetap)
  cosmu = (x_shaded(1)-light_vert%x(1,light_id,time_id))*light_properties(1)+&
          (x_shaded(2)-light_vert%x(2,light_id,time_id))*light_properties(2)+&
          (x_shaded(3)-light_vert%x(3,light_id,time_id))*light_properties(3)
  cosmu = cosmu/sqrt(((x_shaded(1)-light_vert%x(1,light_id,time_id))**2) +&
          ((x_shaded(2)-light_vert%x(2,light_id,time_id))**2) +&
          ((x_shaded(3)-light_vert%x(3,light_id,time_id))**2))
  cospsi = cosmu*light_properties(6) + light_properties(7)*sqrt(1d0-cosmu*cosmu)
  !> compute the directionality function factors
  fact_3 = 5d-1*cospsi*(1d0-cospsi*cospsi)
  cospsi = light_properties(5)*cospsi !< becareful fron now on cospsi = beta*cospsi
  fact3 = fact_3/(1d0-cospsi)
  fact_2 = (1d0-light_properties(5)*light_properties(6)*cosmu)*(((1d0-cospsi)/cospsi)**2)
  fact_1 = (sqrt(((1d0-cospsi)**3)/(5d-1*cospsi)))*(light_properties(4)**3)
  !> compute the directionality function
  if(in_parallel) then
    !$omp taskloop default(shared) private(ii,jj) &
    !$omp firstprivate(fact_1,fact_2,fact_3) collapse(2)
    do ii=1,spectra%n_spectra
      do jj=1,spectra%n_points
        call compute_synchrotron_directionality_funct(spectra%points(jj,ii),&
        light_properties(8),fact_1,fact_2,fact_3,light_properties(9),light_dstb(jj,ii))
      enddo
    enddo
    !$omp end taskloop
  else
    !$omp parallel do default(shared) private(ii,jj) &
    !$omp firstprivate(fact_1,fact_2,fact_3) collapse(2)
    do ii=1,spectra%n_spectra
      do jj=1,spectra%n_points
        call compute_synchrotron_directionality_funct(spectra%points(jj,ii),&
        light_properties(8),fact_1,fact_2,fact_3,light_properties(9),light_dstb(jj,ii))
      enddo
    enddo
    !$omp end parallel do
  endif
end subroutine gyroaverage_synchrotron_directionality_funct

!> computes the spectral irradiance of the gyroaveraged synchrotron radiation
!> using the model of: M. Hoppe, Nucl. Fuson, vol.58, p.026032, 2018
!> emitted towards the shaded point x_shaded
!> inputs:
!>   light_vert:  (gyroaverage_synchrotron_light_dist) 
!>                gyroaverage synchrotron light class
!>   spectra:     (spectrum_base) spectral intervals and integrators
!>   time_id:     (integer) time index
!>   light_id:    (integer) light index
!>   x_shaded:    (real8)(n_x) shaded point position in cartesian coords
!> outputs:
!>   light_vert:           (gyroaverage_synchrotron_light_dist) 
!>                          gyroaverage synchrotron light class
!>   spectra:               (spectrum_base) spectral intervals and integrators
!>   light_spec_irradiance: (real8)(n_points,n_spectra) spectral irradiance
!>                          of the gyroaveraged synchrotron radiation
!>                          spectral distribution
subroutine gyroaverage_synchrotron_spectral_irradiance(light_vert,spectra,&
time_id,light_id,x_shaded,light_spec_irradiance)
  use mod_spectra, only: spectrum_base
  implicit none
  !> Inputs-Outputs:
  class(gyroaverage_synchrotron_light_dist),intent(in) :: light_vert
  class(spectrum_base),intent(inout)                   :: spectra
  !> Inputs:
  integer,intent(in)                          :: time_id,light_id
  real*8,dimension(light_vert%n_x),intent(in) :: x_shaded
  !> Outputs:
  real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_spec_irradiance

  !> compute the directionality function
  call light_vert%directionality_funct(spectra,time_id,light_id,x_shaded,light_spec_irradiance)
  !> denormalise the directionality function
  light_spec_irradiance = light_spec_irradiance*light_vert%properties(10,light_id,time_id)
end subroutine gyroaverage_synchrotron_spectral_irradiance

!> interpolate the JOREK MHD fields required for computing
!> the synchrotron radiation properties
!> inputs:
!>   light_vert:  (gyroaverage_synchrotron_light_dist) 
!>                gyroaverage synchrotron light class
!>   fields:      (fields_base) JOREK MHD fields
!>   particle_in: (particle_base) JOREK particle class
!>   time_id:     (integer) particle simulation time index
!>   mass:        (real8) particle mass
!>   outputs:
!>     mhd_fields: (real8)(n_mhd) JOREK MHD fields in cartesian coordinates
!>                 1-3:   x,y,z components of the electric field
!>                 4-6:   x,y,z components of the magnetic field direction: b
!>                 7:     intensity of the magnetic field: ||B||
!>                 8-10:  x,y,z components of the gradient of ||B||
!>                 11-13: x,y,z components of the curl of b
!>                 14-16: x,y,z components of the b time derivative
subroutine compute_synchrotron_mhd_fields(light_vert,fields,&
particle_in,time_id,mass,mhd_fields)
  use mod_fields,                only: fields_base
  use mod_particle_types,        only: particle_base
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
#ifdef UNIT_TESTS_AFIELDS
  !> used only for unit testing
  use mod_particle_common_test_tools, only: compute_test_E_B_normB_gradB_curlb_Dbdt_fields
#endif
  implicit none
  !> Inputs:
  class(gyroaverage_synchrotron_light_dist),intent(in) :: light_vert
  class(fields_base),intent(in)                        :: fields
  class(particle_base),intent(in)                      :: particle_in
  integer,intent(in)                                   :: time_id
  real*8,intent(in)                                    :: mass
  !> Outputs:
  real*8,dimension(light_vert%n_mhd),intent(out) :: mhd_fields
#ifndef UNIT_TESTS_AFIELDS
  !> compute the JOREK magnetic field
  call fields%calc_EBNormBGradBCurlbDbdt(light_vert%times(time_id),&
  particle_in%i_elm,particle_in%st,particle_in%x(3),mhd_fields(1:3),&
  mhd_fields(4:6),mhd_fields(7),mhd_fields(8:10),mhd_fields(11:13),&
  mhd_fields(14:16))
#else
  !> compute the analytical magnetic field
  call compute_test_E_B_normB_gradB_curlb_Dbdt_fields(particle_in%x,&
  mhd_fields(1:3),mhd_fields(4:6),mhd_fields(7),mhd_fields(8:10),&
  mhd_fields(11:13),mhd_fields(14:16))
#endif
  mhd_fields(1:3)   = vector_cylindrical_to_cartesian(particle_in%x(3),mhd_fields(1:3))
  mhd_fields(4:6)   = vector_cylindrical_to_cartesian(particle_in%x(3),mhd_fields(4:6))
  mhd_fields(8:10)  = vector_cylindrical_to_cartesian(particle_in%x(3),mhd_fields(8:10))
  mhd_fields(11:13) = vector_cylindrical_to_cartesian(particle_in%x(3),mhd_fields(11:13))
  mhd_fields(14:16) = vector_cylindrical_to_cartesian(particle_in%x(3),mhd_fields(14:16))
end subroutine compute_synchrotron_mhd_fields

!> compute the gyroaverage synchrotron light properties from a
!> relativistic guiding center particle.
!> inputs:
!>   light_vert:  (gyroaverage_synchrotron_light_dist) gyroaverage synchrotron lights
!>                with empty synchrotron properties
!>   property_id: (integer) index of the property to be initialised
!>   time_id:     (integer) time index
!>   particle_in: (particle_gc_relativistic) jorek relativistic gc
!>   mass:        (real8) particle mass
!>   mhd_fields:  (real8)(n_mhd) JOREK MHD fields for relativistic guiding centers
!> outputs:
!>   light_vert:  (gyroaverage_synchrotron_light_dist) gyroaverage synchrotron lights
!>                with initialised synchrotron properties:
!>                1:3 -> guiding center velocity direction
!>                4   -> relativitic factor
!>                5   -> beta = velocity/speed of light = v/c
!>                6   -> cos(thetap) cosinus of the pitch angle
!>                7   -> sin(thetap) sinus of the pitch angle
!>                8   -> critical wavelength lambda_c=(4*pi*c*gamma_parallel)/(3*q*B*gamma**2)
!>                9   -> directionality function intensity = 
!>                       (27*q*B*gamma**7)/(128*(pi**2)*m*c*(sin(thetap)**2)*gamma_parallel**4)
!>                10  -> synchrotron power normalisation:
!>                       (((q**2)*B*gamma*gamma_parallel*v_perp)**2)/(6*pi*epsilon0*(m**2)*(c**3))
subroutine compute_gyroaverage_synchrotron_light_properties(light_vert,&
property_id,time_id,particle_in,mass,mhd_fields)
  use constants,           only: PI,EPS_ZERO,EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT     
  use mod_particle_types,  only: particle_base,particle_gc_relativistic
  use mod_gc_relativistic, only: compute_relativistic_factor
  use mod_gc_relativistic, only: compute_relativistic_gc_rhs
  implicit none
  !> inputs-outputs
  class(gyroaverage_synchrotron_light_dist),intent(inout) :: light_vert
  !> inputs
  class(particle_base),intent(in)               :: particle_in
  integer,intent(in)                            :: property_id,time_id
  real*8,intent(in)                             :: mass
  real*8,dimension(light_vert%n_mhd),intent(in) :: mhd_fields
  !> variables
  real*8 :: charge_r8
  real*8 :: beta_perp2,rel_fact_parallel !< (perpendicular velocity/c)**2, parallel relativistic factor

  select type(p_in=>particle_in)
    type is (particle_gc_relativistic)
    !> compute the gc velocity direction
    light_vert%properties(1:4) = compute_relativistic_gc_rhs(integer(p_in%q,kind=4),&
                                 mass,p_in%p(2),p_in%x(1),p_in%p(1),mhd_fields(7),&
                                 mhd_fields(1:3),mhd_fields(4:6),mhd_fields(8:10),&
                                 mhd_fields(11:13),mhd_fields(14:16))
    light_vert%properties(1:3) = light_vert%properties(1:3)/sqrt(light_vert%properties(1)**2 + &
    light_vert%properties(2)**2 + light_vert%properties(3)**2)
    !> compute the relativistic factor
    light_vert%properties(4) = compute_relativistic_factor(p_in,mass,mhd_fields(7))
    !> compute beta beta = v/c
    light_vert%properties(5) = sqrt(1d0 - (1d0/(light_vert%properties(4)**2)))
    !> compute cosinus and sinus of the pitch angle beta_perp2 = (p_perp/(mass*gamma*c))**2
    !> (p_perp/(mass*c))**2 = (2*mu*B)/(mass*(c**2))
    charge_r8 = real(p_in%q,kind=8)
    beta_perp2 = (2d0*p_in%p(2)*mhd_fields(7))/(mass*((SPEED_OF_LIGHT*light_vert%properties(4))**2)) 
    light_vert%properties(6:7) = (/p_in%p(1)/(mass*SPEED_OF_LIGHT*light_vert%properties(4)),&
    sqrt(beta_perp2)/)/light_vert%properties(5)
    !> compute the critical wavelenght
    rel_fact_parallel = sqrt(1d0+((p_in%p(1)/(mass*SPEED_OF_LIGHT))**2))
    light_vert%properties(8) = (4d0*PI*mass*ATOMIC_MASS_UNIT*SPEED_OF_LIGHT*rel_fact_parallel)/&
                             (3d0*charge_r8*EL_CHG*mhd_fields(7)*(light_vert%properties(4)**2))
    !> compute the directionality function intensity
    light_vert%properties(9) = (2.7d1*EL_CHG*charge*mhd_fields(7)*(light_vert%properties(4)**7))/&
    (1.28d2*mass*ATOMIC_MASS_UNIT*SPEED_OF_LIGHT*((PI*light_vert%properties(8)*(rel_fact_parallel**2))**2))
    !> compute the synchrotron power for normalisation
    light_vert%properties(10) = (beta_perp2*((mhd_fields(7)*rel_fact_parallel*&
                                light_vert%properties(4)*((EL_CHG*charge)**2))**2))/&
                                (6d0*PI*EPS_ZERO*SPEED_OF_LIGHT*((mass*ATOMIC_MASS_UNIT)**2))
  end select
end subroutine compute_gyroaverage_synchrotron_light_properties

!> initialise and allocate synchrotron light variables
!> inputs:
!>   light_vert:  (gyroaverage_synchrotron_light_dist) 
!>                gyroaverage synchrotron light class
!>                to be initialised
!> outputs:
!>   light_vert:  (gyroaverage_synchrotron_light_dist) initialised 
!>                gyroaverage synchrotron light class
subroutine setup_gyroaverage_synchrotron_light_class
  use mod_particle_types, only: particle_gc_relativistic_id
  implicit none
  !> Inputs-Outputs:
  class(gyroaverage_synchrotron_light_dist),intent(in) :: light_vert
  !> set-up the gyroaverage synchrotron light variables
  light_vert%n_property_vertex = 10; light_vert%n_mhd = 16;
  light_vert%n_particle_types = 1;
  if(allocated(light_vert%particle_types)) deallocate(light_vert%particle_types)
  light_vert%particle_types = [particle_gc_relativistic_id]
end subroutine setup_gyroaverage_synchrotron_light_class

!> Tools -------------------------------------------------
!> compute the gyroaverage synchrotron radiation directionality function
!> as explained in M. Hoppe, Nucl. Fuson, vol.58, p.026032, 2018
!> inputs:
!>   wavelength_crit: (real8) critical wavelength
!>   fact_1:              (real8) (gamma**3)*sqrt(((1-beta*cospsi)**3)/(0.5*beta*cospsi))
!>   fact_2:              (real8) (((1-beta*cospsi)/(beta*cospsi))**2)*&
!>                                (1-costhetap*cosmu)
!>   fact_3:              (real8) (0.5*cospsi*(sinpsi**2))/(1-beta*cospsi)
!>   dir_funct_intensity: (real8) directionality function intensity
!> outputs:
!>   dir_funct:           (real8) gyroaverage synchrotron radiation directionality funct.
subroutine compute_synchrotron_directionality_funct(wavelength,&
wavelength_crit,fact_1,fact_2,fact_3,dir_funct_intensity,dir_funct)
  use mod_besselk, only: f_besselk
  implicit none
  !> Inputs:
  real*8,intent(in) :: wavelength,wavelength_crit
  real*8,intent(in) :: fact_1,fact_2,fact_3,dir_funct_intensity
  !> Outputs:
  real*8,intent(out) :: dir_funct
  !> Variables
  real*8 :: xi,lambdac_over_lambda,besselk_1,besselk_2
  !> initialisation
  dir_funct = 0d0
  !> compute the directionality function
  lambdac_over_lambda = wavelength_crit/wavelength;
  xi = lambdac_over_lambda*fact_1;
  besselk_1 = f_besselk(onethirds,xi); besselk_2 = f_besselk(twothirds,xi);
  if(isnan(besselk_1)) return; if(isnan(besselk_2)) return;
  dir_funct = dir_funct_intensity*(lambdac_over_lambda**4)*funct_2*(&
              (besselk_2**2) + (besselk_1**2)*fact_3)
end subroutine compute_synchrotron_directionality_funct
!> -------------------------------------------------------
end module mod_gyroaverage_synchrotron_dist_light_vertices
