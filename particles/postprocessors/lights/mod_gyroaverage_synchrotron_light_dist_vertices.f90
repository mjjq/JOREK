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
real*8,parameter :: onethird=1d0/3d0
real*8,parameter :: twothird=2d0/3d0
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
!>                 1-3: x,y,z components of the magnetic field
subroutine compute_synchrotron_mhd_fields(light_vert,fields,&
particle_in,time_id,mass,mhd_fields)
  use mod_fields,                only: fields_base
  use mod_particle_types,        only: particle_base
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
#ifdef UNIT_TESTS_AFIELDS
  !> used only for unit testing
  use mod_particle_common_test_tools, only: compute_test_E_B_fields
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
  !> Variables:
  real*8 :: psi,U
  real*8,dimension(3) :: Efield
#ifndef UNIT_TESTS_AFIELDS
  !> compute the JOREK magnetic field
  call fields%calc_EBpsiU(light_vert%times(time_id),particle_in%i_elm,&
  particle_in%st,particle_in%x(3),Efields,mhd_fields,psi,U)
#else
  !> compute the analytical magnetic field
  call compute_test_E_B_fields(particle_in%x,Efields,mhd_fields)
#endif
  mhd_fields = vector_cylindrical_to_cartesian(particle_in%x(3),mhd_fields)
end subroutine compute_synchrotron_mhd_fields

!> compute the gyroaverage synchrotron light properties from a
!> relativistic guiding center particle.
!> inputs:
!> outputs:
subroutine compute_gyroaverage_synchrotron_light_properties()
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
  light_vert%n_property_vertex = 5; light_vert%n_mhd = 3;
  light_vert%n_particle_types = 1;
  if(allocated(light_vert%particle_types)) deallocate(light_vert%particle_types)
  light_vert%particle_types = [particle_gc_relativistic_id]
end subroutine setup_gyroaverage_synchrotron_light_class


!> Tools -------------------------------------------------
!> -------------------------------------------------------
end module mod_gyroaverage_synchrotron_dist_light_vertices
