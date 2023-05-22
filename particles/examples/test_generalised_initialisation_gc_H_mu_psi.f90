!> The test code is used for emulating the initialise_particles_H_mu_psi
!> method using the initialise_particles_in_phase_space method
!> for guiding centers
!> inputs:
!>   n_x:          (integer) number of variables: must be 7
!>   x:            (real8)(n_x) uniform random numbers in [0,1)
!>   st:           (real8)(2) local particle coordinates 
!>   i_elm:        (integer) jorek element
!>   fields:       (fields_base) JOREK MHD field class
!>   x_min:        (real8)(n_x) lower bound uniform sampling, 1) poloidal flux
!>                 2) poloidal angle, 3) toroidal angle, 4) not used, 5) not used
!>                 6) gyro angle, 7) charge
!>   x_max:        (real8)(n_x) upper bound uniform sampling, 1) poloidal flux
!>                 2) poloidal angle, 3) toroidal angle, 4) not used, 5) not used
!>                 6) gyro angle, 7) charge
!>   n_real_param: (integer) size of the real parameters, must be: 2*n_elements+3
!>   real_param:   (real8)(2*n_elements+3) real gdf_sampler parameters:
!>                 1:n_elements              <- psi_minmax(:,1)
!>                 n_elements+1:2*n_elements <- psi_minmax(:,2)
!>                 2*n_elements+1            <- magnetic axis major radius
!>                 2*n_elements+2            <- magnetic axis vertical position
!>                 2*n_elements+3            <- particle mass
!>   n_int_param:  (integer) size of the integer parameters, must be: 1
!>   int_param:    (integer)(1) integer parameters: 1 <- if 1 paralle velocity is included
!> outputs:
!>   i_elm:        (integer) jorek element
!>   st:           (real8)(2) local particle coordinates 
!>   x:            (real8)(n_x) guiding center variables: 1) R, 2) Z, 3) phi,
!>                 4) energy, 5) magnetic moment, 6) gyro angle, 7) charge
subroutine gdf_psi_H_mu_sampler(n_x,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param)
  use constants,          only: MU_ZERO,EL_CHG,ATOMIC_MASS_UNIT
  use phys_module,        only: central_density
  use mod_model_settings, only: var_T,var_Vpar
  use mod_fields,         only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in) :: n_x,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(n_x),intent(in)            :: x_min,x_max
  real*8,dimension(:),allocatable,intent(:)   :: real_param
  class(fields_base),intent(in)               :: fields
  !> Inputs-Outputs:
  integer,intent(inout)               :: i_elm
  real*8,dimension(2),intent(inout)   :: st
  real*8,dimension(n_x),intent(inout) :: x
  !> Variables:
  real*8 :: psi,U,normB,R,Z,electric_potential,u,temperature_ev
  real*8,dimension(1) :: P
  real*8,dimension(3) :: B,E

  !> compute the psi,theta,phi coordinates
  x(1:3) = x_min(1:3)+(x_max(1:3)-x_min(1:3))*x(1:3)
  call find_theta_psi(fields%node_list,fields%element_list,[real_param(1:fields%element_list%n_elements),\
  real_param(fields%element_list%n_elements+1:2*fields%element_list%n_elements)],\
  x(2),x(1),x(3),real_param(2*fields%element_list%n_elements+1),\
  real_param(2*fields%element_list%n_elements+2),i_elm,st(1),st(2),x(1),x(2))
  if(i_elm.gt.0) then
    call fields%calc_EBpsiU(time,i_elm,st,x(3),E,B,psi,electric_potential); normB = norm2(B);
    !> compute the temperature
    call interp_PRZ(fields%node_list,fields%element_list,i_elm,[var_T],1,st(1),st(2),x(3),P,R,Z)
    temperature_ev = P(1)/(2d0*MU_ZERO*central_density*1d20*EL_CHG) !< temperature [eV]
#ifdef WITH_TiTe
    temperature_ev = 2d0*temperature_ev
#endif
    x(4) = temperature_ev*0.5d0*sample_chi_squared_3(x(4))
    u = 2*mod(x(5),0.5d0)
    x(5) = sign(x(4)*(2d0*u-u**2),x(5)-0.5d0)
    !> include parallel velocity if required
    if(int_param(1).eq.1) then
      call interp_PRZ(fields%node_list,fields%element_list,i_elm,[var_Vpar],1,st(1),st(2),x(3),P,R,Z)
      P(1) = P(1)*sqrt(real_param(2*fields%element_list%n_elements+3)*ATOMIC_MASS_UNIT/EL_CHG)
      x(4) = x(4) + P(1)*(P(1) + 2d0*sign(sqrt(x(4)-x(5)),x(5)))
    endif
      x(5) = x(5)/normB
      !> uniform distribution for both the gyro angle and the charge distribution
      !> very easy to modify for using the coronal charge equilibrium instead
      x(6:7) = x_min(6:7) + (x_max(6:7)-x_min(6:7))*x(6:7) 
  endif
end subroutine gdf_psi_H_mu_sampler
