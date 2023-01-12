!> Basic model-dependend hard-coded run parameters.
module mod_model_settings

  implicit none

! ##################################################################################################
! ####  @USERS: This file should not be modified ###################################################
! ##################################################################################################

  integer, parameter :: jorek_model    = 712       !< JOREK physics model

  logical, parameter :: hydrodynamics   = .false.
  logical, parameter :: reduced_MHD     = .false.
  logical, parameter :: full_MHD        = .true.
  
  logical, parameter :: with_TiTe       = .true.
  logical, parameter :: with_neutrals   = .false.
  logical, parameter :: with_impurities = .false.
  logical, parameter :: with_Vpar       = .false.
  logical, parameter :: with_refluid    = .false.

  integer, parameter :: n_mod_ext            = 0 !< this model is not a model family => no extensions

  integer, parameter :: var_A3   = 1                       ! place of variable psi/mag pot 3               (ps or A3)
  integer, parameter :: var_AR   = 2                       ! place of variable mag pot  1                  (AR)
  integer, parameter :: var_AZ   = 3                       ! place of variable mag pot  2                  (AZ)
  integer, parameter :: var_uR   = 4                       ! place of variable velocity 1                  (UR)
  integer, parameter :: var_uZ   = 5                       ! place of variable velocity 2                  (UZ)
  integer, parameter :: var_up   = 6                       ! place of variable velocity 3                  (Up)
  integer, parameter :: var_rho  = 7                       ! place of variable density                     (r or rho)
  integer, parameter :: var_psi  = 1                       ! place of variable psi/mag pot 3               (ps or A3)  
  integer, parameter :: var_u    = 0                       ! place of variable velocity stream function    (u)
  integer, parameter :: var_zj   = 0                       ! place of variable toroidal current density    (zj)
  integer, parameter :: var_w    = 0                       ! place of variable vorticity                   (w)
  integer, parameter :: var_Vpar = 0                       ! place of variable parallel velocity           (Vpar)
  integer, parameter :: var_jec  = 0                       ! place of variable ECCD current                (jec)
  integer, parameter :: var_jec1 = 0                       ! place of variable ECCD current #1             (jec1)
  integer, parameter :: var_jec2 = 0                       ! place of variable ECCD current #2             (jec2)
  integer, parameter :: var_nre  = 0                       ! place of variable for RE number density       (nre)

! --- variable indices for the model extensions
  integer, parameter :: n_var_base   = var_rho
  integer, parameter :: n_var_T      = merge(n_var_base+2, n_var_base+1, with_TiTe )
  integer, parameter :: if_neutrals  = merge(1, 0, with_neutrals )
  integer, parameter :: if_imp       = merge(1, 0, with_impurities )
  integer, parameter :: var_T        = merge(0, n_var_T, with_TiTe )  ! place of variable temperature                 (T)
  integer, parameter :: var_Ti       = merge(n_var_T-1, 0, with_TiTe )! place of variable ion temperature             (Ti)
  integer, parameter :: var_Te       = merge(n_var_T, 0, with_TiTe )  ! place of variable electron temperature        (Te)
  integer, parameter :: var_rhon     = merge(n_var_T+1, 0, if_neutrals==1 )  ! place of variable neutral density (rhon)
  integer, parameter :: var_rhoimp   = merge( merge(n_var_T+1, 0, if_imp==1 )+1, n_var_T+1, if_neutrals==1)  ! place of variable impurity density (rhoimp)
  integer, parameter :: n_var        = n_var_T+if_neutrals+if_imp

  !> element_matrix and element_matrix_fft combined into a single one?
  logical, parameter :: unified_element_matrix = .true.

end module mod_model_settings
