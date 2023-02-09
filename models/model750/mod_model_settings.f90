!> Basic model-dependend hard-coded run parameters.
module mod_model_settings

  implicit none

  logical, parameter :: with_TiTe       = .true.
  logical, parameter :: with_neutrals   = .false.
  logical, parameter :: with_impurities = .true.
  logical, parameter :: with_Vpar       = .false.
  logical, parameter :: with_refluid    = .false.
  
! ##################################################################################################
! ####  @USERS: This file should not be modified ###################################################
! ##################################################################################################

! The following line is needed by ./util/config.sh:
! #SETTINGS# with_vpar with_TiTe with_neutrals with_impurities

  integer, parameter :: jorek_model    = 750       !< JOREK physics model

  logical, parameter :: hydrodynamics   = .false.
  logical, parameter :: reduced_MHD     = .false.
  logical, parameter :: full_MHD        = .true.
  
  logical, parameter :: model_family      = .true. !< Is this 
  character(len=42)  :: base_mod_descr    = 'Model family for tokamak full MHD'
  
  ! --- extensions to it
  integer, parameter :: n_mod_ext            = 4      !< Number of model extensions
  integer, parameter :: i_ext_TiTe           = 1
  integer, parameter :: i_ext_neutrals       = 2
  integer, parameter :: i_ext_impurities     = 3
  integer, parameter :: i_ext_refluid        = 4
  logical, parameter :: with_ext(n_mod_ext) = (/ with_TiTe, with_neutrals, with_impurities, with_refluid /)

  ! --- variable indices for the base model  
  integer, parameter :: var_A3   = 1                       ! place of variable psi/mag pot 3               (ps or A3)
  integer, parameter :: var_AR   = 2                       ! place of variable mag pot  1                  (AR)
  integer, parameter :: var_AZ   = 3                       ! place of variable mag pot  2                  (AZ)
  integer, parameter :: var_uR   = 4                       ! place of variable velocity 1                  (UR)
  integer, parameter :: var_uZ   = 5                       ! place of variable velocity 2                  (UZ)
  integer, parameter :: var_up   = 6                       ! place of variable velocity 3                  (Up)
  integer, parameter :: var_rho  = 7                       ! place of variable density                     (r or rho)
  ! --- variable indices for the model extensions
  integer, parameter :: n_var_base   = var_rho
  integer, parameter :: n_var_T      = merge(n_var_base+2, n_var_base+1, with_TiTe )  ! place of variable temperature                 (T)
  integer, parameter :: if_neutrals  = merge(1, 0, with_neutrals )
  integer, parameter :: if_imp       = merge(1, 0, with_impurities )
  integer, parameter :: var_T        = merge(0, n_var_T, with_TiTe )  ! place of variable temperature                 (T)
  integer, parameter :: var_Ti       = merge(n_var_T-1, 0, with_TiTe )! place of variable ion temperature             (Ti)
  integer, parameter :: var_Te       = merge(n_var_T, 0, with_TiTe )  ! place of variable electron temperature        (Te)
  integer, parameter :: var_rhon     = merge(n_var_T+1, 0, if_neutrals==1 )  ! place of variable neutral density (rhon)
  integer, parameter :: var_rhoimp   = merge( merge(n_var_T+1, 0, if_imp==1 )+1, n_var_T+1, if_neutrals==1)  ! place of variable impurity density (rhoimp)
  integer, parameter :: n_var        = n_var_T+if_neutrals+if_imp

  integer, parameter :: n_var_TiTe        = sum(merge( (/1/), (/0/), with_TiTe      ))
  integer, parameter :: n_var_neutrals    = sum(merge( (/1/), (/0/), with_neutrals  )) 
  integer, parameter :: n_var_impurities  = sum(merge( (/1/), (/0/), with_impurities)) 
  integer, parameter :: n_var_refluid     = sum(merge( (/1/), (/0/), with_refluid   )) !### not yet
  integer, parameter :: n_var_ext(n_mod_ext) = (/ n_var_TiTe, n_var_neutrals, n_var_impurities, n_var_refluid /)

  ! --- variables not relevant to this model
  integer, parameter :: var_psi  = 1                       ! place of variable psi/mag pot 3               (ps or A3)  
  integer, parameter :: var_u    = 0                       ! place of variable velocity stream function    (u)
  integer, parameter :: var_zj   = 0                       ! place of variable toroidal current density    (zj)
  integer, parameter :: var_w    = 0                       ! place of variable vorticity                   (w)
  integer, parameter :: var_Vpar = 0                       ! place of variable parallel velocity           (Vpar)
  integer, parameter :: var_jec  = 0                       ! place of variable ECCD current                (jec)
  integer, parameter :: var_jec1 = 0                       ! place of variable ECCD current #1             (jec1)
  integer, parameter :: var_jec2 = 0                       ! place of variable ECCD current #2             (jec2)
  integer, parameter :: var_nre  = 0                       ! place of variable for RE number density       (nre)

  !> element_matrix and element_matrix_fft combined into a single one?
  logical, parameter :: unified_element_matrix = .true.

contains
  
!> Is the extension available?
elemental pure logical function ext_available(i_ext)

  implicit none

  ! --- Function parameters
  integer, intent(in) :: i_ext


  ! --- Preset to .true. unless invalid extension number
  if ( (i_ext < 1) .or. (i_ext > n_mod_ext) ) then
    ext_available = .false.
  else
    ext_available = .true.
  end if

  ! --- exceptions for not implemented
  if ( i_ext == i_ext_TiTe ) then
    ext_available = .true.
  else if ( i_ext == i_ext_neutrals ) then
    ext_available = .true.
  else if ( i_ext == i_ext_impurities ) then
    ext_available = .true.
  else if ( i_ext == i_ext_refluid ) then
    ext_available = .false.
  end if

end function ext_available

!> Are two extensions compatible?
pure logical function ext_compatible(i_ext1, i_ext2)

  implicit none

  ! --- Function parameters
  integer, intent(in) :: i_ext1, i_ext2

  ! --- Local variables
  integer :: iext1, iext2

  ! --- sort extension indices, since compatibility is symmetric
  iext1 = min(i_ext1, i_ext2)
  iext2 = max(i_ext1, i_ext2)

  ! --- preset to .true. unless wrong indices were provided
  if ( (iext1 < 1) .or. (iext2 > n_mod_ext) ) then
    ext_compatible = .false.
  else
    ext_compatible = .true.
  end if

  ! --- exceptions for compatibility
  if ( ( iext1 == i_ext_TiTe ) .and. ( iext2 == i_ext_refluid ) ) then
    ext_compatible = .false. ! ### just an example
  end if

end function ext_compatible
  
end module mod_model_settings
