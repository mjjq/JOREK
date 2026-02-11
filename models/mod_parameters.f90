module mod_parameters

use mod_settings
use mod_model_settings

implicit none



contains



! --- Variable names
character(len=64) pure elemental function variable_names(i_var)

implicit none

integer,           intent(in) :: i_var !< Variable index
!logical, optional, intent(in) :: long  !< Long output with description?

character(len=20) :: name
character(len=64) :: descr
logical :: long2

long2 = .false.

if ( (i_var < 1) .or. (i_var > n_var) ) then
  name  = 'NONE'
  descr = 'Illegal variable index!'
  return
else if ( i_var == var_psi ) then
  name  = 'Psi'
  descr = 'poloidal magnetic flux'
else if ( i_var == var_u ) then
  name  = 'u'
  descr = 'velocity stream function'
else if ( i_var == var_zj ) then
  name  = 'zj'
  descr = 'toroidal plasma current'
else if ( i_var == var_w ) then
  name  = 'omega'
  descr = 'toroidal vorticity'
else if ( i_var == var_rho ) then
  name  = 'rho'
  descr = 'particle density'
else if ( i_var == var_T ) then
  name  = 'T'
  descr = 'temperature; sum of electrons and ions'
else if ( i_var == var_Ti ) then
  name  = 'T_i'
  descr = 'ion temperature'
else if ( i_var == var_Te ) then
  name  = 'T_e'
  descr = 'electron temperature'
else if ( i_var == var_Vpar ) then
  name  = 'v_par'
  descr = 'parallel velocity'
else if ( i_var == var_rhon ) then
  name  = 'rho_n'
  descr = 'neutral particle density'
else if ( i_var == var_rhoimp ) then
  name  = 'rho_imp'
  descr = 'impurity density'
else if ( i_var == var_nre ) then
  name  = 'n_RE'
  descr = 'runaway electron particle density'
else if ( i_var == var_AR ) then
  name  = "AR"
  descr = "R component of magnetic vector potential"
else if ( i_var == var_AZ ) then
  name  = "AZ"
  descr = "Z component of magnetic vector potential"
else if ( i_var == var_uR ) then
  name  = "UR"
  descr = "Velocity in R direction"
else if ( i_var == var_uZ ) then
  name  = "UZ"
  descr = "Velocity in Z direction"
else if ( i_var == var_up ) then
  name  = "Uphi"
  descr = "Velocity in phi direction"
else
  write(name, "(I4)") i_var
  descr = "Please check implementation"
end if

if ( .not. long2 ) then
  variable_names = trim(name)
else
  variable_names = trim(name)//' ('//trim(descr)//')'
end if

end function variable_names


character(len=64) pure elemental function variable_names_ext(i_var)

implicit none

integer,           intent(in) :: i_var !< Variable index
!logical, optional, intent(in) :: long  !< Long output with description?

character(len=20) :: name
character(len=64) :: descr
logical :: long2

long2 = .false.
#define n_var_additional 14
if ( (i_var <= n_var) .or. (i_var > (n_var + n_var_additional)) ) then
  name  = 'NONE'
  descr = 'Illegal variable index!'
  return
else if ( i_var == n_var + 1 ) then
  name  = 'BR'
  descr = 'R component of magnetic field'
else if ( i_var == n_var + 2 ) then
  name  = 'BZ'
  descr = 'Z component of magnetic field'
else if ( i_var == n_var + 3 ) then
  name  = 'Bphi'
  descr = 'Phi component of magnetic field'
else if ( i_var == n_var + 4 ) then
  name  = "A_radial"
  descr = "Radial component of magnetic vector potential"
else if ( i_var == n_var + 5 ) then
  name  = "A_poloidal"
  descr = "Poloidal component of magnetic vector potential"
else if ( i_var == n_var + 6 ) then
  name  = "U_radial"
  descr = "Radial component of velocity"
else if ( i_var == n_var + 7 ) then
  name  = "U_poloidal"
  descr = "Poloidal component of velocity"
else if ( i_var == n_var + 8 ) then
  name  = "B_radial"
  descr = "Radial component of magnetic field"
else if ( i_var == n_var + 9 ) then
  name  = "B_poloidal"
  descr = "Poloidal component of magnetic field"
else if ( i_var == n_var + 10 ) then
  name  = "JxB_radial"
  descr = "Radial component of JxB"
else if ( i_var == n_var + 11 ) then
  name  = "gradP_radial"
  descr = "Radial component of pressure gradient"
else if ( i_var == n_var + 12 ) then
  name  = "RUxB_phi"
  descr = "Phi component of R * U x B"
else if ( i_var == n_var + 13 ) then
  name  = "rho_divU"
  descr = "rho * Div U"
else if ( i_var == n_var + 14 ) then
  name  = "U_gradRho"
  descr = "U dot grad(rho)"
else
  write(name, "(I4)") i_var
  descr = "Please check implementation"
end if

if ( .not. long2 ) then
  variable_names_ext = trim(name)
else
  variable_names_ext = trim(name)//' ('//trim(descr)//')'
end if

end function variable_names_ext


#ifdef MODEL_FAMILY

!> Are all selected extensions available and compatible?
logical pure function ext_correct()

implicit none

! --- Local variables
integer :: iext1, iext2

! --- Preset to .true.
ext_correct = .true.

! --- Check binary compatibility of all extensions
do iext1 = 1, n_mod_ext
  if ( .not. with_ext(iext1) ) cycle
  ext_correct = ext_correct .and. ext_available(iext1)
  do iext2 = 1, n_mod_ext
    if ( .not. with_ext(iext2) ) cycle
      ext_correct = ext_correct .and. ext_compatible(iext1, iext2)
  end do
end do

end function ext_correct

#endif



end module mod_parameters
