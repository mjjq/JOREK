module corr_neg

implicit none

contains



!> NUMERICAL IMPROVEMENT FOR CASES WHERE TEMPERATURES CLOSE TO OR BELOW ZERO CAN OCCUR.
!! 
!! PROBLEM:
!! 
!! T^-1.5 is undefined for negative temperatures (resistivity and other quantities).
!! Thus, abs(T) was used so far which may, however, still cause discontinuities
!! in the resistivity.
!! 
!! SOLUTION:
!! 
!! Replace abs(T) by smooth function:
!! 
!!   f(T) = T                             if T>L1+L2
!!   f(T) = L1 + L2 * exp((T-(L2+L1))/L2) otherwise
!!   
!! where L1 and L2 are derived from the input parameter T_1:
!! 
!!   L1 = T_1 * corr_neg_temp_coef(1)
!!   L2 = T_1 * corr_neg_temp_coef(2)
!!   
!! The default values corr_neg_temp_coef(:) = (/ 0.5, 0.5 /) can be
!! changed via the namelist input file. Alternatively, different values
!! can be provided via the optional routine parameter coef.
!!
real*8 function corr_neg_temp(val, coef, val_1)
  
  use phys_module, only: T_1, corr_neg_temp_coef
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  real*8, intent(in), optional :: coef(2)   !< Optional coefficients, if not provided the
                                            !! input parameter corr_neg_temp_coef is used instead.
  real*8, intent(in), optional :: val_1     !< Temperature value floor
  real*8 :: L1, L2, T_floor
 
  if (present(val_1)) then
    T_floor = val_1
  else 
    T_floor = T_1
  end if
 
  if ( present(coef) ) then
    L1 = T_floor * coef(1)
    L2 = T_floor * coef(2)
  else
    L1 = T_floor * corr_neg_temp_coef(1)
    L2 = T_floor * corr_neg_temp_coef(2)
  end if

  corr_neg_temp = val
  if ( val < L1 + L2 ) corr_neg_temp = L1 + L2 * exp( (val-(L1+L2)) / L2 )

end function corr_neg_temp



!> dT_corr/dT
real*8 function dcorr_neg_temp_dT(val, coef, val_1)
  
  use phys_module, only: T_1, corr_neg_temp_coef
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  real*8, intent(in), optional :: coef(2)   !< Optional coefficients, if not provided the
                                            !! input parameter corr_neg_temp_coef is used instead.
  real*8, intent(in), optional :: val_1     !< Temperature value floor
  real*8 :: L1, L2, T_floor

  if (present(val_1)) then
    T_floor = val_1
  else
    T_floor = T_1
  end if
  
  if ( present(coef) ) then
    L1 = T_floor * coef(1)
    L2 = T_floor * coef(2)
  else
    L1 = T_floor * corr_neg_temp_coef(1)
    L2 = T_floor * corr_neg_temp_coef(2)
  end if

  dcorr_neg_temp_dT = 1.d0
  if ( val < L1 + L2 ) dcorr_neg_temp_dT = exp( (val-(L1+L2)) / L2 )

end function dcorr_neg_temp_dT



!> d^2T_corr/dT^2
real*8 function d2corr_neg_temp_dT2(val, coef, val_1)
  
  use phys_module, only: T_1, corr_neg_temp_coef
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  real*8, intent(in), optional :: coef(2)   !< Optional coefficients, if not provided the
                                            !! input parameter corr_neg_temp_coef is used instead.
  real*8, intent(in), optional :: val_1     !< Temperature value floor  
  real*8 :: L1, L2, T_floor

  if (present(val_1)) then
    T_floor = val_1
  else
    T_floor = T_1
  end if
  
  if ( present(coef) ) then
    L1 = T_floor * coef(1)
    L2 = T_floor * coef(2)
  else
    L1 = T_floor * corr_neg_temp_coef(1)
    L2 = T_floor * corr_neg_temp_coef(2)
  end if

  d2corr_neg_temp_dT2 = 0.d0
  if ( val < L1 + L2 ) d2corr_neg_temp_dT2 = exp( (val-(L1+L2)) / L2 ) / L2

end function d2corr_neg_temp_dT2



!> Same for density (so far not used in element_matrix routines).
real*8 function corr_neg_dens(val, coef, val_1)
  
  use phys_module, only: rho_1, corr_neg_dens_coef
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Density value to be "corrected".
  real*8, intent(in), optional :: coef(2)   !< Optional coefficients, if not provided the
                                            !! input parameter corr_neg_temp_coef is used instead.
  real*8, intent(in), optional :: val_1     !< Density value floor  
  real*8 :: L1, L2, rho_floor

  if (present(val_1)) then
    rho_floor = val_1
  else
    rho_floor = rho_1
  end if
  
  if ( present(coef) ) then
    L1 = rho_1 * coef(1)
    L2 = rho_1 * coef(2)
  else
    L1 = rho_1 * corr_neg_dens_coef(1)
    L2 = rho_1 * corr_neg_dens_coef(2)
  end if

  corr_neg_dens = val
  if ( val < L1 + L2 ) corr_neg_dens = L1 + L2 * exp( (val-(L1+L2)) / L2 )

end function corr_neg_dens



!> dT_corr/dT
real*8 function dcorr_neg_dens_drho(val, coef, val_1)
  
  use phys_module, only: rho_1, corr_neg_dens_coef
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Density value to be "corrected".
  real*8, intent(in), optional :: coef(2)   !< Optional coefficients, if not provided the
                                            !! input parameter corr_neg_dens_coef is used instead.
  real*8, intent(in), optional :: val_1     !< Density value floor    
  real*8 :: L1, L2, rho_floor

  if (present(val_1)) then
    rho_floor = val_1
  else
    rho_floor = rho_1
  end if
  
  if ( present(coef) ) then
    L1 = rho_floor * coef(1)
    L2 = rho_floor * coef(2)
  else
    L1 = rho_floor * corr_neg_dens_coef(1)
    L2 = rho_floor * corr_neg_dens_coef(2)
  end if

  dcorr_neg_dens_drho = 1.d0
  if ( val < L1 + L2 ) dcorr_neg_dens_drho = exp( (val-(L1+L2)) / L2 )

end function dcorr_neg_dens_drho



end module corr_neg
