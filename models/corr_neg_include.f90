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

real*8 function corr_neg_temp1(val)
  !$omp declare simd 
  
  use phys_module, only: T_1, corr_neg_temp_coef
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  
  real*8 :: L1, L2
  
  L1 = T_1 * corr_neg_temp_coef(1)
  L2 = T_1 * corr_neg_temp_coef(2)

  corr_neg_temp1 = val
  if ( val < L1 + L2 ) corr_neg_temp1 = L1 + L2 * exp( (val-(L1+L2)) / L2 )

end function corr_neg_temp1

real*8 function corr_neg_temp2(val, coef)
  !$omp declare simd uniform(coef)
! With uniform, we declare thet coeff should be the same for all vector elements.
! Is this correct?
  use phys_module, only: T_1

  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  real*8, intent(in) :: coef(2)   !< Optional coefficients, if not provided the
                                            !! input parameter corr_neg_temp2_coef is used instead.
  
  real*8 :: L1, L2
  
  L1 = T_1 * coef(1)
  L2 = T_1 * coef(2)

  corr_neg_temp2 = val
  if ( val < L1 + L2 ) corr_neg_temp2 = L1 + L2 * exp( (val-(L1+L2)) / L2 )

end function corr_neg_temp2

!> Same for density (so far not used in element_matrix routines).
real*8 function corr_neg_dens1(val)
!$omp declare simd
  use phys_module, only: rho_1, corr_neg_dens_coef
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Density value to be "corrected".
  
  real*8 :: L1, L2
  
  L1 = rho_1 * corr_neg_dens_coef(1)
  L2 = rho_1 * corr_neg_dens_coef(2)

  corr_neg_dens1 = val
  if ( val < L1 + L2 ) corr_neg_dens1 = L1 + L2 * exp( (val-(L1+L2)) / L2 )

end function corr_neg_dens1

!We cannot have optional argument for a vector funct, therefore we overload it
real*8 function corr_neg_dens2(val, coef)
!$omp declare simd uniform(coef)
! With uniform, we declare thet coeff should be the same for all vector elements.
! Is this correct?
  use phys_module, only: rho_1, corr_neg_dens_coef
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Density value to be "corrected".
  real*8, intent(in) :: coef(2)   !< Optional coefficients, if not provided the
                                            !! input parameter corr_neg_temp_coef is used instead.
  
  real*8 :: L1, L2
  
  L1 = rho_1 * coef(1)
  L2 = rho_1 * coef(2)

  corr_neg_dens2 = val
  if ( val < L1 + L2 ) corr_neg_dens2 = L1 + L2 * exp( (val-(L1+L2)) / L2 )

end function corr_neg_dens2
