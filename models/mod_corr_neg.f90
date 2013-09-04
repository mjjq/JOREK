module corr_neg



contains



real*8 function corr_neg_temp(val)
  
  use phys_module, only: T_1, corr_neg_temp_coef
  
  real*8, intent(in) :: val

  real*8 :: L1, L2

  L1 = T_1 * corr_neg_temp_coef(1)
  L2 = T_1 * corr_neg_temp_coef(2)

  corr_neg_temp = val
  if ( val < L1 + L2 ) corr_neg_temp = L1 + L2 * exp( (val-(L1+L2)) / L2 )

end function corr_neg_temp



real*8 function corr_neg_dens(val)
  
  use phys_module, only: rho_1, corr_neg_dens_coef
  
  real*8, intent(in) :: val

  real*8 :: L1, L2

  L1 = rho_1 * corr_neg_dens_coef(1)
  L2 = rho_1 * corr_neg_dens_coef(2)

  corr_neg_dens = val
  if ( val < L1 + L2 ) corr_neg_dens = L1 + L2 * exp( (val-(L1+L2)) / L2 )

end function corr_neg_dens


end module corr_neg
