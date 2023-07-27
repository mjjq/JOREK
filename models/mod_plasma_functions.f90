!> Module containing plasma functions t
module mod_plasma_functions 
  
  use phys_module, only: eta_T_dependent, T_min 
  use mod_model_settings, only:  with_impurities
    
  implicit none
  
  private
  public resistivity
  
  contains
  
  
  !> Determine resistivity
  pure subroutine resistivity(eta_0, T_raw, T_corr, T_max, T0, Z_eff, eta_T,                       & 
                              dZ_eff_dT, dZ_eff_dr0, dZ_eff_drimp0, dr0_corr_dn, drimp0_corr_dn,   & 
                              deta_dT, d2eta_d2T, deta_dr0, deta_drimp0) 

    implicit none
    
    real*8, intent(in)             :: eta_0            ! central resistivity
    real*8, intent(in)             :: T_raw            ! temperature without correction
    real*8, intent(in)             :: T_corr           ! corrected temperature > 0
    real*8, intent(in)             :: T_max            ! max temperature to use in the function
    real*8, intent(in)             :: T0               ! central temperature at equilibrium
    real*8, intent(in)             :: Z_eff            ! effective charge (only used with_impurities at the moment)
    real*8, intent(out)            :: eta_T            ! output resistivity
    real*8, optional, intent(in)   :: dZ_eff_dT        ! Derivative of Zeff w.r.t. the temperature 
    real*8, optional, intent(in)   :: dZ_eff_dr0       ! Derivative of Zeff w.r.t. the total density 
    real*8, optional, intent(in)   :: dZ_eff_drimp0    ! Derivative of Zeff w.r.t. the impurity density
    real*8, optional, intent(in)   :: dr0_corr_dn      ! Derivative of density correction  
    real*8, optional, intent(in)   :: drimp0_corr_dn   ! Derivative of impurity density correction  
    real*8, optional, intent(out)  :: deta_dT          ! 1st derivative with respect to the temperature
    real*8, optional, intent(out)  :: d2eta_d2T        ! 2nd derivative with respect to the temperature
    real*8, optional, intent(out)  :: deta_dr0         ! 1st derivative with respect to the total density
    real*8, optional, intent(out)  :: deta_drimp0      ! 1st derivative with respect to the timpurity density

    !--- Local parameters
    real*8 :: eta_coef, deta_coef_dZeff

    if (present(deta_dT))       deta_dT     = 0.d0
    if (present(d2eta_d2T))     d2eta_d2T   = 0.d0
    if (present(deta_dr0))      deta_dr0    = 0.d0
    if (present(deta_drimp0))   deta_drimp0 = 0.d0

    if ( eta_T_dependent .and. T_corr <= T_max ) then
      eta_T     =   eta_0 * (T_corr/T0)**(-1.5d0)
      if (present(deta_dT))    deta_dT   = - eta_0 * (1.5d0)  * T_corr**(-2.5d0) * T0**(1.5d0)
      if (present(d2eta_d2T))  d2eta_d2T =   eta_0 * (3.75d0) * T_corr**(-3.5d0) * T0**(1.5d0)
    else if (eta_T_dependent .and. T_corr > T_max) then
      eta_T     = eta_0 * (T_max/T0)**(-1.5d0)
    else
      eta_T     = eta_0
    end if

    if ( eta_T_dependent .and. (T_raw .lt. T_min) ) then
      eta_T       = eta_0     * (max(T_raw, T_min)/T0)**(-1.5d0)
      if (present(deta_dT))     deta_dT    = 0.d0
      if (present(d2eta_d2T))   d2eta_d2T  = 0.d0
    endif

    if (with_impurities) then 

      eta_coef     = Z_eff*(1.+1.198*Z_eff+0.222*Z_eff**2)/(1.+2.966*Z_eff+0.753*Z_eff**2)
      eta_coef     = eta_coef / ((1.+1.198+0.222)/(1.+2.966+0.753))

      deta_coef_dZeff = (1.+1.198*Z_eff+0.222*Z_eff**2)/(1.+2.966*Z_eff+0.753*Z_eff**2)
      deta_coef_dZeff = deta_coef_dZeff + Z_eff*(1.198+2.*0.222*Z_eff)/(1.+2.966*Z_eff+0.753*Z_eff**2)
      deta_coef_dZeff = deta_coef_dZeff - Z_eff*(1.+1.198*Z_eff+0.222*Z_eff**2)*(2.966+2.*0.753*Z_eff)/((1.+2.966*Z_eff+0.753*Z_eff**2)**2)
      deta_coef_dZeff = deta_coef_dZeff / ((1.+1.198+0.222)/(1.+2.966+0.753))

      if ( eta_T_dependent ) then
        eta_T       = eta_T * eta_coef
        if (present(deta_dr0) .and. present(dZ_eff_dr0) .and. present(dr0_corr_dn)) then
          deta_dr0    = eta_T * deta_coef_dZeff * dZ_eff_dr0 * dr0_corr_dn
        endif
        if (present(deta_drimp0) .and. present(dZ_eff_drimp0) .and. present(drimp0_corr_dn)) then
            deta_drimp0 = eta_T * deta_coef_dZeff * dZ_eff_drimp0 * drimp0_corr_dn
        endif
        if (present(deta_dT) .and. present(dZ_eff_dT)) then
          deta_dT     = deta_dT * eta_coef + eta_T * deta_coef_dZeff * dZ_eff_dT
        endif
      end if

    endif

  end subroutine resistivity


end module mod_plasma_functions
