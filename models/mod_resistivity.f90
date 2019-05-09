!> Module containing functions to determine particle and heat diffusivities
module mod_resistivity 
  
  use phys_module, only: eta_T_dependent, T_0, T_min, xpoint, eta 
    
  implicit none
  
  private
  public resistivity, dresistivity_dT 
  
 
  contains
  
  
  !> Determine resistivity
  pure function resistivity(T) result(eta_T)
#if _OPENMP >= 201511
    !$omp declare simd
#endif
    implicit none
    
    real*8, intent(in) :: T
    real*8             :: eta_T

    ! --- Temperature dependent resistivity
    if ( eta_T_dependent ) then
      eta_T     = eta   * (T/T_0)**(-1.5d0)
      if ( xpoint .and. (T .lt. T_min) ) then
        eta_T   = eta   * (max(T,T_min)/T_0)**(-1.5d0)
      endif
    else
      eta_T     = eta
    end if

  end function resistivity
    
  !> Determine eta derivative 
  pure function dresistivity_dT(T) result(deta_dT)
#if _OPENMP >= 201511
    !$omp declare simd
#endif
    implicit none
    
    real*8, intent(in) :: T
    real*8             :: deta_dT

    ! --- Temperature dependent resistivity
    if ( eta_T_dependent ) then
      deta_dT   = - eta   * (1.5d0)  * T**(-2.5d0) * T_0**(1.5d0)
      if ( xpoint .and. (T .lt. T_min) ) then
        deta_dT   = 0.d0
      endif
    else
      deta_dT   = 0.d0
    end if
    
   
  end function dresistivity_dT 
  


end module mod_resistivity 
