!> Module containing functions to determine particle and heat diffusivities
module diffusivities
  
  use parameters,  only: jorek_model
  use phys_module, only: num_d_perp, D_perp, num_d_perp_x, num_d_perp_y, num_d_perp_len,           &
    num_zk_perp, ZK_perp, num_zk_perp_x, num_zk_perp_y, num_zk_perp_len
  use profiles,    only: interpolProf
    
  implicit none
  
  private
  public get_dperp, get_zkperp
  
  
  
  contains
  
  
  
  !> Determine perpendicular particle diffusivity, D_perp, as a function of Psi_N
  real*8 function get_dperp(psin)
    
    implicit none
    
    real*8, intent(in) :: psin
    
    if ( num_d_perp ) then
      
      get_dperp = interpolProf(num_d_perp_x, num_d_perp_y, num_d_perp_len, psin)
      
    else
      
      get_dperp = D_perp(1) * ( (1.d0-D_perp(2)) +                                                 &
        D_perp(2)*(0.5d0 - 0.5d0*tanh((psin-D_perp(5))/D_perp(4))) )
      
      if ( jorek_model >= 300 ) then
        
        get_dperp = get_dperp + D_perp(6)*D_perp(2) *                                              &
          ((0.5d0 - 0.5d0*tanh((-psin+D_perp(5)+D_perp(3)) /D_perp(4))))
        
      end if
      
    end if
    
  end function get_dperp
  
  
  
  !> Determine perpendicular heat diffusivity, ZK_perp, as a function of Psi_N
  real*8 function get_zkperp(psin)
    
    implicit none
    
    real*8, intent(in) :: psin
    
    if ( num_zk_perp ) then
      
      get_zkperp = interpolProf(num_zk_perp_x, num_zk_perp_y, num_zk_perp_len, psin)
      
    else
      
      get_zkperp = ZK_perp(1) * ( (1.d0-ZK_perp(2)) +                                              &
        ZK_perp(2) *(0.5d0 - 0.5d0*tanh((psin-ZK_perp(5))/ZK_perp(4))) )
      
      if ( jorek_model >= 300 ) then
        
        get_zkperp = get_zkperp + ZK_perp(6)*ZK_perp(2) *                                          &
          ((0.5d0 - 0.5d0*tanh((-psin+ZK_perp(5)+ZK_perp(3)) /ZK_perp(4))))
        
      end if
      
    end if
    
  end function get_zkperp

end module diffusivities
