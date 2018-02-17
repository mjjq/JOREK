!> Module containing functions to determine particle and heat diffusivities
module diffusivities
  
  use mod_parameters,  only: jorek_model
  use phys_module, only: num_d_perp, D_perp, num_d_perp_x, num_d_perp_y, num_d_perp_len,           &
                         num_zk_perp, ZK_perp, num_zk_perp_x, num_zk_perp_y, num_zk_perp_len,      &
			 xpoint, xcase, rho_0, rho_coef, T_coef
  use profiles,    only: interpolProf
    
  implicit none
  
  private
  public get_dperp, get_zkperp
  
  
  interface get_dperp
    module procedure get_dperp1
    module procedure get_dperp2
  end interface get_dperp
  
  
  interface get_zkperp
    module procedure get_zkperp1
    module procedure get_zkperp2
  end interface get_zkperp
  
  
  
  contains
  
  
  
  !> Determine perpendicular particle diffusivity, D_perp, as a function of Psi_N
  real*8 function get_dperp1(psin)
    
    implicit none
    
    real*8, intent(in) :: psin
    
    if ( num_d_perp ) then
      
      get_dperp1 = interpolProf(num_d_perp_x, num_d_perp_y, num_d_perp_len, psin)
      
    else
      
      get_dperp1 = D_perp(1) * ( (1.d0-D_perp(2)) +                                                 &
        D_perp(2)*(0.5d0 - 0.5d0*tanh((psin-D_perp(5))/D_perp(4))) )
      
      if ( jorek_model >= 300 ) then
        
        get_dperp1 = get_dperp1 + D_perp(6)*D_perp(2) *                                              &
          ((0.5d0 - 0.5d0*tanh((-psin+D_perp(5)+D_perp(3)) /D_perp(4))))
        
      end if
      
    end if
    
  end function get_dperp1
  
  
  
  !> Determine perpendicular heat diffusivity, ZK_perp, as a function of Psi_N
  real*8 function get_zkperp1(psin)
    
    implicit none
    
    real*8, intent(in) :: psin
    
    if ( num_zk_perp ) then
      
      get_zkperp1 = interpolProf(num_zk_perp_x, num_zk_perp_y, num_zk_perp_len, psin)
      
    else
      
      get_zkperp1 = ZK_perp(1) * ( (1.d0-ZK_perp(2)) +                                              &
        ZK_perp(2) *(0.5d0 - 0.5d0*tanh((psin-ZK_perp(5))/ZK_perp(4))) )
      
      if ( jorek_model >= 300 ) then
        
        get_zkperp1 = get_zkperp1 + ZK_perp(6)*ZK_perp(2) *                                          &
          ((0.5d0 - 0.5d0*tanh((-psin+ZK_perp(5)+ZK_perp(3)) /ZK_perp(4))))
        
      end if
      
    end if
    
  end function get_zkperp1
  
  
  
  !> Determine perpendicular particle diffusivity, D_perp, as a function of Psi_N
  real*8 function get_dperp2(psi, psi_norm, psi_axis, psi_bnd, Z, Z_xpoint)
    
    implicit none
    
    real*8, intent(in) :: psi, psi_norm, psi_axis, psi_bnd, Z, Z_xpoint(2)
    real*8	       :: psi_D
    real*8	       :: atn_D, datn_D, atn_D_n, pol_D, dpol_D, D_min
    real*8	       :: Diff(1:5)
    
    ! --- Numerical profile
    if ( num_d_perp ) then
      
      get_dperp2 = interpolProf(num_d_perp_x, num_d_perp_y, num_d_perp_len, psi_norm)
      
    ! --- Normal (analytical) profiles
    else
      
      ! --- Old profile, note, Diff(10) is our switch: =0 => use old profile, =1 => use new profile
      if (D_perp(10) .ne. 1.d0) then
        
	atn_D = (0.5d0 - 0.5d0*tanh((psi_norm-D_perp(5))/D_perp(4)))
        get_dperp2 = D_perp(1) * ( (1.d0-D_perp(2)) + D_perp(2) * atn_D )
        if (jorek_model >= 300) then
          atn_D = ((0.5d0 - 0.5d0*tanh((-psi_norm+D_perp(5)+D_perp(3)) /D_perp(4))))
	  get_dperp2 = get_dperp2 + D_perp(6) * D_perp(2) * atn_D
        endif
	
      ! --- Adapted profile, going like Dperp ~ 1/grad(rho) to obtain constand perp flux in pedestal
      else  
        ! --- Take values from input file (rho_coef, T_coef...) 
      	Diff(1) = rho_coef(1)
      	Diff(2) = rho_coef(2)
      	Diff(3) = rho_coef(3)
      	Diff(4) = rho_coef(4)
	Diff(5) = rho_coef(5)
      	
	! --- Correct for hollow profiles (otherwise D_perp > infinity)
	if (Diff(1) .ge. 0.d0) Diff(1) = -0.1d0

      	! --- First compute the normalisation factor (this is the profile at a given psi)
      	! --- Note that if you choose psi_D=Diff(5) you normalise with the minimal value of D_perp
      	psi_D	 = 0.8       ! At the pedestal top, roughly
      	!psi_D	 = Diff(5)   ! At the point of steepest gradient (ie. will be the smallest D_perp)
      	! --- The tanh part of the initial profile and its derivative 
      	atn_D	 = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(5))/Diff(4))
      	datn_D   =	 - 0.5d0 / cosh((psi_D-Diff(5))/Diff(4))**2.d0 /(Diff(4)*(psi_bnd - psi_axis))
      	! --- The polynomial part of the initial profile and its derivative
      	pol_D	 = 1 + Diff(1)*psi_D    + Diff(2)*psi_D**2.d0	 + Diff(3)*psi_D**3.d0
      	dpol_D   =    (Diff(1)	   + 2.d0*Diff(2)*psi_D     + 3.d0*Diff(3)*psi_D**2.d0)/(psi_bnd - psi_axis)
      	! --- The normalisation factor
      	D_min	 = 1.d0 / ( dpol_D*atn_D + pol_D*datn_D )

      	! --- Then we recompute the profile derivatives for the actual local psi_norm
      	! --- Be careful at core and in SOL, where gradients are ~0 => infinite 1/grad(p)
      	! --- Take a symetric profile at from middle of pedestal outwards, because SOL is too dangerous (gradient too close to zero...)
	if (psi_norm .gt. Diff(5)) then
      	  psi_D = 2.d0*Diff(5) - psi_norm
      	else
      	  psi_D = psi_norm
      	endif
      	if (psi_norm .lt. 0.5d0) psi_D = 0.5d0
      	if (xcase .ne. 2)  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z_xpoint(1)-Z)/0.1d0))
      	if (xcase .ne. 1)  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z-Z_xpoint(2))/0.1d0))
      	! --- The tanh part of the initial profile and its derivative 
      	atn_D	 = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(5))/Diff(4))
      	datn_D   =	 - 0.5d0 / cosh((psi_D-Diff(5))/Diff(4))**2.d0 /(Diff(4)*(psi_bnd - psi_axis))
      	! --- The polynomial part of the initial profile and its derivative
      	pol_D	 = 1 + Diff(1)*psi_D    + Diff(2)*psi_D**2.d0	 + Diff(3)*psi_D**3.d0
      	dpol_D   =    (Diff(1)	   + 2.d0*Diff(2)*psi_D     + 3.d0*Diff(3)*psi_D**2.d0)/(psi_bnd - psi_axis)

      	! --- Build the profile (goes like 1/grad(density_profile) in the pedestal region
      	get_dperp2 = D_perp(1) / (dpol_D*atn_D + pol_D*datn_D) / D_min 
        
      end if
      
    end if
    
  end function get_dperp2
  
  
  
  !> Determine perpendicular heat diffusivity, ZK_perp, as a function of Psi_N
  real*8 function get_zkperp2(psi, psi_norm, psi_axis, psi_bnd, Z, Z_xpoint)
    
    implicit none
    
    real*8, intent(in) :: psi, psi_norm, psi_axis, psi_bnd, Z, Z_xpoint(2)
    real*8	       :: psi_D
    real*8	       :: atn_D, datn_D, atn_D_n, pol_D, dpol_D, D_min
    real*8	       :: Diff(1:5)
    real*8             :: rho_norm
    real*8             :: zn, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz
    
    ! --- Numerical profile
    if ( num_zk_perp ) then
      
      get_zkperp2 = interpolProf(num_zk_perp_x, num_zk_perp_y, num_zk_perp_len, psi_norm)
      
    ! --- Normal (analytical) profiles
    else
      
      ! --- Old profile, note, Diff(10) is our switch: =0 => use old profile, =1 => use new profile
      if (ZK_perp(10) .ne. 1.d0) then
        
	atn_D = (0.5d0 - 0.5d0*tanh((psi_norm-ZK_perp(5))/ZK_perp(4)))
	get_zkperp2 = ZK_perp(1) * ( (1.d0-ZK_perp(2)) + ZK_perp(2) * atn_D )
        if ( jorek_model >= 300 ) then
          atn_D = ((0.5d0 - 0.5d0*tanh((-psi_norm+ZK_perp(5)+ZK_perp(3)) /ZK_perp(4))))
          get_zkperp2 = get_zkperp2 + ZK_perp(6) * ZK_perp(2) * atn_D
        endif

      ! --- Adapted profile, going like Dperp ~ 1/grad(rho) to obtain constand perp flux in pedestal
      else
      	
        ! --- Take values from input file (rho_coef, T_coef...) 
      	Diff(1) = T_coef(1)
      	Diff(2) = T_coef(2)
      	Diff(3) = T_coef(3)
      	Diff(4) = T_coef(4)
	Diff(5) = T_coef(5)
      	
	! --- Correct for hollow profiles (otherwise D_perp > infinity)
	if (Diff(1) .ge. 0.d0) Diff(1) = -0.1d0

      	! --- First compute the normalisation factor (this is the profile at a given psi)
      	! --- Note that if you choose psi_D=Diff(5) you normalise with the minimal value of D_perp
      	psi_D	 = 0.8       ! At the pedestal top, roughly
      	!psi_D	 = Diff(5)   ! At the point of steepest gradient (ie. will be the smallest D_perp)
      	! --- The tanh part of the initial profile and its derivative 
      	atn_D	 = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(5))/Diff(4))
      	datn_D   =	 - 0.5d0 / cosh((psi_D-Diff(5))/Diff(4))**2.d0 /(Diff(4)*(psi_bnd - psi_axis))
      	! --- The polynomial part of the initial profile and its derivative
      	pol_D	 = 1 + Diff(1)*psi_D    + Diff(2)*psi_D**2.d0	 + Diff(3)*psi_D**3.d0
      	dpol_D   =    (Diff(1)	   + 2.d0*Diff(2)*psi_D     + 3.d0*Diff(3)*psi_D**2.d0)/(psi_bnd - psi_axis)
      	! --- The normalisation factor
      	D_min	 = 1.d0 / ( dpol_D*atn_D + pol_D*datn_D )

      	! --- Then we recompute the profile derivatives for the actual local psi_norm
      	! --- Be careful at core and in SOL, where gradients are ~0 => infinite 1/grad(p)
      	! --- Take a symetric profile at from middle of pedestal outwards, because SOL is too dangerous (gradient too close to zero...)
	if (psi_norm .gt. Diff(5)) then
      	  psi_D = 2.d0*Diff(5) - psi_norm
      	else
      	  psi_D = psi_norm
      	endif
      	if (psi_norm .lt. 0.5d0) psi_D = 0.5d0
      	if (xcase .ne. 2)  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z_xpoint(1)-Z)/0.1d0))
      	if (xcase .ne. 1)  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z-Z_xpoint(2))/0.1d0))
      	! --- The tanh part of the initial profile and its derivative 
      	atn_D	 = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(5))/Diff(4))
      	datn_D   =	 - 0.5d0 / cosh((psi_D-Diff(5))/Diff(4))**2.d0 /(Diff(4)*(psi_bnd - psi_axis))
      	! --- The polynomial part of the initial profile and its derivative
      	pol_D	 = 1 + Diff(1)*psi_D    + Diff(2)*psi_D**2.d0	 + Diff(3)*psi_D**3.d0
      	dpol_D   =    (Diff(1)	   + 2.d0*Diff(2)*psi_D     + 3.d0*Diff(3)*psi_D**2.d0)/(psi_bnd - psi_axis)

      	! --- Build the profile (goes like 1/grad(density_profile) in the pedestal region
      	get_zkperp2 = ZK_perp(1) / (dpol_D*atn_D + pol_D*datn_D) / D_min 

        ! --- K-perp is the pressure diffusion, not the temperature diffusion, so need to divide by rho
      	psi_D = psi
      	call density(xpoint, xcase, Z, Z_xpoint, psi_D,psi_axis,psi_bnd, &
      		     zn,dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
      	rho_norm        = zn / rho_0 ! normalise to 1.0 in core
      	get_zkperp2 = get_zkperp2 / rho_norm
	
      end if
      
    end if
    
  end function get_zkperp2





end module diffusivities
