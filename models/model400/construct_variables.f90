!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------- Compute the RZ-coordinates and the Jacobians ------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
subroutine ELM_build_RZ_and_Jacobians(element, nodes, ms, mt)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_RZ_and_Jacobians

  ! --- Modules
  use parameters    
  use basis_at_gaussian
  use equation_variables
  use data_structure
  
  implicit none
  
  ! --- Routine Variables
  type (type_element)	      :: element
  type (type_node)	      :: nodes(n_vertex_max)
  integer		      :: ms, mt
  
  ! --- Internal Variables
  integer		      :: i, j
      
  ! --- Empty before integration
  x_g = 0.d0 ; x_s = 0.d0 ; x_t = 0.d0; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0
  y_g = 0.d0 ; y_s = 0.d0 ; y_t = 0.d0; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0

  ! --- Integrate
  do i=1,n_vertex_max
    do j=1,n_order+1

      x_g  = x_g  + nodes(i)%x(j,1) * element%size(i,j) * H   (i,j,ms,mt)
      x_s  = x_s  + nodes(i)%x(j,1) * element%size(i,j) * H_s (i,j,ms,mt)
      x_t  = x_t  + nodes(i)%x(j,1) * element%size(i,j) * H_t (i,j,ms,mt)

      x_ss = x_ss + nodes(i)%x(j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
      x_st = x_st + nodes(i)%x(j,1) * element%size(i,j) * H_st(i,j,ms,mt)
      x_tt = x_tt + nodes(i)%x(j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

      y_g  = y_g  + nodes(i)%x(j,2) * element%size(i,j) * H   (i,j,ms,mt)
      y_s  = y_s  + nodes(i)%x(j,2) * element%size(i,j) * H_s (i,j,ms,mt)
      y_t  = y_t  + nodes(i)%x(j,2) * element%size(i,j) * H_t (i,j,ms,mt)

      y_ss = y_ss + nodes(i)%x(j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
      y_st = y_st + nodes(i)%x(j,2) * element%size(i,j) * H_st(i,j,ms,mt)
      y_tt = y_tt + nodes(i)%x(j,2) * element%size(i,j) * H_tt(i,j,ms,mt)
    
    enddo
  enddo
  
  BigR    = x_g
  BigR_x  = 1.d0
  
  ! --- Jacobians
  xjac    = x_s*y_t - x_t*y_s

  xjac_x  = (x_ss* y_t**2 - y_ss*x_t*y_t - 2.d0*x_st*y_s*y_t   &       
	   + y_st*(x_s*y_t + x_t*y_s)			       &
	   + x_tt* y_s**2 - y_tt*x_s*y_s) / xjac
	 
  xjac_y  = (y_tt* x_s**2 - x_tt*y_s*x_s - 2.d0*y_st*x_t*x_s   &       
	   + x_st*(y_t*x_s + y_s*x_t)			       &
	   + y_ss* x_t**2 - x_ss*y_t*x_t) / xjac

  return

end subroutine ELM_build_RZ_and_Jacobians








!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------- Compute the variables for the equations --------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
subroutine ELM_build_variables(element, nodes, ms, mt, i_plane)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_variables

  ! --- Modules
  use parameters    
  use basis_at_gaussian
  use equation_variables
  use data_structure
  use phys_module
  
  implicit none
  
  ! --- Routine variables
  type (type_element)	      :: element
  type (type_node)	      :: nodes(n_vertex_max)
  integer		      :: ms, mt, i_plane
  
  ! --- Internal variables
  integer		      :: i, j, k, i_tor
      
  ! --- Empty before integration
  ps0	= 0.d0; ps0_p	= 0.d0; ps0_s	= 0.d0; ps0_t	= 0.d0
  u0	= 0.d0; u0_p	= 0.d0; u0_s	= 0.d0; u0_t	= 0.d0
  zj0	= 0.d0; zj0_p	= 0.d0; zj0_s	= 0.d0; zj0_t	= 0.d0
  w0	= 0.d0; w0_p	= 0.d0; w0_s	= 0.d0; w0_t	= 0.d0; w0_ss	 = 0.d0; w0_tt    = 0.d0; w0_st    = 0.d0
  r0	= 0.d0; r0_p	= 0.d0; r0_s	= 0.d0; r0_t	= 0.d0; r0_ss	 = 0.d0; r0_tt    = 0.d0; r0_st    = 0.d0
  T0	= 0.d0; T0_p	= 0.d0; T0_s	= 0.d0; T0_t	= 0.d0; T0_ss	 = 0.d0; T0_tt    = 0.d0; T0_st    = 0.d0; T0_pp = 0.d0
  T3	= 0.d0; T3_ss	= 0.d0; T3_tt	= 0.d0; T3_st	= 0.d0
  Vpar0 = 0.d0; Vpar0_p = 0.d0; Vpar0_s = 0.d0; Vpar0_t = 0.d0; Vpar0_ss = 0.d0; Vpar0_tt = 0.d0; Vpar0_st = 0.d0
  delta_g = 0.d0 ; delta_s = 0.d0 ; delta_t = 0.d0

  ! --- Integrate
  do i =1,n_vertex_max
    do j=1,n_order+1
      do i_tor =1,n_tor

	! --- Variable 1
	ps0	     = ps0	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	ps0_p	     = ps0_p	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	ps0_s	     = ps0_s	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	ps0_t	     = ps0_t	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)

	! --- Variable 2
	u0	     = u0	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	u0_p	     = u0_p	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	u0_s	     = u0_s	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	u0_t	     = u0_t	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)

	! --- Variable 3
	zj0	     = zj0	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	zj0_p	     = zj0_p	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	zj0_s	     = zj0_s	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	zj0_t	     = zj0_t	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)

	! --- Variable 4
	w0	     = w0	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_p	     = w0_p	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	w0_s	     = w0_s	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_t	     = w0_t	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_ss	     = w0_ss	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_tt	     = w0_tt	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_st	     = w0_st	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)

	! --- Variable 5
	r0	     = r0	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_p	     = r0_p	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	r0_s	     = r0_s	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_t	     = r0_t	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_ss	     = r0_ss	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_tt	     = r0_tt	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_st	     = r0_st	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)

	! --- Variable 6
	T0	     = T0	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_p	     = T0_p	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	T0_s	     = T0_s	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_t	     = T0_t	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_ss	     = T0_ss	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_tt	     = T0_tt	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_st	     = T0_st	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_pp	     = T0_pp	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)
  
	! --- Toroidally localised hyper-parallel-conductivity
	if(i_tor .ne. 1) then
	  T3	     = T3	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	  T3_ss      = T3_ss	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	  T3_tt      = T3_tt	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	  T3_st      = T3_st	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	endif

	! --- Variable 7
	Vpar0	     = Vpar0	  + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_p      = Vpar0_p    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	Vpar0_s      = Vpar0_s    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_t      = Vpar0_t    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_ss     = Vpar0_ss   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_tt     = Vpar0_tt   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_st     = Vpar0_st   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)

	! --- Deltas
	do k=1,n_var
	  delta_g(k) = delta_g(k) + nodes(i)%deltas(i_tor,j,k) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	  delta_s(k) = delta_s(k) + nodes(i)%deltas(i_tor,j,k) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	  delta_t(k) = delta_t(k) + nodes(i)%deltas(i_tor,j,k) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	enddo			
  
      enddo
    enddo
  enddo

  ! --- Variable 1
  ps0_x    = (   y_t * ps0_s - y_s * ps0_t ) / xjac
  ps0_y    = ( - x_t * ps0_s + x_s * ps0_t ) / xjac
  
  ! --- Variable 2
  u0_x     = (   y_t * u0_s - y_s * u0_t ) / xjac
  u0_y     = ( - x_t * u0_s + x_s * u0_t ) / xjac
  vv2	   = BigR**2 *  ( u0_x * u0_x + u0_y *u0_y  )
  
  ! --- Variable 3
  zj0_x    = (   y_t * zj0_s - y_s * zj0_t ) / xjac
  zj0_y    = ( - x_t * zj0_s + x_s * zj0_t ) / xjac
  
  ! --- Variable 4
  w0_x     = (   y_t * w0_s - y_s * w0_t ) / xjac
  w0_y     = ( - x_t * w0_s + x_s * w0_t ) / xjac
  
  ! --- Variable 5
  r0_x     = (   y_t * r0_s - y_s * r0_t ) / xjac
  r0_y     = ( - x_t * r0_s + x_s * r0_t ) / xjac
  r0_hat   = BigR**2 * r0
  r0_x_hat = 2.d0 * BigR * BigR_x  * r0 + BigR**2 * r0_x
  r0_y_hat = BigR**2 * r0_y
  
  ! --- Variable 6
  T0_x     = (   y_t * T0_s  - y_s * T0_t ) / xjac
  T0_y     = ( - x_t * T0_s  + x_s * T0_t ) / xjac
  T0_xx    = (T0_ss * y_t**2 - 2.d0*T0_st * y_s*y_t + T0_tt * y_s**2  & 	    
	      + T0_s * (y_st*y_t - y_tt*y_s )			      &    
	      + T0_t * (y_st*y_s - y_ss*y_t ) )      / xjac**2        & 	
	    - xjac_x * (T0_s * y_t - T0_t * y_s)     / xjac**2
  T0_yy    = (T0_ss * x_t**2 - 2.d0*T0_st * x_s*x_t + T0_tt * x_s**2  & 	    
	      + T0_s * (x_st*x_t - x_tt*x_s )			      &    
	      + T0_t * (x_st*x_s - x_ss*x_t ) )      / xjac**2        & 	
	    - xjac_y * (- T0_s * x_t + T0_t * x_s )  / xjac**2
  T0_xy    = (- T0_ss * y_t*x_t - T0_tt * x_s*y_s		      &
	      + T0_st * (y_s*x_t  + y_t*x_s  )  		      &        
	      - T0_s  * (x_st*y_t - x_tt*y_s )  		      &    
	      - T0_t  * (x_st*y_s - x_ss*y_t )  )    / xjac**2        & 	
	    - xjac_x * (- T0_s * x_t + T0_t * x_s )  / xjac**2
  T0_ps0_x = T0_xx * ps0_y - T0_xy * ps0_x + T0_x * ps0_xy - T0_y * ps0_xx
  T0_ps0_y = T0_xy * ps0_y - T0_yy * ps0_x + T0_x * ps0_yy - T0_y * ps0_xy
  
  ! --- Variable 7
  Vpar0_x  = (   y_t * Vpar0_s - y_s * Vpar0_t ) / xjac
  Vpar0_y  = ( - x_t * Vpar0_s + x_s * Vpar0_t ) / xjac
  
  ! --- Deltas
  delta_u_x= (   y_t * delta_s(2) - y_s * delta_t(2) ) / xjac
  delta_u_y= ( - x_t * delta_s(2) + x_s * delta_t(2) ) / xjac
  
  ! --- Pressure
  P0	   = r0   * T0
  P0_x     = r0_x * T0 + r0 * T0_x
  P0_y     = r0_y * T0 + r0 * T0_y
  P0_s     = r0_s * T0 + r0 * T0_s
  P0_t     = r0_t * T0 + r0 * T0_t
  P0_p     = r0_p * T0 + r0 * T0_p
  
  ! --- Magnetic field amplitude (squared)
  BB2	    = (F0*F0 + ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2
  
  return

end subroutine ELM_build_variables







!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------- Compute the diffusivities and source ---------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
subroutine ELM_build_diffusivities_and_sources(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint, i_plane)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_diffusivities_and_sources

  ! --- Modules
  use parameters    
  use basis_at_gaussian
  use phys_module
  use equation_variables
  use data_structure
  
  implicit none
  
  ! --- Routine variables
  type (type_element)	      :: element
  type (type_node)	      :: nodes(n_vertex_max)
  logical		      :: xpoint2
  integer		      :: xcase2
  integer		      :: i_plane
  real*8		      :: R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint(2)
  
  ! --- Internal variables
  integer		      :: id
  real*8		      :: psi_norm, psi_D
  real*8		      :: atn_D, datn_D, atn_D_n, pol_D, dpol_D, D_min
  real*8		      :: prof(1:2),Diff(1:2,1:10)
      
  ! -------------------------------------
  ! --- Temperature dependent resistivity
  ! -------------------------------------
  if ( eta_T_dependent ) then
    if ( T0 .lt. T_1 ) then
      eta_T	=   eta   * (T_1/T_0)**(-1.5d0)
      deta_dT	= - eta   * (1.5d0)  * T_1**(-2.5d0) * T_0**(1.5d0)
      d2eta_d2T =   eta   * (3.75d0) * T_1**(-3.5d0) * T_0**(1.5d0)
    else
      eta_T	=   eta   * (T0 /T_0)**(-1.5d0)
      deta_dT	= - eta   * (1.5d0)  * T0 **(-2.5d0) * T_0**(1.5d0)
      d2eta_d2T =   eta   * (3.75d0) * T0 **(-3.5d0) * T_0**(1.5d0)
    endif
  else
    eta_T     = eta
    deta_dT   = 0.d0
    d2eta_d2T = 0.d0
  end if
  
  ! -----------------------------------
  ! --- Temperature dependent viscosity
  ! -----------------------------------
  if ( visco_T_dependent ) then
    if ( T0 .lt. T_1 ) then
      visco_T	=   visco * (T_1/T_0)**(-1.5d0)
      dvisco_dT = - visco * (1.5d0)  * T_1**(-2.5d0) * T_0**(1.5d0)
    else
      visco_T	=   visco * (T0 /T_0)**(-1.5d0)
      dvisco_dT = - visco * (1.5d0)  * T0 **(-2.5d0) * T_0**(1.5d0)
    endif
  else
    visco_T   = visco
    dvisco_dT = 0.d0
  end if
  
  ! -------------------------------------------------------------
  ! --- D_perp and K_perp profiles (for fixed pedestal gradients)
  ! -------------------------------------------------------------
  ! --- First need psi_norm
  psi_norm = (ps0 - psi_axis)/(psi_bnd - psi_axis)
  if (xpoint2) then
    if ((psi_norm .lt. 1.d0) .and. (y_g .lt. Z_xpoint(1)) .and. (xcase2 .ne. 2)) then
      psi_norm = 2.d0 - psi_norm
    endif
    if ((psi_norm .lt. 1.d0) .and. (y_g .gt. Z_xpoint(2)) .and. (xcase2 .ne. 1)) then
      psi_norm = 2.d0 - psi_norm
    endif
  endif

  ! --- Take values from input file (rho_coef, T_coef...) 
  do id = 1,10
    if ((id .eq. 7) .or. (id .eq. 8) .or. (id .eq. 9)) then
      Diff(1,id) = rho_coef(id-6) ; Diff(2,id) = T_coef(id-6)
    else
      Diff(1,id) = D_perp(id)	  ; Diff(2,id) = ZK_perp(id)
    endif
  enddo      
  if (Diff(1,10) .eq. 1.d0) then 
    Diff(1,4) = rho_coef(4) ; Diff(1,5) = rho_coef(5) 
  endif
  if (Diff(2,10) .eq. 1.d0) then 
    Diff(2,4) = T_coef(4)   ; Diff(2,5) = T_coef(5) 
  endif
  
  ! --- Build profiles
  do id = 1,2
    if (psi_norm .gt. Diff(id,5)) then
      if (id .eq. 1) then
	psi_D = 2.d0*Diff(id,5) - psi_norm
      else
	psi_D = Diff(id,5)
      endif
    else 
      psi_D = psi_norm
    endif
    if (psi_norm .lt. 0.5d0) then
      psi_D = 0.5d0
    endif 
    if (Diff(id,7) .ge. 0.d0) then
      Diff(id,7) = -0.1d0
    endif 
    if ((xcase2 .ne. 2) .and. (id .eq. 1))  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z_xpoint(1)-y_g)/0.1d0))
    if ((xcase2 .ne. 1) .and. (id .eq. 1))  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((y_g-Z_xpoint(2))/0.1d0))
    atn_D    = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(id,5))/Diff(id,4))
    datn_D   =       - 0.5d0 / cosh((psi_D-Diff(id,5))/Diff(id,4))**2.d0 /(Diff(id,4)*(psi_bnd - psi_axis))
    pol_D    = 1 + Diff(id,7)*psi_D    + Diff(id,8)*psi_D**2.d0      + Diff(id,9)*psi_D**3.d0
    dpol_D   =    (Diff(id,7)	    + 2.d0*Diff(id,8)*psi_D	+ 3.d0*Diff(id,9)*psi_D**2.d0)/(psi_bnd - psi_axis)
    D_min    = 1.d0/( -(1+Diff(id,7)*Diff(id,5)+Diff(id,8)*Diff(id,5)**2.d0+Diff(id,9)*Diff(id,5)**3.d0) * 0.5d0/(Diff(id,4)*(psi_bnd - psi_axis))&
	       + 0.5d0 * (Diff(id,7)	 + 2.d0*Diff(id,8)*Diff(id,5)+ 3.d0*Diff(id,9)*Diff(id,5)**2.d0)/(psi_bnd - psi_axis) )

    prof(id) = (1.d0-Diff(id,10)) * ( Diff(id,1) * (1.d0-Diff(id,2)+Diff(id,2)*(0.5d0 - 0.5d0 * tanh((psi_norm-Diff(id,5))/Diff(id,4)))) &
				    + Diff(id,6) * (0.5d0 - 0.5d0 * tanh((-psi_norm+Diff(id,5)+Diff(id,3))/Diff(id,4)))) &
		    +Diff(id,10)  * ( Diff(id,1) / (dpol_D*atn_D + pol_D*datn_D) / D_min ) &
				  * (1 + Diff(id,6) - Diff(id,6) * tanh(-(psi_norm-(1+4*Diff(id,4)))/Diff(id,4)))   !higher Kperp in SOL
  enddo
  
  ! --- Allocate profiles to corresponding names
  D_prof  = prof(1)
  ZK_prof = prof(2) 
  
  ! --- Avoid negative density 
  if ( r0 .lt. rho_1 ) then
    D_prof = D_prof * 1.d2
  endif       
  if ( T0 .lt. T_1 ) then
    ZK_prof = ZK_prof * 1.d0
  endif

  ! -----------------------------------------------------
  ! --- Parallel conductivity profiles (Braginskii model)
  ! -----------------------------------------------------
  if ( T0 .lt. T_1 ) then
    K_par   = ZK_par	       * T_1**(2.5d0) 
    dZK_par = ZK_par * (2.5d0) * T_1**(1.5d0)
  else
    K_par   = ZK_par	       * T0 **(2.5d0) 
    dZK_par = ZK_par * (2.5d0) * T0 **(1.5d0)
  endif
  
  ! -------------------------
  ! --- Hyper diffusivitities
  ! -------------------------
  eta_numm	 = eta_num			 ! hyper-resistivity
  visco_numm	 = visco_num			 ! hyper-viscosity
  visco_par_numm = visco_par_num		 ! hyper-viscosity
  D_perp_numm	 = D_perp_num			 ! hyper-diffusivity
  K_perp_numm	 = 0.d0!ZK_perp_num		 ! hyper-conductivity
  K_perp_numm2   = ZK_perp_num  		 ! hyper-conductivity
  K_par_num	 = 0.d0!1.d-10  		 ! hyper-parallel-conductivity
  K_par_num2	 = 0.d0!1.d-10  		 ! hyper-parallel-conductivity
  
  if (psi_norm .lt. 0.4d0) eta_numm   = eta_numm   * 1.d2
  if (psi_norm .lt. 0.4d0) visco_numm = visco_numm * 1.d2
  if (psi_norm .lt. 0.2d0) visco_numm = visco_numm * 1.d2
  if (psi_norm .lt. 0.2d0) eta_numm   = eta_numm   * 1.d2

  ! -------------------------------------------------------------------
  ! --- Heating, current and particle source (the same for all i_plane)
  ! -------------------------------------------------------------------
  if (i_plane .eq. 1) then
    call current(xpoint2, xcase2, x_g,y_g, Z_xpoint, ps0,psi_axis,psi_bnd,current_source)

    !call sources(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd,particle_source,heat_source_i,heat_source_e)
    ! --- New source profile: source with exactly the same profile as the initial equilibirum profiles.
    !call density(    xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
    !		     zn,dn_dpsi,  dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)

    !call temperature(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
    !		     zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)

    !particle_source = particlesource * ( zn - r0 )
    !heat_source    = heatsource     * ( zT - T0 )

    !particle_source = particle_source * ( 0.5d0 - 0.5d0 * tanh((psi_norm-0.99)/0.005) )
    !heat_source    = heat_source     * ( 0.5d0 - 0.5d0 * tanh((psi_norm-0.99)/0.005) )
    particle_source = particlesource * (0.5d0 - 0.5d0 * tanh((psi_norm-0.5)/0.005) )
    heat_source     = heatsource     * (0.5d0 - 0.5d0 * tanh((psi_norm-0.5)/0.005) ) 
  endif
  
  return

end subroutine ELM_build_diffusivities_and_sources










!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!-------------------------------- Compute the test functions (ie. the Bezier polynomials) -------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
subroutine ELM_build_test_functions(element, nodes, ms, mt, i_plane, i_vertex, i_order, i_tor, &
				    vv, vv_s,  vv_t,	     vv_p,  vv_x,  vv_y,	       &
        				vv_ss, vv_tt, vv_st, vv_pp, vv_xx, vv_yy, vv_xy        )
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_test_functions

  ! --- Modules
  use parameters    
  use basis_at_gaussian
  use equation_variables
  use data_structure
  
  implicit none
  
  ! --- Routine variables
  type (type_element)	      :: element
  type (type_node)	      :: nodes(n_vertex_max)
  integer		      :: ms, mt, i_plane, i_vertex, i_order, i_tor
  real*8		      :: vv, vv_s,  vv_t,	  vv_p,  vv_x,  vv_y
  real*8		      ::     vv_ss, vv_tt, vv_st, vv_pp, vv_xx, vv_yy, vv_xy
  
  ! --- Internal variables
  real*8, dimension(n_tor,n_plane) :: HHZ, HHZ_p, HHZ_pp
  
  ! --- Toroidal functions	      
  if (n_tor .gt. 3) then
    HHZ   (i_tor,i_plane) = 1.d0
    HHZ_p (i_tor,i_plane) = 1.d0
    HHZ_pp(i_tor,i_plane) = 1.d0
  else
    HHZ   (i_tor,i_plane) = HZ   (i_tor,i_plane)
    HHZ_p (i_tor,i_plane) = HZ_p (i_tor,i_plane)
    HHZ_pp(i_tor,i_plane) = HZ_pp(i_tor,i_plane)
  endif
  
  ! --- Test functions  	      
  vv	= H   (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ	(i_tor,i_plane)
  vv_s  = H_s (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ	(i_tor,i_plane)
  vv_t  = H_t (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ	(i_tor,i_plane)
  vv_p  = H   (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ_p (i_tor,i_plane)
  vv_pp = H   (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ_pp(i_tor,i_plane)

  vv_ss = H_ss(i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ	(i_tor,i_plane)
  vv_tt = H_tt(i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ	(i_tor,i_plane)
  vv_st = H_st(i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ	(i_tor,i_plane)

  vv_x = (  y_t * vv_s - y_s * vv_t ) / xjac
  vv_y = (- x_t * vv_s + x_s * vv_t ) / xjac

  vv_xx = (vv_ss * y_t**2 - 2.d0*vv_st * y_s*y_t + vv_tt * y_s**2   &	      
	 + vv_s * (y_st*y_t - y_tt*y_s )			    &	 
	 + vv_t * (y_st*y_s - y_ss*y_t ) )    / xjac**2 	    &	  
	 - xjac_x * (vv_s * y_t - vv_t * y_s) / xjac**2

  vv_yy = (vv_ss * x_t**2 - 2.d0*vv_st * x_s*x_t + vv_tt * x_s**2   &	      
	 + vv_s * (x_st*x_t - x_tt*x_s )			    &	 
	 + vv_t * (x_st*x_s - x_ss*x_t ) )	/ xjac**2	    &	  
	 - xjac_y * (- vv_s * x_t + vv_t * x_s ) / xjac**2

  vv_xy = (- vv_ss * y_t*x_t - vv_tt * x_s*y_s  		&
	   + vv_st * (y_s*x_t  + y_t*x_s  )			&	 
	   - vv_s  * (x_st*y_t - x_tt*y_s )			&	   
	   - vv_t  * (x_st*y_s - x_ss*y_t )  )       / xjac**2  &		
	   - xjac_x * (- vv_s * x_t + vv_t * x_s )   / xjac**2  	      
  
  return

end subroutine ELM_build_test_functions


