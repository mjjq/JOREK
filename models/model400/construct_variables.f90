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
  use mod_parameters    
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
  
  R    = x_g
  R_x  = 1.d0
  
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
  use mod_parameters    
  use basis_at_gaussian
  use equation_variables
  use data_structure
  use phys_module
  use corr_neg
  
  implicit none
  
  ! --- Routine variables
  type (type_element)	      :: element
  type (type_node)	      :: nodes(n_vertex_max)
  integer		      :: ms, mt, i_plane
  
  ! --- Internal variables
  integer		      :: i, j, k, i_tor
      
  ! --- Empty before integration
  ps0	= 0.d0; ps0_s	= 0.d0; ps0_t	= 0.d0; ps0_ss	 = 0.d0; ps0_tt   = 0.d0; ps0_st   = 0.d0; ps0_p   = 0.d0; ps0_pp   = 0.d0
  u0	= 0.d0; u0_s	= 0.d0; u0_t	= 0.d0; u0_ss	 = 0.d0; u0_tt    = 0.d0; u0_st    = 0.d0; u0_p    = 0.d0; u0_pp    = 0.d0
  zj0	= 0.d0; zj0_s	= 0.d0; zj0_t	= 0.d0; zj0_ss	 = 0.d0; zj0_tt   = 0.d0; zj0_st   = 0.d0; zj0_p   = 0.d0; zj0_pp   = 0.d0
  w0	= 0.d0; w0_s	= 0.d0; w0_t	= 0.d0; w0_ss	 = 0.d0; w0_tt    = 0.d0; w0_st    = 0.d0; w0_p    = 0.d0; w0_pp    = 0.d0
  r0	= 0.d0; r0_s	= 0.d0; r0_t	= 0.d0; r0_ss	 = 0.d0; r0_tt    = 0.d0; r0_st    = 0.d0; r0_p    = 0.d0; r0_pp    = 0.d0
  Ti0	= 0.d0; Ti0_s	= 0.d0; Ti0_t	= 0.d0; Ti0_ss   = 0.d0; Ti0_tt   = 0.d0; Ti0_st   = 0.d0; Ti0_p   = 0.d0; Ti0_pp   = 0.d0
  Te0	= 0.d0; Te0_s	= 0.d0; Te0_t	= 0.d0; Te0_ss   = 0.d0; Te0_tt   = 0.d0; Te0_st   = 0.d0; Te0_p   = 0.d0; Te0_pp   = 0.d0
  T0	= 0.d0; T0_s	= 0.d0; T0_t	= 0.d0; T0_ss    = 0.d0; T0_tt    = 0.d0; T0_st    = 0.d0; T0_p    = 0.d0; T0_pp    = 0.d0
  Vpar0 = 0.d0; Vpar0_s = 0.d0; Vpar0_t = 0.d0; Vpar0_ss = 0.d0; Vpar0_tt = 0.d0; Vpar0_st = 0.d0; Vpar0_p = 0.d0; Vpar0_pp = 0.d0
  delta_g = 0.d0 ; delta_s = 0.d0 ; delta_t = 0.d0

  ! --- Integrate
  do i =1,n_vertex_max
    do j=1,n_order+1
      do i_tor =1,n_tor

	! --- Variable 1
	ps0	     = ps0	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	ps0_s	     = ps0_s	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	ps0_t	     = ps0_t	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	ps0_ss	     = ps0_ss	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	ps0_tt	     = ps0_tt	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	ps0_st	     = ps0_st	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	ps0_p	     = ps0_p	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	ps0_pp	     = ps0_pp	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)

	! --- Variable 2
	u0	     = u0	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	u0_s	     = u0_s	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	u0_t	     = u0_t	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	u0_ss	     = u0_ss	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	u0_tt	     = u0_tt	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	u0_st	     = u0_st	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	u0_p	     = u0_p	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	u0_pp	     = u0_pp	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)

	! --- Variable 3
	zj0	     = zj0	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	zj0_s	     = zj0_s	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	zj0_t	     = zj0_t	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	zj0_ss	     = zj0_ss	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	zj0_tt	     = zj0_tt	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	zj0_st	     = zj0_st	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	zj0_p	     = zj0_p	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	zj0_pp	     = zj0_pp	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)

	! --- Variable 4
	w0	     = w0	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_s	     = w0_s	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_t	     = w0_t	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_ss	     = w0_ss	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_tt	     = w0_tt	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_st	     = w0_st	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_p	     = w0_p	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	w0_pp	     = w0_pp	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)

	! --- Variable 5
	r0	     = r0	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_s	     = r0_s	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_t	     = r0_t	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_ss	     = r0_ss	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_tt	     = r0_tt	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_st	     = r0_st	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_p	     = r0_p	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	r0_pp	     = r0_pp	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)

	! --- Variable 6
	Ti0	     = Ti0	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Ti0_s	     = Ti0_s	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Ti0_t	     = Ti0_t	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Ti0_ss	     = Ti0_ss	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Ti0_tt	     = Ti0_tt	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Ti0_st	     = Ti0_st	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Ti0_p	     = Ti0_p	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	Ti0_pp	     = Ti0_pp	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)
  
	! --- Variable 7
	Vpar0	     = Vpar0	  + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_s      = Vpar0_s    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_t      = Vpar0_t    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_ss     = Vpar0_ss   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_tt     = Vpar0_tt   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_st     = Vpar0_st   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_p      = Vpar0_p    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	Vpar0_pp     = Vpar0_pp	  + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)

	! --- Variable 8
	Te0	     = Te0	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Te0_s	     = Te0_s	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Te0_t	     = Te0_t	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Te0_ss	     = Te0_ss	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Te0_tt	     = Te0_tt	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Te0_st	     = Te0_st	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Te0_p	     = Te0_p	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	Te0_pp	     = Te0_pp	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)
  
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
  ps0_xx   = (ps0_ss * y_t**2 - 2.d0*ps0_st * y_s*y_t + ps0_tt * y_s**2  & 	    
	      + ps0_s * (y_st*y_t - y_tt*y_s )			         &    
	      + ps0_t * (y_st*y_s - y_ss*y_t ) )       / xjac**2         & 	
	    - xjac_x * (ps0_s * y_t - ps0_t * y_s)     / xjac**2
  ps0_yy   = (ps0_ss * x_t**2 - 2.d0*ps0_st * x_s*x_t + ps0_tt * x_s**2  & 	    
	      + ps0_s * (x_st*x_t - x_tt*x_s )			         &    
	      + ps0_t * (x_st*x_s - x_ss*x_t ) )       / xjac**2         & 	
	    - xjac_y * (- ps0_s * x_t + ps0_t * x_s )  / xjac**2
  ps0_xy   = (- ps0_ss * y_t*x_t - ps0_tt * x_s*y_s		         &
	      + ps0_st * (y_s*x_t  + y_t*x_s  )  		         &        
	      - ps0_s  * (x_st*y_t - x_tt*y_s )  		         &    
	      - ps0_t  * (x_st*y_s - x_ss*y_t )  )     / xjac**2         & 	
	    - xjac_x * (- ps0_s * x_t + ps0_t * x_s )  / xjac**2

  ! --- Variable 2
  u0_x     = (   y_t * u0_s - y_s * u0_t ) / xjac
  u0_y     = ( - x_t * u0_s + x_s * u0_t ) / xjac
  u0_xx    = (u0_ss * y_t**2 - 2.d0*u0_st * y_s*y_t + u0_tt * y_s**2  & 	    
	      + u0_s * (y_st*y_t - y_tt*y_s )			      &    
	      + u0_t * (y_st*y_s - y_ss*y_t ) )      / xjac**2        & 	
	    - xjac_x * (u0_s * y_t - u0_t * y_s)     / xjac**2
  u0_yy    = (u0_ss * x_t**2 - 2.d0*u0_st * x_s*x_t + u0_tt * x_s**2  & 	    
	      + u0_s * (x_st*x_t - x_tt*x_s )			      &    
	      + u0_t * (x_st*x_s - x_ss*x_t ) )      / xjac**2        & 	
	    - xjac_y * (- u0_s * x_t + u0_t * x_s )  / xjac**2
  u0_xy    = (- u0_ss * y_t*x_t - u0_tt * x_s*y_s		      &
	      + u0_st * (y_s*x_t  + y_t*x_s  )  		      &        
	      - u0_s  * (x_st*y_t - x_tt*y_s )  		      &    
	      - u0_t  * (x_st*y_s - x_ss*y_t )  )    / xjac**2        & 	
	    - xjac_x * (- u0_s * x_t + u0_t * x_s )  / xjac**2
  vv2	   = R**2 *  ( u0_x * u0_x + u0_y *u0_y  )
  
  ! --- Variable 3
  zj0_x    = (   y_t * zj0_s - y_s * zj0_t ) / xjac
  zj0_y    = ( - x_t * zj0_s + x_s * zj0_t ) / xjac
  zj0_xx   = (zj0_ss * y_t**2 - 2.d0*zj0_st * y_s*y_t + zj0_tt * y_s**2  & 	    
	      + zj0_s * (y_st*y_t - y_tt*y_s )			         &    
	      + zj0_t * (y_st*y_s - y_ss*y_t ) )       / xjac**2         & 	
	    - xjac_x * (zj0_s * y_t - zj0_t * y_s)     / xjac**2
  zj0_yy   = (zj0_ss * x_t**2 - 2.d0*zj0_st * x_s*x_t + zj0_tt * x_s**2  & 	    
	      + zj0_s * (x_st*x_t - x_tt*x_s )			         &    
	      + zj0_t * (x_st*x_s - x_ss*x_t ) )       / xjac**2         & 	
	    - xjac_y * (- zj0_s * x_t + zj0_t * x_s )  / xjac**2
  zj0_xy   = (- zj0_ss * y_t*x_t - zj0_tt * x_s*y_s		         &
	      + zj0_st * (y_s*x_t  + y_t*x_s  )  		         &        
	      - zj0_s  * (x_st*y_t - x_tt*y_s )  		         &    
	      - zj0_t  * (x_st*y_s - x_ss*y_t )  )     / xjac**2         & 	
	    - xjac_x * (- zj0_s * x_t + zj0_t * x_s )  / xjac**2
  
  ! --- Variable 4
  w0_x     = (   y_t * w0_s - y_s * w0_t ) / xjac
  w0_y     = ( - x_t * w0_s + x_s * w0_t ) / xjac
  w0_xx    = (w0_ss * y_t**2 - 2.d0*w0_st * y_s*y_t + w0_tt * y_s**2  & 	    
	      + w0_s * (y_st*y_t - y_tt*y_s )			      &    
	      + w0_t * (y_st*y_s - y_ss*y_t ) )      / xjac**2        & 	
	    - xjac_x * (w0_s * y_t - w0_t * y_s)     / xjac**2
  w0_yy    = (w0_ss * x_t**2 - 2.d0*w0_st * x_s*x_t + w0_tt * x_s**2  & 	    
	      + w0_s * (x_st*x_t - x_tt*x_s )			      &    
	      + w0_t * (x_st*x_s - x_ss*x_t ) )      / xjac**2        & 	
	    - xjac_y * (- w0_s * x_t + w0_t * x_s )  / xjac**2
  w0_xy    = (- w0_ss * y_t*x_t - w0_tt * x_s*y_s		      &
	      + w0_st * (y_s*x_t  + y_t*x_s  )  		      &        
	      - w0_s  * (x_st*y_t - x_tt*y_s )  		      &    
	      - w0_t  * (x_st*y_s - x_ss*y_t )  )    / xjac**2        & 	
	    - xjac_x * (- w0_s * x_t + w0_t * x_s )  / xjac**2
  
  ! --- Variable 5
  rho_corr = corr_neg_dens(r0)
  r0_x     = (   y_t * r0_s - y_s * r0_t ) / xjac
  r0_y     = ( - x_t * r0_s + x_s * r0_t ) / xjac
  r0_xx    = (r0_ss * y_t**2 - 2.d0*r0_st * y_s*y_t + r0_tt * y_s**2  & 	    
	      + r0_s * (y_st*y_t - y_tt*y_s )			      &    
	      + r0_t * (y_st*y_s - y_ss*y_t ) )      / xjac**2        & 	
	    - xjac_x * (r0_s * y_t - r0_t * y_s)     / xjac**2
  r0_yy    = (r0_ss * x_t**2 - 2.d0*r0_st * x_s*x_t + r0_tt * x_s**2  & 	    
	      + r0_s * (x_st*x_t - x_tt*x_s )			      &    
	      + r0_t * (x_st*x_s - x_ss*x_t ) )      / xjac**2        & 	
	    - xjac_y * (- r0_s * x_t + r0_t * x_s )  / xjac**2
  r0_xy    = (- r0_ss * y_t*x_t - r0_tt * x_s*y_s		      &
	      + r0_st * (y_s*x_t  + y_t*x_s  )  		      &        
	      - r0_s  * (x_st*y_t - x_tt*y_s )  		      &    
	      - r0_t  * (x_st*y_s - x_ss*y_t )  )    / xjac**2        & 	
	    - xjac_x * (- r0_s * x_t + r0_t * x_s )  / xjac**2
  r0_hat   = R**2 * r0
  r0_x_hat = 2.d0 * R * R_x  * r0 + R**2 * r0_x
  r0_y_hat = R**2 * r0_y
  
  ! --- Variable 6
  Ti_corr   = corr_neg_temp(Ti0)
  Ti0_x     = (   y_t * Ti0_s  - y_s * Ti0_t ) / xjac
  Ti0_y     = ( - x_t * Ti0_s  + x_s * Ti0_t ) / xjac
  Ti0_xx    = (Ti0_ss * y_t**2 - 2.d0*Ti0_st * y_s*y_t + Ti0_tt * y_s**2 & 	    
	      + Ti0_s * (y_st*y_t - y_tt*y_s )			         &    
	      + Ti0_t * (y_st*y_s - y_ss*y_t ) )        / xjac**2        & 	
	     - xjac_x * (Ti0_s * y_t - Ti0_t * y_s)     / xjac**2
  Ti0_yy    = (Ti0_ss * x_t**2 - 2.d0*Ti0_st * x_s*x_t + Ti0_tt * x_s**2 & 	    
	      + Ti0_s * (x_st*x_t - x_tt*x_s )			         &    
	      + Ti0_t * (x_st*x_s - x_ss*x_t ) )        / xjac**2        & 	
	     - xjac_y * (- Ti0_s * x_t + Ti0_t * x_s )  / xjac**2
  Ti0_xy    = (- Ti0_ss * y_t*x_t - Ti0_tt * x_s*y_s		         &
	       + Ti0_st * (y_s*x_t  + y_t*x_s  )  		         &        
	       - Ti0_s  * (x_st*y_t - x_tt*y_s )  		         &    
	       - Ti0_t  * (x_st*y_s - x_ss*y_t )  )       / xjac**2      & 	
	       - xjac_x * (- Ti0_s * x_t + Ti0_t * x_s )  / xjac**2
  Ti0_ps0_x = Ti0_xx * ps0_y - Ti0_xy * ps0_x + Ti0_x * ps0_xy - Ti0_y * ps0_xx
  Ti0_ps0_y = Ti0_xy * ps0_y - Ti0_yy * ps0_x + Ti0_x * ps0_yy - Ti0_y * ps0_xy

  
  ! --- Variable 7
  Vpar0_x  = (   y_t * Vpar0_s - y_s * Vpar0_t ) / xjac
  Vpar0_y  = ( - x_t * Vpar0_s + x_s * Vpar0_t ) / xjac
  Vpar0_xx = (Vpar0_ss * y_t**2 - 2.d0*Vpar0_st * y_s*y_t + Vpar0_tt * y_s**2  &	 
	      + Vpar0_s * (y_st*y_t - y_tt*y_s )			       &	
	      + Vpar0_t * (y_st*y_s - y_ss*y_t ) )         / xjac**2           &      
	    - xjac_x * (Vpar0_s * y_t - Vpar0_t * y_s)     / xjac**2
  Vpar0_yy = (Vpar0_ss * x_t**2 - 2.d0*Vpar0_st * x_s*x_t + Vpar0_tt * x_s**2  &	 
	      + Vpar0_s * (x_st*x_t - x_tt*x_s )			       &	
	      + Vpar0_t * (x_st*x_s - x_ss*x_t ) )         / xjac**2           &      
	    - xjac_y * (- Vpar0_s * x_t + Vpar0_t * x_s )  / xjac**2
  Vpar0_xy = (- Vpar0_ss * y_t*x_t - Vpar0_tt * x_s*y_s 		       &
	      + Vpar0_st * (y_s*x_t  + y_t*x_s  )  		               &	    
	      - Vpar0_s  * (x_st*y_t - x_tt*y_s )  		               &	
	      - Vpar0_t  * (x_st*y_s - x_ss*y_t )  )       / xjac**2           &      
	    - xjac_x * (- Vpar0_s * x_t + Vpar0_t * x_s )  / xjac**2
  
  ! --- Variable 8
  Te_corr   = corr_neg_temp(Te0) ! For use in eta(T), visco(T), ...
  Te0_x     = (   y_t * Te0_s  - y_s * Te0_t ) / xjac
  Te0_y     = ( - x_t * Te0_s  + x_s * Te0_t ) / xjac
  Te0_xx    = (Te0_ss * y_t**2 - 2.d0*Te0_st * y_s*y_t + Te0_tt * y_s**2 & 	    
	      + Te0_s * (y_st*y_t - y_tt*y_s )			         &    
	      + Te0_t * (y_st*y_s - y_ss*y_t ) )        / xjac**2        & 	
	     - xjac_x * (Te0_s * y_t - Te0_t * y_s)     / xjac**2
  Te0_yy    = (Te0_ss * x_t**2 - 2.d0*Te0_st * x_s*x_t + Te0_tt * x_s**2 & 	    
	      + Te0_s * (x_st*x_t - x_tt*x_s )			         &    
	      + Te0_t * (x_st*x_s - x_ss*x_t ) )        / xjac**2        & 	
	     - xjac_y * (- Te0_s * x_t + Te0_t * x_s )  / xjac**2
  Te0_xy    = (- Te0_ss * y_t*x_t - Te0_tt * x_s*y_s		         &
	       + Te0_st * (y_s*x_t  + y_t*x_s  )  		         &        
	       - Te0_s  * (x_st*y_t - x_tt*y_s )  		         &    
	       - Te0_t  * (x_st*y_s - x_ss*y_t )  )        / xjac**2     & 	
	       - xjac_x * (- Te0_s * x_t + Te0_t * x_s )  / xjac**2
  Te0_ps0_x = Te0_xx * ps0_y - Te0_xy * ps0_x + Te0_x * ps0_xy - Te0_y * ps0_xx
  Te0_ps0_y = Te0_xy * ps0_y - Te0_yy * ps0_x + Te0_x * ps0_yy - Te0_y * ps0_xy
  
  ! --- Deltas
  delta_u_x  = (   y_t * delta_s(2) - y_s * delta_t(2) ) / xjac
  delta_u_y  = ( - x_t * delta_s(2) + x_s * delta_t(2) ) / xjac
  delta_ps_x = (   y_t * delta_s(1) - y_s * delta_t(1) ) / xjac
  delta_ps_y = ( - x_t * delta_s(1) + x_s * delta_t(1) ) / xjac
  
  ! --- Total temperature
  T0	   = (Ti0    + Te0)
  T0_x     = (Ti0_x  + Te0_x)
  T0_y     = (Ti0_y  + Te0_y)
  T0_s     = (Ti0_s  + Te0_s)
  T0_t     = (Ti0_t  + Te0_t)
  T0_p     = (Ti0_p  + Te0_p)
  T0_pp    = (Ti0_pp + Te0_pp)
  T0_xx    = (Ti0_xx + Te0_xx)
  T0_yy    = (Ti0_yy + Te0_yy)
  T0_xy    = (Ti0_xy + Te0_xy)
  
  ! --- Pressure
  P0	   = r0    * (Ti0 + Te0)
  P0_x     = r0_x  * (Ti0 + Te0) + r0 * (Ti0_x + Te0_x)
  P0_y     = r0_y  * (Ti0 + Te0) + r0 * (Ti0_y + Te0_y)
  P0_s     = r0_s  * (Ti0 + Te0) + r0 * (Ti0_s + Te0_s)
  P0_t     = r0_t  * (Ti0 + Te0) + r0 * (Ti0_t + Te0_t)
  P0_p     = r0_p  * (Ti0 + Te0) + r0 * (Ti0_p + Te0_p)
  P0_pp    = r0_pp * (Ti0 + Te0) + r0 * (Ti0_pp + Te0_pp) + 2.d0 * r0_p * (Ti0_p + Te0_p)
  P0_xx    = r0_xx * (Ti0 + Te0) + r0 * (Ti0_xx + Te0_xx) + 2.d0 * r0_x * (Ti0_x + Te0_x)
  P0_yy    = r0_yy * (Ti0 + Te0) + r0 * (Ti0_yy + Te0_yy) + 2.d0 * r0_y * (Ti0_y + Te0_y)
  P0_xy    = r0_xy * (Ti0 + Te0) + r0 * (Ti0_xy + Te0_xy) + r0_x * (Ti0_y + Te0_y) + r0_y * (Ti0_x + Te0_x)
  
  ! --- Ion Pressure
  Pi0	   = r0    * Ti0
  Pi0_x    = r0_x  * Ti0 + r0 * Ti0_x
  Pi0_y    = r0_y  * Ti0 + r0 * Ti0_y
  Pi0_s    = r0_s  * Ti0 + r0 * Ti0_s
  Pi0_t    = r0_t  * Ti0 + r0 * Ti0_t
  Pi0_p    = r0_p  * Ti0 + r0 * Ti0_p
  Pi0_pp   = r0_pp * Ti0 + r0 * Ti0_pp + 2.d0 * r0_p * Ti0_p
  Pi0_xx   = r0_xx * Ti0 + r0 * Ti0_xx + 2.d0 * r0_x * Ti0_x
  Pi0_yy   = r0_yy * Ti0 + r0 * Ti0_yy + 2.d0 * r0_y * Ti0_y
  Pi0_xy   = r0_xy * Ti0 + r0 * Ti0_xy + r0_x * Ti0_y + r0_y * Ti0_x
  
  ! --- Ion Pressure
  Pe0	   = r0    * Te0
  Pe0_x    = r0_x  * Te0 + r0 * Te0_x
  Pe0_y    = r0_y  * Te0 + r0 * Te0_y
  Pe0_s    = r0_s  * Te0 + r0 * Te0_s
  Pe0_t    = r0_t  * Te0 + r0 * Te0_t
  Pe0_p    = r0_p  * Te0 + r0 * Te0_p
  Pe0_pp   = r0_pp * Te0 + r0 * Te0_pp + 2.d0 * r0_p * Te0_p
  Pe0_xx   = r0_xx * Te0 + r0 * Te0_xx + 2.d0 * r0_x * Te0_x
  Pe0_yy   = r0_yy * Te0 + r0 * Te0_yy + 2.d0 * r0_y * Te0_y
  Pe0_xy   = r0_xy * Te0 + r0 * Te0_xy + r0_x * Te0_y + r0_y * Te0_x
  
  ! --- Magnetic field amplitude (squared)
  BB2	    = (F0*F0 + ps0_x * ps0_x + ps0_y * ps0_y )/R**2
  
  return

end subroutine ELM_build_variables







!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------- Compute the diffusivities and source ---------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
subroutine ELM_build_diffusivities_and_sources(element, nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, i_plane)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_diffusivities_and_sources

  ! --- Modules
  use mod_parameters    
  use basis_at_gaussian
  use phys_module
  use equation_variables
  use data_structure
  use mod_bootstrap_functions
  
  implicit none
  
  ! --- Routine variables
  type (type_element)	      :: element
  type (type_node)	      :: nodes(n_vertex_max)
  logical		      :: xpoint2
  integer		      :: xcase2
  integer		      :: i_plane
  real*8		      :: minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
  
  ! --- Internal variables
  integer		      :: id
  real*8		      :: psi_norm, psi_D
  real*8		      :: atn_D, datn_D, atn_D_n, pol_D, dpol_D, D_min
  real*8		      :: prof(1:3),Diff(1:3,1:10)
  real*8                      :: Ti_min_Kpar, Te_min_Kpar
  real*8		      :: V_source, dV_dpsi2, dV_dz2, dV_dpsi_dz, dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz
  real*8		      :: distance_xpoint
  logical, parameter	      :: avoid_xpoint = .true.
  real*8		      :: slope, offset, Z_tmp
  real*8		      :: Rline1, Zline1
  real*8		      :: Rline2, Zline2
  real*8		      :: Rline3, Zline3
  real*8		      :: Rline4, Zline4
  real*8                      :: target_buffer_width, tan_width
  real*8                      :: drho_dpsi, grad_psi
  real*8		      :: zTi_x, zTi_y
  real*8		      :: zTe_x, zTe_y
  real*8		      :: zn_x,  zn_y
  real*8		      :: Jb_0
      
  ! -------------------------------------
  ! --- Temperature dependent resistivity
  ! -------------------------------------
  if ( eta_T_dependent ) then
    if ( Te0 .lt. Te_1 ) then
      eta_Te	 =   eta   * (Te_1/Te_0)**(-1.5d0)
      deta_dTe	 =   0.d0 
      d2eta_d2Te =   0.d0
    else
      eta_Te	 =   eta   * (Te0 /Te_0)**(-1.5d0)
      deta_dTe	 = - eta   * (1.5d0)  * Te0 **(-2.5d0) * Te_0**(1.5d0)
      d2eta_d2Te =   eta   * (3.75d0) * Te0 **(-3.5d0) * Te_0**(1.5d0)
    endif
  else
    eta_Te     = eta
    deta_dTe   = 0.d0
    d2eta_d2Te = 0.d0
  end if
  
  ! -----------------------------------
  ! --- Temperature dependent viscosity
  ! -----------------------------------
  if ( visco_T_dependent ) then
    if ( Te0 .lt. Te_1 ) then
      visco_Te	 =   visco * (Te_1/Te_0)**(-1.5d0)
      dvisco_dTe =   0.d0
    else
      visco_Te	 =   visco * (Te0 /Te_0)**(-1.5d0)
      dvisco_dTe = - visco * (1.5d0)  * Te0 **(-2.5d0) * Te_0**(1.5d0)
    endif
  else
    visco_Te   = visco
    dvisco_dTe = 0.d0
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
      Diff(1,id) = rho_coef(id-6) ; Diff(2,id) = Ti_coef(id-6) ; Diff(3,id) = Te_coef(id-6)
    else
      Diff(1,id) = D_perp(id)	  ; Diff(2,id) = ZK_i_perp(id) ; Diff(3,id) = ZK_e_perp(id)
    endif
  enddo      
  if (Diff(1,10) .eq. 1.d0) then 
    Diff(1,4) = rho_coef(4) ; Diff(1,5) = rho_coef(5) 
  endif
  if (Diff(2,10) .eq. 1.d0) then 
    Diff(2,4) = Ti_coef(4)  ; Diff(2,5) = Ti_coef(5) 
  endif
  if (Diff(3,10) .eq. 1.d0) then 
    Diff(3,4) = Te_coef(4)  ; Diff(3,5) = Te_coef(5) 
  endif
  
  ! --- Build profiles
  do id = 1,3
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
    if (xcase2 .ne. 2)  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z_xpoint(1)-y_g)/0.1d0))
    if (xcase2 .ne. 1)  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((y_g-Z_xpoint(2))/0.1d0))
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
  Ki_prof = prof(2) 
  Ke_prof = prof(3) 
  
  ! -----------------------------------------------------
  ! --- Parallel conductivity profiles (Braginskii model)
  ! -----------------------------------------------------
  Ti_min_Kpar = 5.d0 * Ti_1
  if ( Ti0 .lt. Ti_min_Kpar ) then
    Ki_par  = K_i_par	        * Ti_min_Kpar**(2.5d0) 
    dKi_par = 0.d0 
  else
    Ki_par  = K_i_par	        * Ti0 **(2.5d0) 
    dKi_par = K_i_par * (2.5d0) * Ti0 **(1.5d0)
  endif
  Te_min_Kpar = 5.d0 * Te_1
  if ( Te0 .lt. Te_min_Kpar ) then
    Ke_par  = K_e_par	        * Te_min_Kpar**(2.5d0) 
    dKe_par = 0.d0
  else
    Ke_par  = K_e_par	        * Te0 **(2.5d0) 
    dKe_par = K_e_par * (2.5d0) * Te0 **(1.5d0)
  endif
  
 
  ! -------------------------
  ! --- Hyper diffusivitities
  ! -------------------------
  call temperature_e(xpoint2,xcase2, y_g,Z_xpoint, ps0,psi_axis,psi_bnd, zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2,dTe_dpsi2_dz)
  eta_numm	 = eta_num       * (zTe/Te_0)**(-1.5d0)
  if (eta_numm .gt. 1.d-10) eta_numm = 1.d-10 
  visco_numm	 = visco_num     * (zTe/Te_0)**(-1.5d0)
  if (visco_numm .gt. 1.d-10) visco_numm = 1.d-10 
  visco_par_numm = visco_par_num * (zTe/Te_0)**(-1.5d0)
  if (visco_par_numm .gt. 1.d-10) visco_par_numm = 1.d-10 
  D_perp_numm	 = D_perp_num    * (zTe/Te_0)**(-1.5d0)
  if (D_perp_numm .gt. 1.d-10) D_perp_numm = 1.d-10 
  Ki_perp_numm	 = ZK_perp_num   * (zTe/Te_0)**(-1.5d0)
  if (Ki_perp_numm .gt. 1.d-10) Ki_perp_numm = 1.d-10 
  Ke_perp_numm	 = ZK_perp_num   * (zTe/Te_0)**(-1.5d0)
  if (Ke_perp_numm .gt. 1.d-10) Ke_perp_numm = 1.d-10 
  Ki_par_num	 = 0.d-10
  Ke_par_num	 = 0.d-10		

  ! hyper-resistivity
  eta_numm       = eta_num       
  visco_numm     = visco_num     
  visco_par_numm = visco_par_num 
  D_perp_numm    = D_perp_num    
  Ki_perp_numm   = ZK_perp_num   
  Ke_perp_numm   = ZK_perp_num   
  Ki_par_num     = 0.d-10
  Ke_par_num     = 0.d-10               

  
  ! We need hyper diffusivities mostly at the grid axis
  !if (psi_norm .lt. 0.1d0) eta_numm    = eta_numm   * 1.d2
  !if (psi_norm .lt. 0.1d0) visco_numm  = visco_numm * 1.d2
  !if (psi_norm .lt. 0.05d0) visco_numm = visco_numm * 1.d4
  !if (psi_norm .lt. 0.05d0) eta_numm   = eta_numm   * 1.d4

  ! ------------------------------------------------
  ! --- Taylor Galerkin (TG2) stabilisation switches
  ! ------------------------------------------------
  TG_num1 = tgnum(1);
  TG_num2 = tgnum(2);
  TG_num5 = tgnum(5);
  TG_num6 = tgnum(6);
  TG_num7 = tgnum(7);
  TG_num8 = tgnum(8);
  
  
  ! ---------------------
  ! --- Bootstrap current
  ! ---------------------
  if (bootstrap) then
    ! --- Full Sauter formula
    call bootstrap_current(minRad, R, y_g,                       &
                           R_axis,   Z_axis,   psi_axis,         &
			   R_xpoint, Z_xpoint, psi_bnd, psi_norm,&
			   ps0, ps0_x, ps0_y,                    &
			   r0,  r0_x,  r0_y,                     &
			   Ti0, Ti0_x, Ti0_y,                    &
			   Te0, Te0_x, Te0_y,                  Jb)
    ! --- Full Sauter formula for initial profiles
    call density      (xpoint2,xcase2, y_g,Z_xpoint, ps0,psi_axis,psi_bnd, zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz)
    call temperature_i(xpoint2,xcase2, y_g,Z_xpoint, ps0,psi_axis,psi_bnd, zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2,dTi_dpsi2_dz)
    call temperature_e(xpoint2,xcase2, y_g,Z_xpoint, ps0,psi_axis,psi_bnd, zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2,dTe_dpsi2_dz)
    zTi_x = dTi_dpsi * ps0_x
    zTi_y = dTi_dpsi * ps0_y
    zTe_x = dTe_dpsi * ps0_x
    zTe_y = dTe_dpsi * ps0_y
    zn_x  = dn_dpsi * ps0_x
    zn_y  = dn_dpsi * ps0_y
    call bootstrap_current(minRad, R, y_g,                       &
                           R_axis,   Z_axis,   psi_axis,         &
			   R_xpoint, Z_xpoint, psi_bnd, psi_norm,&
			   ps0, ps0_x, ps0_y,                    &
			   zn,  zn_x,  zn_y,                     &
			   zTi, zTi_x, zTi_y,                    &
			   zTe, zTe_x, zTe_y,                  Jb_0)
    ! --- Subtract the initial equilibrium part
    Jb = Jb - Jb_0
  else
    Jb = 0.d0
  endif
  
  
  ! ------------------------------------------------------
  ! --- Diamagnetic terms, avoid problems at the target...
  ! ------------------------------------------------------
  tau_IC = tauIC
  if (Wdia) W_dia = 1.d0
  if (avoid_xpoint) then
    if (xpoint2 .and.  (xcase2 .ne. 2) ) then
      distance_xpoint = sqrt( (x_g - R_xpoint(1))**2 + (y_g - Z_xpoint(1))**2 )
      tau_IC = tau_IC * (0.5d0 - 0.5d0 * tanh( -(distance_xpoint - 0.05d0)/0.01d0 ) )
    endif
    if (xpoint2 .and.  (xcase2 .ne. 1) ) then
      distance_xpoint = sqrt( (x_g - R_xpoint(2))**2 + (y_g - Z_xpoint(2))**2 )
      tau_IC = tau_IC * (0.5d0 - 0.5d0 * tanh( -(distance_xpoint - 0.15d0)/0.01d0 ) )
    endif
  endif

  ! --- Switch off diamagnetic terms in private region?
  if ( (ps0 .le. psi_bnd) .and. (y_g .le. Z_xpoint(1)) ) tau_IC  = 0.d0
  !if ( (ps0 .le. psi_bnd) .and. (y_g .le. Z_xpoint(1)) ) visco_Te = visco_Te + 1.d3 * visco_Te
  !if ( ( (ps0-psi_axis)/(psi_bnd-psi_axis) .le. 1.01d0) .and. (y_g .le. Z_xpoint(1)) ) tau_IC  = 0.d0
  !if ( ( (ps0-psi_axis)/(psi_bnd-psi_axis) .le. 1.01d0) .and. (y_g .le. Z_xpoint(1)) ) visco_Te = visco_Te + 1.d3 * visco_Te

  ! --- Viscosity buffer (and other buffers) at targets
  if (R_limiter(1) .ne. 0.d0) then
    Rline1 = R_limiter(1)  ; Zline1 = Z_limiter(1)
    Rline2 = R_limiter(2)  ; Zline2 = Z_limiter(2)
    Rline3 = R_limiter(3)  ; Zline3 = Z_limiter(3)
    Rline4 = R_limiter(4)  ; Zline4 = Z_limiter(4)
  else
    Rline1 = 0.d0  ; Zline1 = -10.d0
    Rline2 = 1.d0  ; Zline2 = -11.d0
    Rline3 = 0.d0  ; Zline3 = -10.d0
    Rline4 = 1.d0  ; Zline4 = -11.d0
  endif
  if (R .lt. R_xpoint(1)) then
    slope  = (Zline1 - Zline2) / (Rline1 - Rline2)
    offset = (Rline1*Zline2 - Rline2*Zline1) / (Rline1 - Rline2)
  else
    slope  = (Zline3 - Zline4) / (Rline3 - Rline4)
    offset = (Rline3*Zline4 - Rline4*Zline3) / (Rline3 - Rline4)
  endif
  Z_tmp = offset + slope * R
  target_buffer_width = 0.02d0
  tan_width           = target_buffer_width / 1.d0
  ! --- Choose buffers
  visco_Te   = visco_Te + 1.d3 * visco_Te * (0.5d0 - 0.5d0 * tanh(  (abs(y_g - Z_tmp) - target_buffer_width)/tan_width ))
  !eta_Te     = eta_Te   + 1.d3 * eta_Te   * (0.5d0 - 0.5d0 * tanh(  (abs(y_g - Z_tmp) - target_buffer_width)/tan_width ))
  !dvisco_dTe = dvisco_dTe                 * (0.5d0 - 0.5d0 * tanh( -(abs(y_g - Z_tmp) - target_buffer_width)/tan_width ))
  !deta_dTe   = deta_dTe                   * (0.5d0 - 0.5d0 * tanh( -(abs(y_g - Z_tmp) - target_buffer_width)/tan_width ))
  tau_IC     = tau_IC                     * (0.5d0 - 0.5d0 * tanh( -(abs(y_g - Z_tmp) - target_buffer_width)/tan_width ))

  ! --- Viscosity buffer for the type-2 boundary
  visco_Te   = visco_Te + 1.d3 * visco_Te * (0.5d0 - 0.5d0 * tanh(  -( psi_norm - (rho_coef(5)+4*rho_coef(4)) )/rho_coef(4) ))
  !tau_IC     = tau_IC                     * (0.5d0 - 0.5d0 * tanh(   ( psi_norm - (rho_coef(5)+4*rho_coef(4)) )/rho_coef(4) ))

  ! -------------------------
  ! --- Neoclassical rotation
  ! -------------------------
  epsil   = 1.d-3
  Btheta2 = (ps0_x**2.d0 + ps0_y**2.d0) / R**2
  if ( NEO ) then 
    if (num_neo_file) then
      call neo_coef(xpoint2, xcase2, y_g, Z_xpoint, ps0, psi_axis, psi_bnd, amu_neo_prof, aki_neo_prof)
    else
       amu_neo_prof = amu_neo_const
       aki_neo_prof = aki_neo_const
    endif
  else 
    amu_neo_prof   = 0.d0
    aki_neo_prof   = 0.d0
  endif
  

  ! -------------------------------------------------------------------
  ! --- Heating, current and particle source (the same for all i_plane)
  ! -------------------------------------------------------------------
  if (i_plane .eq. 1) then
    ! --- Current source
    call current(xpoint2, xcase2, x_g,y_g, Z_xpoint, ps0,psi_axis,psi_bnd,current_source)

    ! --- Density source and heating
    call sources(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd,particle_source,heat_source_i,heat_source_e)
    ! --- Old simple uniform sources
    !particle_source   = particlesource * (0.5d0 - 0.5d0 * tanh((psi_norm-0.5)/0.005) )
    !heat_source_i     = heatsource_i   * (0.5d0 - 0.5d0 * tanh((psi_norm-0.5)/0.005) ) 
    !heat_source_e     = heatsource_e   * (0.5d0 - 0.5d0 * tanh((psi_norm-0.5)/0.005) ) 
    ! --- New source profile: depends on initial equilibirum profiles.
    call density      (xpoint2,xcase2, y_g,Z_xpoint, ps0,psi_axis,psi_bnd, zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz)
    call temperature_i(xpoint2,xcase2, y_g,Z_xpoint, ps0,psi_axis,psi_bnd, zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2,dTi_dpsi2_dz)
    call temperature_e(xpoint2,xcase2, y_g,Z_xpoint, ps0,psi_axis,psi_bnd, zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2,dTe_dpsi2_dz)
    particle_source = particlesource * ( zn  - r0 )
    heat_source_i   = heatsource_i   * ( zTi - Ti0 )
    heat_source_e   = heatsource_e   * ( zTe - Te0 )
    particle_source = particle_source * ( 0.5d0 - 0.5d0 * tanh((psi_norm-0.99)/0.005) )
    heat_source_i   = heat_source_i   * ( 0.5d0 - 0.5d0 * tanh((psi_norm-0.99)/0.005) )
    heat_source_e   = heat_source_e   * ( 0.5d0 - 0.5d0 * tanh((psi_norm-0.99)/0.005) )
    if( (xpoint2) .and. (xcase2 .ne. 2) .and. (y_g .le. Z_xpoint(1)) ) particle_source = 0.d0
    if( (xpoint2) .and. (xcase2 .ne. 2) .and. (y_g .le. Z_xpoint(1)) ) heat_source_i   = 0.d0
    if( (xpoint2) .and. (xcase2 .ne. 2) .and. (y_g .le. Z_xpoint(1)) ) heat_source_e   = 0.d0
    if( (xpoint2) .and. (xcase2 .ne. 1) .and. (y_g .ge. Z_xpoint(2)) ) particle_source = 0.d0
    if( (xpoint2) .and. (xcase2 .ne. 1) .and. (y_g .ge. Z_xpoint(2)) ) heat_source_i   = 0.d0
    if( (xpoint2) .and. (xcase2 .ne. 1) .and. (y_g .ge. Z_xpoint(2)) ) heat_source_e   = 0.d0
    if(particle_source .le. 0.d0) particle_source = 0.d0
    if(heat_source_i   .le. 0.d0) heat_source_i   = 0.d0
    if(heat_source_e   .le. 0.d0) heat_source_e   = 0.d0

    ! --- Toroidal momentum source
    if ( abs(V_0) .ge. 1.e-12) then 
      call velocity(xpoint2,xcase2, y_g,Z_xpoint, ps0,psi_axis,psi_bnd, V_source,dV_dpsi_source,dV_dz_source,dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2,dV_dpsi2_dz)
      Vt0_x = dV_dpsi_source * ps0_x
      Vt0_y = dV_dz_source + dV_dpsi_source * ps0_y
    endif
  endif
  
  ! --- Avoid negative density/temperature
  if ( r0 .lt. 1.d-1*rho_1 ) then
    !D_prof  = D_prof  * 1.d2
    !particle_source = 3.d-3
  endif       
  if ( Ti0 .lt. 1.d-0*Ti_1 ) then
    !Ki_prof = Ki_prof * 1.d2
    !Ki_par  = Ki_par  * 1.d2
    !dKi_par = 0.d0
    !heat_source_i = 1.d-5
  endif
  if ( Te0 .lt. 1.d-0*Te_1 ) then
    !Ke_prof = Ke_prof * 1.d2
    !Ke_par  = Ke_par  * 1.d2
    !dKe_par = 0.d0
    !heat_source_e = 1.d-5
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
  use mod_parameters    
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


