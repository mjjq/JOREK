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
  use parameters    
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
	T0	     = T0	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_s	     = T0_s	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_t	     = T0_t	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_ss	     = T0_ss	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_tt	     = T0_tt	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_st	     = T0_st	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_p	     = T0_p	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	T0_pp	     = T0_pp	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)
  
	! --- Variable 7
	Vpar0	     = Vpar0	  + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_s      = Vpar0_s    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_t      = Vpar0_t    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_ss     = Vpar0_ss   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_tt     = Vpar0_tt   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_st     = Vpar0_st   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_p      = Vpar0_p    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	Vpar0_pp     = Vpar0_pp	  + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)

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
  r0_corr = corr_neg_dens(r0)
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
  T0_corr    = corr_neg_temp(T0) ! For use in eta(T), visco(T), ...
  T0_x      = (   y_t * T0_s  - y_s * T0_t ) / xjac
  T0_y      = ( - x_t * T0_s  + x_s * T0_t ) / xjac
  T0_xx     = (T0_ss * y_t**2 - 2.d0*T0_st * y_s*y_t + T0_tt * y_s**2	& 	    
	      + T0_s * (y_st*y_t - y_tt*y_s )			        &    
	      + T0_t * (y_st*y_s - y_ss*y_t ) )        / xjac**2        & 	
	     - xjac_x * (T0_s * y_t - T0_t * y_s)      / xjac**2
  T0_yy     = (T0_ss * x_t**2 - 2.d0*T0_st * x_s*x_t + T0_tt * x_s**2   & 	    
	      + T0_s * (x_st*x_t - x_tt*x_s )			        &    
	      + T0_t * (x_st*x_s - x_ss*x_t ) )        / xjac**2        & 	
	     - xjac_y * (- T0_s * x_t + T0_t * x_s )   / xjac**2
  T0_xy     = (- T0_ss * y_t*x_t - T0_tt * x_s*y_s		        &
	       + T0_st * (y_s*x_t  + y_t*x_s  )  		        &        
	       - T0_s  * (x_st*y_t - x_tt*y_s )  		        &    
	       - T0_t  * (x_st*y_s - x_ss*y_t )  )     / xjac**2        & 	
	       - xjac_x * (- T0_s * x_t + T0_t * x_s ) / xjac**2

  
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
  
  ! --- Deltas
  delta_u_x  = (   y_t * delta_s(2) - y_s * delta_t(2) ) / xjac
  delta_u_y  = ( - x_t * delta_s(2) + x_s * delta_t(2) ) / xjac
  delta_ps_x = (   y_t * delta_s(1) - y_s * delta_t(1) ) / xjac
  delta_ps_y = ( - x_t * delta_s(1) + x_s * delta_t(1) ) / xjac
  
  ! --- Pressure
  P0	   = r0    * T0
  P0_x     = r0_x  * T0 + r0 * T0_x
  P0_y     = r0_y  * T0 + r0 * T0_y
  P0_s     = r0_s  * T0 + r0 * T0_s
  P0_t     = r0_t  * T0 + r0 * T0_t
  P0_p     = r0_p  * T0 + r0 * T0_p
  P0_pp    = r0_pp * T0 + r0 * T0_pp + 2.d0 * r0_p * T0_p
  P0_xx    = r0_xx * T0 + r0 * T0_xx + 2.d0 * r0_x * T0_x
  P0_yy    = r0_yy * T0 + r0 * T0_yy + 2.d0 * r0_y * T0_y
  P0_xy    = r0_xy * T0 + r0 * T0_xy + r0_x * T0_y + r0_y * T0_x
  
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
  use parameters    
  use basis_at_gaussian
  use phys_module
  use equation_variables
  use data_structure
  use diffusivities, only: get_dperp, get_zkperp
  use pellet_module
  use bootstrap_functions
  
  implicit none
  
  ! --- Routine variables
  type (type_element)	      :: element
  type (type_node)	      :: nodes(n_vertex_max)
  logical		      :: xpoint2
  integer		      :: xcase2
  integer		      :: i_plane
  real*8		      :: minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
  
  ! --- Internal variables
  real*8		      :: psi_norm
  real*8		      :: V_source, dV_dpsi2, dV_dz2, dV_dpsi_dz, dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz
  real*8		      :: Ti0, Ti0_x, Ti0_y, Te0, Te0_x, Te0_y
  real*8		      :: zTi, zTi_x, zTi_y, zTe, zTe_x, zTe_y, zn_x, zn_y
  real*8		      :: Jb_0
      
  
  ! -------------------------------------
  ! --- Temperature dependent resistivity
  ! -------------------------------------
  if ( eta_T_dependent ) then
    eta_T     =   eta   * (T0_corr/T_0)**(-1.5d0)
    deta_dT   = - eta	* (1.5d0)  * T0_corr**(-2.5d0) * T_0**(1.5d0)
    d2eta_d2T =   eta	* (3.75d0) * T0_corr**(-3.5d0) * T_0**(1.5d0)
  else
    eta_T     = eta
    deta_dT   = 0.d0
    d2eta_d2T = 0.d0
  end if
  
  
  ! -----------------------------------
  ! --- Temperature dependent viscosity
  ! -----------------------------------
  if ( visco_T_dependent ) then       
    visco_T   =   visco * (T0_corr/T_0)**(-1.5d0)
    dvisco_dT = - visco * (1.5d0)  * T0_corr**(-2.5d0) * T_0**(1.5d0)
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

  ! --- Call Diff functions
  D_prof = get_dperp (psi_norm)
  K_prof = get_zkperp(psi_norm)
  if (xpoint2) then
    if (r0 .lt. 0.d0)  then
      D_prof = D_prof_neg  ! JET : 1.d-4; ITER :  4.d-3
    endif
    if (T0 .lt. 0.d0) then
      K_prof = ZK_prof_neg  ! JET : 1.d-3; ITER : 2.d-2 
    endif
  endif

  
  
  ! -----------------------------------------------------
  ! --- Parallel conductivity profiles (Braginskii model)
  ! -----------------------------------------------------
  if ( ZKpar_T_dependent ) then
    K_par    = ZK_par * (T0_corr/T_0)**(+2.5d0)
    dK_par   = ZK_par * (2.5d0)  * T0_corr**(+1.5d0) * T_0**(-2.5d0)
    if (K_par .gt. ZK_par_max) then
      K_par  = Zk_par_max
      dK_par = 0.d0
    endif
  else
    K_par  = ZK_par
    dK_par = 0.d0
  endif
  
 
  ! -------------------------
  ! --- Hyper diffusivitities
  ! -------------------------
  eta_numm	 = eta_num		! hyper-resistivity
  visco_numm	 = visco_num		! hyper-viscosity
  visco_par_numm = visco_par_num	! hyper-viscosity
  D_perp_numm	 = D_perp_num		! hyper-diffusivity
  K_perp_numm	 = ZK_perp_num		! hyper-conductivity

  
  ! ------------------------------------------------
  ! --- Taylor Galerkin (TG2) stabilisation switches
  ! ------------------------------------------------
  TG_num1 = tgnum(1);
  TG_num2 = tgnum(2);
  TG_num5 = tgnum(5);
  TG_num6 = tgnum(6);
  TG_num7 = tgnum(7);
  
  
  ! ---------------------
  ! --- Bootstrap current
  ! ---------------------
  if (bootstrap) then
    ! --- Full Sauter formula
    Ti0   = T0   / 2.d0 ; Te0	= T0   / 2.d0
    Ti0_x = T0_x / 2.d0 ; Te0_x = T0_x / 2.d0
    Ti0_y = T0_y / 2.d0 ; Te0_y = T0_y / 2.d0
    call bootstrap_current(minRad, R, y_g,                       &
                           R_axis,   Z_axis,   psi_axis,         &
			   R_xpoint, Z_xpoint, psi_bnd, psi_norm,&
			   ps0, ps0_x, ps0_y,                    &
			   r0,  r0_x,  r0_y,                     &
			   Ti0, Ti0_x, Ti0_y,                    &
			   Te0, Te0_x, Te0_y,                  Jb)
    ! --- Full Sauter formula for initial profiles
    call density(    xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
    		     zn,dn_dpsi,  dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
    call temperature(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
    		     zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)
    zTi   = zT / 2.d0             
    zTi_x = dT_dpsi * ps0_x / 2.d0
    zTi_y = dT_dpsi * ps0_y / 2.d0
    zTe   = zTi  
    zTe_x = zTi_x
    zTe_y = zTi_y
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
  
  
  ! -------------------------
  ! --- Neoclassical rotation
  ! -------------------------
  epsil   = 1.d-3
  Btheta2 = (ps0_x**2 + ps0_y**2) / R**2
    amu_neo_prof   = 0.d0
    aki_neo_prof   = 0.d0
  if ( NEO ) then 
    if (num_neo_file) then
      call neo_coef(xpoint2, xcase2, y_g, Z_xpoint, ps0, psi_axis, psi_bnd, amu_neo_prof,          &
        aki_neo_prof)
    else
       amu_neo_prof = amu_neo_const
       aki_neo_prof = aki_neo_const
    endif
  endif
  
  ! -------------------------------------------------------------------
  ! --- Heating, current and particle source (the same for all i_plane)
  ! -------------------------------------------------------------------
  if (i_plane .eq. 1) then
    ! --- Current source
    call current(xpoint2, xcase2, x_g,y_g, Z_xpoint, ps0,psi_axis,psi_bnd,current_source)

    ! --- Density and Temperature source
    call sources(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd,particle_source,heat_source)
    
    ! --- New source profile: source with exactly the same profile as the initial equilibirum profiles.
    call density(    xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
    		     zn,dn_dpsi,  dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
    call temperature(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
    		     zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)
    
    ! --- Toroidal momentum source (NBI)
  dV_dpsi_source = 0.d0
  dV_dz_source   = 0.d0
  if ( ( abs(V_0) .ge. 1.e-12 ) .or. ( num_rot ) ) then
    call velocity(xpoint2, xcase2, y_g, z_xpoint, ps0, psi_axis, psi_bnd, V_source,               &
      dV_dpsi_source, dV_dz_source, dV_dpsi2, dV_dz2, dV_dpsi_dz, dV_dpsi3,dV_dpsi_dz2,           &
      dV_dpsi2_dz)
    if (normalized_velocity_profile) then
      Vt0_x = dV_dpsi_source * ps0_x
      Vt0_y = dV_dz_source + dV_dpsi_source * ps0_y
    else
      Omega_tor0_x = dV_dpsi_source * ps0_x
      Omega_tor0_y = dV_dz_source + dV_dpsi_source * ps0_y
    end if
  end if
    
    ! --- Pellet Source
    source_pellet = 0.d0
    source_volume = 0.d0
    if (use_pellet) then
      call pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
    			  pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta, &
    			  x_g,y_g, ps0, phi, zn, zT, &
    			  central_density, pellet_particles, pellet_density, total_pellet_volume, &
    			  source_pellet, source_volume)
    endif
    
    ! --- Total density source 
    total_rho_source = particle_source + source_pellet
  endif
  
  return

end subroutine ELM_build_diffusivities_and_sources










!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!-------------------------------- Compute the basis functions (ie. the Bezier polynomials) ------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
subroutine ELM_build_basis_functions(element, nodes, ms, mt, i_plane, i_vertex, i_order, i_tor, &
				     vv, vv_s,  vv_t,	     vv_p,  vv_x,  vv_y,	       &
        				 vv_ss, vv_tt, vv_st, vv_pp, vv_xx, vv_yy, vv_xy        )
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_basis_functions

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

end subroutine ELM_build_basis_functions


