!> This routine reorders fluxsurfaces so that pieces are one after the other
subroutine find_wall_crossing_with_flux_surface(node_list, element_list, surface, psi_bnd, Z_xpoint, n_int, R_int, Z_int, index_wall)

  use data_structure
  use phys_module
  use high_resolution_wall
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),    intent(in)		:: node_list
  type (type_element_list), intent(in)		:: element_list
  type (type_surface),      intent(inout)	:: surface
  integer,                  intent(inout)	:: n_int
  real*8,                   intent(inout)	:: R_int(20), Z_int(20)
  real*8,                   intent(inout)	:: psi_bnd,   Z_xpoint(2)
  integer,                  intent(inout)	:: index_wall(20)
  
  ! --- Local variables
  integer			:: i
  integer			:: node1, node2, node3, node4
  integer			:: inside, inside_save
  real*8			:: dl, si, ti
  real*8			:: ss1, dss1
  real*8			:: ss2, dss2
  real*8			:: tt1, dtt1
  real*8			:: tt2, dtt2
  real*8			:: R,R2,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt
  real*8			:: Z,Z2,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt
  real*8			:: surface_accuracy
  integer			:: n_int_tmp, piece(20), wall_piece(20)
  type (type_surface_list)	:: surfaces_tmp
  real*8			:: theta
  real*8			:: s_find(8),t_find(8)
  integer			:: i_elm_find(8),i_find
  
  ! --- Initialise some data
  allocate(surfaces_tmp%psi_values(1))
  surfaces_tmp%psi_values(1) = surface%psi
  surface_accuracy = 1.d-3
  n_int     = 0
  n_int_tmp = 0
  
  ! --- First step along surface with high resolution and check that points are inside wall
  do i=1,surface%n_parts
    do j=surface%parts_index(i),surface%parts_index(i+1)
      i_elm = surface%elm(j)
      
      ! --- Get approximate length of piece
      ss1 = surface%s(1,j)
      tt1 = surface%t(1,j)
      call interp_RZ(node_list,element_list,i_elm,ss1,tt1,    &
                     R,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt, &
    	  	     Z,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt)

      ss2 = surface%s(3,j)
      tt2 = surface%t(3,j)
      call interp_RZ(node_list,element_list,i_elm,ss2,tt2,     &
                     R2,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt, &
      	  	     Z2,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt)
      
      length = sqrt( (R-R2)**2.d0 + (Z-Z2)**2.d0 )
      
      ! --- Deduce number of points on piece required to get accuracy
      n_tmp = length / surface_accuracy
      
      ! --- Get variables for piece extrapolation
      node1 = element_list%element(i_elm)%vertex(1)
      node2 = element_list%element(i_elm)%vertex(2)
      node3 = element_list%element(i_elm)%vertex(3)
      node4 = element_list%element(i_elm)%vertex(4)

      ss1  = surface%s(1,j)
      dss1 = surface%s(2,j)
      ss2  = surface%s(3,j)
      dss2 = surface%s(4,j)

      tt1  = surface%t(1,j)
      dtt1 = surface%t(2,j)
      tt2  = surface%t(3,j)
      dtt2 = surface%t(4,j)

      ! --- Step on each piece point
      do k=1,n_tmp

        dl = -1.d0 + 2.d0 * real(k-1) / real(n_tmp-1)
        call CUB1D(ss1, dss1, ss2, dss2, dl, si, dsi)
        call CUB1D(tt1, dtt1, tt2, dtt2, dl, ti, dti)
        call interp_RZ(node_list,element_list,i_elm,si,ti,      &
	               R,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt, &
      		       Z,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt)
        
	! --- Check if point is inside/outside wall
	call check_point_is_inside_wall(R, Z, inside)
	if (k .eq. 1) inside_save = inside
	
	! --- Have we stepped accross the wall?
	if (inside .ne. inside_save) then
	  if (     ((Z .lt. Z_xpoint(1)) .and. (xcase .ne. 2))
	      .or. ((Z .gt. Z_xpoint(2)) .and. (xcase .ne. 1)) )
	    n_int_tmp = n_int_tmp + 1
	    piece(n_int_tmp) = j
            ! --- Find the wall piece that we've just crossed
	    length_min = 1.d10
	    do l=1,n_limiter-1
	      length = 0.5d0 * (  sqrt( (R-R_limiter(l  ))**2.d0 + (Z-Z_limiter(l  ))**2.d0 ) &
	                        + sqrt( (R-R_limiter(l+1))**2.d0 + (Z-Z_limiter(l+1))**2.d0 ) )
	      if (length .lt. length_min) then
	        length_min = length
		wall_piece(n_int_tmp) = l
	      endif
	    enddo
	    exit
	  endif
	endif
	
      enddo
        
    enddo
  enddo
  
  ! --- If we found just one intersection, simple
  if (n_int .eq. 1) then
    ! --- Get parameters of wall line
    R = R_limiter(wall_piece(1))
    Z = Z_limiter(wall_piece(1))
    theta = atan2( Z_limiter(wall_piece(1)+1) - Z, R_limiter(wall_piece(1)+1) - R )

    surfaces_tmp%flux_surfaces(1)%n_pieces = 1
    surfaces_tmp%flux_surfaces(1)%s(:,1) = surface%s(:,piece(1))
    surfaces_tmp%flux_surfaces(1)%t(:,1) = surface%t(:,piece(1))
    surfaces_tmp%flux_surfaces(1)%elm(1) = surface%elm(piece(1))
    call find_theta_surface(node_list,element_list,surfaces_tmp,1,theta,R,Z,i_elm_find,s_find,t_find,i_find)
    if (i_find .ne. 1) then
      write(*,*)'Warning! Did not find exactly one intersection with the wall:',i_find
      return
    endif
    n_int = n_int + 1
    call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1), &
    		   R,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt, &
    		   Z,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt)
    R_int(n_int) = R
    Z_int(n_int) = Z
    index_wall(n_int) = wall_piece(1)
  endif
  
  ! --- If we found more than one intersection, need to make sure they all belong to the target
  if (n_int .eq. 2) then

  endif
  
  deallocate(surfaces_tmp%psi_values)
  
  return

end subroutine find_wall_crossing_with_flux_surface
  
  




