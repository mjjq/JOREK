!> This routine finds the intersections of a flux surface with the target (ie. wall below Xpoint)
subroutine get_target_flux_surfaces(node_list, element_list, surface_list,       &
                                    psi_bnd, R_axis, Z_axis, R_xpoint, Z_xpoint, &
				    n_int_max, n_int, R_int, Z_int, index_int, ifail)


  use data_structure
  use phys_module
  use high_resolution_wall
  use constants
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),    intent(in)		:: node_list
  type (type_element_list), intent(in)		:: element_list
  type (type_surface_list), intent(inout)	:: surface_list
  integer,                  intent(in)		:: n_int_max
  integer,                  intent(inout)	:: n_int, ifail
  real*8,                   intent(inout)	:: R_int(n_int_max), Z_int(n_int_max)
  real*8,                   intent(inout)	:: R_axis,           Z_axis
  real*8,                   intent(inout)	:: R_xpoint(2),      Z_xpoint(2)
  real*8,                   intent(inout)	:: psi_bnd
  integer,                  intent(inout)	:: index_int(n_int_max,3) ! 1->surface_index, 2->surface_piece_index, 3->wall_piece_index
  
  ! --- Local variables
  type (type_surface_list)	:: surfaces_tmp
  type (type_surface)		:: surface
  integer			:: i, i_surf
  real*8			:: theta_x(2), Z_line(2)
  real*8			:: s_find(8),t_find(8)
  integer			:: i_elm_find(8),i_find
  real*8			:: R,R2,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt
  real*8			:: Z,Z2,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt
  integer			:: n_int_tmp, index_int_tmp(n_int_max,3)
  real*8			:: R_int_tmp(n_int_max), Z_int_tmp(n_int_max)
  logical			:: debug
  character*256			:: filename
  
  write(*,*) '***********************************'
  write(*,*) '*     get_target_flux_surfaces    *'
  write(*,*) '***********************************'
  
  n_int = 0
  debug = .true.
  
  ! --- First get all intersections with wall
  do i_surf = 1,surface_list%n_psi
    call find_wall_crossings_with_flux_surface(node_list, element_list, surface_list%flux_surfaces(i_surf), &
                                               n_int_max, n_int_tmp, R_int_tmp, Z_int_tmp, index_int_tmp, ifail)

    if (ifail .ne. 0) then
      write(*,*) 'Warning! Failed to find all wall intersections for surface',i_surf,ifail
      return
    endif
    
    do i=1,n_int_tmp
      R_int(n_int + i) = R_int_tmp(i)
      Z_int(n_int + i) = Z_int_tmp(i)
      index_int(n_int + i, 1) = i_surf
      index_int(n_int + i, 2) = index_int_tmp(i,2)
      index_int(n_int + i, 3) = index_int_tmp(i,3)
    enddo
    n_int = n_int + n_int_tmp
  enddo
  
  ! --- Some debug plots
  if (debug) then
    write(*,*)'number of total intersections found:',n_int 
    filename = 'plot_wall_intersections.py'
    call print_py_plot_prepare_plot(filename)
    call print_py_plot_ordered_flux_surfaces(filename, node_list, element_list, surface_list)
    call print_py_plot_points(filename,n_int,R_int,Z_int)
    call print_py_plot_wall(filename)
    call print_py_plot_finish_plot(filename)
  endif

  ! --- Copy back into temporary array
  do i=1,n_int
      R_int_tmp(i) = R_int(i)
      Z_int_tmp(i) = Z_int(i)
      index_int_tmp(i,1) = index_int(i,1)
      index_int_tmp(i,2) = index_int(i,2)
      index_int_tmp(i,3) = index_int(i,3)
  enddo
  n_int_tmp = n_int
  
  ! --- Determine angle of X-line
  if (xcase .ne. 2) theta_x(1) = atan2(Z_xpoint(1)-Z_axis,R_xpoint(1)-R_axis) + 0.5d0*PI
  if (xcase .ne. 1) theta_x(2) = atan2(Z_xpoint(2)-Z_axis,R_xpoint(2)-R_axis) + 1.5d0*PI
  if (theta_x(1) .lt. 0.d0)    theta_x(1) = theta_x(1) + 2.d0*PI
  if (theta_x(2) .lt. 0.d0)    theta_x(2) = theta_x(2) + 2.d0*PI
  if (theta_x(1) .gt. 2.d0*PI) theta_x(1) = theta_x(1) - 2.d0*PI
  if (theta_x(2) .gt. 2.d0*PI) theta_x(2) = theta_x(2) - 2.d0*PI
  
  ! --- Now extract target points. First, retain only points below Xpoint-line (or above 2nd Xpoint's)
  n_int = 0
  do i=1,n_int_tmp
    if (xcase .ne. 2) Z_line(1) = Z_xpoint(1) + (R_int_tmp(i) - R_xpoint(1)) * tan(theta_x(1))
    if (xcase .ne. 1) Z_line(2) = Z_xpoint(2) + (R_int_tmp(i) - R_xpoint(2)) * tan(theta_x(2))
    if (     ((xcase .ne. 2) .and. (Z_int_tmp(i) .lt. Z_line(1))) &
        .or. ((xcase .ne. 1) .and. (Z_int_tmp(i) .gt. Z_line(2))) ) then
      n_int = n_int + 1
      R_int(n_int) = R_int_tmp(i)
      Z_int(n_int) = Z_int_tmp(i)
      index_int(n_int,1) = index_int_tmp(i,1)
      index_int(n_int,2) = index_int_tmp(i,2)
      index_int(n_int,3) = index_int_tmp(i,3)
    endif
  enddo
  
  ! --- Some debug plots
  if (debug) then
    write(*,*)'number of total target intersections found:',n_int 
    filename = 'plot_target_intersections.py'
    call print_py_plot_prepare_plot(filename)
    call print_py_plot_ordered_flux_surfaces(filename, node_list, element_list, surface_list)
    call print_py_plot_points(filename,n_int,R_int,Z_int)
    call print_py_plot_wall(filename)
    call print_py_plot_finish_plot(filename)
  endif

  
  ! --- Initialise some data
  !allocate(surfaces_tmp%psi_values(1))
  !surfaces_tmp%psi_values(1) = surface%psi
  
  ! --- If we found just one intersection, simple
  !if (n_int_tmp .eq. 1) then
  !  ! --- Get parameters of wall line
  !  R = R_limiter(index_wall(1))
  !  Z = Z_limiter(index_wall(1))
  !  theta = atan2( Z_limiter(index_wall(1)+1) - Z, R_limiter(index_wall(1)+1) - R )

  !  surfaces_tmp%flux_surfaces(1)%n_pieces = 1
  !  surfaces_tmp%flux_surfaces(1)%s(:,1) = surface%s(:,index_surf(1))
  !  surfaces_tmp%flux_surfaces(1)%t(:,1) = surface%t(:,index_surf(1))
  !  surfaces_tmp%flux_surfaces(1)%elm(1) = surface%elm(index_surf(1))
  !  call find_theta_surface(node_list,element_list,surfaces_tmp,1,theta,R,Z,i_elm_find,s_find,t_find,i_find)
  !  if (i_find .ne. 1) then
  !    write(*,*)'Warning! Did not find exactly one intersection with the wall:',i_find
  !    return
  !  endif
  !  n_int = n_int + 1
  !  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1), &
  !  		   R,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt, &
  !  		   Z,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt)
  !  R_int(n_int) = R
  !  Z_int(n_int) = Z
  !endif
  
  ! --- If we found more than one intersection, need to make sure they all belong to the target
  !if (n_int .eq. 2) then

  !endif
  
  !deallocate(surfaces_tmp%psi_values)


end subroutine get_target_flux_surfaces







!> This routine finds all intersections between a flux surface and the wall
subroutine find_wall_crossings_with_flux_surface(node_list, element_list, surface, n_int_max, n_int, R_int, Z_int, index_int, ifail)

  use data_structure
  use phys_module
  use high_resolution_wall
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),    intent(in)		:: node_list
  type (type_element_list), intent(in)		:: element_list
  type (type_surface),      intent(inout)	:: surface
  integer,                  intent(in)		:: n_int_max
  integer,                  intent(inout)	:: n_int, ifail
  real*8,                   intent(inout)	:: R_int(n_int_max), Z_int(n_int_max)
  integer,                  intent(inout)	:: index_int(n_int_max,3)
  
  ! --- Local variables
  integer			:: i,j,i_elm,n_tmp,k,l
  integer			:: node1, node2, node3, node4
  integer			:: inside, inside_save
  real*8			:: dl, si, ti, dsi, dti, length, length_min
  real*8			:: ss1, dss1, ss_int, dl_int
  real*8			:: ss2, dss2
  real*8			:: tt1, dtt1, tt_int
  real*8			:: tt2, dtt2
  real*8			:: R,R2,R_save,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt
  real*8			:: Z,Z2,Z_save,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt
  real*8			:: R_cub1d(4)
  real*8			:: Z_cub1d(4)
  real*8			:: surface_accuracy
  
  ! --- Initialise some data
  surface_accuracy = 1.d-3
  n_int     = 0
  
  ! --- Step along surface with high resolution and check that points are inside wall
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
	  ! --- Record the surface piece of the intersection
	  n_int = n_int + 1
	  index_int(n_int,2) = j
          ! --- Record (find) the wall piece that we've just crossed
	  do l=1,n_limiter-1
            ! --- Find the intersection itself
            R_cub1d(1) = R_limiter(l    )               ;   	Z_cub1d(1) = Z_limiter(l    )
            R_cub1d(3) = R_limiter(l + 1)               ;   	Z_cub1d(3) = Z_limiter(l + 1)
            R_cub1d(2) = (R_cub1d(3)-R_cub1d(1))/2.d0	;   	Z_cub1d(2) = (Z_cub1d(3)-Z_cub1d(1))/2.d0
            R_cub1d(4) = (R_cub1d(3)-R_cub1d(1))/2.d0	;   	Z_cub1d(4) = (Z_cub1d(3)-Z_cub1d(1))/2.d0
            call find_crossing_on_surface_piece(node_list,element_list,surface,j,R_cub1d,Z_cub1d, &
        	                                R_int(n_int),Z_int(n_int),ss_int,tt_int,dl_int,ifail)
            if (ifail .eq. 0) then
	      index_int(n_int,3) = l
	      exit
	    endif
            if ( (ifail .ne. 0) .and. (l .eq. n_limiter-1) ) then
	      write(*,*)'Warning! Failed to find intersection with wall at',index_int(n_int,2),index_int(n_int,3)
	      return
	    endif
	  enddo
	endif
	
	! --- Save previous point
	inside_save = inside
        R_save = R
	Z_save = Z
	
      enddo
        
    enddo
  enddo
  
  return

end subroutine find_wall_crossings_with_flux_surface
  
  




