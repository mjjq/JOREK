module reorder_surfaces_parameters
  real*8, parameter	:: accuracy = 1.d-7
  integer, parameter	:: n_parts_max = 10
end module reorder_surfaces_parameters


!> This routine reorders fluxsurfaces so that pieces are one after the other
subroutine reorder_flux_surfaces(node_list, element_list, surface_list, ier)

  use data_structure
  use reorder_surfaces_parameters
  
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),    intent(in)		:: node_list
  type (type_element_list), intent(in)		:: element_list
  type (type_surface_list), intent(inout)	:: surface_list
  integer,                  intent(inout)	:: ier
  
  ! --- Local variables
  integer	:: i_surf, i, j, nStart
  integer	:: i_piece, i_piece2, i_piece3
  integer	:: found,   found2,   found3
  integer	:: index1, index2, index_save
  integer	:: n_parts, parts_index(n_parts_max)
  integer	:: n_edge_pieces,     index_edge_pieces(2*n_parts_max)
  integer	:: n_isolated_pieces, index_isolated_pieces(2*n_parts_max)
  logical	:: invert, debug, finished
  integer	:: i_elm
  integer	:: i_elm2
  real*8	:: rr,    ss
  real*8	:: rr2,   ss2
  real*8	:: R, dRR_dr, dRR_ds, dRR_drs, dRR_drr, dRR_dss
  real*8	:: R2,dRR2_dr,dRR2_ds,dRR2_drs,dRR2_drr,dRR2_dss
  real*8	:: Z, dZZ_dr, dZZ_ds, dZZ_drs, dZZ_drr, dZZ_dss
  real*8	:: Z2,dZZ2_dr,dZZ2_ds,dZZ2_drs,dZZ2_drr,dZZ2_dss
  real*8	:: distance
  
  write(*,*) '***********************************'
  write(*,*) '*     reorder_flux_surfaces       *'
  write(*,*) '***********************************'
  

  debug = .true.
  ier   = 0
  
  ! --- Get a plot?
  if (debug) call print_py_plot(node_list, element_list, surface_list)
  
  ! --- Loop over all surfaces
  do i_surf = 1, surface_list%n_psi
    
    ! --- Set zero parts at first
    n_parts = 0
    
    ! --- Find all surface pieces that have no neighbour (ie. end pieces)
    call find_all_edge_pieces(node_list, element_list, surface_list%flux_surfaces(i_surf), &
                              n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)
    
    ! --- If there are isolated pieces, put them at the begining
    do i = 1, n_isolated_pieces
      invert = .false.
      index1 = i
      index2 = index_isolated_pieces(i)
      call swap_surface_pieces(surface_list%flux_surfaces(i_surf), index1, index2, invert, &
                               n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)
      n_parts	     = i
      parts_index(i) = i
    enddo
    
    ! ---------------------------------------------------------
    ! --- In case there are no edge pieces (ie. closed surface)
    ! ---------------------------------------------------------
    if (n_edge_pieces .eq. 0) then
      
      ! --- Start at first piece (check if there are isolated pieces)
      nStart = 1
      if (n_isolated_pieces .gt. 0) nStart = parts_index(n_parts) + 1
      n_parts = n_parts + 1
      parts_index(n_parts) = nStart
      
      ! --- There might be several closed surfaces
      finished = .false.
      do while(.not. finished)
        ! --- Loop over all remaining pieces
	do i_piece = nStart, surface_list%flux_surfaces(i_surf)%n_pieces - 1
	  if (i_piece .eq. surface_list%flux_surfaces(i_surf)%n_pieces - 1) finished = .true.
          call get_next_surface_piece(node_list, element_list, surface_list%flux_surfaces(i_surf), i_piece, &
        			      n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces, found)
      	  ! --- If we didn't find the next piece, this means there could be another closed surface (unlikely but try once at least)
      	  if (found .eq. 0) then
    	    rr    = surface_list%flux_surfaces(i_surf)%s(3,i_piece)
    	    ss    = surface_list%flux_surfaces(i_surf)%t(3,i_piece)
    	    i_elm = surface_list%flux_surfaces(i_surf)%elm(i_piece)
    	    call interp_RZ(node_list,element_list,i_elm,rr,ss,R,dRR_dr,dRR_ds,dRR_drs,dRR_drr,dRR_dss, &
    							      Z,dZZ_dr,dZZ_ds,dZZ_drs,dZZ_drr,dZZ_dss)
	    
            rr2    = surface_list%flux_surfaces(i_surf)%s(1,parts_index(n_parts))
            ss2    = surface_list%flux_surfaces(i_surf)%t(1,parts_index(n_parts))
            i_elm2 = surface_list%flux_surfaces(i_surf)%elm(parts_index(n_parts))
            call interp_RZ(node_list,element_list,i_elm2,rr2,ss2,R2,dRR2_dr,dRR2_ds,dRR2_drs,dRR2_drr,dRR2_dss, &
      								 Z2,dZZ2_dr,dZZ2_ds,dZZ2_drs,dZZ2_drr,dZZ2_dss)
            
            distance = sqrt( (R-R2)**2.d0 + (Z-Z2)**2.d0 )
            if (distance .lt. accuracy) then
      	      found = j
	      nStart = i_piece + 1
              n_parts = n_parts + 1
              parts_index(n_parts) = nStart
      	      exit
            else
      	      write(*,'(A,1i4)') 'Warning! Failed to find parts of the surface',i_surf
      	      ier = 1
	    endif
      	  endif
	  
        enddo
      enddo
    
    ! -----------------------------------------------------
    ! --- In case there are edge pieces (ie. open surfaces)
    ! -----------------------------------------------------
    else
      if (mod(n_edge_pieces,2) .ne. 0) then
        write(*,'(A,1i4)') 'Warning! There are an odd number of edge pieces for surface',i_surf
	ier = 3
	return
      endif
      
      ! --- Start at first piece (check if there are isolated pieces before)
      nStart = 1
      if (n_isolated_pieces .gt. 0) nStart = parts_index(n_parts) + 1
      n_parts = n_parts + 1
      parts_index(n_parts) = nStart
      
      ! --- Loop over end pieces (two per surface)
      do i=1,n_edge_pieces/2
        
	! --- First swap piece so that it's at the begining of our new part.
        invert = .false.
        index1 = index_edge_pieces(i)
        index2 = nStart
        call swap_surface_pieces(surface_list%flux_surfaces(i_surf), index1, index2, invert, &
                                 n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)
        
	! --- Loop over all remaining pieces
        do i_piece = nStart, surface_list%flux_surfaces(i_surf)%n_pieces - 1
          
	  call get_next_surface_piece(node_list, element_list, surface_list%flux_surfaces(i_surf), i_piece, &
                                      n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces, found)
	  
	  ! --- If we didn't find the next piece, check if this should be the end
	  if (found .eq. 0) then
	    
	    ! --- Let's get to the next surface part
	    if (i .ne. n_edge_pieces/2) then
              n_parts = n_parts + 1
              parts_index(n_parts) = i_piece + 1
	      ! --- Make sure that this piece is an edge piece and that its index is > n_edge_pieces/2
              do j=1,n_edge_pieces
                if (index_edge_pieces(j) .eq. i_piece) then
            	  index_save = index_edge_pieces(j)
		  index_edge_pieces(j) = index_edge_pieces(n_edge_pieces/2+i)
		  index_edge_pieces(n_edge_pieces/2+i) = index_save
            	  exit
                endif
              enddo
	      exit
	    endif
	    
	    ! --- We ran out of edge pieces but there is more, meaning that there is probably a closed surface as well
	    finished = .true.
	    if ( (i .eq. n_edge_pieces/2) .and. (i_piece .ne. surface_list%flux_surfaces(i_surf)%n_pieces-1) ) then
	      finished = .false.
              nStart   = i_piece + 1
              n_parts  = n_parts + 1
              parts_index(n_parts) = i_piece + 1
	    endif
    
	  endif
	    
        enddo
	nStart = parts_index(n_parts)
      
        if (.not. finished) then
          ! --- There might be several closed surfaces
          do while(.not. finished)
            ! --- Loop over all remaining pieces
            do i_piece = nStart, surface_list%flux_surfaces(i_surf)%n_pieces - 1
      	      if (i_piece .eq. surface_list%flux_surfaces(i_surf)%n_pieces - 1) finished = .true.
              call get_next_surface_piece(node_list, element_list, surface_list%flux_surfaces(i_surf), i_piece, &
          				  n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces, found)
              ! --- If we didn't find the next piece, this means there could be another closed surface (unlikely but try once at least)
              if (found .eq. 0) then
      	  	rr    = surface_list%flux_surfaces(i_surf)%s(3,i_piece)
      	  	ss    = surface_list%flux_surfaces(i_surf)%t(3,i_piece)
      	  	i_elm = surface_list%flux_surfaces(i_surf)%elm(i_piece)
      	  	call interp_RZ(node_list,element_list,i_elm,rr,ss,R,dRR_dr,dRR_ds,dRR_drs,dRR_drr,dRR_dss, &
      	  							  Z,dZZ_dr,dZZ_ds,dZZ_drs,dZZ_drr,dZZ_dss)
      	  	
          	rr2    = surface_list%flux_surfaces(i_surf)%s(1,parts_index(n_parts))
          	ss2    = surface_list%flux_surfaces(i_surf)%t(1,parts_index(n_parts))
          	i_elm2 = surface_list%flux_surfaces(i_surf)%elm(parts_index(n_parts))
          	call interp_RZ(node_list,element_list,i_elm2,rr2,ss2,R2,dRR2_dr,dRR2_ds,dRR2_drs,dRR2_drr,dRR2_dss, &
          							     Z2,dZZ2_dr,dZZ2_ds,dZZ2_drs,dZZ2_drr,dZZ2_dss)
          	
          	distance = sqrt( (R-R2)**2.d0 + (Z-Z2)**2.d0 )
          	if (distance .lt. accuracy) then
          	  found = j
      	  	  nStart = i_piece + 1
          	  n_parts = n_parts + 1
          	  parts_index(n_parts) = nStart
          	  exit
          	else
          	  write(*,'(A,1i4)') 'Warning! Failed to find parts of the surface',i_surf
          	  ier = 1
      	  	endif
              endif
      	      
            enddo
          enddo
	  
        endif
      
      enddo
    
    endif
    
    ! --- Save parts indexes
    surface_list%flux_surfaces(i_surf)%n_parts = n_parts
    do i=1,n_parts
      surface_list%flux_surfaces(i_surf)%parts_index(i) = parts_index(i)
    enddo
    
    ! --- Print the number of pieces of the last part
    surface_list%flux_surfaces(i_surf)%parts_index(n_parts+1) = surface_list%flux_surfaces(i_surf)%n_pieces + 1
    
  enddo
  
  if (debug) call print_py_plot_ordered(node_list, element_list, surface_list)
  
  return

end subroutine reorder_flux_surfaces
  
  

! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------

!> This routine finds the next piece given a piece index of a fluxsurface
subroutine get_next_surface_piece(node_list, element_list, surface, i_piece, &
                                  n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces, found)
  
  use data_structure
  use reorder_surfaces_parameters
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),    intent(in)		:: node_list
  type (type_element_list), intent(in)		:: element_list
  type (type_surface),      intent(inout)	:: surface
  integer,                  intent(in)		:: i_piece
  integer,                  intent(inout)	:: found
  integer,                  intent(inout)	:: n_edge_pieces,     index_edge_pieces(2*n_parts_max)
  integer,                  intent(inout)	:: n_isolated_pieces, index_isolated_pieces(2*n_parts_max)
  
  ! --- Internal parameters
  integer	:: i, j
  integer	:: index1, index2
  logical	:: invert
  integer	:: i_elm
  integer	:: i_elm2
  real*8	:: rr,    ss
  real*8	:: rr2,   ss2
  real*8	:: R, dRR_dr, dRR_ds, dRR_drs, dRR_drr, dRR_dss
  real*8	:: R2,dRR2_dr,dRR2_ds,dRR2_drs,dRR2_drr,dRR2_dss
  real*8	:: Z, dZZ_dr, dZZ_ds, dZZ_drs, dZZ_drr, dZZ_dss
  real*8	:: Z2,dZZ2_dr,dZZ2_ds,dZZ2_drs,dZZ2_drr,dZZ2_dss
  real*8	:: distance
  
  found = 0
  
  ! --- Get last point of that surface piece
  rr	= surface%s(3,i_piece)
  ss	= surface%t(3,i_piece)
  i_elm = surface%elm(i_piece)
  call interp_RZ(node_list,element_list,i_elm,rr,ss,R,dRR_dr,dRR_ds,dRR_drs,dRR_drr,dRR_dss, &
  						    Z,dZZ_dr,dZZ_ds,dZZ_drs,dZZ_drr,dZZ_dss)

  ! --- Loop over all remaining pieces
  do i=i_piece+1,surface%n_pieces
  
    ! --- Try both ends of the piece
    do j=1,3,2
      rr2    = surface%s(j,i)
      ss2    = surface%t(j,i)
      i_elm2 = surface%elm(i)

      call interp_RZ(node_list,element_list,i_elm2,rr2,ss2,R2,dRR2_dr,dRR2_ds,dRR2_drs,dRR2_drr,dRR2_dss, &
  							   Z2,dZZ2_dr,dZZ2_ds,dZZ2_drs,dZZ2_drr,dZZ2_dss)
      
      distance = sqrt( (R-R2)**2.d0 + (Z-Z2)**2.d0 )
      if (distance .lt. accuracy) then
    	found = j
    	exit
      endif
    enddo
  
    ! --- Have we found the next piece?
    if (found .ne. 0) then
      index1 = i_piece + 1
      index2 = i
      invert = .false.
      if (found .eq. 3) invert = .true.
      call swap_surface_pieces(surface, index1, index2, invert, n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)
      exit
    endif

  enddo
  
  return

end subroutine get_next_surface_piece





! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This routine swaps two pieces of a given fluxsurface
subroutine swap_surface_pieces(surface, index1, index2, invert, n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)

  use data_structure
  use reorder_surfaces_parameters
  implicit none
  
  ! --- Routine parameters
  type (type_surface),      intent(inout)	:: surface
  integer,                  intent(in)		:: index1, index2
  logical,                  intent(in)		:: invert
  integer,                  intent(inout)	:: n_edge_pieces,     index_edge_pieces(2*n_parts_max)
  integer,                  intent(inout)	:: n_isolated_pieces, index_isolated_pieces(2*n_parts_max)
  
  ! --- Internal parameters
  integer	:: i
  integer	:: elm_save
  real*8	:: s_save(4), t_save(4)
  
  ! --- Save firsy piece
  elm_save  = surface%elm(index1)
  s_save(:) = surface%s(:,index1)
  t_save(:) = surface%t(:,index1)
  
  ! --- And swap
  if (invert) then
    ! --- First piece
    surface%elm(index1) = surface%elm(index2)
    surface%s(1,index1) = surface%s(3,index2)
    surface%t(1,index1) = surface%t(3,index2)
    surface%s(2,index1) = surface%s(4,index2)
    surface%t(2,index1) = surface%t(4,index2)
    surface%s(3,index1) = surface%s(1,index2)
    surface%t(3,index1) = surface%t(1,index2)
    surface%s(4,index1) = surface%s(2,index2)
    surface%t(4,index1) = surface%t(2,index2)
    
    ! --- Second piece
    surface%elm(index2) = elm_save
    surface%s(1,index2) = s_save(3)
    surface%t(1,index2) = t_save(3)
    surface%s(2,index2) = s_save(4)
    surface%t(2,index2) = t_save(4)
    surface%s(3,index2) = s_save(1)
    surface%t(3,index2) = t_save(1)
    surface%s(4,index2) = s_save(2)
    surface%t(4,index2) = t_save(2)
    
  else
    ! --- First piece
    surface%elm(index1) = surface%elm(index2)
    surface%s(:,index1) = surface%s(:,index2)
    surface%t(:,index1) = surface%t(:,index2)
    
    ! --- Second piece
    surface%elm(index2) = elm_save
    surface%s(:,index2) = s_save(:)
    surface%t(:,index2) = t_save(:)
  endif
  
  ! --- Then make sure the edge_pieces indices are swapped too
  if(n_edge_pieces .gt. 0) then
    do i=1,n_edge_pieces
      if (index_edge_pieces(i) .eq. index1) then
        index_edge_pieces(i) = index2
        continue
      endif
      if (index_edge_pieces(i) .eq. index2) then
        index_edge_pieces(i) = index1
        continue
      endif
    enddo
  endif
  
  ! --- And make sure the isolated_pieces indices are swapped too
  if(n_isolated_pieces .gt. 0) then
    do i=1,n_isolated_pieces
      if (index_isolated_pieces(i) .eq. index1) then
        index_isolated_pieces(i) = index2
        continue
      endif
      if (index_isolated_pieces(i) .eq. index2) then
        index_isolated_pieces(i) = index1
        continue
      endif
    enddo
  endif
  
  return
end subroutine swap_surface_pieces
  








! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This routine finds all the pieces of a given fluxsurface that have no neighbour
subroutine find_all_edge_pieces(node_list, element_list, surface, n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)



  use data_structure
  use reorder_surfaces_parameters
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),    intent(in)		:: node_list
  type (type_element_list), intent(in)		:: element_list
  type (type_surface),      intent(inout)	:: surface
  integer,                  intent(inout)	:: n_edge_pieces,     index_edge_pieces(2*n_parts_max)
  integer,                  intent(inout)	:: n_isolated_pieces, index_isolated_pieces(2*n_parts_max)
  
  ! --- Internal parameters
  integer	:: i, j, k1, k2
  logical	:: found(2), invert
  integer	:: i_elm
  integer	:: i_elm2
  real*8	:: rr,    ss
  real*8	:: rr2,   ss2
  real*8	:: R, dRR_dr, dRR_ds, dRR_drs, dRR_drr, dRR_dss
  real*8	:: R2,dRR2_dr,dRR2_ds,dRR2_drs,dRR2_drr,dRR2_dss
  real*8	:: Z, dZZ_dr, dZZ_ds, dZZ_drs, dZZ_drr, dZZ_dss
  real*8	:: Z2,dZZ2_dr,dZZ2_ds,dZZ2_drs,dZZ2_drr,dZZ2_dss
  real*8	:: distance
  
  n_isolated_pieces = 0
  n_edge_pieces     = 0
  
  ! --- Check each piece
  do i=1,surface%n_pieces
    found(1) = .false.
    found(2) = .false.
    ! --- Check both sides of the piece
    do k1=1,3,2
      ! --- Get last point of that surface piece
      rr    = surface%s(k1,i)
      ss    = surface%t(k1,i)
      i_elm = surface%elm(i)
      call interp_RZ(node_list,element_list,i_elm,rr,ss,R,dRR_dr,dRR_ds,dRR_drs,dRR_drr,dRR_dss, &
    	  					        Z,dZZ_dr,dZZ_ds,dZZ_drs,dZZ_drr,dZZ_dss)

      ! --- Loop over all other pieces
      do j=1,surface%n_pieces
        if (j .ne. i) then
          
	  ! --- Check both sides of the piece
          do k2=1,3,2
            rr2    = surface%s(k2,j)
            ss2    = surface%t(k2,j)
            i_elm2 = surface%elm(j)

            call interp_RZ(node_list,element_list,i_elm2,rr2,ss2,R2,dRR2_dr,dRR2_ds,dRR2_drs,dRR2_drr,dRR2_dss, &
          							 Z2,dZZ2_dr,dZZ2_ds,dZZ2_drs,dZZ2_drr,dZZ2_dss)
            
            distance = sqrt( (R-R2)**2.d0 + (Z-Z2)**2.d0 )
            if (distance .lt. accuracy) then
              if (k1 .eq. 1) found(1) = .true.
              if (k1 .eq. 3) found(2) = .true.
              exit
            endif
          enddo
	  
	  if ( (k1 .eq. 1) .and. (found(1)) ) exit
	  if ( (k1 .eq. 3) .and. (found(2)) ) exit
	  
        endif
      enddo
      
    enddo
    
    ! --- Have we found an isolated pieces?
    if ( (.not. found(1)) .and. (.not. found(2)) ) then
      n_isolated_pieces = n_isolated_pieces + 1
      index_isolated_pieces(n_isolated_pieces) = i
    else
      ! --- Have we found an edge pieces?
      if ( (.not. found(1)) .or. (.not. found(2)) ) then
        n_edge_pieces = n_edge_pieces + 1
        index_edge_pieces(n_edge_pieces) = i
        
        ! --- The edge pieces need to start at the edge. Invert the piece with itself if needed
        if (.not. found(2)) then
          invert = .true.
	  call swap_surface_pieces(surface, i, i, invert, n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)
        endif
      endif
    endif
      
  enddo
  



end subroutine find_all_edge_pieces





  








! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This plots fluxsurface not using the fact that they are ordered (ie. plotting piece after piece)
subroutine print_py_plot(node_list, element_list, surface_list)

  use data_structure
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),    intent(in)		:: node_list
  type (type_element_list), intent(in)		:: element_list
  type (type_surface_list), intent(inout)	:: surface_list
  
  ! --- Internal variables
  integer	:: i, j
  integer	:: i_elm
  real*8	:: rr,    ss
  real*8	:: R, dRR_dr, dRR_ds, dRR_drs, dRR_drr, dRR_dss
  real*8	:: Z, dZZ_dr, dZZ_ds, dZZ_drs, dZZ_drr, dZZ_dss
  
  open(100,file='plot_surfaces.py')
    write(100,'(A)')			      '#!/usr/bin/env python'
    write(100,'(A)')			      'import numpy as N'
    write(100,'(A)')			      'import pylab'
    write(100,'(A)')			      'def main():'
    write(100,'(A)')			      ' rplot = N.zeros(2)'
    write(100,'(A)')			      ' zplot = N.zeros(2)'
    do i=1,surface_list%n_psi
      do j=1,surface_list%flux_surfaces(i)%n_pieces
  	i_elm = surface_list%flux_surfaces(i)%elm(j)
  	rr    = surface_list%flux_surfaces(i)%s(1,j)
  	ss    = surface_list%flux_surfaces(i)%t(1,j)
  	call interp_RZ(node_list,element_list,i_elm,rr,ss,R,dRR_dr,dRR_ds,dRR_drs,dRR_drr,dRR_dss, &
        						  Z,dZZ_dr,dZZ_ds,dZZ_drs,dZZ_drr,dZZ_dss)
  	write(100,'(A,f15.4)')  	      ' rplot[0] = ',R
  	write(100,'(A,f15.4)')  	      ' zplot[0] = ',Z
        i_elm = surface_list%flux_surfaces(i)%elm(j)
  	rr    = surface_list%flux_surfaces(i)%s(3,j)
  	ss    = surface_list%flux_surfaces(i)%t(3,j)
  	call interp_RZ(node_list,element_list,i_elm,rr,ss,R,dRR_dr,dRR_ds,dRR_drs,dRR_drr,dRR_dss, &
        						  Z,dZZ_dr,dZZ_ds,dZZ_drs,dZZ_drr,dZZ_dss)
  	write(100,'(A,f15.4)')  	      ' rplot[1] = ',R
  	write(100,'(A,f15.4)')  	      ' zplot[1] = ',Z
  	write(100,'(A)')		      ' pylab.plot(rplot,zplot, "r")'
      enddo
    enddo
    write(100,'(A)')			      ' pylab.axis("equal")'
    write(100,'(A)')			      ' pylab.show()'
    write(100,'(A)')			      ' '
    write(100,'(A)')			      'main()'
  close(100)



end subroutine print_py_plot






  








! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This plots fluxsurface using the fact that they are ordered (ie. from the first piece of a part until its last one)
subroutine print_py_plot_ordered(node_list, element_list, surface_list)

  use data_structure
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),    intent(in)		:: node_list
  type (type_element_list), intent(in)		:: element_list
  type (type_surface_list), intent(inout)	:: surface_list
  
  ! --- Internal variables
  integer	:: i, j, k
  integer	:: i_elm
  real*8	:: rr,    ss
  real*8	:: R, dRR_dr, dRR_ds, dRR_drs, dRR_drr, dRR_dss
  real*8	:: Z, dZZ_dr, dZZ_ds, dZZ_drs, dZZ_drr, dZZ_dss
  
    open(101,file='plot_reordered_surfaces.py')
      write(101,'(A)')							'#!/usr/bin/env python'
      write(101,'(A)')							'import numpy as N'
      write(101,'(A)')							'import pylab'
      write(101,'(A)')							'def main():'
      write(101,'(A,i6,A)')						' r = N.zeros(',n_pieces_max,')'
      write(101,'(A,i6,A)')						' z = N.zeros(',n_pieces_max,')'
      do i=1,surface_list%n_psi
        do j=1,surface_list%flux_surfaces(i)%n_parts
          write(101,'(A,i6)')						' n_points = ', surface_list%flux_surfaces(i)%parts_index(j+1) &
	                                                                               -surface_list%flux_surfaces(i)%parts_index(j  ) + 1
          do k = surface_list%flux_surfaces(i)%parts_index(j), surface_list%flux_surfaces(i)%parts_index(j+1)-1
    	    rr	  = surface_list%flux_surfaces(i)%s(1,k)
    	    ss	  = surface_list%flux_surfaces(i)%t(1,k)
    	    i_elm = surface_list%flux_surfaces(i)%elm(k)
    	    call interp_RZ(node_list,element_list,i_elm,rr,ss,&
	                   R,dRR_dr,dRR_ds,dRR_drs,dRR_drr,dRR_dss, &
  			   Z,dZZ_dr,dZZ_ds,dZZ_drs,dZZ_drr,dZZ_dss)
    	    write(101,'(A,i6,A,f15.4)')					' r[',k-surface_list%flux_surfaces(i)%parts_index(j),'] = ',R
    	    write(101,'(A,i6,A,f15.4)')					' z[',k-surface_list%flux_surfaces(i)%parts_index(j),'] = ',Z
          enddo
    	  rr	= surface_list%flux_surfaces(i)%s(3,surface_list%flux_surfaces(i)%parts_index(j+1)-1)
    	  ss	= surface_list%flux_surfaces(i)%t(3,surface_list%flux_surfaces(i)%parts_index(j+1)-1)
    	  i_elm = surface_list%flux_surfaces(i)%elm(surface_list%flux_surfaces(i)%parts_index(j+1)-1)
    	  call interp_RZ(node_list,element_list,i_elm,rr,ss,&
	  		 R,dRR_dr,dRR_ds,dRR_drs,dRR_drr,dRR_dss, &
  	        	 Z,dZZ_dr,dZZ_ds,dZZ_drs,dZZ_drr,dZZ_dss)
    	  write(101,'(A,i6,A,f15.4)')					' r[',surface_list%flux_surfaces(i)%parts_index(j+1)-surface_list%flux_surfaces(i)%parts_index(j),'] = ',R
    	  write(101,'(A,i6,A,f15.4)')					' z[',surface_list%flux_surfaces(i)%parts_index(j+1)-surface_list%flux_surfaces(i)%parts_index(j),'] = ',Z
          write(101,'(A)')						' pylab.plot(r[0:n_points],z[0:n_points], "r")'
        enddo
      enddo
      write(101,'(A)')							' pylab.axis("equal")'
      write(101,'(A)')							' pylab.show()'
      write(101,'(A)')							' '
      write(101,'(A)')							'main()'
    close(101)

end subroutine print_py_plot_ordered




