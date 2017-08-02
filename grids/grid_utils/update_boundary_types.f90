subroutine update_boundary_types(element_list,node_list, across_xpoint)
  
  ! -------------------------------------------------------------------
  ! -------------------------------------------------------------------
  ! --------- IMPORTANT NOTE! THIS ROUTINE ASSUMES THAT NODES ---------
  ! --------- 1-4 ARE XPOINT NODES (AND ALSO 5-8 IF XCASE=3)  ---------
  ! --------- IF YOU ARE BUILDING A PATCH THAT DOES NOT HAVE  ---------
  ! --------- ANY XPOINT IN IT, THEN MAKE SURE THE ELEMENT    ---------
  ! --------- VERTEX INDEXING STARTS AT N_XPOINT+1            ---------
  ! -------------------------------------------------------------------
  ! -------------------------------------------------------------------
  
  
  use mod_parameters
  use data_structure
  use phys_module, only: xcase
  
  implicit none
  
  ! --- Routine variables
  type(type_node_list),    intent(inout) :: node_list	   !< list of grid nodes
  type(type_element_list), intent(inout) :: element_list   !< list of finite elements
  logical,                 intent(in)    :: across_xpoint  !< if true, then we assume we have 4 xpoint nodes, otherwise, only two
                                                           ! eg, you go from 1 to 4, or you go from 1 to 2 (see create_x_node.f90)
  
  ! --- Internal variables
  integer :: i_elm,       i_vertex,  i_node,  i_node_prev
  integer :: i_elm2,      i_vertex2, i_node2
  integer :: elm_sum
  integer :: i_elm_first, i_vertex_first
  integer :: i_elm_now,   i_vertex_now, i_vertex_next, i_vertex_save
  integer :: i_elm_prev
  integer :: iter, n_xpoints
  logical :: found_first, found_next
  logical, parameter :: debug  = .false.
  logical, parameter :: debug2 = .false.
  
  ! --- Some printouts?
  if (debug) then
    write(*,*)'----------------------------------------------------------'
    write(*,*)'---------------- Updating boundary types. ----------------'
    write(*,*)'----------------------------------------------------------'
  endif
  
  if (xcase .eq. 3) then
    n_xpoints = 8
  else
    n_xpoints = 4
  endif
  
  ! --- Initialise first!
  do i_node=1,node_list%n_nodes
    node_list%node(i_node)%boundary = 0
  enddo
  
  ! --- Boundary check
  ! --- Find first element with a boundary
  !write(*,*)'Looking for first boundary element...'
  i_elm_first = 0
  found_first = .false.
  do i_elm=1,element_list%n_elements
    do i_vertex=1,4
      i_node = element_list%element(i_elm)%vertex(i_vertex)
      if ( (i_node .le. 4) .and. (xcase .ne. 3) ) cycle
      if ( (i_node .le. 8) .and. (xcase .eq. 3) ) cycle
      call adjacent_elements(element_list,node_list,i_elm,i_vertex,3,elm_sum)
      ! --- We want a type-3 boundary to start with (note this also make it safer if we have axis nodes on our grid)
      if (elm_sum .eq. 0) then
        ! --- Use one of the adjacent nodes to start from
        i_vertex_next = mod(i_vertex+2,4) + 1
        call adjacent_elements(element_list,node_list,i_elm,i_vertex_next,3,elm_sum)
        if (elm_sum .eq. 1) then
          i_vertex_first = i_vertex_next
          i_elm_first    = i_elm
          found_first = .true.
          exit
        else
          ! --- Try the other direction (this might be the xpoint when we have an individual leg alone)
          i_vertex_next = mod(i_vertex,4) + 1
          call adjacent_elements(element_list,node_list,i_elm,i_vertex_next,3,elm_sum)
          if (elm_sum .eq. 1) then
            i_vertex_first = i_vertex_next
            i_elm_first    = i_elm
            found_first = .true.
            exit
          else
            if (debug) write(*,*)'Something strange when looking for the first element'
            if (debug) write(*,*)'Found a type 3 but its neighbours are not bnd...'
          endif  
        endif  
      endif
    enddo
    if (found_first) exit
  enddo
  if (.not. found_first) then
    write(*,*)'Could not find first boundary element. Aborting...'
    return
  endif
  if (debug) write(*,'(A,i6,2f10.3)')'Found first bnd elm    :',i_node,node_list%node(i_node)%x(1,1:2)
  
  ! --- Make sure our boundary is coherent (should not last longer than the number of elements)
  found_first  = .false.
  i_elm_now    = i_elm_first
  i_vertex_now = i_vertex_first
  i_elm_prev   = 0
  if (debug) write(*,*)'Going around boundary...'
  do iter=1,element_list%n_elements
    !write(*,*)'iter...',iter
    
    ! --- What type is this boundary?
    call adjacent_elements(element_list,node_list,i_elm_now,i_vertex_now,2,elm_sum)
    if ( (elm_sum .eq. 0) .and. (i_elm_prev .eq. 0) ) then
      write(*,*)'Impossible, we should always start on a non-corner node. Aborting...'
      exit
    endif
    i_node = element_list%element(i_elm_now)%vertex(i_vertex_now)
    if (debug2) write(*,'(A,i6,2f10.3)')'Starting on new element:',i_elm_now
    if (debug2) write(*,'(A,i6,2f10.3)')'On new elm, starting at:',i_node,node_list%node(i_node)%x(1,1:2)
    if (elm_sum .eq. 1) node_list%node(i_node)%boundary = 1
    if (elm_sum .eq. 2) node_list%node(i_node)%boundary = 3
    
    ! --- Find the next boundary node on this element
    i_vertex_next = mod(i_vertex_now,4) + 1 ! it should always be the same direction (ie. +1, not -1)... or should it?
    call adjacent_elements(element_list,node_list,i_elm_now,i_vertex_next,3,elm_sum)
    if (debug2) write(*,'(A,i6,2f10.3,i2)')'Trying next bnd node   :',element_list%element(i_elm_now)%vertex(i_vertex_next),&
                                            node_list%node(element_list%element(i_elm_now)%vertex(i_vertex_next))%x(1,1:2),elm_sum
    if (elm_sum .eq. 3) then
      ! --- Something wrong: we went back inside grid. change direction...'
      i_vertex_next = mod(i_vertex_now+2,4) + 1
      call adjacent_elements(element_list,node_list,i_elm_now,i_vertex_next,3,elm_sum)
      if (debug2) write(*,'(A,i6,2f10.3,i2)')'Wrong direction, other :',element_list%element(i_elm_now)%vertex(i_vertex_next),&
                                              node_list%node(element_list%element(i_elm_now)%vertex(i_vertex_next))%x(1,1:2),elm_sum
      if (elm_sum .eq. 3) then
        write(*,*)'Something wrong: we went back inside grid. Aborting...'
        exit
      endif
    endif
    i_node_prev = i_node
    i_node = element_list%element(i_elm_now)%vertex(i_vertex_next)
    if (elm_sum .eq. 0) node_list%node(i_node)%boundary = 3
    if (elm_sum .eq. 1) node_list%node(i_node)%boundary = 1
    if (elm_sum .eq. 2) node_list%node(i_node)%boundary = 3
    ! --- On a corner, we need to repeat (careful, xpoints also have only one elm)
    if ( (elm_sum .eq. 0) .and. (i_node .gt. n_xpoints) ) then
      i_vertex_save = i_vertex_next
      i_vertex_next = mod(i_vertex_next,4) + 1 ! it should always be the same direction (ie. +1, not -1)
      call adjacent_elements(element_list,node_list,i_elm_now,i_vertex_next,3,elm_sum)
      if (debug2) write(*,'(A,i6,2f10.3,i2)')'On corner, next node is:',element_list%element(i_elm_now)%vertex(i_vertex_next),&
                                              node_list%node(element_list%element(i_elm_now)%vertex(i_vertex_next))%x(1,1:2),elm_sum
      if ( (elm_sum .eq. 3) .or. (element_list%element(i_elm_now)%vertex(i_vertex_next) .eq. i_node_prev) ) then
        ! --- Something wrong: we went back inside grid. change direction...'
        i_vertex_next = mod(i_vertex_save+2,4) + 1
        call adjacent_elements(element_list,node_list,i_elm_now,i_vertex_next,3,elm_sum)
        if (debug2) write(*,'(A,i6,2f10.3,i2)')'Wrong direction, other :',element_list%element(i_elm_now)%vertex(i_vertex_next),&
                                                node_list%node(element_list%element(i_elm_now)%vertex(i_vertex_next))%x(1,1:2),elm_sum
        if (elm_sum .eq. 3) then
          write(*,*)'Something wrong at corner: we went back inside grid. Aborting...'
          exit
        endif
      endif
      i_node = element_list%element(i_elm_now)%vertex(i_vertex_next)
      if (elm_sum .eq. 1) node_list%node(i_node)%boundary = 1
      if (elm_sum .eq. 2) node_list%node(i_node)%boundary = 3
      if (elm_sum .eq. 0) then
        write(*,*)'Impossible, you cannot have two corner nodes on the same element. Aborting...'
        exit
      endif
    endif
    if (debug2) write(*,'(A,i6,2f10.3)')'Found next bnd node    :',i_node,node_list%node(i_node)%x(1,1:2)
    i_vertex_now = i_vertex_next
  
    ! --- Find the next element (it should have at least 2 boundary nodes!)
    found_next = .false.
    i_node = element_list%element(i_elm_now)%vertex(i_vertex_now)
    if (i_node .le. n_xpoints) then
      if (across_xpoint) then
        if (mod(i_node,4) .eq. 0) then
          i_node = i_node - 3
        elseif (mod(i_node,4) .eq. 1) then
          i_node = i_node + 3
        else
          if (debug) write(*,*)'Problem: found wrong Xpoint node on our path'
        endif
      else
        if (mod(i_node,2) .eq. 0) then
          i_node = i_node - 1
        else
          i_node = i_node + 1
        endif
      endif
    endif
    do i_elm2=1,element_list%n_elements
      if (i_elm2 .eq. i_elm_now) cycle
      if (i_elm2 .eq. i_elm_prev) cycle
      do i_vertex2=1,4
        i_node2 = element_list%element(i_elm2)%vertex(i_vertex2)
        if (i_node2 .eq. i_node) then
          ! --- If this is one, check that the next node is also a boundary
          i_vertex_next = mod(i_vertex2,4) + 1 ! it should always be the same direction (ie. +1, not -1)
          call adjacent_elements(element_list,node_list,i_elm2,i_vertex_next,3,elm_sum)
          if (elm_sum .le. 2) then
            found_next   = .true.
            i_elm_prev   = i_elm_now
            i_elm_now    = i_elm2
            i_vertex_now = i_vertex2
          else
            ! --- Something wrong: change direction...'
            i_vertex_next = mod(i_vertex2+2,4) + 1
            call adjacent_elements(element_list,node_list,i_elm2,i_vertex_next,3,elm_sum)
            if (elm_sum .le. 2) then
              found_next   = .true.
              i_elm_prev   = i_elm_now
              i_elm_now    = i_elm2
              i_vertex_now = i_vertex2
            endif
          endif
          exit
        endif
      enddo
      if (found_next) exit
    enddo
    if (.not. found_next) then
      write(*,*)'Could not find next element. Aborting...',i_node
      exit
    else
      if (debug2) write(*,'(A,i6,2f10.3)')'Found next bnd elm/node:',i_node2,node_list%node(i_node2)%x(1,1:2)
    endif
    
    ! --- Check if we looped all around
    if (i_elm_now .eq. i_elm_first) then
      if (debug) write(*,*)'Finished boundary'
      if (i_vertex_now .ne. i_vertex_first) write(*,*)'But finished on wrong node!!!'
      exit
    endif
    
  enddo
  
  ! --- Now, we need to check non-corner boundaries to see if they are on the s-side or the t-side
  do i_elm=1,element_list%n_elements
    do i_vertex=1,4
      i_node = element_list%element(i_elm)%vertex(i_vertex)
      if (node_list%node(i_node)%boundary .eq. 1) then
        do i_vertex2=1,4
          if (i_vertex2 .eq. i_vertex) cycle
          i_node2 = element_list%element(i_elm)%vertex(i_vertex2)
          if (node_list%node(i_node)%boundary .gt. 0) then
	    if (i_vertex+i_vertex2 .ne. 5) node_list%node(i_node)%boundary = 2
	    exit
	  endif
        enddo
      endif
    enddo
  enddo
  
  if (debug) then
    open(101,file='plot_bnd_points.py')
      write(101,'(A)')                '#!/usr/bin/env python'
      write(101,'(A)')                'import numpy as N'
      write(101,'(A)')                'import pylab'
      write(101,'(A)')                'def main():'
      elm_sum = 0
      write(101,'(A,i6,A)')           ' r = N.zeros(',node_list%n_nodes,')'
      write(101,'(A,i6,A)')           ' z = N.zeros(',node_list%n_nodes,')'
      do i_elm=1,element_list%n_elements
        do i_vertex=1,4
          i_node = element_list%element(i_elm)%vertex(i_vertex)
          if (node_list%node(i_node)%boundary .eq. 0) cycle
          write(101,'(A,i6,A,f15.4)') ' r[',elm_sum,'] = ',node_list%node(i_node)%x(1,1)
          write(101,'(A,i6,A,f15.4)') ' z[',elm_sum,'] = ',node_list%node(i_node)%x(1,2)
          elm_sum = elm_sum + 1
        enddo
      enddo
      write(101,'(A,i6)')             ' n_points = ',elm_sum
      write(101,'(A)')                ' pylab.plot(r[0:n_points],z[0:n_points], "rx")'
      write(101,'(A)')                ' pylab.axis("equal")'
      write(101,'(A)')                ' pylab.show()'
      write(101,'(A)')                ' '
      write(101,'(A)')                'main()'
    close(101)
  endif
  
  return
  
end subroutine update_boundary_types



subroutine adjacent_elements(element_list,node_list,i_elm,i_vertex,side_max, elm_sum)

  use mod_parameters
  use data_structure
  
  implicit none
  
  ! --- Routine variables
  type(type_node_list),    intent(in)  :: node_list	   !< list of grid nodes
  type(type_element_list), intent(in)  :: element_list   !< list of finite elements
  integer,                 intent(in)  :: i_elm,i_vertex,side_max
  integer,                 intent(out) :: elm_sum
  
  ! --- Internal variables
  integer :: i_elm2, i_vertex2, i_node, i_node2
  
  elm_sum = 0
  i_node = element_list%element(i_elm)%vertex(i_vertex)
  do i_elm2=1,element_list%n_elements
    if (i_elm2 .eq. i_elm) cycle
    do i_vertex2=1,4
      i_node2 = element_list%element(i_elm2)%vertex(i_vertex2)
      if (i_node2 .eq. i_node) then
    	elm_sum = elm_sum + 1
    	exit
      endif
    enddo
    if (elm_sum .eq. side_max) exit
  enddo
  
  return
end subroutine adjacent_elements



