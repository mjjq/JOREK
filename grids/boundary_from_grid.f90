subroutine boundary_from_grid()
! Routine extracts the boundary information (boundary element and node lists)
! from the information stored in the grid (element and node lists).
! 
! Note: Grid nodes with boundary=3 (located at the edges of the boundary in
!       the divertor region) are added twice to the bnd_node_list.
  
  use data_structure
  use nodes_elements
  
  implicit none
  
  integer :: i_elem         ! Element index in element_list.
  integer :: iside          ! Side of the element (1...4).
  integer :: iv1, iv2       ! Numbers of the two nodes on the element side (1...4).
  integer :: inode1, inode2 ! Indices of the nodes in the node_list.
  integer :: b1, b2         ! Boundary properties of the two nodes.
  
  ! --- Empty the boundary node and element lists.
  bnd_node_list%n_bnd_nodes    = 0
  boundary_list%n_bnd_elements = 0
  
  
  ! --- Go through all elements of the grid and find out which are located
  !     at the boundary. Build up the lists of boundary nodes and elements.
  
  do i_elem = 1, element_list%n_elements ! for each element
    
    do iside = 1, 4 ! for each side of the element (1...4)
    
      ! --- Indices of the nodes of the current element, which are located on this side.
      iv1 = iside
      iv2 = mod( iside, 4 ) + 1
      
      ! --- Indices of these nodes in the node_list.
      inode1 = element_list%element(i_elem)%vertex(iv1)
      inode2 = element_list%element(i_elem)%vertex(iv2)
      
      ! --- Boundary information of the two nodes.
      b1 = node_list%node(inode1)%boundary
      b2 = node_list%node(inode2)%boundary
      
      ! --- If the current side of this element is located at the boundary,
      !     add a new boundary element.
      if ( (b1 > 0) .AND. (b2 > 0) ) then
        call add_bnd_elem( i_elem, iv1, iv2, inode1, inode2, iside, b1, b2 )
      end if
      
    end do
    
  end do
  
  call log_bnd_info(.false.)
  !call log_bnd_info(.true.)  ! Output verbose boundary information
  
end subroutine boundary_from_grid







subroutine add_bnd_elem( i_elem, iv1, iv2, inode1, inode2, iside, b1, b2 )
! Adds a boundary element to the boundary element list.
  
  use data_structure
  use nodes_elements
  
  implicit none
  
  integer, intent(in)    :: i_elem         ! Index of the element in the element_list.
  integer, intent(in)    :: iv1, iv2       ! Index of the nodes of the element i_elem.
  integer, intent(in)    :: inode1, inode2 ! Indices of nodes in node_list.
  integer, intent(in)    :: iside          ! Which side of the element is located at the boundary?
  integer, intent(in)    :: b1, b2         ! Boundary information of the two nodes.
  
  integer :: ib1, ib2       ! Indices of the nodes in the bnd_node_list.
  integer :: idir
  
  ! --- Add the two boundary nodes to the list.
  call add_bnd_node( i_elem, iv1, inode1, iside, b1, ib1 )
  call add_bnd_node( i_elem, iv2, inode2, iside, b2, ib2 )
  
  boundary_list%n_bnd_elements = boundary_list%n_bnd_elements + 1
  
  ! --- Store vertex indices belonging to the boundary element.
  boundary_list%bnd_element(boundary_list%n_bnd_elements)%vertex = (/ inode1, inode2 /)
  
  ! --- Store side of the element that is at the boundary.
  boundary_list%bnd_element(boundary_list%n_bnd_elements)%side   = iside
  
  ! --- Store, which element the boundary element belongs to.
  boundary_list%bnd_element(boundary_list%n_bnd_elements)%element = i_elem
  
  ! --- Store direction (directly connected to %side) ############ REMOVE? #############
  boundary_list%bnd_element(boundary_list%n_bnd_elements)%direction(:,1) = 1
  if ( mod(iside,2) == 1 ) then
    idir = 2
  else
    idir = 3
  end if
  boundary_list%bnd_element(boundary_list%n_bnd_elements)%direction(:,2) = idir
  
  ! --- Store element size.
  boundary_list%bnd_element(boundary_list%n_bnd_elements)%size(1,1) = element_list%element(i_elem)%size(iv1,1)
  boundary_list%bnd_element(boundary_list%n_bnd_elements)%size(2,1) = element_list%element(i_elem)%size(iv2,1)
  
  boundary_list%bnd_element(boundary_list%n_bnd_elements)%size(1,2) = element_list%element(i_elem)%size(iv1,idir)
  boundary_list%bnd_element(boundary_list%n_bnd_elements)%size(2,2) = element_list%element(i_elem)%size(iv2,idir)
  
  ! --- Store boundary node indices.
  boundary_list%bnd_element(boundary_list%n_bnd_elements)%bnd_vertex(1) = ib1
  boundary_list%bnd_element(boundary_list%n_bnd_elements)%bnd_vertex(2) = ib2
  
end subroutine add_bnd_elem






subroutine add_bnd_node(i_elem, iv, inode, iside, boundary, bnd_vertex)
! Adds a node to the bnd_node_list avoiding duplicates (except for boundary=3)
! and returns the index of the node in the bnd_node_list as bnd_vertex.
  
  use data_structure
  use nodes_elements
  
  implicit none
  
  integer, intent(in)    :: i_elem         ! Index of the element in the element_list.
  integer, intent(in)    :: iv             ! Index of the node of the element i_elem.
  integer, intent(in)    :: inode          ! Index of the node in the node_list.
  integer, intent(in)    :: iside          ! Which side of the element is located at the boundary?
  integer, intent(in)    :: boundary       ! Boundary property of the node.
  integer, intent(out)   :: bnd_vertex     ! Index of the node in the bnd_node_list.
  
  integer :: i, idir

  ! --- Make sure the node is not in the boundary node list yet (except for boundary=3).
  if ( boundary /= 3 ) then
    do i = 1, bnd_node_list%n_bnd_nodes
      if ( bnd_node_list%bnd_node(i)%index_jorek == inode ) then
        bnd_vertex = i ! Node is already in the list, return its index.
        return
      end if
    end do
  end if
  
  ! --- Add the node.
  bnd_node_list%n_bnd_nodes = bnd_node_list%n_bnd_nodes + 1
  bnd_vertex = bnd_node_list%n_bnd_nodes
  
  bnd_node_list%bnd_node(bnd_vertex)%index_jorek    = inode
  bnd_node_list%bnd_node(bnd_vertex)%index_starwall = bnd_vertex !################ Change later maybe ###################
  
  ! --- Store direction (directly connected to iside)
  bnd_node_list%bnd_node(bnd_vertex)%direction(1) = 1
  if ( mod(iside,2) == 1 ) then
    idir = 2
  else
    idir = 3
  end if
  bnd_node_list%bnd_node(bnd_vertex)%direction(2) = idir
  
end subroutine add_bnd_node






subroutine log_bnd_info(verbose)
! Outputs information about the boundary elements and nodes.
  
  use data_structure
  use nodes_elements
  
  implicit none
  
  logical, intent(in) :: verbose ! Output verbose information?
  integer             :: i
  character(len=20)   :: s
  
  121 format(1X,A,I10,A)
  141 format(3X,A,I10,A)
  161 format(5X,A,I10,A)
  182 format(7X,A,15I10)
  183 format(7X,A,15ES10.2)
  
  write(*,*)
  
  write(*,121) 'BOUNDARY INFORMATION:'
  
  write(*,141) 'n_bnd_elements=', boundary_list%n_bnd_elements
  write(*,141) 'n_bnd_nodes   =', bnd_node_list%n_bnd_nodes
  
  if ( verbose ) then
  
    write(*,*)
    write(*,141) 'BOUNDARY ELEMENTS:'
    do i = 1, boundary_list%n_bnd_elements
      write(s,*) i
      write(*,161) '#'//trim(adjustl(s))//':'
      write(*,182) 'vertex        =', boundary_list%bnd_element(i)%vertex
      write(*,182) 'bnd_vertex    =', boundary_list%bnd_element(i)%bnd_vertex
      write(*,182) 'direction     =', boundary_list%bnd_element(i)%direction
      write(*,182) 'element       =', boundary_list%bnd_element(i)%element
      write(*,182) 'side          =', boundary_list%bnd_element(i)%side
      write(*,183) 'size          =', boundary_list%bnd_element(i)%size
    end do
  
    write(*,*)
    write(*,141) 'BOUNDARY NODES:'
    do i = 1, bnd_node_list%n_bnd_nodes
      write(s,*) i
      write(*,161) '#'//trim(adjustl(s))//':'
      write(*,182) 'index_jorek   =', bnd_node_list%bnd_node(i)%index_jorek
      write(*,182) 'index_starwall=', bnd_node_list%bnd_node(i)%index_starwall
      write(*,182) 'direction     =', bnd_node_list%bnd_node(i)%direction
    end do
  
  end if

  write(*,*)
  
end subroutine log_bnd_info

