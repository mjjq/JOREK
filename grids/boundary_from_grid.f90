subroutine boundary_from_grid()
! Routine extracts the boundary information (boundary element and node lists)
! from the information stored in the grid (element and node lists).
  
  
  
  use data_structure
  use nodes_elements
  
  
  
  implicit none
  
  
  
  integer :: i_elem         ! Element index in element_list.
  integer :: iside          ! Side of the element (1...4).
  integer :: iv1, iv2       ! Numbers of the two nodes on the element side (1...4).
  integer :: inode1, inode2 ! Indices of the nodes in the node_list.
  integer :: b1, b2         ! Boundary properties of the two nodes.
  integer :: ib1, ib2       ! Indices of the nodes in the bnd_node_list.
  

  
  ! --- Empty the boundary node and element lists.
  bnd_node_list%n_bnd_nodes    = 0
  boundary_list%n_bnd_elements = 0
  
  
  ! --- Go through all elements of the grid and find out which are located
  !     at the boundary. Build up the lists of boundary nodes and elements.
  
  do i_elem = 1, element_list%n_elements ! for each element
    
    do iside = 1, 4 ! for each side of the element (1...4)
    
      ! --- Indices of the nodes of the current elements that are located on this side.
      iv1 = iside
      iv2 = mod( iside, 4 ) + 1
      
      ! --- Indices of the nodes in the node_list.
      inode1 = element_list%element(i_elem)%vertex(iv1)
      inode2 = element_list%element(i_elem)%vertex(iv2)
      
      ! --- Boundary information of the two nodes.
      b1 = node_list%node(inode1)%boundary
      b2 = node_list%node(inode2)%boundary
      
      if ( (b1 > 0) .AND. (b2 > 0) ) then
        ! The current side of this element is located at the boundary, so...
        
        ! --- Add the two boundary nodes to the list.
        call add_bnd_node( inode1, b1, ib1 )
        call add_bnd_node( inode2, b2, ib2 )
        
        ! --- Add a new boundary element.
        call add_bnd_elem( i_elem, inode1, inode2, iside, ib1, ib2 )
        
      end if
      
    end do
    
  end do
  
  write(*,*) 'bnd_node_list%n_bnd_nodes    =', bnd_node_list%n_bnd_nodes
  write(*,*) 'boundary_list%n_bnd_elements =', boundary_list%n_bnd_elements

  
  
  
  
  
  
  CONTAINS
  
  
  
  
  
  
  subroutine add_bnd_node(inode, boundary, bnd_vertex)
  ! Adds a node to the bnd_node_list avoiding duplicates (except for boundary=3)
  ! and returns the index of the node in the bnd_node_list as bnd_vertex.
    
    integer, intent(in)    :: inode      ! Index of the node in the node_list.
    integer, intent(in)    :: boundary   ! Boundary property of the node.
    integer, intent(out)   :: bnd_vertex ! Index of the node in the bnd_node_list.
    
    integer :: i
  
    ! Make sure the node is not in the boundary node list yet (except for boundary=3).
    if ( boundary /= 3 ) then
      do i = 1, bnd_node_list%n_bnd_nodes
        if ( bnd_node_list%bnd_node(i)%index_jorek == inode ) then
          bnd_vertex = i
          return
        end if
      end do
    end if
    
    ! Add the node.
    bnd_node_list%n_bnd_nodes = bnd_node_list%n_bnd_nodes + 1
    bnd_vertex = bnd_node_list%n_bnd_nodes
    
    bnd_node_list%bnd_node(bnd_vertex)%index_jorek    = inode
    bnd_node_list%bnd_node(bnd_vertex)%index_starwall = bnd_vertex !################ Change later maybe ###################
    
  end subroutine add_bnd_node



  subroutine add_bnd_elem( i_elem, inode1, inode2, iside, ib1, ib2 )
  ! Adds a boundary element to the boundary element list.
    
    integer, intent(in) :: i_elem         ! Index of the element in the element_list.
    integer, intent(in) :: inode1, inode2 ! Indices of nodes in node_list.
    integer, intent(in) :: iside          ! Which side of the element is located at the boundary?
    integer, intent(in) :: ib1, ib2       ! Indices of nodes in bnd_node_list.
    
    integer :: idir
    
    boundary_list%n_bnd_elements = boundary_list%n_bnd_elements + 1
    
    !   --- Store vertex indices belonging to the boundary element.
    boundary_list%bnd_element(boundary_list%n_bnd_elements)%vertex = (/ inode1, inode2 /)
    
    !   --- Store side of the element that is at the boundary.
    boundary_list%bnd_element(boundary_list%n_bnd_elements)%side   = iside
    
    !   --- Store, which element the boundary element belongs to.
    boundary_list%bnd_element(boundary_list%n_bnd_elements)%element = i_elem
    
    !   --- Store direction (directly connected to %side)
    boundary_list%bnd_element(boundary_list%n_bnd_elements)%direction(:,1) = 1
    if ( mod(iside,2) == 1 ) then
      idir = 2
    else
      idir = 3
    end if
    boundary_list%bnd_element(boundary_list%n_bnd_elements)%direction(:,2) = idir
    
    !   --- Store element size.
    boundary_list%bnd_element(boundary_list%n_bnd_elements)%size(1,1) = element_list%element(i_elem)%size(iv1,1)
    boundary_list%bnd_element(boundary_list%n_bnd_elements)%size(2,1) = element_list%element(i_elem)%size(iv2,1)
    
    boundary_list%bnd_element(boundary_list%n_bnd_elements)%size(1,2) = element_list%element(i_elem)%size(iv1,idir)
    boundary_list%bnd_element(boundary_list%n_bnd_elements)%size(2,2) = element_list%element(i_elem)%size(iv2,idir)
    
    !   --- Store boundary node indices.
    boundary_list%bnd_element(boundary_list%n_bnd_elements)%bnd_vertex(1) = ib1
    boundary_list%bnd_element(boundary_list%n_bnd_elements)%bnd_vertex(2) = ib2
    
  end subroutine add_bnd_elem
  
  
  
end subroutine boundary_from_grid
