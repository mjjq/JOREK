!> Allows to generate the boundary-element and boundary-node data structures for a given grid.
module mod_boundary

  implicit none
  private
  public boundary_from_grid

  contains
  !> Routine extracts the boundary information (boundary element and node lists)
  !! from the information stored in the grid (element and node lists).
  !! 
  !! Note: Grid nodes with boundary=3 (located at the edges of the boundary in
  !!       the divertor region) are intentionally added twice to the bnd_node_list.
  subroutine boundary_from_grid(node_list,element_list,bnd_node_list,bnd_elm_list,infos)

    use data_structure

    implicit none

    type (type_node_list),        intent(inout) :: node_list
    type (type_element_list),     intent(in)    :: element_list
    type (type_bnd_node_list),    intent(inout) :: bnd_node_list
    type (type_bnd_element_list), intent(inout) :: bnd_elm_list
    logical,                      intent(in)    :: infos ! If .true., writes bnd nodes and bnd elements in files 'boundary_nodes.dat' and 'boundaru_elements.dat' 

    integer :: i_elem         ! Element index in element_list.
    integer :: iside          ! Side of the element (1...4).
    integer :: iv1, iv2       ! Numbers of the two nodes on the element side (1...4).
    integer :: inode1, inode2, inode  ! Indices of the nodes in the node_list.
    integer :: b1, b2         ! Boundary properties of the two nodes.
    integer :: i, j, index_bnd, i_dof
    
    ! --- Empty the boundary node and element lists.
    bnd_node_list%n_bnd_nodes   = 0
    bnd_elm_list%n_bnd_elements = 0

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
          call add_bnd_elem( i_elem, iv1, iv2, inode1, inode2, iside, b1, b2, &
                             element_list, bnd_node_list, bnd_elm_list )
        end if

      end do

    end do

    ! --- Sort the boundary elements.
    call sort_bnd_elements( bnd_elm_list )


!=================== definition of boundary_index
    do i=1, bnd_node_list%n_bnd_nodes
       inode = bnd_node_list%bnd_node(i)%index_jorek
       node_list%node(inode)%boundary_index = i          ! the index in the bnd_node_list (not response_matrix)
    end do

!=================== define index in the response matrix (including corners whose nodes appear twice)

    index_bnd = 0

    do i=1, bnd_node_list%n_bnd_nodes

      bnd_node_list%bnd_node(i)%index_starwall = (/index_bnd+1,index_bnd+2/)
      bnd_node_list%bnd_node(i)%n_dof = 2

      index_bnd = index_bnd + bnd_node_list%bnd_node(i)%n_dof

    !  write(*,*)  'bnd index : ',i, bnd_node_list%bnd_node(i)%index_starwall

    enddo

    ! --- Output short (.false.) or verbose (.true.) boundary information.
    call log_bnd_info(infos, node_list, bnd_node_list, bnd_elm_list)

  end subroutine boundary_from_grid

  !> Adds a boundary element to the boundary element list.
  subroutine add_bnd_elem( i_elem, iv1, iv2, inode1, inode2, iside, b1, b2, &
    element_list, bnd_node_list, bnd_elm_list)

    use data_structure

    implicit none

    integer, intent(in)    :: i_elem         ! Index of the element in the element_list.
    integer, intent(in)    :: iv1, iv2       ! Index of the nodes of the element i_elem.
    integer, intent(in)    :: inode1, inode2 ! Indices of nodes in node_list.
    integer, intent(in)    :: iside          ! Which side of the element is located at the boundary?
    integer, intent(in)    :: b1, b2         ! Boundary information of the two nodes.
    type (type_element_list),     intent(in)    :: element_list
    type (type_bnd_node_list),    intent(inout) :: bnd_node_list
    type (type_bnd_element_list), intent(inout) :: bnd_elm_list

    integer :: ib1, ib2       ! Indices of the nodes in the bnd_node_list.
    integer :: idir

    ! --- Add the two boundary nodes to the list.
    call add_bnd_node( i_elem, iv1, inode1, iside, b1, ib1, bnd_node_list )
    call add_bnd_node( i_elem, iv2, inode2, iside, b2, ib2, bnd_node_list )

    bnd_elm_list%n_bnd_elements = bnd_elm_list%n_bnd_elements + 1

    ! --- Store vertex indices belonging to the boundary element.
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%vertex = (/ inode1, inode2 /)

    ! --- Store side of the element that is at the boundary.
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%side   = iside

    ! --- Store, which element the boundary element belongs to.
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%element = i_elem

    ! --- Store direction (directly connected to %side)
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%direction(:,1) = 1
    if ( mod(iside,2) == 1 ) then
      idir = 2
    else
      idir = 3
    end if
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%direction(:,2) = idir

    ! --- Store element size.
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%size(1,1) = element_list%element(i_elem)%size(iv1,1)
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%size(2,1) = element_list%element(i_elem)%size(iv2,1)

    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%size(1,2) = element_list%element(i_elem)%size(iv1,idir)
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%size(2,2) = element_list%element(i_elem)%size(iv2,idir)

    ! --- Store boundary node indices.
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%bnd_vertex(1) = ib1
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%bnd_vertex(2) = ib2

  end subroutine add_bnd_elem

  !> Adds a node to the bnd_node_list avoiding duplicates (except for boundary=3)
  !! and returns the index of the node in the bnd_node_list as bnd_vertex.
  subroutine add_bnd_node(i_elem, iv, inode, iside, boundary, bnd_vertex, bnd_node_list)

    use data_structure

    implicit none

    integer, intent(in)    :: i_elem         ! Index of the element in the element_list.
    integer, intent(in)    :: iv             ! Index of the node of the element i_elem.
    integer, intent(in)    :: inode          ! Index of the node in the node_list.
    integer, intent(in)    :: iside          ! Which side of the element is located at the boundary?
    integer, intent(in)    :: boundary       ! Boundary property of the node.
    integer, intent(inout) :: bnd_vertex     ! Index of the node in the bnd_node_list.
    type (type_bnd_node_list), intent(inout) :: bnd_node_list

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
    bnd_node_list%bnd_node(bnd_vertex)%index_starwall = bnd_vertex

    ! --- Store direction (directly connected to iside)
    bnd_node_list%bnd_node(bnd_vertex)%direction(1) = 1
    if ( mod(iside,2) == 1 ) then
      idir = 2
    else
      idir = 3
    end if
    bnd_node_list%bnd_node(bnd_vertex)%direction(2) = idir

  end subroutine add_bnd_node

  !> Outputs information about the boundary elements and nodes.
  subroutine log_bnd_info(verbose, node_list, bnd_node_list, bnd_elm_list)

    use data_structure

    implicit none

    logical, intent(in) :: verbose ! Output verbose information?
    type (type_node_list),        intent(in) :: node_list
    type (type_bnd_node_list),    intent(in) :: bnd_node_list
    type (type_bnd_element_list), intent(in) :: bnd_elm_list

    integer             :: i
    character(len=20)   :: s
    
    120 format(3x,77('-'))
    121 format(3X,A,I10,A)
    141 format(5X,A,I10,A)
    161 format(7X,A,I10,A)
    182 format(9X,A,15I10)
    183 format(9X,A,15ES10.2)

    write(*,*)
    write(*,120)
    write(*,121) 'BOUNDARY INFORMATION:'
    write(*,120)
    write(*,141) 'n_bnd_elements=', bnd_elm_list%n_bnd_elements
    write(*,141) 'n_bnd_nodes   =', bnd_node_list%n_bnd_nodes
    write(*,120)

    if ( verbose ) then

      write(*,*)
      write(*,120)
      write(*,141) 'BOUNDARY ELEMENTS:'
      write(*,120)
      do i = 1, bnd_elm_list%n_bnd_elements
        write(s,*) i
        write(*,161) '#'//trim(adjustl(s))//':'
        write(*,182) 'vertex        =', bnd_elm_list%bnd_element(i)%vertex
        write(*,182) 'bnd_vertex    =', bnd_elm_list%bnd_element(i)%bnd_vertex
        write(*,182) 'direction     =', bnd_elm_list%bnd_element(i)%direction
        write(*,182) 'element       =', bnd_elm_list%bnd_element(i)%element
        write(*,182) 'side          =', bnd_elm_list%bnd_element(i)%side
        write(*,183) 'size          =', bnd_elm_list%bnd_element(i)%size
      end do
      write(*,120)
      write(*,*)
      
      write(*,120)
      write(*,141) 'BOUNDARY NODES:'
      write(*,120)
      do i = 1, bnd_node_list%n_bnd_nodes
        write(s,*) i
        write(*,161) '#'//trim(adjustl(s))//':'
        write(*,182) 'index_jorek   =', bnd_node_list%bnd_node(i)%index_jorek
        write(*,182) 'index_starwall=', bnd_node_list%bnd_node(i)%index_starwall
        write(*,182) 'direction     =', bnd_node_list%bnd_node(i)%direction
      end do
      write(*,120)
      write(*,*)

      write(*,*) 'Writing boundary elements to "./boundary_elements.dat".'
      open(42, file='./boundary_elements.dat', status='replace', action='write')
      do i = 1, bnd_elm_list%n_bnd_elements
        write(42,*) node_list%node( bnd_elm_list%bnd_element(i)%vertex(1) )%x(1,:)
        write(42,*) node_list%node( bnd_elm_list%bnd_element(i)%vertex(2) )%x(1,:)
        write(42,*)
        write(42,*)
      end do
      close(42)

      write(*,*) 'Writing boundary nodes to "./boundary_nodes.dat".'
      open(42, file='./boundary_nodes.dat', status='replace', action='write')
      do i = 1, bnd_node_list%n_bnd_nodes
        write(42,*) node_list%node( bnd_node_list%bnd_node(i)%index_jorek )%x(1,:)
      end do
      close(42)

    end if

    write(*,*)

  end subroutine log_bnd_info

  !> Sorts the boundary elements.
  subroutine sort_bnd_elements( bnd_elm_list )

    use data_structure

    implicit none

    integer :: current_vertex, iter, iter_max
    integer :: ibnd_elem
    type(type_bnd_element)      :: bnd_elem
    type(type_bnd_element_list) :: sorted_bnd_element_list
    type (type_bnd_element_list), intent(inout) :: bnd_elm_list
    logical ::found_neighbour

    sorted_bnd_element_list%n_bnd_elements = 0

    current_vertex = bnd_elm_list%bnd_element(1)%vertex(1)

    iter = 0
    iter_max = 2*bnd_elm_list%n_bnd_elements
    found_neighbour = .true.

    L_A: do while ( (bnd_elm_list%n_bnd_elements > 0) .and. found_neighbour)

      iter = iter + 1
      if (iter .gt. iter_max) exit

      found_neighbour = .false.

      do ibnd_elem = 1, bnd_elm_list%n_bnd_elements

        bnd_elem = bnd_elm_list%bnd_element(ibnd_elem)

        if ( bnd_elem%vertex(2) == current_vertex ) call reverse_elem( bnd_elem )

        if ( bnd_elem%vertex(1) == current_vertex ) then
          current_vertex = bnd_elem%vertex(2)
          call add_elem( bnd_elem, sorted_bnd_element_list )
          call remove_elem( ibnd_elem, bnd_elm_list )
          found_neighbour = .true.
          exit
        end if

      end do

    end do L_A

    bnd_elm_list = sorted_bnd_element_list

    if (.not. found_neighbour) then
      write(*,*) 'BOUNDARY WARNING : NO NEIGHBOUR BND ELEMENT!',current_vertex
    endif

  end subroutine sort_bnd_elements

  !> Add the given boundary element to the end of the sorted list.
  subroutine add_elem( bnd_elem, sorted_bnd_element_list )

    use data_structure

    implicit none

    type(type_bnd_element),      intent(in)    :: bnd_elem
    type(type_bnd_element_list), intent(inout) :: sorted_bnd_element_list

    integer :: ibnd_elem

    ibnd_elem = sorted_bnd_element_list%n_bnd_elements + 1
    sorted_bnd_element_list%n_bnd_elements = ibnd_elem
    sorted_bnd_element_list%bnd_element(ibnd_elem) = bnd_elem

  end subroutine add_elem

  !> Reverse the given boundary element, i.e., exchange the nodes.
  subroutine reverse_elem( bnd_elem )

    use data_structure

    implicit none

    type(type_bnd_element), intent(inout) :: bnd_elem

    type(type_bnd_element) :: reversed_elem

    !write(*,*) '############################################################################'
    write(*,*) 'REVERSE_ELEM'
    !write(*,*) '############################################################################'

    reversed_elem%vertex(1) = bnd_elem%vertex(2)
    reversed_elem%vertex(2) = bnd_elem%vertex(1)

    reversed_elem%bnd_vertex(1) = bnd_elem%bnd_vertex(2)
    reversed_elem%bnd_vertex(2) = bnd_elem%bnd_vertex(1)

    reversed_elem%direction(1,:) = bnd_elem%direction(2,:)
    reversed_elem%direction(2,:) = bnd_elem%direction(1,:)

    reversed_elem%size(1,:) = bnd_elem%size(2,:)
    reversed_elem%size(2,:) = bnd_elem%size(1,:)
    
    bnd_elem = reversed_elem

  end subroutine reverse_elem

  !> Remove the ibnd_elem-th boundary element from bnd_elm_list.
  subroutine remove_elem( ibnd_elem, bnd_elm_list )

    use data_structure

    implicit none

    integer, intent(in) :: ibnd_elem
    type (type_bnd_element_list), intent(inout) :: bnd_elm_list

    integer :: nbnd_elem

    nbnd_elem = bnd_elm_list%n_bnd_elements

    bnd_elm_list%bnd_element(ibnd_elem:nbnd_elem-1) = bnd_elm_list%bnd_element(ibnd_elem+1:nbnd_elem)

    bnd_elm_list%n_bnd_elements = nbnd_elem - 1

  end subroutine remove_elem

end module mod_boundary
