subroutine update_neighbours_basic(element_list,node_list)
  
  use data_structure
  
  implicit none

  ! --- Routine variables
  type (type_element_list), intent(inout) :: element_list
  type (type_node_list),    intent(in)    :: node_list

  ! --- Internal variables
  integer :: i_elm,   j_elm
  integer :: i_side,  j_side
  integer :: i_node1, j_node1
  integer :: i_node2, j_node2
  real*8  :: progress

  write(*,*)'----------------------------------------------------------'
  write(*,*)'----- Updating Neighbours. This may take a while... ------'
  write(*,*)'----------------------------------------------------------'
  
  ! --- Loop on each element
!$omp parallel do default(private) &
!$omp   shared(element_list,node_list)
  do i_elm = 1, element_list%n_elements
  
    !progress = 1.d2 * float(i_elm) / float(element_list%n_elements)
    !progress = max(0.d0,progress)
    !progress = min(1.d2,progress)
    !write(*,"(' Processing  ... ',I0,'%',A,$)") int(progress),CHAR(13)
    
    ! --- Loop on each side
    do i_side = 1, n_vertex_max
      i_node1 = element_list%element(i_elm)%vertex(i_side)
      i_node2 = element_list%element(i_elm)%vertex(mod(i_side,4)+1)
      ! --- Loop on all other element
      do j_elm = 1, element_list%n_elements
        if (j_elm .eq. i_elm) cycle
        call are_elements_neighbours(element_list%element(i_elm), i_side, &
                                     element_list%element(j_elm), node_list, j_side)
        if (j_side .gt. 0) then
          element_list%element(i_elm)%neighbours(i_side) = j_elm
          element_list%element(j_elm)%neighbours(j_side) = i_elm
          exit
        endif
      enddo
    enddo
  
  enddo
!$omp end parallel do

  !write(*,*) 'Processing  ... 100'
  !write(*,*) 'finished updating neighbours'
  
  return

end subroutine update_neighbours_basic




subroutine are_elements_neighbours(elm1, i_side, elm2, node_list, j_side)

  use data_structure
  use mod_parameters

  implicit none
  
  ! --- Routine variables
  type (type_node_list), intent(in)  :: node_list
  type (type_element),   intent(in)  :: elm1, elm2
  integer,               intent(in)  :: i_side
  integer,               intent(out) :: j_side
  
  ! --- Internal variables
  integer :: i_node1, i_node2
  integer :: j_node1, j_node2
  integer :: i

  ! --- Initialise (not neighbours)
  j_side = 0

  i_node1 = elm1%vertex(i_side)
  i_node2 = elm1%vertex(mod(i_side,4)+1)
  
  ! --- Loop on each side
  do  i = 1, n_vertex_max
    j_node1 = elm2%vertex(i)
    j_node2 = elm2%vertex(mod(i,4)+1)
    if (     ( (j_node1 .eq. i_node1) .and. (j_node2 .eq. i_node2) ) &
        .or. ( (j_node1 .eq. i_node2) .and. (j_node2 .eq. i_node1) ) ) then
      j_side = i
      exit
    endif
  enddo
  
  return

end subroutine are_elements_neighbours
