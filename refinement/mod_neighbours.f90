module mod_neighbours
contains
logical function neighbours(node_list,elm1,elm2,inb1,inb2)
!-------------------------------------------------------
! function to check if the two elements elm1 and elm2
! are neighbours. Two elements are neighbours if they
! share a side, i.e. if two nodes are the same.
! inb1 -> the index of the shared neighbour of elm1
! inb2 -> the index of the shared neighbour of elm2
!   note : does not work for unequal sized neighbours
!-------------------------------------------------------
use data_structure
implicit none

type (type_node_list), intent(in) :: node_list
type (type_element), intent(in)   :: elm1, elm2
integer, intent(out)              :: inb1, inb2
real*8, parameter :: tol=1.d-8

integer :: n1(2), n2(2) ! nodes of side i of elm1 or of side j of elm2
integer :: i, j
neighbours = .false.

! First test by node number
do i=1,4 ! Loop over sides of element 1
  do j=1,4 ! Loop over sides of element 2
    ! node numbers are related to sides as node1=mod(side-1,4)+1, node2=mod(side,4)+1
    n1 = elm1%vertex(mod([i-1,i],4)+1)
    n2 = elm2%vertex(mod([j-1,j],4)+1)
    ! If match cross or straight (i.e. 1->2/2->1 or 1->1/2->2)
    ! Four different cases
    if      (n1(1) .eq. n2(2)) then
      neighbours = (n1(2) .eq. n2(1)) .or. &
        (norm2(node_list%node(n1(2))%x(1,1:2)-node_list%node(n2(1))%x(1,1:2)) .lt. tol)
    else if (n1(2) .eq. n2(1)) then
      neighbours = (n1(1) .eq. n2(2)) .or. &
        (norm2(node_list%node(n1(1))%x(1,1:2)-node_list%node(n2(2))%x(1,1:2)) .lt. tol)
    else if (n1(1) .eq. n2(1)) then
      neighbours = (n1(2) .eq. n2(2)) .or. &
        (norm2(node_list%node(n1(2))%x(1,1:2)-node_list%node(n2(2))%x(1,1:2)) .lt. tol)
    else if (n1(2) .eq. n2(2)) then
      neighbours = (n1(1) .eq. n2(1)) .or. &
        (norm2(node_list%node(n1(1))%x(1,1:2)-node_list%node(n2(1))%x(1,1:2)) .lt. tol)
    end if
    if (neighbours) then
      inb1 = i
      inb2 = j
      return
    end if
  enddo
enddo
end function neighbours
end module mod_neighbours
