pure logical function neighbours(node_list,elm1,elm2,inb1,inb2)
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

type(type_node_list), intent(in) :: node_list
type(type_element), intent(in)   :: elm1, elm2
integer, intent(out), pointer    :: inb1, inb2

integer             :: iv(4,2), i, j, nb
real*8, parameter   :: tol = 1.d-6 !< Distance in meters between nodes to be considered at the same position
real*8 :: dist

neighbours = .false.
nb = 0
do i=1,4
  do j=1,4
    dist = norm2(node_list%node(elm1%vertex(i))%x(1,1:2) - node_list%node(elm2%vertex(j))%x(1,1:2))
    ! If the nodes are the same or are at the same position
    if (elm1%vertex(i) .eq. elm2%vertex(j) .or. dist .lt. tol) then
      nb = nb + 1
      iv(nb,1) = i
      iv(nb,2) = j
    endif
  enddo
enddo
if (nb .gt. 1 ) then
  neighbours=.true.
  inb1 = minval(iv(1:nb,1)) ; if ( abs(iv(1,1)-iv(2,1)) .gt. 1 ) inb1 = 4
  inb2 = minval(iv(1:nb,2)) ; if ( abs(iv(1,2)-iv(2,2)) .gt. 1 ) inb2 = 4
endif
end function
