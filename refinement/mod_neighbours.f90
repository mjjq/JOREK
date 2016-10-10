module mod_neighbours
contains
logical function neighbours(elm1,elm2,inb1,inb2)
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

type (type_element) :: elm1, elm2
integer             :: inb1, inb2, iv(4,2), i, j, nb

neighbours = .false.
nb = 0
! write(*,'(A13,4i6,A2,4i6)') ' NEIGHBOURS? ',elm1%vertex,'  ',elm2%vertex
do i=1,4
  do j=1,4
    if (elm1%vertex(i)==elm2%vertex(j)) then
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
else
endif
return
end
end module mod_neighbours
