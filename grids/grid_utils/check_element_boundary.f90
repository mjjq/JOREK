!> Check if a position is in or outside an element, and if so on which side
!!
!! Uses the entire list of elements and the current element, a new and old position
!! to find the number of the new element j_elm, if it is 
!! lost (outside of the domain) or we need to search for it
subroutine check_element_boundary(element_list,i_elm,x,j_elm,y,changed,lost,search)
use data_structure
implicit none

! -- Input arguments
type (type_element_list), intent(in) :: element_list !< The list of all finite elements
integer, intent(in)                  :: i_elm !< The index of the current element
real*8, intent(in)                   :: x(2) !< The coordinates s and t in the coordinate system of i_elm

! -- Output arguments
real*8, intent(out)                  :: y(2) !< The best guess 
integer, intent(out)                 :: j_elm !< The index of our next best guess of the element containing x
logical, intent(out)                 :: changed !< Whether we changed element and need to recalculate the coordinates
logical, intent(out)                 :: lost !< True if the particle is not in any element we know of
logical, intent(out)                 :: search !< True if there is no neighbour in the direction we found

! -- local variables
integer :: i_side

changed = .false.
lost    = .false.
search  = .false.

if (maxval(abs(x(1:2)-0.5)) .le. 0.5) return ! If s and t are in [0,1] the position is in this element

i_side  = 0.

! Find the direction the position is in from here
if (x(1) .gt. 1.0) then
  i_side = 2
elseif (x(1) .lt. 0.0) then
  i_side = 4
elseif (x(2) .gt. 1.0) then
  i_side = 3
elseif (x(2) .lt. 0.0) then
  i_side = 1
else
  return
endif

! If that does not exist we should search
if (element_list%element(i_elm)%neighbours(i_side) .lt. 0) then
  search = .true.
  return
endif

! If it does exist look up the value of j_elm
if (i_side .gt. 0) then
  j_elm = element_list%element(i_elm)%neighbours(i_side)
  if (j_elm .gt. 0) then
    ! Transform the variable into the new basis functions of that element (first-order approximation)
    y(1)  = element_list%element(i_elm)%transform(i_side,1,1) + element_list%element(i_elm)%transform(i_side,1,2) * x(1) + element_list%element(i_elm)%transform(i_side,1,3) * x(2)
    y(2)  = element_list%element(i_elm)%transform(i_side,2,1) + element_list%element(i_elm)%transform(i_side,2,2) * x(1) + element_list%element(i_elm)%transform(i_side,2,3) * x(2)

    changed = .true.
  else
    lost = .true.
 !   write(*,*) 'particle is lost!'
  endif
endif

return
end
