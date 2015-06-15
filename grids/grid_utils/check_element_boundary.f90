subroutine check_element_boundary(element_list,i_elm,x,x_prev,j_elm,y,delta_x,changed,lost,search)
use data_structure
implicit none

type (type_element_list) :: element_list
integer                  :: i_elm, j_elm, i_side
real*8                   :: x(2), x_prev(2), y(2), delta_x(2), x_in(2)
logical                  :: changed, lost, search

changed = .false.
lost    = .false.
search  = .false.

if (maxval(abs(x(1:2)-0.5)) .le. 0.5) return

delta_x = 0.
i_side  = 0.

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

if (element_list%element(i_elm)%neighbours(i_side) .lt. 0) then
  search = .true.
  return
endif

if (i_side .gt. 0) then

  j_elm = element_list%element(i_elm)%neighbours(i_side)

  if (j_elm .gt. 0) then

    y(1)  = element_list%element(i_elm)%transform(i_side,1,1) + element_list%element(i_elm)%transform(i_side,1,2) * x(1) + element_list%element(i_elm)%transform(i_side,1,3) * x(2)
    y(2)  = element_list%element(i_elm)%transform(i_side,2,1) + element_list%element(i_elm)%transform(i_side,2,2) * x(1) + element_list%element(i_elm)%transform(i_side,2,3) * x(2)

    i_elm = j_elm
    x     = y

    changed = .true.

  else

    lost = .true.
 !   write(*,*) 'particle is lost!'

  endif

endif

return
end