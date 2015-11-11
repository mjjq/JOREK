!> Optimized subroutine to find st coordinates corresponding to x_new=[R_new, Z_new] using
!! The previous values x_old=[R_old, Z_old], st_old(2), i_elm_old and checking adjacent elements first
!! It first checks the current element, using newton iteration to find st_try corresponding to x_new
!! Once this is found it calls check_element_boundary to see if we have crossed an element boundary
!! If this is the case, the search is restarted within the new element.
!! This procedure continues until the required tolerance is reached, or element_try_max is reached.
subroutine find_RZ_nearby(node_list, element_list, x_new, x_old, st_old, st_new, i_elm_old, i_elm_new, ifail)
use data_structure
implicit none
!> Input parameters
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
real*8,                   intent(in)    :: x_new(2), x_old(2), st_old(2)
integer,                  intent(in)    :: i_elm_old
real*8,                   intent(out)   :: st_new(2)
integer,                  intent(out)   :: i_elm_new
integer,                  intent(out)   :: ifail !< if ifail = -1 the position could not be found in the grid

!> Accuracy defaults (tolerances are squared!)
real*8, parameter ::  in_element_tolerance = 1.d-12
real*8, parameter :: out_element_tolerance = 1.d-9
integer :: newton_iter_max = 8
integer :: element_try_max = 6 ! Try newton iterations in at most this many elements
! A loop can occur here!
integer :: num_backtrack_steps = 2 ! Try 0.5**this times the step at a minimum

!> Internal variables
integer :: newton_iter_number
integer :: element_try_index
integer :: backtrack_step
real*8 :: st_jac_det, R_s, R_t, Z_s, Z_t
real*8 :: st_step(2), st_try(2), x_step(2) ! x_step = (R,Z)
real*8 :: err2, err2_old
logical :: changed, lost, search

! Setup initial values
x_step = x_old
i_elm_new = i_elm_old
st_try = st_old
! Find the jacobian at the current s and t position
call try_interp(node_list,element_list,i_elm_new,st_try,x_step,R_s,R_t,Z_s,Z_t,st_jac_det)
err2_old = dot_product(x_step-x_new,x_step-x_new)




! Outer loop tries newton iteration, changes element if necessary
do element_try_index = 1, element_try_max
  ! Inner loop, newton iteration to find s and t in or out of this element
  do newton_iter_number = 1, newton_iter_max
    ! Perform newton iteration by calculating the inverse of the jacobian matrix

    ! Calculate the trial newton step
    st_step(1) = ( Z_t * (x_new(1)-x_step(1)) - R_t * (x_new(2)-x_step(2))) / st_jac_det
    st_step(2) = (-Z_s * (x_new(1)-x_step(1)) + R_s * (x_new(2)-x_step(2))) / st_jac_det

    ! Test different step sizes (backtracking)
    do backtrack_step = 0, num_backtrack_steps
      st_try = (st_new + 0.5**backtrack_step * st_step)

      ! Calculate the x_step corresponding to st_try with the coordinates of element i_elm_new
      call try_interp(node_list,element_list,i_elm_new,st_try,x_step,R_s,R_t,Z_s,Z_t,st_jac_det)
      err2 = dot_product(x_step-x_new,x_step-x_new)
      if (err2 .lt. err2_old) exit
    enddo
    ! Save this step value
    st_new = st_try
    err2_old = err2
    if (isnan(err2)) exit

    if (st_try(1) > 1.d0 .or. st_try(1) < 0.d0 .or. st_try(2) > 1.d0 .or. st_try(2) < 0.d0) then
      if (err2 < out_element_tolerance) exit
    else
      if (err2 < in_element_tolerance) exit
    endif
  enddo

  ! Guards against errors above
  if (isnan(err2)) then
    write(*,*) "DEBUG: NaN encountered after newton iteration, using find_RZ"
    call find_RZ(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,st_new(1),st_new(2),ifail)
    exit
  endif
  if (newton_iter_number .gt. newton_iter_max) then
    write(*,"(A,i4,A,i3,A,2g14.6)") "WARNING: iteration for st did not converge after", newton_iter_max, " tries in element ", i_elm_new, &
    " using find_RZ", x_new
    call find_RZ(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,st_new(1),st_new(2),ifail)
    exit
  endif

  ! Now check to see if we have left the current element
  call check_element_boundary(element_list,i_elm_new,st_try,i_elm_new,st_try,changed,lost,search)

  if (search) then
    ! Use the extremely slow option
    call find_RZ(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,st_new(1),st_new(2),ifail)
    write(*,"(A,i4,A,2g14.6)") "WARNING: check_element_boundary returned search, used find_RZ in element", i_elm_new, ", at position", x_new
    exit ! We have used the nuclear option, stop the loop now
  endif

  if (lost) then
    write(*,*) "WARNING: position ", x_new, "not found in grid"
    ifail = -1
    exit
  endif

  if (.not. changed) then
    ! We have converged!
    exit
  else
    if (element_try_index .eq. element_try_max) then
      write(*,"(A,i4)") "WARNING: insufficient iterations for element change, trying brute force method. start at element", i_elm_new
      call find_RZ(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,st_new(1),st_new(2),ifail)
    else
      ! We have changed element, recalculate st_jac and err2
      call try_interp(node_list,element_list,i_elm_new,st_try,x_step,R_s,R_t,Z_s,Z_t,st_jac_det)
      err2_old = dot_product(x_step-x_new,x_step-x_new)
    endif
  endif
enddo
end subroutine find_RZ_nearby


!> Auxiliary subroutine for find_RZ_nearby
subroutine try_interp(node_list,element_list,i_elm,st,x,R_s,R_t,Z_s,Z_t,st_jac_det)
use data_structure
implicit none
!> Input parameters
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
real*8,                   intent(in)    :: st(2)
integer,                  intent(in)    :: i_elm
real*8,                   intent(out)   :: x(2), R_s, R_t, Z_s, Z_t, st_jac_det

real*8, parameter :: st_jac_det_min = 1.d-6

call interp3_RZ(node_list,element_list,i_elm,st(1),st(2),x(1),R_s,R_t,x(2),Z_s,Z_t)
st_jac_det = R_s * Z_t - R_t * Z_s
if (st_jac_det**2 .lt. st_jac_det_min**2) st_jac_det = sign(st_jac_det_min, st_jac_det)
end subroutine try_interp
