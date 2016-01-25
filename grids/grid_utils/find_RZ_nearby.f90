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
real*8,                   intent(in)    :: x_new(2) !< The new R,Z location
real*8,                   intent(in)    :: x_old(2) !< The old R,Z location
real*8,                   intent(in)    :: st_old(2) !< The old st location (used to compute a guess)
integer,                  intent(in)    :: i_elm_old
real*8,                   intent(out)   :: st_new(2) !< The found new coordinates
integer,                  intent(out)   :: i_elm_new
integer,                  intent(out)   :: ifail !< if ifail = -1 the position could not be found in the grid

!> Accuracy defaults (tolerances are squared!, units of element size)
real*8, parameter ::  element_tolerance = 1.d-12
integer :: newton_iter_max = 4
integer :: element_try_max = 2 ! Try newton iterations in at most this many elements
integer :: num_backtrack_steps = 2 ! Try 0.5**this times the step at a minimum

!> Internal variables
integer :: newton_iter_number
integer :: element_try_index
integer :: backtrack_step
real*8 :: inv_st_jac_det, R_s, R_t, Z_s, Z_t
real*8 :: st_step(2), st_try(2), x_step(2) ! x_step = (R,Z) of trial position
real*8 :: err2, err2_old

!> For output of check_element_boundary
integer(1) :: stat

! Setup initial values
x_step = x_old ! start at the current position
i_elm_new = i_elm_old ! start in the current element
st_try = st_old ! start at the old position
! Find the jacobian at the current s and t position
call try_interp(node_list,element_list,i_elm_new,st_try,x_step,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
err2 = dot_product(x_step-x_new,x_step-x_new)
ifail=0


! Outer loop tries newton iteration, changes element if necessary
do element_try_index = 1, element_try_max
  ! Inner loop, newton iteration to find s and t in or out of this element
  do newton_iter_number = 1, newton_iter_max
    ! Perform newton iteration by calculating the inverse of the jacobian matrix explicitly
    err2_old = err2

    ! Calculate the trial newton step
    st_step(1) = ( Z_t * (x_new(1)-x_step(1)) - R_t * (x_new(2)-x_step(2))) * inv_st_jac_det
    st_step(2) = (-Z_s * (x_new(1)-x_step(1)) + R_s * (x_new(2)-x_step(2))) * inv_st_jac_det


    ! Test different step sizes (backtracking)
    do backtrack_step = 0, num_backtrack_steps
      st_try = (st_new + 0.5**backtrack_step * st_step)

      ! Calculate the x_step corresponding to st_try with the coordinates of element i_elm_new
      call try_interp(node_list,element_list,i_elm_new,st_try,x_step,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
      err2 = dot_product(x_step-x_new,x_step-x_new)
      if (err2 .lt. err2_old) exit
    enddo
    ! Save this trial value
    st_new = st_try

    if (err2 < element_tolerance) exit
  enddo


  if (isnan(err2)) then
    write(*,*) "WARNING: NaN encountered after newton iteration, using find_RZ"
    call find_RZ(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,st_new(1),st_new(2),ifail)
    ifail=2
    return
  endif
  if (newton_iter_number .gt. newton_iter_max) then
    write(*,"(A,i4,A,i5,A,2g14.6,A,3g14.6)") "WARNING: iteration for st did not converge after", newton_iter_max, " tries in element ", i_elm_new, &
    " using find_RZ", x_new, "err2(old)/convergence: ", err2, err2_old, err2_old/err2
      write(*,"(A,2g16.8)") "Find_RZ at ", x_new
    call find_RZ(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,st_new(1),st_new(2),ifail)
    ifail=3
    return
  endif

  ! Now check to see if we have left the current element
  call check_element_boundary(element_list,i_elm_new,st_try,i_elm_new,st_try,stat)

  select case (stat)
  case (1) ! CHANGED
    if (element_try_index .eq. element_try_max) then ! do not do this if it is the last round, we will try find_RZ below
      write(*,"(A,i5)") "WARNING: insufficient iterations for element change, trying brute force method. start at element", i_elm_new
      write(*,"(A,2g16.8)") "Find_RZ at ", x_new
      call find_RZ(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,st_new(1),st_new(2),ifail)
      ifail = 5
      return
    else
      ! We have changed element, recalculate st_jac and err2
      call try_interp(node_list,element_list,i_elm_new,st_try,x_step,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
      !err2 = dot_product(x_step-x_new,x_step-x_new) ! XXX this should be set, logically, but performs better if not
    endif
  case (2) ! LOST
    write(*,*) "WARNING: position ", x_new, "not found in grid"
    ifail = -1
    return
  case (3) ! SEARCH
    call find_RZ(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,st_new(1),st_new(2),ifail)
    write(*,"(A,i5,A,2g14.6)") "WARNING: check_element_boundary returned search, used find_RZ in element", i_elm_new, ", at position", x_new
    ifail=4
    return ! Stop because find_RZ works always
  case default ! SAME == 0
    return ! we are converged anyhow, no problem
  end select

enddo
end subroutine find_RZ_nearby


!> Auxiliary subroutine for find_RZ_nearby
subroutine try_interp(node_list,element_list,i_elm,st,x,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
use data_structure
implicit none
!> Input parameters
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
real*8,                   intent(in)    :: st(2)
integer,                  intent(in)    :: i_elm
real*8,                   intent(out)   :: x(2), R_s, R_t, Z_s, Z_t, inv_st_jac_det

real*8, parameter :: inv_st_jac_det_max = 1.d6

call interp3_RZ(node_list,element_list,i_elm,st(1),st(2),x(1),R_s,R_t,x(2),Z_s,Z_t)
inv_st_jac_det = 1.d0/(R_s * Z_t - R_t * Z_s)
if (inv_st_jac_det**2 .gt. inv_st_jac_det_max**2) inv_st_jac_det = sign(inv_st_jac_det_max, inv_st_jac_det)
end subroutine try_interp
