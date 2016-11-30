!> Find the element and position at this psi, with poloidal angle theta.
!> Root finding on the line from R_axis, Z_axis, parametrised by u.
!> R = R_axis + u cos(theta); Z = Z_axis + u sin(theta).
!> The objective function is (psi - psi(u)) where psi(u) = psi(st(RZ(u)))
!> and we calculate st with find_RZ_nearby and RZ as above.
!> Use Newton's method to find the root. This code could be upgraded to use a
!> higher-order method, as we also have the second derivatives.
!>
!> A reasonable initial guess is to assume psi to be linear in C*u/a
subroutine find_theta_psi(node_list,element_list,theta,psi,phi,u_guess,R_axis,Z_axis,i_elm,s,t,R,Z)
use constants
use data_structure
implicit none

! --- Routine parameters
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
real*8,                   intent(in)    :: theta
real*8,                   intent(in)    :: psi, phi
real*8,                   intent(in)    :: u_guess !< an estimated initial position
real*8,                   intent(in)    :: R_axis
real*8,                   intent(in)    :: Z_axis
integer,                  intent(out)   :: i_elm
real*8,                   intent(out)   :: s, t, R, Z

! --- Internal variables
real*8  :: u, du, R_try, Z_try, s_out, t_out, err, st_out(2)
real*8, dimension(1) :: P, P_s, P_t, P_phi
real*8  :: R_s, R_t, Z_s, Z_t
real*8  :: inv_st_jac, inv_u_jac, psi_R, psi_Z, psi_u
integer :: i_elm_out, ifail, backtrack_step, newton_iter_number
integer, parameter :: num_backtrack_steps = 6
real*8, parameter  :: backtrack_factor = 0.8d0
integer, parameter :: newton_iter_max = 20
integer, parameter :: i_var(1) = [1], n_ivar=1
real*8, parameter  :: tol = 1d-4

! 0. Preparation for the loop
! Keep trying with smaller u until we find an element
u = u_guess
do backtrack_step = 0, num_backtrack_steps
  R_try = R_axis + u*cos(theta)
  Z_try = Z_axis + u*sin(theta)
  call find_RZ(node_list,element_list,R_try,Z_try,R,Z,i_elm_out,s_out,t_out,ifail)
  if (ifail .eq. 0) then
    exit
  else
    u = u*backtrack_factor
  end if
end do
inv_u_jac = 1.d0/(sin(theta)*cos(theta))

! 1. Iterate to find the right value of psi(u)
do newton_iter_number = 1, newton_iter_max
  ! Calculate psi and the derivative in u direction
  call interp_PRZ(node_list,element_list,i_elm_out,i_var,n_ivar,s_out,t_out,phi, &
      P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
  inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
  psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
  psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
  psi_u    = (psi_R*sin(theta) + psi_Z*cos(theta)) * inv_u_jac

  ! Swap variables
  i_elm = i_elm_out
  s = s_out
  t = t_out
  ! Calculate the error, exit if it is low enough
  err = (P(1) - psi)
  if (abs(err) .lt. tol) return

  ! Formulate a new guess (u) based on psi and psi' 
  du = -err/psi_u
  write(*,"(12g14.6)") P(1)-psi, u, du, i_elm_out, s, t
  ! Backtrack if we go out of the domain
  ! Input here:
  ! u_guess, i_elm, s, t (consistent with eachother)
  do backtrack_step = 0, num_backtrack_steps
    ! Calculate R and Z for this guess
    R_try = R_axis + (u+du)*cos(theta)
    Z_try = Z_axis + (u+du)*sin(theta)
    ! Find st based on the old value and the size of this step
    ! strange argument order: new old (RZ), old new (st), old new (i_elm)
    call find_RZ_nearby(node_list, element_list, [R_try,Z_try], [R,Z], &
                    [s,t], st_out, i_elm, i_elm_out, ifail)
                ! watch out for array temporary
    s_out = st_out(1)
    t_out = st_out(2)
    if (ifail .ne. 0) then
      du = du*backtrack_factor
    else
      u = u+du
      exit ! the backtrack loop
    end if
  end do
  if (ifail .ne. 0 .or. isnan(P(1)) .or. abs(P(1)) .lt. tol) then
    write(*,"(A,g12.6,A,g12.6,A,3g12.6,A,g12.6,A)") "Position not found for psi=", psi, ", theta=",theta, " last guess=", u, R_try, Z_try, " (psi=",P(1),")"
    i_elm = 0
    return
  end if
end do

! 2. Test if finding was successful
if (newton_iter_number .gt. newton_iter_max) then
  ! Indicate this by setting i_elm
  write(*,*) "Too many iterations"
  i_elm = 0
end if
end subroutine find_theta_psi
