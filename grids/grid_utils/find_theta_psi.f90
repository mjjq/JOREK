!> Find the element and position at this psi, with poloidal angle theta.
!> Root finding on the line from R_axis, Z_axis, parametrised by u.
!> R = R_axis + u cos(theta); Z = Z_axis + u sin(theta).
!> The objective function is (psi - psi(u)) where psi(u) = psi(st(RZ(u)))
!> and we calculate st with find_RZ_nearby and RZ as above.
!> Use Newton's method to find the root. This code could be upgraded to use a
!> higher-order method, as we also have the second derivatives.
subroutine find_theta_psi(node_list,element_list,psi_minmax,theta,psi,phi,R_axis,Z_axis,i_elm,s,t,R,Z)
use constants
use data_structure
implicit none

! --- Routine parameters
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
real*8, dimension(element_list%n_elements,2), intent(in) :: psi_minmax !< list of minima and maxima of psi in these elements
real*8,                   intent(in)    :: theta
real*8,                   intent(in)    :: psi, phi
real*8,                   intent(in)    :: R_axis
real*8,                   intent(in)    :: Z_axis
integer,                  intent(out)   :: i_elm
real*8,                   intent(out)   :: s, t, R, Z

! --- Internal variables
real*8  :: u, u_init(3), P_init(3), du, R_try, Z_try, s_out, t_out, err, st_out(2)
real*8, dimension(1) :: P, P_s, P_t, P_phi
real*8  :: R_s, R_t, Z_s, Z_t, theta_normalized
real*8  :: inv_st_jac, psi_R, psi_Z, psi_u
integer :: i_elm_out, ifail, backtrack_step, newton_iter_number
integer, parameter :: num_backtrack_steps = 10 ! to prevent going over the border of the domain
real*8, parameter  :: backtrack_factor = 0.99d0
integer, parameter :: newton_iter_max = 6
integer, parameter :: i_var(1) = [1], n_ivar=1
real*8, parameter  :: tol = 1d-8

integer :: i, j
logical :: out_of_domain
logical, dimension(element_list%n_elements) :: psi_right
real*8, dimension(n_vertex_max) :: angle
real*8, dimension(n_dim) :: x
real*8 :: amax, amin

! 0. Preparation
! Normalize theta to -pi,pi
theta_normalized = modulo(theta+PI,TWOPI)-PI
out_of_domain = .false.

! 1. Find the right elements (all elements where the theta and psi lines enter)
psi_right = (psi_minmax(:,1) .lt. psi) .and. (psi_minmax(:,2) .gt. psi)
do i=1,element_list%n_elements
  if (psi_right(i)) then
    ! Calculate angles
    do j=1,n_vertex_max
      x = node_list%node(element_list%element(i)%vertex(j))%x(1,:)
      angle(j) = atan2(x(2)-Z_axis,x(1)-R_axis)
    end do
    ! workaround for the one element at the boundary
    amax = maxval(angle)
    amin = minval(angle)
    ! separate case if the discontinuity (at theta=-pi) is in this element
    if (amax .gt. PI*0.5d0 .and. amin .lt. -PI*0.5d0) then
      ! We need now to compare the minimum on the positive side and the maximum on the negative side
      if (minval(angle,angle .gt. PI*0.5d0) .lt. theta .and. maxval(angle,angle .lt. -PI*0.5d0)+TWOPI .gt. theta) then
        i_elm = i
      end if
    else
      if (amax .gt. theta_normalized .and. amin .lt. theta_normalized) i_elm = i
    end if
  end if
end do
if (i_elm .eq. 0) then
  write(*,"(A,g12.6,A,g12.6,A)") "WARNING: no suitable elements found, skipping for psi=", psi, " theta=", theta/PI, "pi"
  return
end if

  
  

! 2. Preparation for the loop
! Calculate three different guesses for u
do i=1,3
  call interp_PRZ(node_list, element_list,i_elm,i_var,n_ivar,&
      0.25d0*real(i,8),0.5d0,phi,& ! this order works best for most elements
      P,P_s,P_t,P_phi,&
      R,R_s,R_t,Z,Z_s,Z_t)
  u_init(i) = norm2([R-R_axis,Z-Z_axis]) ! approximation
  angle(i)  = atan2(Z-Z_axis,R-R_axis) ! reuse angle even though the size differs, careful
  P_init(i) = P(1)
end do
i = minloc(abs(P_init-psi),1)
u = u_init(i)/cos(angle(i)-theta)
write(*,"(6g16.8)") P_init-psi, u_init
! Keep trying with smaller u until we find an element
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
if (backtrack_step .gt. num_backtrack_steps) then
  write(*,*) "Cannot find initial position after ", backtrack_step
  return
end if

! 3. Iterate to find the right value of psi(u)
do newton_iter_number = 1, newton_iter_max
  ! Calculate psi and the derivative in u direction
  call interp_PRZ(node_list,element_list,i_elm_out,i_var,n_ivar,s_out,t_out,phi, &
      P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
  inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
  psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
  psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
  psi_u    = (psi_R*cos(theta) + psi_Z*sin(theta))

  ! Swap variables
  i_elm = i_elm_out
  s = s_out
  t = t_out
  ! Calculate the error, exit if it is low enough
  err = (P(1) - psi)
  if (abs(err) .lt. tol) return

  ! Formulate a new guess (u) based on psi and psi' 
  du = -err/psi_u
  write(*,"(12g14.6)") P(1)-psi, u, du, i_elm_out, s, t, out_of_domain
  if (out_of_domain .and. du .gt. 0.d0) du = -du ! u must go down in this case
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
  u = u+du
  if (ifail .ne. 0) then
    out_of_domain = .true.
    i_elm_out = i_elm ! keep going in the fields of the old element
  else
    out_of_domain = .false.
  end if
  if (isnan(P(1)) .or. abs(P(1)) .lt. tol) then
    write(*,"(A,g12.6,A,g12.6,A,3g12.6,A,g12.6,A)") "Position not found for psi=", psi, ", theta=",theta, " last guess=", u, R_try, Z_try, " (psi=",P(1),")"
    i_elm = 0
    return
  end if
end do

! 4. Warn if failed
if (newton_iter_number .gt. newton_iter_max) then
  ! Indicate this by setting i_elm
  write(*,*) "Too many iterations, skipping"
  i_elm = 0
end if
end subroutine find_theta_psi
