recursive subroutine find_RZ(node_list,element_list,R_find,Z_find,R_out,Z_out,ielm_out,s_out,t_out,ifail)
!-------------------------------------------------------------------------
! solves two non-linear equations using Newtons method (from numerical recipes)
! LU decomposition replaced by explicit solution of 2x2 matrix.
!
! finds the crossing of two coordinate lines given as a series of cubics
!-------------------------------------------------------------------------
use data_structure

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list

integer :: i, k, ifail, ntrial, istart
integer :: ielm_out
real*8  :: R_find, Z_find, R_out,Z_out,s_out,t_out, Rmin, Rmax, Zmin, Zmax
real*8  :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8  :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8  :: tolx, tolf, errx, errf, temp, dis
real*8  :: x(2), FVEC(2), FJAC(2,2), p(2)

logical              :: minmax_initialised
real*8, allocatable  :: elements_minmax(:,:)

save minmax_initialised, elements_minmax

ielm_out = 0

if (.not. allocated(elements_minmax)) then

  allocate(elements_minmax(4,element_list%n_elements))
  minmax_initialised = .false.
!  write(*,*) ' *** FIND_RZ : initialising ***'

elseif (size(elements_minmax,2) .ne. element_list%n_elements) then

  deallocate(elements_minmax)
  allocate(elements_minmax(4,element_list%n_elements))
  minmax_initialised = .false.
!  write(*,*) ' *** FIND_RZ : re-initialising ***'

endif

if (.not. minmax_initialised) then

  do k=1,element_list%n_elements

    call RZ_minmax(node_list,element_list,k,Rmin,Rmax,Zmin,Zmax)

    elements_minmax(:,k) = (/ Rmin, Rmax, Zmin, Zmax /)

  enddo

  minmax_initialised = .true.
!  write(*,*) ' *** FIND_RZ : initialised ***'

endif


do k=1,element_list%n_elements

!  call RZ_minmax(node_list,element_list,k,Rmin,Rmax,Zmin,Zmax)

  Rmin =  elements_minmax(1,k)
  Rmax =  elements_minmax(2,k)
  Zmin =  elements_minmax(3,k)
  Zmax =  elements_minmax(4,k)

!  Rmin = Rmin - 0.05
!  Rmax = Rmax + 0.05
!  Zmin = Zmin - 0.05
!  Zmax = Zmax + 0.05

  if ((R_find .ge. Rmin) .and. (R_find .le. Rmax) .and. &
      (Z_find .ge. Zmin) .and. (Z_find .le. Zmax) ) then

    ntrial = 20
    tolx = 1.d-8
    tolf = 1.d-15

    do istart = 1,5

      if (istart .eq. 1) then
        x(1) = 0.5d0
        x(2) = 0.5d0
      elseif (istart .eq. 2) then
        x(1) = 0.75d0
        x(2) = 0.75d0
      elseif (istart .eq. 3) then
        x(1) = 0.75d0
        x(2) = 0.25d0
      elseif (istart .eq. 4) then
        x(1) = 0.25d0
        x(2) = 0.75d0
      elseif (istart .eq. 5) then
        x(1) = 0.25d0
        x(2) = 0.25d0
      endif

      ifail = 999

      do i=1,ntrial

        call interp_RZ(node_list,element_list,k,x(1),x(2),RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                        ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

        FVEC(1)   = RRg1 - R_find
        FVEC(2)   = ZZg1 - Z_find
        FJAC(1,1) = dRRg1_dr
        FJAC(1,2) = dRRg1_ds
        FJAC(2,1) = dZZg1_dr
        FJAC(2,2) = dZZg1_ds

        errf=abs(fvec(1))+abs(fvec(2))

!      write(*,'(A,i3,8e16.8)') ' newton   : ',i,errf,errx,x,RRg1,R_find,ZZg1,Z_find
!      write(*,'(A,i3,8e16.8)') ' newton   : ',i,dRRg1_dr,dRRg1_ds,dZZg1_dr,dZZg1_ds

        if (errf .le. tolf) then

          s_out     = x(1)
          t_out     = x(2)

          ielm_out  = k
          R_out     = RRg1
          Z_out     = ZZg1

!        write(*,'(A,i3,4e16.8)') ' newton (1) : ',i,errf,errx,x

          ifail = 0
          return
        endif

        p = -fvec

        temp = p(1)
        dis  = fjac(2,2)*fjac(1,1)-fjac(1,2)*fjac(2,1)

        if (dis .ne. 0.d0) then
          p(1) = (fjac(2,2)*p(1)-fjac(1,2)*p(2))/dis
          p(2) = (fjac(1,1)*p(2)-fjac(2,1)*temp)/dis
        else
          exit
        endif

        errx=abs(p(1)) + abs(p(2))

        p = min(p,+0.25d0)
        p = max(p,-0.25d0)

        x = x + p

        x = max(x,+0.d0)
        x = min(x,+1.d0)

        if (errx .le. tolx) then

          s_out     = x(1)
          t_out     = x(2)

          ielm_out  = k
          R_out     = RRg1
          Z_out     = ZZg1

!        write(*,'(A,i3,4e16.8)') ' newton (2) : ',i,errf,errx,x

          ifail = 0
          return
        endif

      enddo
    enddo

  endif

enddo

if (ielm_out .eq. 0) ifail = 99

!write(*,'(A,8e16.8)') ' find_RZ wrong exit ',x,errx,errf

return
end subroutine find_RZ





!> Optimized version of find_RZ.
recursive subroutine find_RZ2(node_list, element_list, R_find, Z_find, R_out, Z_out, ielm_out,      &
  s_out, t_out, ifail)

use data_structure

implicit none

! --- Constants
integer, parameter :: niter     = 20                   !< Maximum number of Newton iterations
real*8,  parameter :: tolf      = 1.d-6                !< Tolerance for spatial distance
real*8,  parameter :: tolx      = 1.d-15               !< Tolerance for iteration step width
real*8,  parameter :: delta     = 0.05d0               !< Maximum number of Newton iterations

! --- Routine parameters
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
real*8,                   intent(in)    :: R_find, Z_find
real*8,                   intent(out)   :: R_out, Z_out, s_out, t_out
integer,                  intent(out)   :: ielm_out, ifail

integer :: i, j, k, iv, istart
real*8  :: Rmin, Rmax, Zmin, Zmax, temp, dis, RR, RR_s, RR_t, ZZ, ZZ_s, ZZ_t, x(2), fvec(2), p(2)

ielm_out = 0
ifail    = 99

L_EL: do k = 1, element_list%n_elements

  call RZ_minmax(node_list, element_list, k, Rmin, Rmax, Zmin, Zmax) ! <<< most expensive call!!!
  
  if ( (R_find > Rmin - delta) .and. (R_find < Rmax + delta) .and. (Z_find > Zmin - delta) .and.   &
     (Z_find < Zmax + delta) ) then ! (If the element could be relevant, proceed:)
    
    L_ST: do istart = 1, 5
      
      ! Try up to five different starting positions inside the element:
      if (istart == 1) then
        x(:) = (/ 0.50d0, 0.50d0 /)
      else if (istart == 2) then
        x(:) = (/ 0.23d0, 0.23d0 /)
      else if (istart == 3) then
        x(:) = (/ 0.77d0, 0.77d0 /)
      else if (istart == 4) then
        x(:) = (/ 0.77d0, 0.23d0 /)
      else if (istart == 5) then
        x(:) = (/ 0.23d0, 0.77d0 /)
      end if
      
      do i = 1, niter
        
        call interp_RZ2(node_list, element_list, k, x(1), x(2), RR, RR_s, RR_t, ZZ, ZZ_s, ZZ_t)
        
        fvec(:) = (/ RR - R_find, ZZ - Z_find /)
        
        if (sqrt(sum(fvec**2)) <= tolf) then
          ielm_out = k
          exit L_EL
        endif
        
        dis  = ZZ_t * RR_s - RR_t * ZZ_s
        if (dis == 0.d0) exit L_ST
        
        p(:) = (/ RR_t * fvec(2) - ZZ_t * fvec(1), ZZ_s * fvec(1) - RR_s * fvec(2) /) / dis
        
        p(:) = max( min(p(:),     +0.25d0), -0.25d0 ) ! (limit iteration step size)
        x(:) = max( min(x(:)+p(:),+1.00d0), -0.00d0 ) ! (restict s and t to valid range)
        
        if (sqrt(sum(p**2)) <= tolx) then
          ielm_out  = k
          exit L_EL
        end if
        
      end do
    end do L_ST
    
  end if
  
end do L_EL

if ( ielm_out /= 0 ) then
  s_out     = x(1)
  t_out     = x(2)
  R_out     = RR
  Z_out     = ZZ
  ifail     = 0
end if

end subroutine find_RZ2


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
real*8, parameter :: out_element_tolerance = 1.d-12
integer :: newton_iter_max = 8
integer :: element_try_max = 6 ! Try newton iterations in at most this many elements
! A loop can occur here!

!> Internal variables
integer :: newton_iter_number
integer :: element_try_index
real*8 :: st_jac_det, R_s, R_t, Z_s, Z_t
real*8 :: st_step(2), st_try(2), x_step(2) ! x_step = (R,Z)
real*8 :: err2
logical :: changed, lost, search

! Setup initial values
x_step = x_old
i_elm_new = i_elm_old
! Find the jacobian at the current s and t position
call interp3_RZ(node_list,element_list,i_elm_old,st_old(1),st_old(2),x_step(1),R_s,R_t,&
                                                                     x_step(2),Z_s,Z_t)
st_jac_det = R_s * Z_t - R_t * Z_s



! Outer loop tries newton iteration, changes element if necessary
do element_try_index = 1, element_try_max
  ! Inner loop, newton iteration to find s and t in or out of this element
  do newton_iter_number = 1, newton_iter_max
    ! Perform newton iteration by calculating the inverse of the jacobian matrix

    ! Calculate the trial newton step
    st_step(1) = ( Z_t * (x_new(1)-x_step(1)) - R_t * (x_new(2)-x_step(2))) / st_jac_det
    st_step(2) = (-Z_s * (x_new(1)-x_step(1)) + R_s * (x_new(2)-x_step(2))) / st_jac_det

    ! Test this step
    st_try = st_new + st_step
    ! Calculate the x_step corresponding to st_try with the coordinates of element i_elm_new
    call interp3_RZ(node_list,element_list,i_elm_new,st_try(1),st_try(2),x_step(1),R_s,R_t,x_step(2),Z_s,Z_t)
    err2 = dot_product(x_step-x_new,x_step-x_new)
    if (st_try(1) > 1.d0 .or. st_try(1) < 0.d0 .or. st_try(2) > 1.d0 .or. st_try(2) < 0.d0) then
      if (err2 < out_element_tolerance) exit ! This loop is done!
    else
      if (err2 < in_element_tolerance) exit ! This loop is done!
    endif
    write(*,*) newton_iter_number, err2, st_jac_det
    ! Save this step value
    st_new = st_try
    ! Calculate the new determinant
    st_jac_det = R_s * Z_t - R_t * Z_s
    if (isnan(err2)) call exit(1)
  enddo
  if (newton_iter_number .gt. newton_iter_max) then
    write(*,"(A,i3,A)") "WARNING: newton iteration for st position did not converge after", newton_iter_max, " tries"
  endif

  ! Now check to see if we have left the current element
  call check_element_boundary(element_list,i_elm_new,st_try,i_elm_new,st_try,changed,lost,search)

  if (search) then
    ! Use the extremely slow option
    call find_RZ(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,st_new(1),st_new(2),ifail)
    write(*,*) "WARNING: used find_RZ"
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
    write(*,*) "changed", element_try_index, i_elm_new
    ! We have changed element, recalculate st_jac
    call interp3_RZ(node_list,element_list,i_elm_new,st_try(1),st_try(2),x_step(1),R_s,R_t,&
                                                                         x_step(2),Z_s,Z_t)
    st_jac_det = R_s * Z_t - R_t * Z_s
  endif

  st_new = st_try
enddo

if (changed) then
  write(*,*) "WARNING: insufficient iterations for element change, trying brute force method"
  call find_RZ(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,st_new(1),st_new(2),ifail)
else
  ! Save the final values
  st_new = st_try
endif
end subroutine find_RZ_nearby
