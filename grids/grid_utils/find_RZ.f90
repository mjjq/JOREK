recursive subroutine find_RZ2(node_list,element_list,R_find,Z_find,R_out,Z_out,ielm_out,s_out,t_out,ifail)
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

ielm_out = 0

!write(*,'(A,2e16.8)') ' find_RZ : ',R_find,Z_find

do k=1,element_list%n_elements

  call RZ_minmax(node_list,element_list,k,Rmin,Rmax,Zmin,Zmax)

  Rmin = Rmin - 0.05
  Rmax = Rmax + 0.05
  Zmin = Zmin - 0.05
  Zmax = Zmax + 0.05

  if ((R_find .ge. Rmin) .and. (R_find .le. Rmax) .and. &
      (Z_find .ge. Zmin) .and. (Z_find .le. Zmax) ) then

    ntrial = 20
    tolx = 1.d-6
    tolf = 1.d-15

    do istart = 1,5

      if (istart .eq. 1) then
        x(1) = 0.5d0
        x(2) = 0.5d0
      elseif (istart .eq. 2) then
        x(1) = -0.71d0
        x(2) = -0.71d0
      elseif (istart .eq. 3) then
        x(1) =  0.71d0
        x(2) = -0.71d0
      elseif (istart .eq. 4) then
        x(1) =  0.71d0
        x(2) =  0.71d0
      elseif (istart .eq. 5) then
        x(1) = -0.71d0
        x(2) =  0.71d0
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

        x = max(x,-0.000d0)
        x = min(x,+1.000d0)

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
end subroutine find_RZ2





!> Optimized version of find_RZ.
recursive subroutine find_RZ(node_list, element_list, R_find, Z_find, R_out, Z_out, ielm_out,      &
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

end subroutine find_RZ
