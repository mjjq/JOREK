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
end
