subroutine find_RZ(node_list,element_list,R_find,Z_find,R_out,Z_out,ielm_out,s_out,t_out,ifail)
!-------------------------------------------------------------------------
!< Find all elements for which minmax is correct and run find_RZ_single on those.
!< Return the first result.
!-------------------------------------------------------------------------
use data_structure
use mod_element_rtree
implicit none

type (type_node_list), intent(in)    :: node_list
type (type_element_list), intent(in) :: element_list
real*8, intent(in)     :: R_find, Z_find
real*8, intent(out)    :: R_out,Z_out,s_out,t_out
integer, intent(inout) :: ielm_out
integer, intent(out)   :: ifail

logical, save        :: rtree_initialised = .false.
integer, save        :: n_elements = 0 !< Number of elements in the existing rtree. Used to check if the grid has changed.
real*8, save         :: r_median = 0.d0 !< Position of the median node. Used to check if the grid has changed.

integer :: k
integer, dimension(:), allocatable :: i_elms

ielm_out = 0
if (element_list%n_elements .ne. n_elements) rtree_initialised = .false.
if (node_list%n_nodes .gt. 0 .and. abs(node_list%node((node_list%n_nodes+1)/2)%x(1,1) - r_median) .gt. 1d-80) rtree_initialised = .false.
if (.not. rtree_initialised) then
  call populate_element_rtree(node_list, element_list) ! not OMP safe, call once outside of openmp
  rtree_initialised = .true.
  n_elements = element_list%n_elements
  if (node_list%n_nodes .gt. 0) r_median = node_list%node((node_list%n_nodes+1)/2)%x(1,1)
end if

call elements_containing_point(node_list, element_list, R_find, Z_find, i_elms)

! then loop through all
do k=1,size(i_elms)
  call find_RZ_single(node_list,element_list,i_elms(k),R_find,Z_find,R_out,Z_out,ielm_out,s_out,t_out,ifail)
  if (ifail .eq. 0) return
enddo

if (ielm_out .eq. 0) ifail = 99
end subroutine find_RZ


subroutine find_RZ_single(node_list,element_list,i_elm,R_find,Z_find,R_out,Z_out,ielm_out,s_out,t_out,ifail)
!-------------------------------------------------------------------------
!< solves two non-linear equations using Newtons method (from numerical recipes)
!< LU decomposition replaced by explicit solution of 2x2 matrix.
!<
!< finds the crossing of two coordinate lines given as a series of cubics in element
!< i_elm
!-------------------------------------------------------------------------
use data_structure
implicit none

type (type_node_list), intent(in)    :: node_list
type (type_element_list), intent(in) :: element_list
integer, intent(in)    :: i_elm
real*8, intent(in)     :: R_find, Z_find
real*8, intent(out)    :: R_out,Z_out,s_out,t_out
integer, intent(out)   :: ielm_out
integer, intent(out)   :: ifail

integer :: i, ntrial, istart
real*8  :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8  :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8  :: tolx, tolf, errx, errf, temp, dis
real*8  :: x(2), FVEC(2), FJAC(2,2), p(2)

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

    call interp_RZ(node_list,element_list,i_elm,x(1),x(2),RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
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

      ielm_out  = i_elm
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

      ielm_out  = i_elm
      R_out     = RRg1
      Z_out     = ZZg1

!        write(*,'(A,i3,4e16.8)') ' newton (2) : ',i,errf,errx,x

      ifail = 0
      return
    endif

  enddo
enddo
end subroutine find_RZ_single
