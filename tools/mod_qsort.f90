! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing
! copied from http://www.fortran.com/qsort_c.f95

! Made F conformant by Walt Brainerd
! Modified for integers by Daan van Vugt
module mod_qsort
implicit none
private
public :: qsort

interface qsort
  module procedure qsort_int
end interface qsort

contains

recursive subroutine qsort_int(A)
  integer, intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call partition_int(A, iq)
     call qsort_int(A(:iq-1))
     call qsort_int(A(iq:))
  endif
end subroutine qsort_int

subroutine partition_int(A, marker)
  integer, intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
  integer :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine partition_int

end module mod_qsort
