module mod_quicksort
  implicit none
  private
  public :: quicksort
contains
!< A simple recursive quicksort implementation
pure recursive subroutine quicksort(a, first, last)
  real*8, intent(inout) :: a(:) !< Array to be sorted
  integer, intent(in)   :: first, last !< First and last element of subsection to sort
  real*8 :: x, t
  integer :: i, j

  x = a((first+last)/2)
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     ! Swap a(i) and a(j)
     t = a(i)
     a(i) = a(j)
     a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort
end module mod_quicksort
