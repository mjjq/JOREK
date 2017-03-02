!> Module to find roots of functions
module mod_rootfinding
  implicit none
  private
  public :: newtons_method
  public :: halleys_method
contains
  !> Use newton's method to solve f(x) == y0, starting at x0
  pure subroutine newtons_method(f, df, y0, x0, x, ierr)
    !real*8, external     :: f     !< Function value at x
    !real*8, external     :: df    !< Derivative at x
    real*8, intent(in)   :: y0    !< Intersection to find
    real*8, intent(in)   :: x0    !< Initial value
    real*8, intent(out)  :: x     !< Result value
    integer, intent(out) :: ierr  !< Status code. If == 0 we found a result
    interface
      pure function f(x)
        real*8, intent(in) :: x
      end function f
      pure function df(x)
        real*8, intent(in) :: x
      end function df
    end interface

    integer, parameter :: n_iter = 10
    real*8, parameter  :: tolerance = 1d-10

    integer :: i
    real*8  :: y
    
    ierr = 0
    x = x0
    do i=1,n_iter
      y = f(x) - y0
      x = x - y/df(x)

      if (abs(y) .le. tolerance) return
    end do
    ierr = 1 ! We did not find a root
  end subroutine newtons_method

  !> Use Halley's method to solve f(x) == y0, starting at x0
  pure subroutine halleys_method(f, df, ddf, y0, x0, x, ierr)
    !real*8, external     :: f     !< Function value at x
    !real*8, external     :: df    !< First derivative at x
    !real*8, external     :: ddf   !< Second derivative at x
    real*8, intent(in)   :: y0    !< Intersection to find
    real*8, intent(in)   :: x0    !< Initial value
    real*8, intent(out)  :: x     !< Result value
    integer, intent(out) :: ierr  !< Status code. If == 0 we found a result
    interface
      pure function f(x)
        real*8, intent(in) :: x
      end function f
      pure function df(x)
        real*8, intent(in) :: x
      end function df
      pure function ddf(x)
        real*8, intent(in) :: x
      end function ddf
    end interface

    integer, parameter :: n_iter = 10
    real*8, parameter  :: tolerance = 1d-10

    integer :: i
    real*8  :: y, dy, ddy
    
    ierr = 0
    x = x0
    do i=1,n_iter
      y   = f(x) - y0
      dy  = df(x)
      ddy = ddf(x)
      x = x - 2.d0*y*dy/(2.d0*dy*dy - y*ddy)

      if (abs(y) .le. tolerance) return
    end do
    ierr = 1 ! We did not find a root
  end subroutine halleys_method
end module mod_rootfinding
