!> Hermite-Birkhoff interpolation of vectors and matrices
!> See https://arxiv.org/pdf/1704.08955.pdf page 7 for a description.
module mod_hermite_birkhoff
public :: HB_interp, HB_interp_dt
contains

pure subroutine HB_interp(t0, t1, n, y0, y1, dy0, dy1, t, y)
  real*8, intent(in)                :: t0, t1
  integer, intent(in)               :: n
  real*8, intent(in), dimension(n)  :: y0, y1, dy0, dy1
  real*8, intent(in)                :: t
  real*8, intent(out), dimension(n) :: y

  real*8, dimension(2) :: li, dli, A, B

  ! Prepare needed variables
  li  = [(t - t1), (t0 - t)]/(t0 - t1)
  dli = [1.d0/(t0 - t1), 1.d0/(t1 - t0)]
  A   = [(1.d0 - 2.d0*(t - t0)*dli(1))*li(1)*li(1), &
         (1.d0 - 2.d0*(t - t1)*dli(2))*li(2)*li(2)]
  B   = [(t - t0)*li(1)*li(1), (t - t1)*li(2)*li(2)]

  y = y0 * A(1) + y1 * A(2) + dy0 * B(1) + dy1 * B(2)
end subroutine HB_interp

pure subroutine HB_interp_dt(t0, t1, n, y0, y1, dy0, dy1, t, y)
  real*8, intent(in)                :: t0, t1
  integer, intent(in)               :: n
  real*8, intent(in), dimension(n)  :: y0, y1, dy0, dy1
  real*8, intent(in)                :: t
  real*8, intent(out), dimension(n) :: y
  real*8 :: t_inv
  real*8, dimension(2) :: A, B

  ! Prepare needed variables
  t_inv = 1.d0/(t1-t0)
  A   = 6.d0*(t1-t)*(t0-t)*[1.d0,-1.d0]*(t_inv*t_inv*t_inv)
  B   = [(t1-t)*(t1+2.d0*t0-3.d0*t),&
       (t-t0)*(3.d0*t-2.d0*t1-t0)]*(t_inv*t_inv)

  y = y0 * A(1) + y1 * A(2) + dy0 * B(1) + dy1 * B(2)
end subroutine HB_interp_dt
end module mod_hermite_birkhoff
