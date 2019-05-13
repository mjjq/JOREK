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
  A   = [(1.d0 - 2.d0*(t - t0)*dli(1))*li(1)**2, &
         (1.d0 - 2.d0*(t - t1)*dli(2))*li(2)**2]
  B   = [(t - t0)*li(1)**2, (t - t1)*li(2)**2]

  y = y0 * A(1) + y1 * A(2) + dy0 * B(1) + dy1 * B(2)
end subroutine HB_interp

pure subroutine HB_interp_dt(t0, t1, n, y0, y1, dy0, dy1, t, y)
  real*8, intent(in)                :: t0, t1
  integer, intent(in)               :: n
  real*8, intent(in), dimension(n)  :: y0, y1, dy0, dy1
  real*8, intent(in)                :: t
  real*8, intent(out), dimension(n) :: y

  real*8, dimension(2) :: ti, li, dli, A, B

  ! Prepare needed variables
  ti  = [t0,t1]
  li  = [(t - t1), (t0 - t)]/(t0 - t1)
  dli = [1.d0/(t0 - t1), 1.d0/(t1 - t0)]
  A   = -2.d0*dli*li**2 + (1.d0-2.d0*(t-ti)*dli)*2.d0*li*dli
  B   = li**2 + (t - ti)*2.d0*li*dli

  y = y0 * A(1) + y1 * A(2) + dy0 * B(1) + dy1 * B(2)
end subroutine HB_interp_dt
end module mod_hermite_birkhoff
