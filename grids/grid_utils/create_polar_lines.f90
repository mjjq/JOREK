subroutine create_polar_lines(n_polar, R_beg, Z_beg, R_mid, Z_mid, R_end, Z_end, delta, R_polar, Z_polar)
! --- This routine creates a cubic spline for the polar lines (the ones that are
! --- across flux surfaces). They have 3 control points, two at the end, and one
! --- in the middle, at the turning point, typically the separatrix

  implicit none
  
  ! --- Routine variables
  integer, intent(in)  :: n_polar
  real*8,  intent(in)  :: R_beg(n_polar), Z_beg(n_polar)
  real*8,  intent(in)  :: R_mid(n_polar), Z_mid(n_polar)
  real*8,  intent(in)  :: R_end(n_polar), Z_end(n_polar)
  real*8,  intent(in)  :: delta(n_polar)
  real*8,  intent(out) :: R_polar(3,4,n_polar), Z_polar(3,4,n_polar)
  
  ! --- Internal variables
  integer :: i

  !------------------------------ Central part (without legs)
  do i=1,n_polar

    R_polar(1,1,i) = R_beg(i)
    R_polar(1,4,i) = delta(i) * R_beg(i) + (1.d0 - delta(i)) * R_mid(i)
    R_polar(1,2,i) = ( 2.d0 * R_polar(1,1,i)  +         R_polar(1,4,i) ) / 3.d0
    R_polar(1,3,i) = (        R_polar(1,1,i)  +  2.d0 * R_polar(1,4,i) ) / 3.d0

    Z_polar(1,1,i) = Z_beg(i)
    Z_polar(1,4,i) = delta(i) * Z_beg(i) + (1.d0 - delta(i)) * Z_mid(i)
    Z_polar(1,2,i) = ( 2.d0 * Z_polar(1,1,i)  +         Z_polar(1,4,i) ) / 3.d0
    Z_polar(1,3,i) = (        Z_polar(1,1,i)  +  2.d0 * Z_polar(1,4,i) ) / 3.d0

    R_polar(3,1,i) = R_end(i)
    R_polar(3,4,i) = delta(i) * R_end(i) + (1.d0 - delta(i)) * R_mid(i)
    R_polar(3,2,i) = ( 2.d0 * R_polar(3,1,i)  +         R_polar(3,4,i) ) / 3.d0
    R_polar(3,3,i) = (        R_polar(3,1,i)  +  2.d0 * R_polar(3,4,i) ) / 3.d0

    Z_polar(3,1,i) = Z_end(i)
    Z_polar(3,4,i) = delta(i) * Z_end(i) + (1.d0 - delta(i)) * Z_mid(i)
    Z_polar(3,2,i) = ( 2.d0 * Z_polar(3,1,i)  +         Z_polar(3,4,i) ) / 3.d0
    Z_polar(3,3,i) = (        Z_polar(3,1,i)  +  2.d0 * Z_polar(3,4,i) ) / 3.d0

    R_polar(2,1,i) = R_polar(1,4,i)
    R_polar(2,4,i) = R_polar(3,4,i)
    R_polar(2,2,i) = ( R_polar(2,1,i) +  2.d0 * R_mid(i) ) / 3.d0
    R_polar(2,3,i) = ( R_polar(2,4,i) +  2.d0 * R_mid(i) ) / 3.d0

    Z_polar(2,1,i) = Z_polar(1,4,i)
    Z_polar(2,4,i) = Z_polar(3,4,i)
    Z_polar(2,2,i) = ( Z_polar(2,1,i) +  2.d0 * Z_mid(i) ) / 3.d0
    Z_polar(2,3,i) = ( Z_polar(2,4,i) +  2.d0 * Z_mid(i) ) / 3.d0

  enddo

  return

end subroutine create_polar_lines







subroutine create_polar_lines_4(n_polar, R_beg, Z_beg, R_mid1, Z_mid1, R_mid2, Z_mid2, R_end, Z_end, delta1, delta2, R_polar, Z_polar)

  implicit none
  
  ! --- Routine variables
  integer, intent(in)  :: n_polar
  real*8,  intent(in)  :: R_beg (n_polar), Z_beg (n_polar)
  real*8,  intent(in)  :: R_mid1(n_polar), Z_mid1(n_polar)
  real*8,  intent(in)  :: R_mid2(n_polar), Z_mid2(n_polar)
  real*8,  intent(in)  :: R_end (n_polar), Z_end (n_polar)
  real*8,  intent(in)  :: delta1(n_polar), delta2(n_polar)
  real*8,  intent(out) :: R_polar(5,4,n_polar), Z_polar(5,4,n_polar)
  
  ! --- Internal variables
  integer :: i

  !------------------------------ Central part (without legs)
  do i=1,n_polar

    R_polar(1,1,i) = R_beg(i)
    R_polar(1,4,i) = delta1(i) * R_beg(i) + (1.d0 - delta1(i)) * R_mid1(i)
    R_polar(1,2,i) = ( 2.d0 * R_polar(1,1,i)  +         R_polar(1,4,i) ) / 3.d0
    R_polar(1,3,i) = (        R_polar(1,1,i)  +  2.d0 * R_polar(1,4,i) ) / 3.d0

    Z_polar(1,1,i) = Z_beg(i)
    Z_polar(1,4,i) = delta1(i) * Z_beg(i) + (1.d0 - delta1(i)) * Z_mid1(i)
    Z_polar(1,2,i) = ( 2.d0 * Z_polar(1,1,i)  +         Z_polar(1,4,i) ) / 3.d0
    Z_polar(1,3,i) = (        Z_polar(1,1,i)  +  2.d0 * Z_polar(1,4,i) ) / 3.d0

    R_polar(3,1,i) = delta1(i) * R_mid2(i) + (1.d0 - delta1(i)) * R_mid1(i)
    R_polar(3,4,i) = delta2(i) * R_mid1(i) + (1.d0 - delta2(i)) * R_mid2(i)
    R_polar(3,2,i) = ( 2.d0 * R_polar(3,1,i)  +         R_polar(3,4,i) ) / 3.d0
    R_polar(3,3,i) = (        R_polar(3,1,i)  +  2.d0 * R_polar(3,4,i) ) / 3.d0

    Z_polar(3,1,i) = delta1(i) * Z_mid2(i) + (1.d0 - delta1(i)) * Z_mid1(i)
    Z_polar(3,4,i) = delta2(i) * Z_mid1(i) + (1.d0 - delta2(i)) * Z_mid2(i)
    Z_polar(3,2,i) = ( 2.d0 * Z_polar(3,1,i)  +         Z_polar(3,4,i) ) / 3.d0
    Z_polar(3,3,i) = (        Z_polar(3,1,i)  +  2.d0 * Z_polar(3,4,i) ) / 3.d0

    R_polar(5,1,i) = R_end(i)
    R_polar(5,4,i) = delta2(i) * R_end(i) + (1.d0 - delta2(i)) * R_mid2(i)
    R_polar(5,2,i) = ( 2.d0 * R_polar(5,1,i)  +         R_polar(5,4,i) ) / 3.d0
    R_polar(5,3,i) = (        R_polar(5,1,i)  +  2.d0 * R_polar(5,4,i) ) / 3.d0

    Z_polar(5,1,i) = Z_end(i)
    Z_polar(5,4,i) = delta2(i) * Z_end(i) + (1.d0 - delta2(i)) * Z_mid2(i)
    Z_polar(5,2,i) = ( 2.d0 * Z_polar(5,1,i)  +         Z_polar(5,4,i) ) / 3.d0
    Z_polar(5,3,i) = (        Z_polar(5,1,i)  +  2.d0 * Z_polar(5,4,i) ) / 3.d0

    R_polar(2,1,i) = R_polar(1,4,i)
    R_polar(2,4,i) = R_polar(3,1,i)
    R_polar(2,2,i) = ( R_polar(2,1,i) +  2.d0 * R_mid1(i) ) / 3.d0
    R_polar(2,3,i) = ( R_polar(2,4,i) +  2.d0 * R_mid1(i) ) / 3.d0

    Z_polar(2,1,i) = Z_polar(1,4,i)
    Z_polar(2,4,i) = Z_polar(3,1,i)
    Z_polar(2,2,i) = ( Z_polar(2,1,i) +  2.d0 * Z_mid1(i) ) / 3.d0
    Z_polar(2,3,i) = ( Z_polar(2,4,i) +  2.d0 * Z_mid1(i) ) / 3.d0

    R_polar(4,1,i) = R_polar(5,4,i)
    R_polar(4,4,i) = R_polar(3,4,i)
    R_polar(4,2,i) = ( R_polar(4,1,i) +  2.d0 * R_mid2(i) ) / 3.d0
    R_polar(4,3,i) = ( R_polar(4,4,i) +  2.d0 * R_mid2(i) ) / 3.d0

    Z_polar(4,1,i) = Z_polar(5,4,i)
    Z_polar(4,4,i) = Z_polar(3,4,i)
    Z_polar(4,2,i) = ( Z_polar(4,1,i) +  2.d0 * Z_mid2(i) ) / 3.d0
    Z_polar(4,3,i) = ( Z_polar(4,4,i) +  2.d0 * Z_mid2(i) ) / 3.d0

  enddo

  return

end subroutine create_polar_lines_4



subroutine create_polar_lines_5(n_polar, R_mid, Z_mid, delta, R_polar, Z_polar)

  implicit none
  
  ! --- Routine variables
  integer, intent(in)  :: n_polar
  real*8,  intent(in)  :: R_mid(n_polar,5), Z_mid(n_polar,5)
  real*8,  intent(in)  :: delta(n_polar,3)
  real*8,  intent(out) :: R_polar(7,4,n_polar), Z_polar(7,4,n_polar)
  
  ! --- Internal variables
  integer :: i

  !------------------------------ Central part (without legs)
  do i=1,n_polar

    R_polar(1,1,i) = R_mid(i,1)
    R_polar(1,4,i) = delta(i,1) * R_mid(i,1) + (1.d0 - delta(i,1)) * R_mid(i,2)
    R_polar(1,2,i) = ( 2.d0 * R_polar(1,1,i)  +         R_polar(1,4,i) ) / 3.d0
    R_polar(1,3,i) = (        R_polar(1,1,i)  +  2.d0 * R_polar(1,4,i) ) / 3.d0

    Z_polar(1,1,i) = Z_mid(i,1)
    Z_polar(1,4,i) = delta(i,1) * Z_mid(i,1) + (1.d0 - delta(i,1)) * Z_mid(i,2)
    Z_polar(1,2,i) = ( 2.d0 * Z_polar(1,1,i)  +         Z_polar(1,4,i) ) / 3.d0
    Z_polar(1,3,i) = (        Z_polar(1,1,i)  +  2.d0 * Z_polar(1,4,i) ) / 3.d0

    R_polar(3,1,i) = delta(i,1) * R_mid(i,3) + (1.d0 - delta(i,1)) * R_mid(i,2)
    R_polar(3,4,i) = delta(i,2) * R_mid(i,2) + (1.d0 - delta(i,2)) * R_mid(i,3)
    R_polar(3,2,i) = ( 2.d0 * R_polar(3,1,i)  +         R_polar(3,4,i) ) / 3.d0
    R_polar(3,3,i) = (        R_polar(3,1,i)  +  2.d0 * R_polar(3,4,i) ) / 3.d0

    Z_polar(3,1,i) = delta(i,1) * Z_mid(i,3) + (1.d0 - delta(i,1)) * Z_mid(i,2)
    Z_polar(3,4,i) = delta(i,2) * Z_mid(i,2) + (1.d0 - delta(i,2)) * Z_mid(i,3)
    Z_polar(3,2,i) = ( 2.d0 * Z_polar(3,1,i)  +         Z_polar(3,4,i) ) / 3.d0
    Z_polar(3,3,i) = (        Z_polar(3,1,i)  +  2.d0 * Z_polar(3,4,i) ) / 3.d0

    R_polar(5,1,i) = delta(i,2) * R_mid(i,4) + (1.d0 - delta(i,2)) * R_mid(i,3)
    R_polar(5,4,i) = delta(i,3) * R_mid(i,3) + (1.d0 - delta(i,3)) * R_mid(i,4)
    R_polar(5,2,i) = ( 2.d0 * R_polar(5,1,i)  +         R_polar(5,4,i) ) / 3.d0
    R_polar(5,3,i) = (        R_polar(5,1,i)  +  2.d0 * R_polar(5,4,i) ) / 3.d0

    Z_polar(5,1,i) = delta(i,2) * Z_mid(i,4) + (1.d0 - delta(i,2)) * Z_mid(i,3)
    Z_polar(5,4,i) = delta(i,3) * Z_mid(i,3) + (1.d0 - delta(i,3)) * Z_mid(i,4)
    Z_polar(5,2,i) = ( 2.d0 * Z_polar(5,1,i)  +         Z_polar(5,4,i) ) / 3.d0
    Z_polar(5,3,i) = (        Z_polar(5,1,i)  +  2.d0 * Z_polar(5,4,i) ) / 3.d0

    R_polar(7,1,i) = R_mid(i,5)
    R_polar(7,4,i) = delta(i,3) * R_mid(i,5) + (1.d0 - delta(i,3)) * R_mid(i,4)
    R_polar(7,2,i) = ( 2.d0 * R_polar(7,1,i)  +         R_polar(7,4,i) ) / 3.d0
    R_polar(7,3,i) = (        R_polar(7,1,i)  +  2.d0 * R_polar(7,4,i) ) / 3.d0

    Z_polar(7,1,i) = Z_mid(i,5)
    Z_polar(7,4,i) = delta(i,3) * Z_mid(i,5) + (1.d0 - delta(i,3)) * Z_mid(i,4)
    Z_polar(7,2,i) = ( 2.d0 * Z_polar(7,1,i)  +         Z_polar(7,4,i) ) / 3.d0
    Z_polar(7,3,i) = (        Z_polar(7,1,i)  +  2.d0 * Z_polar(7,4,i) ) / 3.d0

    R_polar(2,1,i) = R_polar(1,4,i)
    R_polar(2,4,i) = R_polar(3,1,i)
    R_polar(2,2,i) = ( R_polar(2,1,i) +  2.d0 * R_mid(i,2) ) / 3.d0
    R_polar(2,3,i) = ( R_polar(2,4,i) +  2.d0 * R_mid(i,2) ) / 3.d0

    Z_polar(2,1,i) = Z_polar(1,4,i)
    Z_polar(2,4,i) = Z_polar(3,1,i)
    Z_polar(2,2,i) = ( Z_polar(2,1,i) +  2.d0 * Z_mid(i,2) ) / 3.d0
    Z_polar(2,3,i) = ( Z_polar(2,4,i) +  2.d0 * Z_mid(i,2) ) / 3.d0

    R_polar(4,1,i) = R_polar(3,4,i)
    R_polar(4,4,i) = R_polar(5,4,i)
    R_polar(4,2,i) = ( R_polar(4,1,i) +  2.d0 * R_mid(i,3) ) / 3.d0
    R_polar(4,3,i) = ( R_polar(4,4,i) +  2.d0 * R_mid(i,3) ) / 3.d0

    Z_polar(4,1,i) = Z_polar(3,4,i)
    Z_polar(4,4,i) = Z_polar(5,4,i)
    Z_polar(4,2,i) = ( Z_polar(4,1,i) +  2.d0 * Z_mid(i,3) ) / 3.d0
    Z_polar(4,3,i) = ( Z_polar(4,4,i) +  2.d0 * Z_mid(i,3) ) / 3.d0

    R_polar(6,1,i) = R_polar(5,4,i)
    R_polar(6,4,i) = R_polar(7,4,i)
    R_polar(6,2,i) = ( R_polar(4,1,i) +  2.d0 * R_mid(i,4) ) / 3.d0
    R_polar(6,3,i) = ( R_polar(4,4,i) +  2.d0 * R_mid(i,4) ) / 3.d0

    Z_polar(6,1,i) = Z_polar(5,4,i)
    Z_polar(6,4,i) = Z_polar(7,4,i)
    Z_polar(6,2,i) = ( Z_polar(4,1,i) +  2.d0 * Z_mid(i,4) ) / 3.d0
    Z_polar(6,3,i) = ( Z_polar(4,4,i) +  2.d0 * Z_mid(i,4) ) / 3.d0

  enddo

  return

end subroutine create_polar_lines_5







subroutine create_polar_lines_simple(n_points, R_points, Z_points, R_polar, Z_polar)

  implicit none
  
  ! --- Routine variables
  integer, intent(in)  :: n_points
  real*8,  intent(in)  :: R_points(n_points), Z_points(n_points)
  real*8,  intent(out) :: R_polar(n_points-1,4), Z_polar(n_points-1,4)
  
  ! --- Internal variables
  integer :: i

  do i=1,n_points-1

    R_polar(i,1) = R_points(i)
    R_polar(i,4) = R_points(i+1)
    R_polar(i,2) = ( 2.d0 * R_polar(i,1)  +         R_polar(i,4) ) / 3.d0
    R_polar(i,3) = (        R_polar(i,1)  +  2.d0 * R_polar(i,4) ) / 3.d0

    Z_polar(i,1) = Z_points(i)
    Z_polar(i,4) = Z_points(i+1)
    Z_polar(i,2) = ( 2.d0 * Z_polar(i,1)  +         Z_polar(i,4) ) / 3.d0
    Z_polar(i,3) = (        Z_polar(i,1)  +  2.d0 * Z_polar(i,4) ) / 3.d0

  enddo

  return

end subroutine create_polar_lines_simple







