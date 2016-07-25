!> Subroutine which defines the basis functions derived from
!! a mixed Bezier/Cubic finite element representation.
!!
!! - index 1 : counts the vertex
!! - index 2 : counts the variables (p,u,v,w)
!! - the functions are defined on the interval [0,1][0,1]
!! \see ::basisfunctions and ::basisfunctions1
pure subroutine basisfunctions3(s, t, H, H_s, H_t)

implicit none

! --- Routine parameters
real*8, intent(in)  :: s          !< s-coordinate in the element
real*8, intent(in)  :: t          !< t-coordinate in the element
real*8, intent(out) :: H(4,4)     !< Basis functions
real*8, intent(out) :: H_s(4,4)   !< Basis functions derived with respect to s
real*8, intent(out) :: H_t(4,4)   !< Basis functions derived with respect to t

!---------------------------------------------------------- vertex (1)
H(1,1)   =(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_s(1,1) =6.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(1,1) =6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t

H(1,2)   =3.d0*(-1.d0 + s)**2*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_s(1,2) =3.d0*(-1.d0 + s)*(-1.d0 + 3.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(1,2) =18.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t

H(1,3)   =3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)**2*t
H_s(1,3) =18.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*t
H_t(1,3) =3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*(-1.d0 + 3.d0*t)

H(1,4)   =9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)**2*t
H_s(1,4) =9.d0*(-1.d0 + s)*(-1.d0 + 3.d0*s)*(-1.d0 + t)**2*t
H_t(1,4) =9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*(-1.d0 + 3.d0*t)

!---------------------------------------------------------- vertex (2)
H(2,1)   =-(s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t))
H_s(2,1) =-6.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(2,1) =-6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t

H(2,2)   =-3.d0*(-1.d0 + s)*s**2*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_s(2,2) =-3.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(2,2) =-18.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t

H(2,3)   =-3.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)**2*t
H_s(2,3) =-18.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*t
H_t(2,3) =3.d0*s**2*(-3.d0 + 2.d0*s)*(1.d0 - 3.d0*t)*(-1.d0 + t)

H(2,4)   =-9.d0*(-1.d0 + s)*s**2*(-1.d0 + t)**2*t
H_s(2,4) =-9.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*t
H_t(2,4) =9.d0*(-1.d0 + s)*s**2*(1.d0 - 3.d0*t)*(-1.d0 + t)

!---------------------------------------------------------- vertex (3)
H(3,1)   =s**2*(-3.d0 + 2.d0*s)*t**2*(-3.d0 + 2.d0*t)
H_s(3,1) =6.d0*(-1.d0 + s)*s*t**2*(-3.d0 + 2.d0*t)
H_t(3,1) =6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t

H(3,2)   =3.d0*(-1.d0 + s)*s**2*t**2*(-3.d0 + 2.d0*t)
H_s(3,2) =3.d0*s*(-2 + 3.d0*s)*t**2*(-3.d0 + 2.d0*t)
H_t(3,2) =18.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t

H(3,3)   =3.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t**2
H_s(3,3) =18.d0*(-1.d0 + s)*s*(-1.d0 + t)*t**2
H_t(3,3) =3.d0*s**2*(-3.d0 + 2.d0*s)*t*(-2.d0 + 3.d0*t)

H(3,4)   =9.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t**2
H_s(3,4) =9.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)*t**2
H_t(3,4) =9.d0*(-1.d0 + s)*s**2*t*(-2.d0 + 3.d0*t)

!---------------------------------------------------------- vertex (4)
H(4,1)   =-((-1.d0 + s)**2*(1.d0 + 2.d0*s)*t**2*(-3.d0 + 2.d0*t))
H_s(4,1) =-6.d0*(-1.d0 + s)*s*t**2*(-3.d0 + 2.d0*t)
H_t(4,1) =-6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t

H(4,2)   =-3.d0*(-1.d0 + s)**2*s*t**2*(-3.d0 + 2.d0*t)
H_s(4,2) =3.d0*(1.d0 - 3*s)*(-1.d0 + s)*t**2*(-3.d0 + 2.d0*t)
H_t(4,2) =-18.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t

H(4,3)   =-3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t**2
H_s(4,3) =-18.d0*(-1.d0 + s)*s*(-1.d0 + t)*t**2
H_t(4,3) =-3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*t*(-2.d0 + 3.d0*t)

H(4,4)   =-9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t**2
H_s(4,4) =9.d0*(1.d0 - 3.d0*s)*(-1.d0 + s)*(-1.d0 + t)*t**2
H_t(4,4) =-9.d0*(-1.d0 + s)**2*s*t*(-2.d0 + 3.d0*t)

return
end subroutine basisfunctions3
