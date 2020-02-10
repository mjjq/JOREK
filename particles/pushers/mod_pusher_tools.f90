!> This modules contains tools needed by the different
!> particle pushers. Generally, they are mathematical 
!> stand-alone routines which can be used by multiple
!> independent algorithm
module mod_pusher_tools
implicit none
private
!> public procedures
public get_orthonormals
public cayley_transform,approximated_cayley_transform
contains

!---------------------------------------------------------------------------

!> Get two vectors orthogonal to a given vector.
!> This is the RZPhi-version, the right-handed version will have different directions but will also work.
pure subroutine get_orthonormals(b, e1, e2)
  use mod_math_operators, only: cross_product
  implicit none

  real*8, dimension(3), intent(in)  :: b !< Does not need to be normalized
  real*8, dimension(3), intent(out) :: e1, e2

  ! Take r as a reference vector (this will fail therefore if B is purely in r direction)
  if (b(2) .eq. 0.d0 .and. b(3) .eq. 0.d0) then
    e1 = [0.d0, 1.d0, 0.d0]
    e2 = [0.d0, 0.d0, 1.d0]
  else
    e1 = cross_product(b, [1.d0, 0.d0, 0.d0])
    ! Normalize
    e1 = e1/norm2(e1)
    ! Obtain a second reference vector
    e2 = cross_product(b, e1)
    ! Normalize
    e2 = e2/norm2(e2)
  end if
end subroutine get_orthonormals

!---------------------------------------------------------------------------

!> This function computes the right-handed Cayley transform of a vector vec multiplied 
!> by a scalar alpha. The Cayley transform is defined as:
!> cayley(alpha*B) = (I-alpha*B)^(-1) * (I+alpha*B)
!> where B is the vector product skew symmetric matrix of the vector vec
!> and I is the
pure function cayley_transform(alpha,vec)
  ! defining input variables
  real(kind=8),intent(in) :: alpha !< multiplicative constant
  real(kind=8),dimension(3),intent(in) :: vec !< vector to be transformed
  ! defining output variables
  real(kind=8),dimension(3,3) :: cayley_transform !< Cayley transform of vec
  ! defining internal variables
  real(kind=8),dimension(3,3) :: A,B !< (I-alpha*B)^(-1) and (I+alpha*B)

! computing (I+alpha*B)
  B(1:3,1) = (/1.d0,(-alpha*vec(3)),alpha*vec(2)/)
  B(1:3,2) = (/alpha*vec(3),1.d0,(-alpha*vec(1))/)
  B(1:3,3) = (/(-alpha*vec(2)),alpha*vec(1),1.d0/)

! computing (I-alpha*B)^(-1)
  A(1:3,1) = (/(1.0 + alpha*alpha*vec(1)*vec(1)),(alpha*(alpha*vec(1)*vec(2) - vec(3))),&
            (alpha*(alpha*vec(3)*vec(1) + vec(2)))/)
  A(1:3,2) = (/(alpha*(alpha*vec(2)*vec(1) + vec(3))),(1.0 + alpha*alpha*vec(2)*vec(2)),&
            (alpha*(alpha*vec(3)*vec(2) - vec(1)))/)
  A(1:3,3) = (/(alpha*(alpha*vec(3)*vec(1) - vec(2))),(alpha*(alpha*vec(2)*vec(3) + vec(1))),&
            (1.0 + alpha*alpha*vec(3)*vec(3))/)

! computing the Cayley transform
  cayley_transform = (matmul(A,B))/(1.d0 + alpha*alpha*(dot_product(vec,vec)))

end function cayley_transform

!---------------------------------------------------------------------------

!> This function computes the approximated Cayley transform a vector vec 
!> multiplied by a constant alpha as reported in R. Zhang et al., Phys. Plasmas 22 (2015) 044501.
pure function approximated_cayley_transform(alpha,vec)
  ! defining input variables
  real(kind=8),intent(in) :: alpha !< multiplicative constant
  real(kind=8),dimension(3),intent(in) :: vec !< vector to be transformed
  ! defining output variables
  real(kind=8),dimension(3,3) :: approximated_cayley_transform !< approximated Cayley transform of V
  ! internal variables
  real(kind=8) :: coefficient

  ! compute the approximated transfrom coefficient
  coefficient = (2.0*alpha)/(1.0 + alpha*alpha*(dot_product(vec,vec)))

  ! compute the approximated Cayley transform
  approximated_cayley_transform(1:3,1) = &
    (/1.0-coefficient*alpha*(vec(3)*vec(3) + vec(2)*vec(2)),&
    coefficient*(vec(3) + alpha*vec(2)*vec(1)),&
    coefficient*(alpha*vec(3)*vec(1) - vec(2))/)
  approximated_cayley_transform(1:3,2) = &
    (/coefficient*(alpha*vec(1)*vec(2) - vec(3)),&
    1.0-coefficient*alpha*(vec(3)*vec(3) + vec(1)*vec(1)),&
    coefficient*(alpha*vec(3)*vec(2) + vec(1))/)
  approximated_cayley_transform(1:3,3) = &
    (/coefficient*(alpha*vec(1)*vec(3) + vec(2)),&
    coefficient*(alpha*vec(2)*vec(3) - vec(1)),1.0-&
    coefficient*alpha*(vec(2)*vec(2) + vec(1)*vec(1))/)
end function approximated_cayley_transform

!---------------------------------------------------------------------------

end module mod_pusher_tools
