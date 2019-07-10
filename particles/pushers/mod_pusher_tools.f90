!> This modules contains tools needed by the different
!> particle pushers. Generally, they are mathematical 
!> stand-alone routines which can be used by multiple
!> independent algorithm
module mod_pusher_tools
implicit none
private
!> public procedures
public get_orthonormals
public left_handed_cross_product, right_handed_cross_product
public vector_transform_RZPHI_to_XYZ,vector_transform_XYZ_to_RZPHI
public cayley_transform,approximated_cayley_transform
contains

!---------------------------------------------------------------------------

!> Get two vectors orthogonal to a given vector.
!> This is the RZPhi-version, the right-handed version will have different directions but will also work.
pure subroutine get_orthonormals(b, e1, e2)
  real*8, dimension(3), intent(in)  :: b !< Does not need to be normalized
  real*8, dimension(3), intent(out) :: e1, e2

  ! Take r as a reference vector (this will fail therefore if B is purely in r direction)
  if (b(2) .eq. 0.d0 .and. b(3) .eq. 0.d0) then
    e1 = [0.d0, 1.d0, 0.d0]
    e2 = [0.d0, 0.d0, 1.d0]
  else
    e1 = left_handed_cross_product(b, [1.d0, 0.d0, 0.d0])
    ! Normalize
    e1 = e1/norm2(e1)
    ! Obtain a second reference vector
    e2 = left_handed_cross_product(b, e1)
    ! Normalize
    e2 = e2/norm2(e2)
  end if
end subroutine get_orthonormals

!---------------------------------------------------------------------------

!> The cross product in a left-handed coordinate system (e.g. RZPhi)
pure function left_handed_cross_product(a, b)
  real*8, dimension(3) :: left_handed_cross_product
  real*8, dimension(3), intent(in) :: a, b

  left_handed_cross_product(1) = a(2) * b(3) - a(3) * b(2)
  left_handed_cross_product(2) = a(3) * b(1) - a(1) * b(3)
  left_handed_cross_product(3) = a(1) * b(2) - a(2) * b(1)
end function left_handed_cross_product

!---------------------------------------------------------------------------

!> The cross product in a right-handed coordinate system (e.g. XYZ or RPhiZ)
pure function right_handed_cross_product(a, b)
  real*8, dimension(3) :: right_handed_cross_product
  real*8, dimension(3), intent(in) :: a, b

  right_handed_cross_product(1) = a(2) * b(3) - a(3) * b(2)
  right_handed_cross_product(2) = a(3) * b(1) - a(1) * b(3)
  right_handed_cross_product(3) = a(1) * b(2) - a(2) * b(1)
end function right_handed_cross_product

!---------------------------------------------------------------------------

!> This function rotates a vector from a \{R,Z,phi\} basis to a \{X,Y,Z\}
!> This is done performing a first swap and reflaction of the Y axis
!> and then, a rotation of angle phi (clockwise)
pure function vector_transform_RZPHI_to_XYZ(phi,a) result(b)
  ! declare input variables
  real(kind=8), intent(in) :: phi
  real(kind=8), dimension(3), intent(in) :: a
  ! declare output variables
  real(kind=8), dimension(3) :: b
  ! declare internal variables
  real(kind=8),dimension(2) :: sincosphi

  ! computing sinus and cosinus
  sincosphi = (/sin(phi),cos(phi)/)

  b(1) = a(1)*sincosphi(2) - a(3)*sincosphi(1) 
  b(2) = -1.d0*(a(1)*sincosphi(1) + a(3)*sincosphi(2))
  b(3) = a(2)

end function vector_transform_RZPHI_to_XYZ

!---------------------------------------------------------------------------

!> This function rotates a vector from a \{X,Y,Z\} basis to a \{R,Z,phi\}
pure function vector_transform_XYZ_to_RZPHI(phi,a) result(b)
  ! declare input variables
  real(kind=8), intent(in) :: phi
  real(kind=8), dimension(3), intent(in) :: a
  ! declare output variables
  real(kind=8), dimension(3) :: b
  ! declare internal variables
  real(kind=8),dimension(2) :: sincosphi

  ! computing sinus and cosinus
  sincosphi = (/sin(phi),cos(phi)/)

  b(1) = a(1)*sincosphi(2) - a(2)*sincosphi(1) 
  b(2) = a(3)
  b(3) = -1.d0*(a(1)*sincosphi(1) + a(2)*sincosphi(2))

end function vector_transform_XYZ_to_RZPHI

!---------------------------------------------------------------------------

!> This function computes the right handed Cayley transform of a vector vec multiplied 
!> by a constant \alpha. The Cayley transform is defined as:
!> cayley(alpha\cdot B) = \(I-\alpha \cdot B\\)^\{-1\} \cdot \(I+\alpha \cdot b\)
!> Where B is the vector product skew symmetric matrix of the vector B,
!> \alpha is a constant and I is the identity matrix.
!> inputs:
!>   alpha: (real8) multiplicative constant
!>   vec:     (real8)(3) vector to be transformed
!> outputs:
!>   cayley_transofrm: (real8)(3,3) the Cayley transfrom of vec
pure function cayley_transform(alpha,vec)
! defining input variables
real(kind=8),intent(in) :: alpha !< multiplicative constant
real(kind=8),dimension(3),intent(in) :: vec !< vector to be transformed
! defining output variables
real(kind=8),dimension(3,3) :: cayley_transform !< Cayley transform of vec
! defining internal variables
real(kind=8),dimension(3,3) :: A,B !< \(I-\alpha \cdot B\\)^\{-1\} and \(I+\alpha \cdot b\)

! computing \(I+\alpha \cdot b\)
 B(1:3,1) = (/1.d0,(-alpha*vec(3)),alpha*vec(2)/)
 B(1:3,2) = (/alpha*vec(3),1.d0,(-alpha*vec(1))/)
 B(1:3,3) = (/(-alpha*vec(2)),alpha*vec(1),1.d0/)

!computing \(I-\alpha \cdot B\\)^\{-1\}
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
!> multiplied by a constant \alpha as reported in R. Zhang, Phys. of Plasmas,
!> vol.22, p.044501, 2015.
!> inputs:
!>   alpha: (real8) multiplicative constant
!>   vec:     (real8)(3) vector to be transformed
!> outputs:
!>   approximated_cayley_transofrm: (real8)(3,3) the approximated Cayley 
!>                                               transfrom of vec
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
