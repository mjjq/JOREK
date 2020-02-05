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
public coordinate_transform_RZPHI_to_XYZ, coordinate_transform_XYZ_to_RZPHI
public transform_derivatives_st_to_RZ,transform_derivatives_RZ_to_st
public transform_first_derivatives_st_to_RZ
public transform_second_derivatives_st_to_RZ
!> interfaces

!> overload transform derivatives from local to global coordinates
interface transform_derivatives_st_to_RZ
   module procedure transform_first_derivatives_st_to_RZ ,&
        transform_second_derivatives_st_to_RZ
end interface transform_derivatives_st_to_RZ

!> overload transfrom derivatives from global to local coorsinates
interface transform_derivatives_RZ_to_st
   module procedure transform_first_derivatives_RZ_to_st,&
        transform_second_derivatives_RZ_to_st
end interface transform_derivatives_RZ_to_st
  
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

!> This function computes the cartesian coordinate position
!> \{X,Y,Z\} from cylindrical one \{R,Z,\phi\}.
!> Be careful: positive angle in \{X,Y,Z\} are counter-clockwise
!> while in \{R,Z,\phi\} are clockwise hence -phi is used.
!> inputs:
!>   a: (real8)(3) position in \{R,Z,\phi\} coordinates 
!> outputs:
!>   b: (real8)(3) position in \{X,Y,Z} coordinates
pure function coordinate_transform_RZPHI_to_XYZ(a) result(b)
  ! declare input variables
  real(kind=8),dimension(3),intent(in) :: a
  ! declare output variables
  real(kind=8),dimension(3) :: b

  b(1) = a(1)*cos(a(3))
  b(2) = -1.d0*a(1)*sin(a(3)) !< negative sign: clockwise -> counter-clockwise
  b(3) = a(2)

end function coordinate_transform_RZPHI_to_XYZ

!---------------------------------------------------------------------------

!> This function computes the cartesian coordinate position
!> \{R,Z,\phi\} from cylindrical one \{X,Y,Z}.
!> Be careful: positive angle in \{X,Y,Z\} are counter-clockwise
!> while in \{R,Z,\phi\} are clockwise hence -phi is used.
!> inputs:
!>   a: (real8)(3) position in \{X,Y,Z} coordinates 
!> outputs:
!>   b: (real8)(3) position in \{R,Z,\phi\} coordinates
pure function coordinate_transform_XYZ_to_RZPHI(a) result(b)
  ! declare input variables
  real(kind=8),dimension(3),intent(in) :: a
  ! declare output variables
  real(kind=8),dimension(3) :: b

  b(1) = sqrt(a(1)*a(1)+a(2)*a(2))
  b(2) = a(3)
  b(3) = atan2(-a(2),a(1)) !< negative sign: counter-clockwise -> clockwise

end function coordinate_transform_XYZ_to_RZPHI

!---------------------------------------------------------------------------

!> This function rotates a vector from a \{R,Z,\phi\} basis to a \{X,Y,Z\}
!> This is done performing a first swap and reflaction of the Y axis
!> and then, a rotation of angle phi (clockwise)
!> inputs:
!>   phi: (real8) angle positive-clockwise
!>   a:   (real8)(3) vector in the reference \{R,Z,\phi\} 
!> outputs:
!>   b: (real8)(3) vector in the reference \{X,Y,Z\}
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

!> This function rotates a vector from a \{X,Y,Z\} basis to a \{R,Z,\phi\}
!> inputs:
!>   phi: (real8) angle positive-clockwise
!>   a:   (real8)(3) vector in the reference \{X,Y,Z\} 
!> outputs:
!>   b: (real8)(3) vector in the reference \{R,Z,\phi\}
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

!--------------------------------------------------------------------------
!> This procedure expresses first order derivatives from the local (s,t)
!> to the global (R,Z) reference system.
!> inputs:
!>   n_v: (integer) number of physical quantities
!>   P_s: (real8)(n_v) physcial quantitiy first derivative in s
!>   P_t: (real8)(n_v) physcial quantity first derivatives in t
!>   R_s: (real8) major radius first derivative in s
!>   R_z: (real8) major radius first derivative in t
!>   Z_s: (real8) vertical position first derivative in s
!>   Z_t: (real8) vertical position first derivative in t
!> outputs:
!>   P_R: (real8)(n_v) physical quantity first derivative in R
!>   P_Z: (real8)(n_v) physcial quantity first derivative in Z
pure subroutine transform_first_derivatives_st_to_RZ(P_R,P_Z,n_v,&
     P_s,P_t,R_s,R_t,Z_s,Z_t)
  implicit none
  !> define input variables
  integer,intent(in) :: n_v
  real(kind=8),dimension(n_v),intent(in) :: P_s,P_t
  real(kind=8),intent(in) :: R_s,R_t,Z_s,Z_t
  !> define output variables
  real(kind=8),dimension(n_v),intent(out) :: P_R,P_Z
  !> internal variables
  real(kind=8) :: inverse_jacobian

  !> compute the inverse of the matrix jacobian
  inverse_jacobian = 1.d0/(R_s*Z_t-R_t*Z_s)
  !> compute transform derivatives in global coordinates
  P_R = (P_s*Z_t-P_t*Z_s)*inverse_jacobian!< R
  P_Z = (P_t*R_s-P_s*R_t)*inverse_jacobian!< Z
  
end subroutine transform_first_derivatives_st_to_RZ

!---------------------------------------------------------------------------

!> This procedure expresses second order derivatives from the local (s,t)
!> to the global (R,Z) reference system. Functions are assumed to be Hessian
!> inputs:
!>   n_v:  (integer) number of physical quantities
!>   P_ss: (real8)(n_v) physical quantity second derivatives in s
!>   P_st: (real8)(n_v) physcal quanitiy cross derivatives in s,t
!>   P_tt: (real8)(n_v) physical quantitiy second derivative in t
!>   P_R:  (real8)(n_v) physcal quantitiy first derivative in R
!>   P_Z:  (real8)(n_v) phyisical quantitiy first derivative in Z
!>   R_s:  (real8) major radius first derivative in s
!>   R_t:  (real8) major radius first derivative in t
!>   R_ss: (real8) major radius second derivative in s
!>   R_st: (real8) major radius cross derivative in s,t
!>   R_tt: (real8) major radius second derivative in t)
!>   Z_s:  (real8) vertical position first derivative in s
!>   Z_t:  (real8) vertical position first derivative in t
!>   Z_ss: (real8) vertical position second derivative in s
!>   Z_st: (real8) vertical position cross derivative in s,t
!>   Z_tt: (real8) vertical position second derivative in t
!> outputs:
!>   P_RR: (real8)(n_v) physical quantitiy second derivatives in R
!>   P_RZ: (real8)(n_v) physical quantity cross derivatives in R,Z
!>   P_ZZ: (real8)(n_v) physical quantitiy second derivatives in Z
pure subroutine transform_second_derivatives_st_to_RZ(P_RR,P_RZ,P_ZZ,n_v,&
     P_ss,P_st,P_tt,P_R,P_Z,R_s,R_t,R_ss,R_st,R_tt,Z_s,Z_t,Z_ss,Z_st,Z_tt)
  implicit none
  !> declare input variables
  integer,intent(in) :: n_v
  real(kind=8),dimension(n_v),intent(in) :: P_ss,P_st,P_tt,P_R,P_Z
  real(kind=8),intent(in) :: R_s,R_t,R_ss,R_st,R_tt
  real(kind=8),intent(in) :: Z_s,Z_t,Z_ss,Z_st,Z_tt
  !> declare output variables
  real(kind=8),dimension(n_v),intent(out) :: P_RR,P_RZ,P_ZZ
  !> declare input variables
  integer :: i !< index
  real(kind=8),dimension(10) :: transformation_matrix !< 10:jacobian
  real(kind=8),dimension(n_v) :: RHS_RR,RHS_RZ,RHS_ZZ

  !< compute matrix elements
  transformation_matrix(1:9) = [R_s*R_s,2.d0*R_s*Z_s,Z_s*Z_s,&
       R_t*R_s,R_s*Z_t+Z_s*R_t,Z_s*Z_t,R_t*R_t,&
       2.d0*R_t*Z_t,Z_t*Z_t]
  !< compute the inverse of the matrix jacobian
  transformation_matrix(10) = 1.d0/(transformation_matrix(1)*(&
       transformation_matrix(5)*transformation_matrix(9)-&
       transformation_matrix(8)*transformation_matrix(6))+&
       transformation_matrix(2)*(transformation_matrix(7)*&
       transformation_matrix(6)-transformation_matrix(4)*&
       transformation_matrix(9))+transformation_matrix(3)*(&
       transformation_matrix(4)*transformation_matrix(8)-&
       transformation_matrix(5)*transformation_matrix(7)))
  !> compute RHS
  RHS_RR = P_ss - R_ss*P_R - Z_ss*P_Z !< RR
  RHS_RZ = P_st - R_st*P_R - Z_st*P_Z !< RZ
  RHS_ZZ = P_tt - R_tt*P_R - Z_tt*P_Z !< ZZ
  !> compute second order derivatives
  P_RR = ((transformation_matrix(5)*transformation_matrix(9)-&
       transformation_matrix(8)*transformation_matrix(6))*RHS_RR +&
       (transformation_matrix(8)*transformation_matrix(3)-&
       transformation_matrix(2)*transformation_matrix(9))*RHS_RZ +&
       (transformation_matrix(2)*transformation_matrix(6)-&
       transformation_matrix(5)*transformation_matrix(3))*&
       RHS_ZZ)*transformation_matrix(10) !< RR
  P_RZ = ((transformation_matrix(7)*transformation_matrix(6)-&
       transformation_matrix(4)*transformation_matrix(9))*RHS_RR +&
       (transformation_matrix(1)*transformation_matrix(9)-&
       transformation_matrix(7)*transformation_matrix(3))*RHS_RZ +&
       (transformation_matrix(4)*transformation_matrix(3)-&
       transformation_matrix(1)*transformation_matrix(6))*&
       RHS_ZZ)*transformation_matrix(10) !< RZ
  P_ZZ = ((transformation_matrix(4)*transformation_matrix(8)-&
       transformation_matrix(5)*transformation_matrix(7))*RHS_RR +&
       (transformation_matrix(2)*transformation_matrix(7)-&
       transformation_matrix(1)*transformation_matrix(8))*RHS_RZ +&
       (transformation_matrix(1)*transformation_matrix(5)-&
       transformation_matrix(2)*transformation_matrix(4))*&
       RHS_ZZ)*transformation_matrix(10) !< ZZ
  
end subroutine transform_second_derivatives_st_to_RZ

!---------------------------------------------------------------------------------------

!> This procedure transform back the first order derivatives from global RZ coordinates
!> to local s,t coordinates. This method is used mainly for testing.
!> inputs:
!>   P_R,P_Z: (real8)(n_v) physcial quantity first derivatives in R and t
!>   R_s,R_t: (real8) major radius first derivative in s and t
!>   Z_s,Z_t: (real) vertical position first and second order derivatives in s and t
!> outputs:
!>   P_s,P_t: (real8)(n_v) physical quantitiy first derivatives in s and t
pure subroutine transform_first_derivatives_RZ_to_st(P_s,P_t,n_v,P_R,P_Z,&
     R_s,R_t,Z_s,Z_t)
  implicit none
  !> define input variables
  integer,intent(in) :: n_v
  real(kind=8),dimension(n_v),intent(in) :: P_R,P_Z
  real(kind=8),intent(in) :: R_s,R_t,Z_s,Z_t
  !> define output varibales
  real(kind=8),dimension(n_v),intent(out) :: P_s,P_t

  !> compute derivatives
  P_s = P_R*R_s + P_Z*Z_s !< ds
  P_t = P_R*R_t + P_Z*Z_t !< dt
  
end subroutine transform_first_derivatives_RZ_to_st
!---------------------------------------------------------------------------------------

!> This procedure transform back the second order derivatives from globa RZ coordinates
!> to local s,t coordinates. This method is used mainly for testing.
!> inputs:
!>   n_v:            (integer) number fo physical quantities
!>   P_R,P_Z:        (real8)(n_v) physical quantitiy first derivatives
!>                   in R and Z
!>   P_RR,P_RZ,P_ZZ: (real8)(n_v) physical quantitiy second and cross
!>                   dervatives in R and Z
!>   R_s,R_t:        (real8) major radius first derivative in s and t
!>   Z_s,Z_t:        (real8) vertical position first derivatives in s and t
!>   R_ss,R_st,R_tt: (real8) major radius second and cross derivatives in s and t
!>   Z_ss,Z_st,Z_tt: (real8) vertical position second and cross derivatives in s and t
!> outputs:
!>   P_ss,P_st,P_tt: (real)(n_v) physical quantity second and cross
!>                   derivatives in s and t
pure subroutine transform_second_derivatives_RZ_to_st(P_ss,P_st,P_tt,n_v,P_R,P_Z,&
     P_RR,P_RZ,P_ZZ,R_s,R_t,R_ss,R_st,R_tt,Z_s,Z_t,Z_ss,Z_st,Z_tt)
  implicit none
  !> declare input variables
  integer,intent(in) :: n_v
  real(kind=8),dimension(n_v),intent(in) :: P_RR,P_RZ,P_ZZ,P_R,P_Z
  real(kind=8),intent(in) :: R_s,R_t,R_ss,R_st,R_tt
  real(kind=8),intent(in) :: Z_s,Z_t,Z_ss,Z_st,Z_tt
  !> delcare output varibales
  real(kind=8),dimension(n_v),intent(out) :: P_ss,P_st,P_tt

  !> compute second order derivatives
  P_ss = P_RR*R_s*R_s + 2.d0*P_RZ*R_s*Z_s + P_ZZ*Z_s*Z_s + P_R*R_ss + P_Z*Z_ss !< dsds
  P_st = P_RR*R_t*R_s + P_RZ*(Z_t*R_s+R_t*Z_s) + P_ZZ*Z_t*Z_s + P_R*R_st + P_Z*Z_st !<dsdt
  P_tt = P_RR*R_t*R_t + 2.d0*P_RZ*R_t*Z_t + P_ZZ*Z_t*Z_t + P_R*R_tt + P_Z*Z_tt !< dtdt
  
end subroutine transform_second_derivatives_RZ_to_st


!---------------------------------------------------------------------------------------

end module mod_pusher_tools
