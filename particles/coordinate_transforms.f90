module coordinate_transforms
  private

  public :: XYZtoRZPhi
  public :: RZPhiToXYZ
contains
  !> This function converts a vector in xyz coordinates to RZPhi coordinates
  pure function XYZtoRZPhi(xyz, origin) result(RZPhi)
    real*8, intent(in)           :: xyz(3) !< The vector components in xyz coordinates
    real*8                       :: RZPhi(3) !< The vector components in RZPhi coordinates
    real*8, intent(in), optional :: origin(3) !< The xyz coordinates of the base of the vector, assumed 0 if omitted

    real*8 :: phi, dp(2,2) ! Angle and dot products
    if (present(origin) .and. origin(2) .ne. 0.d0 .and. origin(1) .ne. 0.d0) then ! Catch atan2(0,0)
      phi = atan2(origin(2),origin(1)) ! Calculate angle to find inner products between e_x . e_r etc
    else ! If we are at the origin then lim x -> 0 atan2(0,x) = 0
      phi = 0.d0
    end if
    dp = dot_products(phi)

    RZPhi(1) = xyz(1)*dp(1,1)+xyz(2)*dp(1,2)
    RZPhi(2) = xyz(3)
    RZPhi(3) = xyz(1)*dp(2,1)+xyz(2)*dp(2,2)
  end function XYZtoRZPhi
  !> This function converts a position in RZPhi coordinates to xyz coordinates
  pure function RZPhiToXYZ(RZPhi, origin) result(xyz)
    real*8, intent(in)           :: RZPhi(3) !< The vector components in RZPhi coordinates
    real*8                       :: xyz(3) !< The vector components in xyz coordinates
    real*8, intent(in), optional :: origin(3) !< The RZPhi coordinates of the base of the vector, assumed 0 if omitted

    real*8 :: phi, dp(2,2) ! Angle and dot products
    if (present(origin)) then
      phi = origin(3)
    else ! Else we are at the origin so phi = 0
      phi = 0.d0
    end if
    dp = dot_products(phi)

    xyz(1) = RZPhi(1)*dp(1,1)+RZPhi(3)*dp(1,2)
    xyz(2) = RZPhi(1)*dp(2,1)+RZPhi(3)*dp(2,2)
    xyz(3) = RZPhi(2)
  end function RZPhiToXYZ
  !> This function contains the dot products in the JOREK coordinate
  !! system between the basis vectors e_x,e_y and e_r,e_theta
  pure function dot_products(phi) result(dp)
    real*8, intent(in) :: phi
    real*8             :: dp(2,2), cp, sp
    cp = cos(phi)
    sp = sin(phi)
    dp = reshape((/cp,-sp,-sp,-cp/),(/2,2/)) ! Since this is symmetric row/column major does not matter
  end function dot_products
end module coordinate_transforms
