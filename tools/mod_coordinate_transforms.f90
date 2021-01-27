!> Module to calculate coordinate transforms between XYZ and RZPhi coordinate systems.
module mod_coordinate_transforms
  implicit none
  private
  public cartesian_to_cylindrical
  public cylindrical_to_cartesian
  public vector_cartesian_to_cylindrical
  public vector_cylindrical_to_cartesian
  public vector_rotation
contains
  !> convert a position in xyz coordinates to RZPhi coordinates
  pure function cartesian_to_cylindrical(xyz) result(cyl)
    real*8, intent(in)           :: xyz(3) !< The position in xyz coordinates
    real*8                       :: cyl(3) !< The position in RZPhi coordinates

    cyl(1) = sqrt(xyz(1)*xyz(1) + xyz(2)*xyz(2))
    cyl(2) = xyz(3)
    cyl(3) = atan2(-xyz(2), xyz(1))
  end function cartesian_to_cylindrical

  !> convert a position in RZPhi coordinates to xyz coordinates
  pure function cylindrical_to_cartesian(cyl) result(xyz)
    real*8, intent(in)           :: cyl(3) !< The position in RZPhi coordinates
    real*8                       :: xyz(3) !< The position in xyz coordinates

    xyz(1) =  cyl(1)*cos(cyl(3))
    xyz(2) = -cyl(1)*sin(cyl(3))
    xyz(3) =  cyl(2)
  end function cylindrical_to_cartesian

  !> convert a vector in (ex,ey,ez) basis into (eR,eZ,ephi) basis
  pure function vector_cartesian_to_cylindrical(phi,a) result(b)
    real*8, intent(in)               :: phi !< The local toroidal angle
    real*8, dimension(3), intent(in) :: a   !< The vector components in (ex,ey,ez) basis
    real*8, dimension(3)             :: b   !< The vector components in (eR,eZ,ephi) basis
    real*8, dimension(2)             :: sincosphi

    sincosphi = (/sin(phi),cos(phi)/)
   
    b(1) = a(1)*sincosphi(2) - a(2)*sincosphi(1) 
    b(2) = a(3)
    b(3) = -1.d0*(a(1)*sincosphi(1) + a(2)*sincosphi(2))
  end function vector_cartesian_to_cylindrical  

  !> convert a vector in (eR,eZ,ephi) basis into (ex,ey,ez)  basis
  pure function vector_cylindrical_to_cartesian(phi,a) result(b)
    real*8, intent(in)               :: phi !< The local toroidal angle
    real*8, dimension(3), intent(in) :: a   !< The vector components in (eR,eZ,ephi) basis
    real*8, dimension(3)             :: b   !< The vector components in (ex,ey,ez) basis
    real*8, dimension(2)             :: sincosphi

    sincosphi = (/sin(phi),cos(phi)/)

    b(1) = a(1)*sincosphi(2) - a(3)*sincosphi(1) 
    b(2) = -1.d0*(a(1)*sincosphi(1) + a(3)*sincosphi(2))
    b(3) = a(2)
  end function vector_cylindrical_to_cartesian

  !> multiply a vector with the rotation matrix for angle phi around the z-axis in RZPhi coordinates
  pure function vector_rotation(in, phi) result(out)
    real*8, intent(in) :: in(3) !< Input vector in RZPhi coordinates
    real*8, intent(in) :: phi
    real*8             :: out(3) !< Output vector in RZPhi coordinates

    out(1) = cos(phi) * in(1) - sin(phi) * in(3)
    out(2) = in(2)
    out(3) = sin(phi) * in(1) + cos(phi) * in(3)
  end function vector_rotation
end module mod_coordinate_transforms
