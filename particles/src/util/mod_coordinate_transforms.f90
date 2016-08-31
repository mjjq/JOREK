!> Module to calculate coordinate transforms between XYZ and cyl coordinate systems.
module mod_coordinate_transforms
  private

  public :: cartesian_to_cylindrical
  public :: cylindrical_to_cartesian
  public :: vector_rotation
contains
  !> convert a position in xyz coordinates to RPhiZ coordinates
  pure function cartesian_to_cylindrical(xyz) result(cyl)
    real*8, intent(in)           :: xyz(3) !< The vector components in xyz coordinates
    real*8                       :: cyl(3) !< The vector components in RPhiZ coordinates

    cyl(1) = sqrt(xyz(1)**2 + xyz(2)**2)
    cyl(2) = atan2(xyz(2), xyz(1))
    cyl(3) = xyz(3)
  end function cartesian_to_cylindrical

  !> converts a position in RPhiZ coordinates to xyz coordinates
  pure function cylindrical_to_cartesian(cyl) result(xyz)
    real*8, intent(in)           :: cyl(3) !< The vector components in RPhiZ coordinates
    real*8                       :: xyz(3) !< The vector components in xyz coordinates

    xyz(1) = cyl(1)*cos(cyl(2))
    xyz(2) = cyl(1)*sin(cyl(2))
    xyz(3) = cyl(3)
  end function cylindrical_to_cartesian

  !> multiply a vector with the rotation matrix for angle phi
  pure function vector_rotation(in, phi) result(out)
    real*8, intent(in) :: in(3)
    real*8, intent(in) :: phi
    real*8             :: out(3)

    out(1:2) = [cos(phi) * in(1) + sin(phi) * in(2), &
               -sin(phi) * in(1) + cos(phi) * in(2)]
    out(3) = in(3)
  end function vector_rotation
end module mod_coordinate_transforms
