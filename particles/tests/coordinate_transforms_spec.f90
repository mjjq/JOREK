!> Test the coordinate transform functions between cylindrical and
!> cartesian coordinates
module coordinate_transforms_spec
use mod_coordinate_transforms, only: cylindrical_to_cartesian, &
                                     cartesian_to_cylindrical, &
                                     vector_rotation
use mod_constants, only: PI
use fruit
implicit none

contains

subroutine test_cylindrical_to_cartesian
  real*8, parameter :: tolerance = 1.d-15
  real*8, dimension(3) :: xyz, cyl

  cyl = [1.d0, 0.d0, 2.d0]
  xyz = cylindrical_to_cartesian(cyl)
  call assert_equals(1.d0, xyz(1), tolerance, "r nonzero phi zero")
  call assert_equals(0.d0, xyz(2), "r nonzero phi zero => y zero")
  call assert_equals(2.d0, xyz(3), "z unchanged")

  cyl = [1.d0, PI/2, 2.d0]
  xyz = cylindrical_to_cartesian(cyl)
  call assert_equals(0.d0, xyz(1), tolerance, "r nonzero phi pi/2")
  call assert_equals(1.d0, xyz(2), tolerance, "r nonzero phi pi => y zero")
  call assert_equals(2.d0, xyz(3), "z unchanged")

  cyl = [1.d0, PI, 2.d0]
  xyz = cylindrical_to_cartesian(cyl)
  call assert_equals(-1.d0, xyz(1), tolerance, "r nonzero phi pi")
  call assert_equals(0.d0, xyz(2), tolerance, "r nonzero phi pi => y zero")
  call assert_equals(2.d0, xyz(3), "z unchanged")

  cyl = [0.d0, 0.d0, 2.d0]
  xyz = cylindrical_to_cartesian(cyl)
  call assert_equals(0.d0, xyz(1), tolerance, "r zero phi zero")
  call assert_equals(0.d0, xyz(2), tolerance, "r zero phi zero => y zero")
  call assert_equals(2.d0, xyz(3), "z unchanged")

  cyl = [0.d0, PI, 2.d0]
  xyz = cylindrical_to_cartesian(cyl)
  call assert_equals(0.d0, xyz(1), tolerance, "r zero phi pi")
  call assert_equals(0.d0, xyz(2), "r zero phi pi => y zero")
  call assert_equals(2.d0, xyz(3), "z unchanged")
end subroutine test_cylindrical_to_cartesian

subroutine test_cartesian_to_cylindrical
  real*8, parameter :: tolerance = 1.d-15
  real*8, dimension(3) :: xyz, cyl

  xyz = [1.d0, 0.d0, 2.d0]
  cyl = cartesian_to_cylindrical(xyz)
  call assert_equals(1.d0, cyl(1), tolerance, "correct radius")
  call assert_equals(0.d0, cyl(2), "correct angle treatment")
  call assert_equals(2.d0, cyl(3), "z unchanged")

  xyz = [0.d0, 0.d0, 2.d0]
  cyl = cartesian_to_cylindrical(xyz)
  call assert_equals(0.d0, cyl(1), tolerance, "r zero phi pi")
  call assert_equals(0.d0, cyl(2), tolerance, "r zero => phi zero")
  call assert_equals(2.d0, cyl(3), "z unchanged")

  xyz = [1.d0, 1.d0, 2.d0]
  cyl = cartesian_to_cylindrical(xyz)
  call assert_equals(sqrt(2.d0), cyl(1), tolerance, "pythagoras => sqrt(2)")
  call assert_equals(PI/4, cyl(2), "angle pi/4")
  call assert_equals(2.d0, cyl(3), "z unchanged")

  xyz = [-1.d0, 0.d0, 2.d0]
  cyl = cartesian_to_cylindrical(xyz)
  call assert_equals(1.d0, cyl(1), tolerance, "r correct for negative values")
  call assert_equals(PI, cyl(2), "angle should be PI or -PI")
  call assert_equals(2.d0, cyl(3), "z unchanged")
end subroutine test_cartesian_to_cylindrical

subroutine test_vector_rotation
  real*8, parameter :: tolerance = 1d-15
  real*8, dimension(3) :: out

  out = vector_rotation([1.d0, 3.d0, 3.d0], 0.d0)
  call assert_equals(1.d0, out(1), tolerance, "no rotation keeps vector x")
  call assert_equals(3.d0, out(2), tolerance, "no rotation keeps vector y")
  call assert_equals(3.d0, out(3), tolerance, "keeps vector z")

  out = vector_rotation([1.d0, 3.d0, 3.d0], PI)
  call assert_equals(-1.d0, out(1), tolerance, "reverses vector x")
  call assert_equals(-3.d0, out(2), tolerance, "reverses vector y")
  call assert_equals(3.d0, out(3), tolerance, "keeps vector z")

  out = vector_rotation([1.d0, 3.d0, -2.d0], PI/2)
  call assert_equals(3.d0, out(1), tolerance, "rotates y to x")
  call assert_equals(-1.d0, out(2), tolerance, "rotates x to -y")
  call assert_equals(-2.d0, out(3), tolerance, "keeps z")

  out = vector_rotation([1.d0, 3.d0, -2.d0], 5*PI/2)
  call assert_equals(3.d0, out(1), tolerance, "rotates y to x")
  call assert_equals(-1.d0, out(2), tolerance, "rotates x to -y")
  call assert_equals(-2.d0, out(3), tolerance, "keeps z")
end subroutine test_vector_rotation
end module coordinate_transforms_spec
