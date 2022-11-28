!> mod_coordinate_transformrs_test contains variables and
!> procedures for testing the mod_coordinate_transforms
!> module procedures
module mod_coordinate_transforms_test
use constants, only: PI,TWOPI
use fruit
implicit none

private
public :: run_fruit_coordinate_transforms

!> Variables--------------------------------------------
integer,parameter :: n_points=4                !< number of test positions
integer,parameter :: n_origins=4               !< number of sphere origins
real*8,parameter  :: tol_r8=2.5d-14            !< tolerance double
real*4,parameter  :: tol_r4=real(7.5d-6,kind=4) !< tolerance float
!> tolerance for calculatons
real*4,parameter  :: tol_calc_r4=real(5.0d-4,kind=4)
real*8,parameter  :: tol_calc_r8=5.0d-12
real*4,parameter  :: zero_r4=real(0.d0,kind=4)
real*4,parameter  :: one_r4=real(1.d0,kind=4)
real*4,dimension(3),parameter :: zeros_r4=real((/0.d0,0.d0,0.d0/),kind=4)
real*8,dimension(3),parameter :: zeros_r8=(/0.d0,0.d0,0.d0/)
!> intervals for randomly chosing the first, second and third
!> position components
real*8,dimension(3),parameter :: x_lowbnd=(/-2.3d1,-3.2d2,-9.d-1/)
real*8,dimension(3),parameter :: x_uppbnd=(/4.23d2,1.45d1,7.50d1/)
real*8,dimension(2),parameter :: r_interval=(/3.23d-3,4.53d2/)
real*8,dimension(2),parameter :: theta_interval=(/0d0,PI/)
real*8,dimension(2),parameter :: phi_interval=(/0.d0,TWOPI/)
!> intervals for randomly chosing the first, second and third
!> vector components
real*8,dimension(3),parameter :: a_lowbnd=(/-3.41d2,-4.67d1,-9.35d1/)
real*8,dimension(3),parameter :: a_uppbnd=(/6.75d1,8.70d1,2.43d2/)
real*4,dimension(3,n_points)  :: x_r4           !< set o positions
real*4,dimension(3,n_origins) :: origin_r4      !< set of origins
real*4,dimension(3,n_origins) :: v1_r4,v2_r4    !< random vectors
real*4,dimension(3,n_origins) :: T_r4,N_r4,B_r4 !< sphere directions 
real*4,dimension(3,n_points)  :: rthetaphi_r4   !< spherical coordinates
real*8,dimension(3,n_points)  :: x_r8           !< set o positions
real*8,dimension(n_points)    :: phi_r8         !< set of toroidal angles
real*8,dimension(3,n_origins) :: origin_r8      !< set of origins
real*8,dimension(3,n_origins) :: v1_r8,v2_r8    !< random vectors
real*8,dimension(3,n_origins) :: T_r8,N_r8,B_r8 !< sphere directioins
real*8,dimension(3,n_points)  :: rthetaphi_r8   !< spherical coordinates

!> Interfaces ------------------------------------------
!> function for testing basis orthonormality
interface test_orthonormality_basis
  module procedure test_orthonormality_basis_r4
  module procedure test_orthonormality_basis_r8
end interface test_orthonormality_basis

contains

!> Fruit test basket -----------------------------------
!> Test basket for executing set-up, tests and tear-down
subroutine run_fruit_coordinate_transforms()
  implicit none
  write(*,'(/A)') "  ... setting-up: coordinate transfroms tests"
  call setup
  write(*,'(/A)') "  ... running: coordinate transforms tests"
  call test_cartesian_tofrom_cylindrical_transform
  call test_cartesian_tofrom_spherical_latitude_transform
  call test_cartesian_tofrom_spherical_colatitude_std_transform
  call test_cartesian_tofrom_cylindrical_vector_rotation
  call test_vectors_to_orthonormal_basis
  call test_cartesian_tofrom_spherical
  write(*,'(/A)') "  ... tearing-down: coordinate transforms tests"
  call teardown
end subroutine run_fruit_coordinate_transforms

!> Set-up and tear-down --------------------------------
!> Set-up test features common to all tests
subroutine setup()
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_math_operators, only: cross_product
  use mod_sampling,       only: sample_uniform_sphere_corona_rthetaphi
  implicit none
  !> variables
  integer             :: ii
  real*8,dimension(2) :: r3_interval,costheta_interval
  real*8,dimension(3) :: rnd_r8

  !> generate set of random toroidal angles
  call gnu_rng_interval(n_points,phi_interval,phi_r8)
  !> generate random spherical points
  r3_interval = r_interval**3; costheta_interval = cos(theta_interval);
  do ii=1,n_points
    call random_number(rnd_r8)
    rthetaphi_r8(:,ii) = sample_uniform_sphere_corona_rthetaphi(&
    r3_interval,costheta_interval,phi_interval,rnd_r8)
  enddo
  !> generate random positions (assume cartesian coord.)
  do ii=1,n_points
    call gnu_rng_interval(n_points,x_lowbnd,x_uppbnd,x_r8(:,ii))
  enddo
  !> generate random origins (assume cartesian coord.)
  do ii=1,n_origins
    call gnu_rng_interval(n_points,x_lowbnd,x_uppbnd,origin_r8(:,ii))
    !> generate random orthonormal directions (assume cartesian coord.)
    call gnu_rng_interval(n_points,a_lowbnd,a_uppbnd,v1_r8(:,ii))
    call gnu_rng_interval(n_points,a_lowbnd,a_uppbnd,v2_r8(:,ii))
    !> compute the vector for the ith origin
    v1_r8(:,ii) = v1_r8(:,ii)-origin_r8(:,ii)
    v2_r8(:,ii) = v2_r8(:,ii)-origin_r8(:,ii)
    !> compute basis
    T_r8(:,ii) = v1_r8(:,ii)/norm2(v1_r8(:,ii))
    N_r8(:,ii) = v2_r8(:,ii) - (dot_product(T_r8(:,ii),v2_r8(:,ii)))*T_r8(:,ii)
    N_r8(:,ii) = N_r8(:,ii)/norm2(N_r8(:,ii))
    B_r8(:,ii) = cross_product(T_r8(:,ii),N_r8(:,ii))
    B_r8(:,ii) = B_r8(:,ii)/norm2(B_r8(:,ii))
    !> test the correct generation of the orthonormal basis
    call test_orthonormality_basis(T_r8(:,ii),N_r8(:,ii),B_r8(:,ii))
  enddo

  !> convert to float precision
  x_r4 = real(x_r8,kind=4); origin_r4 = real(origin_r8,kind=8); 
  v1_r4 = real(v1_r8,kind=4); v2_r4 = real(v2_r8,kind=4); T_r4 = real(T_r8,kind=4); 
  N_r4 = real(N_r8,kind=4); B_r4 = real(B_r8,kind=4); 
  rthetaphi_r4 = real(rthetaphi_r8,kind=4);
  do ii=1,n_origins
    call test_orthonormality_basis(T_r4(:,ii),N_r4(:,ii),B_r4(:,ii))
  enddo
end subroutine setup

!> Clean-up all common test features
subroutine teardown()
  implicit none
  !> set all variables to 0
  x_r4 = zero_r4; origin_r4 = zero_r4; v1_r4 = zero_r4;
  v2_r4 = zero_r4; T_r4 = zero_r4; N_r4 = zero_r4; B_r4 = zero_r4;
  rthetaphi_r4 = zero_r4; x_r8 = 0.d0; origin_r8 = 0.d0; 
  v1_r8 = 0.d0; v2_r8 = 0.e0; T_r8 = 0.d0; N_r8 = 0.d0; B_r8 = 0.d0;
  phi_r8 = 0.d0; rthetaphi_r8 = 0.d0;
end subroutine teardown
!> Tests -----------------------------------------------

!> Test cartesian to spherical and spherical to cartesian
!> for both sing and double precision
subroutine test_cartesian_tofrom_spherical()
  use mod_coordinate_transforms, only: vectors_cartesian_to_spherical
  use mod_coordinate_transforms, only: vectors_spherical_to_cartesian
  implicit none
  !> variables
  integer                      :: ii
  real*4,dimension(3,n_points) :: x_cart_new_r4,zeros_3nv_r4
  real*8,dimension(3,n_points) :: x_cart_new_r8,zeros_3nv_r8
  !> test single precision
  zeros_3nv_r4 = real(0d0,kind=4);
  do ii=1,n_points
    x_cart_new_r4 = x_r4;
    call vectors_cartesian_to_spherical(n_points,rthetaphi_r4(:,ii),x_cart_new_r4)
    call vectors_spherical_to_cartesian(n_points,rthetaphi_r4(:,ii),x_cart_new_r4)
    call assert_equals(x_r4-x_cart_new_r4,zeros_3nv_r4,3,n_points,tol_calc_r4,&
    "Error test cartesian to/from spherical (float): x-cartesian mismatch!")
  enddo
  zeros_3nv_r8 = 0d0;
  do ii=1,n_points
    x_cart_new_r8 = x_r8;
    call vectors_cartesian_to_spherical(n_points,rthetaphi_r8(:,ii),x_cart_new_r8)
    call vectors_spherical_to_cartesian(n_points,rthetaphi_r8(:,ii),x_cart_new_r8)
    call assert_equals(x_r8-x_cart_new_r8,zeros_3nv_r8,3,n_points,tol_calc_r8,&
    "Error test cartesian to/from spherical (double): x-cartesian mismatch!")
  enddo
end subroutine test_cartesian_tofrom_spherical

!> Test cartesian to cylindrical and cylindrical to 
!> cartesian transformations for single and double
!> precision functions
subroutine test_cartesian_tofrom_cylindrical_transform()
  use mod_coordinate_transforms, only: cartesian_to_cylindrical
  use mod_coordinate_transforms, only: cylindrical_to_cartesian
  implicit none
  !> variables
  integer             :: ii
  real*4,dimension(3) :: x_cart_new_r4,x_cyl_r4
  real*8,dimension(3) :: x_cart_new_r8,x_cyl_r8
  !> test single procision
  do ii=1,n_points
    x_cyl_r4 = cartesian_to_cylindrical(x_r4(:,ii))
    x_cart_new_r4 = cylindrical_to_cartesian(x_cyl_r4)   
    call assert_equals(x_r4(:,ii)-x_cart_new_r4,zeros_r4,3,tol_calc_r4,&
    "Error test cartesian to/from cylindrical (float): x-cartesian mismatch!")
  enddo
  !> test double procision
  do ii=1,n_points
    x_cyl_r8 = cartesian_to_cylindrical(x_r8(:,ii))
    x_cart_new_r8 = cylindrical_to_cartesian(x_cyl_r8)
    call assert_equals(x_r8(:,ii)-x_cart_new_r8,zeros_r8,3,tol_calc_r8,&
    "Error test cartesian to/from cylindrical (double): x-cartesian mismatch!")
  enddo
end subroutine test_cartesian_tofrom_cylindrical_transform

!> Test cartesian to spherical (latitude) and spherical
!> (latitude) transformations for single and double
!> precision functions
subroutine test_cartesian_tofrom_spherical_latitude_transform()
  use mod_coordinate_transforms, only: cartesian_to_spherical_latitude
  use mod_coordinate_transforms, only: spherical_latitude_to_cartesian
  implicit none
  !> variables
  integer             :: ii,jj
  real*4,dimension(3) :: x_cart_new_r4,rpsichi_r4
  real*8,dimension(3) :: x_cart_new_r8,rpsichi_r8
  !> test single precision
  do jj=1,n_origins
    do ii=1,n_points
      rpsichi_r4 = cartesian_to_spherical_latitude(x_r4(:,ii),&
      origin_r4(:,jj),T_r4(:,jj),N_r4(:,jj),B_r4(:,jj))
      x_cart_new_r4 = spherical_latitude_to_cartesian(rpsichi_r4,&
      origin_r4(:,jj),T_r4(:,jj),N_r4(:,jj),B_r4(:,jj))
      call assert_equals(x_r4(:,ii)-x_cart_new_r4,zeros_r4,3,tol_calc_r4,&
      "Error test cartesian to/from spherical-latitude (float): x-cartesian mismatch!") 
    enddo
  end do
  !> test double precision
  do jj=1,n_origins
    do ii=1,n_points
      rpsichi_r8 = cartesian_to_spherical_latitude(x_r8(:,ii),&
      origin_r8(:,jj),T_r8(:,jj),N_r8(:,jj),B_r8(:,jj))
      x_cart_new_r8 = spherical_latitude_to_cartesian(rpsichi_r8,&
      origin_r8(:,jj),T_r8(:,jj),N_r8(:,jj),B_r8(:,jj))
      call assert_equals(x_r8(:,ii)-x_cart_new_r8,zeros_r8,3,tol_calc_r8,&
      "Error test cartesian to/from spherical-latitude (double): x-cartesian mismatch!") 
    enddo
  end do
end subroutine test_cartesian_tofrom_spherical_latitude_transform

!> Test cartesian to spherical (colatitude) and spherical
!> (latitude) transformations for single and double
!> precision functions. Standard direction version
subroutine test_cartesian_tofrom_spherical_colatitude_std_transform()
  use mod_coordinate_transforms, only: cartesian_to_spherical_colatitude
  use mod_coordinate_transforms, only: spherical_colatitude_to_cartesian
  implicit none
  !> variables
  integer             :: ii,jj
  real*4,dimension(3) :: x_cart_new_r4,rthetachi_r4
  real*8,dimension(3) :: x_cart_new_r8,rthetachi_r8
  !> test single precision
  do jj=1,n_origins
    do ii=1,n_points
      rthetachi_r4 = cartesian_to_spherical_colatitude(x_r4(:,ii),origin_r4(:,jj))
      x_cart_new_r4 = spherical_colatitude_to_cartesian(rthetachi_r4,origin_r4(:,jj))
      call assert_equals(x_r4(:,ii)-x_cart_new_r4,zeros_r4,3,tol_calc_r4,&
      "Error test cartesian to/from spherical-colatitude standard (float): x-cartesian mismatch!") 
    enddo
  end do
  !> test double precision
  do jj=1,n_origins
    do ii=1,n_points
      rthetachi_r8 = cartesian_to_spherical_colatitude(x_r8(:,ii),origin_r8(:,jj))
      x_cart_new_r8 = spherical_colatitude_to_cartesian(rthetachi_r8,origin_r8(:,jj))
      call assert_equals(x_r8(:,ii)-x_cart_new_r8,zeros_r8,3,tol_calc_r8,&
      "Error test cartesian to/from spherical-colatitude standard (double): x-cartesian mismatch!") 
    enddo
  end do
end subroutine test_cartesian_tofrom_spherical_colatitude_std_transform

!> test generation of orthonormal basis for 3d cartesian coordinates
subroutine test_vectors_to_orthonormal_basis()
  use mod_coordinate_transforms, only: vectors_to_orthonormal_basis
  implicit none
  integer :: ii
  real*4,dimension(3) :: T_new_r4,N_new_r4,B_new_r4
  real*8,dimension(3) :: T_new_r8,N_new_r8,B_new_r8

  !> test orthonormal basis single precision
  do ii=1,n_origins
    call vectors_to_orthonormal_basis(v1_r4(:,ii),v2_r4(:,ii),T_new_r4,N_new_r4,B_new_r4)
    call test_orthonormality_basis(T_new_r4,N_new_r4,B_new_r4)
    call assert_equals(T_r4(:,ii),T_new_r4,3,tol_r4,&
    "Error vectors to orhtonormal basis (float): T basis mismatch!)")
    call assert_equals(N_r4(:,ii),N_new_r4,3,tol_r4,&
    "Error vectors to orhtonormal basis (float): T basis mismatch!)")
    call assert_equals(B_r4(:,ii),B_new_r4,3,tol_r4,&
    "Error vectors to orhtonormal basis (float): B basis mismatch!)")
  enddo
  !> test orthonormal basis double precision
  do ii=1,n_origins
    call vectors_to_orthonormal_basis(v1_r8(:,ii),v2_r8(:,ii),T_new_r8,N_new_r8,B_new_r8)
    call test_orthonormality_basis(T_new_r8,N_new_r8,B_new_r8)
    call assert_equals(T_r8(:,ii),T_new_r8,3,tol_r8,&
    "Error vectors to orhtonormal basis (double): T basis mismatch!)")
    call assert_equals(N_r8(:,ii),N_new_r8,3,tol_r8,&
    "Error vectors to orhtonormal basis (double): T basis mismatch!)")
    call assert_equals(B_r8(:,ii),B_new_r8,3,tol_r8,&
    "Error vectors to orhtonormal basis (double): B basis mismatch!)")
  enddo 

end subroutine test_vectors_to_orthonormal_basis

!> Test vector transformation from cartesian to cylindrical and back
subroutine test_cartesian_tofrom_cylindrical_vector_rotation()
  use mod_coordinate_transforms, only: vector_cartesian_to_cylindrical
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  implicit none
  integer :: ii,jj
  real*8,dimension(3) :: T_cyl_r8,N_cyl_r8,B_cyl_r8
  real*8,dimension(3) :: T_cart_r8,N_cart_r8,B_cart_r8

  do jj=1,n_origins
    do ii=1,n_points
      !> transform cartesian basis to cylindrical 
      !> and check that it is still a basis
      T_cyl_r8 = vector_cartesian_to_cylindrical(phi_r8(ii),T_r8(:,jj))
      N_cyl_r8 = vector_cartesian_to_cylindrical(phi_r8(ii),N_r8(:,jj))
      B_cyl_r8 = vector_cartesian_to_cylindrical(phi_r8(ii),B_r8(:,jj))
      call test_orthonormality_basis(T_cyl_r8,N_cyl_r8,B_cyl_r8)
      !> transform back and check consistency
      T_cart_r8 = vector_cylindrical_to_cartesian(phi_r8(ii),T_cyl_r8)
      N_cart_r8 = vector_cylindrical_to_cartesian(phi_r8(ii),N_cyl_r8)
      B_cart_r8 = vector_cylindrical_to_cartesian(phi_r8(ii),B_cyl_r8)
      !> check correctness
      call assert_equals(T_cart_r8,T_r8(:,jj),3,tol_r8,&
      "Error test vector cartesian to/from cylindrical: T vector mismatch!")
      call assert_equals(N_cart_r8,N_r8(:,jj),3,tol_r8,&
      "Error test vector cartesian to/from cylindrical: N vector mismatch!")
      call assert_equals(B_cart_r8,B_r8(:,jj),3,tol_r8,&
      "Error test vector cartesian to/from cylindrical: B vector mismatch!")
    enddo
  enddo
end subroutine test_cartesian_tofrom_cylindrical_vector_rotation

!> test vector rotation of a toroidal angle
subroutine test_vector_rotation_toroidal_angle()
  use mod_coordinate_transforms, only: vector_rotation
  use mod_coordinate_transforms, only: cartesian_velocity_to_cylindrical
  implicit none
  integer :: ii,jj
  real*8,dimension(3) :: T_rot_fwd_r8,N_rot_fwd_r8,B_rot_fwd_r8
  real*8,dimension(3) :: T_rot_bck_r8,N_rot_bck_r8,B_rot_bck_r8

  do jj=1,n_origins
    do ii=1,n_points
      !> forward rotation of a toroidal angle phi and check orthonormality
      T_rot_fwd_r8 = vector_rotation(T_r8(:,jj),phi_r8(ii))
      N_rot_fwd_r8 = vector_rotation(N_r8(:,jj),phi_r8(ii))
      B_rot_fwd_r8 = vector_rotation(B_r8(:,jj),phi_r8(ii))
      call test_orthonormality_basis(T_rot_fwd_r8,N_rot_fwd_r8,B_rot_fwd_r8)
      !> backward rotation and check
      T_rot_bck_r8 = cartesian_velocity_to_cylindrical(T_rot_fwd_r8,phi_r8(ii))
      N_rot_bck_r8 = cartesian_velocity_to_cylindrical(N_rot_fwd_r8,phi_r8(ii))
      B_rot_bck_r8 = cartesian_velocity_to_cylindrical(B_rot_fwd_r8,phi_r8(ii))
      call assert_equals(T_rot_bck_r8,T_r8(:,jj),3,tol_r8,&
      "Error test vector rotation (toroidal): T vector mismatch!")
      call assert_equals(N_rot_bck_r8,N_r8(:,jj),3,tol_r8,&
      "Error test vector rotation (toroidal): N vector mismatch!")
      call assert_equals(B_rot_bck_r8,B_r8(:,jj),3,tol_r8,&
      "Error test vector rotation (toroidal): B vector mismatch!")
    enddo
  enddo
end subroutine test_vector_rotation_toroidal_angle

!> Tools -----------------------------------------------
!> method for testing the orthonormality 
!> of a basis function. Single precision.
subroutine test_orthonormality_basis_r4(v1,v2,v3)
  implicit none
  real*4,dimension(3),intent(in) :: v1,v2,v3
  call assert_equals(dot_product(v1,v1),one_r4,tol_r4,&
  "Error basis orthonormality (float): v1 is not normalized!")
  call assert_equals(dot_product(v2,v2),one_r4,tol_r4,&
  "Error basis orthonormality (float): v2 is not normalized!")
  call assert_equals(dot_product(v3,v3),one_r4,tol_r4,&
  "Error basis orthonormality (float): v3 is not normalized!")
  call assert_equals(dot_product(v1,v2),zero_r4,tol_r4,&
  "Error basis orthonormality (float): v1 and v2 are not orthogonal!")
  call assert_equals(dot_product(v1,v3),zero_r4,tol_r4,&
  "Error basis orthonormality (float): v1 and v3 are not orthogonal!")
  call assert_equals(dot_product(v2,v3),zero_r4,tol_r4,&
  "Error basis orthonormality (float): v2 and v3 are not orthogonal!")
end subroutine test_orthonormality_basis_r4

!> method for testing the orthonormality 
!> of a basis function. Double precision.
subroutine test_orthonormality_basis_r8(v1,v2,v3)
  implicit none
  real*8,dimension(3),intent(in) :: v1,v2,v3
  call assert_equals(dot_product(v1,v1),1.d0,tol_r8,&
  "Error basis orthonormality (double): v1 is not normalized!")
  call assert_equals(dot_product(v2,v2),1.d0,tol_r8,&
  "Error basis orthonormality (double): v2 is not normalized!")
  call assert_equals(dot_product(v3,v3),1.d0,tol_r8,&
  "Error basis orthonormality (double): v3 is not normalized!")
  call assert_equals(dot_product(v1,v2),0.d0,tol_r8,&
  "Error basis orthonormality (double): v1 and v2 are not orthogonal!")
  call assert_equals(dot_product(v1,v3),0.d0,tol_r8,&
  "Error basis orthonormality (double): v1 and v3 are not orthogonal!")
  call assert_equals(dot_product(v2,v3),0.d0,tol_r8,&
  "Error basis orthonormality (double): v2 and v3 are not orthogonal!")
end subroutine test_orthonormality_basis_r8

!>------------------------------------------------------
end module mod_coordinate_transforms_test
