!> the module mod_lights_assert_equals_tools contains procedures
!> for comparing various data structures
module mod_assert_equals_tools
use fruit
implicit none

private
public :: assert_equals_rel_error

!> Variables -----------------------------------------------------
!> Interfaces ----------------------------------------------------
interface assert_equals_rel_error
  module procedure assert_equals_rel_error_r8
  module procedure assert_equals_rel_error_1d_r8
  module procedure assert_equals_rel_error_2d_r8
  module procedure assert_equals_rel_error_3d_r8
end interface assert_equals_rel_error

contains

!> Procedures ----------------------------------------------------
!> assert equals for 0D-real8 arrays with relative error
subroutine assert_equals_rel_error_r8(val_1,val_2,tol,message)
  implicit none
  !> inputs
  real*8,intent(in) :: val_1,val_2,tol
  character(len=*),intent(in) :: message
  !> variables
  integer :: ii
  real*8  :: error
  !> compare with relative error
  if(val_2.ne.0.d0) error = abs((val_1-val_2)/val_2)
  call assert_equals(error,0.d0,tol,message)
end subroutine assert_equals_rel_error_r8

!> assert equals for 1D-real8 arrays with relative error
subroutine assert_equals_rel_error_1d_r8(size_1,arr_1,&
arr_2,tol,message)
  implicit none
  !> inputs
  integer,intent(in) :: size_1
  real*8,intent(in) :: tol
  real*8,dimension(size_1),intent(in) :: arr_1,arr_2
  character(len=*),intent(in) :: message
  !> variables
  integer :: ii
  real*8,dimension(size_1) :: error,zeros
  !> compare with relative error
  zeros = 0.d0; error = arr_1-arr_2;
  where(arr_2.ne.0.d0) error = abs(error/arr_2)
  call assert_equals(error,zeros,size_1,tol,message)
end subroutine assert_equals_rel_error_1d_r8

!> assert equals for 2D-real8 arrays with relative error
subroutine assert_equals_rel_error_2d_r8(size_1,size_2,&
arr_1,arr_2,tol,message)
  implicit none
  !> inputs
  integer,intent(in) :: size_1,size_2
  real*8,intent(in) :: tol
  real*8,dimension(size_1,size_2),intent(in) :: arr_1,arr_2
  character(len=*),intent(in) :: message
  !> variables
  integer :: ii
  real*8,dimension(size_1,size_2) :: error,zeros
  !> compare with relative error
  zeros = 0.d0; error = arr_1-arr_2;
  where(arr_2.ne.0.d0) error = abs(error/arr_2)
  call assert_equals(error,zeros,size_1,size_2,tol,message)
end subroutine assert_equals_rel_error_2d_r8

!> assert equals for 3D-real8 arrays with relative error
subroutine assert_equals_rel_error_3d_r8(size_1,size_2,size_3,&
arr_1,arr_2,tol,message)
  implicit none
  !> inputs
  integer,intent(in) :: size_1,size_2,size_3
  real*8,intent(in) :: tol
  real*8,dimension(size_1,size_2,size_3),intent(in) :: arr_1,arr_2
  character(len=*),intent(in) :: message
  !> variables
  integer :: ii
  real*8,dimension(size_1,size_2) :: error,zeros
  !> compare with relative error
  zeros = 0.d0
  do ii=1,size_3
    error = 0.d0; error = arr_1(:,:,ii)-arr_2(:,:,ii);
    where(arr_2(:,:,ii).ne.0.d0) error = abs(error/arr_2(:,:,ii))
    call assert_equals(error,zeros,size_1,size_2,tol,message)
  enddo
end subroutine assert_equals_rel_error_3d_r8

!>----------------------------------------------------------------
end module mod_assert_equals_tools

