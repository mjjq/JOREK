!> the module mod_lights_assert_equals_tools contains procedures
!> for comparing various data structures
module mod_assert_equals_tools
use fruit
implicit none

private
public :: assert_equals_rel_error
public :: assert_equals_allocatable_arrays

!> Variables -----------------------------------------------------
!> Interfaces ----------------------------------------------------
interface assert_equals_rel_error
  module procedure assert_equals_rel_error_r8
  module procedure assert_equals_rel_error_1d_r8
  module procedure assert_equals_rel_error_2d_r8
  module procedure assert_equals_rel_error_3d_r8
end interface assert_equals_rel_error

interface assert_equals_allocatable_arrays
  module procedure assert_equals_allocatable_array_value_1d_int
  module procedure assert_equals_allocatable_array_value_1d_r8
  module procedure assert_equals_allocatable_array_value_2d_r8
  module procedure assert_equals_allocatable_array_value_3d_r8
  module procedure assert_equals_allocatable_array_value_4d_r8
  module procedure assert_equals_allocatable_arrays_1d_r8
  module procedure assert_equals_allocatable_arrays_2d_r8
end interface assert_equals_allocatable_arrays

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

!> assert equals allocatable arrays 1d real 8
subroutine assert_equals_allocatable_arrays_1d_r8(n_values,&
array_test,array_sol,tol,message)
  implicit none
  !> inputs:
  integer,intent(in)                         :: n_values
  real*8,intent(in)                          :: tol
  real*8,dimension(:),allocatable,intent(in) :: array_test
  real*8,dimension(n_values),intent(in)      :: array_sol
  character(len=*),intent(in)                :: message
  !> tests
  call assert_true(allocated(array_test),trim(message//" not allocated!"))
  call assert_equals(size(array_test),n_values,trim(message//" size mismatch!"))
  if(allocated(array_test)) &
  call assert_equals(array_test,array_sol,n_values,tol,&
  trim(message//" mismatch!"))
end subroutine assert_equals_allocatable_arrays_1d_r8

!> assert equals allocatble array 2d real 8
subroutine assert_equals_allocatable_arrays_2d_r8(n_values_1,&
n_values_2,array_test,array_sol,tol,message) 
  implicit none
  !> inputs:
  integer,intent(in)                                 :: n_values_1,n_values_2
  real*8,intent(in)                                  :: tol
  real*8,dimension(:,:),allocatable,intent(in)       :: array_test
  real*8,dimension(n_values_1,n_values_2),intent(in) :: array_sol
  character(len=*),intent(in)                        :: message
  !> tests
  call assert_true(allocated(array_test),trim(message//" not allocated!"))
  call assert_equals(shape(array_test),(/n_values_1,n_values_2/),&
  2,trim(message//" size mismatch!"))
  if(allocated(array_test)) &
  call assert_equals(array_test,array_sol,n_values_1,n_values_2,tol,&
  trim(message//" mismatch!"))
end subroutine assert_equals_allocatable_arrays_2d_r8

!> assert allocatable array equal to value 1d integer
subroutine assert_equals_allocatable_array_value_1d_int(n_values,&
array_test,value_sol,message)
  implicit none
  !> inputs:
  integer,intent(in) :: n_values
  integer,intent(in) :: value_sol
  integer,dimension(:),allocatable,intent(in) :: array_test
  character(len=*),intent(in) :: message
  !> test
  call assert_true(allocated(array_test),trim(message//" not allocated!"))
  call assert_equals(size(array_test),n_values,trim(message//" size mismatch!")) 
  if(allocated(array_test)) &
  call assert_true(all(array_test.eq.value_sol),trim(message//" mismatch!"))
end subroutine assert_equals_allocatable_array_value_1d_int

!> assert allocatable array equal to value 1d real8
subroutine assert_equals_allocatable_array_value_1d_r8(n_values,&
array_test,value_sol,message)
  implicit none
  !> inputs:
  integer,intent(in) :: n_values
  real*8,intent(in) :: value_sol
  real*8,dimension(:),allocatable,intent(in) :: array_test
  character(len=*),intent(in) :: message
  !> test
  call assert_true(allocated(array_test),trim(message//" not allocated!"))
  call assert_equals(size(array_test),n_values,trim(message//" size mismatch!")) 
  if(allocated(array_test)) &
  call assert_true(all(array_test.eq.value_sol),trim(message//" mismatch!"))
end subroutine assert_equals_allocatable_array_value_1d_r8

!> assert allocatable array equal to value 2d real8
subroutine assert_equals_allocatable_array_value_2d_r8(n_values_1,&
n_values_2,array_test,value_sol,message)
  implicit none
  !> inputs:
  integer,intent(in) :: n_values_1,n_values_2
  real*8,intent(in) :: value_sol
  real*8,dimension(:,:),allocatable,intent(in) :: array_test
  character(len=*),intent(in) :: message
  !> test
  call assert_true(allocated(array_test),trim(message//" not allocated!"))
  call assert_equals(shape(array_test),(/n_values_1,n_values_2/),&
  2,trim(message//" size mismatch!")) 
  if(allocated(array_test)) &
  call assert_true(all(array_test.eq.value_sol),trim(message//" mismatch!"))
end subroutine assert_equals_allocatable_array_value_2d_r8

!> assert allocatable array equal to value 3d real8
subroutine assert_equals_allocatable_array_value_3d_r8(n_values_1,&
n_values_2,n_values_3,array_test,value_sol,message)
  implicit none
  !> inputs:
  integer,intent(in) :: n_values_1,n_values_2,n_values_3
  real*8,intent(in) :: value_sol
  real*8,dimension(:,:,:),allocatable,intent(in) :: array_test
  character(len=*),intent(in) :: message
  !> test
  call assert_true(allocated(array_test),trim(message//" not allocated!"))
  call assert_equals(shape(array_test),(/n_values_1,n_values_2,n_values_3/),&
  3,trim(message//" size mismatch!")) 
  if(allocated(array_test)) &
  call assert_true(all(array_test.eq.value_sol),trim(message//" mismatch!"))
end subroutine assert_equals_allocatable_array_value_3d_r8

!> assert allocatable array equal to value 4d real8
subroutine assert_equals_allocatable_array_value_4d_r8(n_values_1,&
n_values_2,n_values_3,n_values_4,array_test,value_sol,message)
  implicit none
  !> inputs:
  integer,intent(in) :: n_values_1,n_values_2,n_values_3,n_values_4
  real*8,intent(in) :: value_sol
  real*8,dimension(:,:,:,:),allocatable,intent(in) :: array_test
  character(len=*),intent(in) :: message
  !> test
  call assert_true(allocated(array_test),trim(message//" not allocated!"))
  call assert_equals(shape(array_test),(/n_values_1,n_values_2,n_values_3,&
  n_values_4/),4,trim(message//" size mismatch!")) 
  if(allocated(array_test)) &
  call assert_true(all(array_test.eq.value_sol),trim(message//" mismatch!"))
end subroutine assert_equals_allocatable_array_value_4d_r8

!>----------------------------------------------------------------
end module mod_assert_equals_tools

