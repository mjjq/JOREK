! The mod_tiles_test module contains variables and procedures used
! for testing the good functionality of the mod_tiles modules
module mod_tiles_test
use fruit
implicit none

private
public :: run_fruit_tiles

! Module variables ------------------------------------------------
integer,parameter              :: N_rows=100 !< number of tile rows
integer,parameter              :: N_cols=100 !< number of tile columns
integer,parameter              :: N_rows_new=75
integer,parameter              :: N_cols_new=80
integer,parameter              :: N_rows_data=50
integer,parameter              :: N_cols_data=30
real*8,parameter               :: tol_real8=1.d-16 !< tolerance for error check
integer,dimension(2),parameter :: interval_int_1d=(/-1000,1000/)
real*8,dimension(2),parameter  :: interval_real8_1d=(/-1.d2,1.d2/)
! define interger/double 1d and 2d test arrays
integer,dimension(:),allocatable   :: data_int_1d
integer,dimension(:,:),allocatable :: data_int_2d
real*8,dimension(:),allocatable    :: data_real8_1d
real*8,dimension(:,:),allocatable  :: data_real8_2d

contains

! Tests basketes --------------------------------------------------

! run_fruit_tiles executes the tiles test set-up, tear-down
! and run the tests
subroutine run_fruit_tiles()
  implicit none
  
  ! execute setup -> tests -> teardown
  write(*,'(/A)') "  ... setting-up: tiles tests"
  call setup 
  write(*,'(/A)') "  ... running: tiles tests"
  call test_alloc_dealloc_noinit !< test tile de-allocation, init=0
  call test_alloc_dealloc_init   !< test tile de-allocation, init=data
  call test_tile_resize          !< test resize tile array
  write(*,'(/A)') "  ... tearing-down: tiles test"
  call teardown

end subroutine run_fruit_tiles

! Set-up and tear-down --------------------------------------------

! setup initiliases the module variables
subroutine setup()
  use mod_dynamic_array_tools, only: allocate_check
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none

  ! allocate integer test arrays
  call allocate_check(N_rows,data_int_1d)
  call allocate_check(N_rows,N_cols,data_int_2d)
  ! allocate double test arryas
  call allocate_check(N_rows,data_real8_1d)
  call allocate_check(N_rows,N_cols,data_real8_2d)

  ! generate integer random number within range
  call gnu_rng_interval(N_rows,interval_int_1d,data_int_1d)
  call gnu_rng_interval(N_rows,N_cols,interval_int_1d,data_int_2d)
  ! generate double random number within range
  call gnu_rng_interval(N_rows,interval_real8_1d,data_real8_1d)
  call gnu_rng_interval(N_rows,N_cols,interval_real8_1d,data_real8_2d)

end subroutine setup

! teardown cleans-up the module variables
subroutine teardown()
  use mod_dynamic_array_tools, only: deallocate_check
  implicit none

  ! deallocate integer test arrays
  call deallocate_check(data_int_1d)
  call deallocate_check(data_int_2d)
  ! deallocate double test arrays
  call deallocate_check(data_real8_1d)
  call deallocate_check(data_real8_2d)

end subroutine teardown

! Tests -----------------------------------------------------------

! test allocation and deallocation of every type of tiles
! initialization is set to zero
subroutine test_alloc_dealloc_noinit()
  use mod_tiles,only: tile_int_1d,tile_int_2d
  use mod_tiles,only: tile_real8_1d,tile_real8_2d
  implicit none

  ! variables
  type(tile_int_1d)   :: int_tile_1d
  type(tile_int_2d)   :: int_tile_2d
  type(tile_real8_1d) :: real8_tile_1d
  type(tile_real8_2d) :: real8_tile_2d
  integer,dimension(N_rows)        :: int_zero_array_1d
  integer,dimension(N_rows,N_cols) :: int_zero_array_2d
  real*8,dimension(N_rows)         :: real8_zero_array_1d 
  real*8,dimension(N_rows,N_cols)  :: real8_zero_array_2d

  ! init all arrays to zero
  int_zero_array_1d = 0
  int_zero_array_2d = 0
  real8_zero_array_1d = 0.d0
  real8_zero_array_2d = 0.d0

  ! check allocation and allocation for each tile type
  !> int_tile_1d
  call int_tile_1d%allocate_tile(N_rows)
  call assert_equals(int_tile_1d%data_array,int_zero_array_1d,N_rows,&
  "Error: allocation and init to 0 of tile interger-1D failed!")
  call int_tile_1d%deallocate_tile()
  call assert_false(allocated(int_tile_2d%data_array),&
  "Error: deallocation of tile integer-1D failed!")
  !> int_tile_2d
  call int_tile_2d%allocate_tile(N_rows,N_cols)
  call assert_equals(int_tile_2d%data_array,int_zero_array_2d,N_rows,&
  N_cols,"Error: allocation and init to 0 of tile interger-1D failed!")
  call int_tile_2d%deallocate_tile()
  call assert_false(allocated(int_tile_2d%data_array),&
  "Error: deallocation of tile integer-2D failed!")
  !> real8_tile_1d
  call real8_tile_1d%allocate_tile(N_rows)
  call assert_equals(real8_tile_1d%data_array,real8_zero_array_1d,N_rows,&
  tol_real8,"Error: allocation and double to 0 of tile double-1D failed!")
  call real8_tile_1d%deallocate_tile()
  call assert_false(allocated(real8_tile_1d%data_array),&
  "Error: deallocation of tile double-1D failed!")
  !> real8_tile_2d
  call real8_tile_2d%allocate_tile(N_rows,N_cols)
  call assert_equals(real8_tile_2d%data_array,real8_zero_array_2d,N_rows,N_cols,&
  tol_real8,"Error: allocation and double to 0 of tile double-1D failed!")
  call real8_tile_2d%deallocate_tile
  call assert_false(allocated(real8_tile_2d%data_array),&
  "Error: deallocation of tile double-2D failed!")

end subroutine test_alloc_dealloc_noinit

!> test allocation and deallocation of tiles with value initialisation
subroutine test_alloc_dealloc_init()
  use mod_tiles,only: tile_int_1d,tile_int_2d
  use mod_tiles,only: tile_real8_1d,tile_real8_2d
  implicit none

  ! variables
  type(tile_int_1d)   :: int_tile_1d
  type(tile_int_2d)   :: int_tile_2d
  type(tile_real8_1d) :: real8_tile_1d
  type(tile_real8_2d) :: real8_tile_2d

  ! check allocation with data initialisation and deallocation
  !> int_tile_1d
  call int_tile_1d%allocate_tile(N_rows,data_int_1d)
  call assert_equals(int_tile_1d%data_array,data_int_1d,N_rows,&
  "Error: allocation and init to data of tile interger-1D failed!")
  call int_tile_1d%deallocate_tile()
  call assert_false(allocated(int_tile_2d%data_array),&
  "Error: deallocation of tile integer-1D failed!")
  !> int_tile_2d
  call int_tile_2d%allocate_tile(N_rows,N_cols,data_int_2d)
  call assert_equals(int_tile_2d%data_array,data_int_2d,N_rows,&
  N_cols,"Error: allocation and init to data of tile interger-1D failed!")
  call int_tile_2d%deallocate_tile()
  call assert_false(allocated(int_tile_2d%data_array),&
  "Error: deallocation of tile integer-2D failed!")
  !> real8_tile_1d
  call real8_tile_1d%allocate_tile(N_rows,data_real8_1d)
  call assert_equals(real8_tile_1d%data_array,data_real8_1d,N_rows,&
  tol_real8,"Error: allocation and double to data of tile double-1D failed!")
  call real8_tile_1d%deallocate_tile()
  call assert_false(allocated(real8_tile_1d%data_array),&
  "Error: deallocation of tile double-1D failed!")
  !> real8_tile_2d
  call real8_tile_2d%allocate_tile(N_rows,N_cols,data_real8_2d)
  call assert_equals(real8_tile_2d%data_array,data_real8_2d,N_rows,N_cols,&
  tol_real8,"Error: allocation and double to data of tile double-1D failed!")
  call real8_tile_2d%deallocate_tile
  call assert_false(allocated(real8_tile_2d%data_array),&
  "Error: deallocation of tile double-2D failed!")
end subroutine test_alloc_dealloc_init

!> test the ability to resize the tile array
subroutine test_tile_resize()
  use mod_tiles,only: tile_int_1d,tile_int_2d
  use mod_tiles,only: tile_real8_1d,tile_real8_2d
  implicit none

  ! variables
  type(tile_int_1d)   :: int_tile_1d
  type(tile_int_2d)   :: int_tile_2d
  type(tile_real8_1d) :: real8_tile_1d
  type(tile_real8_2d) :: real8_tile_2d
  integer,dimension(N_rows_new)            :: resize_data_int_1d
  integer,dimension(N_rows_new,N_cols_new) :: resize_data_int_2d
  real*8,dimension(N_rows_new)             :: resize_data_real8_1d
  real*8,dimension(N_rows_new,N_cols_new)  :: resize_data_real8_2d

  ! initialise resized data array
  !> resized integer 1D
  resize_data_int_1d = 0
  resize_data_int_1d(1:N_rows_data) = data_int_1d(1:N_rows_data)
  !> resized interger 2D
  resize_data_int_2d = 0
  resize_data_int_2d(1:N_rows_data,1:N_cols_data) = &
  data_int_2d(1:N_rows_data,1:N_cols_data)
  !> resized double 1D
  resize_data_real8_1d = 0.d0
  resize_data_real8_1d(1:N_rows_data) = data_real8_1d(1:N_rows_data)
  !> resized double 2D
  resize_data_real8_2d = 0.d0
  resize_data_real8_2d(1:N_rows_data,1:N_cols_data) = &
  data_real8_2d(1:N_rows_data,1:N_cols_data) 

  ! initialise tile array
  call int_tile_1d%allocate_tile(N_rows,data_int_1d)
  call int_tile_2d%allocate_tile(N_rows,N_cols,data_int_2d)
  call real8_tile_1d%allocate_tile(N_rows,data_real8_1d)
  call real8_tile_2d%allocate_tile(N_rows,N_cols,data_real8_2d)

  ! resize tiles
  call int_tile_1d%resize_tile(N_rows_new,N_rows_data)
  call int_tile_2d%resize_tile(N_rows_new,N_cols_new,N_rows_data,N_cols_data)
  call real8_tile_1d%resize_tile(N_rows_new,N_rows_data)
  call real8_tile_2d%resize_tile(N_rows_new,N_cols_new,N_rows_data,N_cols_data)

  ! check procedures
  call assert_equals(int_tile_1d%data_array,resize_data_int_1d,N_rows_new,&
  "Error: resize tile array interger-1D failed!")
  call assert_equals(int_tile_2d%data_array,resize_data_int_2d,N_rows_new,&
  N_cols_new,"Error: resize tile array interger-2D failed!")
  call assert_equals(real8_tile_1d%data_array,resize_data_real8_1d,N_rows_new,&
  tol_real8,"Error: resize tile array double-1D failed!")
  call assert_equals(real8_tile_2d%data_array,resize_data_real8_2d,N_rows_new,&
  N_cols_new,tol_real8,"Error: resize tile array double-2D failed!")

  ! deallocate tiles
  call int_tile_1d%deallocate_tile()
  call int_tile_2d%deallocate_tile()
  call real8_tile_1d%deallocate_tile()
  call real8_tile_2d%deallocate_tile()
  
end subroutine test_tile_resize

!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------

end module mod_tiles_test
