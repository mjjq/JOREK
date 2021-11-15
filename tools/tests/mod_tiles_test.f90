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
real*8,parameter               :: tol_real8=1.d-16 !< tolerance for error check
integer,dimension(2),parameter :: interval_int_1d_1=(/-1000,1000/)
integer,dimension(2),parameter :: interval_int_1d_2=(/7500,50000/)
real*8,dimension(2),parameter  :: interval_real8_1d_1=(/-1.d2,1.d2/)
real*8,dimension(2),parameter  :: interval_real8_1d_2=(/1.5d2,3.d3/)
! define interger/double 1d and 2d test arrays
integer,dimension(:),allocatable   :: data_int_1d_1,data_int_1d_2
integer,dimension(:,:),allocatable :: data_int_2d_1,data_int_2d_2
real*8,dimension(:),allocatable    :: data_real8_1d_1,data_real8_1d_2
real*8,dimension(:,:),allocatable  :: data_real8_2d_1,data_real8_2d_2

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
  call allocate_check(N_rows,data_int_1d_1)
  call allocate_check(N_rows,data_int_1d_2)
  call allocate_check(N_rows,N_cols,data_int_2d_1)
  call allocate_check(N_rows,N_cols,data_int_2d_2)
  ! allocate double test arryas
  call allocate_check(N_rows,data_real8_1d_1)
  call allocate_check(N_rows,data_real8_1d_2)
  call allocate_check(N_rows,N_cols,data_real8_2d_1)
  call allocate_check(N_rows,N_cols,data_real8_2d_2)

  ! generate integer random number within range
  call gnu_rng_interval(N_rows,interval_int_1d_1,data_int_1d_1)
  call gnu_rng_interval(N_rows,interval_int_1d_2,data_int_1d_2)
  call gnu_rng_interval(N_rows,N_cols,interval_int_1d_1,data_int_2d_1)
  call gnu_rng_interval(N_rows,N_cols,interval_int_1d_2,data_int_2d_2)
  ! generate double random number within range
  call gnu_rng_interval(N_rows,interval_real8_1d_1,data_real8_1d_1)
  call gnu_rng_interval(N_rows,interval_real8_1d_2,data_real8_1d_2)
  call gnu_rng_interval(N_rows,N_cols,interval_real8_1d_1,data_real8_2d_1)
  call gnu_rng_interval(N_rows,N_cols,interval_real8_1d_2,data_real8_2d_2)

end subroutine setup

! teardown cleans-up the module variables
subroutine teardown()
  use mod_dynamic_array_tools, only: deallocate_check
  implicit none

  ! deallocate integer test arrays
  call deallocate_check(data_int_1d_1)
  call deallocate_check(data_int_1d_2)
  call deallocate_check(data_int_2d_1)
  call deallocate_check(data_int_2d_2)
  ! deallocate double test arrays
  call deallocate_check(data_real8_1d_1)
  call deallocate_check(data_real8_1d_2)
  call deallocate_check(data_real8_2d_1)
  call deallocate_check(data_real8_2d_2)

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
  tol_real8,"Error: allocation and real8 to 0 of tile real8-1D failed!")
  call real8_tile_1d%deallocate_tile()
  call assert_false(allocated(real8_tile_1d%data_array),&
  "Error: deallocation of tile real8-1D failed!")
  !> real8_tile_2d
  call real8_tile_2d%allocate_tile(N_rows,N_cols)
  call assert_equals(real8_tile_2d%data_array,real8_zero_array_2d,N_rows,N_cols,&
  tol_real8,"Error: allocation and real8 to 0 of tile real8-1D failed!")
  call real8_tile_2d%deallocate_tile
  call assert_false(allocated(real8_tile_2d%data_array),&
  "Error: deallocation of tile real8-2D failed!")

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
  call int_tile_1d%allocate_tile(N_rows,data_int_1d_1)
  call assert_equals(int_tile_1d%data_array,data_int_1d_1,N_rows,&
  "Error: allocation and init to data of tile interger-1D failed!")
  call int_tile_1d%deallocate_tile()
  call assert_false(allocated(int_tile_2d%data_array),&
  "Error: deallocation of tile integer-1D failed!")
  !> int_tile_2d
  call int_tile_2d%allocate_tile(N_rows,N_cols,data_int_2d_1)
  call assert_equals(int_tile_2d%data_array,data_int_2d_1,N_rows,&
  N_cols,"Error: allocation and init to data of tile interger-1D failed!")
  call int_tile_2d%deallocate_tile()
  call assert_false(allocated(int_tile_2d%data_array),&
  "Error: deallocation of tile integer-2D failed!")
  !> real8_tile_1d
  call real8_tile_1d%allocate_tile(N_rows,data_real8_1d_1)
  call assert_equals(real8_tile_1d%data_array,data_real8_1d_1,N_rows,&
  tol_real8,"Error: allocation and real8 to data of tile real8-1D failed!")
  call real8_tile_1d%deallocate_tile()
  call assert_false(allocated(real8_tile_1d%data_array),&
  "Error: deallocation of tile real8-1D failed!")
  !> real8_tile_2d
  call real8_tile_2d%allocate_tile(N_rows,N_cols,data_real8_2d_1)
  call assert_equals(real8_tile_2d%data_array,data_real8_2d_1,N_rows,N_cols,&
  tol_real8,"Error: allocation and real8 to data of tile real8-1D failed!")
  call real8_tile_2d%deallocate_tile
  call assert_false(allocated(real8_tile_2d%data_array),&
  "Error: deallocation of tile real8-2D failed!")
end subroutine test_alloc_dealloc_init

!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------

end module mod_tiles_test
