! The mod_tiles_test module contains variables and procedures used
! for testing the good functionality of the mod_tiles modules
module mod_tiles_test
use fruit
implicit none

private
public :: run_fruit_tiles

! Module variables ------------------------------------------------
integer,parameter    :: N_rows=100 !< number of tile rows
integer,parameter    :: N_cols=100 !< number of tile columns
integer,dimension(2) :: interval_int_1d_1=(/-1000,1000/)
integer,dimension(2) :: interval_int_1d_2=(/7500,50000/)
real*8,dimension(2)  :: interval_real8_1d_1=(/-1.d2,1.d2/)
real*8,dimension(2)  :: interval_real8_1d_2=(/1.5d2,3.d3/)
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
!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------

end module mod_tiles_test
