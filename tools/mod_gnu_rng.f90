! the mod_gnu_rng modules contains procedures
! helping the generation of random numbers
module mod_gnu_rng
implicit none

private
public :: gnu_rng_interval

! gnu_rng_interval: template interface for the
! gnu_rng_interval implementation
interface gnu_rng_interval
  module procedure gnu_rng_interval_1d
  module procedure gnu_rng_interval_2d
end interface gnu_rng_interval

contains

! gnu_rng_interval_1d generates a 1D random number array
! using the gnu-fortran intrinsic function within a 
! predefined interval (uniform distribution)
! inputs:
!   N:            (int) length of the array
!   interval:     (2)(real8) interval within the value
!                            are generated
!   rng_array_1d: (N)(real8) 1D array of uniform random
!                            numbers within an interval
! outputs:
!   rng_array_1d: (N)(real8) 1D array of uniform random
!                            numbers within an interval
subroutine gnu_rng_interval_1d(N,interval,rng_array_1d)
  implicit none

  ! inputs:
  integer,intent(in)                :: N
  real*8,dimension(2),intent(in)    :: interval
  ! inputs-outputs:
  real*8,dimension(N),intent(inout) :: rng_array_1d

  ! generate random number array
  call random_number(rng_array_1d)
  rng_array_1d = interval(1) + (interval(2)-&
  interval(1))*rng_array_1d

end subroutine gnu_rng_interval_1d

! gnu_rng_interval_2d generates a 2D random number array
! using the gnu-fortran intrinsic function within a 
! predefined interval (uniform distribution)
! inputs:
!   N_rows:       (int) number of array rows
!   N_cols:       (int) number of array cols
!   interval:     (2)(real8) interval within the value
!                            are generated
!   rng_array_2d: (N_rows,N_cols)(real8) 2D array of 
!                            uniform random numbers 
!                            within an interval
! outputs:
!   rng_array_2d: (N_rows,N_cols)(real8) 2D array of 
!                            uniform random numbers 
!                            within an interval
subroutine gnu_rng_interval_2d(N_rows,N_cols,interval,rng_array_2d)
  implicit none

  ! inputs:
  integer,intent(in)                            :: N_rows,N_cols
  real*8,dimension(2),intent(in)                :: interval
  ! inputs-outputs:
  real*8,dimension(N_rows,N_cols),intent(inout) :: rng_array_2d

  ! generate random number array
  call random_number(rng_array_2d)
  rng_array_2d = interval(1) + (interval(2)-&
  interval(1))*rng_array_2d

end subroutine gnu_rng_interval_2d

end module mod_gnu_rng
