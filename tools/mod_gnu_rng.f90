! the mod_gnu_rng modules contains procedures
! helping the generation of random numbers
module mod_gnu_rng
implicit none

private
public :: gnu_rng_interval
public :: gnu_rng_array_norep

! gnu_rng_interval: template interface for the
! gnu_rng_interval implementation
interface gnu_rng_interval
  module procedure gnu_rng_int_interval_single
  module procedure gnu_rng_int_interval_1d
  module procedure gnu_rng_real8_interval_1d
  module procedure gnu_rng_real8_interval_2d
end interface gnu_rng_interval

! gnu_rng_array_norep: template function for 
! procedures implementing the generation of
! random arrays without repetitions
interface gnu_rng_array_norep
  module procedure gnu_rng_int_array_1d_norep
end interface

contains

! Uniform RNG ------------------------------------------------

! gnu_rng_int_interval_single computes a
! random integer from uniform distribution
! inputs:
!   interval: (2)(integer) minimum and maximum values
! outputs:
!   id: (index) index
subroutine gnu_rng_int_interval_single(interval,id)
  implicit none

  ! inputs
  integer,dimension(2),intent(in) :: interval
  ! outputs
  integer,intent(out) :: id
  ! variables
  real*8 :: rnd

   ! compute random integer 
   call random_number(rnd)
   id = floor(interval(1)+(interval(2)-interval(1)+1)*rnd)

end subroutine gnu_rng_int_interval_single

! gnu_rng_int_interval_1dcomputes an array of
! random integers from a uniform distribution
! inputs:
!   N:        (intrger) length of the array
!   interval: (integer)(2) minimum and maximum values
! outputs:
!   ids: (index)(N) array of random integers
subroutine gnu_rng_int_interval_1d(N,interval,ids)
  implicit none

  ! inputs
  integer,intent(in) :: N
  integer,dimension(2),intent(in) :: interval
  ! outputs
  integer,dimension(N),intent(out) :: ids
  ! variables
  real*8,dimension(N) :: rnds

  ! compute random integer
  call random_number(rnds)
  ids = floor(interval(1)+(interval(2)-interval(1)+1)*rnds)
end subroutine gnu_rng_int_interval_1d

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
subroutine gnu_rng_real8_interval_1d(N,interval,rng_array_1d)
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

end subroutine gnu_rng_real8_interval_1d

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
subroutine gnu_rng_real8_interval_2d(N_rows,N_cols,interval,rng_array_2d)
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

end subroutine gnu_rng_real8_interval_2d

! Additional tools -------------------------------------------

! gnu_rng_int_array_1d_norep generates an array of random
! integer number without repetitions within an interval.
! inputs:
!   N:         (integer) number of elements
!   interval:  (2)(integer) minimum and maximum value
!   n_max:     (integer) maximum value
!   ierr:      (integer) error value
!   max_it_in: (integer)(optional) maximum number of iterations
! outputs:
!   ids:   (integer)(N) array of non repeated random integers
subroutine gnu_rng_int_array_1d_norep(N,interval,ids,ierr,max_it_in)
  implicit none

  ! inputs
  integer,intent(in) :: N
  integer,dimension(2),intent(in) :: interval
  integer,intent(in),optional :: max_it_in
  ! outputs
  integer,dimension(N),intent(out) :: ids
  ! inputs-outputs
  integer,intent(inout) :: ierr
  ! variables
  integer :: ii,jj,max_it
  real*8 :: rnd
  
  ! check for optional inputs
  max_it = 1000000
  if(present(max_it_in)) max_it = max_it_in 
  ! generate sequence of random numbers
  call gnu_rng_int_interval_1d(N,interval,ids)
  ! try to correct for repeated ids
  do jj=1,N
      ii=1
      do while((count(ids==ids(jj)).gt.1).and.(ii.le.max_it))
      call random_number(rnd)
      ids(jj) = floor(interval(1)+(interval(2)-interval(1)+1)*rnd)
      ii=ii+1
    enddo
    if(ii.gt.max_it) ierr=1
  enddo
end subroutine gnu_rng_int_array_1d_norep

end module mod_gnu_rng
