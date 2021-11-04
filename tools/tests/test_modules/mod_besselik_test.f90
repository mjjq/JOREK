! This module contains procedures for testing procedured
! computing the modified bessel functions of fractional order
module mod_besselik_test
use fruit 
implicit none

!> variables common to all tests
integer,parameter               :: max_it=1000000 !< maximum number of iterations
real*8,parameter                :: eps=1.d-10     !< variation values
real*8,parameter                :: tol=1.d-16     !< tolerance for success
integer                         :: Nx             !< number of elements
real*8                          :: x_thre         !< threshold value
real*8,dimension(:),allocatable :: x              !< input array x

! procedures for checking and allocating arrays
interface allocate_check
  module procedure allocate_check_integer
	module procedure allocate_check_double
end interface

! procedures for checking and deallocating arrays
interface deallocate_check
  module procedure deallocate_check_integer
	module procedure deallocate_check_double
end interface
contains

! setup initialise the module variables
! inputs:
!   Nx_val:    (integer)(optional) size of the x array
!   Nx_lt_val: (integer)(optional) number of x lower than threshold
!   Nx_ge_val: (integer)(optional) number of x higher than threshold
!   x_val:     (double)(optional) threshold value
!   x_min:     (double)(optional) lower bound of x
!   x_max:     (double)(optional) upper bound of x
subroutine setup(Nx_val,Nx_lt_val,Nx_ge_val,x_val,x_min,x_max)
  implicit none
  
	! inputs
  integer,intent(in),optional :: Nx_val,Nx_lt_val,Nx_ge_val
  real*8,intent(in),optional :: x_val,x_min,x_max
	! variables
	integer :: Nx_lt_ge,Nx_lt_minx,Nx_ge_minx
	integer,dimension(:),allocatable :: ids
	real*8,dimension(:),allocatable  :: rnd
	real*8 :: x_val_loc,x_min_loc,x_max_loc

  ! set default inputs
  Nx = 1000
	Nx_lt_minx = 300
	Nx_ge_minx = 600
	Nx_lt_ge = Nx_lt_minx+Nx_ge_minx
	x_min_loc = -1.d2
  x_max_loc = 1.d2
	x_thre = 2.d0
	! check for inputs
	if(present(Nx_val)) Nx = Nx_val
	if(present(Nx_lt_val)) Nx_lt_minx = Nx_lt_val
  if(present(Nx_ge_val)) Nx_ge_minx = Nx_ge_val
	if(present(x_val)) x_thre = x_val
	if(present(x_min)) x_min_loc = x_min
	if(present(x_max)) x_max_loc = x_max
  Nx_lt_ge = Nx_lt_minx+Nx_ge_minx

	! allocate arrays for tests
	call allocate_check(Nx_lt_ge,ids)
	call allocate_check(Nx,x)
	call allocate_check(Nx,rnd)
  ! allocate internal arrays
	x=0.d0; rnd=0.d0;

	! generate random sequences of xs and ids (we use gnu rng for the tests)
  call generate_int_rnd_array(Nx_lt_ge,1,Nx,ids)
	call random_number(rnd)
  x(ids(1:Nx_lt_minx)) = x_min_loc + (x_thre-eps-x_min_loc)*rnd(1:Nx_lt_minx)
	x(ids(Nx_lt_minx+1:Nx_lt_ge)) = x_thre + (x_max_loc-x_thre)*rnd(Nx_lt_minx+1:Nx_lt_ge)

	! deallocate arrays
	call deallocate_check(ids)
	call deallocate_check(rnd)

end subroutine setup

! teardown of the unit tests
subroutine teardown()
  implicit none
	call deallocate_check_double(x)
end subroutine teardown

! test_split_array tests the split_array_value routine
subroutine test_split_array_value
use mod_besselik, only: split_array_value
use fruit
  implicit none

	! variable for running the test
	integer :: Nx_lt_minx_sol,Nx_ge_minx_sol !< values from routine
	real*8,dimension(Nx) :: x_sol,x_lt_minx_sol,x_ge_minx_sol !< solution from routine
  integer,dimension(Nx) :: ids_lt_minx_sol,ids_ge_minx_sol
  
	! init array
	x_lt_minx_sol=0.d0; x_ge_minx_sol=0.d0
  ! apply solution
	call split_array_value(Nx,x,x_thre,Nx_lt_minx_sol,&
	Nx_ge_minx_sol,x_lt_minx_sol,x_ge_minx_sol,&
	ids_lt_minx_sol,ids_ge_minx_sol)
	
	! reconstruct solution
	x_sol(ids_lt_minx_sol) = x_lt_minx_sol(1:Nx_lt_minx_sol)
	x_sol(ids_ge_minx_sol) = x_ge_minx_sol(1:Nx_ge_minx_sol)

	! check results
	call assert_true(all(x_lt_minx_sol(Nx_lt_minx_sol+1:Nx).le.tol),&
	"Error in x_lt_minx order: non-zero values found in tail")
	call assert_true(all((x_ge_minx_sol(Nx_ge_minx_sol+1:Nx)).le.tol),&
	"Error in x_lt_minx order: non-zero values found in tail")
	call assert_equals(x_sol,x,Nx,tol,&
	"Constructed and reconstructed x array must be the same")

end subroutine test_split_array_value

! allocate integer array if not allocated
subroutine allocate_check_integer(N,array)
  implicit none
	integer,intent(in) :: N
	integer,dimension(:),allocatable,intent(inout) :: array
	if(.not.allocated(array)) allocate(array(N))
end subroutine allocate_check_integer

! allocate double array if not allocated
subroutine allocate_check_double(N,array)
  implicit none
	integer,intent(in) :: N
	real*8,dimension(:),allocatable,intent(inout) :: array
	if(.not.allocated(array)) allocate(array(N))
end subroutine allocate_check_double

! deallocate integer array if allocated
subroutine deallocate_check_integer(array)
  implicit none
	integer,dimension(:),allocatable,intent(inout) :: array
  if(allocated(array)) deallocate(array)
end subroutine deallocate_check_integer

! deallocate array if allocated
subroutine deallocate_check_double(array)
  implicit none
	real*8,dimension(:),allocatable,intent(inout) :: array
  if(allocated(array)) deallocate(array)
end subroutine deallocate_check_double

! generate_int_rnd_array generates an array of random
! integer number without repetitions within an interval.
! inputs:
!   N: (integer) number of elements
!   n_min: (integer) minimum value
!   n_max: (integer) maximum value
! outputs:
!   ids:   (integer)(N) array of non repeated random integers
subroutine generate_int_rnd_array(N,n_min,n_max,ids)
  implicit none
	
	integer,intent(in) :: N,n_min,n_max
	integer,dimension(N),intent(out) :: ids
	integer :: ii,jj,N_repeated
	logical,dimension(N) :: mask
	real*8 :: rnd
	real*8,dimension(N) :: rnd_1d
  
	! generate sequence of random numbers
	call random_number(rnd_1d) 
	ids = floor(n_min+(n_max-n_min+1)*rnd_1d)
	! try to correct for repeated ids
  do jj=1,N
	  ii=1
		do while((count(ids==ids(jj)).gt.1).and.(ii.le.max_it))
		  call random_number(rnd)
		  ids(jj) = floor(n_min+(n_max-n_min+1)*rnd)
			ii=ii+1
		enddo
		call assert_false(ii.gt.max_it,'Failed to generate unique index array')
	enddo
end subroutine generate_int_rnd_array

end module mod_besselik_test
