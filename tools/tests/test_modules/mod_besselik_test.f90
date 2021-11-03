! This module contains procedures for testing procedured
! computing the modified bessel functions of fractional order
module mod_besselik_test
use fruit 
implicit none

!> variables common to all tests
integer,parameter :: max_it=1000000 !< maximum number of iterations
real*8,parameter :: eps=1.d-10
real*8,parameter :: tol=1.d-16 !< tolerance for success
integer :: Nx,Nx_lt_minx,Nx_ge_minx
real*8 :: x_thre
real*8,dimension(:),allocatable :: x,x_lt_minx,x_ge_minx

contains

!> This function initialise the module variables
subroutine setup(Nx_val,Nx_lt_val,Nx_ge_val,x_val,x_min,x_max)
  implicit none

  integer,intent(in),optional :: Nx_val,Nx_lt_val,Nx_ge_val
  real*8,intent(in),optional :: x_val,x_min,x_max
	integer :: Nx_lt_ge
	integer,dimension(:),allocatable :: ids
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
	call allocate_check_double(Nx,x)
	call allocate_check_double(Nx,x_lt_minx)
  call allocate_check_double(Nx,x_ge_minx)
  ! allocate internal arrays
	allocate(ids(Nx_lt_ge))
	x=0.d0; x_lt_minx=0.d0; x_ge_minx=0.d0;
	! generate random sequences of xs and ids (we use gnu rng for the tests)
  call generate_int_rnd_array(Nx_lt_ge,1,Nx,ids)
	call random_number(x_lt_minx(1:Nx_lt_minx))
  x(ids(1:Nx_lt_minx)) = x_min_loc + (x_thre-eps-x_min_loc)*x_lt_minx(1:Nx_lt_minx)
	call random_number(x_ge_minx(1:Nx_ge_minx))
	x(ids(Nx_lt_minx+1:Nx_lt_ge)) = x_thre + (x_max_loc-x_thre)*x_ge_minx(1:Nx_ge_minx)
	! reorder arrays for comparison
	x_lt_minx = 0.d0; x_ge_minx=0.d0
	x_lt_minx(ids(1:Nx_lt_minx)) = x(ids(1:Nx_lt_minx))
	x_ge_minx(ids(Nx_lt_minx+1:Nx_ge_minx)) = x(ids(Nx_lt_minx+1:Nx_lt_ge))
	! deallocate internal arrays
  deallocate(ids)

end subroutine setup

! Teardown of the unit tests
subroutine teardown()
  implicit none
	call deallocate_check_double(x)
  call deallocate_check_double(x_lt_minx)
	call deallocate_check_double(x_ge_minx)
end subroutine teardown

! This function testst the split_array_value routine
subroutine test_split_array_value
use mod_besselik, only: split_array_value
use fruit
  implicit none

	! variable for running the test
	integer :: ii,N_zeros,Nx_lt_minx_sol,Nx_ge_minx_sol !< values from routine
	real*8,dimension(Nx) :: x_lt_minx_sol,x_ge_minx_sol !< solution from routine

  ! init array
	x_lt_minx_sol=0.d0; x_ge_minx_sol=0.d0
  ! apply solution
	call split_array_value(Nx,x,x_thre,Nx_lt_minx_sol,&
	Nx_ge_minx_sol,x_lt_minx_sol,x_ge_minx_sol)
  N_zeros = count(x.eq.0.d0)
	if(0.d0.lt.x_thre) then
	  Nx_lt_minx = Nx_lt_minx+N_zeros
	else
	  Nx_ge_minx = Nx_ge_minx+N_zeros
	endif
	! check results
	call assert_equals(Nx_lt_minx_sol,Nx_lt_minx,&
	"Number of x-values smaller than x_threshold must be the same")
	call assert_equals(Nx_ge_minx_sol,Nx_ge_minx,&
	"Number of x-values greater than x_threshold must be the same")
	call assert_equals(x_lt_minx_sol,x_lt_minx,Nx,tol,&
	"x-values smaller than x_threshold must be the same")
	call assert_equals(x_ge_minx_sol,x_ge_minx,Nx,tol,&
	"x-values grater than x_threshold must be the same")

end subroutine test_split_array_value

! allocate array if not allocated
subroutine allocate_check_double(N,array)
  implicit none
	integer,intent(in) :: N
	real*8,dimension(:),allocatable,intent(inout) :: array
	if(.not.allocated(array)) allocate(array(N))
end subroutine allocate_check_double

! deallocate array if allocated
subroutine deallocate_check_double(array)
  implicit none
	real*8,dimension(:),allocatable,intent(inout) :: array
  if(allocated(array)) deallocate(array)
end subroutine deallocate_check_double

! This function generates an array of random integer number
! without repetitions within an interval.
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
