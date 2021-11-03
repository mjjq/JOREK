! This module contains procedures for testing procedured
! computing the modified bessel functions of fractional order
module mod_besselik_test
use fruit 
implicit none

!> variables common to all tests
real*8,parameter :: eps=1.d-10
integer :: Nx,Nx_lt_minx,Nx_ge_minx
real*8 :: x_thre
real*8,dimension(:),allocatable :: x,x_lt_minx,x_ge_minx

contains

!> This function initialise the module variables
subroutine setup(Nx_val,Nx_lt_val,Nx_ge_val,x_val,x_min,x_max)
  implicit none

  integer,intent(in),optional :: Nx_val,Nx_lt_val,Nx_ge_val
  real*8,intent(in),optional :: x_val,x_min,x_max
	integer,dimension(:),allocatable :: ids
	real*8 :: x_val_loc,x_min_loc,x_max_loc


  ! set default inputs
  Nx = 1000
	Nx_lt_minx = 300
	Nx_ge_minx = 600
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

	! allocate arrays for tests
	call allocate_check_double(Nx,x)
	call allocate_check_double(Nx,x_lt_minx)
  call allocate_check_double(Nx,x_ge_minx)
  ! allocate internal arrays
	allocate(ids(Nx))
	! generate random sequences of xs and ids (we use gnu rng for the tests)
	call random_number(x_lt_minx(1:Nx_lt_minx))
  x_lt_minx(1:Nx_lt_minx) = x_min_loc + (x_thre-eps-x_min_loc)*x_lt_minx(1:Nx_lt_minx)
	call random_number(x_ge_minx(1:Nx_ge_minx))
	x_lt_minx(1:Nx_ge_minx) = x_thre + (x_max_loc-x_thre)*x_ge_minx(1:Nx_ge_minx)
	call generate_int_rnd_array(Nx,1,Nx,ids)
	! assemble x array
	x(ids(1:Nx_lt_minx)) = x_lt_minx(1:Nx_lt_minx)
	x(ids(Nx_lt_minx+1:Nx_ge_minx)) = x_ge_minx(1:Nx_ge_minx)
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
	real*8,dimension(N) :: rnd
  
	! generate sequence of random numbers
	call random_number(rnd) 
	ids = floor(n_min+(n_max-n_min+1)*rnd)
	! try to correct for repeated ids
  do jj=1,N
	  ii=n_min
		do while((count(ids==ids(jj)).gt.1).and.(ii.ge.n_max))
		  ids(jj) = ii
			ii = ii+1
		enddo
		call assert_false(ii.gt.n_max,'Failed to generate unique index array')
	enddo
	
end subroutine generate_int_rnd_array

end module mod_besselik_test
