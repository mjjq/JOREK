! This module contains procedures for testing procedured
! computing the modified bessel functions of fractional order
module mod_besselik_test
use fruit 
implicit none

!> variables common to all tests
integer,parameter               :: max_it=1000000 !< maximum number of iterations
real*8,parameter                :: eps=1.d-10     !< variation values
real*8,parameter                :: tol=1.d-16     !< tolerance for success
! order of the bessel function for testing
integer :: N_nu=3
real*8,dimension(3),parameter :: nu=(/1.d0/3.d0,2.d0/3.d0,5.d0/3.d0/)
! values of x used for testing the bessel function for x>=x_val=2
integer                             :: Nx_ge_2=10
real*8,dimension(10),parameter :: x_ge_2=(/7.12d0,8.17d0,2.22d0,19.2d0,&
                                       4.12d0,3.7d0,6.05d0,21.18d0,14.4d0,2.d0/)
! values of x used for testing the bessel function for x<x_val=2
integer                             :: Nx_lt_2=10
real*8,dimension(10),parameter :: x_lt_2=(/1.d-4,1.82d0,0.36d0,5.28d-3,&
                                       2.9d-6,0.27d0,1.74d0,1.d0,1.d-8,1.99d0/)
! values of the bessel function 2nd kind for x_ge_2 and nu=1/3,2/3,5/2
! as computed by matlab (mat) and python (py)
real*8,dimension(30),parameter :: knu_ge_2_mat = (/&
     (/3.764170850828786d-4,1.231043600701636d-4,8.899156804745512d-2,&
     1.307444784232888d-9,9.879203508537923d-3,1.584159892509702d-2,&
     1.188697476014876d-3,1.719285416494302d-10,1.832364180921524d-7,&
     1.165449612961646d-1/),(/3.847728624956395d-4,1.255008926911841d-4,&
     9.476954484002706d-2,1.318559733980535d-9,1.024573759899036d-2,&
     1.649057589808570d-2,1.219489768014882d-3,1.732560154686813d-10,&
     1.852997068247446d-7,1.248389274881278d-1/),(/4.484719282468561d-4,&
     1.435859419047102d-4,1.459102135970209d-1,1.399011432425981d-9,&
     1.319497619429532d-2,2.178414879828106d-2,1.457455551610443d-3,&
     1.828354359187838d-10,2.003937983537028d-7,1.997709129549173d-1/)/)
real*8,dimension(30),parameter :: knu_ge_2_py = (/&
     (/3.7641708508287861d-4,1.2310436007016362d-4,8.8991568047455116d-2,&
       1.3074447842328883d-9,9.8792035085379232d-3,1.5841598925097022d-2,&
       1.1886974760148762d-3,1.7192854164943015d-10,1.8323641809215239d-7,&
       1.1654496129616169d-1,3.8477286249563949d-4,1.2550089269118409d-4,&
       9.4769544840027059d-2,1.3185597339805353d-9,1.0245737598990364d-2,&
       1.6490575898085696d-2,1.2194897680148818d-3,1.7325601546868132d-10,&
       1.8529970682474454d-7,1.2483892748813481d-1/),(/4.4847192824685610d-4,&
       1.4358594190471018d-4,1.4591021359702089d-1,1.3990114324259808d-9,&
       1.3194976194295318d-2,2.1784148798281056d-2,1.4574555516104425d-3,&
       1.8283543591878372d-10,2.0039379835370282d-7,1.9977091295491711d-1/)/)
! values of the bessel function 2nd kind for x_lt_2 and nu=1/3,2/3,5/2
! as computed by matlab (mat) and python (py)
real*8,dimension(30),parameter :: knu_lt_2_mat = (/&
     (/3.628396070100760d1,1.459197599310183d-1,1.313842449151210d0,&
     9.411032704078833d0,1.183202403880544d2,1.626582338799409d0,&
     1.614737084812082d-1,4.384306334415338d-1,7.833229062430851d2,&
     1.179972353239948d-1/),(/4.988585910043111d2,1.571933403240230d-1,&
     1.679684673610099d0,3.540754405283199d1,5.285034257170110d3,&
     2.179872025239006d0,1.744438023694717d-1,4.944750621042081d-1,&
     2.315509105323807d5,1.264314529006558d-1/),(/6.651484164018186d+6,&
     2.610797162123548d-1,7.534896795855283d0,8.950710035944485d3,&
     2.429900926214546d9,1.239138246343648d1,2.951471202585817d-1,&
     1.097730716247145d0,3.087345473843414d13,2.027084265472155d-1/)/)
real*8,dimension(30),parameter :: knu_lt_2_py = (/&
     (/3.6283960701007608d1,1.4591975993101533d-1,1.3138424491512091d0,&
     9.4110327040788331d0,1.1832024038805436d2,1.6265823387994081d0,&
     1.6147370848120510d-1,4.3843063344153255d-1,7.8332290624308507d2,&
     1.1799723532399106d-1/),(/4.9885859100431156d2,1.5719334032402899d-1,&
     1.6796846736101023d0,3.5407544052832023d1,5.2850342571701158d3,&
     2.1798720252390091d0,1.7444380236947751d-1,4.9447506210421149d-1,&
     2.3155091053238098d5,1.2643145290066365d-1/),(/6.6514841640181877d6,&
     2.6107971621235415d-1,7.5348967958552855d0,8.9507100359444867d3,&
     2.4299009262145467d9,1.2391382463436482d1,2.9514712025858114d-1,&
     1.0977307162471459d0,3.0873454738434148d13,2.0270842654721480d-1/)/)
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
