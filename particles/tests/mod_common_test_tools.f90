!> mod_common_test_tools contains and procedure common
!> to multiple unit tests
module mod_common_test_tools
implicit none

private
public :: omp_initialize_rngs

contains
!> Procedures -------------------------------------------------------------
!> initialise the random number generators in omp loops
!> inputs:
!>   n_points_loc: (integer) number of points for rngs init
subroutine omp_initialize_rngs(n_points_loc,n_rngs,rngs)
  use mod_rng, only: type_rng
  use mod_random_seed,only: random_seed
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: n_points_loc,n_rngs  
  !> inputs-outputs
  class(type_rng),dimension(n_rngs),intent(inout) :: rngs !< random number generators
  !> variables
  integer :: n_threads,n_points_per_thread,thread_id,ifail
  !> initialise the rngs using the pcg32
  n_threads = 1
  thread_id = 1
  !$omp parallel default(private) shared(rngs,n_threads,ifail)
  !$ thread_id = omp_get_thread_num()+1
  !$ n_threads = omp_get_num_threads()
  n_points_per_thread = n_points_loc/n_threads
  call rngs(thread_id)%initialize(n_dims=n_points_per_thread,&
  seed=random_seed(),n_streams=n_threads,i_stream=thread_id,ierr=ifail)
  !$omp end parallel 
end subroutine omp_initialize_rngs
!>-------------------------------------------------------------------------
end module mod_common_test_tools
