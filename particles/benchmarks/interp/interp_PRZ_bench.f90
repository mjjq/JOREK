!> Benchmark program to test interpolation speeds of mod_interp_PRZ
!> Set n to switch element and position every n steps
program interp_PRZ_bench
use data_structure
use projection_helpers, only: default_polar_grid, f_R, project_f
use mpi
use mod_pcg32_rng
use mod_random_seed
use mod_interp_PRZ
implicit none

integer, parameter :: n_switch = 100 !< Go to a new element every 100 steps (mimic jorek particle pattern)
! we go to a new position inside the element every step.
integer, parameter :: N_tries = 1000000

integer :: ierr, provided, i, k, i_elm
type(type_node_list) :: node_list
type(type_element_list) :: element_list
real*8 :: t0, t1
type(pcg32_rng) :: rng
real*8 :: s, t, phi
real*8 :: R, R_s, R_t, Z, Z_s, Z_t
real*8, dimension(1) :: P, P_s, P_t, P_phi, u

! Setup the grid
call MPI_INIT_THREAD(MPI_THREAD_SINGLE, provided, ierr)
call initialise_basis

call default_polar_grid(node_list, element_list, 80)
call project_f(node_list, element_list, f_R)
! Copy these (from the first harmonic) to every harmonic
do i=1,node_list%n_nodes
  do k=1,n_order+1
    node_list%node(i)%values(1,k,2:) = node_list%node(i)%values(1,k,1)
    node_list%node(i)%deltas(1,k,2:) = node_list%node(i)%deltas(1,k,1)
  end do
enddo

! Wait after this omp part for a while
call sleep(5)

call rng%initialize(1, random_seed(), 1, 1)

call cpu_time(t0)
! Start the test
i_elm = 1
do i=1,N_tries
  if (mod(i,n_switch) == 0 .and. n_switch > 0) then
    call rng%next(u)
    i_elm = nint(1.d0+u(1)*real(element_list%n_elements-1))
  end if
  call rng%next(u); s = u(1)
  call rng%next(u); t = u(1)
  call rng%next(u); phi = u(1)
  call interp_PRZ(node_list, element_list, i_elm, [1], 1, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
  call rng%next(u)
  if (u(1) .gt. 0.999999) write(*,"(10g16.7)") P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t
end do
call cpu_time(t1)

write(*,*) "Time for ", N_tries, " iterations: ", t1-t0, "s"
write(*,*) "Time per iteration: ", ((t1-t0)/real(N_tries)) *1d6, " ms"

end program interp_PRZ_bench
