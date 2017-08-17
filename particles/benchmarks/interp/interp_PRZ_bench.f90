!> Benchmark program to test interpolation speeds of mod_interp_PRZ
!> Set n to switch element and position every n steps
program interp_PRZ_bench
use data_structure
use projection_helpers, only: default_polar_grid, f_R, project_f
use mpi
use mod_pcg32_rng
use mod_random_seed
use mod_interp_PRZ
use phys_module
implicit none

integer, parameter :: n_switch = 100 !< Go to a new element every 100 steps (mimic jorek particle pattern)
! we go to a new position inside the element every step.
integer, parameter :: N_loops = 100000
integer, parameter :: N_vector = 16
integer, parameter :: n_v = 1, i_var(n_v) = [1]

integer :: ierr, provided, i, j, k, i_elm
type(type_node_list) :: node_list
type(type_element_list) :: element_list
real*8 :: t0, t1
type(pcg32_rng) :: rng
real*8, dimension(n_vector) :: R, R_s, R_t, Z, Z_s, Z_t
real*8, dimension(n_vector,n_v) :: P, P_s, P_t, P_phi
real*8, dimension(n_vector) :: u, s, t, phi

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

call rng%initialize(n_vector, random_seed(), 1, 1)

call cpu_time(t0)
! Start the test
i_elm = 1
do i=1,N_loops
  if (mod(i,n_switch) == 0 .and. n_switch > 0) then
    call rng%next(u)
    i_elm = nint(1.d0+u(1)*real(element_list%n_elements-1))
  end if
  call rng%next(u); s = u
  call rng%next(u); t = u
  call rng%next(u); phi = u
  !call interp_PRZ_vec(node_list, element_list, i_elm, i_var, n_v, n_vector, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
  do j=1,n_vector
    call interp_PRZ(node_list, element_list, i_elm, i_var, n_v, s(j), t(j), phi(j), &
        P(j,:), P_s(j,:), P_t(j,:), P_phi(j,:), R(j), R_s(j), R_t(j), Z(j), Z_s(j), Z_t(j))
  end do
  call rng%next(u)
  do j=1,n_vector
    if (u(j) .gt. 0.999999) write(*,"(10g16.7)") P(j,1), P_s(j,1), P_t(j,1), P_phi(j,1), R(j), R_s(j), R_t(j), Z(j), Z_s(j), Z_t(j)
  end do
end do
call cpu_time(t1)

write(*,*) "Time for ", N_loops*N_vector, " iterations: ", t1-t0, "s"
write(*,*) "Time per iteration: ", ((t1-t0)/real(N_loops*N_vector)) *1d6, " microseconds"

end program interp_PRZ_bench
