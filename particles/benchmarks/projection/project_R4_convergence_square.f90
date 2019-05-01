!> Benchmark the projection routines by projecting R^4 onto a square grid
!> and calculating the error
program project_R4_convergence_square
use data_structure
use basis_at_gaussian, only: initialise_basis
use projection_helpers, only: project_f, elements_mean_rms, F_R4, &
    default_square_grid
use mpi
implicit none

integer :: i, j, u, n
integer :: ierr, provided
type(type_node_list) :: node_list
type(type_element_list) :: element_list
real*8 :: m, e
integer, parameter :: n_exp = 8

call MPI_INIT_THREAD(MPI_THREAD_SINGLE, provided, ierr)
call initialise_basis

open(newunit=u, file='project_R4_square.csv', status='replace')
write(u,"(A)") "n, epsilon"

do i=0,nint(1.8d0*real(n_exp))
  n = nint(3d0*10.d0**(real(i)/real(n_exp)))
  call default_square_grid(node_list, element_list, n)
  write(u,"(i3)",advance="no") n-1
  call project_f(node_list, element_list, f_R4)
  call elements_mean_rms(node_list, element_list, f_R4, m, e)
  write(u,"(A, g14.8)") ", ", e
end do
close(u)
call MPI_Finalize(ierr)
end program project_R4_convergence_square
