!> Benchmark the projection function by projection 0, 1, R, RZ and R^4
!> onto a square, polar and flux aligned grid and write the error to a file
program project_f_convergence
use data_structure
use basis_at_gaussian, only: initialise_basis
use projection_helpers, only: project_f, elements_mean_rms, f_0, f_1, f_R, f_RZ, F_R4, &
    default_square_grid, default_polar_grid, default_flux_grid
use mpi
implicit none

abstract interface
  function f(R,Z)
    real*8, intent(in) :: R, Z
    real*8 :: f
  end function f
end interface
procedure (f), pointer :: f_ptr

integer :: i, j, u
integer :: ierr, provided
type(type_node_list) :: node_list
type(type_element_list) :: element_list
real*8 :: m, e

call MPI_INIT_THREAD(MPI_THREAD_SINGLE, provided, ierr)
call initialise_basis

open(newunit=u, file='project_f.csv', status='replace')
write(u,"(A)") "Grid; f(R,Z)=0; f(R,Z)=1; f(R,Z)=R; f(R,Z)=RZ; f(R,Z)=R^4"

do i=1,5
  select case(i)
  case(1)
    node_list%n_nodes = 0
    call default_square_grid(node_list, element_list, 10)
    write(u,"(A)",advance="no") "square_10_10"
  case(2)
    call default_polar_grid(node_list, element_list, 31)
    write(u,"(A)",advance="no") "polar_30_31"
  case(3)
    call default_polar_grid(node_list, element_list, 32)
    write(u,"(A)",advance="no") "polar_30_32"
  case(4)
    call default_flux_grid(node_list, element_list, 31)
    write(u,"(A)",advance="no") "flux_30_31"
  case(5)
    call default_flux_grid(node_list, element_list, 32)
    write(u,"(A)",advance="no") "flux_30_32"
  end select

  do j=1,5
    select case(j)
    case(1)
      f_ptr => f_0
    case(2)
      f_ptr => f_1
    case(3)
      f_ptr => f_R
    case(4)
      f_ptr => f_RZ
    case(5)
      f_ptr => f_R4
    end select

    call project_f(node_list, element_list, f_ptr)
    call elements_mean_rms(node_list, element_list, f_ptr, m, e)
    write(u,"(A, g14.8)", advance="no") "; ", e
  end do
  write(u,*) ""
end do
close(u)
call MPI_Finalize(ierr)
end program project_f_convergence
