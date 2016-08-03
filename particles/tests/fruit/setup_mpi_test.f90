module setup_mpi_test
contains
subroutine setup
  use mpi
  integer :: ierr
  call MPI_init(ierr)
end subroutine setup

subroutine teardown
  use mpi
  integer :: ierr
  call MPI_finalize(ierr)
end subroutine teardown
end module setup_mpi_test
