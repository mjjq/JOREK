module mpi_mod

#ifndef LAHEY
#ifdef MPI_F08
  use mpi_f08
#else
  use mpi
#endif
#endif

  implicit none

! --- For Lahey.
#ifdef LAHEY
  include 'mpif.h'
#endif

end module mpi_mod
