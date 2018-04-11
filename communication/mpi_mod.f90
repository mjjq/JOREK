module mpi_mod

! backwards compatibility
#ifdef LAHEY
#define INCLUDE_MPIFH 1
#endif

#ifndef INCLUDE_MPIFH
#ifdef MPI_F08
  use mpi_f08
#else
  use mpi
#endif
#endif

  implicit none

  ! --- If 
#ifdef INCLUDE_MPIFH
  include 'mpif.h'
#endif

end module mpi_mod
