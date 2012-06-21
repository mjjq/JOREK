module mpi_mod

! --- For forchk
#ifdef FORCHECK
  use mpi
#endif
!
  implicit none

! --- For Lahey.
#ifndef FORCHECK
  include 'mpif.h'
#endif

end module mpi_mod
