!> Definitions of integer types for solvers to switch between short and long ints
module mod_integer_types
  use mpi
  use iso_c_binding
  use iso_fortran_env

! --- Generic integers, valid for all solvers
#ifdef INTSIZE64
  integer, parameter :: int_all = int64
  integer, parameter :: MPI_INTEGER_ALL = MPI_INTEGER8
  integer, parameter :: C_INT_ALL = C_INT64_T
#else
  integer, parameter :: int_all = int32
  integer, parameter :: MPI_INTEGER_ALL = MPI_INTEGER
  integer, parameter :: C_INT_ALL = C_INT
#endif

end module mod_integer_types


