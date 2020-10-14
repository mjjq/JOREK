!> Definitions of integer types for solvers to switch between short and long ints
module mod_integer_types
  use mpi
  use iso_c_binding
  use iso_fortran_env

! --- Generic integers, valid for all solvers
#ifdef INTSIZE64
  integer, parameter :: int_all = int64!selected_int_kind (8) ! apparently, integer*8 and integer(kind=selected_int_kind(8)) is not the same!
  integer, parameter :: MPI_INTEGER_ALL = MPI_INTEGER8
  integer, parameter :: C_INT_ALL = C_LONG
#else
  integer, parameter :: int_all = int32!selected_int_kind (4)
  integer, parameter :: MPI_INTEGER_ALL = MPI_INTEGER
  integer, parameter :: C_INT_ALL = C_INT
#endif

end module mod_integer_types


