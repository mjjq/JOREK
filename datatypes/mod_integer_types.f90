!> Definitions of integer types for solvers to switch between short and long ints
module mod_integer_types

! --- Generic integers, valid for all solvers
#ifdef INTSIZE64
  integer, parameter :: int_all = selected_int_kind (8)
#else
  integer, parameter :: int_all = selected_int_kind (4)
#endif

! --- Solver-specific integers, that may change between solvers
#ifdef INTSIZE64
#ifdef USE_MUMPS
  integer, parameter :: int_spec1 = selected_int_kind (4)
#else
  integer, parameter :: int_spec1 = selected_int_kind (8)
#endif
#else
  integer, parameter :: int_spec1 = selected_int_kind (4)
#endif

end module mod_integer_types


