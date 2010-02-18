module pastix_module             ! PastiX specific variables
  integer(kind=8)       :: pastix_data
  integer               :: pastix_iparm(64)
  real*8                :: pastix_dparm(64)
  integer,allocatable   :: pastix_perm_vars(:), pastix_iperm_vars(:)
  integer,allocatable   :: sparskit_work(:)
  integer,allocatable   :: ihwb(:),iwk(:)

  logical               :: use_pastix, pastix_initialised, pastix_analysed, pastix_smp_only

  integer, parameter    :: pastix_facto    = 2
  integer, parameter    :: pastix_rhs      = 0
  integer, parameter    :: pastix_sym      = 1
  integer               :: pastix_nthrd    = 1
  integer, parameter    :: pastix_iter     = 255
  integer, parameter    :: pastix_ricar    = 0
  integer, parameter    :: pastix_iluk     = 3
  integer, parameter    :: pastix_amalg    = 5
  integer, parameter    :: pastix_endsolve = 6
  real*8,  parameter    :: pastix_epsilon  = 1.d-15
  real*8,  parameter    :: pastix_pivot    = 1.d-64
endmodule pastix_module
