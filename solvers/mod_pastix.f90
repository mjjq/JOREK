module pastix_module             ! PastiX specific variables
  integer(kind=8)       :: pastix_data
  integer               :: pastix_iparm(64)
  real*8                :: pastix_dparm(64)
  integer,allocatable   :: pastix_perm_vars(:), pastix_iperm_vars(:)
  integer,allocatable   :: sparskit_work(:)
  integer,allocatable   :: ihwb(:),iwk(:)
  integer               :: n_block, nnz_block, block_size, block_size2

  logical               :: use_pastix, pastix_initialised, pastix_analysed, pastix_smp_only, no_zeros_pastix

  integer, parameter    :: pastix_facto    = 2
  integer, parameter    :: pastix_rhs      = 0
  integer, parameter    :: pastix_sym      = 1
  integer               :: pastix_nthrd    = 1
  integer, parameter    :: pastix_iter     = 250
  integer, parameter    :: pastix_verb     = 1
  integer, parameter    :: pastix_ricar    = 0
  integer, parameter    :: pastix_iluk     = 3
  integer, parameter    :: pastix_amalg    = 5 
  integer, parameter    :: pastix_endsolve = 5 ! 5: direct solve, 6: refinement step (not required usually)
  real*8,  parameter    :: pastix_epsilon  = 1.d-12
  real*8,  parameter    :: pastix_pivot    = 1.d-64
end module pastix_module
