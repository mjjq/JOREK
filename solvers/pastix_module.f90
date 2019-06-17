module pastix_module             ! PastiX specific variables
#ifndef USE_PASTIX6
#ifdef USE_PASTIX 
#include "pastix_fortran.h"
#else
#include "no_pastix_fortran.h"
#endif
#endif

#ifdef USE_PASTIX6
  use iso_c_binding
  use pastixf
  use pastix_enums
  use spmf
#endif

  integer,allocatable   :: sparskit_work(:)
  integer,allocatable   :: ihwb(:),iwk(:)
  integer               :: n_block, nnz_block, block_size, block_size2

  logical               :: use_pastix, pastix_initialised, pastix_analysed, pastix_smp_only, no_zeros_pastix

#ifndef USE_PASTIX6
  integer(kind=8)       :: pastix_data
  integer               :: pastix_iparm(IPARM_SIZE)
  real*8                :: pastix_dparm(DPARM_SIZE)
  integer,allocatable   :: pastix_perm_vars(:), pastix_iperm_vars(:)
  integer, parameter    :: pastix_facto    = API_FACT_LU
  integer, parameter    :: pastix_sym      = API_SYM_NO
  integer, parameter    :: pastix_verb     = API_VERBOSE_NO
  integer, parameter    :: pastix_endsolve = API_TASK_SOLVE
  integer, parameter    :: pastix_rhs      = 0
#else
  type(pastix_data_t), pointer     :: pastix_data
  type(spmatrix_t),    allocatable :: pastix_spm
  type(spmatrix_t),    allocatable :: pastix_spm_check
  integer(kind=pastix_int_t), target   :: pastix_iparm(iparm_size)
  real(kind=c_double),        target   :: pastix_dparm(dparm_size)
  integer, parameter    :: pastix_facto    = PastixFactLU
  integer, parameter    :: pastix_sym      = PastixGeneral
  integer, parameter    :: pastix_verb     = PastixVerboseNo 
#endif

  integer               :: pastix_nthrd    = 1
  integer, parameter    :: pastix_iter     = 250
  integer, parameter    :: pastix_ricar    = 0
  integer, parameter    :: pastix_iluk     = 3
  integer, parameter    :: pastix_amalg    = 5 
  real*8,  parameter    :: pastix_epsilon  = 1.d-12
  real*8                :: pastix_pivot    = 1.d-64
  !> Sometimes PaStiX is faster if we limit the number of threads it is allowed to use.
  !! These limits apply to the total number of threads/node = threads/mpi_task * mpi_tasks/node.
  !! By default, the number of threads equals to OMP_NUM_THREADS, and the high limits here
  !! do not change this. 
  integer               :: pastix_maxthrd = 1024

  contains

  subroutine pastix_init_num_threads(my_id)
    use mpi_mod
!$  use omp_lib
    implicit none
    integer, intent(in) :: my_id
    !$omp parallel default(none) shared(pastix_nthrd)
    !$omp master
!$      pastix_nthrd = omp_get_num_threads()
    !$omp end master
    !$omp end parallel
    if (pastix_nthrd * get_tasks_per_node() > pastix_maxthrd) then
      pastix_nthrd = max(pastix_maxthrd / get_tasks_per_node(), 1)
    endif
    if (my_id .eq. 0) then
      write(*,'(i5,A,i5)') my_id,' PastiX n_threads : ', pastix_nthrd
    end if
  end subroutine
end module pastix_module
