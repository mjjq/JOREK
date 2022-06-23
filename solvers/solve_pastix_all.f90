#if defined(USE_PASTIX_NEW)
subroutine solve_pastix_all(ptss, n_cpu, my_id, ad_mat, rhs_vec)
#include "pastix_fortran.h"
  use tr_module 
  use mod_parameters, only: n_tor, n_var
  use mpi_mod
  use mod_clock
  use mod_integer_types
  use data_structure, only: type_SP_MATRIX, type_RHS
  use mod_pastix, only: type_PASTIX_SOLVER, pastix_initialize, scale_by_coulmns, pastix_analyze, pastix_set_mat, &
                        pastix_factorize, pastix_solve
  use mod_coicsr

  implicit none

  type(type_PASTIX_SOLVER) :: ptss
  integer                  :: my_id, n_cpu, ierr
  integer(kind=C_INT_ALL)  :: n, nnz
  
  type(type_SP_MATRIX)     :: ad_mat, ac_mat
  type(type_RHS)           :: rhs_vec
  
  integer                            :: index_min, index_max
  integer(kind=int_all)              :: m_loc  
  type(clcktype)                     :: t_itstart, t0, t1, t2, t3
  real*8                             :: tsecond
  integer                            :: block_size, block_size2
  integer(kind=int_all)              :: n_block, nnz_block
  integer(kind=int_all)              :: i
  integer(kind=int_all), allocatable :: sparskit_work(:)

  index_min = ad_mat%index_min
  index_max = ad_mat%index_max
  m_loc = (index_max - index_min + 1) * n_tor * n_var

  call MPI_Allreduce(m_loc,ac_mat%ng,1,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(ad_mat%nnz,ac_mat%nnz,1,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)
  
  n = ac_mat%ng
  nnz = ac_mat%nnz
  
  call clck_time(t0)
  call scale_by_coulmns(ad_mat)
  call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time scale :', tsecond  

  allocate(ac_mat%irn(ac_mat%nnz))
  allocate(ac_mat%jcn(ac_mat%nnz))
  allocate(ac_mat%val(ac_mat%nnz))
  
  call clck_time(t0)
  call split_allgathersolve(n_cpu,my_id,ad_mat,ac_mat)
  call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time mpi_gather :', tsecond
  
#ifdef USE_BLOCK
  block_size  = n_tor * n_var
#else
  block_size  = 1
#endif
  block_size2 = block_size**2  

  n_block   = ac_mat%ng/block_size
  nnz_block = ac_mat%nnz/block_size2
  ac_mat%nblock = n_block
  ac_mat%nzblock = nnz_block  

  if (block_size>1) then
    do i=1,nnz_block  
      ac_mat%irn(i) = (ac_mat%irn((i - 1)*block_size2 + 1) - 1)/block_size + 1 
      ac_mat%jcn(i) = (ac_mat%jcn((i - 1)*block_size2 + 1) - 1)/block_size + 1 
    enddo
  endif

  allocate(sparskit_work(n_block+1))
  
  call clck_time(t0)  
  call coicsr2(n_block,nnz_block,ac_mat%val,ac_mat%irn,ac_mat%jcn,block_size,sparskit_work)
  call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time coicsr2 :', tsecond
  
  deallocate(sparskit_work)
  
  if (.not. ptss%initialized) then
    call pastix_initialize(ptss,MPI_COMM_WORLD)
  endif  
 
  if (.not.ptss%analyzed) then
! -- For PaStiX solver before version 6.x
    allocate(ptss%perm_vars(n_block))
    allocate(ptss%iperm_vars(n_block))
    
    ptss%iparm(IPARM_DOF_NBR) = block_size
    
    call pastix_set_mat(ptss,ac_mat,rhs_vec)
    
    call clck_time(t0)    
    call pastix_analyze(ptss,ac_mat)
    call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time reordering :', tsecond
  endif
  
  call pastix_factorize(ptss,ac_mat)
  call pastix_solve(ptss,ac_mat,rhs_vec)
  
  return
  
  
end subroutine solve_pastix_all
#endif










#ifdef USE_PASTIX
subroutine solve_pastix_all(ptss,n_cpu,my_id, ad_mat, rhs_vec)
!---------------------------------------------------------------------
! subroutine solves the complete system of equation using pastix with
! distributed matrix on the main group mpi_comm_world
!---------------------------------------------------------------------
use tr_module 
use mod_parameters
!use mumps_module
use pastix_module
use global_distributed_matrix
use mpi_mod
use mod_clock
use mod_coicsr
use phys_module, only: use_BLR_compression, epsilon_BLR, just_in_time_BLR, pastix_blr_abs_tol
use mod_integer_types
use data_structure, only: type_SP_MATRIX, type_RHS
use mod_pastix, only: type_PASTIX_SOLVER, scale_by_coulmns, pastix_init_nthreads
 
implicit none

#ifdef USE_PASTIX
#include "pastix_fortran.h"
#else
#include "no_pastix_fortran.h"
#endif

! --- Routine variables
integer                           :: n_cpu, my_id
! --- Local variables
real*8,               allocatable :: column_local(:)
type(clcktype)                    :: t_itstart, t0, t1, t2, t3
real*8                            :: tsecond
integer                           :: i, k, j, ierr
integer(kind=int_all)             :: m_loc
integer(kind=int_all),allocatable :: counts(:), displacements(:)
integer(kind=int_all), parameter  :: Int1=1

type(type_SP_MATRIX)     :: ad_mat, ac_mat
type(type_RHS)           :: rhs_vec
integer                            :: index_min, index_max
type(type_PASTIX_SOLVER) :: ptss
integer(kind=int_all) :: n, nnz
integer(kind=int_all),allocatable, target :: perm_vars(:), iperm_vars(:)

!write(*,*) my_id,'*********************************'
!write(*,*) my_id,'*  solve global matrix (PastiX) *'
!write(*,*) my_id,'*********************************'

index_min = ad_mat%index_min
index_max = ad_mat%index_max
m_loc = (index_max - index_min + 1) * n_tor * n_var
!mumps_par%nz_loc = nz_glob

call MPI_Allreduce(m_loc,n,1,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_Allreduce(ad_mat%nnz,nnz,1,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)

!mumps_par%n = n
!mumps_par%nz = nnz

ac_mat%ng  = n
ac_mat%nnz = nnz

allocate(ac_mat%irn(ac_mat%nnz))
allocate(ac_mat%jcn(ac_mat%nnz))
allocate(ac_mat%val(ac_mat%nnz))

call scale_by_coulmns(ad_mat)

call clck_time(t0)

call split_allgathersolve(n_cpu,my_id,ad_mat,ac_mat)

call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time mpi_gather :', tsecond

call clck_time(t0)

#ifdef USE_BLOCK
block_size  = n_tor * n_var
#else
block_size = 1
#endif
block_size2 = block_size**2

n_block   = ac_mat%ng/block_size
nnz_block = ac_mat%nnz/block_size2

if (block_size > 1) then
  do i=1,nnz_block  
    ac_mat%irn(i) = (ac_mat%irn((i-1)*block_size2+1) - 1) / block_size + 1 
    ac_mat%jcn(i) = (ac_mat%jcn((i-1)*block_size2+1) - 1) / block_size + 1 
  enddo
endif

if (allocated(sparskit_work)) deallocate(sparskit_work)
allocate(sparskit_work(n_block+1))

call coicsr2(n_block,nnz_block,ac_mat%val,ac_mat%irn(1:nnz_block),ac_mat%jcn(1:nnz_block),block_size,sparskit_work)

if (allocated(sparskit_work)) deallocate(sparskit_work)

call clck_time(t1)
call clck_ldiff(t0,t1,tsecond)
if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time coicsr :', tsecond

! End of matrix preparation

allocate(perm_vars(n_block)); perm_vars = 0
allocate(iperm_vars(n_block)); iperm_vars = 0

ptss%perm_vars  => perm_vars
ptss%iperm_vars => iperm_vars

call pastix_init_nthreads(ptss)
if (my_id .eq. 0) write(*,'(i5,A,i5)') my_id,' PastiX n_threads : ', ptss%iparm(IPARM_THREAD_NBR)

if (.not. ptss%initialized) then
  
  ptss%iparm(IPARM_MODIFY_PARAMETER)  = API_NO          ! insert default values
  ptss%iparm(IPARM_START_TASK)        = API_TASK_INIT   ! initializse
  ptss%iparm(IPARM_END_TASK)          = API_TASK_INIT  

  if (my_id .eq. 0) then
    write(*,*) '***********************************'
    write(*,*) '* initialise PastiX               *'
    write(*,*) '***********************************'
  endif
                        
  call pastix_fortran(ptss%idata,MPI_COMM_WORLD,n_block,ac_mat%jcn,ac_mat%irn,ac_mat%val, &
                        ptss%perm_vars,ptss%iperm_vars,rhs_vec%val,Int1,ptss%iparm,ptss%dparm)
 
  ptss%iparm(IPARM_VERBOSE)               = pastix_verb              
  ptss%iparm(IPARM_ITERMAX)               = pastix_iter                ! refinement : max number of iterations

  ptss%iparm(IPARM_FACTORIZATION)         = pastix_facto
  ptss%iparm(IPARM_THREAD_NBR)            = pastix_nthrd               ! number of threads
  ptss%iparm(IPARM_INCOMPLETE)            = pastix_ricar
  ptss%iparm(IPARM_LEVEL_OF_FILL)         = pastix_iluk
  ptss%dparm(DPARM_EPSILON_REFINEMENT)    = pastix_epsilon             ! error level refinement
  ptss%dparm(DPARM_EPSILON_MAGN_CTRL)     = pastix_pivot               ! pivot threshold
  ptss%iparm(IPARM_DOF_NBR)               = block_size                 ! block size
  ptss%iparm(IPARM_RHS_MAKING)            = pastix_rhs                 ! right hand side (0 : use RHS)
  ptss%iparm(IPARM_SYM)                   = pastix_sym
  ptss%iparm(IPARM_AMALGAMATION_LEVEL)    = pastix_amalg
#ifdef FUNNELED
  ptss%iparm(IPARM_THREAD_COMM_MODE)      = API_THREAD_FUNNELED
#else
  ptss%iparm(IPARM_THREAD_COMM_MODE)      = API_THREAD_MULTIPLE
#endif

  ptss%initialized = .true.  
  
endif

if (.not. ptss%analyzed) then
  
  ptss%iparm(IPARM_THREAD_NBR) = pastix_nthrd
  ptss%iparm(IPARM_START_TASK) = API_TASK_ORDERING
  ptss%iparm(IPARM_END_TASK)   = API_TASK_ANALYSE

  call clck_time(t0)
  
  if (my_id .eq. 0) then
    write(*,*) '***********************************'
    write(*,*) '* analyse PastiX                  *'
    write(*,*) '***********************************'
  endif

  call pastix_fortran(ptss%idata,MPI_COMM_WORLD, n_block, ac_mat%jcn, ac_mat%irn, ac_mat%val, &
                      ptss%perm_vars,ptss%iperm_vars,rhs_vec%val,Int1,ptss%iparm,ptss%dparm)
 
  call clck_time(t1)
  call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time analysis :', tsecond

  ptss%analyzed = .true.

endif
call clck_time(t0)

ptss%iparm(IPARM_START_TASK) = API_TASK_NUMFACT
ptss%iparm(IPARM_END_TASK)   = pastix_endsolve

if (my_id .eq. 0) then
  write(*,*) '***********************************'
  write(*,*) '* solve PastiX                     *'
  write(*,*) '***********************************'
endif

call pastix_fortran(ptss%idata,MPI_COMM_WORLD, n_block, ac_mat%jcn, ac_mat%irn, ac_mat%val, &
                    ptss%perm_vars,ptss%iperm_vars,rhs_vec%val,Int1,ptss%iparm,ptss%dparm)                    

call clck_time(t1)
call clck_ldiff(t0,t1,tsecond)
if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed facto/solve :', tsecond

do k=1,ac_mat%ng
  rhs_vec%val(k) =  rhs_vec%val(k)/ad_mat%column_scaling(k)
enddo


return
end
#endif

#ifdef USE_PASTIX6
!> subroutine solves the complete system of equation using PaStiX 6.2
subroutine solve_pastix_all(n_cpu,my_id,index_min,index_max)
  use mod_pastix

  use tr_module 
  use mod_parameters
  use mumps_module
  use global_distributed_matrix, only: column_scaling, deltas, nz_glob, &
                                       irn_glob, jcn_glob, a_glob, rhs_glob  
  use mpi_mod
  use mod_clock
  use mod_integer_types
  use data_structure, only: type_SP_MATRIX  

!$ use omp_lib

  implicit none

! --- Routine parameters
  integer,               intent(in) :: n_cpu, my_id
  integer,               intent(in) :: index_min, index_max

! --- Local variables
  real*8,               allocatable :: column_local(:)
  type(clcktype)                    :: t_itstart, t0, t1, t2, t3
  real*8                            :: tsecond
  integer                           :: i, k, j, ierr
  integer(kind=int_all)             :: m_loc
  integer(kind=int_all),allocatable :: counts(:), displacements(:)
  integer(kind=int_all), parameter  :: Int1=1
  type(type_SP_MATRIX)   :: ad_mat, ac_mat
  
!write(*,*) my_id,'*****************************************'
!write(*,*) my_id,'*  solve global matrix using PaStiX 6.2 *'
!write(*,*) my_id,'*****************************************'

  m_loc = (index_max - index_min + 1) * n_tor * n_var
  mumps_par%nz_loc = nz_glob

  call MPI_Allreduce(m_loc,mumps_par%n,1,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(mumps_par%nz_loc,mumps_par%nz,1,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)
  
  call clck_time(t0)

!------------------------------------------------------ collect the distributed matrix onto all procs
  if (allocated(counts))        call tr_deallocate(counts,"counts",CAT_DMATRIX)
  if (allocated(displacements)) call tr_deallocate(displacements,"displacements",CAT_DMATRIX)

  call tr_allocate(counts,1,n_cpu,"counts",CAT_DMATRIX)
  call tr_allocate(displacements,1,n_cpu,"displacements",CAT_DMATRIX)

  call MPI_Allgather(mumps_par%nz_loc,1,MPI_INTEGER_ALL,counts,1,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)

  displacements(1) = 0
  do i=2,n_cpu
    displacements(i) = displacements(i-1) + counts(i-1)
  enddo

  if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"mumps_par%IRN",CAT_DMATRIX)
  if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"mumps_par%JCN",CAT_DMATRIX)
  if (associated(mumps_par%a) )  call tr_deallocatep(mumps_par%a,"mumps_par%A",CAT_DMATRIX)
  if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"mumps_par%rhs",CAT_DMATRIX)

  call tr_allocatep(mumps_par%irn,Int1,mumps_par%nz,"mumps_par%IRN",CAT_DMATRIX)
  call tr_allocatep(mumps_par%jcn,Int1,mumps_par%nz,"mumps_par%JCN",CAT_DMATRIX)
  call tr_allocatep(mumps_par%a,Int1,mumps_par%nz,"mumps_par%A",CAT_DMATRIX)
  call tr_allocatep(mumps_par%rhs,Int1,mumps_par%n,"mumps_par%rhs",CAT_DMATRIX)

  call split_allgathersolve(n_cpu,my_id,counts,displacements,ad_mat,ac_mat)

  call MPI_AllReduce(rhs_glob,mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  
  call clck_time(t1)
  call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time mpi_gather :', tsecond
  
  call clck_time(t0)
  
  if (.not. spm_initialized) then
    call pastix_init(MPI_COMM_WORLD)
    spm_initialized = .true.
  endif

  call pastix_set_mat(mumps_par%n,mumps_par%nz,mumps_par%irn,mumps_par%jcn,mumps_par%a,1,MPI_COMM_WORLD,&
                         UPDATE=spm_analyzed,DISTRIBUTED=.false.,EQUILIBRIUM=.false.)
  
  if (.not. spm_analyzed) then
    call clck_time(t0)
    call pastix_analyze()    
    spm_analyzed = .true.
    call clck_time(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time analysis :', tsecond    
  endif
  
  call clck_time(t0)

  call pastix_factorize()   
  call pastix_solve(mumps_par%n,mumps_par%rhs,REFINE=.true.)
 
  call clck_time(t1)
  call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed facto/solve :', tsecond

  deltas(1:mumps_par%n) =  mumps_par%rhs(1:mumps_par%n)

  if (allocated(counts))        call tr_deallocate(counts,"counts",CAT_DMATRIX)
  if (allocated(displacements)) call tr_deallocate(displacements,"displacements",CAT_DMATRIX)  
  
  return
end subroutine solve_pastix_all
#endif
