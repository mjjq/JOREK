subroutine solve_pastix_all(n_cpu,my_id,index_min,index_max)
!---------------------------------------------------------------------
! subroutine solves the complete system of equation using pastix with
! distributed matrix on the main group mpi_comm_world
!---------------------------------------------------------------------
use tr_module 
use mod_parameters
use mumps_module
use pastix_module
use global_distributed_matrix
use mpi_mod
use mod_clock
use mod_coicsr
use phys_module, only: use_BLR_compression, epsilon_BLR, just_in_time_BLR

!$ use omp_lib

#ifdef USE_PASTIX6
! -- For PaStiX solver version 6.x
use iso_c_binding
use pastixf
use pastix_enums
use spmf
#endif
 
implicit none

#ifndef USE_PASTIX6
#ifdef USE_PASTIX
#include "pastix_fortran.h"
#else
#include "no_pastix_fortran.h"
#endif
#endif


integer                  :: n_cpu, index_min, index_max       ! global index_min, index_max for this cpu
real*8,allocatable       :: column_local(:)
type(clcktype)           :: t_itstart, t0, t1, t2, t3
real*8                   :: tsecond
integer                  :: i, k, j, ierr, my_id, m_loc
integer,allocatable      :: counts(:), displacements(:)
#ifdef USE_PASTIX6
! -- For PaStiX solver version 6.x
integer(c_int)     :: pastix_info
type(c_ptr)        :: pastix_rhs_ptr
integer(kind=spm_int_t), dimension(:), pointer       :: pastix_colptr
integer(kind=spm_int_t), dimension(:), pointer       :: pastix_rowptr
real(kind=c_double)    , dimension(:), pointer       :: pastix_values
#endif


!write(*,*) my_id,'*********************************'
!write(*,*) my_id,'*  solve global matrix (PastiX) *'
!write(*,*) my_id,'*********************************'

m_loc = (index_max - index_min + 1) * n_tor * n_var
mumps_par%nz_loc = nz_glob

call MPI_Allreduce(m_loc,mumps_par%N,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_Allreduce(mumps_par%NZ_loc,mumps_par%nz,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

!------------------------------------------------------- colunm scaling of global distributed matrix
call clck_time(t0)

if (allocated(column_scaling))  call tr_deallocate(column_scaling,"column_scaling",CAT_DMATRIX)
if (allocated(column_local))    call tr_deallocate(column_local,"column_local",CAT_DMATRIX)
call tr_allocate(column_scaling,1,mumps_par%N,"column_scaling",CAT_DMATRIX)
call tr_allocate(column_local,1,mumps_par%N,"column_local",CAT_DMATRIX)

column_local = 1.d-20;   column_scaling = 1.d-20
do k=1,nz_glob
  j = jcn_glob(k)
  column_local(j) = max(column_local(j),abs(A_glob(k)))
enddo

call MPI_AllReduce(column_local,column_scaling,mumps_par%N,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

do k = 1, nz_glob
  j = jcn_glob(k)
  A_glob(k) = A_glob(k) / column_scaling(j)
enddo

call clck_time(t1)
call clck_ldiff(t0,t1,tsecond)
if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time scale :', tsecond

call clck_time(t0)

!------------------------------------------------------ collect the distributed matrix onto all procs
if (allocated(counts))        call tr_deallocate(counts,"counts",CAT_DMATRIX)
if (allocated(displacements)) call tr_deallocate(displacements,"displacements",CAT_DMATRIX)

call tr_allocate(counts,1,n_cpu,"counts",CAT_DMATRIX)
call tr_allocate(displacements,1,n_cpu,"displacements",CAT_DMATRIX)

call MPI_Allgather(mumps_par%nz_loc,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

displacements(1) = 0
do i=2,n_cpu
  displacements(i) = displacements(i-1) + counts(i-1)
enddo

if (associated(mumps_par%IRN)) call tr_deallocatep(mumps_par%IRN,"mumps_par%IRN",CAT_DMATRIX)
if (associated(mumps_par%JCN)) call tr_deallocatep(mumps_par%JCN,"mumps_par%JCN",CAT_DMATRIX)
if (associated(mumps_par%A) )  call tr_deallocatep(mumps_par%A,"mumps_par%A",CAT_DMATRIX)
if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"mumps_par%rhs",CAT_DMATRIX)

call tr_allocatep(mumps_par%IRN,1,mumps_par%nz,"mumps_par%IRN",CAT_DMATRIX)
call tr_allocatep(mumps_par%JCN,1,mumps_par%nz,"mumps_par%JCN",CAT_DMATRIX)
call tr_allocatep(mumps_par%A,1,mumps_par%nz,"mumps_par%A",CAT_DMATRIX)
call tr_allocatep(mumps_par%rhs,1,mumps_par%n,"mumps_par%rhs",CAT_DMATRIX)

call MPI_AllgatherV(IRN_glob,mumps_par%nz_loc,MPI_INTEGER,mumps_par%IRN, &
                    counts,displacements,MPI_INTEGER,MPI_COMM_WORLD,ierr)

call MPI_AllgatherV(JCN_glob,mumps_par%nz_loc,MPI_INTEGER,mumps_par%JCN, &
                    counts,displacements,MPI_INTEGER,MPI_COMM_WORLD,ierr)

call MPI_AllgatherV(A_glob,mumps_par%nz_loc,MPI_DOUBLE_PRECISION,mumps_par%A, &
                    counts,displacements,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

call MPI_AllReduce(RHS_glob,mumps_par%RHS,mumps_par%N,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

call clck_time(t1)
call clck_ldiff(t0,t1,tsecond)
if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time mpi_gather :', tsecond
call clck_time(t0)

#ifdef USE_BLOCK
!---------------------------- reduce IRN,JCN to make use of blocksize ntor*nvar
!                             temporary solution before using blocks everywhere

block_size  = n_tor * n_var
block_size2 = block_size**2
!---------------------------- reduce IRN,JCN to make use of blocksize ntor*nvar
n_block   = mumps_par%n  / block_size
nnz_block = mumps_par%nz / block_size2

do i=1,nnz_block  
  mumps_par%irn(i) = (mumps_par%irn((i-1)*block_size2+1) - 1) / block_size + 1 
  mumps_par%jcn(i) = (mumps_par%jcn((i-1)*block_size2+1) - 1) / block_size + 1 
enddo

if (allocated(sparskit_work)) deallocate(sparskit_work)
allocate(sparskit_work(n_block+1))

call coicsr2(n_block,nnz_block,mumps_par%A,mumps_par%IRN(1:nnz_block),mumps_par%JCN(1:nnz_block),block_size,sparskit_work)

#ifndef USE_PASTIX6
! -- For PaStiX solver before version 6.x
if (.not. allocated(pastix_perm_vars))  call tr_allocate(pastix_perm_vars,1,n_block,"pastix_perm_vars",CAT_UNKNOWN)
if (.not. allocated(pastix_iperm_vars)) call tr_allocate(pastix_iperm_vars,1,n_block,"pastix_iperm_vars",CAT_UNKNOWN)
#endif

#else

if (allocated(sparskit_work)) deallocate(sparskit_work)
allocate(sparskit_work(mumps_par%N + 1))

call coicsr(mumps_par%N,mumps_par%NZ,1,mumps_par%A,mumps_par%IRN,mumps_par%JCN,sparskit_work)
#endif

if (allocated(sparskit_work)) deallocate(sparskit_work)
call clck_time(t1)
call clck_ldiff(t0,t1,tsecond)
if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time coicsr :', tsecond


#ifndef USE_PASTIX6
! -- For PaStiX solver before version 6.x
if (.not. allocated(pastix_perm_vars))  call tr_allocate(pastix_perm_vars,1,mumps_par%n,"pastix_perm_vars",CAT_UNKNOWN)
if (.not. allocated(pastix_iperm_vars)) call tr_allocate(pastix_iperm_vars,1,mumps_par%n,"pastix_iperm_vars",CAT_UNKNOWN)
#else
! -- For PaStiX solver version 6.x
allocate(pastix_spm) ! Replace by tr_allocate etc.?!
call spmInit(pastix_spm)

#ifdef USE_BLOCK
pastix_spm%n           =  n_block
pastix_spm%nnz         =  nnz_block
pastix_spm%dof         =  block_size
#else
pastix_spm%n           =  mumps_par%n
pastix_spm%nnz         =  mumps_par%nz
pastix_spm%dof         =  1
#endif
call spmUpdateComputedFields(pastix_spm)
call spmAlloc(pastix_spm)

call c_f_pointer(pastix_spm%colptr,pastix_colptr, [pastix_spm%n+1])
call c_f_pointer(pastix_spm%rowptr,pastix_rowptr, [pastix_spm%nnz])
call c_f_pointer(pastix_spm%values,pastix_values, [mumps_par%nz])

pastix_colptr      = mumps_par%jcn(1:pastix_spm%n+1)
pastix_rowptr      = mumps_par%irn(1:pastix_spm%nnz)
pastix_values      = mumps_par%A(1:mumps_par%nz)
    
!pastix_colptr      = mumps_par%jcn
!pastix_rowptr      = mumps_par%irn
!pastix_values      = mumps_par%A
#endif

call pastix_init_num_threads(my_id)

if (.not. pastix_initialised) then

#ifndef USE_PASTIX6
  ! -- For PaStiX solver before version 6.x
  pastix_iparm(IPARM_MODIFY_PARAMETER)  = API_NO          ! insert default values
  pastix_iparm(IPARM_START_TASK)        = API_TASK_INIT   ! initializse
  pastix_iparm(IPARM_END_TASK)          = API_TASK_INIT
#else
  ! -- For PaStiX solver version 6.x
  call pastixInitParam(pastix_iparm, pastix_dparm)
#endif
  if (my_id .eq. 0) then
    write(*,*) '***********************************'
    write(*,*) '* initialise PastiX               *'
    write(*,*) '***********************************'
  endif
  
#ifndef USE_PASTIX6
  ! -- For PaStiX solver before version 6.x
#ifdef USE_BLOCK
  call pastix_fortran(pastix_data,MPI_COMM_WORLD,n_block,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                        pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#else
  call pastix_fortran(pastix_data,MPI_COMM_WORLD,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                        pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#endif
#endif

  ! pastix input parameters working in Pastix5 and Pastix6
  pastix_iparm(IPARM_VERBOSE)               = pastix_verb              
  pastix_iparm(IPARM_ITERMAX)               = pastix_iter                ! refinement : max number of iterations

  pastix_iparm(IPARM_FACTORIZATION)         = pastix_facto
  pastix_iparm(IPARM_THREAD_NBR)            = pastix_nthrd               ! number of threads
  pastix_iparm(IPARM_INCOMPLETE)            = pastix_ricar
  pastix_iparm(IPARM_LEVEL_OF_FILL)         = pastix_iluk
  pastix_dparm(DPARM_EPSILON_REFINEMENT)    = pastix_epsilon             ! error level refinement
  pastix_dparm(DPARM_EPSILON_MAGN_CTRL)     = pastix_pivot               ! pivot threshold
#ifdef USE_BLOCK
  pastix_iparm(IPARM_DOF_NBR)               = block_size                 ! block size
#else
  pastix_iparm(IPARM_DOF_NBR)               = 1
#endif



#ifndef USE_PASTIX6
  ! -- For PaStiX solver before version 6.x
  pastix_iparm(IPARM_RHS_MAKING)            = pastix_rhs                 ! right hand side (0 : use RHS)
  pastix_iparm(IPARM_SYM)                   = pastix_sym
  pastix_iparm(IPARM_AMALGAMATION_LEVEL)    = pastix_amalg
#ifdef WORLDWAR2
#ifdef FUNNELED
  pastix_iparm(IPARM_THREAD_COMM_MODE)      = API_THREAD_FUNNELED
#else
  pastix_iparm(IPARM_THREAD_COMM_MODE)      = API_THREAD_MULTIPLE
#endif
#endif


#else
  ! -- For PaStiX solver version 6.x
  pastix_iparm(IPARM_MTX_TYPE)              = pastix_sym
  pastix_iparm(IPARM_AMALGAMATION_LVLCBLK)  = pastix_amalg

! TEMPORARY: not yet relevant for Pastix6 as MPI parallelisation is not implemented
!#ifdef WORLDWAR2
!#ifdef FUNNELED
!   pastix_iparm(IPARM_THREAD_COMM_MODE)      = PastixThreadFunneled
!#else
!   pastix_iparm(IPARM_THREAD_COMM_MODE)      = PastixThreadMultiple
!#endif
!#endif

  ! BLR Compression
  if (use_BLR_compression) then
    if (just_in_time_BLR) then
      pastix_iparm(IPARM_COMPRESS_WHEN)     = PastixCompressWhenEnd ! Just-in-Time (speed optimal)
    else 
      pastix_iparm(IPARM_COMPRESS_WHEN)     = PastixCompressWhenBegin ! Minimal-memory (default)
    endif
    pastix_dparm(DPARM_COMPRESS_TOLERANCE)  = epsilon_BLR

!!   Additional PaStiX compression parameters (currently set to their default values)
!    pastix_iparm(IPARM_COMPRESS_ORTHO)      = PastixCompressOrthoCGS
!    pastix_iparm(IPARM_COMPRESS_METHOD)     = PastixCompressMethodPQRCP
!    pastix_iparm(IPARM_COMPRESS_MIN_WIDTH)  = 120
!    pastix_iparm(IPARM_COMPRESS_MIN_HEIGHT) = 20
!    pastix_dparm(DPARM_COMPRESS_MIN_RATIO)  = 1.0
  endif

  call pastixInit(pastix_data, 0, pastix_iparm, pastix_dparm)    ! TEMPORARY: 0 should be pastix_comm but pastix6 is not yet MPI parallelised!
#endif
  pastix_initialised = .true.
endif

if (.not. pastix_analysed) then

#ifndef USE_PASTIX6
  ! -- For PaStiX solver before version 6.x
  pastix_iparm(IPARM_THREAD_NBR) = pastix_nthrd
  pastix_iparm(IPARM_START_TASK) = API_TASK_ORDERING
  pastix_iparm(IPARM_END_TASK)   = API_TASK_ANALYSE
#endif 
  call clck_time(t0)
  
  if (my_id .eq. 0) then
    write(*,*) '***********************************'
    write(*,*) '* analyse PastiX                  *'
    write(*,*) '***********************************'
  endif

#ifndef USE_PASTIX6
  ! -- For PaStiX solver before version 6.x
#ifdef USE_BLOCK
  
  call pastix_fortran(pastix_data,MPI_COMM_WORLD, n_block, &
                      mumps_par%jcn(1:n_block+1), mumps_par%irn(1:nnz_block), mumps_par%A(1:mumps_par%nz), &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

#else

  call pastix_fortran(pastix_data,MPI_COMM_WORLD, mumps_par%n, &
                      mumps_par%jcn(1:mumps_par%n+1), mumps_par%irn(1:mumps_par%nz), mumps_par%A(1:mumps_par%nz), &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

#endif

#else
  ! -- For PaStiX solver version 6.x
#ifdef USE_BLOCK
  call pastix_subtask_order(pastix_data,pastix_spm,pastix_myorder,pastix_info)
  call pastix_subtask_symbfact(pastix_data,pastix_info)
  call pastix_subtask_reordering(pastix_data,pastix_info)

  ! Expand spm matrix and pastix analysis substructures because rest of Pastix6 cannot handle multiple dofs (yet)
  call pastixExpand(pastix_data,pastix_spm)

  call pastix_subtask_blend(pastix_data,pastix_info)
#else
  call pastix_task_analyze(pastix_data,pastix_spm,pastix_info)
#endif

#endif
 
  call clck_time(t1)
  call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time analysis :', tsecond

  pastix_analysed = .true.
#if (defined(USE_PASTIX6) && defined(USE_BLOCK))
  pastix_analysed = .false. ! Necessary for now such that the spm expansion is done every time step. 
                            ! Can be removed once the PaStiX team has implemented multi-dof for all pastix_subtasks.
#endif
endif ! .not. pastix_analysed

call clck_time(t0)

#ifndef USE_PASTIX6
! -- For PaStiX solver before version 6.x
pastix_iparm(IPARM_THREAD_NBR) = pastix_nthrd
pastix_iparm(IPARM_START_TASK) = API_TASK_NUMFACT
pastix_iparm(IPARM_END_TASK)   = pastix_endsolve
#ifdef WORLDWAR2
#ifdef FUNNELED
  pastix_iparm(IPARM_THREAD_COMM_MODE)  = API_THREAD_FUNNELED
#else
  pastix_iparm(IPARM_THREAD_COMM_MODE)  = API_THREAD_MULTIPLE
#endif
#endif
#endif

if (my_id .eq. 0) then
  write(*,*) '***********************************'
  write(*,*) '* call PastiX                     *'
  write(*,*) '***********************************'
endif

#ifndef USE_PASTIX6
! -- For PaStiX solver before version 6.x
#ifdef USE_BLOCK

pastix_iparm(IPARM_DOF_NBR)            = block_size

call pastix_fortran(pastix_data,MPI_COMM_WORLD, n_block,                                                 &
                    mumps_par%jcn(1:n_block+1), mumps_par%irn(1:nnz_block), mumps_par%A(1:mumps_par%nz), &
                    pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

#else

call pastix_fortran(pastix_data,MPI_COMM_WORLD, mumps_par%n, mumps_par%jcn, mumps_par%irn, mumps_par%A, &
                    pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

#endif

#else
! -- For PaStiX solver version 6.x

call pastix_task_numfact(pastix_data,pastix_spm,pastix_info)
!call pastix_subtask_spm2bcsc(pastix_data,pastix_spm,pastix_info )
!call pastix_subtask_bcsc2ctab(pastix_data,pastix_info )
!call pastix_subtask_sopalin(pastix_data,pastix_info )

pastix_rhs_ptr = c_loc(mumps_par%rhs)
call pastix_task_solve(pastix_data,1,pastix_rhs_ptr,pastix_spm%n,pastix_info)
call spmExit(pastix_spm)
deallocate(pastix_spm)
#endif


call clck_time(t1)
call clck_ldiff(t0,t1,tsecond)
if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed facto/solve :', tsecond

do k=1,mumps_par%n
  deltas(k) =  mumps_par%rhs(k)  / column_scaling(k)
enddo

if (allocated(column_local))  call tr_deallocate(column_local,"column_local",CAT_DMATRIX)
if (allocated(counts))        call tr_deallocate(counts,"counts",CAT_DMATRIX)
if (allocated(displacements)) call tr_deallocate(displacements,"displacements",CAT_DMATRIX)

return
end
