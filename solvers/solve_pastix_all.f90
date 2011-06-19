subroutine solve_pastix_all(n_cpu,my_id,index_min,index_max)
!---------------------------------------------------------------------
! subroutine solves the complete system of equation using pastix with
! distributed matrix on the main group mpi_comm_world
!---------------------------------------------------------------------
use tr_module 
use parameters
use mumps_module
use pastix_module
use global_distributed_matrix
implicit none
include 'mpif.h'
#include "pastix_fortran.h"

integer                  :: n_cpu, index_min, index_max       ! global index_min, index_max for this cpu
real*8,allocatable       :: column_local(:)
real*8                   :: t_analysis_0, t_analysis_1, t_fact_0, t_fact_1, t_comm_0, t_comm_1
real*8                   :: t_scale_0, t_scale_1
integer                  :: i, k, j, ierr, my_id, m_loc
integer,allocatable      :: counts(:), displacements(:)

integer, external :: omp_get_num_threads, omp_get_thread_num
  
!write(*,*) my_id,'*********************************'
!write(*,*) my_id,'*  solve global matrix (PastiX) *'
!write(*,*) my_id,'*********************************'

m_loc = (index_max - index_min + 1) * n_tor * n_var
mumps_par%nz_loc = nz_glob

call MPI_Allreduce(m_loc,mumps_par%N,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_Allreduce(mumps_par%NZ_loc,mumps_par%nz,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

!------------------------------------------------------- colunm scaling of global distributed matrix
call cpu_time(t_scale_0)

if (allocated(column_scaling))  call tr_deallocate(column_scaling,"column_scaling")
if (allocated(column_local))    call tr_deallocate(column_local,"column_local")
call tr_allocate(column_scaling,1,mumps_par%N,"column_scaling")
call tr_allocate(column_local,1,mumps_par%N,"column_local")

column_local = 1.d-20;   column_scaling = 1.d-20
do k=1,nz_glob
  j = jcn_glob(k)
  column_local(j) = max(column_local(j),abs(A_glob(k)))
enddo

call MPI_AllReduce(column_local,column_scaling,mumps_par%N,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
do k=1,nz_glob
  j = jcn_glob(k)
  A_glob(k) = A_glob(k) / column_scaling(j)
enddo

call cpu_time(t_scale_1)

if (my_id .eq. 0)  write(*,'(A,f8.3)') ' PASTIX, scale     : ',t_scale_1-t_scale_0


call cpu_time(t_comm_0)

!------------------------------------------------------ collect the distributed matrix onto all procs
if (allocated(counts))        call tr_deallocate(counts,"counts")
if (allocated(displacements)) call tr_deallocate(displacements,"displacements")

call tr_allocate(counts,1,n_cpu,"counts")
call tr_allocate(displacements,1,n_cpu,"displacements")

call MPI_Allgather(mumps_par%nz_loc,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

displacements(1) = 0
do i=2,n_cpu
  displacements(i) = displacements(i-1) + counts(i-1)
enddo

if (associated(mumps_par%IRN)) call tr_deallocatep(mumps_par%IRN,"mumps_par%IRN")
if (associated(mumps_par%JCN)) call tr_deallocatep(mumps_par%JCN,"mumps_par%JCN")
if (associated(mumps_par%A) )  call tr_deallocatep(mumps_par%A,"mumps_par%A")
if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"mumps_par%rhs")

call tr_allocatep(mumps_par%IRN,1,mumps_par%nz,"mumps_par%IRN")
call tr_allocatep(mumps_par%JCN,1,mumps_par%nz,"mumps_par%JCN")
call tr_allocatep(mumps_par%A,1,mumps_par%nz,"mumps_par%A")
call tr_allocatep(mumps_par%rhs,1,mumps_par%n,"mumps_par%rhs")

call MPI_AllgatherV(IRN_glob,mumps_par%nz_loc,MPI_INTEGER,mumps_par%IRN, &
                    counts,displacements,MPI_INTEGER,MPI_COMM_WORLD,ierr)

call MPI_AllgatherV(JCN_glob,mumps_par%nz_loc,MPI_INTEGER,mumps_par%JCN, &
                    counts,displacements,MPI_INTEGER,MPI_COMM_WORLD,ierr)

call MPI_AllgatherV(A_glob,mumps_par%nz_loc,MPI_DOUBLE_PRECISION,mumps_par%A, &
                    counts,displacements,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

call MPI_AllReduce(RHS_glob,mumps_par%RHS,mumps_par%N,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)


#IFDEF USE_BLOCK
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

if (.not. allocated(pastix_perm_vars))  call tr_allocate(pastix_perm_vars,1,n_block,"pastix_perm_vars")
if (.not. allocated(pastix_iperm_vars)) call tr_allocate(pastix_iperm_vars,1,n_block,"pastix_iperm_vars")

#ELSE

if (allocated(sparskit_work)) deallocate(sparskit_work)
allocate(sparskit_work(mumps_par%N + 1))

call coicsr(mumps_par%N,mumps_par%NZ,1,mumps_par%A,mumps_par%IRN,mumps_par%JCN,sparskit_work)

if (.not. allocated(pastix_perm_vars))  call tr_allocate(pastix_perm_vars,1,mumps_par%n,"pastix_perm_vars")
if (.not. allocated(pastix_iperm_vars)) call tr_allocate(pastix_iperm_vars,1,mumps_par%n,"pastix_iperm_vars")

#ENDIF

call cpu_time(t_comm_1)

if (my_id .eq. 0)  write(*,'(A,f8.3)') ' PASTIX, comm      : ',t_comm_1-t_comm_0


!$omp parallel default(none) shared(pastix_nthrd)    
pastix_nthrd = omp_get_num_threads()
!$omp end parallel

if (my_id .eq. 0) write(*,'(i5,A,i5)') my_id,' PastiX n_threads : ',pastix_nthrd 

if (.not. pastix_initialised) then

  pastix_iparm(IPARM_MODIFY_PARAMETER)  = API_NO          ! insert default values
  pastix_iparm(IPARM_START_TASK)        = API_TASK_INIT   ! initializse
  pastix_iparm(IPARM_END_TASK)          = API_TASK_INIT

  if (my_id .eq. 0) then
    write(*,*) '***********************************'
    write(*,*) '* initialise PastiX               *'
    write(*,*) '***********************************'
  endif
  
#IFDEF USE_BLOCK
  call pastix_fortran(pastix_data,MPI_COMM_WORLD,n_block,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#ELSE
  call pastix_fortran(pastix_data,MPI_COMM_WORLD,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#ENDIF

  pastix_iparm(IPARM_VERBOSE)            = pastix_verb
  pastix_iparm(IPARM_ITERMAX)            = pastix_iter                 ! refinement : max number of iterations

  pastix_iparm(IPARM_FACTORIZATION)      = pastix_facto
  pastix_iparm(IPARM_THREAD_NBR)         = pastix_nthrd               ! number of threads
  pastix_iparm(IPARM_RHS_MAKING)         = pastix_rhs                 ! right hand side (0 : use RHS)
  
  pastix_iparm(IPARM_SYM)                = pastix_sym
  
  pastix_iparm(IPARM_INCOMPLETE)         = pastix_ricar
  pastix_iparm(IPARM_LEVEL_OF_FILL)      = pastix_iluk
  pastix_iparm(IPARM_AMALGAMATION_LEVEL) = pastix_amalg
  
  pastix_dparm(DPARM_EPSILON_REFINEMENT) = pastix_epsilon             ! error level refinement
  pastix_dparm(DPARM_EPSILON_MAGN_CTRL)  = pastix_pivot               ! pivot threshold?
  pastix_initialised = .true.

#IFDEF USE_BLOCK
  pastix_iparm(IPARM_DOF_NBR)            = block_size
#ELSE
  pastix_iparm(IPARM_DOF_NBR)            = 1
#ENDIF

endif

if (.not. pastix_analysed) then

  pastix_iparm(IPARM_START_TASK) = API_TASK_ORDERING
  pastix_iparm(IPARM_END_TASK)   = API_TASK_ANALYSE
 
  call cpu_time(t_analysis_0)
  
  if (my_id .eq. 0) then
    write(*,*) '***********************************'
    write(*,*) '* analyse PastiX                  *'
    write(*,*) '***********************************'
  endif

#IFDEF USE_BLOCK
  
  call pastix_fortran(pastix_data,MPI_COMM_WORLD, n_block, &
                      mumps_par%jcn(1:n_block+1), mumps_par%irn(1:nnz_block), mumps_par%A(1:mumps_par%nz), &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

#ELSE

  call pastix_fortran(pastix_data,MPI_COMM_WORLD, mumps_par%n, &
                      mumps_par%jcn(1:mumps_par%n+1), mumps_par%irn(1:mumps_par%nz), mumps_par%A(1:mumps_par%nz), &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

#ENDIF

  call cpu_time(t_analysis_1)

  pastix_analysed = .true.

  if (my_id .eq. 0)  write(*,'(A,f8.3)') ' PASTIX, analysis  : ',t_analysis_1-t_analysis_0

endif

call cpu_time(t_fact_0)

pastix_iparm(IPARM_START_TASK) = API_TASK_NUMFACT
pastix_iparm(IPARM_END_TASK)   = pastix_endsolve

if (my_id .eq. 0) then
  write(*,*) '***********************************'
  write(*,*) '* call PastiX                     *'
  write(*,*) '***********************************'
endif

#IFDEF USE_BLOCK

pastix_iparm(IPARM_DOF_NBR)            = block_size

call pastix_fortran(pastix_data,MPI_COMM_WORLD, n_block,                                                 &
                    mumps_par%jcn(1:n_block+1), mumps_par%irn(1:nnz_block), mumps_par%A(1:mumps_par%nz), &
                    pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

#ELSE

call pastix_fortran(pastix_data,MPI_COMM_WORLD, mumps_par%n, mumps_par%jcn, mumps_par%irn, mumps_par%A, &
                    pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

#ENDIF

call cpu_time(t_fact_1)

if (my_id .eq. 0) write(*,'(A,f8.3)')  ' PASTIX, fact/solv : ',t_fact_1-t_fact_0

do k=1,mumps_par%n
  deltas(k) =  mumps_par%rhs(k)  / column_scaling(k)
enddo

if (allocated(column_local))  call tr_deallocate(column_local,"column_local")
if (allocated(counts))        call tr_deallocate(counts,"counts")
if (allocated(displacements)) call tr_deallocate(displacements,"displacements")

return
end
