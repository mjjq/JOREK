subroutine solve_pastix_all(n_cpu,my_id,index_min,index_max)
!---------------------------------------------------------------------
! subroutine solves the complete system of equation using pastix with
! distributed matrix on the main group mpi_comm_world
!---------------------------------------------------------------------
use parameters
use mumps_module
use pastix_module
use global_distributed_matrix
implicit none
include 'mpif.h'

integer                  :: n_cpu, index_min, index_max       ! global index_min, index_max for this cpu
real*8,allocatable       :: column_local(:)
integer, allocatable     :: pastix_loc2glb(:)
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

if (allocated(column_scaling))  deallocate(column_scaling)
if (allocated(column_local))    deallocate(column_local)
allocate(column_scaling(mumps_par%N),column_local(mumps_par%N))

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
if (allocated(counts))        deallocate(counts)
if (allocated(displacements)) deallocate(displacements)

allocate(counts(n_cpu),displacements(n_cpu))

call MPI_Allgather(mumps_par%nz_loc,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

displacements(1) = 0
do i=2,n_cpu
  displacements(i) = displacements(i-1) + counts(i-1)
enddo

if (associated(mumps_par%IRN)) deallocate(mumps_par%IRN)
if (associated(mumps_par%JCN)) deallocate(mumps_par%JCN)
if (associated(mumps_par%A) )  deallocate(mumps_par%A)
if (associated(mumps_par%rhs)) deallocate(mumps_par%rhs)

allocate(mumps_par%IRN(mumps_par%nz),mumps_par%JCN(mumps_par%nz),mumps_par%A(mumps_par%nz))
allocate(mumps_par%rhs(mumps_par%n))

call MPI_AllgatherV(IRN_glob,mumps_par%nz_loc,MPI_INTEGER,mumps_par%IRN, &
                    counts,displacements,MPI_INTEGER,MPI_COMM_WORLD,ierr)

call MPI_AllgatherV(JCN_glob,mumps_par%nz_loc,MPI_INTEGER,mumps_par%JCN, &
                    counts,displacements,MPI_INTEGER,MPI_COMM_WORLD,ierr)

call MPI_AllgatherV(A_glob,mumps_par%nz_loc,MPI_DOUBLE_PRECISION,mumps_par%A, &
                    counts,displacements,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

call MPI_AllReduce(RHS_glob,mumps_par%RHS,mumps_par%N,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

if (allocated(sparskit_work)) deallocate(sparskit_work)
allocate(sparskit_work(mumps_par%N + 1))

call coicsr(mumps_par%N,mumps_par%NZ,1,mumps_par%A,mumps_par%IRN,mumps_par%JCN,sparskit_work)

call cpu_time(t_comm_1)

if (my_id .eq. 0)  write(*,'(A,f8.3)') ' PASTIX, comm      : ',t_comm_1-t_comm_0


if (.not. allocated(pastix_perm_vars))  allocate(pastix_perm_vars(mumps_par%n))
if (.not. allocated(pastix_iperm_vars)) allocate(pastix_iperm_vars(mumps_par%n))

!$omp parallel default(none) shared(pastix_nthrd)    
        pastix_nthrd = omp_get_num_threads()
!$omp end parallel

write(*,'(i5,A,i5)') my_id,' PastiX n_threads : ',pastix_nthrd 

if (.not. pastix_initialised) then

  pastix_iparm(1)  = 0          ! insert default values
  pastix_iparm(2)  = 0          ! initializse
  pastix_iparm(3)  = 0

write(*,*) '***********************************'
write(*,*) '* initialise PastiX               *'
write(*,*) '***********************************'

  call pastix_fortran(pastix_data,MPI_COMM_WORLD,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

  pastix_iparm(4) = pastix_verb
  pastix_initialised = .true.

endif

if (.not. pastix_analysed) then

  pastix_iparm(2) = 1
  pastix_iparm(3) = 3
  pastix_iparm(6) = pastix_iter          ! refinement : max number of iterations

 ! pastix_iparm(5) = n_tor * n_var   ! degrees of freedom per node (not correct)

  pastix_iparm(31) = pastix_facto
  pastix_iparm(35) = pastix_nthrd         ! numthreads : number of threads
  pastix_iparm(39) = pastix_rhs         ! right hand side (0 : use RHS)
  pastix_iparm(37) = pastix_iluk 
  pastix_iparm(41) = pastix_sym

  pastix_iparm(42) = pastix_ricar
  pastix_iparm(37) = pastix_iluk
  pastix_iparm(14) = pastix_amalg

  pastix_dparm(6)  = pastix_epsilon    ! error level refinement
  pastix_dparm(11) = pastix_pivot    ! pivot threshold?

  call cpu_time(t_analysis_0)

write(*,*) '***********************************'
write(*,*) '* analyse PastiX                  *'
write(*,*) '***********************************'

  call pastix_fortran(pastix_data,MPI_COMM_WORLD, mumps_par%n, &
                      mumps_par%jcn(1:mumps_par%n+1), mumps_par%irn(1:mumps_par%nz), mumps_par%A(1:mumps_par%nz), &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

  call cpu_time(t_analysis_1)

  pastix_analysed = .true.

  if (my_id .eq. 0)  write(*,'(A,f8.3)') ' PASTIX, analysis  : ',t_analysis_1-t_analysis_0

endif

call cpu_time(t_fact_0)

pastix_iparm(2) = 4
pastix_iparm(3) = pastix_endsolve
pastix_iparm(6) = pastix_iter          ! refinement : max number of iterations

pastix_iparm(31) = pastix_facto
pastix_iparm(35) = pastix_nthrd         !   numthreads : number of threads
pastix_iparm(39) = pastix_rhs        ! right hand side (0 : use RHS)
pastix_iparm(37) = pastix_iluk
pastix_iparm(41) = pastix_sym

pastix_iparm(42) = pastix_ricar
pastix_iparm(37) = pastix_iluk
pastix_iparm(14) = pastix_amalg

pastix_dparm(6)  = pastix_epsilon    ! error level refinement
pastix_dparm(11) = pastix_pivot    ! pivot threshold?

write(*,*) '***********************************'
write(*,*) '* call PastiX                     *'
write(*,*) '***********************************'

call pastix_fortran(pastix_data,MPI_COMM_WORLD, mumps_par%n, mumps_par%jcn, mumps_par%irn, mumps_par%A, &
                    pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

call cpu_time(t_fact_1)

if (my_id .eq. 0) write(*,'(A,f8.3)')  ' PASTIX, fact/solv : ',t_fact_1-t_fact_0

do k=1,mumps_par%n
  deltas(k) =  mumps_par%rhs(k)  / column_scaling(k)
!  write(*,*) k,deltas(k)
enddo

return
end
