#ifdef USE_PASTIX
!> subroutine solves the complete system of equation using pastix with
!  distributed matrix ad_mat on the main group mpi_comm_world.
!  For pastix5 solver matrix is centralized into ac_mat
subroutine solve_pastix_all(ptss, ad_mat, rhs_vec)
  use tr_module 
  use mod_parameters, only: n_tor, n_var
  use global_distributed_matrix
  use mpi_mod
  use mod_clock
  use mod_coicsr
  
  use mod_integer_types
  use data_structure, only: type_SP_MATRIX, type_RHS
  use mod_pastix, only:     type_PASTIX_SOLVER, scale_by_cols, pastix_solve, pastix_factorize, pastix_analyze, pastix_initialize
  
#include "pastix_fortran.h"  
   
  implicit none

  type(clcktype)                    :: t_itstart, t0, t1, t2, t3
  real*8                            :: tsecond
  integer                           :: n_cpu, my_id, ierr, comm
  integer                           :: i, k, j
    
  type(type_SP_MATRIX)               :: ad_mat, ac_mat
  type(type_RHS)                     :: rhs_vec
  integer                            :: index_min, index_max
  integer                            :: block_size2
  integer(kind=int_all)              :: n_block, nnz_block
  type(type_PASTIX_SOLVER) :: ptss
  
  integer(kind=int_all), allocatable         :: sparskit_work(:)
  
  comm = ad_mat%comm
  
  call MPI_COMM_RANK(comm, my_id, ierr)
  call MPI_COMM_SIZE(comm, n_cpu, ierr)  

  !write(*,*) my_id,'*********************************'
  !write(*,*) my_id,'*  solve global matrix (PaStiX) *'
  !write(*,*) my_id,'*********************************'
  
  call MPI_Allreduce(ad_mat%nnz,ac_mat%nnz,1,MPI_INTEGER_ALL,MPI_SUM,comm,ierr)
  
  ac_mat%ng = ad_mat%ng
  ac_mat%block_size = ad_mat%block_size
  ac_mat%comm = ad_mat%comm  
  
  allocate(ac_mat%irn(ac_mat%nnz))
  allocate(ac_mat%jcn(ac_mat%nnz))
  allocate(ac_mat%val(ac_mat%nnz))
  
  call scale_by_cols(ad_mat)
  
  call clck_time(t0)
  
  call split_allgathersolve(n_cpu,my_id,ad_mat,ac_mat)
  
  call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time mpi_gather :', tsecond
  
  call clck_time(t0)

  block_size2 = ac_mat%block_size**2
  
  n_block   = ac_mat%ng/ac_mat%block_size
  nnz_block = ac_mat%nnz/block_size2
  
  ac_mat%nblock = n_block
  ac_mat%nzblock = nnz_block
  
  if (ac_mat%block_size > 1) then
    do i=1,nnz_block  
      ac_mat%irn(i) = (ac_mat%irn((i-1)*block_size2+1) - 1)/ac_mat%block_size + 1 
      ac_mat%jcn(i) = (ac_mat%jcn((i-1)*block_size2+1) - 1)/ac_mat%block_size + 1 
    enddo
  endif
  
  allocate(sparskit_work(n_block+1))
  
  call coicsr2(n_block,nnz_block,ac_mat%val,ac_mat%irn(1:nnz_block),ac_mat%jcn(1:nnz_block),ac_mat%block_size,sparskit_work)
  
  if (allocated(sparskit_work)) deallocate(sparskit_work)
  
  call clck_time(t1)
  call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time coicsr :', tsecond
  
  ! End of matrix preparation
  
  if (.not. ptss%initialized) then
  
    call pastix_initialize(ptss, comm)
    
  endif
  
  if (.not. ptss%analyzed) then
    
    ptss%iparm(IPARM_DOF_NBR)    = ac_mat%block_size
    ptss%nblock = n_block
    
    call pastix_analyze(ptss,ac_mat)
  
  endif
  
  call pastix_factorize(ptss,ac_mat)

  deallocate(ac_mat%irn)
  deallocate(ac_mat%jcn)
  deallocate(ac_mat%val)                      
  
  call pastix_solve(ptss,rhs_vec)
  
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
