#ifdef USE_STRUMPACK      
!> subroutine solves the complete system of equation using STRUMPACK
! takes distributed matrix ad_mat, centralize it and solve, placing the solution into the rhs_vec
subroutine solve_strumpack_all(spss1, ad_mat, rhs_vec)
  use strumpack_module
  use mod_strumpack

  use tr_module 
  use mod_parameters, only: n_tor, n_var
  use mpi_mod
  use mod_clock
  use mod_integer_types
  use data_structure, only: type_SP_MATRIX, type_RHS

  implicit none

! --- Local variables
  integer                           :: index_min, index_max
  integer(kind=int_all)             :: m_loc  
  type(clcktype)                    :: t_itstart, t0, t1, t2, t3
  real*8                            :: tsecond
  integer                           :: i, k, j, ierr
  integer(kind=int_all), parameter  :: Int1=1
  integer                           :: my_id, n_cpu
  
  integer(kind=C_INT_ALL)     :: n, nnz
  
  type(type_SP_MATRIX)        :: ad_mat, ac_mat
  type(type_RHS)              :: rhs_vec
  type(type_STRUMPACK_SOLVER) :: spss1
  
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)  

!write(*,*) my_id,'*********************************'
!write(*,*) my_id,'*  solve global matrix using STRUMPACK *'
!write(*,*) my_id,'*********************************'

  index_min = ad_mat%index_min
  index_max = ad_mat%index_max
  m_loc = (index_max - index_min + 1) * n_tor * n_var

  call MPI_Allreduce(m_loc,ac_mat%ng,1,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(ad_mat%nnz,ac_mat%nnz,1,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)
  
  n = ac_mat%ng
  nnz = ac_mat%nnz
  
  call clck_time(t0)

  allocate(ac_mat%irn(nnz))
  allocate(ac_mat%jcn(nnz))
  allocate(ac_mat%val(nnz))
  ac_mat%nnz = nnz

  call split_allgathersolve(n_cpu,my_id,ad_mat,ac_mat)
  
#ifdef USE_BLOCK
  ac_mat%block_size  = n_tor * n_var
#else
  ac_mat%block_size = 1
#endif  
  
  call clck_time(t1)
  call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time mpi_gather :', tsecond
  
  if (.not. spss1%initialized) then
    call strumpack_init_core(spss1, MPI_COMM_WORLD)
    spss1%initialized = .true.
    spss1%equilibrium = .false.
  endif
  
  call strumpack_set_mat_core(spss1, ac_mat)
  
  if (.not. spss1%analyzed) then
    call clck_time(t0)
    
    call strumpack_analyze_core(spss1)
    spss1%analyzed = .true.
    
    call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time analysis :', tsecond    
  endif
  
  call clck_time(t0)
  
  call strumpack_factorize_core(spss1)
  
  call strumpack_solve_core(spss1, rhs_vec)
  
  call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed facto/solve :', tsecond    
  
  return
  
  
  if (.not. spss_initialized) then
    call strumpack_init(MPI_COMM_WORLD)
    spss_initialized = .true.
  endif  
  call strumpack_set_mat(n,nnz,ac_mat%irn,ac_mat%jcn,ac_mat%val,1, MPI_COMM_WORLD,&
                         UPDATE=spss1%analyzed, DISTRIBUTED=.false.,EQUILIBRIUM=.false.)

  if (.not. spss_analyzed) then
    call clck_time(t0)
    call strumpack_analyze(MPI_COMM_WORLD)    
    spss_analyzed = .true.
    call clck_time(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time analysis :', tsecond    
  endif
  
  call clck_time(t0)

  call strumpack_factorize(MPI_COMM_WORLD)   
  call strumpack_solve(n,rhs_vec%val,MPI_COMM_WORLD)
 
  call clck_time(t1)
  call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed facto/solve :', tsecond  
 
  deallocate(ac_mat%irn)
  deallocate(ac_mat%jcn)
  deallocate(ac_mat%val)  

  return
end subroutine solve_strumpack_all
#endif  
