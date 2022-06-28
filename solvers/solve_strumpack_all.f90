#ifdef USE_STRUMPACK      
!> subroutine solves the complete system of equation using STRUMPACK
! takes distributed matrix ad_mat, centralize it and solve, placing the solution into the rhs_vec
subroutine solve_strumpack_all(spss, ad_mat, rhs_vec)
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
  integer                           :: my_id, n_cpu, comm
  
  integer(kind=C_INT_ALL)     :: n, nnz
  
  type(type_SP_MATRIX)        :: ad_mat, ac_mat
  type(type_RHS)              :: rhs_vec
  type(type_STRUMPACK_SOLVER) :: spss
  
  comm = ad_mat%comm
  
  call MPI_COMM_RANK(comm, my_id, ierr)
  call MPI_COMM_SIZE(comm, n_cpu, ierr)  

!write(*,*) my_id,'***************************************'
!write(*,*) my_id,'* solve global matrix using STRUMPACK *'
!write(*,*) my_id,'***************************************'

  index_min = ad_mat%index_min
  index_max = ad_mat%index_max
  m_loc = (index_max - index_min + 1) * n_tor * n_var

  call MPI_Allreduce(m_loc,ac_mat%ng,1,MPI_INTEGER_ALL,MPI_SUM,comm,ierr)
  call MPI_Allreduce(ad_mat%nnz,ac_mat%nnz,1,MPI_INTEGER_ALL,MPI_SUM,comm,ierr)
  
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
  
  if (.not. spss%initialized) then
    call strumpack_init_core(spss, comm)
    spss%initialized = .true.
    spss%equilibrium = .false.
  endif
  
  call strumpack_set_mat_core(spss, ac_mat)
  
  if (.not. spss%analyzed) then
    call clck_time(t0)
    
    call strumpack_analyze_core(spss)
    spss%analyzed = .true.
    
    call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time analysis :', tsecond    
  endif
  
  call clck_time(t0)
  
  call strumpack_factorize_core(spss)
  
  call strumpack_solve_core(spss, rhs_vec)
  
  call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed facto/solve :', tsecond    

  deallocate(ac_mat%irn)
  deallocate(ac_mat%jcn)
  deallocate(ac_mat%val)  

  return
end subroutine solve_strumpack_all
#endif  
