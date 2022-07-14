#ifdef USE_STRUMPACK      
!> subroutine solves the complete system of equation using STRUMPACK
! takes distributed matrix ad_mat, centralize it and solve, placing the solution into the rhs_vec
subroutine solve_strumpack_all(spss, ad_mat, rhs_vec, solve_only)
  use mod_strumpack

  use tr_module 
  use mpi_mod
  use mod_clock
  use mod_integer_types
  use data_structure, only: type_SP_MATRIX, type_RHS

  implicit none

  type(type_SP_MATRIX)        :: ad_mat, ac_mat
  type(type_RHS)              :: rhs_vec
  type(type_STRUMPACK_SOLVER) :: spss
  
  
! --- Local variables
  type(clcktype)              :: t_itstart, t0, t1, t2, t3
  real*8                      :: tsecond
  integer                     :: my_id, n_cpu, comm, ierr
  logical                     :: solve_only
  logical                     :: centralize = .true.
  
  if (.not.solve_only) then
    
    comm = ad_mat%comm
    
    call MPI_COMM_RANK(comm, my_id, ierr)
    call MPI_COMM_SIZE(comm, n_cpu, ierr)  
  
  !write(*,*) my_id,'***************************************'
  !write(*,*) my_id,'* solve sparse matrix using STRUMPACK *'
  !write(*,*) my_id,'***************************************'
    
    centralize = (n_cpu>1).and.(.not.ad_mat%row_distributed)
  
    if (centralize) then
  ! centralize distributed matrix ad_mat into ac_mat
    
      call MPI_Allreduce(ad_mat%nnz,ac_mat%nnz,1,MPI_INTEGER_ALL,MPI_SUM,comm,ierr)
    
      ac_mat%ng = ad_mat%ng
      ac_mat%block_size = ad_mat%block_size
      ac_mat%comm = ad_mat%comm
      
      call clck_time(t0)
    
      allocate(ac_mat%irn(ac_mat%nnz))
      allocate(ac_mat%jcn(ac_mat%nnz))
      allocate(ac_mat%val(ac_mat%nnz))
    
      call split_allgathersolve(n_cpu,my_id,ad_mat,ac_mat)
      
      call clck_time(t1)
      call clck_ldiff(t0,t1,tsecond)
      if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time mpi_gather :', tsecond
      
    else
      ac_mat = ad_mat
    endif
    
    if (.not. spss%initialized) then
      call strumpack_init_core(spss, comm)
      spss%initialized = .true.
    endif
    
    call strumpack_set_mat_core(spss, ac_mat)
    
    if (.not. spss%analyzed) then
      call clck_time(t0)
      
      call strumpack_analyze_core(spss)
      spss%analyzed = .true.
      
      call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
      if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time analysis:', tsecond    
    endif
    
    call clck_time(t0)
    
    call strumpack_factorize_core(spss)
    
    call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time factorize:', tsecond
    
    if (centralize) then
      deallocate(ac_mat%irn)
      deallocate(ac_mat%jcn)
      deallocate(ac_mat%val)
    endif    
    
  endif ! .not.solve_only
  
  call clck_time(t0)
  
  call strumpack_solve_core(spss, rhs_vec)
  
  call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time solve:', tsecond    



  return
end subroutine solve_strumpack_all
#endif  
