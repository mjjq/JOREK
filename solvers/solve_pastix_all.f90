#ifdef USE_PASTIX
!> subroutine solves the complete system of equation using pastix with
!  distributed matrix ad_mat on the main group mpi_comm_world.
!  For pastix5 solver matrix is centralized into ac_mat
subroutine solve_pastix_all(ptss, ad_mat, rhs_vec, solve_only)
  use tr_module 
  use mod_parameters, only: n_tor, n_var
  use mpi_mod
  use mod_clock
  
  use mod_integer_types
  use data_structure, only: type_SP_MATRIX, type_RHS
  use mod_pastix, only:     type_PASTIX_SOLVER, pastix_set_mat, pastix_solve, pastix_factorize, pastix_analyze, pastix_initialize
   
  implicit none

  type(clcktype)                    :: t_itstart, t0, t1, t2, t3
  real*8                            :: tsecond
  integer                           :: n_cpu, my_id, ierr, comm
  type(type_SP_MATRIX)              :: ad_mat, ac_mat
  type(type_RHS)                    :: rhs_vec
  type(type_PASTIX_SOLVER)          :: ptss
  logical                           :: solve_only
  
  comm = ad_mat%comm
  
  call MPI_COMM_RANK(comm, my_id, ierr)
  call MPI_COMM_SIZE(comm, n_cpu, ierr)
  
  if (.not.solve_only) then
  
    !write(*,*) my_id,'*********************************'
    !write(*,*) my_id,'*  solve global matrix (PaStiX) *'
    !write(*,*) my_id,'*********************************'
    
    if (.not. ptss%initialized) then
    
      call pastix_initialize(ptss)
      
    endif
    
    call pastix_set_mat(ptss, ad_mat, ac_mat)
    
    ptss%rhs_val => rhs_vec%val
    
    
    if (.not. ptss%analyzed) then
    
      call pastix_analyze(ptss)
    
    endif
    
    call pastix_factorize(ptss)
    
    if (n_cpu>1) then
      deallocate(ac_mat%irn)
      deallocate(ac_mat%jcn)
      deallocate(ac_mat%val)
    endif
  endif ! .not.solve_only
  
  call clck_time(t0)
  
  call pastix_solve(ptss,rhs_vec)
  
  call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time solve:', tsecond    
  
  return
end




#elif USE_PASTIX6
!> subroutine solves the complete system of equation using PaStiX 6.2
subroutine solve_pastix_all(ptss, ad_mat, rhs_vec, solve_only)

  use tr_module 
  use mod_parameters, only: n_tor, n_var
  use mpi_mod
  use mod_clock
  use mod_coicsr
  
  use mod_integer_types
  use data_structure, only: type_SP_MATRIX, type_RHS
  use mod_pastix, only:     type_PASTIX_SOLVER, pastix_solve, pastix_factorize, pastix_analyze, pastix_initialize
  
  
  implicit none

  type(clcktype)                    :: t_itstart, t0, t1, t2, t3
  real*8                            :: tsecond
  integer                           :: n_cpu, my_id, ierr, comm
  integer                           :: i, j
  integer(kind=int_all)             :: k, nnz
  integer*8 :: check_data
    
  type(type_SP_MATRIX)               :: ad_mat, ac_mat
  type(type_RHS)                     :: rhs_vec
  integer                            :: index_min, index_max
  integer                            :: block_size2
  integer(kind=int_all)              :: nblock, nnz_block
  type(type_PASTIX_SOLVER)           :: ptss
  logical                            :: solve_only
  
  integer(kind=int_all), allocatable         :: sparskit_work(:)
  
  comm = ad_mat%comm
  
  call MPI_COMM_RANK(comm, my_id, ierr)
  call MPI_COMM_SIZE(comm, n_cpu, ierr)
  
  if (.not.solve_only) then
  
    !write(*,*) my_id,'*********************************'
    !write(*,*) my_id,'*  solve global matrix (PaStiX) *'
    !write(*,*) my_id,'*********************************'
    
    if (.not.ptss%equilibrium) then
    
      call scale_by_cols(ad_mat)
      if (associated(ptss%solution_scaling)) then
        deallocate(ptss%solution_scaling); ptss%solution_scaling => Null()
      endif
      allocate(ptss%solution_scaling(ad_mat%ng))
      do k = 1, ad_mat%ng
        ptss%solution_scaling(k) = ad_mat%column_scaling(k)
      enddo
      ptss%scaled = .true.
      
    endif
    
    if (n_cpu>1) then
    
      call MPI_Allreduce(ad_mat%nnz,ac_mat%nnz,1,MPI_INTEGER_ALL,MPI_SUM,comm,ierr)
      
      ac_mat%ng = ad_mat%ng
      ac_mat%block_size = ad_mat%block_size
      ac_mat%comm = ad_mat%comm  
      
      allocate(ac_mat%irn(ac_mat%nnz))
      allocate(ac_mat%jcn(ac_mat%nnz))
      allocate(ac_mat%val(ac_mat%nnz))
      
      call clck_time(t0)
      
      call split_allgathersolve(n_cpu,my_id,ad_mat,ac_mat)
      
      call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
      if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time mpi_gather :', tsecond
      
    else
      ac_mat = ad_mat
    endif
    
  endif
    
  call exit(0)


  
  return
  
end subroutine solve_pastix_all
#endif
