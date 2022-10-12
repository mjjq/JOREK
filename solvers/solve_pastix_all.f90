#ifdef USE_PASTIX
!> subroutine solves the complete system of equation using pastix with
!  distributed matrix ad_mat on the main group mpi_comm_world.
!  For pastix5 solver matrix is centralized into ac_mat
subroutine solve_pastix_all(ptss, ad_mat, rhs_vec, solve_only)
  use tr_module 
  use mod_parameters, only: n_tor, n_var
  use mpi_mod
  use mod_clock
  use mod_coicsr
  
  use mod_integer_types
  use data_structure, only: type_SP_MATRIX, type_RHS
  use mod_pastix, only:     type_PASTIX_SOLVER, pastix_solve, pastix_factorize, pastix_analyze, pastix_initialize
  
#include "pastix_fortran.h"  
   
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
    
    call clck_time(t0)
  
    block_size2 = ac_mat%block_size**2
    
    nblock   = ac_mat%ng/ac_mat%block_size
    nnz_block = ac_mat%nnz/block_size2
    
    ac_mat%nblock = nblock
    ac_mat%nzblock = nnz_block
    
    if (ac_mat%block_size > 1) then
      do i = 1,nnz_block  
        ac_mat%irn(i) = (ac_mat%irn((i-1)*block_size2+1) - 1)/ac_mat%block_size + 1 
        ac_mat%jcn(i) = (ac_mat%jcn((i-1)*block_size2+1) - 1)/ac_mat%block_size + 1 
      enddo
    endif
    
    allocate(sparskit_work(nblock+1))
    
    call coicsr2(nblock,nnz_block,ac_mat%val,ac_mat%irn(1:nnz_block),ac_mat%jcn(1:nnz_block),ac_mat%block_size,sparskit_work)
    
    deallocate(sparskit_work)
    
    call clck_time(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time coicsr :', tsecond
    
    ! End of matrix preparation
    ptss%irn => ac_mat%irn
    ptss%jcn => ac_mat%jcn
    ptss%val => ac_mat%val
    ptss%comm = ac_mat%comm
    ptss%rhs_val => rhs_vec%val
    ptss%nblock = ac_mat%nblock
    
    
    if (.not. ptss%initialized) then
    
      call pastix_initialize(ptss)
      
    endif
    
    if (ptss%equilibrium) then
    ! combine duplicated values
      nnz = ac_mat%jcn(ac_mat%ng + 1) - 1
      call pastix_fortran_checkmatrix(check_data, ac_mat%comm, &
       Int1, ptss%sym, Int1, ac_mat%ng, ac_mat%jcn, ac_mat%irn, ac_mat%val, -Int1, Int1)

      ac_mat%nnz = ac_mat%jcn(ac_mat%ng+1) - 1
      if (ac_mat%nnz < nnz) then
         call pastix_fortran_checkmatrix_end(check_data, Int1, ac_mat%irn, ac_mat%val, Int1)
      endif
    endif
    
    if (.not. ptss%analyzed) then
    
      ptss%iparm(IPARM_DOF_NBR) = ac_mat%block_size           
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
