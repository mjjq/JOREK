module mod_sparse
  use iso_c_binding
  use mpi
  use mod_sparse_data
  use phys_module, only:    use_pastix, use_mumps, use_strumpack, iter_precon, gmres_max_iter
  

  private
  public :: solve_sparse_system

  contains

!> solve Ax = rhs
!! sol_vec contains the initial guess
!! sol_vec, rhs_vec are broadcasted
!! solve_type - type of system, e.g. GS equilibrium, MHD system with preconditioner, etc.
  subroutine solve_sparse_system(a_mat, rhs_vec, sol_vec, solver, solve_type)

    use data_structure, only: type_SP_MATRIX, type_PRECOND, type_RHS
    use mod_integer_types
    use mod_clock
    use mod_sparse_data, only: type_SP_SOLVER
#ifdef USE_PASTIX
    use mod_pastix, only: type_PASTIX_SOLVER, pastix_finalize
#endif
    use mod_preconditioner, only: initialize_preconditioner, reset_preconditioner
    use mod_distribute_preconditioner_core, only: update_pc_mat, update_pc_rhs, gather_solution
    !use mod_gmres_core, only: gmres_driver
    use mod_bicgstab_core, only: bicgstab_driver
    
    implicit none
    
    type(type_SP_MATRIX)     :: a_mat
    type(type_RHS)           :: rhs_vec, sol_vec
    integer                  :: solve_type
    integer                  :: my_id, n_cpu, ierr
    
    type(clcktype)           :: t_itstart, t0, t1, t2, t3
    real*8                   :: tsecond
    type(type_SP_SOLVER)     :: solver
    integer(kind=int_all)    :: i
    integer :: max_it = 10
    real(kind=8) :: tol = 1.e-7
    
    call MPI_COMM_SIZE(a_mat%comm, n_cpu, ierr)
    call MPI_COMM_RANK(a_mat%comm, my_id, ierr)
    sol_vec%n = rhs_vec%n
    
    if (solve_type.eq.MHD_EQUILI) then
    
      write(*,*) solve_type
      
    elseif (solve_type.eq.MHD_DIRECT) then
      
      if (my_id.eq.0) write(*,*) "Solving MHD system using direct solver"
    
      if (use_mumps) then
      
        call solve_mumps_all(my_id)
#ifdef USE_STRUMPACK        
      elseif (use_strumpack) then
      
        call solve_strumpack_all(solver%spss, a_mat, rhs_vec)
#endif
#ifdef USE_PASTIX
      elseif (use_pastix) then
      
        call solve_pastix_all(solver%ptss, a_mat, rhs_vec)
        !call pastix_finalize(solver%ptss)
#endif
      endif
      
      do i=1,rhs_vec%n
        sol_vec%val(i) =  rhs_vec%val(i)
      enddo      

    elseif (solve_type.eq.MHD_PRECON) then

      if (my_id.eq.0) write(*,*) "Solving MHD system using iterative solver"
      
! Finding PC solution
      
      if (.not.solver%pc%initialized) call initialize_preconditioner(solver%pc,a_mat%comm)
      
      call update_pc_mat(solver%pc,a_mat)
      
      call update_pc_rhs(solver%pc,rhs_vec)
      
      call solve_strumpack_all(solver%spss, solver%pc%mat, solver%pc%rhs)
      
      call gather_solution(solver%pc,sol_vec)
   
      
! iterative part

      solver%iter_prev  = iter_precon
      solver%iter_gmres = gmres_max_iter
      
      !call gmres_driver(a_mat, rhs_vec, sol_vec, solver)
      call bicgstab_driver(a_mat, sol_vec%val, rhs_vec%val, max_it, tol, solver%pc%comm, solver%pc%MPI_COMM_N, solver%pc%MPI_COMM_MASTER, solver)
      
     
      !call reset_preconditioner(solver%pc)

    endif
    
  end subroutine solve_sparse_system
 
  
end module mod_sparse
