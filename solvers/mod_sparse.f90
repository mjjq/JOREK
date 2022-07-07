module mod_sparse
  use mod_sparse_data
  use phys_module, only: use_mumps, use_pastix, use_strumpack, &
                         use_mumps_eq, use_pastix_eq, use_strumpack_eq
  

  private
  public :: solve_sparse_system, solver_finalize

  contains

!> solve Ax = rhs
!! sol_vec contains the initial guess
!! sol_vec, rhs_vec are broadcasted
!! solve_type - type of system, e.g. GS equilibrium, MHD system with preconditioner, etc.
  subroutine solve_sparse_system(a_mat, rhs_vec, sol_vec, solver)

    use data_structure, only: type_SP_MATRIX, type_PRECOND, type_RHS
    use mod_integer_types
    use mod_clock
    use mod_sparse_data, only: type_SP_SOLVER
#ifdef USE_PASTIX
    use mod_pastix, only: type_PASTIX_SOLVER
#endif
    use mod_preconditioner, only: initialize_preconditioner, reset_preconditioner
    use mod_distribute_preconditioner_core, only: update_pc_mat, update_pc_rhs, gather_solution
    !
#ifdef USE_BICGSTAB    
    use mod_bicgstab_core, only: bicgstab_driver
#else
    use mod_gmres_core, only: gmres_driver
#endif
    
    implicit none
    
    type(type_SP_MATRIX)     :: a_mat
    type(type_RHS)           :: rhs_vec, sol_vec
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

    if (.not.solver%iterative) then
    
      if (solver%equilibrium) then
        if (my_id.eq.0) write(*,*) "Solving MHD equilibrium system"
      else
        if (my_id.eq.0) write(*,*) "Solving MHD system using direct solver"
      endif
    
      if ((use_mumps.and..not.solver%equilibrium).or.(use_mumps_eq.and.solver%equilibrium)) then
#ifdef USE_MUMPS
        if (my_id.eq.0) write(*,*) "Using MUMPS solver"
        solver%mmss%equilibrium = solver%equilibrium
        !call solve_mumps_all_core(solver%mmss, a_mat, rhs_vec, solver%solve_only)
        call solve_mumps_all(solver%mmss, a_mat, rhs_vec, solver%solve_only)
#endif        
      elseif ((use_strumpack.and..not.solver%equilibrium).or.(use_strumpack_eq.and.solver%equilibrium)) then        
#ifdef USE_STRUMPACK
        if (my_id.eq.0) write(*,*) "Using STRUMPACK solver"
        solver%spss%equilibrium = solver%equilibrium
        call solve_strumpack_all(solver%spss, a_mat, rhs_vec, solver%solve_only)
#endif
      elseif ((use_pastix.and..not.solver%equilibrium).or.(use_pastix_eq.and.solver%equilibrium)) then
#ifdef USE_PASTIX
        if (my_id.eq.0) write(*,*) "Using PaStiX solver"
        solver%ptss%equilibrium = solver%equilibrium
        solver%ptss%refine = .true.
        call solve_pastix_all(solver%ptss, a_mat, rhs_vec, solver%solve_only)
#endif
      endif
      
      do i=1,rhs_vec%n
        sol_vec%val(i) =  rhs_vec%val(i)
      enddo
      
      solver%step_success = .true.

    elseif (solver%iterative) then

      if (my_id.eq.0) write(*,*) "Solving MHD system using iterative solver"

      
      ! condition for no PC update
      solver%solve_only = (solver%istep > 1) .and. ((solver%iter_gmres + solver%iter_prev <= 2*solver%iter_precon) &
                                             .and.  (solver%n_since_update < solver%max_steps_noUpdate))      
      
      if (solver%solve_only) then 
        solver%n_since_update = solver%n_since_update + 1
      else
        solver%n_since_update = 0
      endif         
      
      if (.not.solver%pc%initialized) call initialize_preconditioner(solver%pc,a_mat%comm)
      
! Finding PC solution
      if (.not.solver%solve_only) then
        call update_pc_mat(solver%pc,a_mat)
      endif
      
      call update_pc_rhs(solver%pc,rhs_vec)
        
      if (use_mumps) then
      
        !call solve_mumps_all(solver%mpss, solver%pc%mat, solver%pc%rhs, solver%solve_only)
        
#ifdef USE_STRUMPACK        
      elseif (use_strumpack) then
      
        call solve_strumpack_all(solver%spss, solver%pc%mat, solver%pc%rhs, solver%solve_only)
#endif
#ifdef USE_PASTIX
      elseif (use_pastix) then
      
        call solve_pastix_all(solver%ptss, solver%pc%mat, solver%pc%rhs, solver%solve_only)

#endif
      endif      
      
      call gather_solution(solver%pc,sol_vec)
      
! iterative part
      solver%iter_prev  = solver%iter_gmres
      solver%iter_gmres = solver%iter_max

#ifdef USE_BICGSTAB      
      call bicgstab_driver(a_mat, rhs_vec, sol_vec, solver)
#else
      call gmres_driver(a_mat, rhs_vec, sol_vec, solver)
#endif

      if (my_id.eq.0) write(*,'(A32,I5)') 'Number of iterations: ', solver%iter_gmres
      
      solver%step_success = (solver%iter_gmres .lt. solver%iter_max)

    endif
    
  end subroutine solve_sparse_system
  
  !call pastix_finalize(solver%ptss)     
  !call reset_preconditioner(solver%pc)
  
  subroutine solver_finalize(solver)
#ifdef USE_PASTIX
    use mod_pastix, only: pastix_finalize
#endif
#ifdef USE_MUMPS
    use mod_mumps, only: mumps_finalize
#endif
    implicit none
      
    type(type_SP_SOLVER)     :: solver
    
    write(*,*) "Finalizing solver"    
    
#ifdef USE_PASTIX
    if (solver%ptss%initialized) call pastix_finalize(solver%ptss)
#endif
#ifdef USE_MUMPS
    if (solver%mmss%initialized) call mumps_finalize(solver%mmss)
#endif
    
    solver%solve_only   = .false.
    solver%step_success = .false.
    solver%iterative    = .false.
    solver%equilibrium  = .false.


    return
  end subroutine solver_finalize
 
  
end module mod_sparse
