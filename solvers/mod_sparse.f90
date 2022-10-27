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
  subroutine solve_sparse_system(a_mat, rhs_vec, sol_vec, solver, sim)
    use mod_integer_types
    use mod_clock
    use data_structure, only: type_SP_MATRIX, type_PRECOND, type_RHS, type_MHD_SIM
    use mod_sparse_data, only: type_SP_SOLVER
    use mod_preconditioner, only: initialize_preconditioner, reset_preconditioner, update_pc_rhs, gather_solution
    use mod_distribute_preconditioner, only: update_pc_mat
#ifdef USE_BICGSTAB
    use mod_bicgstab, only: bicgstab_driver
#else
    use mod_gmres, only: gmres_driver
#endif

    implicit none

    type(type_SP_SOLVER)     :: solver
    type(type_SP_MATRIX)     :: a_mat
    type(type_RHS)           :: rhs_vec, sol_vec
    type(type_MHD_SIM), optional  :: sim
    
    integer                  :: my_id, n_cpu, ierr
    type(clcktype)           :: t_itstart, t0, t1, t2, t3
    real*8                   :: tsecond
    integer(kind=int_all)    :: i

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
        call solve_mumps_all(solver%mmss, a_mat, rhs_vec, solver%solve_only)
#endif
      elseif ((use_strumpack.and..not.solver%equilibrium).or.(use_strumpack_eq.and.solver%equilibrium)) then
#ifdef USE_STRUMPACK
        if (my_id.eq.0) write(*,*) "Using STRUMPACK solver"
        solver%spss%equilibrium = solver%equilibrium
        call solve_strumpack_all(solver%spss, a_mat, rhs_vec, solver%solve_only)
#endif
      elseif ((use_pastix.and..not.solver%equilibrium).or.(use_pastix_eq.and.solver%equilibrium)) then
#if (defined USE_PASTIX) || (defined USE_PASTIX6)
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
        call update_pc_mat(solver%pc,a_mat,sim)
      endif

      call update_pc_rhs(solver%pc,rhs_vec)

      if (use_mumps) then
#ifdef USE_MUMPS
        call solve_mumps_all(solver%mmss, solver%pc%mat, solver%pc%rhs, solver%solve_only)
#endif
      elseif (use_strumpack) then
#ifdef USE_STRUMPACK
        call solve_strumpack_all(solver%spss, solver%pc%mat, solver%pc%rhs, solver%solve_only)
#endif
      elseif (use_pastix) then
#if (defined USE_PASTIX) || (defined USE_PASTIX6)
        call solve_pastix_all(solver%ptss, solver%pc%mat, solver%pc%rhs, solver%solve_only)
#endif
      endif

      call MPI_Barrier(a_mat%comm, ierr)
      
      call gather_solution(solver%pc,sol_vec)
      
      !write(*,*) sol_vec%val(1:2), sol_vec%val(sol_vec%n-1:sol_vec%n)
      !call MPI_Barrier(MPI_COMM_WORLD,ierr); call MPI_Finalize(ierr); call exit(0)

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
#if (defined USE_PASTIX) || (defined USE_PASTIX6)
    use mod_pastix, only: pastix_finalize
#endif
#ifdef USE_MUMPS
    use mod_mumps, only: mumps_finalize
#endif
    implicit none

    type(type_SP_SOLVER)     :: solver

    write(*,*) "Finalizing solver"

#if (defined USE_PASTIX) || (defined USE_PASTIX6)
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
