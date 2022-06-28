module mod_sparse
  use iso_c_binding
  use mpi
  use mumps_module, only:   mumps_par
  use phys_module, only:    use_pastix, use_mumps, use_strumpack
#ifdef USE_PASTIX
  use mod_pastix, only:     type_PASTIX_SOLVER
#endif
#ifdef USE_STRUMPACK
  use mod_strumpack, only:  type_STRUMPACK_SOLVER
#endif
  use data_structure, only: type_PRECOND

  type type_SP_SOLVER
#ifdef USE_PASTIX
    type(type_PASTIX_SOLVER)    :: ptss
#endif
#ifdef USE_STRUMPACK
    type(type_STRUMPACK_SOLVER) :: spss
#endif
    type(type_PRECOND)          :: pc
  end type type_SP_SOLVER
  
  integer, parameter :: MHD_EQUILI = 0
  integer, parameter :: MHD_DIRECT = 1
  integer, parameter :: MHD_PRECON = 2
   
  private
  public :: solve_sparse_system, MHD_EQUILI, MHD_DIRECT, MHD_PRECON, type_SP_SOLVER

  contains

!> solve Ax = rhs
!! sol_vec contains the initial guess
!! sol_vec, rhs_vec are broadcasted
!! solve_type - type of system, e.g. GS equilibrium, MHD system with preconditioner, etc.
  subroutine solve_sparse_system(a_mat, rhs_vec, solver, solve_type)

    use data_structure, only: type_SP_MATRIX, type_PRECOND, type_RHS
    use mod_clock
#ifdef USE_PASTIX
    use mod_pastix, only: type_PASTIX_SOLVER, pastix_finalize
#endif
    use mod_preconditioner, only: initialize_preconditioner
    
    implicit none
    
    type(type_SP_MATRIX)     :: a_mat
    type(type_RHS)           :: rhs_vec
    integer                  :: solve_type
    integer                  :: my_id, n_cpu, ierr
    
    type(clcktype)           :: t_itstart, t0, t1, t2, t3
    real*8                   :: tsecond
    type(type_SP_SOLVER)     :: solver

    
    if (solve_type.eq.MHD_EQUILI) then
    
      write(*,*) solve_type
      
    elseif (solve_type.eq.MHD_DIRECT) then
    
      call MPI_COMM_SIZE(a_mat%comm, n_cpu, ierr)
      call MPI_COMM_RANK(a_mat%comm, my_id, ierr)
      
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

    elseif (solve_type.eq.MHD_PRECON) then
      call MPI_COMM_SIZE(a_mat%comm, n_cpu, ierr)
      call MPI_COMM_RANK(a_mat%comm, my_id, ierr)    
      if (my_id.eq.0) write(*,*) "Solving MHD system using iterative solver"
      call initialize_preconditioner(solver%pc,a_mat%comm)
    else
      write(*,*) solve_type
    endif
    
  end subroutine solve_sparse_system
  
  
end module mod_sparse
