module mod_newton_solver

implicit none

  type type_newton_solver
    integer                       :: maxNewton    = 20
    real(kind=8)                  :: gamma_Newton = 0.5
    real(kind=8)                  :: alpha_Newton = 2.d0

    contains
    procedure :: setup
  end type type_newton_solver

contains

!> Solve sparse system using inexact Newton iterative method
  subroutine solve_newton(a_mat, rhs_vec, deltas, solver, mhd_sim, tag)
    use mod_sparse_data,      only: type_SP_SOLVER
    use mod_sparse,           only: solve_sparse_system
    use construct_matrix_mod, only: construct_matrix
    use mod_simulation_data,  only: type_MHD_SIM
    use data_structure,       only: type_RHS, node_list_save, node_list_restore
    use mod_sparse_matrix,    only: type_SP_MATRIX

    type(type_SP_SOLVER)          :: solver
    type(type_SP_MATRIX)          :: a_mat
    type(type_RHS)                :: rhs_vec, deltas
    type(type_MHD_SIM)            :: mhd_sim
    integer                       :: tag

    type(type_newton_solver)      :: newton_solver
    type(type_SP_MATRIX)          :: a_matk
    real(kind=8), allocatable     :: store_value(:,:,:,:), store_delta(:,:,:,:)
    integer                       :: inewton
    real(kind=8)                  :: normRHSn, normRHSprevious, normRHScurrent
    real(kind=8)                  :: tol_Newton, tol_Gmres
    type(type_RHS)                :: deln, rhsk, residue, dum_vec
    logical                       :: converged
    real(kind=8)                  :: solve_tol
    logical                       :: verbose = .false.

    external                      :: matrix_vector
    real(kind=8), external        :: dnrm2
    !real(kind=C_DOUBLE), external :: dnrm2

    verbose = solver%verbose.and.(tag.eq.0)
    call newton_solver%setup()

    deln%n = deltas%n
    if (associated(deln%val)) deallocate(deln%val)
    allocate(deln%val(deln%n)); deln%val(1:deln%n) = 0.d0

    rhsk%n = deltas%n
    if (associated(rhsk%val)) deallocate(rhsk%val)
    allocate(rhsk%val(rhsk%n)); rhsk%val(1:rhsk%n) = 0.d0

    dum_vec%n = deltas%n

    converged = .false.

    call a_mat%copy_to(a_matk)
    call node_list_save(mhd_sim%node_list, store_value, store_delta)

    normRHSn = dnrm2(rhs_vec%n, rhs_vec%val, 1)
    normRHScurrent = normRHSn

    solve_tol = solver%iter_tol

    newton_loop: do inewton=1, newton_solver%maxNewton

      if (inewton.gt.1) call matrix_vector(deln, a_mat, rhsk)
      normRHSprevious = normRHScurrent
      rhsk%val(1:rhsk%n) = rhs_vec%val(1:rhs_vec%n) - rhsk%val(1:rhsk%n)
      normRHScurrent = dnrm2(rhsk%n, rhsk%val, 1)

      tol_Newton = normRHScurrent/normRHSn
      tol_Gmres = newton_solver%gamma_Newton*(normRHScurrent/normRHSprevious)**newton_solver%alpha_Newton

      if (verbose) write(*,*) "Newton:", inewton, tol_Newton, tol_Gmres

      if (tol_Newton.le.(solve_tol)) then
        converged = .true.
        if (verbose) then
          write(*,*) "Newton converged:", inewton - 1, tol_Newton
        endif

        call node_list_restore(mhd_sim%node_list, store_value, store_delta)
        deltas%val(1:deltas%n) = deln%val(1:deln%n)
        solver%iter_tol = solve_tol ! restore tolerance

        call a_matk%reset()
        deallocate(deln%val,rhsk%val)

        exit newton_loop

      elseif (inewton.gt.1) then

        call update_values(mhd_sim%element_list,mhd_sim%node_list, deltas)
        call update_deltas(mhd_sim%node_list, deltas)

        allocate(dum_vec%val(dum_vec%n))
        call construct_matrix(mhd_sim, mhd_sim%local_elms, mhd_sim%n_local_elms, a_matk, dum_vec, harmonic_matrix=.false.)
        deallocate(dum_vec%val)

        if (inewton.eq.newton_solver%maxNewton) then
          write(*,*) "Newton failed to converge:", inewton - 1, tol_Newton
          solver%step_success = .false.
        endif

      endif

      solver%iter_tol = tol_Gmres
      call solve_sparse_system(a_matk, rhsk, deltas, solver, mhd_sim)
      solver%solve_only = .true.

      if (.not.solver%step_success) return

      deln%val(1:deln%n) = deln%val(1:deln%n) + deltas%val(1:deltas%n)

    enddo newton_loop

    return

  end subroutine solve_newton

!> Set Newton parameters
  subroutine setup(self)
    use phys_module, only: maxNewton, gamma_Newton, alpha_Newton
    class(type_newton_solver)     :: self

    self%maxNewton = maxNewton
    self%gamma_Newton = gamma_Newton
    self%alpha_Newton = alpha_Newton
    return
  end subroutine setup

end module mod_newton_solver