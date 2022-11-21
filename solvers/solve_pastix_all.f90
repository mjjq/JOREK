#if defined(USE_PASTIX)||defined(USE_PASTIX6)
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

    ptss%rhs_val => rhs_vec%val

#ifdef USE_PASTIX
    call pastix_set_mat(ptss, ad_mat, ac_mat)

    if (.not. ptss%initialized) then
      ptss%comm = ad_mat%comm
      call pastix_initialize(ptss)
    endif
#elif USE_PASTIX6
    if (.not. ptss%initialized) then
      ptss%comm = ad_mat%comm
      call pastix_initialize(ptss)
    endif

    call pastix_set_mat(ptss, ad_mat, ac_mat)
#endif

    if (.not. ptss%analyzed) then
      if (my_id .eq. 0) write(*,*) "PaStiX: analyzing matrix"
      call clck_time(t0)

      call pastix_analyze(ptss)

      call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
      if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time analysis:', tsecond
    endif

    if (my_id .eq. 0) write(*,*) "PaStiX: factorizing matrix"
    call clck_time(t0)

    call pastix_factorize(ptss)

    call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time factorization:', tsecond

    if (n_cpu>1) then
      deallocate(ac_mat%irn)
      deallocate(ac_mat%jcn)
      deallocate(ac_mat%val)
    endif

  endif

  call clck_time(t0)

  call pastix_solve(ptss,rhs_vec)

  call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time solve:', tsecond

  return
end
#endif
