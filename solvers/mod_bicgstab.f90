module mod_bicgstab
#ifdef USE_BICGSTAB
  use iso_c_binding
  use mpi
  use mumps_module, only: mumps_par
  !use mod_settings, only: n_tor
  use mod_parameters, only : n_tor, n_var
  use phys_module, only: use_pastix, use_mumps, use_strumpack
  use global_distributed_matrix, only: ndof_glob, nz_glob, local_index_start


  implicit none

  type SPARSE_MATRIX_T
    integer(kind=C_INT), pointer :: irn(:), jcn(:)
    real(kind=C_DOUBLE), pointer :: val(:)
    integer                      :: indexing
    integer(kind=C_INT)          :: n
    integer(kind=C_INT)          :: nnz
  end type SPARSE_MATRIX_T

  type(SPARSE_MATRIX_T)          :: cooA

  integer                        :: my_id, my_id_n
  integer                        :: MPI_GLOB, MPI_COMM_N, MPI_COMM_MASTER


  private
  public :: bicgstab_driver

  contains

!> solve Ax = b using iterative preconditioned BiCGStab method
!! x=sol contains the initial guess
!! max_it - maximum number of iterations
!! tol - iteration tolerance
  subroutine bicgstab_driver(irn, jcn, val, x, b, max_it, tol, comm_glob, comm_n, comm_master)
    implicit none

    integer(kind=C_INT), pointer, intent(in)    :: irn(:), jcn(:)
    real(kind=C_DOUBLE), pointer, intent(in)    :: val(:)
    real(kind=C_DOUBLE), allocatable            :: b(:)
    real(kind=C_DOUBLE), allocatable            :: x(:)
    !real(kind=8), allocatable, intent(inout)    :: sol(:)
    integer, intent(in)                         :: comm_glob, comm_n, comm_master
    integer, intent(inout)                      :: max_it

    !real(kind=C_DOUBLE), allocatable :: x(:) ! local copy of sol, not necessary, but sol is not a pointer
    real(kind=C_DOUBLE), allocatable, target :: r(:), r_tld(:), s(:), s_hat(:), tmp(:), p(:), p_hat(:), v(:), t(:)
    integer                          :: iter, flag
    integer                          :: ierr
    real(kind=C_DOUBLE)              :: tol, error, alpha, beta, omega, bnrm2, rho, rho_1, resid, snrm2
    integer                          :: n, nnz

    real(kind=C_DOUBLE), external :: dnrm2, ddot

    ! set module values
    n = ndof_glob
    nnz = nz_glob
    MPI_GLOB = comm_glob
    MPI_COMM_N = comm_n
    MPI_COMM_MASTER = comm_master
    call MPI_COMM_RANK(MPI_GLOB, my_id, ierr)
    call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr)

    cooA%irn => irn
    cooA%jcn => jcn
    cooA%val => val
    cooA%indexing = 1
    cooA%n = n
    cooA%nnz = nnz

    iter = 0
    flag = 0

    allocate(r(n), r_tld(n), s(n), s_hat(n), tmp(n), p(n), p_hat(n), v(n), t(n))

    bnrm2 = dnrm2(n, b, 1)

    if (bnrm2 == 0.0) bnrm2 = 1.0

    call matv(x,tmp)
    r = b - tmp

    error = dnrm2(n, r, 1)/bnrm2;
    if (error < tol) then
      if (my_id.eq.0) write(*,*) "bicgstab exiting, initial relative error:", error
      max_it = 1
      return
    endif

    omega  = 1.0
    r_tld = r

    do iter = 1, max_it
      rho = ddot(n,r_tld,1,r,1)

      if (rho == 0.0 ) exit

      if ( iter > 1 ) then
        beta  = (rho/rho_1)*(alpha/omega)
        p = r + beta*(p - omega*v)
      else
        p = r
      endif

      call pc(p,p_hat)
      call matv(p_hat,v)
      alpha = rho/ddot(n,r_tld,1,v,1)
      s = r - alpha*v
      snrm2 = dnrm2(n, s, 1)

      if (snrm2 < tol) then
        x = x + alpha*p_hat
        resid = snrm2/bnrm2
        exit
      endif

      call pc(s,s_hat)
      call matv(s_hat,t)
      omega = ddot(n,t,1,s,1)/ddot(n,t,1,t,1)
      x = x + alpha*p_hat + omega*s_hat
      r = s - omega*t
      error = dnrm2(n, r, 1)/bnrm2

      if (error <= tol) exit
      if (omega == 0.0) exit

      rho_1 = rho
    enddo

    if ((error <= tol).or.(snrm2 <= tol)) then ! converged
     if ( snrm2 <= tol ) error = snrm2/bnrm2;
     flag =  0
     if (my_id.eq.0) write(*,*) "bicgstab completed successfully with n_iter: ", iter
    elseif (omega == 0.0) then ! breakdown
     flag = -2
     if (my_id.eq.0) write(*,*) "bicgstab fails, flag: ", flag
    elseif (rho == 0.0) then
     flag = -1
     write(*,*) "bicgstab fails, flag: ", flag
    else ! no convergence
     flag = 1
     if (my_id.eq.0) write(*,*) "bicgstab failed to converge, relative error: ", error
    endif

    !write(*,*) "solution", x(1), x(n)
    max_it = iter - 1

    deallocate(r, r_tld, s, s_hat, tmp, p, p_hat, v, t)

  end subroutine bicgstab_driver

!> get matrix-vector product b=Ax
  subroutine matv(x,b)
    implicit none

    real(kind=C_DOUBLE), allocatable  :: x(:), b(:)
    integer                           :: i, j, ir, jc
    integer                           :: ierr
    integer                           :: iA_start, ix_start, iy_start, index_offset
    integer                           :: n_blocksize, n_blocks
    real(kind=C_DOUBLE)               :: b_tmp_block(n_tor*n_var)

    b = 0.0

    do i=1,cooA%nnz
      ir = cooA%irn(i)
      jc = cooA%jcn(i)
      b(ir) = b(ir) + cooA%val(i) * x(jc)
    enddo

!    n_blocksize  = n_tor * n_var
!    n_blocks     = nz_glob/(n_blocksize*n_blocksize)
!    index_offset = (local_index_start(my_id+1)-1) * n_blocksize
!
!!$omp parallel default(none)                                      &
!!$omp shared(cooA, x, n_blocks, n_blocksize, index_offset)        &
!!$omp private(i,iA_start,ix_start, iy_start, b_tmp_block) &
!!$omp reduction(+:b)
!
!!$omp do schedule(guided)
!    do i = 1, n_blocks
!
!      iA_start = (i - 1)*n_blocksize*n_blocksize
!      ix_start = cooA%jcn(iA_start + 1)
!      iy_start = cooA%irn(iA_start + 1) - index_offset
!
!      call dgemv('T',n_blocksize,n_blocksize,1.d0,cooA%val(iA_start + 1),n_blocksize,x(ix_start),1,0.d0,b_tmp_block,1)
!
!      b(iy_start:iy_start + n_blocksize - 1) = b(iy_start:iy_start + n_blocksize - 1) + b_tmp_block(1:n_blocksize)
!
!    enddo
!!$omp end do
!!$omp end parallel

    call MPI_AllReduce(MPI_IN_PLACE,b,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_GLOB,ierr)


  end subroutine matv

!> apply preconditioner b = M\x
  subroutine pc(x,b)
    use preconditioner_module, only: my_row_index, my_row_factor
#ifdef USE_STRUMPACK
    use strumpack_module, only: strumpack_solve
#endif
    implicit none

    real(kind=C_DOUBLE), allocatable, target :: x(:), b(:)
    integer :: i
    integer :: ierr

    if (my_id_n.eq.0) then
      call MPI_BCAST(x,ndof_glob,MPI_DOUBLE_PRECISION,0,MPI_COMM_MASTER,ierr)
      do i = 1, mumps_par%n
        mumps_par%rhs(i) = x(my_row_index(i))
      enddo
    endif

#ifdef USE_STRUMPACK
    if (use_strumpack) then
      call MPI_BCAST(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
      call strumpack_solve(mumps_par%n,mumps_par%rhs,MPI_COMM_N)
    endif
#endif

    if (my_id_n.eq.0) then
      b = 0.0
      do i = 1, mumps_par%n
        b(my_row_index(i)) = mumps_par%rhs(i)*my_row_factor
      enddo
      call MPI_AllReduce(MPI_IN_PLACE,b,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_MASTER,ierr)
      ! now all my_id_n==0 ranks have the global solution vector
    endif
    call MPI_BCAST(b,ndof_glob,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
    ! now all ranks have the global solution vector

  end subroutine pc

#endif
end module mod_bicgstab
