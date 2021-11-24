module mod_bicgstab
#ifdef USE_BICGSTAB
  use iso_c_binding

  implicit none

  type SPARSE_MATRIX_T
    integer(kind=C_INT), pointer :: irn(:), jcn(:)
    real(kind=C_DOUBLE), pointer :: val(:)
    integer :: indexing
    integer(kind=C_INT) :: n
    integer(kind=C_INT) :: nnz
  end type SPARSE_MATRIX_T
  
  type(SPARSE_MATRIX_T) :: cooA
 
  private
  public :: bicgstab_solve

  contains
  
  subroutine bicgstab_solve(irn, jcn, val, indexing, n, nnz, x, b, max_it, tol)
    implicit none
    
    integer(kind=C_INT), pointer :: irn(:), jcn(:)
    real(kind=C_DOUBLE), pointer :: val(:)    
    real(kind=C_DOUBLE), pointer :: x(:), b(:)
    real(kind=C_DOUBLE), pointer :: r(:), r_tld(:), s(:), s_hat(:), tmp(:), p(:), p_hat(:), v(:), t(:)
    integer(kind=C_INT) :: n, nnz, indexing
    integer :: max_it, iter, flag
    real(kind=C_DOUBLE) :: tol, error, alpha, beta, omega, bnrm2, rho, rho_1, resid, snrm2
    
    real(kind=C_DOUBLE), external :: dnrm2, ddot
    
    cooA%irn => irn
    cooA%jcn => jcn
    cooA%val => val
    cooA%indexing = indexing
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
    if (error < tol) return
    
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
     write(*,*) "bicgstab completed successfully with", iter, "iterations"
    elseif (omega == 0.0) then ! breakdown
     flag = -2
    elseif (rho == 0.0) then
     flag = -1
    else ! no convergence
     flag = 1
     write(*,*) "bicgstab failed to converge, error:", error
    endif
    
    write(*,*) "solution", x(1), x(n)

    deallocate(r, r_tld, s, s_hat, tmp, p, p_hat, v, t)
    
  end subroutine bicgstab_solve

  subroutine matv(x,b)
    implicit none

    real(kind=C_DOUBLE), pointer :: x(:), b(:)
    integer :: stt, i, j, ir, jc
    real(kind=c_double) alpha, beta, eps

    b = 0.0
    do i=1,cooA%nnz
      ir = cooA%irn(i)
      jc = cooA%jcn(i)
      b(ir) = b(ir) + cooA%val(i) * x(jc)
    enddo

  end subroutine matv
  
  subroutine pc(x,b)
    implicit none

    real(kind=C_DOUBLE), pointer :: x(:), b(:)
    
    b = x ! returns the same value
 
  end subroutine pc  
  
#endif
end module mod_bicgstab
