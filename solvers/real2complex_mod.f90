!> converting real harmonic blocks into the complex ones
module real2complex_mod 

contains

subroutine real2complex(my_id)
 
  use mod_parameters, only: n_tor, n_var
  use mumps_module
  use global_distributed_matrix
  use mpi_mod 
   
  implicit none

  integer, intent(in) :: my_id

  integer :: i, j, k, l, m

  if (my_id .eq. 0) then
    write(*,*) my_id,'*****************************************************'
    write(*,*) my_id,'* converting real harmonic block into the complex one  *'
    write(*,*) my_id,'*****************************************************'
  endif

  if (associated(A_cmplx))  deallocate(A_cmplx) 
  if (associated(rhs_cmplx))deallocate(rhs_cmplx)
  if (associated(irn_cmplx))deallocate(irn_cmplx) 
  if (associated(jcn_cmplx))deallocate(jcn_cmplx) 

  if(my_id .eq. 0) then 
    n_cmplx = mumps_par%n
    nz_cmplx = mumps_par%nz
  else
    n_cmplx = mumps_par%n/2
    nz_cmplx = mumps_par%nz/4
  endif

  allocate(A_cmplx(1:nz_cmplx))
  allocate(irn_cmplx(1:nz_cmplx))
  allocate(jcn_cmplx(1:nz_cmplx))
  allocate(rhs_cmplx(1:n_cmplx))

  if(my_id .eq. 0) then
    do i = 1, mumps_par%n
      rhs_cmplx(i) = CMPLX(mumps_par%rhs(i))
    enddo  
    do i = 1, mumps_par%nz
      A_cmplx(i) = CMPLX(mumps_par%A(i)) 
      irn_cmplx(i) = mumps_par%irn(i)
      jcn_cmplx(i) = mumps_par%jcn(i)
    enddo 
  else !-- if my_id .ne. 0 
    do i = 1, mumps_par%n, 2
      rhs_cmplx((i+1)/2) = CMPLX(mumps_par%rhs(i),mumps_par%rhs(i+1))
    enddo  
    l = mumps_par%nz/(4*n_var) 
    do j = 1, l
      do i = 1, 2*n_var, 2
        m = i  + 4*(j-1)*n_var
        k = (i + 2*(j-1)*n_var + 1)/2
        irn_cmplx(k) = (mumps_par%irn(m)+1)/2    
        jcn_cmplx(k) = (mumps_par%jcn(m)+1)/2
        A_cmplx(k) = CMPLX(mumps_par%A(m), -mumps_par%A(m+1))
      enddo
    enddo
  endif      

  return
end subroutine real2complex

end module real2complex_mod
