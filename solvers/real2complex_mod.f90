!> converting real harmonic blocks into the complex ones
module real2complex_mod 

contains

#ifdef USE_COMPLEX_PRECOND
  !-- converting matrix entries into the complex form
  subroutine real2complex_a(my_id, my_id_n)
 
    use mod_parameters, only: n_var
    use mumps_module
    use global_distributed_matrix
    use mpi_mod 
   
    implicit none

    integer, intent(in) :: my_id, my_id_n
    integer             :: i, j, k, m

    !if (my_id .eq. 0) then
    !  write(*,*) my_id,'********************************************************'
    !  write(*,*) my_id,'* converting real harmonic block into the complex one  *'
    !  write(*,*) my_id,'********************************************************'
    !endif

    if(my_id_n .eq. 0) then 
      if (allocated(A_cmplx))  deallocate(A_cmplx) 
      if (allocated(irn_cmplx))deallocate(irn_cmplx) 
      if (allocated(jcn_cmplx))deallocate(jcn_cmplx) 

      if(my_id .eq. 0) then
        nz_cmplx = mumps_par%nz
      else
        nz_cmplx = mumps_par%nz/4
      endif

      allocate(A_cmplx(1:nz_cmplx))
      allocate(irn_cmplx(1:nz_cmplx))
      allocate(jcn_cmplx(1:nz_cmplx))

      if(my_id .eq. 0) then

        do i = 1, mumps_par%nz
          A_cmplx(i) = CMPLX(mumps_par%A(i)) 
          irn_cmplx(i) = mumps_par%irn(i)
          jcn_cmplx(i) = mumps_par%jcn(i)
        enddo 

      else !-- if my_id .ne. 0 

        do j = 1, nz_cmplx/n_var
          do i = 1, 2*n_var-1, 2
            m = i  + 4*(j-1)*n_var
            k = (i + 2*(j-1)*n_var + 1)/2
            irn_cmplx(k) = (mumps_par%irn(m)+1)/2    
            jcn_cmplx(k) = (mumps_par%jcn(m)+1)/2
            !-- Uncomment the line below to remove enforced symmetry 
            !A_cmplx(k) = CMPLX(mumps_par%A(m), -mumps_par%A(m+1)) 
            !-- Comment out the line below to remove enforced symmetry 
            A_cmplx(k) = CMPLX((mumps_par%A(m) + mumps_par%A(m+2*n_var+1))*0.5d0, -((mumps_par%A(m+1))-(mumps_par%A(m+2*n_var)))*0.5d0) 
          enddo
        enddo 
  
      endif !-- closing if(my_id .eq. 0)     
    endif !-- closing if(my_id_n .eq. 0)     
  
    return
  end subroutine real2complex_a

  !-- converting RHS into the complex form
  subroutine real2complex_rhs(my_id, my_id_n, rhs_cmplx)
 
    use mumps_module
    use global_distributed_matrix, only: n_cmplx 
    use mpi_mod 
   
    implicit none

    integer,                     intent(in)    :: my_id, my_id_n
    double complex, allocatable, intent(inout) :: rhs_cmplx(:) 
    integer                                    :: i

    !if (my_id .eq. 0) then
    !  write(*,*) my_id,'*****************************************'
    !  write(*,*) my_id,'* converting RHS into the complex form  *'
    !  write(*,*) my_id,'*****************************************'
    !endif

    if(my_id_n .eq. 0) then 
      if (allocated(rhs_cmplx))deallocate(rhs_cmplx)
      if(my_id .eq. 0) then
        n_cmplx = mumps_par%n
      else
        n_cmplx = mumps_par%n/2
      endif
      allocate(rhs_cmplx(1:n_cmplx))
  
      if(my_id .eq. 0) then
        do i = 1, mumps_par%n
          rhs_cmplx(i) = CMPLX(mumps_par%rhs(i))
        enddo  
      else !-- if my_id .ne. 0 
        do i = 1, mumps_par%n-1, 2
          rhs_cmplx((i+1)/2) = CMPLX(mumps_par%rhs(i),mumps_par%rhs(i+1))
        enddo  
      endif      
    endif      

    return
  end subroutine real2complex_rhs

#endif

end module real2complex_mod
