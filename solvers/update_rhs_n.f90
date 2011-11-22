subroutine update_rhs_n(my_id,my_id_n, i_tor, MPI_COMM_MASTER)
!***********************************************************************
!* subroutine updates the RHS with the explicit part of the matrix     *
!* multiplied with the previous solution (deltas)                      *
!***********************************************************************
use tr_module 
use parameters
use mumps_module
use global_distributed_matrix
implicit none
include 'mpif.h'

integer             :: i_tor(:)
integer             :: MPI_COMM_MASTER
integer             :: my_id, my_id_n, i, ir, jc
real*8,allocatable  :: rhs_delta(:), rhs_delta_n(:), deltas_n(:)

interface 
   subroutine gmres_matrix_vector(x,y,my_id,my_id_n, i_tor, MPI_COMM_MASTER)      
     integer             :: i_tor(:), MPI_COMM_MASTER
     real*8              :: x(:), y(:)
     integer             :: my_id, my_id_n
   end subroutine gmres_matrix_vector
end interface

!write(*,*) my_id,my_id_n,' starting update_rhs'

call tr_allocate(rhs_delta,1,ndof_glob,"rhs_delta",CAT_PRECOND)             ! size could be zero when my_id .gt.0 but needs to be allocated

rhs_delta(1:ndof_glob) = 0.d0

if (my_id_n .eq.0) then
  call tr_allocate(rhs_delta_n,1,mumps_par%n,"rhs_delta_n",CAT_PRECOND)
  call tr_allocate(deltas_n,1,mumps_par%n,"deltas_n",CAT_PRECOND)
  rhs_delta_n(1:mumps_par%n) = 0.d0
  deltas_n(1:mumps_par%n)    = 0.d0
endif

call gmres_matrix_vector(deltas,rhs_delta,my_id,my_id_n, i_tor, MPI_COMM_MASTER)

call distribute_vector(my_id,rhs_delta,rhs_delta_n)
call distribute_vector(my_id,deltas,deltas_n)

if (my_id_n .eq. 0) then              ! rhs lives only on masters

  do i=1,mumps_par%n
    mumps_par%rhs(i) = mumps_par%rhs(i) + rhs_delta_n(i)                  ! positive sign is on purpose (for stability reasons)
  enddo

  do i=1,mumps_par%nz
    ir = mumps_par%irn(i)
    jc = mumps_par%jcn(i)
    mumps_par%rhs(ir) = mumps_par%rhs(ir) - mumps_par%A(i) * deltas_n(jc) ! negative sign is on purpose (for stability reasons)
  enddo

endif

call tr_deallocate(rhs_delta,"rhs_delta",CAT_PRECOND)
if (my_id_n .eq.0) then
   call tr_deallocate(rhs_delta_n,"rhs_delta_n",CAT_PRECOND)
   call tr_deallocate(deltas_n,"deltas_n",CAT_PRECOND)
end if

return
end
