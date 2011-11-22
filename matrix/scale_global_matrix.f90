subroutine scale_global_matrix
!***********************************************************************
!* column scaling of the global matrix                                 *
!***********************************************************************
use tr_module 
use global_distributed_matrix
implicit none
include 'mpif.h'
real*8, allocatable :: column_local(:)
integer :: k,j,ierr

if (allocated(column_scaling))  call tr_deallocate(column_scaling,"column_scaling",CAT_DMATRIX)
if (allocated(column_local))    call tr_deallocate(column_local,"column_local",CAT_DMATRIX)
call tr_allocate(column_scaling,1,ndof_glob,"column_scaling",CAT_DMATRIX)
call tr_allocate(column_local,1,ndof_glob,"column_local",CAT_DMATRIX)

column_local = 1.d-20;   column_scaling = 1.d-20

do k=1,nz_glob
  j = jcn_glob(k)
  column_local(j) = max(column_local(j),abs(A_glob(k)))
enddo

call MPI_AllReduce(column_local,column_scaling,ndof_glob,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

do k=1,nz_glob
  j = jcn_glob(k)
  A_glob(k) = A_glob(k) / column_scaling(j)
enddo

call tr_deallocate(column_local,"column_local",CAT_DMATRIX)

return
end
