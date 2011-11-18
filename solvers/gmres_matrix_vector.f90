subroutine gmres_matrix_vector(x,y,my_id,my_id_n, i_tor, MPI_COMM_MASTER)
!-----------------------------------------------------------------------
! sparse matrix vector product using coordinate scheme
! to be called only on all of MPI_COMM_WORLD
!
! result (y) is only known on id=0
!-----------------------------------------------------------------------
use tr_module 
use parameters
use global_distributed_matrix
use murge_module
implicit none
include 'mpif.h'

integer             :: i_tor(:), MPI_COMM_MASTER
real*8              :: x(:), y(:), t1, t2, t3, t4, t5
real*8, allocatable :: y_tmp(:), y_tmp2(:), x_tmp(:)
real*8              :: y_tmp_block(n_tor*n_var)
integer,allocatable :: recv_counts(:), recv_disp(:)
integer             :: n, i, ir, jc, ierr, my_id, my_id_n, n_cpu, index_offset
integer             :: n_blocksize, n_blocks, iA_start, ix_start, iy_start, ndof_local

integer index_ytmp_min,index_ytmp_max
logical found_value

!write(*,*) my_id,my_id_n,' GMRES matrix_vector ',ndof_glob
call cpu_time(t1)

call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)     ! the number of cpus

call cpu_time(t2)

call MPI_BCAST(x,ndof_glob,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

call cpu_time(t3)

y = 0
if ((use_murge) .and. (use_murge_element)) then
#ifdef USE_MURGE
   CALL MURGE_SetGlobalRHS(murge_id_prod, x, -1, MURGE_ASSEMBLY_OVW , ierr)
   CALL MURGE_GetGlobalProduct(murge_id_prod, y, 0, ierr)
#else
   print *, "Binary built without murge"
   call abort()
#endif
else
   n_blocksize = n_tor * n_var
   n_blocks    = nz_glob/n_blocksize**2
   ndof_local   = (local_index_end(my_id+1) - local_index_start(my_id+1) + 1) * n_blocksize
   index_offset = (local_index_start(my_id+1)-1) * n_blocksize 
   call tr_allocate(y_tmp,1,ndof_local,"y_tmp")
   y_tmp(1:ndof_local) = 0.d0
!$omp parallel default(none)                                                                                         &
!$omp   shared(y_tmp, A_glob, jcn_glob, irn_glob,x, n_blocks,n_blocksize, nz_glob, local_index_start, index_offset)  &
!$omp   private(i,iA_start,ix_start, iy_start, ir, jc, y_tmp_block)

!$omp do

   do i=1, n_blocks
      
      iA_start = (i-1) * n_blocksize**2
      ix_start = jcn_glob(iA_start+1)
      iy_start = irn_glob(iA_start+1) - index_offset
      
      call dgemv('T',n_blocksize,n_blocksize,1.d0,A_glob(iA_start+1),n_blocksize,x(ix_start),1,0.d0,y_tmp_block,1)
      
      !$omp critical
      y_tmp(iy_start:iy_start+n_blocksize-1) = y_tmp(iy_start:iy_start+n_blocksize-1) + y_tmp_block(1:n_blocksize)
      !$omp end critical
      
   end do
!$omp end do
!$omp end parallel

   !do i=1,nz_glob
   !  ir = irn_glob(i)
   !  jc = jcn_glob(i)
   !  y_tmp(ir) = y_tmp(ir) + A_glob(i) * x(jc)
   !enddo

   call cpu_time(t4)
   
   y(1:ndof_glob) = 0.d0
   
   call tr_allocate(recv_counts,1,n_cpu,"recv_counts")
   call tr_allocate(recv_disp,1,n_cpu,"recv_disp")

   do i=1,n_cpu
      recv_counts(i) = (local_index_end(i) - local_index_start(i) + 1) * n_blocksize
   enddo

   recv_disp(1) = 0
   do i=2,n_cpu
      recv_disp(i) = recv_disp(i-1) + recv_counts(i-1)
   enddo

   call mpi_gatherv(y_tmp,ndof_local,MPI_DOUBLE_PRECISION,y,recv_counts,recv_disp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call tr_deallocate(y_tmp,"y_tmp")
   call tr_deallocate(recv_counts,"recv_counts")
   call tr_deallocate(recv_disp,"recv_disp")
end if

call cpu_time(t5)


!write(*,'(A,i3,3f14.6)') ' M-V timing  barrier: ',my_id,t2-t1
!write(*,'(A,i3,3f14.6)') ' M-V timing  bcast  : ',my_id,t3-t2
!write(*,'(A,i3,3f14.6)') ' M-V timing  dgemv  : ',my_id,t4-t3
!write(*,'(A,i3,3f14.6)') ' M-V timing  reduce : ',my_id,t5-t4
!write(*,'(A,i3,3f14.6)') ' M-V timing  TOTAL  : ',my_id,t5-t1

return
end
