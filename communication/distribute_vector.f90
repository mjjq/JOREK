subroutine distribute_vector(my_id,rhs,rhs_dis)
!-----------------------------------------------------------------------
! distribute vector rhs to the MASTERS of each toroidal harmonic
!-----------------------------------------------------------------------
use tr_module 
use parameters
use global_distributed_matrix
use mumps_module

implicit none
include 'mpif.h'

real*8               :: rhs(*), rhs_dis(*)
integer              :: my_id, my_id_n, in, j, M_cpu, n_cpu, ifactor, idisp, n_loc_n, nz_loc_n, n_send, ierr, n_recv
integer              :: index_snd
real*8,  allocatable :: Asnd_buffer(:), Rsnd_buffer(:)
integer, allocatable :: send_counts(:), send_disp(:), recv_counts(:), recv_disp(:), sizes(:)

call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)

n_loc_n  = ndof_glob  / n_tor
M_cpu    = n_cpu / ((n_tor+1)/2)

call tr_allocate(send_counts,1,n_cpu    ,"dv_send_counts",CAT_DMATRIX)
call tr_allocate(send_disp,1,n_cpu      ,"dv_send_disp"  ,CAT_DMATRIX)
call tr_allocate(recv_counts,1,n_cpu    ,"dv_recv_counts",CAT_DMATRIX)
call tr_allocate(recv_disp,1,n_cpu      ,"dv_recv_disp"  ,CAT_DMATRIX)
call tr_allocate(Rsnd_buffer,1,ndof_glob,"dv_Rsnd_buffer",CAT_DMATRIX)

if (my_id .eq. 0) then

  Rsnd_buffer(1:n_loc_n) = rhs(1:ndof_glob:n_tor)

  do in=2, (n_tor+1)/2

    index_snd = n_loc_n + (in-2)*2*n_loc_n

    Rsnd_buffer(index_snd+1:index_snd+2*n_loc_n:2) = rhs(2*in-2:ndof_glob:n_tor)
    Rsnd_buffer(index_snd+2:index_snd+2*n_loc_n:2) = rhs(2*in-1:ndof_glob:n_tor)

  enddo

  send_counts(1) = n_loc_n
  send_disp(1)   = 0
  idisp          = send_counts(1)

  do j=2,n_cpu

    if (mod(j-1,M_cpu) .eq. 0) then
      send_counts(j) = 2*n_loc_n
      send_disp(j)   = idisp
      idisp          = idisp + send_counts(j)
    else
      send_counts(j) = 0
      send_disp(j)   = 0
    endif

  enddo

endif

ifactor = 2
if (my_id .eq. 0)            ifactor = 1
if (mod(my_id,M_cpu) .ne. 0) ifactor = 0

n_send =  ifactor*n_loc_n
mumps_par%n =  ifactor*n_loc_n

nz_loc_n = nz_glob    / n_tor**2
call tr_allocate(sizes,1,n_cpu,"dh_sizes",CAT_DMATRIX)

call mpi_allgather(nz_loc_n,1,MPI_INTEGER,sizes,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

recv_counts = 0
recv_disp   = 0

ifactor = 4                                   ! WHY 4, I have forgotten
if (my_id .eq. 0)            ifactor = 1
if (mod(my_id,M_cpu) .ne. 0) ifactor = 0

if (mod(my_id,M_cpu) .eq. 0) then
  do j=1,N_cpu
    recv_counts(j) = ifactor*sizes(j)
  enddo
  recv_disp(1)   = 0
  do j=2,N_cpu
    recv_disp(j) = recv_disp(j-1) + ifactor*sizes(j-1)
  enddo
endif

N_recv = sum(recv_counts)

mumps_par%nz = N_recv


call mpi_scatterv(Rsnd_buffer,send_counts,send_disp,MPI_DOUBLE_PRECISION, &
                  rhs_dis,n_send,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

call tr_deallocate(Rsnd_buffer  ,"dv_Rsnd_buffer" ,CAT_DMATRIX)
call tr_deallocate(send_counts  ,"dv_send_counts" ,CAT_DMATRIX)
call tr_deallocate(send_disp    ,"dv_send_disp"   ,CAT_DMATRIX)
call tr_deallocate(recv_counts  ,"dv_recv_counts" ,CAT_DMATRIX)
call tr_deallocate(recv_disp    ,"dv_recv_disp"   ,CAT_DMATRIX)
call tr_deallocate(sizes        ,"dh_sizes"       ,CAT_DMATRIX)
return
end
