subroutine distribute_harmonics(my_id,my_id_n,n_cpu)
!---------------------------------------------------------------------
! extracts the reduced local matrices for each toroidal harmonic from
! the global matrix (distributed) :
!
!    A_glob(1:nz_glob), rhs_glob(1:ndof_glob)
!    irn_glob(1:nz_glob)
!    jcn_glob(1:nz_glob)
!
! sends the reduced local matrices to the masters only
!  (i.e. centralised matrices) :
!
!    mumps_par%A(1:mumps_par%nz), mumsp_par%rhs(1:mumps_par%n)
!    mumps_par%irn(1:mumps_par%nz)
!    mumps_par%jcn(1:mumps_par%nz)
!
!---------------------------------------------------------------------
use tr_module
use parameters
use global_distributed_matrix
use mumps_module

implicit none
include 'mpif.h'

integer :: my_id, my_id_n, n_cpu, m_cpu, i, j, idisp, nz_loc_n, n_loc_n, n_recv, in, ierr, ifactor
integer :: index(n_tor+1), index_snd, n_i, n_j, ibufsize
integer :: i_reduced, j_reduced

real*8,  allocatable :: Asnd_buffer(:), Rsnd_buffer(:)
integer, allocatable :: isnd_buffer(:), jsnd_buffer(:)
integer, allocatable :: send_counts(:), send_disp(:), recv_counts(:), recv_disp(:), sizes(:)

write(*,*) my_id,'*********************************'
write(*,*) my_id,'*      distributing matrix      *'
write(*,*) my_id,'*********************************'

ibufsize=0
do i=1,nz_glob                                    ! determine buffersize
  n_i = (mod(irn_glob(i)-1,n_tor) + 1) / 2        ! the toroidal modenumbers for this row-index
  n_j = (mod(jcn_glob(i)-1,n_tor) + 1) / 2        ! the toroidal modenumbers for this column-index
  if (n_i .eq. n_j) then                          ! select only the contributions from each toroidal harmonic
    ibufsize = ibufsize + 1
  endif
enddo
write(*,*) my_id,' ibufsize : ',ibufsize,nz_glob

call tr_allocate(Asnd_buffer,1,ibufsize,"dh_Asnd_buffer")
call tr_allocate(isnd_buffer,1,ibufsize,"dh_isnd_buffer")
call tr_allocate(jsnd_buffer,1,ibufsize,"dh_jsnd_buffer")
call tr_allocate(send_counts,1,n_cpu,"dh_send_counts")
call tr_allocate(send_disp,1,n_cpu,"dh_send_disp")
call tr_allocate(recv_counts,1,n_cpu,"dh_recv_counts")
call tr_allocate(recv_disp,1,n_cpu,"dh_recv_disp")

nz_loc_n = nz_glob    / n_tor**2
n_loc_n  = ndof_glob  / n_tor

M_cpu = n_cpu / ((n_tor+1)/2)

write(*,*) my_id,' nz_glob, nz_loc_n : ',nz_glob, nz_loc_n
write(*,*) my_id,' n, n_loc_n       : ',ndof_glob , n_loc_n
write(*,*) my_id,' n_cpu, M_cpu     : ',n_cpu, M_cpu


index(1) = 0
if (n_tor .gt. 1) index(2) = nz_loc_n

do i=3,(n_tor+1)/2
  index(i) = index(i-1) + 4*nz_loc_n                   ! offset for each harmonic in the send buffer (factor 4 because the number
enddo                                                  ! matrix elements of harmonic n is 4 times the size of the n=0 block

do i=1,nz_glob

  n_i = (mod(irn_glob(i)-1,n_tor) + 1) / 2        ! the toroidal modenumbers for this row-index
  n_j = (mod(jcn_glob(i)-1,n_tor) + 1) / 2        ! the toroidal modenumbers for this column-index

  if (n_i .eq. n_j) then                          ! select only the contributions from each toroidal harmonic

    index(n_i+1) = index(n_i+1) + 1

    Asnd_buffer(index(n_i+1)) = A_glob(i)
    isnd_buffer(index(n_i+1)) = irn_glob(i)
    jsnd_buffer(index(n_i+1)) = jcn_glob(i)

  endif

enddo

send_counts(1) = index(1)
send_disp(1)   = 0
idisp          = send_counts(1)

do j=2,n_cpu

  n_i = (j-1) / M_cpu

!  write(*,*) my_id,' j, n_i, index(ni+1) : ', j, n_i, index(n_i+1)

  if (mod(j-1,M_cpu) .eq. 0) then
    send_counts(j) = index(n_i+1) - index(n_i)
    send_disp(j)   = idisp
    idisp          = idisp + send_counts(j)
  else
    send_counts(j) = 0
    send_disp(j)   = 0
  endif

enddo

!write(*,'(i3,A,12i8)') my_id,' send_counts : ',send_counts
!write(*,'(i3,A,12i8)') my_id,' send_disp   : ',send_disp

call tr_allocate(sizes,1,n_cpu,"dh_sizes")

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

!write(*,*) my_id,' N_recv      : ',N_recv
!write(*,*) my_id,' ifactor     : ',ifactor
!write(*,'(i3,A,12i8)') my_id,' recv_counts : ',recv_counts
!write(*,'(i3,A,12i8)') my_id,' recv_disp   : ',recv_disp
!write(*,'(i3,A,12i8)') my_id,' sizes       : ',sizes

if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,"dh_mumps_par%A")
if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"dh_mumps_par%irn")
if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"dh_mumps_par%jcn")

call tr_allocatep(mumps_par%A,1,N_recv,"dh_mumps_par%A")
call tr_allocatep(mumps_par%irn,1,N_recv,"dh_mumps_par%irn")
call tr_allocatep(mumps_par%jcn,1,N_recv,"dh_mumps_par%jcn")

mumps_par%A = 0.d0
mumps_par%irn = 0
mumps_par%jcn = 0

call mpi_alltoallv(Asnd_buffer,send_counts,send_disp,MPI_DOUBLE_PRECISION, &
                   mumps_par%A,recv_counts,recv_disp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

call mpi_alltoallv(isnd_buffer,send_counts,send_disp,MPI_INTEGER, &
                   mumps_par%irn,recv_counts,recv_disp,MPI_INTEGER,MPI_COMM_WORLD,ierr)

call mpi_alltoallv(jsnd_buffer,send_counts,send_disp,MPI_INTEGER, &
                   mumps_par%jcn,recv_counts,recv_disp,MPI_INTEGER,MPI_COMM_WORLD,ierr)

!----------------------------------- change indices of the local matrices to local indices

if (my_id_n .eq. 0) then

  do i=1,mumps_par%nz

    n_i = (mod(mumps_par%irn(i)-1,n_tor) + 1) / 2
    n_j = (mod(mumps_par%jcn(i)-1,n_tor) + 1) / 2

    if (n_j .eq. 0) then
      j_reduced = (mumps_par%jcn(i)-1) / n_tor + 1
    else
      j_reduced = 2 * int((mumps_par%jcn(i)-1) / n_tor) + mod(mod(mumps_par%jcn(i)-1,n_tor)+1,2) + 1
    endif

    if (n_i .eq. 0) then
      i_reduced = (mumps_par%irn(i)-1) / n_tor + 1
    else
      i_reduced = 2 * int((mumps_par%irn(i)-1) / n_tor) + mod(mod(mumps_par%irn(i)-1,n_tor)+1,2) + 1
    endif

    mumps_par%irn(i) = i_reduced
    mumps_par%jcn(i) = j_reduced

  enddo

!  write(*,*) my_id,' distribute_harm check irn : ',minval(mumps_par%irn),maxval(mumps_par%irn)
!  write(*,*) my_id,' distribute_harm check jcn : ',minval(mumps_par%jcn),maxval(mumps_par%jcn)

endif


call tr_deallocate(Asnd_buffer,"dh_Asnd_buffer")
call tr_deallocate(isnd_buffer,"dh_isnd_buffer")
call tr_deallocate(jsnd_buffer,"dh_jsnd_buffer")
call tr_deallocate(send_counts,"dh_send_counts")
call tr_deallocate(send_disp  ,"dh_send_disp  ")
call tr_deallocate(recv_counts,"dh_recv_counts")
call tr_deallocate(recv_disp  ,"dh_recv_disp  ")
call tr_deallocate(sizes      ,"dh_sizes      ")

ifactor = 2
if (my_id .eq. 0)            ifactor = 1
if (mod(my_id,M_cpu) .ne. 0) ifactor = 0

mumps_par%n =  ifactor*n_loc_n
if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"dh_mumps_par%rhs")
call tr_allocatep(mumps_par%rhs,1,mumps_par%n,"dh_mumps_par%rhs")

call distribute_vector(my_id,rhs_glob,mumps_par%rhs)

return
end
