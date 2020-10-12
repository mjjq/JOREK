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
use mod_parameters
use global_distributed_matrix
use mumps_module
use mpi_mod
use mod_integer_types

implicit none

interface 
   subroutine distribute_vector(my_id,rhs,rhs_dis,again)
     real*8               :: rhs(:), rhs_dis(:)
     integer              :: my_id
     logical              :: again
   end subroutine distribute_vector
end interface

integer                :: my_id, my_id_n, n_cpu, m_cpu, idisp, in, ierr, ifactor
integer(kind=int_all)  :: Int1=1, i, j, nz_loc_n, n_loc_n, nrecv_max, nsend_max
integer(kind=int_all)  :: index((n_tor+1)/2), index_snd, n_i, n_j, ibufsize
integer(kind=int_all)  :: i_reduced, j_reduced

real*8,                allocatable :: Asnd_buffer(:)
integer(kind=int_all), allocatable :: isnd_buffer(:), jsnd_buffer(:)
integer(kind=int_all), allocatable :: sizes_nz(:), sizes_buff(:)
integer(kind=int_all), allocatable :: long_recv_counts(:), long_recv_disp(:)
integer(kind=int_all)              :: mod_frac, mod_arg_i, mod_arg_j
integer,               allocatable :: send_counts(:), send_disp(:), recv_counts(:), recv_disp(:)

integer(kind=int_all)  :: INT_MAX
logical                :: need_to_split
integer                :: n_split, i_split, count_all

real*8,                allocatable :: Arecv_buffer(:)
integer(kind=int_all), allocatable :: irecv_buffer(:), jrecv_buffer(:), n_recv_prev(:)
integer,               allocatable :: index_split_local(:,:), index_split(:,:,:)

if (my_id .eq. 0) then
  write(*,*) my_id,'*********************************'
  write(*,*) my_id,'*      distributing matrix      *'
  write(*,*) my_id,'*********************************'
endif

! --- integer limit for short integers (normally 2147483647)
#ifdef INTSIZE64
  ! --- Not sure why, but it seems MPI fails even with counters below the long-int limit
  ! --- Maybe MPI has some internal working arrays that need to be larger than the counters? half seems to work well...
  INT_MAX = 1000000000 !1000000000
#else
  ! --- If we're not using long-ints, then there is nothing to split anyway
  INT_MAX = 2147000000
#endif


! --- Allocate MPI send/recv counters
call tr_allocate(send_counts     ,1,n_cpu,"dh_send_counts"     ,CAT_DMATRIX)
call tr_allocate(send_disp       ,1,n_cpu,"dh_send_disp"       ,CAT_DMATRIX)
call tr_allocate(recv_counts     ,1,n_cpu,"dh_recv_counts"     ,CAT_DMATRIX)
call tr_allocate(recv_disp       ,1,n_cpu,"dh_recv_disp"       ,CAT_DMATRIX)
call tr_allocate(long_recv_counts,1,n_cpu,"dh_long_recv_counts",CAT_DMATRIX)
call tr_allocate(long_recv_disp  ,1,n_cpu,"dh_long_recv_disp"  ,CAT_DMATRIX)
call tr_allocate(n_recv_prev     ,1,n_cpu,"dh_n_recv_prev"     ,CAT_DMATRIX)

! --- Get the size of harmonics matrix on this process
ibufsize=0
do i=1,nz_glob                                    ! determine buffersize
  ! --- just to keep safe, because fortran modulo is a short-int function...
  mod_frac  = ( irn_glob(i)-1 ) / n_tor
  mod_arg_i = ( irn_glob(i)-1 ) - mod_frac * n_tor
  mod_frac  = ( jcn_glob(i)-1 ) / n_tor
  mod_arg_j = ( jcn_glob(i)-1 ) - mod_frac * n_tor
  n_i = (mod(mod_arg_i,n_tor) + 1) / 2            ! the toroidal modenumbers for this row-index
  n_j = (mod(mod_arg_j,n_tor) + 1) / 2            ! the toroidal modenumbers for this column-index
  if (n_i .eq. n_j) then                          ! select only the contributions from each toroidal harmonic
    ibufsize = ibufsize + 1
  endif
enddo

! --- Get the size of the matrix corresponding to a single i_tor block (note, a non-zero harmonic has two blocks)
nz_loc_n = nz_glob    / n_tor**2
n_loc_n  = ndof_glob  / n_tor

! --- Number of cpu's per harmonic
M_cpu = n_cpu / ((n_tor+1)/2)

! --- The size of the harmonic matrices on each process
call tr_allocate(sizes_nz  ,1,n_cpu,"dh_sizes_nz"  ,CAT_DMATRIX)
call tr_allocate(sizes_buff,1,n_cpu,"dh_sizes_buff",CAT_DMATRIX)
call mpi_allgather(nz_loc_n,1,MPI_INTEGER_ALL,sizes_nz  ,1,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)
call mpi_allgather(ibufsize,1,MPI_INTEGER_ALL,sizes_buff,1,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)

! --- The starting index of each harmonic block in the full Aglob (without splits)
index(1) = 0
if (n_tor .gt. 1) then
  index(2) = nz_loc_n
  do i=3,(n_tor+1)/2
    index(i) = index(i-1) + 4*nz_loc_n  ! offset for each harmonic in the send buffer (factor 4 because the number
  enddo                                 ! matrix elements of harmonic n is 4 times the size of the n=0 block
endif

! --- Figure out the final size of the harmonic matrix this process is responsible for
long_recv_counts = 0
if (my_id .eq. 0) then
  ifactor = 1 ! the first set  of cpu's is  for n=0 (only one  block  of nz_loc_n)
else
  ifactor = 4 ! the other sets of cpu's are for n>0 (each four blocks of nz_loc_n)
endif
! --- We receive the matrix only on the main process of each harmonic
if (mod(my_id,M_cpu) .eq. 0) then
  do j=1,N_cpu
    long_recv_counts(j) = ifactor*sizes_nz(j)
  enddo
endif
mumps_par%nz = sum(long_recv_counts)

! --- Allocate harmonic matrix (one per process, for each harmonics, each process on that harmonic has a copy of the matrix)
if (associated(mumps_par%A  )) call tr_deallocatep(mumps_par%A  ,"dh_mumps_par%A"  ,CAT_DMATRIX)
if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"dh_mumps_par%irn",CAT_DMATRIX)
if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"dh_mumps_par%jcn",CAT_DMATRIX)
call tr_allocatep(mumps_par%A  ,Int1,mumps_par%nz,"dh_mumps_par%A"  ,CAT_DMATRIX)
call tr_allocatep(mumps_par%irn,Int1,mumps_par%nz,"dh_mumps_par%irn",CAT_DMATRIX)
call tr_allocatep(mumps_par%jcn,Int1,mumps_par%nz,"dh_mumps_par%jcn",CAT_DMATRIX)
mumps_par%A = 0.d0
mumps_par%irn = 0
mumps_par%jcn = 0

! --- Figure out the size of the largest harmonic matrix received on any of all processes
long_recv_counts = 0
ifactor = 4 ! for n>0 each matrix has four blocks of nz_loc_n
do j=1,N_cpu
  long_recv_counts(j) = ifactor*sizes_nz(j)
enddo
nrecv_max = sum(long_recv_counts)

! --- Figure out the size of the largest sent data on any of all processes
nsend_max = maxval(sizes_buff)

! --- Make sure we are below the short-int limit (necessary for MPI communication counters)
need_to_split = .false.
if (max(nrecv_max,nsend_max) .gt. INT_MAX) need_to_split = .true.


! -----------------------
! --- START MPI-SPLIT ---
! -----------------------
! --- Splitting the matrix, and doing mpi_alltoallv in chunks
if (need_to_split) then

  ! --- Allocate buffers where distributed matrix will be copied before being sent to corresponding harmonic processes
  call tr_allocate(Asnd_buffer,Int1,INT_MAX,"dh_Asnd_buffer",CAT_DMATRIX)
  call tr_allocate(isnd_buffer,Int1,INT_MAX,"dh_isnd_buffer",CAT_DMATRIX)
  call tr_allocate(jsnd_buffer,Int1,INT_MAX,"dh_jsnd_buffer",CAT_DMATRIX)

  ! --- Allocate buffers where harmonics matrix will be sent to before being copied into each harmonic matrix
  call tr_allocate(Arecv_buffer,Int1,INT_MAX,"dh_Arecv_buffer",CAT_DMATRIX)
  call tr_allocate(irecv_buffer,Int1,INT_MAX,"dh_irecv_buffer",CAT_DMATRIX)
  call tr_allocate(jrecv_buffer,Int1,INT_MAX,"dh_jrecv_buffer",CAT_DMATRIX)

  ! --- Record displacement in main matrix to copy buffers into correct location
  long_recv_disp   = 0
  ! --- We receive the matrix only on the main process of each harmonic
  if (mod(my_id,M_cpu) .eq. 0) then
    if (my_id .eq. 0) then
      ifactor = 1 ! the first set  of cpu's is  for n=0 (only one  block  of nz_loc_n)
    else
      ifactor = 4 ! the other sets of cpu's are for n>0 (each four blocks of nz_loc_n)
    endif
    long_recv_disp(1)   = 0
    do j=2,N_cpu
      long_recv_disp(j) = long_recv_disp(j-1) + ifactor*sizes_nz(j-1)
    enddo
  endif

  ! --- Split respective to the max send/recv
  if (nrecv_max .gt. nsend_max) then
    ! --- If it's because of what we send, easy, just split what you send.
    ! --- However, if it's because of what we receive, then we assume that the matrix 
    ! --- is ~evenly spread out between all cpu's, so we derive a new INT_MAX as:
    INT_MAX = real(INT_MAX) / real(n_cpu)
    ! --- But, note that each sender has a block made of (1 + 4(n_tor-1)/2)*nz_loc_n
    ! --- sub-blocks, ie. (1/4 + (n_tor-1)/2) blocks of (*4*nz_loc_n)
    INT_MAX = real(INT_MAX) * (real((n_tor-1)/2) + 0.25)
    ! --- To keep safe, just take 80% of all that, to account for the uneven distribution
    ! --- of the matrix among all cpu's
    INT_MAX = 0.80 * real(INT_MAX)
  endif
  n_split = nsend_max / INT_MAX + 1

  ! --- Before starting, get indexing for each cpu, for each harmonic, for each split
  call tr_allocate(index_split_local,        1,n_split,1,(n_tor+1)/2+1,"dh_index_split_local",CAT_DMATRIX)
  call tr_allocate(index_split      ,1,n_cpu,1,n_split,1,(n_tor+1)/2+1,"dh_index_split"      ,CAT_DMATRIX)

  index_split_local = 0
  do i_split=1,n_split

    ! --- How do we split the matrix? We need new indices for each harmonic
    ibufsize  = 0
    count_all = 0
    do i=1,nz_glob                                    ! determine buffersize
      ! --- just to keep safe, because fortran modulo is a short-int function...
      mod_frac  = ( irn_glob(i)-1 ) / n_tor
      mod_arg_i = ( irn_glob(i)-1 ) - mod_frac * n_tor
      mod_frac  = ( jcn_glob(i)-1 ) / n_tor
      mod_arg_j = ( jcn_glob(i)-1 ) - mod_frac * n_tor
      n_i = (mod(mod_arg_i,n_tor) + 1) / 2            ! the toroidal modenumbers for this row-index
      n_j = (mod(mod_arg_j,n_tor) + 1) / 2            ! the toroidal modenumbers for this column-index
      if (n_i .eq. n_j) then                          ! select only the contributions from each toroidal harmonic
        count_all = count_all + 1
        ! --- Start counting at the relevant split
        if (count_all .gt. (i_split-1)*INT_MAX) then
          index_split_local(i_split,n_i+1) = index_split_local(i_split,n_i+1) + 1
          ibufsize = ibufsize + 1
          if (ibufsize .ge. INT_MAX) exit
        endif
      endif
    enddo
    ! --- The start of each harmonic block is at the end of the previous block, obviously...
    do i=(n_tor+1)/2+1,2,-1
      index_split_local(i_split,i) = 0
      do j=1,i-1
        index_split_local(i_split,i) = index_split_local(i_split,i) + index_split_local(i_split,j)
      enddo
    enddo
    index_split_local(i_split,1) = 0

  enddo

  ! --- Gather indices from all cpu's
  do i_split=1,n_split
    do i=1,(n_tor+1)/2+1
      call mpi_allgather(index_split_local(i_split,i),1,MPI_INTEGER,index_split(:,i_split,i),1,MPI_INTEGER ,MPI_COMM_WORLD,ierr)
    enddo
  enddo

  ! --- Loop over each mpi-chunk
  n_recv_prev = 0
  do i_split=1,n_split

    ! --- We already know the index of those we send
    do i=1,(n_tor+1)/2
      index(i) = index_split_local(i_split,i)
    enddo

    ! --- Copy matrix into send-buffer
    ibufsize  = 0
    count_all = 0
    do i=1,nz_glob
      ! --- just to keep safe, because fortran modulo is a short-int function...
      mod_frac  = ( irn_glob(i)-1 ) / n_tor
      mod_arg_i = ( irn_glob(i)-1 ) - mod_frac * n_tor
      mod_frac  = ( jcn_glob(i)-1 ) / n_tor
      mod_arg_j = ( jcn_glob(i)-1 ) - mod_frac * n_tor
      n_i = (mod(mod_arg_i,n_tor) + 1) / 2            ! the toroidal modenumbers for this row-index
      n_j = (mod(mod_arg_j,n_tor) + 1) / 2            ! the toroidal modenumbers for this column-index
      if (n_i .eq. n_j) then                          ! select only the contributions from each toroidal harmonic
        count_all = count_all + 1
        ! --- Start counting at the relevant split
        if (count_all .gt. (i_split-1)*INT_MAX) then
          index(n_i+1) = index(n_i+1) + 1             ! index reference for each harmonic block
          Asnd_buffer(index(n_i+1)) = A_glob(i)
          isnd_buffer(index(n_i+1)) = irn_glob(i)
          jsnd_buffer(index(n_i+1)) = jcn_glob(i)
          ibufsize = ibufsize + 1
          if (ibufsize .ge. INT_MAX) exit
        endif
      endif
    enddo

    ! --- MPI send-counters for the location of each block
    send_counts(1) = index(1)
    send_disp(1)   = 0
    idisp          = send_counts(1)
    do j=2,n_cpu
      n_i = (j-1) / M_cpu
      ! --- We send the matrix only to the main process of each harmonic
      if (mod(j-1,M_cpu) .eq. 0) then
        send_counts(j) = index(n_i+1) - index(n_i)
        send_disp(j)   = idisp
        idisp          = idisp + send_counts(j)
      else
        send_counts(j) = 0
        send_disp(j)   = 0
      endif
    enddo

    ! --- MPI receive-counters for the location of each block
    recv_counts = 0
    recv_disp   = 0
    ! --- We receive the matrix only on the main process of each harmonic
    if (mod(my_id,M_cpu) .eq. 0) then
      i = my_id / M_cpu + 1 ! the corresponding toroidal n
      do j=1,N_cpu
        recv_counts(j) = index_split(j,i_split,i+1) - index_split(j,i_split,i)
      enddo
      recv_disp(1)   = 0
      do j=2,N_cpu
        recv_disp(j) = recv_disp(j-1) + recv_counts(j-1)
      enddo
    endif

    ! --- Send the data to the buffers
    call mpi_alltoallv(Asnd_buffer ,send_counts,send_disp,MPI_DOUBLE_PRECISION, &
                       Arecv_buffer,recv_counts,recv_disp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

    call mpi_alltoallv(isnd_buffer ,send_counts,send_disp,MPI_INTEGER_ALL, &
                       irecv_buffer,recv_counts,recv_disp,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)

    call mpi_alltoallv(jsnd_buffer ,send_counts,send_disp,MPI_INTEGER_ALL, &
                       jrecv_buffer,recv_counts,recv_disp,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)

    ! --- Copy from buffer into local harmonic matrix
    do j=1,n_cpu
      do i=1,recv_counts(j)
        mumps_par%A  (long_recv_disp(j)+n_recv_prev(j)+i) = Arecv_buffer(recv_disp(j)+i)
        mumps_par%irn(long_recv_disp(j)+n_recv_prev(j)+i) = irecv_buffer(recv_disp(j)+i)
        mumps_par%jcn(long_recv_disp(j)+n_recv_prev(j)+i) = jrecv_buffer(recv_disp(j)+i)
      enddo
      n_recv_prev(j) = n_recv_prev(j) + recv_counts(j)
    enddo

  enddo ! n_split

  ! --- Deallocate split buffers
  call tr_deallocate(Arecv_buffer,"dh_Arecv_buffer",CAT_DMATRIX)
  call tr_deallocate(irecv_buffer,"dh_irecv_buffer",CAT_DMATRIX)
  call tr_deallocate(jrecv_buffer,"dh_jrecv_buffer",CAT_DMATRIX)

  ! --- Deallocate split indices
  call tr_deallocate(index_split      ,"dh_index_split"      ,CAT_DMATRIX)
  call tr_deallocate(index_split_local,"dh_index_split_local",CAT_DMATRIX)

! --- Not splitting the matrix, can just mpi_alltoallv with everything right away
else

  ! --- Allocate buffers where distributed matrix will be copied before being sent to corresponding harmonic processes
  call tr_allocate(Asnd_buffer,Int1,ibufsize,"dh_Asnd_buffer",CAT_DMATRIX)
  call tr_allocate(isnd_buffer,Int1,ibufsize,"dh_isnd_buffer",CAT_DMATRIX)
  call tr_allocate(jsnd_buffer,Int1,ibufsize,"dh_jsnd_buffer",CAT_DMATRIX)

  ! --- Copy matrix into send-buffer
  do i=1,nz_glob
    n_i = (mod(irn_glob(i)-1,n_tor) + 1) / 2        ! the toroidal modenumbers for this row-index
    n_j = (mod(jcn_glob(i)-1,n_tor) + 1) / 2        ! the toroidal modenumbers for this column-index
    if (n_i .eq. n_j) then                          ! select only the contributions from each toroidal harmonic
      index(n_i+1) = index(n_i+1) + 1               ! index reference for each harmonic block
      Asnd_buffer(index(n_i+1)) = A_glob(i)
      isnd_buffer(index(n_i+1)) = irn_glob(i)
      jsnd_buffer(index(n_i+1)) = jcn_glob(i)
    endif
  enddo

  ! --- MPI send-counters for the location of each block
  send_counts(1) = index(1)
  send_disp(1)   = 0
  idisp          = send_counts(1)
  do j=2,n_cpu
    n_i = (j-1) / M_cpu
    ! --- We send the matrix only to the main process of each harmonic
    if (mod(j-1,M_cpu) .eq. 0) then
      send_counts(j) = index(n_i+1) - index(n_i)
      send_disp(j)   = idisp
      idisp          = idisp + send_counts(j)
    else
      send_counts(j) = 0
      send_disp(j)   = 0
    endif
  enddo

  ! --- MPI receive-counters for the location of each block
  recv_counts = 0
  recv_disp   = 0
  ! --- We receive the matrix only on the main process of each harmonic
  if (mod(my_id,M_cpu) .eq. 0) then
    if (my_id .eq. 0) then
      ifactor = 1 ! the first set  of cpu's is  for n=0 (only one  block  of nz_loc_n)
    else
      ifactor = 4 ! the other sets of cpu's are for n>0 (each four blocks of nz_loc_n)
    endif
    do j=1,N_cpu
      recv_counts(j) = ifactor*sizes_nz(j)
    enddo
    recv_disp(1)   = 0
    do j=2,N_cpu
      recv_disp(j) = recv_disp(j-1) + ifactor*sizes_nz(j-1)
    enddo
  endif

  call mpi_alltoallv(Asnd_buffer,send_counts,send_disp,MPI_DOUBLE_PRECISION, &
                     mumps_par%A,recv_counts,recv_disp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call mpi_alltoallv(isnd_buffer,send_counts,send_disp,MPI_INTEGER_ALL, &
                     mumps_par%irn,recv_counts,recv_disp,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)

  call mpi_alltoallv(jsnd_buffer,send_counts,send_disp,MPI_INTEGER_ALL, &
                     mumps_par%jcn,recv_counts,recv_disp,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)

endif ! split





! --- Deallocate global buffers
call tr_deallocate(Asnd_buffer,"dh_Asnd_buffer",CAT_DMATRIX)
call tr_deallocate(isnd_buffer,"dh_isnd_buffer",CAT_DMATRIX)
call tr_deallocate(jsnd_buffer,"dh_jsnd_buffer",CAT_DMATRIX)

! --- Deallocate MPI send/recv counters
call tr_deallocate(send_counts     ,"dh_send_counts"     ,CAT_DMATRIX)
call tr_deallocate(send_disp       ,"dh_send_disp"       ,CAT_DMATRIX)
call tr_deallocate(recv_counts     ,"dh_recv_counts"     ,CAT_DMATRIX)
call tr_deallocate(recv_disp       ,"dh_recv_disp"       ,CAT_DMATRIX)
call tr_deallocate(long_recv_counts,"dh_long_recv_counts",CAT_DMATRIX)
call tr_deallocate(long_recv_disp  ,"dh_long_recv_disp"  ,CAT_DMATRIX)
call tr_deallocate(n_recv_prev     ,"dh_n_recv_prev"     ,CAT_DMATRIX)
call tr_deallocate(sizes_nz        ,"dh_sizes_nz"        ,CAT_DMATRIX)
call tr_deallocate(sizes_buff      ,"dh_sizes_buff"      ,CAT_DMATRIX)







! --- change indices of the local matrices to local indices
if (my_id_n .eq. 0) then

  do i=1,mumps_par%nz

    ! --- just to keep safe, because fortran modulo is a short-int function...
    mod_frac  = ( mumps_par%irn(i)-1 ) / n_tor
    mod_arg_i = ( mumps_par%irn(i)-1 ) - mod_frac * n_tor
    mod_frac  = ( mumps_par%jcn(i)-1 ) / n_tor
    mod_arg_j = ( mumps_par%jcn(i)-1 ) - mod_frac * n_tor
    n_i = (mod(mod_arg_i,n_tor) + 1) / 2
    n_j = (mod(mod_arg_j,n_tor) + 1) / 2

    if (n_j .eq. 0) then
      j_reduced = (mumps_par%jcn(i)-1) / n_tor + 1
    else
      ! --- just to keep safe, because fortran modulo is a short-int function...
      mod_frac  = ( mumps_par%jcn(i)-1 ) / n_tor
      mod_arg_j = ( mumps_par%jcn(i)-1 ) - mod_frac * n_tor
#ifdef INTSIZE64
      j_reduced = 2 * int8((mumps_par%jcn(i)-1) / n_tor) + mod(mod(mumps_par%jcn(i)-1,n_tor)+1,2) + 1
#else
      j_reduced = 2 * int((mumps_par%jcn(i)-1) / n_tor) + mod(mod(mumps_par%jcn(i)-1,n_tor)+1,2) + 1
#endif
    endif

    if (n_i .eq. 0) then
      i_reduced = (mumps_par%irn(i)-1) / n_tor + 1
    else
      ! --- just to keep safe, because fortran modulo is a short-int function...
      mod_frac  = ( mumps_par%irn(i)-1 ) / n_tor
      mod_arg_i = ( mumps_par%irn(i)-1 ) - mod_frac * n_tor
#ifdef INTSIZE64
      i_reduced = 2 * int8((mumps_par%irn(i)-1) / n_tor) + mod(mod(mumps_par%irn(i)-1,n_tor)+1,2) + 1
#else
      i_reduced = 2 * int((mumps_par%irn(i)-1) / n_tor) + mod(mod(mumps_par%irn(i)-1,n_tor)+1,2) + 1
#endif
    endif

    mumps_par%irn(i) = i_reduced
    mumps_par%jcn(i) = j_reduced

  enddo

endif







! --- Now do the same thing for the RHS matrix (done inside distribute_vector)
if (my_id .lt. M_cpu) then
   ifactor = 1 ! n=0 has only one block
else 
   ifactor = 2 ! n>0 has sine and cosine
end if

mumps_par%n = ifactor*n_loc_n
if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"dh_mumps_par%rhs",CAT_DMATRIX)
call tr_allocatep(mumps_par%rhs,Int1,mumps_par%n,"dh_mumps_par%rhs",CAT_DMATRIX)

call distribute_vector(my_id,rhs_glob,mumps_par%rhs,.false.)






return
end
