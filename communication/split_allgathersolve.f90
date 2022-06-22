subroutine split_allgathersolve(n_cpu,my_id,ad_mat,ac_mat)
!Split MPI_ALLGATHERV if MPI counts beyond 64-int

  use mpi_mod
  use mod_integer_types
  use data_structure, only: type_SP_MATRIX
  use tr_module

  implicit none

  ! --- Routine variables
  integer               :: n_cpu,my_id
  !integer(kind=int_all) :: counts(n_cpu),displacements(n_cpu)

  ! --- Local variables
  integer               :: i, i_cpu, ierr
  integer, allocatable  :: counts_int(:),displacements_int(:)
  integer(kind=int_all),allocatable :: counts(:), displacements(:) ! should be placed inside split_allgathersolve
  
  integer(kind=int_all) :: i_long

  integer(kind=int_all)  :: INT_MAX=1000000000 ! very conservative, could be up to ~2147000000
  logical                :: need_to_split
  integer                :: n_split, i_split
  integer(kind=int_all)  :: count_split

  integer(kind=int_all), parameter   :: Int1=1

  real*8,                allocatable :: Asend_buffer(:)
  integer(kind=int_all), allocatable :: isend_buffer(:), jsend_buffer(:)
  real*8,                allocatable :: Arecv_buffer(:)
  integer(kind=int_all), allocatable :: irecv_buffer(:), jrecv_buffer(:)
  integer(kind=int_all), allocatable :: index_buffer(:), index_target(:)
  
  type(type_SP_MATRIX)   :: ad_mat, ac_mat

  
  if (allocated(counts))        call tr_deallocate(counts,"counts",CAT_DMATRIX)
  if (allocated(displacements)) call tr_deallocate(displacements,"displacements",CAT_DMATRIX)

  call tr_allocate(counts,1,n_cpu,"counts",CAT_DMATRIX)
  call tr_allocate(displacements,1,n_cpu,"displacements",CAT_DMATRIX)  
  
  call MPI_Allgather(ad_mat%nnz,1,MPI_INTEGER_ALL,counts,1,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)

  displacements(1) = 0
  do i=2,n_cpu
    displacements(i) = displacements(i-1) + counts(i-1)
  enddo

  ! --- If we're using short integers, then this is just a simple wrapper, no need to split
#ifndef INTSIZE64
  call MPI_AllgatherV(ad_mat%irn,ad_mat%nnz,MPI_INTEGER_ALL,ac_mat%irn, &
                      counts,displacements,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)
  call MPI_AllgatherV(ad_mat%jcn,ad_mat%nnz,MPI_INTEGER_ALL,ac_mat%jcn, &
                      counts,displacements,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)
  call MPI_AllgatherV(ad_mat%val,ad_mat%nnz,MPI_DOUBLE_PRECISION,ac_mat%val, &
                      counts,displacements,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
                      
  call tr_deallocate(counts,"counts",CAT_DMATRIX)
  call tr_deallocate(displacements,"displacements",CAT_DMATRIX)
                      
  return
#else

  ! --- Otherwise, we might need to split the MPI communications

  ! --- Check if we need to split
  need_to_split = .false.
  if (ac_mat%nnz .gt. INT_MAX) need_to_split = .true.

  ! --- Allocate short-integer counts and displacements for MPI calls
  ! --- Counts still need to be copied because MPI count types are always short ints
  call tr_allocate(counts_int,1,n_cpu,"SPL_GATH_counts",CAT_DMATRIX)
  call tr_allocate(displacements_int,1,n_cpu,"SPL_GATH_displacements",CAT_DMATRIX)

  ! --- Split MPI calls
  if (need_to_split) then

    ! --- Allocate buffers where split matrix will be communicated before being copied into centralised matrix
    call tr_allocate(Arecv_buffer,Int1,INT_MAX,"SPL_GATH_Arecv_buffer",CAT_DMATRIX)
    call tr_allocate(irecv_buffer,Int1,INT_MAX,"SPL_GATH_irecv_buffer",CAT_DMATRIX)
    call tr_allocate(jrecv_buffer,Int1,INT_MAX,"SPL_GATH_jrecv_buffer",CAT_DMATRIX)
    call tr_allocate(Asend_buffer,Int1,INT_MAX,"SPL_GATH_Asend_buffer",CAT_DMATRIX)
    call tr_allocate(isend_buffer,Int1,INT_MAX,"SPL_GATH_isend_buffer",CAT_DMATRIX)
    call tr_allocate(jsend_buffer,Int1,INT_MAX,"SPL_GATH_jsend_buffer",CAT_DMATRIX)

    call tr_allocate(index_buffer,1,n_cpu,"SPL_GATH_index_buffer",CAT_DMATRIX)
    call tr_allocate(index_target,1,n_cpu,"SPL_GATH_index_target",CAT_DMATRIX)

    ! --- Split respective to the max send/recv
    n_split = ac_mat%nnz / INT_MAX + 1
    if (my_id .eq. 0) write(*,*) 'Warning: Splitting Main Matrix MPI centralisation',n_split

    !write(*,*)'******* BEFORE SPLIT:'
    !write(*,'(A,2i20)') 'nnz   ',my_id,mumps_par%nz
    !write(*,'(A,10i20)')'counts',my_id,counts
    !write(*,'(A,10i20)')'disps ',my_id,displacements
    do i_split=1,n_split
 
      ! --- Split counts for each MPI chunk
      do i_cpu=1,n_cpu
        count_split = counts(i_cpu) / n_split
        if (i_split .gt. 1) then
          index_buffer(i_cpu) = (i_split-1)*count_split
        else
          index_buffer(i_cpu) = 0
        endif
        counts_int(i_cpu) = count_split
        if (i_split .eq. n_split) then
          counts_int(i_cpu) = counts(i_cpu) - (n_split-1)*count_split
        endif
      enddo

      ! --- Split displacements for each MPI chunk
      displacements_int = 0
      index_target = 0
      count_split = counts(1) / n_split
      index_target(1) = displacements(1) + (i_split-1)*count_split
      do i_cpu=2,n_cpu
        displacements_int(i_cpu) = displacements_int(i_cpu-1) + counts_int(i_cpu-1)
        count_split = counts(i_cpu) / n_split
        index_target(i_cpu) = displacements(i_cpu) + (i_split-1)*count_split
      enddo

      ! --- Copy distributed matrix into send/recv buffers
      i_cpu = my_id+1
      isend_buffer(1:counts_int(i_cpu)) = ad_mat%irn(index_buffer(i_cpu)+1:index_buffer(i_cpu)+counts_int(i_cpu))
      jsend_buffer(1:counts_int(i_cpu)) = ad_mat%jcn(index_buffer(i_cpu)+1:index_buffer(i_cpu)+counts_int(i_cpu))
      Asend_buffer(1:counts_int(i_cpu)) = ad_mat%val(index_buffer(i_cpu)+1:index_buffer(i_cpu)+counts_int(i_cpu))

      ! --- Gather data on buffers
      call MPI_AllgatherV(isend_buffer,counts_int(my_id+1),MPI_INTEGER_ALL,irecv_buffer, &
                          counts_int,displacements_int,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(jsend_buffer,counts_int(my_id+1),MPI_INTEGER_ALL,jrecv_buffer, &
                          counts_int,displacements_int,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)
      call MPI_AllgatherV(Asend_buffer,counts_int(my_id+1),MPI_DOUBLE_PRECISION,Arecv_buffer, &
                          counts_int,displacements_int,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      ! --- Copy from buffer into centralised matrix
      do i_cpu=1,n_cpu
        ac_mat%val(index_target(i_cpu)+1:index_target(i_cpu)+counts_int(i_cpu)) &
          = Arecv_buffer(displacements_int(i_cpu)+1:displacements_int(i_cpu)+counts_int(i_cpu))
        ac_mat%irn(index_target(i_cpu)+1:index_target(i_cpu)+counts_int(i_cpu)) &
          = irecv_buffer(displacements_int(i_cpu)+1:displacements_int(i_cpu)+counts_int(i_cpu))
        ac_mat%jcn(index_target(i_cpu)+1:index_target(i_cpu)+counts_int(i_cpu)) &
          = jrecv_buffer(displacements_int(i_cpu)+1:displacements_int(i_cpu)+counts_int(i_cpu))
      enddo

    !write(*,*)'******* NEW SPLIT:',i_split,n_split
    !write(*,'(A,10i20)')'counts',my_id,counts_int
    !write(*,'(A,10i20)')'index ',my_id,index_buffer
    !write(*,'(A,10i20)')'disps ',my_id,displacements_int
    !write(*,'(A,10i20)')'target',my_id,index_target

    enddo
 
    ! --- Deallocate buffers of split matrix
    call tr_deallocate(Arecv_buffer,"SPL_GATH_Arecv_buffer",CAT_DMATRIX)
    call tr_deallocate(irecv_buffer,"SPL_GATH_irecv_buffer",CAT_DMATRIX)
    call tr_deallocate(jrecv_buffer,"SPL_GATH_jrecv_buffer",CAT_DMATRIX)
    call tr_deallocate(Asend_buffer,"SPL_GATH_Asend_buffer",CAT_DMATRIX)
    call tr_deallocate(isend_buffer,"SPL_GATH_isend_buffer",CAT_DMATRIX)
    call tr_deallocate(jsend_buffer,"SPL_GATH_jsend_buffer",CAT_DMATRIX)

    call tr_deallocate(index_buffer,"SPL_GATH_index_buffer",CAT_DMATRIX)
    call tr_deallocate(index_target,"SPL_GATH_index_target",CAT_DMATRIX)

  ! --- Don't split MPI calls
  else
    ! --- Counts still need to be copied because MPI count types are always short ints
    counts_int(:) = counts(:)
    displacements_int(:) = displacements(:)

    call MPI_AllgatherV(ad_mat%irn,ad_mat%nnz,MPI_INTEGER_ALL,ac_mat%irn, &
                        counts_int,displacements_int,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)
    call MPI_AllgatherV(ad_mat%jcn,ad_mat%nnz,MPI_INTEGER_ALL,ac_mat%jcn, &
                        counts_int,displacements_int,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)
    call MPI_AllgatherV(ad_mat%val,ad_mat%nnz,MPI_DOUBLE_PRECISION,ac_mat%val, &
                        counts_int,displacements_int,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  endif

  ! --- Deallocate short-integer counts and displacements
  call tr_deallocate(counts_int,"SPL_GATH_counts",CAT_DMATRIX)
  call tr_deallocate(displacements_int,"SPL_GATH_displacements",CAT_DMATRIX)
  call tr_deallocate(counts,"counts",CAT_DMATRIX)
  call tr_deallocate(displacements,"displacements",CAT_DMATRIX)
  
  return
#endif  

end subroutine split_allgathersolve
