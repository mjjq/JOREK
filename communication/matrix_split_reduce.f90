!> Reduces distributed matrix ad_mat into centralized matrix ac_mat
! using split MPI_Allgatherv if MPI counts beyond INT_MAX
subroutine matrix_split_reduce(ad_mat, ac_mat)

  use mpi_mod
  use mod_integer_types
  use data_structure, only: type_SP_MATRIX
  use tr_module
  implicit none

  type(type_SP_MATRIX), intent(inout)  :: ad_mat
  type(type_SP_MATRIX), intent(inout)  :: ac_mat


  integer               :: my_id, n_cpu, comm
  integer               :: i, i_cpu, ierr
  integer(kind=int_all) :: is, ie
  integer, allocatable  :: counts_int(:),displacements_int(:)
  integer(kind=int_all),allocatable :: counts(:), displacements(:) ! should be placed inside split_allgathersolve

  integer(kind=int_all) :: i_long

  logical                :: need_to_split
  integer                :: n_split, i_split
  integer(kind=int_all)  :: count_split

  real*8,                allocatable :: Asend_buffer(:)
  integer(kind=int_all), allocatable :: isend_buffer(:), jsend_buffer(:)
  real*8,                allocatable :: Arecv_buffer(:)
  integer(kind=int_all), allocatable :: irecv_buffer(:), jrecv_buffer(:)
  integer(kind=int_all), allocatable :: index_buffer(:), index_target(:)

  comm = ad_mat%comm
  call MPI_COMM_RANK(comm, my_id, ierr)
  call MPI_COMM_SIZE(comm, n_cpu, ierr)

  if (allocated(counts))        call tr_deallocate(counts,"counts",CAT_DMATRIX)
  if (allocated(displacements)) call tr_deallocate(displacements,"displacements",CAT_DMATRIX)

  call tr_allocate(counts,1,n_cpu,"counts",CAT_DMATRIX)
  call tr_allocate(displacements,1,n_cpu,"displacements",CAT_DMATRIX)

  call MPI_Allgather(ad_mat%nnz,1,MPI_INTEGER_ALL,counts,1,MPI_INTEGER_ALL,comm,ierr)

  displacements(1) = 0
  do i=2,n_cpu
    displacements(i) = displacements(i-1) + counts(i-1)
  enddo
  
  call ad_mat%copy_to(ac_mat, with_data=.false.) ! copy matrix parameters except arrays
  ac_mat%reduced = .true.

  ! Allocate centralized matrix
  if (associated(ac_mat%irn)) call tr_deallocatep(ac_mat%irn,"RMatrix",CAT_DMATRIX)
  if (associated(ac_mat%jcn)) call tr_deallocatep(ac_mat%jcn,"RMatrix",CAT_DMATRIX)
  if (associated(ac_mat%val)) call tr_deallocatep(ac_mat%val,"RMatrix",CAT_DMATRIX)

  ac_mat%nnz = sum(counts(1:n_cpu))

  call tr_allocatep(ac_mat%irn,Int1,ac_mat%nnz,"RMatrix",CAT_DMATRIX)
  call tr_allocatep(ac_mat%jcn,Int1,ac_mat%nnz,"RMatrix",CAT_DMATRIX)
  call tr_allocatep(ac_mat%val,Int1,ac_mat%nnz,"RMatrix",CAT_DMATRIX)

  ! --- Check if we need to split
  need_to_split = .false.
  if (ac_mat%nnz .gt. INT_MAX) need_to_split = .true.

  ! --- Allocate short-integer counts and displacements for MPI calls
  ! --- Counts still need to be copied because MPI count types are always short ints
  call tr_allocate(counts_int,1,n_cpu,"SPL_GATH_counts",CAT_DMATRIX)
  call tr_allocate(displacements_int,1,n_cpu,"SPL_GATH_displacements",CAT_DMATRIX)

  ! --- Split MPI calls
  if (need_to_split) then

    call tr_allocate(index_buffer,1,n_cpu,"SPL_GATH_index_buffer",CAT_DMATRIX)
    call tr_allocate(index_target,1,n_cpu,"SPL_GATH_index_target",CAT_DMATRIX)

    ! --- Split respective to the max send/recv
    n_split = ac_mat%nnz / INT_MAX + 1
    if (my_id .eq. 0) write(*,*) 'Warning: splitting matrix MPI centralisation', n_split

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

      i_cpu = my_id+1

      is = index_buffer(i_cpu) + 1
      ie = index_buffer(i_cpu) + counts_int(i_cpu)

      call MPI_AllgatherV(ad_mat%irn(is:ie),counts_int(my_id+1),MPI_INTEGER_ALL,ac_mat%irn, &
                          counts_int,index_target,MPI_INTEGER_ALL,comm,ierr)

      call MPI_AllgatherV(ad_mat%jcn(is:ie),counts_int(my_id+1),MPI_INTEGER_ALL,ac_mat%jcn, &
                          counts_int,index_target,MPI_INTEGER_ALL,comm,ierr)

      call MPI_AllgatherV(ad_mat%val(is:ie),counts_int(my_id+1),MPI_DOUBLE_PRECISION,ac_mat%val, &
                          counts_int,index_target,MPI_DOUBLE_PRECISION,comm,ierr)

    enddo

    call tr_deallocate(index_buffer,"SPL_GATH_index_buffer",CAT_DMATRIX)
    call tr_deallocate(index_target,"SPL_GATH_index_target",CAT_DMATRIX)

  ! --- Don't split MPI calls
  else
    ! --- Counts still need to be copied because MPI count types are always short ints
    counts_int(:) = counts(:)
    displacements_int(:) = displacements(:)

    call MPI_AllgatherV(ad_mat%irn,ad_mat%nnz,MPI_INTEGER_ALL,ac_mat%irn, &
                        counts_int,displacements_int,MPI_INTEGER_ALL,comm,ierr)
    call MPI_AllgatherV(ad_mat%jcn,ad_mat%nnz,MPI_INTEGER_ALL,ac_mat%jcn, &
                        counts_int,displacements_int,MPI_INTEGER_ALL,comm,ierr)
    call MPI_AllgatherV(ad_mat%val,ad_mat%nnz,MPI_DOUBLE_PRECISION,ac_mat%val, &
                        counts_int,displacements_int,MPI_DOUBLE_PRECISION,comm,ierr)
  endif

  ! --- Deallocate short-integer counts and displacements
  call tr_deallocate(counts_int,"SPL_GATH_counts",CAT_DMATRIX)
  call tr_deallocate(displacements_int,"SPL_GATH_displacements",CAT_DMATRIX)
  call tr_deallocate(counts,"counts",CAT_DMATRIX)
  call tr_deallocate(displacements,"displacements",CAT_DMATRIX)

  return
end subroutine matrix_split_reduce