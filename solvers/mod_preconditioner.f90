module mod_preconditioner

  private
  public initialize_preconditioner, reset_preconditioner

  contains

  subroutine initialize_preconditioner(pc,comm_glob)
    use phys_module, only: autodistribute_modes, n_mode_families, autodistribute_ranks, centralize_harm_mat
    use mod_parameters, only: n_tor
    use data_structure, only: type_PRECOND
    use mpi_mod
    implicit none

    type(type_PRECOND) :: pc
    integer            :: comm_glob, my_id, n_cpu, ierr
    integer            :: i
    character(len=256) :: s

    pc%comm = comm_glob

    call MPI_COMM_RANK(pc%comm, my_id, ierr)
    call MPI_COMM_SIZE(pc%comm, n_cpu, ierr)

    pc%my_id = my_id
    pc%n_cpu = n_cpu
    if (pc%my_id.eq.0) write(*,*) "Initializing preconditioner"

    pc%autodistribute_ranks = autodistribute_ranks
    pc%autodistribute_modes = autodistribute_modes
    pc%mat%row_distributed  = .not.centralize_harm_mat

    if (pc%autodistribute_modes) then
      pc%n_mode_families = (n_tor + 1)/2
    else
      pc%n_mode_families = n_mode_families
    endif

    call distribute_ranks(n_cpu, pc)

    call create_communicators(pc)
    pc%mat%comm = pc%MPI_COMM_N ! communicator for PC matrix distribution

    call distribute_modes(pc)

    pc%initialized = .true.

    if (my_id.eq.0) then
      do i=1, pc%n_mode_families
        write(s,'(A17,i4,A12)') " mode_family_id: ", i, " MPI ranks: "
        write(*,*) trim(s), pc%mode_families_ranks(i,1:pc%ranks_per_family(i))
      enddo
      do i=1, pc%n_mode_families
        write(s,'(A17,i4,A9,f6.2)') " mode_family_id: ", i, " weight: ", pc%row_factor
        write(*,*) trim(s), " modes:", pc%mode_families_modes(i,1:pc%modes_per_family(i))
      enddo
    endif

    return

  end subroutine initialize_preconditioner

!> Set up MPI communicators for mode families and corresponding masters
  subroutine create_communicators(pc)
    use data_structure, only: type_PRECOND
    use mpi_mod
    implicit none

    type(type_PRECOND) :: pc

    integer, allocatable :: i_tor(:), ranks_tmp(:)
    integer :: my_id, n_cpu, ierr

    my_id = pc%my_id
    n_cpu = pc%n_cpu

    call MPI_COMM_SPLIT(pc%comm, pc%family_id, my_id, pc%MPI_COMM_N, ierr)
    if (ierr.ne.0) then
      write(*,*) "Error in creating MPI_COMM_N"
      call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
    endif
    call MPI_COMM_RANK(pc%MPI_COMM_N, pc%my_id_n, ierr)
    call MPI_COMM_SIZE(pc%MPI_COMM_N, pc%n_cpu_n, ierr)

    allocate(i_tor(n_cpu)); i_tor = 0
    i_tor(my_id + 1) = pc%my_id_n
    call MPI_Allreduce(MPI_IN_PLACE,i_tor,n_cpu,MPI_INT,MPI_SUM,pc%comm,ierr)
    call MPI_COMM_SPLIT(pc%comm,i_tor(my_id+1),my_id,pc%MPI_COMM_TRANS,ierr)

    pc%n_masters = pc%n_mode_families
    allocate(ranks_tmp(pc%n_masters)); ranks_tmp=0;

    if (pc%my_id_n.eq.0) ranks_tmp(pc%family_id) = my_id
    call MPI_AllReduce(MPI_IN_PLACE,ranks_tmp,pc%n_masters,MPI_INT,MPI_SUM,pc%comm,ierr)
    call MPI_COMM_GROUP(pc%comm,pc%MPI_GROUP_WORLD,ierr)
    call MPI_GROUP_INCL(pc%MPI_GROUP_WORLD,pc%n_masters,ranks_tmp,pc%MPI_GROUP_MASTER,ierr)
    call MPI_COMM_CREATE(pc%comm,pc%MPI_GROUP_MASTER,pc%MPI_COMM_MASTER,ierr)

    if (pc%my_id_n .eq. 0) then
     call MPI_COMM_RANK(pc%MPI_COMM_MASTER, pc%my_id_master, ierr)
    endif

    if ((my_id.eq.0).and.(pc%my_id_n.ne.0)) then
      write(*,*) "Error in creating communicators: my_id==0 must have my_id_n==0"
      call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
    endif

    deallocate(i_tor, ranks_tmp)

    return

  end subroutine create_communicators

  !> Distribute toroidal modes among mode families
  subroutine distribute_modes(pc)
    use data_structure, only: type_PRECOND
    use phys_module, only: modes_per_family, mode_families_modes, weights_per_family
    implicit none

    type(type_PRECOND) :: pc
    integer            :: i, j, n_fam_max

    allocate(pc%modes_per_family(pc%n_mode_families))

    if (pc%autodistribute_modes) then
      pc%row_factor = 1.0
      pc%modes_per_family(1) = 1
      if (pc%n_mode_families>1) pc%modes_per_family(2:pc%n_mode_families) = 2
    else
      do i = 1, pc%n_mode_families
        pc%row_factor = weights_per_family(i)
        pc%modes_per_family(i) = modes_per_family(i)
      enddo
    endif

    n_fam_max = 1
    do i = 1, pc%n_mode_families
      n_fam_max = max(n_fam_max,pc%modes_per_family(i))
    enddo

    allocate(pc%mode_families_modes(pc%n_mode_families,n_fam_max))
    pc%mode_families_modes(:,:) = -1

    if (pc%autodistribute_modes) then
      pc%mode_families_modes(1,1) = 1
      if (pc%n_mode_families.gt.1) then
        do i = 2, pc%n_mode_families
          pc%mode_families_modes(i,1) =  (i - 1)*2
          pc%mode_families_modes(i,2) =  (i - 1)*2 + 1
        enddo
      endif
    else
      do i = 1, pc%n_mode_families
        do j = 1, modes_per_family(i)
          pc%mode_families_modes(i,j) = mode_families_modes(i,j)
        enddo
      enddo
    endif

    pc%mode_set_n = pc%modes_per_family(pc%family_id)
    allocate(pc%mode_set(pc%mode_set_n))
    pc%mode_set(1:pc%mode_set_n) = pc%mode_families_modes(pc%family_id,1:pc%mode_set_n)

  end subroutine distribute_modes

  !> Distribute MPI ranks among mode families
  subroutine distribute_ranks(n_cpu,pc)
    use data_structure, only: type_PRECOND
    use phys_module, only: ranks_per_family
    implicit none

    type(type_PRECOND) :: pc
    integer, intent(in)  :: n_cpu
    integer :: mcpu, r, i, j

    allocate(pc%rank_range(pc%n_mode_families + 1))
    allocate(pc%rank_id(n_cpu))
    allocate(pc%ranks_per_family(pc%n_mode_families))
    allocate(pc%mode_families_ranks(pc%n_mode_families,n_cpu))

    do i = 1, pc%n_mode_families
      do j = 1, n_cpu
        pc%mode_families_ranks(i,j) = -1
      enddo
    enddo

    mcpu = n_cpu/pc%n_mode_families
    r = mod(n_cpu,pc%n_mode_families)

    if (pc%autodistribute_ranks) then
      do i=1, pc%n_mode_families
        pc%ranks_per_family(i) = mcpu
        if ((r.gt.0).and.(i.le.r))  pc%ranks_per_family(i) = pc%ranks_per_family(i) + 1 ! add extra rank if avaiable
      enddo
    else
      do i = 1, pc%n_mode_families
        pc%ranks_per_family(i) = ranks_per_family(i)
      enddo
    endif

    pc%rank_range(1) = 1
    do i = 2, pc%n_mode_families+1
      pc%rank_range(i) = pc%rank_range(i-1) + pc%ranks_per_family(i-1)
    enddo

    ! check for consistency
    r = 0
    do i= 2, pc%n_mode_families + 1
      r = r + pc%rank_range(i) - pc%rank_range(i-1)
    enddo
    if (r.ne.n_cpu) then
      write(*,*) "Error in distribution of ranks"
      call exit(0)
    endif

    do i = 1, n_cpu
      do j = 2, pc%n_mode_families + 1
        if ((i.ge.pc%rank_range(j-1)).and.(i.lt.pc%rank_range(j))) then
          pc%rank_id(i) = j - 1
          exit
        endif
      enddo
    enddo

    do j=1,pc%n_mode_families
      r = 0
      do i = 1, n_cpu
        if (pc%rank_id(i).eq.j) then
          r = r + 1
          pc%mode_families_ranks(j,r) = i - 1
        endif
      enddo
    enddo

    pc%family_id = pc%rank_id(pc%my_id + 1)

    return

  end subroutine distribute_ranks

!> Deallocate arrays and reset to the default values
  subroutine reset_preconditioner(pc)
    use data_structure, only: type_PRECOND
    implicit none

    type(type_PRECOND) :: pc !, pc_def

    if (.not.pc%initialized) then

      write(*,*) "Preconditioner is not initialized"

    else

      if (pc%analyzed) then

        deallocate(pc%rhs%val)
        deallocate(pc%row_index)
        deallocate(pc%send_counts, pc%recv_counts)
        deallocate(pc%send_disp, pc%recv_disp)
        deallocate(pc%istart, pc%ifinish)
        deallocate(pc%n_per_rank)

        deallocate(pc%mat%val)
        deallocate(pc%mat%irn)
        deallocate(pc%mat%jcn)
        if (pc%mat%scaled) deallocate(pc%mat%column_scaling)
        pc%mat%scaled = .false.
        pc%mat%row_distributed = .false.
        pc%mat%col_distributed = .false.
        pc%mat%indexing = 1
        pc%mat%block_size = 1
        pc%analyzed = .false.

      endif

      deallocate(pc%mode_families_ranks)
      deallocate(pc%mode_families_modes)
      deallocate(pc%rank_id)
      deallocate(pc%mode_set)

      pc%initialized = .false.

    endif

    return

  end subroutine reset_preconditioner


end module mod_preconditioner











!
!
!
!
!
! !> Determine mapping from local to globar row index for the RHS
!  subroutine map_row_index(ndof)
!
!    use mod_parameters, only: n_tor
!    use tr_module
!    use mod_integer_types
!
!    implicit none
!
!    integer(kind=int_all), intent(in) :: ndof
!    integer(kind=int_all) :: i, ndof_family
!    integer(kind=int_all), parameter   :: Int1=1
!    integer :: im
!
!    ndof_family = my_mode_set_n*ndof/n_tor
!
!    if (allocated(my_row_index)) call tr_deallocate(my_row_index,"my_row_index",CAT_DMATRIX)
!    call tr_allocate(my_row_index,Int1,ndof_family,"my_row_index",CAT_DMATRIX)
!    !if (allocated(my_row_factor)) call tr_deallocate(my_row_factor,"my_row_index",CAT_DMATRIX)
!    !call tr_allocate(my_row_factor,Int1,ndof_family,"my_row_factor",CAT_DMATRIX)
!    my_row_factor = weights_per_family(my_family_id)
!
!    do i = 0, ndof/n_tor - 1
!      do im=1,my_mode_set_n
!        my_row_index(im+i*my_mode_set_n) =  my_mode_set(im) + i*n_tor
!      enddo
!    enddo
!
!    if (autodistribute_modes) my_row_factor = 1.0
!
!    write(*,*) "my_family_id my_row_factor", my_family_id, weights_per_family(my_family_id)
!
!    return
!  end subroutine map_row_index
!
!  !> Check mode families consistency
!  subroutine check_preconditioner_consistency
!
!    use mpi
!    use mod_parameters, only: n_tor
!
!    implicit none
!
!    integer :: i, j, my_id, n_cpu, ierr
!    logical :: cck
!
!    call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
!    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)
!
!    if (n_mode_families.gt.n_cpu) then
!      if (my_id.eq.0) write(*,*) "Error: number of cpu must be >= number of mode families", n_mode_families
!      call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
!    endif
!
!    if (.not.autodistribute_modes) then
!      ! check if mode families are specified
!      if (n_mode_families.le.0) then
!        if (my_id.eq.0) write(*,*) "Error: number of mode families is not specified or <=0"
!        call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
!      endif
!
!      ! check if mode families consist of correct modes
!      do i = 1, n_mode_families
!        do j = 1, modes_per_family(i)
!          if ((mode_families_modes(i,j).lt.1).or.(mode_families_modes(i,j).gt.n_tor)) then
!            if (my_id.eq.0) write(*,*) "Error: incorrect specification of mode_families_modes", i
!            call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
!          endif
!        enddo
!      enddo
!
!    endif
!
!    if (.not.autodistribute_ranks) then
!      ! check if numbers of ranks are provided
!      do i = 1, n_mode_families
!        if (ranks_per_family(i).le.0) then
!          if (my_id.eq.0) write(*,*) "Error: ranks_per_family is not correctly specified", i
!          call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
!        endif
!      enddo
!      ! check if sum of all ranks per family is equal to the total number of ranks
!      if (sum(ranks_per_family(1:n_mode_families)).ne.n_cpu) then
!        if (my_id.eq.0) write(*,*) "Error: sum of all ranks_per_family must be equal to total number of MPI ranks"
!        call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
!      endif
!    endif
!
!#ifndef USE_PASTIX6
!    if ((.not.use_strumpack).and.(.not.centralize_harm_mat)) then
!      if (my_id.eq.0) write(*,*) "Warning: PC matrix distribution is only supported by STRUMPACK"
!      centralize_harm_mat = .true.
!    endif
!#endif
!
!    if ((use_strumpack).and.(centralize_harm_mat)) then
!      if (my_id.eq.0) write(*,*) "Warning: centralization of PC matrix is not advized when using STRUMPACK"
!    endif
!
!  end subroutine check_preconditioner_consistency
!!#endif
!
!
!

!end module preconditioner_module
