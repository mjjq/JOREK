!> Module allowing to exchange indices (node%index) of the finite element grid in order to optimize
!! the parallelization of the free boundary part of the code. The exchanging of the indices makes
!! sure that all MPI tasks contribute to the calculation of the boundary integral, i.e., feel
!! responsible for some of the boundary elements.
module mod_exchange_indices

implicit none

logical, save              :: initialized       = .false. !< Has the module been initialized?
logical, save              :: indices_exchanged = .false. !< Have the indices been exchanged w.r.t.
                                                          !! their normal order?
integer, save, allocatable :: exchange_table(:,:)         !< Table with the indices that should be
                                                          !! exchanged
integer, save              :: len_exchange                !< How many entries in the exchange table?


contains



!> Exchange some indices of grid nodes in order to parallelize the vacuum boundary integral.
subroutine exchange_indices(node_list, my_id, n_cpu, back)
  
  use data_structure
  
  ! --- Routine parameters
  type(type_node_list), intent(inout) :: node_list
  integer,              intent(in)    :: my_id
  integer,              intent(in)    :: n_cpu
  logical,              intent(in)    :: back   !< Change the indices back (would not strictly
                                                !! be needed but is a good way to check that we
                                                !! don't exchange the wrong number of times)
  
  ! --- Local variables
  logical, parameter   :: DEBUG_OUTPUT = .false.
  integer, allocatable :: first_index_usable(:), mm(:)
  integer :: i, j, k, l, ind_max, n_bnd, ind_bnd, ind1, ind2
  logical :: skip
  
  if ( indices_exchanged .neqv. back ) then
    write(*,*) 'ERROR: Somewhere in the code you call exchange_indices too often or not often enough.'
    stop
  else if ( indices_exchanged .and. ( .not. initialized ) ) then
    write(*,*) 'ERROR: Indices have already been exchanged but has not been initialized? Internal bug!'
    stop
  end if
  
  if ( .not. initialized ) then
    
    ! --- Determine maximum index in the grid and number of boundary nodes
    ind_max = -1
    n_bnd   = 0
    do i = 1, node_list%n_nodes
      ind_max = max(ind_max, maxval(node_list%node(i)%index(:)))
      if ( node_list%node(i)%boundary > 0 ) n_bnd = n_bnd + 1
    end do
    allocate(exchange_table(n_bnd*8,2))
    allocate(first_index_usable(n_cpu))
    allocate(mm(n_cpu))
    if ( DEBUG_OUTPUT .and. (my_id == 0) ) then
      write(*,*) 'ind_max =', ind_max
      write(*,*) 'n_bnd   =', n_bnd
    end if
    
    ! --- Find out which grid nodes are usable for "exchanging indices"; store the number of the
    !     first index for each MPI rank that can be used; skip the grid center nodes as these might
    !     not have four independent grid indices
    ii: do i = 1, n_cpu
      do j = 1, node_list%n_nodes
        if ( (node_list%node(j)%index(1) > (real(i-1)/real(n_cpu))*ind_max + 1) .and. (.not. (node_list%node(j)%axis_node)) ) then
          first_index_usable(i) = node_list%node(j)%index(1)
          cycle ii
        end if
      end do
    end do ii
    
    ! --- Prepare a "table" of indices to be exchanged
    ind_bnd = 0
    j       = 1
    mm(:)   = first_index_usable(:)
    if ( DEBUG_OUTPUT ) write(*,*) 'mm before:', mm(:)
    do i = 1, node_list%n_nodes
      if ( node_list%node(i)%boundary > 0 ) then
        k = (real(j)/real(8*n_bnd))*n_cpu + 1 ! with which MPI rank to we want to exchange this index?
        if ( k == n_cpu ) cycle ! the last MPI rank doesn't need to exchange with itself
        do l = 1, 4 ! the four dofs of one node
          ind1    = node_list%node(i)%index(l) ! exchange this index
          ind2    = mm(k) + l - 1              ! with this one for which MPI rank k is responsible
          exchange_table(j,:) = (/ind1, ind2/)
          j = j + 1
          exchange_table(j,:) = (/ind2, ind1/)
          j = j + 1
          if ( DEBUG_OUTPUT ) write(*,*) 'ex: ', ind1, '<->', ind2
        end do
        mm(k) = mm(k) + 4
      end if
    end do
    len_exchange = j - 1
    if ( DEBUG_OUTPUT .and. (my_id == 0) ) write (*,*) 'len_exchange ', len_exchange
    
    deallocate(first_index_usable, mm)
    
    initialized = .true.
    
  end if ! (.not. initialized)
  
  ! --- Exchange the indices
  l = 0
  do i = 1, node_list%n_nodes
    do j = 1, 4
      do k = 1, len_exchange
        if ( node_list%node(i)%index(j) == exchange_table(k,1) ) then
          node_list%node(i)%index(j) = exchange_table(k,2)
          l = l + 1
          exit
        end if
      end do
    end do
  end do
  if ( DEBUG_OUTPUT .and. (my_id == 0) ) write (*,*) 'num exchanged ', l
  
  indices_exchanged = .not. indices_exchanged ! switch the state
  
end subroutine exchange_indices



end module mod_exchange_indices
