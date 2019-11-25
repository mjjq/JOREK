!> Exchange some indices of grid nodes in order to parallelize the vacuum boundary integral.
subroutine exchange_indices_for_vacuum(node_list, my_id, n_cpu)
  
  use data_structure
  
  ! --- Routine parameters
  type(type_node_list), intent(inout) :: node_list
  integer,              intent(in)    :: my_id
  integer,              intent(in)    :: n_cpu
  
  ! --- Local variables
  logical, parameter   :: DEBUG_OUTPUT = .false.
  integer, allocatable :: exchange_table(:,:)
  integer :: i, j, k, l, ind_max, n_bnd, ind_bnd, ind1, ind2, len_exchange
  logical :: skip
  
  ! --- Determine maximum index in the grid and number of boundary nodes
  ind_max = -1
  n_bnd   = 0
  do i = 1, node_list%n_nodes
    ind_max = max(ind_max, maxval(node_list%node(i)%index(:)))
    if ( node_list%node(i)%boundary > 0 ) n_bnd = n_bnd + 1
  end do
  allocate(exchange_table(n_bnd*8,2))
  if ( DEBUG_OUTPUT .and. (my_id == 0) ) then
    write(*,*) 'ind_max =', ind_max
    write(*,*) 'n_bnd   =', n_bnd
  end if
  
  ! --- Prepare exchange table
  ind_bnd = 0
  k       = 0
  do i = 1, node_list%n_nodes
    if ( node_list%node(i)%boundary > 0 ) then
      do j = 1, 4
        ind_bnd = ind_bnd + 1
        ind1    = node_list%node(i)%index(j)
        ind2    = real(ind_bnd)/real(4*n_bnd)*ind_max
        
        skip = (ind1==ind2)
        do l = 1, k
          skip = skip .or. ( exchange_table(l,1) == ind2 )
        end do
        if ( skip ) cycle

        if ( DEBUG_OUTPUT .and. (my_id == 0) ) write(*,*) 'ex:', ind1, ind2
        k = k + 1
        exchange_table(k,:) = (/ ind1, ind2 /)
        k = k + 1
        exchange_table(k,:) = (/ ind2, ind1 /)
        
      end do
    end if
  end do
  len_exchange = k
  if ( DEBUG_OUTPUT .and. (my_id == 0) ) write (*,*) 'len_exchange ', len_exchange
  
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
  deallocate(exchange_table)
  
end subroutine exchange_indices_for_vacuum
