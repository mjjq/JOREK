!> Program to convert a JOREK2 restart file into binary VTK format
program jorek_convert_grid_order

use mod_parameters
use data_structure
use phys_module
use mod_import_restart
use mod_export_restart
use basis_at_gaussian
use mod_grid_conversions
use mod_poiss
use mod_boundary
use mpi_mod
use live_data

implicit none

type (type_node_list)   ,     pointer :: node_list
type (type_element_list),     pointer :: element_list
type (type_bnd_element_list), pointer :: bnd_elm_list    
type (type_bnd_node_list),    pointer :: bnd_node_list 

type (type_node_list)   ,     pointer :: newnode_list
type (type_element_list),     pointer :: newelement_list

integer :: my_id, k_tor, i, ierr, k, iv, j, index
integer :: i_elm, i_vertex, i_degrees, i_node, i_var, i_tor, i_dim

integer              :: n_axis_nodes
integer, allocatable :: axis_nodes(:)
logical, allocatable :: i_am_axis_node (:)
integer              :: first_xpoint_nodes(2)
integer, allocatable :: n_parents(:)         ! for each node, want the number of parent elements
integer, allocatable :: node_parents(:,:)    ! for each node, want to know the 4 parent elements
integer, allocatable :: parent_elm_node(:,:) ! for each node, want to know the corresonding vertex for the 4 parent elements

integer :: required, provided, StatInfo
integer :: rank, comm_size, n_cpu
integer :: MPI_COMM_N, MPI_GROUP_MASTER, MPI_GROUP_WORLD, MPI_COMM_MASTER, MPI_COMM_TRANS

write(*,*) '***************************************'
write(*,*) '*      jorek_convert_grid_order       *'
write(*,*) '***************************************'

#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif
call MPI_Init_thread(required, provided, StatInfo)
call init_threads()  ! on some systems init_threads needs to come after mpi_init_thread
call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
n_cpu = comm_size
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
my_id = rank
if ( my_id == 0 ) call init_live_data()

allocate(node_list)
allocate(element_list)
allocate(newnode_list)
allocate(newelement_list)
allocate(bnd_elm_list)
allocate(bnd_node_list)

! --- Initialisation
my_id     = 0
call initialise_and_broadcast_parameters(my_id, "__NO_FILENAME__")
do k_tor=1, n_tor
  mode(k_tor) = + int(k_tor / 2) * n_period
enddo
call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)
call initialise_basis                              ! define the basis functions at the Gaussian points
call tr_meminit(my_id, n_cpu)

! --- Copy nodes/elements
newnode_list%n_nodes = node_list%n_nodes
newelement_list%n_elements = element_list%n_elements

do i_elm = 1,newelement_list%n_elements
  newelement_list%element(i_elm)%vertex(1:n_vertex_max)           = element_list%element(i_elm)%vertex(1:n_vertex_max)        
  newelement_list%element(i_elm)%neighbours(1:n_vertex_max)       = element_list%element(i_elm)%neighbours(1:n_vertex_max)    
  do i_degrees = 1,4
    newelement_list%element(i_elm)%size(1:n_vertex_max,i_degrees) = element_list%element(i_elm)%size(1:n_vertex_max,i_degrees)
  enddo
  newelement_list%element(i_elm)%father                           = element_list%element(i_elm)%father                      
  newelement_list%element(i_elm)%n_sons                           = element_list%element(i_elm)%n_sons                      
  newelement_list%element(i_elm)%n_gen                            = element_list%element(i_elm)%n_gen                       
  newelement_list%element(i_elm)%sons(1:4)                        = element_list%element(i_elm)%sons(1:4)                     
  newelement_list%element(i_elm)%contain_node(1:5)                = element_list%element(i_elm)%contain_node(1:5)             
  newelement_list%element(i_elm)%nref                             = element_list%element(i_elm)%nref                        
enddo

newnode_list%n_dof = node_list%n_dof
do i_node = 1,newnode_list%n_nodes
  do i_degrees = 1,4
    do i_dim = 1,n_dim
      newnode_list%node(i_node)%x(1:n_coord_tor,i_degrees,i_dim) = node_list%node(i_node)%x(1:n_coord_tor,i_degrees,i_dim)
    enddo
    do i_var = 1,n_var
      newnode_list%node(i_node)%values(1:n_tor,i_degrees,i_var)  = node_list%node(i_node)%values(1:n_tor,i_degrees,i_var) 
      newnode_list%node(i_node)%deltas(1:n_tor,i_degrees,i_var)  = 0.d0 !node_list%node(i_node)%deltas(1:n_tor,i_degrees,i_var) 
    enddo
#ifdef fullmhd
    newnode_list%node(i_node)%psi_eq(i_degrees)   = node_list%node(i_node)%psi_eq(i_degrees)   
    newnode_list%node(i_node)%Fprof_eq(i_degrees) = node_list%node(i_node)%Fprof_eq(i_degrees) 
#elif altcs
    newnode_list%node(i_node)%psi_eq(i_degrees)   = node_list%node(i_node)%psi_eq(i_degrees)
#endif
    newnode_list%node(i_node)%index(i_degrees)    = node_list%node(i_node)%index(i_degrees)
  enddo
  newnode_list%node(i_node)%boundary       = node_list%node(i_node)%boundary       
  newnode_list%node(i_node)%boundary_index = node_list%node(i_node)%boundary_index 
  newnode_list%node(i_node)%axis_node      = node_list%node(i_node)%axis_node      
  newnode_list%node(i_node)%parents(1:2)   = node_list%node(i_node)%parents(1:2)     
  newnode_list%node(i_node)%parent_elem    = node_list%node(i_node)%parent_elem    
  newnode_list%node(i_node)%ref_lambda     = node_list%node(i_node)%ref_lambda     
  newnode_list%node(i_node)%ref_mu         = node_list%node(i_node)%ref_mu         
  newnode_list%node(i_node)%constrained    = node_list%node(i_node)%constrained    
enddo

! --- Reset vector sizes
do i_node = 1,newnode_list%n_nodes
  newnode_list%node(i_node)%values(1,5:n_degrees,1) = 0.d0
  do i_degrees = 2,4
    do i_dim = 1,n_dim
      newnode_list%node(i_node)%x(1:n_coord_tor,i_degrees,i_dim) = &
        newnode_list%node(i_node)%x(1:n_coord_tor,i_degrees,i_dim) * float(3)/float(n_order)
    enddo
    newnode_list%node(i_node)%values(1,i_degrees,1) = &
      newnode_list%node(i_node)%values(1,i_degrees,1) * float(3)/float(n_order)
  enddo
enddo
! --- Convert to higher order
call set_high_order_sizes(newelement_list)
call approximate_2nd_derivatives(newnode_list,newelement_list)
do i=1,newnode_list%n_nodes
  newnode_list%node(i)%x(1,7:n_degrees,:) = 0.d0
enddo

! --- Get parent elements
allocate( axis_nodes     (  newnode_list%n_nodes) )
allocate( i_am_axis_node (  newnode_list%n_nodes) )
allocate( n_parents      (  newnode_list%n_nodes) )
allocate( node_parents   (8,newnode_list%n_nodes) )
allocate( parent_elm_node(8,newnode_list%n_nodes) )
n_parents       = 0
node_parents    = 0
parent_elm_node = 0

! --- Find parent elements
n_axis_nodes = 0
first_xpoint_nodes = 0
do i_node = 1, newnode_list%n_nodes
  n_parents(i_node) = 0
  do i_elm = 1, newelement_list%n_elements
    do i_vertex = 1, n_vertex_max
      if (newelement_list%element(i_elm)%vertex(i_vertex) .eq. i_node) then
        n_parents(i_node) = n_parents(i_node) + 1
        node_parents   (n_parents(i_node),i_node) = i_elm
        parent_elm_node(n_parents(i_node),i_node) = i_vertex
        exit
      endif
    enddo
  enddo
  ! --- The axis nodes
  i_am_axis_node(i_node) = .false.
  if (n_parents(i_node) .gt. 10) then
    n_axis_nodes = n_axis_nodes + 1
    axis_nodes(n_axis_nodes) = i_node
    i_am_axis_node(i_node) = .true.
  endif
  ! --- The axis nodes
  if (xpoint) then
    if ( (n_parents(i_node) .eq. 2) .and. (newnode_list%node(i_node)%boundary .eq. 0) .and. (first_xpoint_nodes(1) .eq. 0) ) then
      first_xpoint_nodes(1) = i_node
      if (xcase .eq. 3) first_xpoint_nodes(2) = i_node + 4
    endif
  endif
enddo

! --- Redefine indices
index = 0
do i=1,newnode_list%n_nodes

  do k=1,n_degrees

    index = index + 1
    newnode_list%node(i)%index(k) = index

    ! Remove all but one node at axis
    if (force_central_node) then
      if ( i_am_axis_node(i) .and. (k.eq.1)) then
        newnode_list%node(i)%index(k) = newnode_list%node(axis_nodes(1))%index(1)
        index = index - 1
      endif
    endif
     
    ! Remove all but one node at first Xpoint
    if (xpoint) then
      if (i .eq. first_xpoint_nodes(1)+1) then
        if ( (k.eq.1) .or. (k.eq.3) .or. (k.eq.6) ) then
          newnode_list%node(i)%index(k) = newnode_list%node(first_xpoint_nodes(1))%index(k)
          index = index - 1
        endif
      endif
      if (i .eq. first_xpoint_nodes(1)+2) then
        if ( (k.eq.1) .or. (k.eq.2) .or. (k.eq.5) ) then
          newnode_list%node(i)%index(k) = newnode_list%node(first_xpoint_nodes(1)+1)%index(k)
          index = index - 1
        endif
      endif
      if (i .eq. first_xpoint_nodes(1)+3) then
        if ( (k.eq.1) .or. (k.eq.2) .or. (k.eq.5) ) then
          newnode_list%node(i)%index(k) = newnode_list%node(first_xpoint_nodes(1))%index(k)
          index = index - 1
        endif
        if ( (k.eq.3) .or. (k.eq.6) ) then
          newnode_list%node(i)%index(k) = newnode_list%node(first_xpoint_nodes(1)+2)%index(k)
          index = index - 1
        endif
      endif
    endif

    ! Remove all but one node at second Xpoint
    if ( xpoint .and. (xcase .eq. 3) ) then
      if (i .eq. first_xpoint_nodes(2)+1) then
        if ( (k.eq.1) .or. (k.eq.3) .or. (k.eq.6) ) then
          newnode_list%node(i)%index(k) = newnode_list%node(first_xpoint_nodes(2))%index(k)
          index = index - 1
        endif
      endif
      if (i .eq. first_xpoint_nodes(2)+2) then
        if ( (k.eq.1) .or. (k.eq.2) .or. (k.eq.5) ) then
          newnode_list%node(i)%index(k) = newnode_list%node(first_xpoint_nodes(2)+1)%index(k)
          index = index - 1
        endif
      endif
      if (i .eq. first_xpoint_nodes(2)+3) then
        if ( (k.eq.1) .or. (k.eq.2) .or. (k.eq.5) ) then
          newnode_list%node(i)%index(k) = newnode_list%node(first_xpoint_nodes(2))%index(k)
          index = index - 1
        endif
        if ( (k.eq.3) .or. (k.eq.6) ) then
          newnode_list%node(i)%index(k) = newnode_list%node(first_xpoint_nodes(2)+2)%index(k)
          index = index - 1
        endif
      endif
    endif
  
  enddo  
enddo

! --- Axis indices
if (fix_axis_nodes) then
  do k=1, newelement_list%n_elements
    do iv=1,4
      j = newelement_list%element(k)%vertex(iv)
      if (newnode_list%node(j)%axis_node) then
        newelement_list%element(k)%size(iv,3) = 0.d0
        newelement_list%element(k)%size(iv,4) = 0.d0
      endif
    enddo
  enddo
  if (n_order .ge. 5) call set_high_order_sizes_on_axis(newnode_list,newelement_list)
endif

newnode_list%node(first_xpoint_nodes(1)  )%x(1,5:n_degrees,:) = 0.d0
newnode_list%node(first_xpoint_nodes(1)+1)%x(1,5:n_degrees,:) = 0.d0
newnode_list%node(first_xpoint_nodes(1)+2)%x(1,5:n_degrees,:) = 0.d0
newnode_list%node(first_xpoint_nodes(1)+3)%x(1,5:n_degrees,:) = 0.d0
if (xcase .eq. 3) then
  newnode_list%node(first_xpoint_nodes(2)  )%x(1,5:n_degrees,:) = 0.d0
  newnode_list%node(first_xpoint_nodes(2)+1)%x(1,5:n_degrees,:) = 0.d0
  newnode_list%node(first_xpoint_nodes(2)+2)%x(1,5:n_degrees,:) = 0.d0
  newnode_list%node(first_xpoint_nodes(2)+3)%x(1,5:n_degrees,:) = 0.d0
endif

! --- Copy new nodes/elements
node_list%n_nodes = newnode_list%n_nodes
node_list%node(1:node_list%n_nodes) = newnode_list%node(1:node_list%n_nodes)
element_list%n_elements = newelement_list%n_elements
element_list%element(1:element_list%n_elements) = newelement_list%element(1:element_list%n_elements)

! --- Redo equilibrium
call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.) 
call export_boundary(node_list, bnd_elm_list, bnd_node_list)
call broadcast_boundary(my_id,bnd_elm_list,bnd_node_list) 
call equilibrium(my_id, node_list, element_list, bnd_node_list, bnd_elm_list, xpoint,xcase, .false.)

! --- Redo initial conditions
call initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint,xcase)

! --- Export equilibrium
call export_restart(node_list, element_list, 'jorek_restart')

write(*,*)'done'

deallocate(node_list)
deallocate(element_list)
deallocate(newnode_list)
deallocate(newelement_list)
deallocate(bnd_elm_list)
deallocate(bnd_node_list)
deallocate(axis_nodes)
deallocate(i_am_axis_node)
deallocate(n_parents)
deallocate(node_parents)
deallocate(parent_elm_node)

call MPI_FINALIZE(IERR)                                ! clean up MPI

end program jorek_convert_grid_order
