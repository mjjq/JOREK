!>  Module that contais functions to fill IDSs from JOREK data
module mod_jorek2IMAS

#ifdef USE_IMAS
  use ids_schemas !, only: ids_equilibrium
  use ids_routines, only: imas_open_env, &
     imas_create_env, imas_close, ids_get, ids_put
#endif

  use mod_parameters 
  use data_structure
  
  implicit none
  
 
  public
 
  contains



  ! --- Fills Bezier coefficients in GGD
  subroutine fill_Bezier_coefficients( ggd_scalar, node_list, var_index, grid_ind, grid_sub_ind, res_fact )
  
    implicit none
  
    ! --- External parameters
    type(ids_generic_grid_scalar),  intent(inout) ::  ggd_scalar
    type (type_node_list), intent(in)  :: node_list
    integer,               intent(in)  :: var_index, grid_ind, grid_sub_ind
    real*8,                intent(in)  :: res_fact
  
    ! --- Local parameters
    integer :: num_nodes, idof, inode, icoeff
    
    num_nodes = node_list%n_nodes
  
    if ( associated(ggd_scalar%coefficients) ) then
      deallocate( ggd_scalar%coefficients )
      allocate( ggd_scalar%coefficients(n_tor, num_nodes*(n_order+1)) ) 
    else
      allocate( ggd_scalar%coefficients(n_tor, num_nodes*(n_order+1)) ) 
    endif
  
    do inode=1, num_nodes
      do idof=1, n_order+1
  
        icoeff = inode + (idof-1)*num_nodes
  
        ggd_scalar%coefficients(:,icoeff)=node_list%node(inode)%values(:,idof,var_index)
  
      enddo
    enddo
  
    ggd_scalar%grid_index = grid_ind
    ggd_scalar%grid_subset_index = grid_sub_ind
 
   ! --- Re-scale coefficients
    ggd_scalar%coefficients = ggd_scalar%coefficients * res_fact
  
  end subroutine fill_Bezier_coefficients
  
  
  
  
  !< Fills JOREK grid into GGD
  subroutine grid2ggd( grid, node_list, element_list )
  
    implicit none
  
    ! --- External parameters
    type(ids_generic_grid_aos3_root),     pointer   :: grid
    type (type_node_list),    intent(in)            :: node_list
    type (type_element_list), intent(in)            :: element_list
  
  
    ! --- Local parameters
    type(ids_generic_grid_dynamic_space), pointer   ::  space_RZ
    type(ids_generic_grid_dynamic_space), pointer   ::  space_fourier
    type(ids_generic_grid_dynamic_space_dimension), pointer :: ids_cells
    
    integer:: idx, shot_number, run_number, num_nodes, num_cells   
    integer :: gs_index, i, j
    integer, allocatable :: vertex_elm_array(:,:)
    real*8,  allocatable :: RZ(:,:)
    
    integer :: itor, idof, n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub
  
    ! --- create vertex - elements array
    allocate(  vertex_elm_array(n_vertex_max, element_list%n_elements)  )
    do i=1, element_list%n_elements
      vertex_elm_array(:,i) = element_list%element(i)%vertex
    enddo 
  
    ! Get values of R and Z at the nodes
    allocate( RZ(node_list%n_nodes, 2) )
    do i=1, node_list%n_nodes
      RZ(i,:) = node_list%node(i)%x(1,1,:)  
    enddo
  
    ! Write grid geometry
    allocate(  grid%space(2)                           )
    allocate(  grid%space(1)%objects_per_dimension(4)  )
    allocate(  grid%space(1)%coordinates_type(2)       )
  
    ! Set coordinates type to [R, Z]
    grid%space(1)%coordinates_type = (/ 4, 3 /)
  
    allocate(grid%identifier%description(1))
    allocate(grid%identifier%name(1))
    grid%identifier%description(1) = "Mesh JOREK output HDF5 file grid with quantities"
    grid%identifier%name = "JOREK output HDF5 file grid with quantities"
    grid%identifier%index = 1
  
    num_nodes = size(RZ,1)
    space_RZ  => grid%space(1)
  
    ! Fill simplified grid nodes (uses geometry instead of geometry_2D and misses Bezier representation) 
    allocate( space_RZ%objects_per_dimension(1)%object(num_nodes) )  ! Allocate to number of nodes
    do i=1, num_nodes 
      allocate( space_RZ%objects_per_dimension(1)%object(i)%geometry(2) ) ! Allocate dimensions per each node
      space_RZ%objects_per_dimension(1)%object(i)%geometry(:) = RZ(i,:)
    enddo
  
    ! Fill dummy variables for 1D elements (edges)
    allocate( space_RZ%objects_per_dimension(2)%object(1) )          ! Allocate just one edge    
    allocate( space_RZ%objects_per_dimension(2)%object(1)%nodes(1))  
    space_RZ%objects_per_dimension(2)%object(1)%nodes(1) = 0
  
    ! Fill JOREK 2D elements (or cells)
    num_cells = element_list%n_elements
    ids_cells => space_RZ%objects_per_dimension(3)
    allocate(    ids_cells%object(num_cells) )
    do i=1, num_cells
      allocate(  ids_cells%object(i)%nodes(n_vertex_max)  )
      ids_cells%object(i)%nodes(:) = vertex_elm_array(:,i) 
    enddo
  
    ! Writing grid_subsets
    allocate(grid%grid_subset(2))  ! 2 grid subsets 
  
    ! Subset for points
    gs_index = 1
  
    allocate( grid%grid_subset(gs_index)%identifier%name(1)         )
    allocate( grid%grid_subset(gs_index)%identifier%description(1)  )
    grid%grid_subset(gs_index)%identifier%name(1)        = "Nodes"
    grid%grid_subset(gs_index)%identifier%index          = gs_index 
    grid%grid_subset(gs_index)%identifier%description(1) = "All points/nodes/vertices/0D objects in the domain."
    grid%grid_subset(gs_index)%dimension                 = 1
   
    allocate( grid%grid_subset(gs_index)%element(num_nodes) )
   
    do i=1, num_nodes  
      allocate(  grid%grid_subset(gs_index)%element(i)%object(1)   )
      grid%grid_subset(gs_index)%element(i)%object(1)%space     = 1
      grid%grid_subset(gs_index)%element(i)%object(1)%index     = i 
      grid%grid_subset(gs_index)%element(i)%object(1)%dimension = 1
    enddo 
  
    ! Subset for cells
    gs_index = 2
  
    allocate(grid%grid_subset(gs_index)%identifier%name(1))
    allocate(grid%grid_subset(gs_index)%identifier%description(1))
  
    grid%grid_subset(gs_index)%identifier%name(1)        = "2D cells "
    grid%grid_subset(gs_index)%identifier%index          = gs_index 
    grid%grid_subset(gs_index)%identifier%description(1) = "All points/nodes/vertices/0D objects in the domain."
    grid%grid_subset(gs_index)%dimension                 = 3   ! A bit confusing, but 1 means point, 2 line and 3 surface
   
    allocate( grid%grid_subset(gs_index)%element(num_cells) )
  
    do i=1, num_cells
      allocate(  grid%grid_subset(gs_index)%element(i)%object(1)   )
      grid%grid_subset(gs_index)%element(i)%object(1)%space     = 1
      grid%grid_subset(gs_index)%element(i)%object(1)%index     = i 
      grid%grid_subset(gs_index)%element(i)%object(1)%dimension = 3
    enddo 
  
    ! Fill toroidal space (must be adapted for multiple time slices?)
    space_fourier  => grid%space(2)
    allocate(    space_fourier%coordinates_type(1)    )
    allocate(    space_fourier%identifier%description(1)  )
    space_fourier%coordinates_type(1) = 5          ! The coordinate type is 5, phi angle
    space_fourier%geometry_type%index = n_period   ! Fourier periodicity
    space_fourier%identifier%description(1) = "Toroidal space"             
  
    allocate(  space_fourier%objects_per_dimension(1)               )  ! We have only one dimension of
    allocate(  space_fourier%objects_per_dimension(1)%object(n_tor) )  ! toroidal harmonics
  
    do i=1, n_tor
      allocate(  space_fourier%objects_per_dimension(1)%object(i)%geometry(1) )  ! toroidal harmonics
      space_fourier%objects_per_dimension(1)%object(i)%geometry(1) = i
    enddo
  
    ! Fill in grid Bezier coefficients
    ! Needs generalization for JOREK 3D STELLERATOR EXTENSION!!
    space_RZ%geometry_type%index = 0  ! Standard geometry (non Fourier)
    do i=1, num_nodes
      allocate( space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(2,n_order+1) )
      space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(1,:) = node_list%node(i)%x(1,:,1 )   ! R dofs
      space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(2,:) = node_list%node(i)%x(1,:,2 )   ! Z dofs
    enddo
  
    ! JOREK element sizes
    do i=1, num_cells
      allocate( space_RZ%objects_per_dimension(3)%object(i)%geometry_2d(n_order+1, n_vertex_max) )
      do j=1, n_vertex_max
        space_RZ%objects_per_dimension(3)%object(i)%geometry_2d(:,j) = element_list%element(i)%size(j,:)   
      enddo
    enddo
  
  end subroutine grid2ggd




end module mod_jorek2IMAS
