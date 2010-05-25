module data_structure                                 ! contains definition of the data structure

  use parameters

  type type_node                                      ! type definition of a node (i.e. a vertex)
    real*8     :: x(n_order+1,n_dim)                  ! x,y,z coordinates of points and additional nodal geometry
    real*8     :: values(n_tor,n_order+1,n_var)
    real*8     :: deltas(n_tor,n_order+1,n_var)
    integer    :: index(n_order+1)                    ! the index in the main matrix
    integer    :: boundary                            ! = 1, 2 or 3 for boundary nodes
    integer    :: parents(2)			      ! Parent nodes (used if node is constrained)"refinement"
    integer    :: parent_elem			      ! which element do parent nodes belong to ? "refinement"
    real*8     :: ref_lambda, ref_mu		      ! Local coordinates of node inside the parent element. "refinement"
    logical    :: constrained			      ! Constrained node or not..."refinement"
  endtype type_node                                   ! x(:,1) : position, x(:,2) : vector u, x(:,3) : vector v, x(4) : vector w

  type type_node_list                                 ! type definition of a list of nodes
    integer :: n_nodes                                ! the number of nodes in the list
    integer :: n_dof                                  ! the total number of degrees of freedom
    type (type_node)     :: node(n_nodes_max)         ! an allocatable list of nodes
  endtype type_node_list

  type type_element                                   ! type definition for one elements
    integer :: vertex(n_vertex_max)                   ! the nodes of the corners
    integer :: neighbours(n_vertex_max)               ! the neighbouring elements
    real*8  :: size(n_vertex_max,n_order+1)           ! the size of the vectors at each vertex of the element
    integer :: father				      ! index of father element (0 if no father)"refinement"
    integer :: n_sons				      ! Number of sons elements"refinement"
    integer :: n_gen				      ! Generation rank of the element"refinement"
    integer :: sons(4)				      ! Sons of the element (=0 if no son)"refinement"
    integer :: contain_node(5)                        ! nodes belonging to the element"refinement"
    integer :: nref				      ! How the element has been refined (if so)"refinement"
  endtype type_element

  type type_element_list                              ! type definition for a list of elements
    integer :: n_elements                             ! the number of elements in the list
    type (type_element)  :: element(n_elements_max)   ! the list of elements
  endtype type_element_list

  type type_boundary                                  ! type definition for one boundary element (1D element)
    integer :: vertex(2)                              ! the nodes of the corners
    integer :: direction(2,2)                         ! indicates which direction of the nodes is along the boundary (2 or 3)
    integer :: element                                ! boundary element is part of this element
    integer :: side                                   ! boundary element corresponds to this side of the originating element
    real*8  :: size(2,2)                              ! the size of the vectors at each vertex of the element : size(vertex,order)
  endtype type_boundary

  type type_boundary_list                             ! type definition for a list of boundary elements
    integer :: n_boundary                             ! the number of boundary elements in the list
    type (type_boundary) :: boundary(n_boundary_max)  ! the list of boundary elements
  endtype type_boundary_list

  type type_surface                                   ! type definition for a fluxsurface (in 2D)
    integer :: n_pieces                               ! the number of pieces (each piece is a 3rd order polynomial)
    integer :: elm(n_pieces_max)                      ! the element containg the current piece
    real*8  :: s(4,n_pieces_max), t(4,n_pieces_max)   ! 4 variables per line piece of the flux surface
   endtype

  type  type_surface_list                               ! type definition for a list of surfaces
    integer                         :: n_psi            ! the number of surfaces
    real*8,allocatable              :: psi_values(:)    ! the values of the poloidal flux at the surfaces
    type (type_surface),allocatable :: flux_surfaces(:) ! the list of surfaces
  endtype

endmodule

