module parameters

  implicit none

  integer      :: n_var, n_dim, n_order, n_tor, n_period, n_plane, n_vertex_max
  integer      :: n_nodes_max ,n_elements_max, n_boundary_max, n_pieces_max,n_degrees
  integer      :: nref_max ,n_ref_list                 !  ref
  parameter (n_var          = 6)                       ! the number of variables
  parameter (n_dim          = 2)                       ! the number of dimensions
  parameter (n_order        = 3)                       ! order of the polynomial basis
  parameter (n_tor          = 3)                       ! the number of toroidal harmonics
  parameter (n_period       = 1)                       ! periodicity in toroidal direction
  parameter (n_plane        = 4)                       ! the number of toroidal angles
  parameter (n_vertex_max   = 4)                       ! the maximum number of corners of an element
  parameter (n_nodes_max    = 15001)                   ! the maximum number of nodes
  parameter (n_elements_max = 15001)                   ! the maximum number of elements
  parameter (n_boundary_max = 1001)                    ! the maximum number of boundary elements
  parameter (n_pieces_max   = 10001)                   ! the maximum number of line pieces describing a flux surface
  parameter (n_degrees      = n_order+1)               ! degrees of freedom per variable per node
  parameter (nref_max = 10000, n_ref_list=10000)       !  ref                                                         
  character*11 :: variable_names(n_var)
 
  parameter (variable_names = (/ 'Flux       ','Potential  ','Current    ','Vorticity  ','Density    ','Temperature' /))

endmodule
