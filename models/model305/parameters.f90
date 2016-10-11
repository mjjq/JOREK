!> Basic model-dependend hard-coded run parameters.
module parameters

  implicit none

  integer, parameter :: jorek_model    = 305       !< JOREK physics model

  integer, parameter :: n_var          = 8         !< number of variables
  integer, parameter :: n_dim          = 2         !< number of dimensions
  integer, parameter :: n_order        = 3         !< order of the polynomial basis
  integer, parameter :: n_tor          = 3         !< number of toroidal harmonics
  integer, parameter :: n_period       = 1         !< periodicity in toroidal direction
  integer, parameter :: n_plane        = 8         !< number of toroidal angles
  integer, parameter :: n_vertex_max   = 4         !< maximum number of corners of an element
  integer, parameter :: n_nodes_max    = 100001    !< maximum number of nodes
  integer, parameter :: n_elements_max = 100001    !< maximum number of elements
  integer, parameter :: n_boundary_max = 100001    !< maximum number of boundary elements
  integer, parameter :: n_pieces_max   = 100001    !< maximum number of line pieces describing a flux surface
  integer, parameter :: n_degrees      = n_order+1 !< degrees of freedom per variable per node
  integer, parameter :: nref_max       = 100000    !< (refinement)
  integer, parameter :: n_ref_list     = 100000    !< (refinement)
 
  !> Names of the physical variables
  character(len=11) :: variable_names(n_var) =                       &
    (/ Trim('Flux       '),Trim('Potential  '),Trim('Current    '),Trim('Vorticity  '),      &
       Trim('Density    '),Trim('Temperature'),Trim('V_parallel '), Trim('J_ECCD    ') /)

end module parameters
