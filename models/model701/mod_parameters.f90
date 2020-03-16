!> Basic model-dependend hard-coded run parameters.
module mod_parameters

  implicit none

  integer, parameter :: jorek_model    = 701       !< JOREK physics model
  
  parameter (var_u1         = 1)                       ! place of variable velocity 1  (u1)
  parameter (var_u2         = 2)                       ! place of variable velocity 2  (u2)
  parameter (var_u3         = 3)                       ! place of variable velocity 3  (u3)
  parameter (var_A1         = 4)                       ! place of variable mag pot  1  (A1)
  parameter (var_A2         = 5)                       ! place of variable mag pot  2  (A2)
  parameter (var_A3         = 6)                       ! place of variable psi/mag pot 3(A3)
  parameter (var_T          = 7)                       ! place of variable temperature (T)

  integer, parameter :: n_var          = 7         !< number of variables
  integer, parameter :: n_dim          = 2         !< number of dimensions
  integer, parameter :: n_order        = 3         !< order of the polynomial basis
  integer, parameter :: n_tor          = 1         !< number of toroidal harmonics
  integer, parameter :: n_period       = 1         !< periodicity in toroidal direction
  integer, parameter :: n_plane        = 4         !< number of toroidal angles
  integer, parameter :: n_vertex_max   = 4         !< maximum number of corners of an element
  integer, parameter :: n_nodes_max    = 20001     !< maximum number of nodes
  integer, parameter :: n_elements_max = 20001     !< maximum number of elements
  integer, parameter :: n_boundary_max = 1001      !< maximum number of boundary elements
  integer, parameter :: n_pieces_max   = 6001      !< maximum number of line pieces describing a flux surface
  integer, parameter :: n_degrees      = n_order+1 !< degrees of freedom per variable per node
  integer, parameter :: nref_max       = 10000     !< (refinement)
  integer, parameter :: n_ref_list     = 10000     !< (refinement)
  
  !> Names of the physical variables
  character(len=11) :: variable_names(n_var) =                       &
    (/ 'V_R        ','V_Z        ','V_phi      ','A_R        ',      &
       'A_Z        ','A_3        ','Temperature' /)

end module mod_parameters
