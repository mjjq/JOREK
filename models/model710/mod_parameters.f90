!> Basic model-dependend hard-coded run parameters.
module mod_parameters

  implicit none

  integer, parameter :: jorek_model    = 710       !< JOREK physics model

  integer, parameter :: var_A3 = 1                       ! place of variable psi/mag pot 3(A3)
  integer, parameter :: var_AR = 2                       ! place of variable mag pot  1  (AR)
  integer, parameter :: var_AZ = 3                       ! place of variable mag pot  2  (AZ)
  integer, parameter :: var_uR = 4                       ! place of variable velocity 1  (uR)
  integer, parameter :: var_uZ = 5                       ! place of variable velocity 2  (uZ)
  integer, parameter :: var_up = 6                       ! place of variable velocity 3  (up)
  integer, parameter :: var_r  = 7                       ! place of variable density     (r)
  integer, parameter :: var_T  = 8                       ! place of variable temperature (T)

  integer, parameter :: n_var          = 8         !< number of variables
  integer, parameter :: n_dim          = 2         !< number of dimensions
  integer, parameter :: n_order        = 3         !< order of the polynomial basis
  integer, parameter :: n_tor          = 1         !< number of toroidal harmonics
  integer, parameter :: n_period       = 1         !< periodicity in toroidal direction
  integer, parameter :: n_plane        = 1         !< number of toroidal angles
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
    (/ 'A_3        ','A_R        ','A_Z        ','u_R        ',      &
       'u_Z        ','u_phi      ','Density    ','Temperature' /)

end module mod_parameters
