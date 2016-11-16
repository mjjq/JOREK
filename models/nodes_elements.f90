!> Data structures representing the grid elements and nodes.
module nodes_elements
  
  use data_structure
  use mod_parameters
  
  type (type_node_list)        :: node_list      !< List of grid nodes.
  type (type_element_list)     :: element_list   !< List of grid elements.
  type (type_bnd_element_list) :: bnd_elm_list   !< List of boundary elements.
  type (type_bnd_node_list)    :: bnd_node_list  !< List of boundary nodes.

#ifdef IMPORT_PERTURBATIONS
  type (type_node_list)        :: node_list2      !< List of grid nodes.
  type (type_element_list)     :: element_list2   !< List of grid elements.
#endif  
  
end module nodes_elements
