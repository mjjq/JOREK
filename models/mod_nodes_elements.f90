!> Data structures representing the grid elements and nodes.
module nodes_elements
  
  use data_structure
  use parameters
  
  type (type_node_list)        :: node_list      !< List of grid nodes.
  type (type_element_list)     :: element_list   !< List of grid elements.
  type (type_bnd_element_list) :: bnd_elm_list   !< List of boundary elements.
  type (type_bnd_node_list)    :: bnd_node_list  !< List of boundary nodes.
  
end module nodes_elements
