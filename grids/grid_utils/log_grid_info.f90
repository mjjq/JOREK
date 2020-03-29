subroutine log_grid_info(verbose, node_list, element_list)

use data_structure

implicit none

! --- Routine parameters
logical,                 intent(in)  :: verbose
type(type_node_list),    intent(in)  :: node_list
type(type_element_list), intent(in)  :: element_list

! --- Local variables
integer, parameter :: maxbnd = 200
integer :: i, boundary_types(0:maxbnd), n_axis
real*8  :: Rmin, Rmax, Zmin, Zmax

write(*,*)
write(*,*) 'Number of grid nodes    / Max: ', node_list%n_nodes, size(node_list%node,1)
write(*,*) 'Number of grid elements / Max: ', element_list%n_elements, size(element_list%element,1)
write(*,*)
if ( verbose ) then
  write(*,*) 'Writing out grid node coordinates sorted by boundary type to fort.4?? files'
  write(*,*)
end if

boundary_types(:) = 0
Rmin = 1.d99
Rmax = -1.d99
Zmin = 1.d99
Zmax = -1.d99
n_axis = 0

do i = 1, node_list%n_nodes
  
  ! --- Detect how many nodes of each boundary type exist
  if ( (node_list%node(i)%boundary < 0) .or. (node_list%node(i)%boundary > maxbnd) ) then
    write(*,*) 'Skipping node ', i, ' with unusual boundary value ', node_list%node(i)%boundary
  else
    boundary_types(node_list%node(i)%boundary) = boundary_types(node_list%node(i)%boundary) + 1
    
    if ( verbose ) write(400+node_list%node(i)%boundary,*) node_list%node(i)%x(1,1:2)
  end if
  
  ! --- Determine geometrical region covered with nodes
  Rmin = min( Rmin, node_list%node(i)%x(1,1) )
  Rmax = max( Rmax, node_list%node(i)%x(1,1) )
  Zmin = min( Zmin, node_list%node(i)%x(1,2) )
  Zmax = max( Zmax, node_list%node(i)%x(1,2) )
  
  ! --- Count nodes on axis
  if ( node_list%node(i)%axis_node ) n_axis = n_axis + 1
end do

write(*,*) 'R range of grid nodes: ', Rmin, Rmax
write(*,*) 'Z range of grid nodes: ', Zmin, Zmax
write(*,*) 'Number of nodes per boundary type:'
do i = 0, maxbnd
  if (boundary_types(i)/=0) write(*,*) 'type', i,': ',boundary_types(i),' nodes'
end do
write(*,*)

end subroutine log_grid_info
