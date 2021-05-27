!> Demonstration of the diagnostic framework mod_position / mod_expression / mod_four_filter / mod_straight_field_line / ...
program test_rhs_diagno 
  
  use mod_parameters
  use data_structure
  use phys_module
  use mod_boundary
  use mod_new_diag
  use basis_at_gaussian
  use mod_import_restart
  use equil_info
  
  implicit none
  
  type(type_node_list),         pointer :: node_list
  type(type_element_list),      pointer :: element_list
  type (type_bnd_element_list), pointer :: bnd_elm_list
  type (type_bnd_node_list),    pointer :: bnd_node_list
  type(t_pol_pos_list) :: pol_pos_list
  type(t_tor_pos_list) :: tor_pos_list
  type(t_four_filter)  :: filter
  type(t_expr_list)    :: expr_list
  integer :: my_id, ierr, k_tor, i, j, k, n(4)
  real*8, allocatable :: result(:,:,:,:), res0d(:), res1d(:,:), res2d(:,:,:)
  complex*16, allocatable :: cp(:,:,:,:)
  
  
  ! --- Normal initialization
  allocate(node_list)
  allocate(element_list)
  allocate(bnd_elm_list)
  allocate(bnd_node_list)
  my_id = 0
  call initialise_parameters(my_id, "__NO_FILENAME__")
  call det_modes()
  call import_restart(node_list, element_list, 'jorek_restart',  rst_format, ierr, .true.)
  call initialise_basis()
  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
  
  ! --- Initialize the plasma equilibrium data structure
  call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
  call print_equil_state(.false.)
  
  ! --- Initialize the new_diag framework and print some information (.true.)
  call init_new_diag(.true.)
  
  
  expr_list = exprs((/'R         ', 'Z         ', 'BR        ', &
    'BZ        ', 'Psi       ', 'T         '/), 6, 2)

  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, ES, grid=.true., nsub=7)
  call create_tor_pos(tor_pos_list, ierr, nphi=2)
  
  call eval_expr(ES, JOREK_UNITS, expr_list, pol_pos_list, tor_pos_list, result, ierr)
  
  call reduce_result_to_2d(ierr, result, res2d, i1=1)
  call write_vtk_2d(ierr, expr_list, res2d, 'test_all.vtk', (/1,2/), close1=.true.)
  
 
end program test_rhs_diagno
