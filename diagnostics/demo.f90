!> Demonstration of the diagnostic framework mod_position / mod_expression
program demo
  
  use parameters
  use data_structure
  use phys_module
  use mod_position
  use mod_expression
  use boundary
  
  implicit none
  
  type(type_node_list)   ,      pointer :: node_list
  type(type_element_list),      pointer :: element_list
  type (type_bnd_element_list), pointer :: bnd_elm_list
  type (type_bnd_node_list),    pointer :: bnd_node_list
  type(t_equil_state) :: equil_state
  type(t_pol_pos_list) :: pol_pos_list
  type(t_tor_pos_list) :: tor_pos_list
  integer :: my_id, ierr, k_tor, i
  real*8, allocatable :: result(:,:,:)
  
  allocate(node_list)
  allocate(element_list)
  allocate(bnd_elm_list)
  allocate(bnd_node_list)
  my_id = 0
  call initialise_parameters(my_id, "__NO_FILENAME__")
  
  call det_modes()
  
  call import_restart(node_list, element_list, 'jorek_restart.rst', rst_format, ierr)
  
  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
  
  call update_equil_state(node_list, element_list, bnd_elm_list, xpoint, xcase, equil_state)
  
  ! --- Evaluate several expressions simultaneously at one position.
  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, equil_state, R=10.5, Z=0.5)
  
  call create_tor_pos(tor_pos_list, ierr, phi=0.)
  
  call eval_expr(equil_state, .false., (/ EXPR_R, EXPR_Z, EXPR_PHI, EXPR_XJAC, EXPR_PSI, EXPR_U, &
    EXPR_ZJ, EXPR_W, EXPR_RHO, EXPR_T, EXPR_VPAR, EXPR_PRES, EXPR_BABS, EXPR_BTOR, EXPR_BR,      &
    EXPR_BZ, EXPR_CURRDENS /), 17,    &
    pol_pos_list, tor_pos_list, result, ierr)
  
  write(*,*) result
  
  ! --- Evaluate several expressions on the midplane.
  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, equil_state, Rstart=9.,         &
    Rend=11., Zstart=0., Zend=0., n=300)
  
  call create_tor_pos(tor_pos_list, ierr, phi=0.)
  
  call eval_expr(equil_state, .false., (/ EXPR_R, EXPR_PSIN, EXPR_ETAT, EXPR_BABS, EXPR_T,         &
    EXPR_RHO, EXPR_ZJ, EXPR_CURRDENS /), 8, pol_pos_list, tor_pos_list, result, ierr)
  
  do i = 1, 300
    write(42,'(99es20.12)') result(:,i,1)
  end do
  
  ! --- Evaluate expressions along toroidal direction.
  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, equil_state, R=10.5, Z=0.5)
  
  call create_tor_pos(tor_pos_list, ierr, nphi=128)
  
  call eval_expr(equil_state, .false., (/ EXPR_PHI, EXPR_T /), 2, pol_pos_list, tor_pos_list, result, ierr)
  
  do i = 1, 128
    write(43,'(99es20.12)') result(:,1,i)
  end do
  
end program demo
