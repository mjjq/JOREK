!> Demonstration of the diagnostic framework mod_position / mod_expression
program demo
  
  use parameters
  use data_structure
  use phys_module
  use mod_position
  use mod_expression
  use boundary
  use mod_four_filter
  
  implicit none
  
  type(type_node_list)   ,      pointer :: node_list
  type(type_element_list),      pointer :: element_list
  type (type_bnd_element_list), pointer :: bnd_elm_list
  type (type_bnd_node_list),    pointer :: bnd_node_list
  type(t_equil_state) :: equil_state
  type(t_pol_pos_list) :: pol_pos_list
  type(t_tor_pos_list) :: tor_pos_list
  type(t_four_filter)  :: filter
  type(t_expr_list)    :: expr_list
  integer :: my_id, ierr, k_tor, i, j, k, n(4)
  real*8, allocatable :: result(:,:,:,:), temp(:,:,:,:)
  
  
  real*8, allocatable :: test0(:,:,:), test1(:,:,:), test2(:,:,:)
  real*8 :: phi, tht
  integer*8 :: plan, rank, dims(2), howmany, inembed(2), istride, idist, onembed(2), ostride, odist
  integer*8, parameter :: N0=1, N1=33, N2=19
  
  
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
  call print_equil_state(equil_state, .false.)
  
  
  
  ! --- Initialize expression module and print all currently available expressions
  call init_expr()
  call print_expr(exprs_all)
  
  
  
  ! --- How to select expressions, some examples:
  expr_list = exprs_all
  expr_list = exprs_magfield
  expr_list = join_exprs(exprs_basicvar,exprs_magfield,exprs_cylcoord) !<<< does not work yet, to be done later
  expr_list = exprs((/'B_R', 'xjac', 'T', 'rho', 'zj'/), 5)
  
  
  
  ! --- Evaluate several expressions at one single position.
  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, equil_state, R=10.5, Z=0.5)
  call create_tor_pos(tor_pos_list, ierr, phi=0.)
  expr_list = exprs((/'B_R', 'xjac', 'T', 'rho', 'zj'/), 5)
  call eval_expr(equil_state, JOREK_UNITS, expr_list, pol_pos_list, tor_pos_list, result, ierr)
  
  write(*,*)
  write(*,'(1x,a,99f12.5)') 'Result at single point=', result
  write(*,*)
  
  
  
  ! --- Or as a single call:
  call eval_expr(equil_state, JOREK_UNITS, exprs((/'B_R', 'xjac', 'T', 'rho', 'zj'/), 5),          &
    pol_pos(node_list,element_list,equil_state,R=10.5,Z=0.5), tor_pos(phi=0.), result, ierr)
  
  
  
  ! --- Evaluate several expressions on the midplane.
  call eval_expr(equil_state, JOREK_UNITS, expr_list,                                              &
    pol_pos(node_list,element_list,equil_state,Rstart=9.,Rend=11.,Zstart=0.,Zend=0.,n=300),        &
    tor_pos(phi=0.), result, ierr)
  
  do i = 1, 300
    write(42,'(99es20.12)') result(:,1,i,1)
  end do
  
  
  
  ! --- Evaluate expressions along toroidal direction.
  call eval_expr(equil_state, JOREK_UNITS, expr_list,                                              &
    pol_pos(node_list,element_list,equil_state,R=10.5,Z=0.5), tor_pos(nphi=128), result, ierr)
  
  do i = 1, 128
    write(43,'(99es20.12)') result(:,1,1,i)
  end do
  
  
  
  ! --- Evaluate expressions on flux surfaces using straight field line theta.
  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, equil_state, PsiNmin=0.01,      &
    PsiNmax=0.99, nPsiN=5, nTht=16)
  call create_tor_pos(tor_pos_list, ierr, nphi=16)
  call eval_expr(equil_state, JOREK_UNITS, expr_list, pol_pos_list, tor_pos_list, result, ierr)
  
  
  
  ! --- Check that Fourier transforming forward and backward recovers the original data
  n(:) = (/ size(result,1), size(result,2), size(result,3), size(result,4) /)
  allocate( temp(n(1),n(2),n(3),n(4)) )
  temp(:,:,:,:) = result(:,:,:,:)
  call perform_four_trafo(result, .true.)
  call perform_four_trafo(result, .false.)
  write(*,'(1x,a,es20.10)') 'Relative difference after forward and backward transform: ',          &
    maxval(abs(result-temp)) / maxval(abs(temp))
  
  
  
  ! --- Apply Fourier filter: Keep harmonics (m,n) = (0...1,0...5) and (m,n)=(<any>,16)
  do i = 1, n(2)
    write(37,*) result(1,i,3,1)
  end do
  
  call perform_four_trafo(result, .true.)
  call init_four_filter(filter)
  call filter_add(filter, ierr, m_start=0, m_end=1, n_end=5)
  call filter_add(filter, ierr, n=16)
  call print_filter(filter)
  call apply_four_filter(result, filter)
  call perform_four_trafo(result, .false.)
  
  ! --- Or simple filtering is possible as a single call:
  call transform_and_filter(result, simple_filter(m=0,n=0)) ! Keep only (m,n)=(0,0)
  do i = 1, n(2)
    write(39,*) result(1,i,3,1)
  end do
  
end program demo
