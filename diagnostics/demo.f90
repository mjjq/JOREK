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
  
  
  
  ! --- Initialize expression module and print available expressions
  call init_expr()
  call print_expr()
  
  
  
  ! --- Evaluate several expressions simultaneously at one position.
  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, equil_state, R=10.5, Z=0.5)
  
  call create_tor_pos(tor_pos_list, ierr, phi=0.)
  
  call eval_expr(equil_state, .false., (/ 'R           ', 'Z           ', 'phi         ', 'xjac        ' /), 4,    &
    pol_pos_list, tor_pos_list, result, ierr)
  
  write(*,*) result
  
  
  
  ! --- Evaluate several expressions on the midplane.
  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, equil_state, Rstart=9.,         &
    Rend=11., Zstart=0., Zend=0., n=300)
  
  call create_tor_pos(tor_pos_list, ierr, phi=0.)
  
  call eval_expr(equil_state, .false., (/ 'R           ', 'Psi_N       ', 'u           ', 'T           ' /), 4, pol_pos_list, tor_pos_list, result, ierr)
  
  do i = 1, 300
    write(42,'(99es20.12)') result(:,1,i,1)
  end do
  
  
  
  ! --- Evaluate expressions along toroidal direction.
  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, equil_state, R=10.5, Z=0.5)
  
  call create_tor_pos(tor_pos_list, ierr, nphi=128)
  
  call eval_expr(equil_state, .false., (/ 'phi         ', 'T           ' /), 2, pol_pos_list, tor_pos_list, result, ierr)
  
  do i = 1, 128
    write(43,'(99es20.12)') result(:,1,1,i)
  end do
  
  
  
  ! --- Evaluate expressions on flux surfaces using straight field line theta.
  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, equil_state, PsiNmin=0.01, PsiNmax=0.99, nPsiN=5, nTht=16)
  call create_tor_pos(tor_pos_list, ierr, nphi=16)
  
  call eval_expr(equil_state, .false., (/ 'B_tor       ' /), 1, pol_pos_list, tor_pos_list, result, ierr)
  
  
  
  ! --- Check that Fourier transforming forward and backward recovers the original data
  n(:) = (/ size(result,1), size(result,2), size(result,3), size(result,4) /)
  allocate( temp(n(1),n(2),n(3),n(4)) )
  temp(:,:,:,:) = result(:,:,:,:)
  call perform_four_trafo(result, .true.)
  call perform_four_trafo(result, .false.)
  write(*,*) 'rel. diff=',maxval(abs(result-temp)) / maxval(abs(temp))
  
  
  
  ! --- Apply a simple Fourier filter and output some data before and afterwards
  do i = 1, n(2)
    write(37,*) result(1,i,3,1)
  end do
  
  !   --- Most general variant allowing to use complicated filters
  call perform_four_trafo(result, .true.)
  call init_four_filter(filter)
  call filter_add(filter, ierr, m_start=0, m_end=1, n_end=5) ! Keep harmonics (m,n) = (0...1,0...5)
  call print_filter(filter)
  call apply_four_filter(result, filter)
  call perform_four_trafo(result, .false.)
  do i = 1, n(2)
    write(38,*) result(1,i,3,1)
  end do
  
  !   --- Short form for simple filter
  call transform_and_filter(result, simple_filter(m=0,n=0))
  do i = 1, n(2)
    write(39,*) result(1,i,3,1)
  end do
  
end program demo
