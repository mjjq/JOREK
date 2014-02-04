!> Demonstration of the diagnostic framework mod_position / mod_expression / mod_four_filter / mod_straight_field_line / ...
program demo
  
  use parameters
  use data_structure
  use phys_module
  use boundary
  use mod_new_diag
  use basis_at_gaussian
  
  implicit none
  
  type(type_node_list),         pointer :: node_list
  type(type_element_list),      pointer :: element_list
  type (type_bnd_element_list), pointer :: bnd_elm_list
  type (type_bnd_node_list),    pointer :: bnd_node_list
  type(t_equil_state)  :: equil_state
  type(t_pol_pos_list) :: pol_pos_list
  type(t_tor_pos_list) :: tor_pos_list
  type(t_four_filter)  :: filter
  type(t_expr_list)    :: expr_list
  integer :: my_id, ierr, k_tor, i, j, k, n(4)
  real*8, allocatable :: result(:,:,:,:), res0d(:), res1d(:,:), res2d(:,:,:), temp(:,:,:,:)
  
  allocate(node_list)
  allocate(element_list)
  allocate(bnd_elm_list)
  allocate(bnd_node_list)
  my_id = 0
  call initialise_parameters(my_id, "__NO_FILENAME__")
  
  call det_modes()
  
  call import_restart(node_list, element_list, 'jorek_restart.rst', rst_format, ierr)

  call initialise_basis                                       ! define the basis functions at the Gaussian points
  
  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
  call update_equil_state(node_list, element_list, bnd_elm_list, xpoint, xcase, equil_state)
  call print_equil_state(equil_state, .false.)
  
  
  
  ! --- Initialize expression module and print all currently available expressions
  call init_expr()
  call print_exprs(exprs_all)
  
  
  
  ! --- How to select expressions, some examples:
  expr_list = exprs_all
  expr_list = exprs_magfield
  expr_list = join_exprs(exprs_basicvar,exprs_magfield,exprs_cylcoord) !<<< does not work yet, to be done later
  expr_list = exprs((/'Psi ', 'B_R ', 'xjac', 'T   ', 'rho ', 'zj  '/), 6)
  
  
  
  ! --- Evaluate several expressions at one single position.
  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, equil_state, R=3., Z=0.1)
  call create_tor_pos(tor_pos_list, ierr, phi=0.)
  expr_list = exprs((/'B_R ', 'xjac', 'T   ', 'rho ', 'zj  '/), 5)
  call eval_expr(equil_state, JOREK_UNITS, expr_list, pol_pos_list, tor_pos_list, result, ierr)
  
  
  
  ! --- Print results in two different ways to the screen
  call reduce_result_to_0d(ierr, result, res0d, 1, 1, 1)
  call write_ascii_0d(ierr, equil_state, expr_list, res0d, FORM_TABLE, header=.true.)
  call write_ascii_0d(ierr, equil_state, expr_list, res0d, FORM_LIST)
  
  
  
  ! --- Or as a single call:
  call eval_expr(equil_state, JOREK_UNITS, exprs((/'B_R ', 'xjac', 'T   ', 'rho ', &
       'zj  '/), 5),          &
    pol_pos(node_list,element_list,equil_state,R=3.,Z=0.1), tor_pos(phi=0.), result, ierr)
  
  
  
  ! --- Evaluate several expressions on the outboard midplane and write to file.
  expr_list = exprs((/'R    ', 'Z    ', 'Psi_N', 'Psi  ', 'theta', 'x    ', 'y    ', &
       'phi  ', 'B_R  ', 'xjac ', 'T    ', 'rho  ', 'zj   ', 'omega', 'u    '/), 15)
  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, equil_state,                    &
      Rstart=equil_state%R_axis, Rend=equil_state%R_midpl(2)-1.d-3, Zstart=equil_state%Z_axis,     &
      Zend=equil_state%Z_axis, n=500)
  call eval_expr(equil_state, JOREK_UNITS, expr_list, pol_pos_list, tor_pos(phi=0.), result, ierr)
  call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)
  call write_ascii_1d(ierr, equil_state, expr_list, res1d, FORM_TABLE, header=.true.,              &
    filename='midplane_profiles.dat', append=.false., comment='Various profiles')
  
  
  
  ! --- Evaluate expressions along toroidal direction.
  call eval_expr(equil_state, JOREK_UNITS, expr_list,                                              &
    pol_pos(node_list,element_list,equil_state,R=3.,Z=0.1), tor_pos(nphi=128), result, ierr)
  
  
  
  ! --- Evaluate expressions on flux surfaces using straight field line theta.
  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, equil_state, PsiNmin=0.01,      &
    PsiNmax=0.99, nPsiN=16, nTht=64)
  call create_tor_pos(tor_pos_list, ierr, nphi=8)
  call eval_expr(equil_state, JOREK_UNITS, expr_list, pol_pos_list, tor_pos_list, result, ierr)
  n(:) = (/ size(result,1), size(result,2), size(result,3), size(result,4) /)
  
  
  
  ! --- Output 2D VTK files
  call reduce_result_to_2d(ierr, result, res2d, i1=1)
  call write_vtk_2d(ierr, expr_list, res2d, 'testA.vtk', (/1,2/), close1=.true.) ! vs R,Z
  res2d(:,:,5)=res2d(:,:,5)/(2.d0*PI)
  call write_vtk_2d(ierr, expr_list, res2d, 'testB.vtk', (/3,5/)) ! vs PsiN,theta
  call write_vtk_2d(ierr, expr_list, res2d, 'testC.vtk', (/6,2/)) ! vs x,z
  call reduce_result_to_2d(ierr, result, res2d, i3=16)
  call write_vtk_2d(ierr, expr_list, res2d, 'testC.vtk', (/5,8/)) ! vs theta,phi
  
  
  
  ! --- Check that 2D Fourier transform forward and backward recovers the original data
  allocate( temp(n(1),n(2),n(3),n(4)) )
  temp(:,:,:,:) = result(:,:,:,:)
  call perform_four_trafo(result, POLTOR_TRAFO, FORWARD_TRAFO)
  call perform_four_trafo(result, POLTOR_TRAFO, BACKWARD_TRAFO)
  write(*,'(1x,a,es20.10)') 'Relative difference after forward and backward transform: ',          &
    maxval(abs(result-temp)) / maxval(abs(temp))
  
  
  
  ! --- Check that 1D toroidal Fourier transform forward and backward recovers the original data
  call perform_four_trafo(result, TOROIDAL_TRAFO, FORWARD_TRAFO)
  call perform_four_trafo(result, TOROIDAL_TRAFO, BACKWARD_TRAFO)
  write(*,'(1x,a,es20.10)') 'Relative difference after forward and backward transform: ',          &
    maxval(abs(result-temp)) / maxval(abs(temp))
  
  
  
  ! --- Check that 1D poloidal Fourier transform forward and backward recovers the original data
  call perform_pol_trafo(result, FORWARD_TRAFO)
  call perform_pol_trafo(result, BACKWARD_TRAFO)
  write(*,'(1x,a,es20.10)') 'Relative difference after forward and backward transform: ',          &
    maxval(abs(result-temp)) / maxval(abs(temp))
  
  
  
  ! --- Apply Fourier filter: Keep harmonics (m,n) = (0...1,0...5) and (m,n)=(<any>,16)
  call perform_four_trafo(result, POLTOR_TRAFO, FORWARD_TRAFO)
  call init_four_filter(filter)
  call filter_add(filter, ierr, m_start=0, m_end=1, n_end=5)
  call filter_add(filter, ierr, n=16)
  call print_filter(filter)
  call apply_four_filter(result, filter)
  call perform_four_trafo(result, POLTOR_TRAFO, BACKWARD_TRAFO)
  
  
  
  ! --- Or simple filtering is possible as a single call:
  call transform_and_filter(result, simple_filter(m=0,n=0)) ! Keep only (m,n)=(0,0)
  
  
  
  ! --- Write result to a file = poloidally and toroidally averaged profiles
  call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)
  call write_ascii_1d(ierr, equil_state, expr_list, res1d, FORM_TABLE, header=.true.,              &
    filename='average_profiles.dat')
  
  
  
  ! --- Toroidally averaged expressions on the outboard midplane.
  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, equil_state,                    &
      Rstart=equil_state%R_axis+1d-3, Rend=equil_state%R_midpl(2)-1d-3, Zstart=equil_state%Z_axis, &
      Zend=equil_state%Z_axis, n=200)
  call eval_expr(equil_state, JOREK_UNITS, expr_list, pol_pos_list, tor_pos(nphi=16), result, ierr)
  call transform_and_filter(result, simple_filter(n=0))
  call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)
  call write_ascii_1d(ierr, equil_state, expr_list, res1d, FORM_TABLE, header=.true.,              &
    filename='toroidally_averaged_midplane_profiles.dat', append=.false.)
  
  
  
  ! --- Toroidally averaged midplane profiles (hfs and lfs)
  call midplane_profiles(node_list, element_list, equil_state, JOREK_UNITS, expr_list, res1d,      &
    LOWFIELD_SIDE, 200, ierr)
  call write_ascii_1d(ierr, equil_state, expr_list, res1d, FORM_TABLE, header=.true.,              &
    filename='lfs_profiles.dat', append=.false.)
  call midplane_profiles(node_list, element_list, equil_state, JOREK_UNITS, expr_list, res1d,      &
    HIGHFIELD_SIDE, 200, ierr)
  call write_ascii_1d(ierr, equil_state, expr_list, res1d, FORM_TABLE, header=.true.,              &
    filename='hfs_profiles.dat', append=.false.)
  
  
  
end program demo
