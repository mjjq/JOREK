!> Program to convert a JOREK2 restart file into binary VTK format
program jorek2_stan

  use parameters, only: n_var, variable_names
  use data_structure
  use phys_module
  use basis_at_gaussian
  use nodes_elements
  use high_resolution_wall

  implicit none
  
  type (type_surface_list)	:: flux_list
  integer			:: i, k_tor, my_id, ierr, ifail
  integer			:: i_elm_xpoint(2), i_elm_axis
  real*8			:: R_xpoint(2), Z_xpoint(2), psi_xpoint(2), s_xpoint(2), t_xpoint(2)
  real*8			:: R_axis,	Z_axis,      psi_axis,      s_axis,	 t_axis
  real*8			:: psi_bnd
  integer			:: n_target
  integer, allocatable		:: index_target(:,:)
  real*8,  allocatable		:: R_target(:), Z_target(:)
  
  ! --- Initialise input parameters and read the input namelist.
  my_id     = 0
  call initialise_parameters(my_id, "__NO_FILENAME__")
  
  ! --- Define mode.
  do k_tor=1, n_tor
    mode(k_tor) = + int(k_tor / 2) * n_period
  enddo
  
  ! --- Import restart and define bases.
  call import_restart(node_list,element_list, 'jorek_restart.rst', rst_format, ierr)
  call initialise_basis
  
  ! --- Find axis
  call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
  
  ! --- Find Xpoint
  if (xpoint) then
    call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
    psi_bnd = psi_xpoint(1)
    if( (xcase .eq. 2) .or. ((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
      psi_bnd = psi_xpoint(2)
    endif
  else
    psi_bnd = 0.d0
  endif
  
  ! --- Define flux values
  flux_list%n_psi = 100
  allocate(flux_list%psi_values(flux_list%n_psi))
  do i=1,flux_list%n_psi-1
    flux_list%psi_values(i) = psi_axis + 1.7 * (psi_bnd - psi_axis) * real(i) / real(flux_list%n_psi-1)
  enddo
  flux_list%psi_values(flux_list%n_psi) = psi_bnd

  ! --- Find flux surfaces
  call find_flux_surfaces(xpoint,xcase,node_list,element_list,flux_list)
  do i=1,flux_list%n_psi
    flux_list%flux_surfaces(i)%psi = flux_list%psi_values(i)
  enddo
  
  ! --- Get rid of far private surfaces
  call clean_private(node_list,element_list,flux_list,psi_axis,psi_bnd,Z_xpoint)
  
  ! --- Order flux surfaces
  call reorder_flux_surfaces(node_list, element_list, flux_list, ierr)
  if (ierr .ne. 0) write(*,*)'Warning! reorder_flux_surfaces failed:',ierr
  
  ! --- Get high resolution wall
  call get_high_resolution_wall(node_list, element_list)
  
  ! --- Get flux surfaces that hit the target (assuming maximum 20 intersections per surface should be enough)
  allocate(R_target(20*flux_list%n_psi), Z_target(20*flux_list%n_psi), index_target(20*flux_list%n_psi,3))
  call get_target_flux_surfaces(node_list, element_list, flux_list, &
                                psi_bnd, R_axis, Z_axis, R_xpoint, Z_xpoint, &
				20*flux_list%n_psi, n_target, R_target, Z_target, index_target, ifail)
  
  
  
  
  deallocate(flux_list%psi_values)
  deallocate(R_target, Z_target, index_target)
  return
  
   
end program jorek2_stan





