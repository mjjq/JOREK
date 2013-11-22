!> Program to convert a JOREK2 restart file into binary VTK format
program jorek2_stan

  use parameters, only: n_var, variable_names
  use data_structure
  use phys_module
  use basis_at_gaussian
  use nodes_elements
  use high_resolution_wall
  use constants
  use tr_module 
  use grid_xpoint_data

  implicit none
  
  ! --- local variables
  type (type_surface_list) :: flux_list, sep_list

  type (type_strategic_points) , pointer     :: stpts
  type (type_new_points)       , pointer     :: nwpts

  real*8              :: psi_axis,      R_axis,      Z_axis,      s_axis,      t_axis
  real*8              :: psi_xpoint(2), R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2)
  integer             :: n_psi
  integer             :: i_elm_axis, i_elm_xpoint(2), i_elm_find(8), ifail
  integer             :: my_id
  real*8              :: psi_bnd, psi_bnd2
  real*8              :: sigmas(16)
  integer             :: n_grids(10)

  integer			:: i, k_tor, ierr
  integer			:: n_target
  integer, allocatable		:: index_target(:,:)
  real*8,  allocatable		:: R_target(:), Z_target(:)
  
  ! --- Initialise input parameters and read the input namelist.
  my_id = 0
  call initialise_parameters(my_id, "__NO_FILENAME__")
  
  ! --- Define mode.
  do k_tor=1, n_tor
    mode(k_tor) = + int(k_tor / 2) * n_period
  enddo
  
  ! --- Import restart and define bases.
  call import_restart(node_list,element_list, 'jorek_restart.rst', rst_format, ierr)
  call initialise_basis
  
  allocate(stpts)
  allocate(nwpts)
  call tr_register_mem(sizeof(stpts),"stpts",CAT_GRID)
  call tr_register_mem(sizeof(nwpts),"nwpts",CAT_GRID)
  !-------------------------------------------------------------------------------------------!
  !---------------------------- Initialise some internal data --------------------------------!
  !-------------------------------------------------------------------------------------------!
  
  !-------------------------------- Reset some parameters if they are inconsistent with XCASE
  if (xcase .eq. 1) then
    n_outer   = 0
    n_inner   = 0
    n_up_priv = 0
    n_up_leg  = 0
  endif
  if (xcase .eq. 2) then
    n_outer   = 0
    n_inner   = 0
    n_private = 0
    n_leg     = 0
  endif
  if ( (xcase .eq. 3) .and. (mod(n_tht,2) .ne. 0) )  n_tht = n_tht + 1
  if ( (xcase .ne. 3) .and. (mod(n_tht,2) .eq. 0) )  n_tht = n_tht + 1

  !-------------------------------- Build up some arrays to send as routine parameters (avoid long lists...)
  sigmas(1)  = SIG_closed  ; sigmas(2)  = SIG_theta
  sigmas(3)  = SIG_open    ; sigmas(4)  = SIG_outer   ; sigmas(5)  = SIG_inner
  sigmas(6)  = SIG_private ; sigmas(7)  = SIG_up_priv
  sigmas(8)  = SIG_leg_0   ; sigmas(9)  = SIG_leg_1
  sigmas(10) = SIG_up_leg_0; sigmas(11) = SIG_up_leg_1
  sigmas(12) = dPSI_open   ; sigmas(13) = dPSI_outer  ; sigmas(14) = dPSI_inner
  sigmas(15) = dPSI_private; sigmas(16) = dPSI_up_priv

  n_grids(1) = n_flux   ; n_grids(2) = n_tht
  n_grids(3) = n_open   ; n_grids(4) = n_outer  ; n_grids(5) = n_inner
  n_grids(6) = n_private; n_grids(7) = n_up_priv
  n_grids(8) = n_leg    ; n_grids(9) = n_up_leg
  n_grids(10)= 0 ! keep for n_tht_outer, which determines the angle of second Xpoint

  !-------------------------------------------------------------------------------------------!
  !----------------------------- Find MagAxis and Xpoint -------------------------------------!
  !-------------------------------------------------------------------------------------------!

  my_id = 1 ! Just don't want the printout...
  call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

  psi_bnd  = 0.d0
  psi_bnd2 = 0.d0
  call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
  if(xcase .eq. 1) psi_bnd = psi_xpoint(1)
  if(xcase .eq. 2) psi_bnd = psi_xpoint(2)
  if(xcase .eq. 3) then
    if(psi_xpoint(2) .lt. psi_xpoint(1)) then
      psi_bnd  = psi_xpoint(2)
      psi_bnd2 = psi_xpoint(1)
    else
      psi_bnd  = psi_xpoint(1)
      psi_bnd2 = psi_xpoint(2)  
    endif
    ! If we have a symmetric double-null, force the single separatrix
    if (abs(psi_xpoint(1)-psi_xpoint(2)) .lt. symmetric_threshold) then
      psi_xpoint(1) = (psi_xpoint(1)+psi_xpoint(2))/2.d0
      psi_xpoint(2) = psi_xpoint(1)
      psi_bnd  = psi_xpoint(1)
      psi_bnd2 = psi_bnd  
      n_grids(3) = 0
    endif
  endif




  !-------------------------------------------------------------------------------------------!
  !--------------- Define the flux values on which grid will be aligned ----------------------!
  !-------------------------------------------------------------------------------------------!

  !-------------------------------- Define number of psi values and allocate flux_list structure
  n_psi	        = n_flux   + n_open   + n_outer   + n_inner   + n_private   + n_up_priv + 1   ! this includes the magnetic axis
  flux_list%n_psi = n_psi - 1
  call tr_allocate(flux_list%psi_values,1,flux_list%n_psi,"flux_list%psi_values",CAT_GRID)

  !-------------------------------- Allocate sep_list structure (that's for plotting only)
  sep_list%n_psi =3
  if(xcase .eq. 3) sep_list%n_psi =6
  if (allocated(sep_list%psi_values)) call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
  call tr_allocate(sep_list%psi_values,1,sep_list%n_psi,"sep_list%psi_values",CAT_GRID)

  !-------------------------------- Call the routine
  call define_flux_values(node_list, element_list, flux_list, sep_list, &
                          xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_axis, n_grids, sigmas)

  !call plot_flux_surfaces(node_list,element_list,flux_list,.true.,1,psi_xpoint,R_xpoint,Z_xpoint,.true.,xcase)
  !call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,psi_xpoint,R_xpoint,Z_xpoint,.true.,xcase)

  if (allocated(sep_list%flux_surfaces))     deallocate(sep_list%flux_surfaces)

  ! --- Record psi-values of flux-surfaces individually
  do i=1,flux_list%n_psi
    flux_list%flux_surfaces(i)%psi = flux_list%psi_values(i)
  enddo
  
  ! --- Get rid of far private surfaces
  call clean_surfaces(node_list,element_list,flux_list,n_grids,psi_xpoint,R_xpoint,Z_xpoint)
  
  ! --- Order flux surfaces
  call reorder_flux_surfaces(node_list, element_list, flux_list, ierr)
  if (ierr .ne. 0) write(*,*)'Warning! reorder_flux_surfaces failed:',ierr
  
  ! --- Get high resolution wall
  call get_high_resolution_wall(node_list, element_list)
  
  !-------------------------------------------------------------------------------------------!
  !-------- Find all strategic points (leg corners, strike points and private middles) -------!
  !-------------------------------------------------------------------------------------------!

  !-------------------------------- Call the routine
  call find_strategic_points(node_list, element_list, flux_list, xcase, force_horizontal_Xline, &
                             R_xpoint, Z_xpoint, psi_xpoint, R_axis, Z_axis, n_grids, stpts)


  ! --- Get flux surfaces that hit the target (assuming maximum 20 intersections per surface should be enough)
  allocate(R_target(20*flux_list%n_psi), Z_target(20*flux_list%n_psi), index_target(20*flux_list%n_psi,3))
  call get_target_flux_surfaces(node_list, element_list, flux_list, stpts, &
                                psi_bnd, R_axis, Z_axis, R_xpoint, Z_xpoint, &
				20*flux_list%n_psi, n_target, R_target, Z_target, index_target, ifail)
  
  
  
  
  deallocate(flux_list%psi_values)
  deallocate(R_target, Z_target, index_target)
  
   
end program jorek2_stan





