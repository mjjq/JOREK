!> Program to convert a JOREK restart file into a radiation IDS
!> AT THE MOMENT ONLY WORKS WHEM THE RADIATION IS GIVEN IN THE AUX NODES
program jorek2_radiation_IDS

#ifdef USE_IMAS
  use ids_schemas !, only: ids_equilibrium
  use ids_routines, only: imas_open_env, &
     imas_create_env, imas_close, ids_get, ids_put
#endif

  use mod_parameters
  use data_structure
  use phys_module
  use constants
  use mpi_mod
  use mod_import_restart
  use mod_jorek2IMAS 
  
  implicit none
  
  type (type_node_list)   ,     pointer :: node_list
  type (type_node_list)   ,     pointer :: aux_node_list
  type (type_element_list),     pointer :: element_list
  
  integer    :: i, j, k, m, etype, irst, int, i_var, i_tor, k_tor, i_plane, index, index_node, my_id, ierr
  integer    :: var_rad
  character  :: buffer*80, lf*1, str1*12, str2*12
  real*8     :: fact_time, fact_rad, rho0
  
  ! **********************************************************************************
  ! ******************************* IMAS **********************************************
  ! **********************************************************************************
  type(ids_radiation),                      target  :: radiation_ids
  type(ids_generic_grid_scalar),      pointer :: ggd_scalar
  type(ids_generic_grid_aos3_root),   pointer :: grid
  
  character(len=200):: user, database
  integer:: idx, shot_number, run_number, num_nodes, stat
  
  integer :: n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub, n_grid
  ! **********************************************************************************


  namelist /imas_params/ shot_number, run_number, user, database 

  call flush_it(6)

  allocate(aux_node_list)
  allocate(node_list)
  allocate(element_list)

  ! --- Initialise input parameters and read the input namelist.
  my_id     = 0
  call initialise_parameters(my_id, "__NO_FILENAME__")
  
  ! --- Preset parameters
  database = 'test'
  shot_number=111112;   run_number=1
  call getenv('USER',user)
  
  ! --- Read parameters from namelist file 'vtk.nml' if it exists
  open(42, file='imas.nml', action='read', status='old', iostat=ierr)
  if ( ierr == 0 ) then
    write(*,*) 'Reading parameters from imas.nml namelist.'
    read(42, imas_params)
    close(42)
  end if

  call flush_it(6)

  ! --- Import JOREK restart file
  do k_tor=1, n_tor
    mode(k_tor) = + int(k_tor / 2) * n_period
  enddo

  ! --- Import radiation from the auxiliary nodes
  call import_hdf5_restart_aux(aux_node_list, 'aux_node_list_restart.h5', rst_format, ierr)

  ! --- Import grid
  call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)

  ! --- Number of grids and grid subsets
  n_grid       = 1
  n_grid_sub   = 1
  grid_ind     = 1  ! Index
  grid_sub_ind = 1  ! Index

  ! --- Put the grid in GGD
  allocate( radiation_ids%grid_ggd(n_grid) )
  grid => radiation_ids%grid_ggd(grid_ind)
  call grid2ggd( grid, node_list, element_list )

  ! --- Number of time slices (restarts)
  n_slice = 1  

  radiation_ids%ids_properties%homogeneous_time = 1
  allocate(  radiation_ids%time(n_slice)  )

  ! --- Normalization factors for IMAS
  rho0               = central_density * 1.d20 * central_mass * mass_proton
  sqrt_mu0_rho0      = sqrt( mu_zero * rho0 )
  sqrt_mu0_over_rho0 = sqrt( mu_zero / rho0 )

  fact_rad = 1.d0 / ( (gamma-1.d0) * MU_ZERO * sqrt_mu0_rho0 )
  fact_time =  sqrt_mu0_rho0 


  ! --- Some information of the radiation type
  allocate( radiation_ids%process(1))   ! --- 1 type of radiation
  allocate( radiation_ids%process(1)%identifier%name(1) )

  radiation_ids%process(1)%identifier%name  = "Line radiation"
  radiation_ids%process(1)%identifier%index = 10

  ! --- Set times
  i_slice = 1
  radiation_ids%time(i_slice) = t_now * fact_time
  allocate( radiation_ids%process(1)%ggd(n_slice) )

  ! --- Export radiation in Bezier coefficients
  var_rad = 2
  
  allocate( radiation_ids%process(1)%ggd(i_slice)%ion(1))
  allocate( radiation_ids%process(1)%ggd(i_slice)%ion(1)%emissivity(n_grid_sub))

  allocate( radiation_ids%process(1)%ggd(i_slice)%ion(1)%label(1) )  
  radiation_ids%process(1)%ggd(i_slice)%ion(1)%label = imp_type(index_main_imp) 

  ggd_scalar => radiation_ids%process(1)%ggd(i_slice)%ion(1)%emissivity(grid_sub_ind)
  call fill_Bezier_coefficients( ggd_scalar, aux_node_list, var_rad, grid_ind, grid_sub_ind, fact_rad )

  ! --- Put data into local database
  ! --- Try to open shot and number if it exists
  write(*,*) '  Adding radiation IDS to shot run = ', shot_number, run_number
  call imas_open_env( 'ids', shot_number,run_number,idx,user,database,'3',stat)! 3 is the database version  

  if (stat /= 0) then  ! --- Create a new shot if it doesn't exist
    write(*,*) '  Shot/run number did not exist, creating new one...'
    call imas_create_env('ids',shot_number,run_number, 0,0,idx,user,database,'3') 
  endif

  call ids_put(idx,'radiation',radiation_ids)
  call imas_close(idx)


end program jorek2_radiation_IDS
