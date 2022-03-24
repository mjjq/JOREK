!> Program to convert a JOREK restart file into an MHD IDS
program jorek2_MHD_IDS

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
  type (type_element_list),     pointer :: element_list
  
  integer    :: i, j, k, m, etype, irst, int, i_var, i_tor, k_tor, i_plane, index, index_node, my_id, ierr
  character  :: buffer*80, lf*1, str1*12, str2*12
  real*8     :: fact_T, fact_time, fact_v, fact_zj, fact_psi, rho0, fact_phi, fact_rho
  
  
  ! **********************************************************************************
  ! ******************************* IMAS **********************************************
  ! **********************************************************************************
  type(ids_mhd),                      target  :: mhd_ids
  type(ids_generic_grid_scalar),      pointer :: ggd_scalar
  type(ids_generic_grid_aos3_root),   pointer :: grid
  
  character(len=200):: user, database
  integer:: idx, shot_number, run_number, num_nodes
  
  integer :: n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub, n_grid
  ! **********************************************************************************


  namelist /imas_params/ shot_number, run_number, user, database 

  call flush_it(6)

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
  call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)

  ! --- Number of grids and grid subsets
  n_grid       = 1
  n_grid_sub   = 1
  grid_ind     = 1  ! Index
  grid_sub_ind = 1  ! Index

  ! --- Put the grid in GGD
  allocate( mhd_ids%grid_ggd(n_grid) )
  grid => mhd_ids%grid_ggd(grid_ind)
  call grid2ggd( grid, node_list, element_list )

  ! --- Number of time slices (restarts)
  n_slice = 1  

  mhd_ids%ids_properties%homogeneous_time = 1
  allocate(  mhd_ids%time(1)  )
  mhd_ids%time = 0.d0 
  allocate( mhd_ids%ggd(n_slice) )

  ! --- Normalization factors for IMAS
  rho0               = central_density * 1.d20 * central_mass * mass_proton
  sqrt_mu0_rho0      = sqrt( mu_zero * rho0 )
  sqrt_mu0_over_rho0 = sqrt( mu_zero / rho0 )

  fact_psi  = -2.d0 * PI                  ! Transform to COCOS convention 8 --> 11
  fact_time =  sqrt_mu0_rho0 
  fact_v    =  1.d0 /  sqrt_mu0_rho0 
  fact_phi  =  1.d0 /  sqrt_mu0_rho0 * F0 
  fact_zj   = -1.d0 / mu_zero * (-1.d0)   ! Last sign due to COCOS transformation
  fact_rho  =  rho0 
  fact_T    =  1.d0 / ( EL_CHG * mu_zero * central_density * 1.d20 )   


  ! --- Set times
  i_slice = 1
  mhd_ids%ggd(i_slice)%time = t_now * fact_time

  num_nodes = node_list%n_nodes

  ! --- Fill MHD data
  do i=1, n_var 

    ! --- Poloidal magnetic flux
    if (variable_names(i) == 'Psi') then      
      allocate( mhd_ids%ggd(i_slice)%psi(n_grid_sub))
      ggd_scalar => mhd_ids%ggd(i_slice)%psi(grid_sub_ind)
      call fill_Bezier_coefficients( ggd_scalar, node_list, var_psi, grid_ind, grid_sub_ind, fact_psi )
    endif

    ! --- Electrostatic potential 
    if (variable_names(i) == 'u') then      
      allocate( mhd_ids%ggd(i_slice)%phi_potential(n_grid_sub))
      ggd_scalar => mhd_ids%ggd(i_slice)%phi_potential(grid_sub_ind)
      call fill_Bezier_coefficients( ggd_scalar, node_list, var_u, grid_ind, grid_sub_ind, fact_phi )
    endif

    ! --- Toroidal current density. THIS IS NOT CORRECT!!!  A FACTOR 1/R IS MISSING!!!!!!!!!!
    if (variable_names(i) == 'zj') then      
      allocate( mhd_ids%ggd(i_slice)%j_tor(n_grid_sub))
      ggd_scalar => mhd_ids%ggd(i_slice)%j_tor(grid_sub_ind)
      call fill_Bezier_coefficients( ggd_scalar, node_list, var_zj, grid_ind, grid_sub_ind, fact_zj )
    endif

    ! --- Toroidal vorticity  THIS IS NOT CORRECT!!!  A FACTOR R IS MISSING!!!!!!!!!!
    if (variable_names(i) == 'omega') then      
      allocate( mhd_ids%ggd(i_slice)%vorticity(n_grid_sub))
      ggd_scalar => mhd_ids%ggd(i_slice)%vorticity(grid_sub_ind)
      call fill_Bezier_coefficients( ggd_scalar, node_list, var_w, grid_ind, grid_sub_ind, fact_v )
    endif

    ! --- Mass density
    if (variable_names(i) == 'rho') then      
      allocate( mhd_ids%ggd(i_slice)%mass_density(n_grid_sub))
      ggd_scalar => mhd_ids%ggd(i_slice)%mass_density(grid_sub_ind)
      call fill_Bezier_coefficients( ggd_scalar, node_list, var_rho, grid_ind, grid_sub_ind, fact_rho )
    endif

    ! --- Total temperature 
    if (variable_names(i) == 'T') then      
      ! --- Te
      allocate( mhd_ids%ggd(i_slice)%electrons%temperature(n_grid_sub))
      ggd_scalar => mhd_ids%ggd(i_slice)%electrons%temperature(grid_sub_ind)
      call fill_Bezier_coefficients( ggd_scalar, node_list, var_T, grid_ind, grid_sub_ind, fact_T*0.5d0 )

      ! --- Ti
      allocate( mhd_ids%ggd(i_slice)%t_i_average(n_grid_sub))
      ggd_scalar => mhd_ids%ggd(i_slice)%t_i_average(grid_sub_ind)
      call fill_Bezier_coefficients( ggd_scalar, node_list, var_T, grid_ind, grid_sub_ind, fact_T*0.5d0 )
    endif

    ! --- Ion temperature
    if (variable_names(i) == 'T_i') then      
      allocate( mhd_ids%ggd(i_slice)%t_i_average(n_grid_sub))
      ggd_scalar => mhd_ids%ggd(i_slice)%t_i_average(grid_sub_ind)
      call fill_Bezier_coefficients( ggd_scalar, node_list, var_Ti, grid_ind, grid_sub_ind, fact_T )
    endif

    ! --- Electron temperature
    if (variable_names(i) == 'T_e') then      
      allocate( mhd_ids%ggd(i_slice)%electrons%temperature(n_grid_sub))
      ggd_scalar => mhd_ids%ggd(i_slice)%electrons%temperature(grid_sub_ind)
      call fill_Bezier_coefficients( ggd_scalar, node_list, var_Te, grid_ind, grid_sub_ind, fact_T )
    endif

    ! --- Parallel velocity. NOT CORRECT!!!!!!! A FACTOR B_TOT IS MISSING!!!!!!!!!!
    if (variable_names(i) == 'v_par') then      
      allocate( mhd_ids%ggd(i_slice)%velocity_parallel(n_grid_sub))
      ggd_scalar => mhd_ids%ggd(i_slice)%velocity_parallel(grid_sub_ind)
      call fill_Bezier_coefficients( ggd_scalar, node_list, var_vpar, grid_ind, grid_sub_ind, fact_v )
    endif

  enddo

  ! --- Put data into local database
  call imas_create_env('ids',shot_number,run_number, 0,0,idx,user,database,'3') ! 3 is the database version  
  call ids_put(idx,'mhd',mhd_ids)
  call imas_close(idx)


end program jorek2_MHD_IDS
