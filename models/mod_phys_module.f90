!> Input parameters and physical variables.
module phys_module
  
  use parameters
  use constants
  
  implicit none
  
  !> @name Various parameters
  real*8  :: eta                  !< Resistivity
  real*8  :: eta_T_0              !< Initial resistivity
  logical :: eta_T_dependent      !< Resistivity dependent on temperature? Otherwise constant.
  real*8  :: visco                !< Viscosity

  ! varivale needed by rst_bin2hdf5 and rst_hdf52bin
  real*8  :: visco_rst
  real*8  :: visco_par_rst
  real*8  :: eta_rst
  !

  real*8  :: visco2               !< Second coefficient of viscosity
  logical :: visco_T_dependent    !< Viscosity dependent on temperature? Otherwise constant.
  real*8  :: visco_par            !< Parallel viscosity
  real*8  :: F0                   !< Determines fixed toroidal magnetic field: \f$ B_\phi = F_0/R \f$
  real*8  :: central_density      !< particle density at the magnetic axis (in units of \f$10^{20} m^{-3}\f$)
  real*8  :: central_mass         !< average mass (assumed to be constant in space for the moment)
  real*8  :: gamma                !< ratio of specific heat (=5/3)
  real*8  :: Q_bar                !< (model400)
  real*8  :: sigma                !< (model400)
  real*8  :: tauIC                !< (diamagnetic terms)
  logical :: Wdia                 !< (diamagnetic vorticity)
  real*8  :: gamma_sheath         !< sheath boundary condition open fieldlines (model303)
  real*8  :: density_reflection   !< density reflection coeeficient open fieldlines (model303)
  integer :: mode(n_tor)          !< Toroidal mode number corresponding to the JOREK modes, e.g., for n_period=8 and n_tor=3, mode(:)=0,8,8
  integer :: nout                 !< Output a restart file every nout timesteps.
  integer :: xcase                !< DoubleNull: 1->LowerXpoint. 2->UpperXpoint. 3->doubleNull.
  integer :: rst_format           !< 0 == olf format, 1 == new format for restart file.
  logical :: restart              !< Restart a code run from the restart file jorek_restart.rst?
  logical :: regrid               !< Re-generate the flux-aligned grid (does not work currently)?
  logical :: import_equil         
  logical :: xpoint               !< X-point geometry?
  logical :: bootstrap            !< Bootstrap-current?
  logical :: refinement           !< Use mesh refinement?
  logical :: bc_natural_flux      !< boundary conditions for flux surface boundaries (2 and 3)
  logical :: bc_natural_open      !< use natural boundary conditions on the open fieldlines
  logical :: produce_live_data    !< Write data to 'energies.dat', 'growth_rates.dat', and 'times.dat' during the code run?
  logical :: grid_to_wall         !< extend the grid to a physical wall
  logical :: adaptive_time        !< requires no_mpi for Pastix library
  logical :: equil                !< compute equilibrium
  logical :: bench_without_plot   !< .true. for benchmark (mesuring elapsed time without plot phases) 
  logical :: gmres                !< Use iterative GMRES solver
  integer :: gmres_max_iter       !< Maximum number of GMRES iterations
  logical :: linear_run           !< Perform a linear run where the equilibrium quantities (i_tor=1) do not change with time?
  logical :: export_for_nemec     !< Export data such that the NEMEC Code can reconstruct the same equilibrium?

#ifdef USE_HDF5
  ! for HDF5 diagnostics
  logical :: save_diagnostics_HDF5  !< Export data in HDF5 format
#endif
  ! Number of digits that ends the filenames of diagnostics
  integer             :: nbdigits   = 5
  character(20)       :: numfmt     = "'_d',i5.5"
  character(20)       :: numfmt_rst = "'_r',i3.3"
  ! Identity of the processor
  integer  :: pglobal_id
  
  real*8, allocatable :: energies(:,:,:)  !< Magnetic and kinetic mode energies at timesteps.
  real*8, allocatable :: energies2(:,:,:)  !< global density and temperature at timesteps.
  real*8, allocatable :: energies3(:,:,:)  !< global currents (general and total eccd) at timesteps.
  real*8, allocatable :: energies4(:,:,:)  !< global applied eccd currents j1 and j2 at timesteps.

  character(len=3)    :: mode_type(n_tor) !< 'cos' or 'sin'
  
  !> Points used as limiters (see routine find_limiter)
  integer, parameter :: max_limiter = 1000 !< Maximum number of limiter points
  integer :: n_limiter                     !< Number of limiter points
  real*8  :: R_limiter(max_limiter)        !< R-positions of the limiter points
  real*8  :: Z_limiter(max_limiter)        !< Z-positions of the limiter points
  integer :: first_target_point		   !< index of the first target point on the limiter (for xpoint_grid_wall)
  integer :: last_target_point		   !< index of the last  target point on the limiter (does NOT need to be > first_target_point)
  
  !> @name Define X-point geometry by geometrical properties
  !!
  !! \f[
  !! \Psi(\theta) =
  !!        -x_{shift}\sin(\theta)
  !!        +x_{left}\cos(\theta)
  !!        +x_{ampl}\left[
  !!            \left(\frac{x_{width}\cdot(\theta-x_{theta})}{x_{sig}}\right)^2-1
  !!          \right]exp\left[-\left(\frac{\theta-x_{theta}}{x_{sig}}\right)^2\right]
  !! \f]
  real*8  :: xampl, xwidth, xsig, xtheta, xshift, xleft
  
  !> @name Heat and particle sources
  !!
  !! \f[
  !! S(\Psi_N) = S_0 \cdot \left[0.5 - 0.5 \tanh\left(\frac{\Psi_N - \Psi_{N,0}}{\sigma}\right) \right]
  !! \f]
  !!
  !! The following parameters can be set via the namelist input file:
  !! - \f$ S_0 \f$ denotes the source strenght (e.g., heatsource)
  !! - \f$ \Psi_{N,0} \f$ denotes the position around which the source is ramped down (e.g., heatsource_psin)
  !! - \f$ \sigma \f$ denotes the width over which the source is ramped down (e.g., heatsource_sig)
  !!
  real*8  :: particlesource      !< Particle source strength
  real*8  :: particlesource_psin !< Position around which source is ramped down
  real*8  :: particlesource_sig  !< Width over which source is ramped down
  real*8  :: heatsource          !< Heat source strength
  real*8  :: heatsource_psin     !< Position around which source is ramped down
  real*8  :: heatsource_sig      !< Width over which source is ramped down
  real*8  :: heatsource_i        !< Heat source strength (ions), model4xx only
  real*8  :: heatsource_e        !< Heat source strength (electrons), model4xx only
  
  !> @name Hyper-resistivity, -viscosity and -diffusivities
  real*8  :: eta_num, visco_num, visco_par_num, D_perp_num, Zk_perp_num
  
  !> @name Timestepping parameters
  real*8  :: tstep             		!< Size of the timesteps (\f$ \Delta t \f$)
  real*8  :: tstep_n(10)       		!< Alternative to tstep: Up to ten values may be given
  integer :: nstep             		!< Number of timesteps to perform
  integer :: nstep_n(10)       		!< Alternative to nstep: Up to ten values may be given
  real*8  :: t_start           		!< Time value at the start of the code run (zero or from restart file)
  real*8  :: t_now             		!< Current time value in the simulation
  integer :: index_start       		!< Time step index at the beginning of the code run (zero or from restart file)
  integer :: index_now         		!< Current time step index
  real*8, allocatable :: xtime(:) 	!< Time values corresponding to the timesteps.
  character(len=80) :: time_evol_scheme !< Time evolution scheme to use. This determines the values
                               		!! used for time_evol_theta and time_evol_zeta.
  real*8  :: time_evol_theta   		!< Time evolution parameter theta (see documentation)
  real*8  :: time_evol_zeta    		!< Time evolution parameter zeta (see documentation)

  integer :: rst_hdf5

#ifdef USE_HDF5
  real*8  :: h5_diag_nbtime    		!< the HDF5 diagnostics are saved every "h5_diag_nbtime" Alven times
  integer :: h5_nbsave_all     		!< number of HDF5 files written [or # of times the HDF5 saving has been called]
#endif
  
  !> @name Machine name
  character(len=512) :: tokamak_device 	!< Name of the tokamak device we are simulating

  !> @name Analytical boundary of initial grid
  !!
  !! Analytical definition of the boundary of the non flux-aligned initial polar grid.
  !!
  !! - \f$ Z=Z_{geo} + a_{min} \epsilon \sin(\theta) \f$
  !!
  !! - for \f$ \theta < \pi \f$:
  !!   \f$ R=R_{geo} + a_{min} \cos\left[\theta+T_u\sin(\theta)+Q_u\sin(2\theta)\right] \f$
  !!
  !! - for \f$ \theta \ge \pi \f$:
  !!   \f$ R=R_{geo} + a_{min} \cos\left[\theta+T_l\sin(\theta)+Q_l\sin(2\theta)\right] \f$
  !!
  real*8  :: amin              !< Minor radius
  real*8  :: ellip             !< Ellipticity
  real*8  :: tria_u            !< Upper triangularity
  real*8  :: tria_l            !< Lower triangularity
  real*8  :: quad_u            !< Upper quadrangularity
  real*8  :: quad_l            !< Lower quadrangularity
  
  !> @name Fourier expanded boundary of initial grid
  !! Boundary of the non flux-aligned initial polar grid given as Fourier series
  integer, parameter :: n_bnd_max = 3000 	!< Max number of entries in boundary points
  integer :: mf              		 	!< Number of entries in fbnd and fpsi
  real*8  :: fbnd(n_bnd_max)        		!< Fourier expansion of boundary
  real*8  :: fpsi(n_bnd_max)        		!< Fourier expansion of the poloidal flux at the boundary
  
  !> @name Numerical boundary of initial grid
  !! Numerical definition of the boundary of the non flux-aligned initial polar grid.
  integer :: n_boundary       			!< Number of points in R_boundary, Z_boundary, psi_boundary.
  real*8  :: R_boundary  (n_bnd_max)		!< Numerical R values defining the boundary
  real*8  :: Z_boundary  (n_bnd_max)		!< Numerical Z values defining the boundary
  real*8  :: psi_boundary(n_bnd_max)		!< Numerical values giving the poloidal flux at the boundary
  
  !> @name PF coils definition for initial equilibrium (MAST)
  !! Numerical definition of the PF coils definition for initial equilibrium (MAST)
  integer :: n_pfc            !< Number of coils
  real*8  :: Rmin_pfc(20)     !< Minimum R of coil
  real*8  :: Rmax_pfc(20)     !< Maximum R of coil
  real*8  :: Zmin_pfc(20)     !< Minimum Z of coil
  real*8  :: Zmax_pfc(20)     !< Maximum Z of coil
  real*8  :: current_pfc(20)  !< Current density in the coil
  
  !> @name Pellet-related input parameters
  real*8  :: pellet_amplitude  !< amplitude of density source (when pellet modelled as density source)
  real*8  :: pellet_R          !< major radius position pellet
  real*8  :: pellet_Z          !< Z position pellet
  real*8  :: pellet_phi        !< width of the pellet cloud (density source) in toroidal angle
  real*8  :: pellet_ellipse    !< the ellipticity of the pellet source
  real*8  :: pellet_radius     !< radius of the simulation pellet
  real*8  :: pellet_sig        !< width of smoothing of density source (arctan((r-pellet_radius)/pellet_sig))
  real*8  :: pellet_length     !< width of smoothing of density source in toroidal angle
  real*8  :: pellet_theta      !< orientation of the pellet ellipse
  real*8  :: pellet_psi        !< pellet_width in poloidal flux
  real*8  :: pellet_delta_psi  !< width of smoothing in poloidal flux
  real*8  :: pellet_velocity_R !< pellet velocity component radial direction
  real*8  :: pellet_velocity_Z !< pellet velocity component Z direction
  real*8  :: pellet_density    !< pellet density (in units 10^20 m^-3)
  real*8  :: pellet_particles  !< the number of particles in the pellet (in units of 10^20)
  logical :: use_pellet
  
  !> @name Massive gas injection-related input parameters
  real*8  :: mgi_amplitude      !< amplitude of neutral density source
  real*8  :: mgi_R             !< major radius position of neutral density source
  real*8  :: mgi_Z             !< Z position of neutral density source
  real*8  :: mgi_phi           !< width of the neutral density source in toroidal direction
  real*8  :: mgi_radius        !< radius of the neutral density source in poloidal plane
  real*8  :: mgi_sig           !< width of smoothing of the neutral density source in poloidal plane
  real*8  :: mgi_deltaphi      !< width of smoothing of the neutral density source in toroidal direction
  real*8  :: ksi_ion           !< energy cost of each ionization
  real*8  :: A_Dmv             !< Cross sectional area of DMV (Disruption mitigation valve) pipe
  real*8  :: K_Dmv             !< Correction parameter describing the gas expansion near the pipe orifice
  real*8  :: L_tube            !< Pipe length
  real*8  :: V_Dmv             !< Volume of the DMV reservoir
  real*8  :: P_Dmv             !< Pressure in the DMV reservoir
  real*8  :: t_mgi             !< Beginning of the MGI
  real*8  :: delta_n_convection!< Switch to activate the convection term for neutrals (at the plasma velocity)
  logical :: JET_MGI !< Switch to have a real time dependent MGI or a constant injection
  logical :: ASDEX_MGI
  
  !> @name Free boundary extension
  !! Input parameters related to the free boundary extension (folder vacuum/).
  logical :: freeboundary_equil!< use a free or fixed boundary equilibrium?
  logical :: freeboundary      !< use free or fixed boundary conditions in time-evolution?
  logical :: resistive_wall    !< use a resistive or ideal wall?    (free boundary only)
  
  !> @name Rectangular Grid
  !! Parameters defining a rectangular grid in R- and Z-directions in the poloidal plane.
  integer :: n_R               !< Number of grid points in R-direction
  integer :: n_Z               !< Number of grid points in Z-direction
  real*8  :: R_begin           !< Left boundary of grid in R-direction
  real*8  :: R_end             !< Right boundary of grid in R-direction
  real*8  :: Z_begin           !< Lower boundary of grid in Z-direction
  real*8  :: Z_end             !< Upper boundary of grid in Z-direction
  
  !> @name Polar Grid
  !! Parameters defining a non flux-aligned polar grid in the poloidal plane.
  logical :: force_horizontal_Xline     !< Force the grid line through Xpoint to be horizontal
  					!  (instead of perpendicular to line between Xpoint and Axis)
  integer :: n_radial          		!< Number of radial grid points
  integer :: n_pol             		!< Number of poloidal grid points
  real*8  :: R_geo             		!< Center of the grid
  real*8  :: Z_geo             		!< Center of the grid
  real*8  :: psi_axis_init     		!< Initial guess for Psi at the magnetic axis
  real*8  :: XR_r(2)           		!< Psi_N position of radial grid accumulation (two positions)
  real*8  :: SIG_r(2)          		!< Width of grid accumulation (two positions)
  real*8  :: XR_tht(2)         		!< Position of poloidal grid accumulation (0...1, two positions)
  real*8  :: SIG_tht(2)        		!< Width of grid accumulation (two positions)
  
  !> @name Flux surface grid
  !! Parameters defining a flux-aligned grid without X-point in the poloidal plane.
  integer :: n_flux            !< Number of radial grid points
  integer :: n_tht             !< Number of poloidal grid points
  real*8  :: xr1               !< Grid accumulation parameter
  real*8  :: xr2               !< Grid accumulation parameter
  real*8  :: sig1              !< Grid accumulation parameter
  real*8  :: sig2              !< Grid accumulation parameter
  
  !> @name Flux surface grid with X-point
  !! Parameters defining a flux-aligned grid with X-point in the poloidal plane.
  integer :: n_open            !< Number of 'radial' grid points in the open flux region - in between the two separatrices in case of double-null
  integer :: n_outer           !< Number of 'radial' grid points in the open flux region on the outer side (LFS) in case of double-null
  integer :: n_inner           !< Number of 'radial' grid points in the open flux region on the inner side (HFS) in case of double-null
  integer :: n_private         !< Number of 'radial' grid points in the private flux region on the lower side
  integer :: n_leg             !< Number of 'poloidal' grid points along the divertor legs on the lower side
  integer :: n_up_priv         !< Number of 'radial' grid points in the private flux region on the upper side (upper Xpoint or double-null)
  integer :: n_up_leg          !< Number of 'poloidal' grid points along the divertor legs on the upper side (upper Xpoint or double-null)
  integer :: n_ext             !< Number of 'radial' grid points from the outermost flux surface to wall)
  real*8  :: SIG_closed        !< Width with grid accumulation
  real*8  :: SIG_open          !< Width with grid accumulation
  real*8  :: SIG_outer         !< Width with grid accumulation
  real*8  :: SIG_inner         !< Width with grid accumulation
  real*8  :: SIG_private       !< Width with grid accumulation
  real*8  :: SIG_up_priv       !< Width with grid accumulation
  real*8  :: SIG_theta         !< Width with grid accumulation
  real*8  :: SIG_leg_0         !< Width with grid accumulation
  real*8  :: SIG_leg_1         !< Width with grid accumulation
  real*8  :: SIG_up_leg_0      !< Width with grid accumulation
  real*8  :: SIG_up_leg_1      !< Width with grid accumulation
  real*8  :: dPSI_open         !< Delta Psi grid extends into the open flux region
  real*8  :: dPSI_outer        !< Delta Psi grid extends into the open flux region
  real*8  :: dPSI_inner        !< Delta Psi grid extends into the open flux region
  real*8  :: dPSI_private      !< Delta Psi grid extends into the private flux region
  real*8  :: dPSI_up_priv      !< Delta Psi grid extends into the private flux region
  
  !> @name Analytical heat, particle and neutral particles diffusivity parameters
  real*8  :: D_perp(10), D_par
  real*8  :: ZK_perp(10), ZK_par, ZK_par_max, ZK_i_perp(10), ZK_e_perp(10), K_i_par, K_e_par
  real*8  :: D_neutral_x, D_neutral_y, D_neutral_p
  logical :: ZKpar_T_dependent

  !> @name Numerical heat and particle diffusivity profiles
  character(len=512)  :: d_perp_file        !< ASCII file the profile is read from
  character(len=512)  :: zk_perp_file       !< ASCII file the profile is read from
  character(len=512)  :: zk_e_perp_file     !< ASCII file the profile is read from (model400)
  character(len=512)  :: zk_i_perp_file     !< ASCII file the profile is read from (model400)
  logical             :: num_d_perp         !< is set true if d_perp_file /= 'none'
  logical             :: num_zk_perp        !< is set true if zk_perp_file /= 'none'
  logical             :: num_zk_e_perp      !< is set true if zk_e_perp_file /= 'none' (model400)
  logical             :: num_zk_i_perp      !< is set true if zk_i_perp_file /= 'none' (model400)
  integer             :: num_d_perp_len     !< Number of datapoints in the profile
  integer             :: num_zk_perp_len    !< Number of datapoints in the profile
  integer             :: num_zk_e_perp_len  !< Number of datapoints in the profile (model400)
  integer             :: num_zk_i_perp_len  !< Number of datapoints in the profile (model400)
  real*8, allocatable :: num_d_perp_x(:)    !< Psi_N values of the profile
  real*8, allocatable :: num_d_perp_y(:)    !< D_perp values of the profile
  real*8, allocatable :: num_zk_perp_x(:)   !< Psi_N values of the profile
  real*8, allocatable :: num_zk_perp_y(:)   !< ZK_perp values of the profile
  real*8, allocatable :: num_zk_e_perp_x(:) !< Psi_N values of the profile (model400)
  real*8, allocatable :: num_zk_e_perp_y(:) !< ZK_perp values of the profile (model400)
  real*8, allocatable :: num_zk_i_perp_x(:) !< Psi_N values of the profile (model400)
  real*8, allocatable :: num_zk_i_perp_y(:) !< ZK_perp values of the profile (model400)
  
  !> @name Analytical input profile for the density
  real*8  :: rho_0, rho_1,  rho_coef(10)
  
  !> @name Numerical input profile for the density
  character(len=512)  :: rho_file        !< ASCII file the profile is read from.
  logical             :: num_rho         !< is set true if rho_file /= 'none'
  integer             :: num_rho_len     !< Number of points in profile
  real*8, allocatable :: num_rho_x(:)    !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_rho_y0(:)   !< Values of density profile
  real*8, allocatable :: num_rho_y1(:)   !< First derivatives of density profile (\f$ d\rho/d\Psi_N \f$)
  real*8, allocatable :: num_rho_y2(:)   !< Second derivatives of density profile (\f$ d^2\rho/d\Psi_N^2 \f$)
  real*8, allocatable :: num_rho_y3(:)   !< Third derivatives of density profile (\f$ d^3\rho/d\Psi_N^3 \f$)

  !> @name Analytical input profile for the temperature
  real*8  :: T_0,   T_1,    T_coef(10)
  real*8  :: Ti_0,  Ti_1,   Ti_coef(10)
  real*8  :: Te_0,  Te_1,   Te_coef(10)
  
  !> @name Numerical input profile for the temperature
  character(len=512)  :: T_file          !< ASCII file the profile is read from.
  logical             :: num_T           !< is set true if T_file /= 'none'
  integer             :: num_T_len       !< Number of points in profile
  real*8, allocatable :: num_T_x(:)      !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_T_y0(:)     !< Values of temperature profile
  real*8, allocatable :: num_T_y1(:)     !< First derivatives of temperature profile (\f$ dT/d\Psi_N \f$)
  real*8, allocatable :: num_T_y2(:)     !< Second derivatives of temperature profile (\f$ d^2T/d\Psi_N^2 \f$)
  real*8, allocatable :: num_T_y3(:)     !< Third derivatives of temperature profile (\f$ d^3T/d\Psi_N^3 \f$)
  
  !> @name Numerical input profile for the ion temperature (model400)
  character(len=512)  :: Ti_file         !< ASCII file the profile is read from.
  logical             :: num_Ti          !< is set true if T_file /= 'none'
  integer             :: num_Ti_len      !< Number of points in profile
  real*8, allocatable :: num_Ti_x(:)     !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_Ti_y0(:)    !< Values of temperature profile
  real*8, allocatable :: num_Ti_y1(:)    !< First derivatives of temperature profile (\f$ dT/d\Psi_N \f$)
  real*8, allocatable :: num_Ti_y2(:)    !< Second derivatives of temperature profile (\f$ d^2T/d\Psi_N^2 \f$)
  real*8, allocatable :: num_Ti_y3(:)    !< Third derivatives of temperature profile (\f$ d^3T/d\Psi_N^3 \f$)
  
  !> @name Numerical input profile for the electron temperature (model400)
  character(len=512)  :: Te_file         !< ASCII file the profile is read from.
  logical             :: num_Te          !< is set true if T_file /= 'none'
  integer             :: num_Te_len      !< Number of points in profile
  real*8, allocatable :: num_Te_x(:)     !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_Te_y0(:)    !< Values of temperature profile
  real*8, allocatable :: num_Te_y1(:)    !< First derivatives of temperature profile (\f$ dT/d\Psi_N \f$)
  real*8, allocatable :: num_Te_y2(:)    !< Second derivatives of temperature profile (\f$ d^2T/d\Psi_N^2 \f$)
  real*8, allocatable :: num_Te_y3(:)    !< Third derivatives of temperature profile (\f$ d^3T/d\Psi_N^3 \f$)  
  
  !> @name Analytical input profile for the neutral density (model 500)
  real*8  :: rhon_0, rhon_1,  rhon_coef(10)
  
  !> @name Numerical input profile for the neutral density (model 500)
  character(len=512)  :: rhon_file        !< ASCII file the profile is read from.
  logical             :: num_rhon         !< is set true if rho_file /= 'none'
  integer             :: num_rhon_len     !< Number of points in profile
  real*8, allocatable :: num_rhon_x(:)    !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_rhon_y0(:)   !< Values of neutral density profile
  real*8, allocatable :: num_rhon_y1(:)   !< First derivatives of neutral density profile (\f$ d\rhon/d\Psi_N \f$)
  real*8, allocatable :: num_rhon_y2(:)   !< Second derivatives of neutral density profile (\f$ d^2\rhon/d\Psi_N^2 \f$)
  real*8, allocatable :: num_rhon_y3(:)   !< Third derivatives of neutral density profile (\f$ d^3\rhon/d\Psi_N^3 \f$)
  
  !> @name Analytical input profile for FFprime
  real*8  :: FF_0,  FF_1,   FF_coef(10)
  
  !> @name Numerical input profile for FFprime
  character(len=512)  :: ffprime_file    !< ASCII file the profile is read from.
  logical             :: num_ffprime     !< is set true if ffprime_file /= 'none'
  integer             :: num_ffprime_len !< Number of points in profile
  real*8, allocatable :: num_ffprime_x(:)!< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_ffprime_y0(:) !< Values of FFprime profile
  real*8, allocatable :: num_ffprime_y1(:) !< First derivatives of FFprime profile (\f$ dFF'/d\Psi_N \f$)
  real*8, allocatable :: num_ffprime_y2(:) !< Second derivatives of FFprime profile (\f$ d^2FF'/d\Psi_N^2 \f$)

  !> --- Numerical input profiles for neoclassical coefficients 
  logical             :: NEO         ! If (NEO == .t.) neoclassical effects are considered
  character(len=512)  :: neo_file    ! ASCII file the aki and amu profiles is read from.
  logical             :: num_neo_file    ! is set true if neo_file /= 'none'
  integer             :: num_neo_len     ! Number of points in aki_neo, mu_neo profiles
  real*8, allocatable :: num_neo_psi(:)  ! Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_aki_value(:)  !numerical aki profile (PsiN values)
  real*8, allocatable :: num_amu_value(:)!numerical amu profile (PsiN values)
  real*8  :: aki_neo_const, amu_neo_const !if ( (NEO) .and. (neo_file=='none')), aki_neo and amu_neo are constants

  !> @name RMP profiles
  logical :: output_bnd_elements !< If .true., writes bnd nodes and bnd elements in files 'boundary_nodes.dat' and 'boundary_elements.dat' 
  logical :: RMP_on            !< Activates RMPs on boundary if .true.
  character(len=512)  :: RMP_psi_cos_file  !< ASCII file the profiles of psi_RMP_cos and derivatives are read from
  character(len=512)  :: RMP_psi_sin_file  !< ASCII file the profiles of psi_RMP_sin and derivatives are read from
  real*8  :: lambda, tset      !< parameters for time dependence of psi_RMP (sigmoid)
  real*8  :: RMP_start_time    !< time when RMP coils have been activated (RMP_on = .t.)
  real*8, allocatable :: psi_RMP_cos(:)
  real*8, allocatable :: dpsi_RMP_cos_dR(:)
  real*8, allocatable :: dpsi_RMP_cos_dZ(:)
  real*8, allocatable :: psi_RMP_sin(:)
  real*8, allocatable :: dpsi_RMP_sin_dR(:)
  real*8, allocatable :: dpsi_RMP_sin_dZ(:) 
  integer             :: RMP_har_cos,RMP_har_sin ! Harmoics numbers for RMP-cos and RMP-sin(for ex. ntor=3, nperiod=2,RMP_har_cos=2, RMP_har_sin=3)

  !> @name toroidal rotation profile
  real*8              :: V_0,   V_1,    V_coef(10)! analytical // rotation profile similar to temperature and density in model 303
  character(len=512)  :: R_Z_psi_bnd_file !< ASCII file for R_boundary,Z_boundary, psi_boundary, with n_boundary size.
  character(len=512)  :: wall_file        !< ASCII file for external wall geometry, if n_ext is greater than zero.
  
  !> @name Numerical input profile for the toroidal rotation
  character(len=512)  :: rot_file        !< ASCII file the profile is read from.
  logical             :: num_rot         !< is set true if rot_file /= 'none'
  integer             :: num_rot_len     !< Number of points in profile
  real*8, allocatable :: num_rot_x(:)    !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_rot_y0(:)   !< Values of toroidal rotation profile
  real*8, allocatable :: num_rot_y1(:)   !< First derivatives of toroidal rotation profile with respect to $\Psi_{N}$
  real*8, allocatable :: num_rot_y2(:)   !< Second derivatives of toroidal rotation profile with respect to $\Psi_{N}$ 
  real*8, allocatable :: num_rot_y3(:)   !< Third derivatives of toroidal rotation profile with respect to $\Psi_{N}$ 

  !> @name gmres parameters
  integer             :: iter_precon    !< if number of gmres iterations > iter_precon, the preconditioner is updated
  integer             :: gmres_m        !< gmres restart (dimension)
  real*8              :: gmres_4        !< see gmres manual (error ratio between preconditioned and non-preconditioned error
  real*8              :: gmres_tol      !< the tolerance for the gmres iterations

  !> @name Taylor-Galerkin Stabilisation coefficients
  real*8              :: tgnum(n_var)
  
  !> @name Numerical parameters
  real*8              :: D_prof_neg     !< Diffusion coefficient in regions with negative density
  real*8              :: ZK_prof_neg    !< Diffusion coefficient in regions with negative temperature
  real*8              :: T_min          !< minimum temperature (limits on the temperature dependence of resistivity etc.
  integer             :: n_tor_fft_thresh !< If n_tor >= n_tor_fft_thresh, element_matrix_fft will be used
  integer*8           :: fftw_plan      !< Required for FFTW library
  real*8              :: corr_neg_temp_coef(2) !< Parameters used in models/corr_neg.f90
  real*8              :: corr_neg_dens_coef(2) !< Parameters used in models/corr_neg.f90

  !> @name ECCD current sources
  real*8  :: jecamp             ! parameter, not to be confused with jec_source in element_matrix.f90
  real*8  :: jec_pos1, jec_pos2, jec_pos3, jec_pos4
  real*8  :: jec_width, jec_width2 
  real*8  :: nu_jec_fast         ! 1/collision frequency
  real*8  :: nu_jec1_fast,nu_jec2_fast         ! 1/collision frequency
  real*8  :: mod_jec            ! extra parameters for ECCD
  real*8  :: JJ_par             ! velocity of resonent electrons
  real*8  :: jw1,jw2,jw3        ! parameters to determine current source

  !> @name (Currently unused)
  real*8  :: zjz_0, zjz_1,  zj_coef(10)
  real*8  :: D_neutral
  
end module phys_module
