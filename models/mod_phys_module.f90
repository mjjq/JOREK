!> Input parameters and physical variables.
module phys_module

  use parameters
  
  implicit none
  
  !> @name Various parameters
  real*8  :: eta               !< Resistivity
  logical :: eta_T_dependent   !< Resistivity dependent on temperature? Otherwise constant.
  real*8  :: visco             !< Viscosity
  logical :: visco_T_dependent !< Viscosity dependent on temperature? Otherwise constant.
  real*8  :: visco_par         !< Parallel viscosity
  real*8  :: F0                !< Determines fixed toroidal magnetic field: \f$ B_\phi = F_0/R \f$
  real*8  :: central_density   !< particle density at the magnetic axis (in units of \f$10^{20} m^{-3}\f$)
  real*8  :: gamma             !< ratio of specific heat (=5/3)
  real*8  :: Q_bar             !< (model400)
  real*8  :: sigma             !< (model400)
  real*8  :: tauIC             !< (model302 and 701)
  real*8  :: gamma_sheath      !< sheath boundary condition open fieldlines (model303)
  real*8  :: gamma_sheath_i    !< sheath boundary condition open fieldlines (model400)
  real*8  :: gamma_sheath_e    !< sheath boundary condition open fieldlines (model400)
  real*8  :: density_reflection!< density reflection coeeficient open fieldlines (model303)
  integer :: mode(n_tor)       !< Toroidal mode number corresponding to the JOREK modes, e.g., for n_period=8 and n_tor=3, mode(:)=0,8,8
  integer :: nout              !< Output a restart file every nout timesteps.
  logical :: restart           !< Restart a code run from the restart file jorek_restart.rst?
  logical :: regrid            !< Re-generate the flux-aligned grid (does not work currently)?
  logical :: import_equil
  logical :: xpoint            !< X-point geometry?
  logical :: refinement        !< Use mesh refinement?
  logical :: bc_natural_open   !< use natural boundary conditions on the open fieldlines
  logical :: produce_live_data !< Write data to 'energies.dat', 'growth_rates.dat', and 'times.dat' during the code run?
  logical :: grid_to_wall      !< extend the grid to a physical wall
  logical :: adaptive_time     !< requires no_mpi for Pastix library
  logical :: equil             !< compute equilibrium
  logical :: bench_without_plot!< .true. for benchmark (mesuring elapsed time without plot phases) 
  logical :: gmres             !< Use iterative GMRES solver
  integer :: gmres_max_iter    !< Maximum number of GMRES iterations
  logical :: linear_run        !< Perform a linear run where the equilibrium quantities (i_tor=1) do not change with time?
  
  real*8, allocatable :: energies(:,:,:)  !< Magnetic and kinetic mode energies at timesteps.
  character(len=3)    :: mode_type(n_tor) !< 'cos' or 'sin'
  
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
  
  !> @name Heat and particle diffusivity parameters
  real*8  :: D_perp(10), D_par
  real*8  :: ZK_perp(10), ZK_par, ZK_i_perp(10), ZK_e_perp(10), K_i_par, K_e_par
  
  !> @name Numerical resistivity, viscosity and diffusivities
  real*8  :: eta_num, visco_num, visco_par_num, D_perp_num, Zk_perp_num
  
  !> @name Timestepping parameters
  real*8  :: tstep             !< Size of the timesteps (\f$ \Delta t \f$)
  real*8  :: tstep_n(10)       !< Alternative to tstep: Up to ten values may be given
  integer :: nstep             !< Number of timesteps to perform
  integer :: nstep_n(10)       !< Alternative to nstep: Up to ten values may be given
  real*8  :: t_start           !< Time value at the start of the code run (zero or from restart file)
  real*8  :: t_now             !< Current time value in the simulation
  integer :: index_start       !< Time step index at the beginning of the code run (zero or from restart file)
  integer :: index_now         !< Current time step index
  real*8, allocatable :: xtime(:) !< Time values corresponding to the timesteps.
  
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
  integer :: mf                !< Number of entries in fbnd and fpsi
  real*8  :: fbnd(1026)        !< Fourier expansion of boundary
  real*8  :: fpsi(1026)        !< Fourier expansion of the poloidal flux at the boundary
  
  !> @name Numerical boundary of initial grid
  !! Numerical definition of the boundary of the non flux-aligned initial polar grid.
  integer :: n_boundary        !< Number of points in R_boundary, Z_boundary, psi_boundary.
  real*8  :: R_boundary(1026)  !< Numerical R values defining the boundary
  real*8  :: Z_boundary(1026)  !< Numerical Z values defining the boundary
  real*8  :: psi_boundary(1026)!< Numerical values giving the poloidal flux at the boundary
  
  !> @name Pellet-related input parameters
  real*8  :: pellet_amplitude  !< amplitude of density source (when pellet modelled as density source)
  real*8  :: pellet_R          !< major radius position pellet
  real*8  :: pellet_Z          !< Z position pellet
  real*8  :: pellet_phi        !< width of the pellet cloud (densioty source) in toroidal 
  real*8  :: pellet_radius     !< radius of the simulation pellet
  real*8  :: pellet_sig        !< width of smoothing of density source (arctan((r-pellet_radius)/pellet_sig))
  real*8  :: pellet_length     !< width of smoothing of density source in toroidal angle
  real*8  :: pellet_psi        !< pellet_width in poloidal flux
  real*8  :: pellet_delta_psi  !< width of smoothing in poloidal flux
  real*8  :: pellet_velocity_R !< pellet velocity component radial direction
  real*8  :: pellet_velocity_Z !< pellet velocity component Z direction
  real*8  :: pellet_density    !< pellet density (in units 10^20 m^-3)
  real*8  :: pellet_particles  !< the number of particles in the pellet (in units of 10^20)
  logical :: use_pellet
  
  !> @name Free boundary extension
  !! Input parameters related to the free boundary extension (folder vacuum/).
  logical :: freeboundary_equil!< use a free or fixed boundary equilibrium?
  logical :: freeboundary      !< use free or fixed boundary conditions in time-evolution?
  logical :: use_starwall      !< use the STARWALL vacuum solution? (free boundary only)
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
  integer :: n_radial          !< Number of radial grid points
  integer :: n_pol             !< Number of poloidal grid points
  real*8  :: R_geo             !< Center of the grid
  real*8  :: Z_geo             !< Center of the grid
  
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
  integer :: n_open            !< Number of 'radial' grid points in the open flux region
  integer :: n_private         !< Number of 'radial' grid points in the private flux region
  integer :: n_leg             !< Number of 'poloidal' grid points along the divertor legs
  real*8  :: SIG_closed        !< Width with grid accumulation
  real*8  :: SIG_open          !< Width with grid accumulation
  real*8  :: SIG_private       !< Width with grid accumulation
  real*8  :: SIG_theta         !< Width with grid accumulation
  real*8  :: SIG_leg_0         !< Width with grid accumulation
  real*8  :: SIG_leg_1         !< Width with grid accumulation
  real*8  :: dPSI_open         !< Delta Psi grid extends into the open flux region
  real*8  :: dPSI_private      !< Delta Psi grid extends into the private flux region
  
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
  
  !> @name Numerical input profile for the ion temperature
  character(len=512)  :: Ti_file         !< ASCII file the profile is read from.
  logical             :: num_Ti          !< is set true if T_file /= 'none'
  integer             :: num_Ti_len      !< Number of points in profile
  real*8, allocatable :: num_Ti_x(:)     !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_Ti_y0(:)    !< Values of temperature profile
  real*8, allocatable :: num_Ti_y1(:)    !< First derivatives of temperature profile (\f$ dT/d\Psi_N \f$)
  real*8, allocatable :: num_Ti_y2(:)    !< Second derivatives of temperature profile (\f$ d^2T/d\Psi_N^2 \f$)
  real*8, allocatable :: num_Ti_y3(:)    !< Third derivatives of temperature profile (\f$ d^3T/d\Psi_N^3 \f$)
  
  !> @name Numerical input profile for the electron temperature
  character(len=512)  :: Te_file         !< ASCII file the profile is read from.
  logical             :: num_Te          !< is set true if T_file /= 'none'
  integer             :: num_Te_len      !< Number of points in profile
  real*8, allocatable :: num_Te_x(:)     !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_Te_y0(:)    !< Values of temperature profile
  real*8, allocatable :: num_Te_y1(:)    !< First derivatives of temperature profile (\f$ dT/d\Psi_N \f$)
  real*8, allocatable :: num_Te_y2(:)    !< Second derivatives of temperature profile (\f$ d^2T/d\Psi_N^2 \f$)
  real*8, allocatable :: num_Te_y3(:)    !< Third derivatives of temperature profile (\f$ d^3T/d\Psi_N^3 \f$)
  
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
  
  !> @name (Currently unused)
  real*8  :: zjz_0, zjz_1,  zj_coef(10)
  real*8  :: D_neutral
  
end module phys_module
