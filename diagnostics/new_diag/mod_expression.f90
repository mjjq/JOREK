!> This module calculates "arbitrary" expressions at arbitrary positions for diagnostic purposes.
!!
!! Structure of the result(:,:,:,:) array:
!! - 1st index corresponds to toroidal positions
!! - 2nd and 3rd index correspond to positions in poloidal plane (2nd poloidal, 3rd radial)
!! - 4th index corresponds to the expressions
!!
!! To obtain good performance, evaluate several expressions at several positions in a single call.
!!
!! Note: Add a new expressions to init_expr() and eval_expr() consistently.
module mod_expression
  
  
  
  
  
  use constants
  use mod_parameters
  use mod_position
  use phys_module
  use diffusivities
  use corr_neg
  use mod_basisfunctions
  use mod_bootstrap_functions
  use mod_poloidal_currents
  
  
  
  
  implicit none
  
  
  
  
  
  public
  private add
  
  
  
  
  
  ! --- Constants
  character(len=14), parameter, private :: THIS_MOD_NAME = 'mod_expression'
  integer,           parameter, private :: LEN_NAME      = 12
  integer,           parameter, private :: LEN_DESCR     = 54
  integer,           parameter, private :: N_EXPR_MAX    = 2000
  !   --- Constants selecting a unit system for expression output
  integer,           parameter          :: JOREK_UNITS   = 0 !< Output expressions in JOREK units
  integer,           parameter          :: SI_UNITS      = 1 !< Output expressions in SI units
  
  
  
  
  
  ! --- Externally visible variables.
  !> Expressions evaluate to this value for positions outside the JOREK domain.
  !! This value may be changed by programs using the new_diag package!
  real*8, save :: expr_outside_value = 0.d0
  !> Coordinate expressions evaluate to this value for positions outside the JOREK domain.
  !! This value may be changed by programs using the new_diag package!
  real*8, save :: expr_outside_coord = 0.d0
  
  
  
  
  
  !> Datatype containing information on a single expression.
  type :: t_expr
    character(len=LEN_NAME)  :: name  !< Short name for the expression
    character(len=LEN_DESCR) :: descr !< Brief explanation of the expression
  end type t_expr
  
  ! > List of expressions.
  type :: t_expr_list
    integer                   :: n_expr  = 0      !< Number of expressions in this list
    integer                   :: n_coord = 0      !< Treat the first n_coord expr.s as coordinates
    type(t_expr)              :: expr(N_EXPR_MAX) !< The expressions
  end type t_expr_list
  
  
  
  
  
  ! --- Standard lists of expressions (require init_expr call first!).
  type(t_expr_list), save :: exprs_all, exprs_all_int     !< All available expressions.
  
  
  
  
  
  contains
  
  
  
  
  
  !> Initialization for the module. Should be called before any other 
  subroutine init_expr()
    
    exprs_all%n_expr     = 0
    exprs_all_int%n_expr = 0
    call add(exprs_all, 'R           ', 'Cylindrical Coordinate R (== Major Radius)            ')
    call add(exprs_all, 'Z           ', 'Cylindrical Coordinate Z                              ')
    call add(exprs_all, 'phi         ', 'Cylindrical Coordinate phi                            ')
    call add(exprs_all, 'theta       ', 'Poloidal Angle With Respect to Magnetic Axis          ')
    call add(exprs_all, 'theta_star  ', 'Poloidal Straight Field Line Angle (for flux surfaces)')
    call add(exprs_all, 'length      ', 'Length Along Poloidal Line (for poloidal lines)       ')
    call add(exprs_all, 'r_minor     ', 'Minor Radius From A = r_minor^2 pi (for flux surfaces)')
    call add(exprs_all, 'x           ', 'Cartesian Coordinate x                                ')
    call add(exprs_all, 'y           ', 'Cartesian Coordinate y                                ')
    call add(exprs_all, 'z           ', 'Cartesian Coordinate z (== Cylindrical Z)             ')
    call add(exprs_all, 'Psi_N       ', 'Normalized Poloidal Magnetic Flux                     ')
    call add(exprs_all, 'xjac        ', '2D Jacobian in the Poloidal Plane                     ')
    call add(exprs_all, 't           ', 'Simulation time                                       ')
    call add(exprs_all, 'Psi         ', 'Poloidal Magnetic Flux                                ')
    call add(exprs_all, 'u           ', 'Velocity Stream Function                              ')
    call add(exprs_all, 'Phi         ', 'Electric Potential Phi                                ')
    call add(exprs_all, 'zj          ', 'Toroidal Current Density Multiplied by 1/R            ')
    call add(exprs_all, 'currdens    ', 'Physical Toroidal Current Density (== zj/R)           ')
    call add(exprs_all, 'FFprime_loc ', 'Local FFprime value, calculated from 3D JxB=\grad p   ')
    call add(exprs_all, 'Jpol        ', 'Poloidal current value in the poloidal field direction')
    call add(exprs_all, 'omega       ', 'Toroidal Vorticity Component                          ')
    call add(exprs_all, 'rho         ', 'Mass Density                                          ')
    call add(exprs_all, 'ne          ', 'Electron Density                                      ')
    call add(exprs_all, 'T           ', 'Temperature (Electrons plus Ions)                     ')
    call add(exprs_all, 'Te          ', 'Electron temperature (assuming Ti=Te)                 ')
    call add(exprs_all, 'vpar        ', 'Parallel Velocity (along magnetic field lines)        ')
    call add(exprs_all, 'eta_T       ', 'Temperature Dependent Resistivity                     ')
    call add(exprs_all, 'visco_T     ', 'Temperature Dependent Viscosity                       ')
    call add(exprs_all, 'zkpar_T     ', 'Temperature Dependent Parallel Heat Diffusivity       ')
    call add(exprs_all, 'dprof       ', 'Particle Diffusivity                                  ')
    call add(exprs_all, 'zkprof      ', 'Perpendicular Heat Diffusivity                        ')
    call add(exprs_all, 'pres        ', 'Total Pressure                                        ')
    call add(exprs_all, 'B_abs       ', 'Norm of the Magnetic Field Vector                     ')
    call add(exprs_all, 'B_tor       ', 'Toroidal Magnetic Field Component                     ')
    call add(exprs_all, 'B_R         ', 'Magnetic Field Component Along R                      ')
    call add(exprs_all, 'B_Z         ', 'Vertical Magnetic Field Component                     ')
    call add(exprs_all, 'B_theta     ', 'Poloidal Magnetic Field Component                     ')
    call add(exprs_all, 'Er          ', 'Radial Electric Field                                 ')
    call add(exprs_all, 'Vtheta_i    ', 'Ion Poloidal Velocity                                 ')
    call add(exprs_all, 'Mach_par    ', 'Parallel Mach Number                                  ')
    call add(exprs_all, 'Mach_pol    ', 'Poloidal Mach Number                                  ')
    call add(exprs_all, 'V_sound     ', 'Sound Speed                                           ')
    call add(exprs_all, 'V_neo       ', 'Neoclassical Velocity                                 ')
    call add(exprs_all, 'Vperp_e     ', 'Electron Perpendicular Velocity                       ')
    call add(exprs_all, 'Vperp_i     ', 'Ion Perpendicular Velocity                            ')
    call add(exprs_all, 'V_ExB       ', 'ExB Velocity                                          ')
    call add(exprs_all, 'Vstar_e     ', 'Electron Diamagnetic Velocity                         ')
    call add(exprs_all, 'Vstar_i     ', 'Ion Diamagnetic Velocity                              ')
    call add(exprs_all, 'ki_neo      ', 'Neoclassical Heat Diffusivity                         ')
    call add(exprs_all, 'mu_neo      ', 'Neoclassical Friction Coefficient                     ')
    call add(exprs_all, 'T_e         ', 'Electron temperature                                  ')
    call add(exprs_all, 'T_i         ', 'Ion temperature                                       ')
    call add(exprs_all, 'E_||        ', 'E_|| for RE acceleration                              ')
    call add(exprs_all, 'E_crit      ', 'E_crit for RE acceleration                            ')
#if JOREK_MODEL == 303 || JOREK_MODEL == 333 || JOREK_MODEL == 400 || JOREK_MODEL == 500
    call add(exprs_all, 'J_bootstrap ', 'Bootstrap Current                                     ')
#endif
#if JOREK_MODEL == 500
    call add(exprs_all, 'radiation   ', 'Radiation terms for bolometry diagnostic              ')
    call add(exprs_all, 'brem        ', 'Brem terms for bolometry diagnostic                   ')
#endif
    ! --- List of volume and boundary integrals
    call add(exprs_all_int, 'E_tot       ', 'Total energy                                          ')
    call add(exprs_all_int, 'Wmag_tot    ', 'Total magnetic energy                                 ')
    call add(exprs_all_int, 'Ohmic_tot   ', 'Total ohmic heating                                   ')
    call add(exprs_all_int, 'Thermal_tot ', 'Total thermal energy                                  ')
    call add(exprs_all_int, 'Helicity_tot', 'Total magnetic helicity                               ')
    call add(exprs_all_int, 'Ip_tot      ', 'Total toroidal plasma current                         ')
    call add(exprs_all_int, 'Kin_par_tot ', 'Total parallel kinetic energy                         ')
    call add(exprs_all_int, 'Kin_perp_tot', 'Total perpendicular kinetic energy                    ')
    call add(exprs_all_int, 'Mag_work_tot', 'Total magnetic work = -\int v\cdot(JxB) dV            ')
    call add(exprs_all_int, 'Thm_work_tot', 'Total thermal work  = \int vpar\cdot\nabla p dV       ')
    call add(exprs_all_int, 'Part_src_tot', 'Total particle source                                 ')
    call add(exprs_all_int, 'Heat_src_tot', 'Total heat source                                     ')
    call add(exprs_all_int, 'Viscpar_diss', 'Total parallel viscosity dissipation                  ')
    call add(exprs_all_int, 'Wmag_src_tot', 'Total magnetic energy source (from current source)    ')
    call add(exprs_all_int, 'li3         ', 'Internal inductance inside LCFS, li(3)                ')
    call add(exprs_all_int, 'li3_tot     ', 'Internal inductance inside grid, li(3)                ')
    call add(exprs_all_int, 'betap       ', 'Poloidal beta, of the plasma inside LCFS              ')
    call add(exprs_all_int, 'area        ', 'Poloidal cross section area inside LCFS               ')
    call add(exprs_all_int, 'volume      ', 'Plasma volume, inside LCFS                            ')
    call add(exprs_all_int, 'P_vn        ', 'Boundary flux of outgoing pressure                    ')
    call add(exprs_all_int, 'qn_par      ', 'Boundary flux of the parallel thermal conduction      ')
    call add(exprs_all_int, 'qn_perp     ', 'Boundary flux of the perpendicular thermal conduction ')
    call add(exprs_all_int, 'kinpar_flux ', 'Boundary flux of parallel kinetic energy              ')
 
  end subroutine init_expr
  
  
  
  
  
  !> [Private] Auxilliary routine for init_expr.
  subroutine add(expr_list, name, descr)
    
    ! --- Routine parameters
    type(t_expr_list),        intent(inout) :: expr_list
    character(len=LEN_NAME),  intent(in)    :: name
    character(len=LEN_DESCR), intent(in)    :: descr
    
    expr_list%n_expr = expr_list%n_expr + 1
    expr_list%expr(expr_list%n_expr)%name  = name
    expr_list%expr(expr_list%n_expr)%descr = descr
    
  end subroutine add
  
  
  
  
  
  !> Creates a subset of all available expressions.
  function exprs(name, n_expr, n_coord) result(expr_list)
    type(t_expr_list) :: expr_list
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':exprs'
    
    ! --- Routine parameters
    character(len=*),  intent(in) :: name(n_expr)
    integer,           intent(in) :: n_expr
    integer, optional, intent(in) :: n_coord
    
    ! --- Local variables
    integer :: i, j, k
    
    k = 0
    do i = 1, n_expr
      j = get_expr_num(exprs_all, trim(name(i)))
      if ( j < 1 ) then
        write(*,*) 'WARNING in '//trim(THIS_ROUTINE_NAME)//': Unknown expression "'//trim(name(i)) &
          //'" ignored.'
        cycle
      end if
      k = k + 1
      expr_list%expr(k) = exprs_all%expr(j)
    end do
    expr_list%n_expr = k
    
    expr_list%n_coord = 0
    if ( present(n_coord) ) expr_list%n_coord = n_coord
    
  end function exprs
  
  
  !> Creates a subset of all available expressions.
  function exprs_int(name, n_expr, n_coord) result(expr_list)
    type(t_expr_list) :: expr_list
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':exprs'
    
    ! --- Routine parameters
    character(len=*),  intent(in) :: name(n_expr)
    integer,           intent(in) :: n_expr
    integer, optional, intent(in) :: n_coord
    
    ! --- Local variables
    integer :: i, j, k
    
    k = 0
    do i = 1, n_expr
      j = get_expr_num_int(exprs_all_int, trim(name(i)))
      if ( j < 1 ) then
        write(*,*) 'WARNING in '//trim(THIS_ROUTINE_NAME)//': Unknown expression "'//trim(name(i)) &
          //'" ignored.'
        cycle
      end if
      k = k + 1
      expr_list%expr(k) = exprs_all_int%expr(j)
    end do
    expr_list%n_expr = k
    
    expr_list%n_coord = 0
    if ( present(n_coord) ) expr_list%n_coord = n_coord
    
  end function exprs_int
  
  
  
  !> Merge several expression lists.    ### NOT USED AT PRESENT
  recursive function join_exprs(list1, list2, list3, list4, list5, list6, list7, list8, list9)     &
    result(expr_list)
    type(t_expr_list) :: expr_list
    
    ! --- Routine parameters
    type(t_expr_list),           intent(in) :: list1, list2
    type(t_expr_list), optional, intent(in) :: list3, list4, list5, list6, list7, list8, list9
    
    integer :: i, j
    logical :: found
    
    expr_list = list1
    
    do i = 1, list2%n_expr
      
      ! --- Already in the list? (avoid duplicates)
      found = .false.
      do j = 1, expr_list%n_expr
        if ( expr_list%expr(j)%name == list2%expr(i)%name ) then
          found = .true.
          exit
        end if
      end do
      
      ! --- Add to expr_list
      if ( .not. found ) then
        expr_list%n_expr = expr_list%n_expr + 1
        expr_list%expr(expr_list%n_expr) = list2%expr(i)
      end if
      
    end do
    
    ! --- Recursive joins of the remaining lists
    if ( present(list3) ) expr_list = join_exprs(expr_list, list3)
    if ( present(list4) ) expr_list = join_exprs(expr_list, list4)
    if ( present(list5) ) expr_list = join_exprs(expr_list, list5)
    if ( present(list6) ) expr_list = join_exprs(expr_list, list6)
    if ( present(list7) ) expr_list = join_exprs(expr_list, list7)
    if ( present(list8) ) expr_list = join_exprs(expr_list, list8)
    if ( present(list9) ) expr_list = join_exprs(expr_list, list9)
    
  end function join_exprs
  
  
  
  
  
  !> Prints an expression list in a readable table format.
  subroutine print_exprs(expr_list, short)
    
    ! --- Routine parameters
    type(t_expr_list),  intent(in) :: expr_list
    logical, optional,  intent(in) :: short
    
    ! --- Local variables
    integer :: i
    
    if ( present(short) .and. (short) ) then
      
      905 format(1x,a,',')
      do i = 1, expr_list%n_expr
        write(*,905,advance='no') trim(expr_list%expr(i)%name)
      end do
      write(*,*)
      
    else
      
      900 format(1x,a)
      901 format(1x,i6.6,' | ',a,' | ',a)
      902 format(1x,80('-'))
      
      write(*,*)
      write(*,*) 'List of Diagnostic Expressions:'
      write(*,*)
      
      write(*,902)
      write(*,900) 'Number | Name         | Description'
      write(*,902)
      do i = 1, expr_list%n_expr
        write(*,901) i, expr_list%expr(i)%name, expr_list%expr(i)%descr
      end do
      write(*,902)
      write(*,*)
      
    end if
    
  end subroutine print_exprs
  
  
  
  
  !> Find out expression number in an expression list.
  integer function get_expr_num(expr_list, name) result(num)
    
    ! --- Routine parameters
    type(t_expr_list),       intent(in) :: expr_list
    character(len=*), intent(in) :: name
    
    ! --- Local variables
    integer :: i
    
    num = -99
    do i = 1, exprs_all%n_expr
      if ( trim(exprs_all%expr(i)%name) == trim(name) ) then
        num = i
        exit
      end if
    end do
    
  end function get_expr_num
  
  
   !> Find out expression number in an expression list.
  integer function get_expr_num_int(expr_list, name) result(num)
    
    ! --- Routine parameters
    type(t_expr_list),       intent(in) :: expr_list
    character(len=*), intent(in) :: name
    
    ! --- Local variables
    integer :: i
    
    num = -99
    do i = 1, exprs_all_int%n_expr
      if ( trim(exprs_all_int%expr(i)%name) == trim(name) ) then
        num = i
        exit
      end if
    end do
    
  end function get_expr_num_int
  
  
  
  !> Check that an expression exists.
  logical function expr_exists(expr_list, name) result(exists)
    
    ! --- Routine parameters
    type(t_expr_list),       intent(in) :: expr_list
    character(len=LEN_NAME), intent(in) :: name
    
    exists = ( get_expr_num(expr_list, name) > 0 )
    
  end function expr_exists
  
  
  
  
  
  !> Get the description of an expression.
  character(len=LEN_DESCR) function expr_descr(name) result(descr)
    
    ! --- Routine parameters
    character(len=LEN_NAME) :: name
    
    integer :: num
    
    num = get_expr_num(exprs_all, name)
    
    if ( num > 0 ) then
      descr = exprs_all%expr(num)%descr
    else
      descr = '<UNKNOWN EXPRESSION>'
    end if
    
  end function expr_descr
  
  
  
  
  
  !> Evaluate one/several expressions at one/several poloidal and one/several toroidal positions.
  subroutine eval_expr(eq, units, expr_list, pol_pos_list, tor_pos_list, result, ierr)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':eval_expr'
    
    ! --- Routine parameters
    type(t_equil_state),          intent(in)    :: eq
    integer,                      intent(in)    :: units !< Output in which units?
    type(t_expr_list),            intent(in)    :: expr_list
    type(t_pol_pos_list), target, intent(in)    :: pol_pos_list
    type(t_tor_pos_list), target, intent(in)    :: tor_pos_list
    real*8, allocatable,          intent(inout) :: result(:,:,:,:)
    integer,                      intent(out)   :: ierr
    
    ! --- Local variables
    type(t_pol_pos), pointer :: pol_pos
    type(t_tor_pos), pointer :: tor_pos
    type(type_element)       :: element
    type(type_node)          :: nodes(n_vertex_max)
    integer :: ipolpos, jpolpos, itorpos, iexpr, ielm, i, j, k, i_tor
    real*8  :: xjac, xjac_R, xjac_Z, R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, &
      s, t, H(n_vertex_max,n_order+1), H_s(n_vertex_max,n_order+1), H_t(n_vertex_max,n_order+1),   &
      H_st(n_vertex_max,n_order+1), H_ss(n_vertex_max,n_order+1), H_tt(n_vertex_max,n_order+1),    &
      HZ(n_tor), HZ_p(n_tor), HZ_pp(n_tor), phi, res, BigR, BigR_R, x_cart, y_cart, theta
    real*8  :: ps0, ps0_s, ps0_t, ps0_ss, ps0_tt, ps0_st, ps0_p, ps0_pp, u0, u0_s, u0_t, u0_ss,    &
      u0_tt, u0_st, u0_p, u0_pp, zj0, zj0_s, zj0_t, zj0_ss, zj0_tt, zj0_st, zj0_p, zj0_pp, w0,     &
      w0_s, w0_t, w0_ss, w0_tt, w0_st, w0_p, w0_pp, r0, r0_s, r0_t, r0_ss, r0_tt, r0_st, r0_p,     &
      r0_pp, T0, T0_s, T0_t, T0_ss, T0_tt, T0_st, T0_p, T0_pp, Vpar0, Vpar0_s, Vpar0_t, Vpar0_ss,  &
      Vpar0_tt, Vpar0_st, Vpar0_p, Vpar0_pp, psi_norm
    real*8  :: ps0_R, ps0_Z, ps0_RR, ps0_ZZ, ps0_RZ, u0_R, u0_Z, u0_RR, u0_ZZ, u0_RZ, vv2, zj0_R,  &
      zj0_Z, zj0_RR, zj0_ZZ, zj0_RZ, w0_R, w0_Z, w0_RR, w0_ZZ, w0_RZ, r0_R, r0_Z, r0_RR, r0_ZZ,    &
      r0_RZ, r0_hat, r0_R_hat, r0_Z_hat, T0_R, T0_Z, T0_RR, T0_ZZ, T0_RZ, T0_ps0_R, T0_ps0_Z,      &
      Vpar0_R, Vpar0_Z, Vpar0_RR, Vpar0_ZZ, Vpar0_RZ, P0, P0_R, P0_Z, P0_s, P0_t, P0_p, P0_pp,     &
      P0_RR, P0_ZZ, P0_RZ, BB2, B_tor, B_R, B_Z, Btheta, psi_abs, E_par, E_crit, ln_Lambda
    real*8  :: eta_T, deta_dT, d2eta_d2T, visco_T, dvisco_dT, ZKpar_T, dZKpar_dT, D_prof, ZK_prof
    real*8 :: Ti0, Ti0_s, Ti0_t, Ti0_st, Ti0_ss, Ti0_tt, Ti0_p, Ti0_pp, Te0, Te0_s, Te0_t, Te0_st, &
      Te0_ss, Te0_tt, Te0_p, Te0_pp, Ti0_R, Ti0_Z, Te0_R, Te0_Z, Er, Vtheta, Mach_par, Mach_pol,   &
      Vsound, Vneo, Vperp_e, Vperp_i, V_ExB, Vstar_e, Vstar_i, mu_neo, ki_neo, J_boot 
    real*8 :: FFprime_loc, Jpol
    real*8 :: hh, hh_s, hh_t, hh_ss, hh_tt, hh_st, hhz, hhz_p, hhz_pp, sz, vv(n_var)
    real*8 :: delta_g(n_var), delta_s(n_var), delta_t(n_var)
    ! --- Normalization factors
    real*8  :: rho_norm, fact_time, fact_mu_zero, fact_ne, fact_rho, fact_T, fact_vpar,            &
      fact_resistiv, fact_Er
#if JOREK_MODEL == 500
    real*8  :: coef_rad_1
    real*8  :: T_rad, LradDrays_T, LradDcont_T
    real*8  :: rn0, rn0_s, rn0_t, rn0_ss, rn0_tt, rn0_st, rn0_p, rn0_pp
    real*8  :: Arad_bg, Brad_bg, Crad_bg, frad_bg, dfrad_bg_dT
#endif
    
    ierr = 0
    
    ! --- Some sanity checks
    if ( expr_list%n_expr < 1 ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': expr_list%n_expr < 1 encountered.'
      ierr = -101
      return
    else if ( minval(pol_pos_list%n_pos) < 1 ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': pol_pos_list%n_pos < 1 encountered.'
      ierr = -102
      return
    else if ( .not. allocated(pol_pos_list%pos) ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': pol_pos_list%pos not allocated.'
      ierr = -103
      return
    else if ( tor_pos_list%n_pos < 1 ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': tor_pos_list%n_pos < 1 encountered.'
      ierr = -102
      return
    else if ( .not. allocated(tor_pos_list%pos) ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': tor_pos_list%pos not allocated.'
      ierr = -103
      return
    end if
    
    if ( allocated(result) ) deallocate(result)
    allocate( result(tor_pos_list%n_pos, pol_pos_list%n_pos(1), pol_pos_list%n_pos(2),             &
      expr_list%n_expr) )
    
    ! --- Loop over positions in the poloidal plane
    loop_pol1: do ipolpos = 1, pol_pos_list%n_pos(1)
      loop_pol2: do jpolpos = 1, pol_pos_list%n_pos(2)
        pol_pos => pol_pos_list%pos(ipolpos,jpolpos)
        
        if ( pol_pos%outside ) then
          result(:, ipolpos, jpolpos, :expr_list%n_coord  ) = expr_outside_coord
          result(:, ipolpos, jpolpos, expr_list%n_coord+1:) = expr_outside_value
          cycle
        end if
        
        ! --- Poloidal Coordinates and Basis Functions
        R    = pol_pos%R
        R_s  = pol_pos%R_s
        R_t  = pol_pos%R_t
        R_st = pol_pos%R_st
        R_ss = pol_pos%R_ss
        R_tt = pol_pos%R_tt
        Z    = pol_pos%Z
        Z_s  = pol_pos%Z_s
        Z_t  = pol_pos%Z_t
        Z_st = pol_pos%Z_st
        Z_ss = pol_pos%Z_ss
        Z_tt = pol_pos%Z_tt
        s    = pol_pos%s
        t    = pol_pos%t
        ielm = pol_pos%ielm
        BigR   = R    ! Just two different names for R
        BigR_R = 1.d0 ! Trivial derivative
        call basisfunctions(s, t, H, H_s, H_t, H_st, H_ss, H_tt)
        
        ! --- Poloidal angle theta
        theta = atan2(Z-eq%Z_axis, R-eq%R_axis)
        if ( theta < 0.d0 ) theta = theta + 2.d0*PI
        
        ! --- 2D Jacobian
        xjac   = R_s * Z_t - R_t * Z_s
        xjac_R = ( R_ss * Z_t**2 - 2*R_st * Z_s*Z_t + R_tt * Z_s**2 + R_s * (Z_st*Z_t - Z_tt*Z_s )   &
                 + R_t * (Z_st*Z_s - Z_ss*Z_t ) ) / xjac
        xjac_Z = ( Z_ss * R_t**2 - 2*Z_st * R_s*R_t + Z_tt * R_s**2 + Z_s * (R_st*R_t - R_tt*R_s )   &
                 + Z_t * (R_st*R_s - R_ss*R_t ) ) / xjac
        
        ! --- Elements and nodes
        element  = pol_pos%element
        nodes(:) = pol_pos%nodes(:)
        
        ! --- Loop over toroidal positions
        loop_tor: do itorpos = 1, tor_pos_list%n_pos
          tor_pos => tor_pos_list%pos(itorpos)
          
          ! --- Toroidal Coordinate and Basis Functions
          phi     = tor_pos%phi
          HZ(1)   = 1.d0
          HZ_p(1) = 0.d0
          do i = 1, (n_tor-1) / 2
            HZ(2*i)      =                           cos(mode(2*i)  *phi)
            HZ_p(2*i)    = - float(mode(2*i))      * sin(mode(2*i)  *phi)
            HZ_pp(2*i)   = - float(mode(2*i))**2   * cos(mode(2*i)  *phi)
            HZ(2*i+1)    =                           sin(mode(2*i+1)*phi)
            HZ_p(2*i+1)  = + float(mode(2*i+1))    * cos(mode(2*i+1)*phi)
            HZ_pp(2*i+1) = - float(mode(2*i+1))**2 * sin(mode(2*i+1)*phi)
          end do
          
          ! --- Cartesian Coordinates
          x_cart = + R*cos(phi)
          y_cart = - R*sin(phi)
          
          ps0   = 0.d0; ps0_s   = 0.d0; ps0_t   = 0.d0; ps0_ss   = 0.d0; ps0_tt   = 0.d0; ps0_st   = 0.d0; ps0_p   = 0.d0; ps0_pp   = 0.d0
          u0    = 0.d0; u0_s    = 0.d0; u0_t    = 0.d0; u0_ss    = 0.d0; u0_tt    = 0.d0; u0_st    = 0.d0; u0_p    = 0.d0; u0_pp    = 0.d0
          zj0   = 0.d0; zj0_s   = 0.d0; zj0_t   = 0.d0; zj0_ss   = 0.d0; zj0_tt   = 0.d0; zj0_st   = 0.d0; zj0_p   = 0.d0; zj0_pp   = 0.d0
          w0    = 0.d0; w0_s    = 0.d0; w0_t    = 0.d0; w0_ss    = 0.d0; w0_tt    = 0.d0; w0_st    = 0.d0; w0_p    = 0.d0; w0_pp    = 0.d0
          r0    = 0.d0; r0_s    = 0.d0; r0_t    = 0.d0; r0_ss    = 0.d0; r0_tt    = 0.d0; r0_st    = 0.d0; r0_p    = 0.d0; r0_pp    = 0.d0
          T0    = 0.d0; T0_s    = 0.d0; T0_t    = 0.d0; T0_ss    = 0.d0; T0_tt    = 0.d0; T0_st    = 0.d0; T0_p    = 0.d0; T0_pp    = 0.d0
          Ti0   = 0.d0; Ti0_s   = 0.d0; Ti0_t   = 0.d0; Ti0_ss   = 0.d0; Ti0_tt   = 0.d0; Ti0_st   = 0.d0; Ti0_p   = 0.d0; Ti0_pp   = 0.d0
          Te0   = 0.d0; Te0_s   = 0.d0; Te0_t   = 0.d0; Te0_ss   = 0.d0; Te0_tt   = 0.d0; Te0_st   = 0.d0; Te0_p   = 0.d0; Te0_pp   = 0.d0
          Vpar0 = 0.d0; Vpar0_s = 0.d0; Vpar0_t = 0.d0; Vpar0_ss = 0.d0; Vpar0_tt = 0.d0; Vpar0_st = 0.d0; Vpar0_p = 0.d0; Vpar0_pp = 0.d0
          delta_g(:) = 0.d0; delta_s(:) = 0.d0; delta_t(:) = 0.d0
#if JOREK_MODEL == 500
          rn0 = 0.d0
          rn0_s = 0.0
          rn0_t = 0.0
          rn0_ss = 0.0
          rn0_tt = 0.0
          rn0_st = 0.0
          rn0_p = 0.0
          rn0_pp = 0.0
#endif
          
          ! --- Reconstruct variables
          do i = 1, n_vertex_max
            do j = 1, n_order+1
              
              sz    = element%size(i,j)
              hh    = H   (i,j)
              hh_s  = H_s (i,j)
              hh_t  = H_t (i,j)
              hh_ss = H_ss(i,j)
              hh_tt = H_tt(i,j)
              hh_st = H_st(i,j)
              
              do i_tor = 1, n_tor
                
                hhz    = HZ   (i_tor)
                hhz_p  = HZ_p (i_tor)
                hhz_pp = HZ_pp(i_tor)
                vv(:)  = nodes(i)%values(i_tor,j,:)
                
                ! --- Poloidal Flux
                ps0      = ps0      + vv(1) * sz * hh    * hhz
                ps0_s    = ps0_s    + vv(1) * sz * hh_s  * hhz
                ps0_t    = ps0_t    + vv(1) * sz * hh_t  * hhz
                ps0_ss   = ps0_ss   + vv(1) * sz * hh_ss * hhz
                ps0_tt   = ps0_tt   + vv(1) * sz * hh_tt * hhz
                ps0_st   = ps0_st   + vv(1) * sz * hh_st * hhz
                ps0_p    = ps0_p    + vv(1) * sz * hh    * hhz_p
                ps0_pp   = ps0_pp   + vv(1) * sz * hh    * hhz_pp
                
                ! --- Stream Function
                u0       = u0       + vv(2) * sz * hh    * hhz
                u0_s     = u0_s     + vv(2) * sz * hh_s  * hhz
                u0_t     = u0_t     + vv(2) * sz * hh_t  * hhz
                u0_ss    = u0_ss    + vv(2) * sz * hh_ss * hhz
                u0_tt    = u0_tt    + vv(2) * sz * hh_tt * hhz
                u0_st    = u0_st    + vv(2) * sz * hh_st * hhz
                u0_p     = u0_p     + vv(2) * sz * hh    * hhz_p
                u0_pp    = u0_pp    + vv(2) * sz * hh    * hhz_pp
                
                ! --- Current
                zj0      = zj0      + vv(3) * sz * hh    * hhz
                zj0_s    = zj0_s    + vv(3) * sz * hh_s  * hhz
                zj0_t    = zj0_t    + vv(3) * sz * hh_t  * hhz
                zj0_ss   = zj0_ss   + vv(3) * sz * hh_ss * hhz
                zj0_tt   = zj0_tt   + vv(3) * sz * hh_tt * hhz
                zj0_st   = zj0_st   + vv(3) * sz * hh_st * hhz
                zj0_p    = zj0_p    + vv(3) * sz * hh    * hhz_p
                zj0_pp   = zj0_pp   + vv(3) * sz * hh    * hhz_pp
                
                ! --- Vorticity
                w0       = w0       + vv(4) * sz * hh    * hhz
                w0_s     = w0_s     + vv(4) * sz * hh_s  * hhz
                w0_t     = w0_t     + vv(4) * sz * hh_t  * hhz
                w0_ss    = w0_ss    + vv(4) * sz * hh_ss * hhz
                w0_tt    = w0_tt    + vv(4) * sz * hh_tt * hhz
                w0_st    = w0_st    + vv(4) * sz * hh_st * hhz
                w0_p     = w0_p     + vv(4) * sz * hh    * hhz_p
                w0_pp    = w0_pp    + vv(4) * sz * hh    * hhz_pp
                
                ! --- Density
                r0       = r0       + vv(5) * sz * hh    * hhz
                r0_s     = r0_s     + vv(5) * sz * hh_s  * hhz
                r0_t     = r0_t     + vv(5) * sz * hh_t  * hhz
                r0_ss    = r0_ss    + vv(5) * sz * hh_ss * hhz
                r0_tt    = r0_tt    + vv(5) * sz * hh_tt * hhz
                r0_st    = r0_st    + vv(5) * sz * hh_st * hhz
                r0_p     = r0_p     + vv(5) * sz * hh    * hhz_p
                r0_pp    = r0_pp    + vv(5) * sz * hh    * hhz_pp
                
#if JOREK_MODEL == 400
                ! --- Ion temperature
                Ti0       = Ti0       + vv(6) * sz * hh    * hhz
                Ti0_s     = Ti0_s     + vv(6) * sz * hh_s  * hhz
                Ti0_t     = Ti0_t     + vv(6) * sz * hh_t  * hhz
                Ti0_ss    = Ti0_ss    + vv(6) * sz * hh_ss * hhz
                Ti0_tt    = Ti0_tt    + vv(6) * sz * hh_tt * hhz
                Ti0_st    = Ti0_st    + vv(6) * sz * hh_st * hhz
                Ti0_p     = Ti0_p     + vv(6) * sz * hh    * hhz_p
                Ti0_pp    = Ti0_pp    + vv(6) * sz * hh    * hhz_pp
                
                ! --- Electron temperature
                Te0       = Te0       + vv(8) * sz * hh    * hhz
                Te0_s     = Te0_s     + vv(8) * sz * hh_s  * hhz
                Te0_t     = Te0_t     + vv(8) * sz * hh_t  * hhz
                Te0_ss    = Te0_ss    + vv(8) * sz * hh_ss * hhz
                Te0_tt    = Te0_tt    + vv(8) * sz * hh_tt * hhz
                Te0_st    = Te0_st    + vv(8) * sz * hh_st * hhz
                Te0_p     = Te0_p     + vv(8) * sz * hh    * hhz_p
                Te0_pp    = Te0_pp    + vv(8) * sz * hh    * hhz_pp
#else
                ! --- Temperature (ion + electron) in models .ne. 400
                T0       = T0       + vv(6) * sz * hh    * hhz
                T0_s     = T0_s     + vv(6) * sz * hh_s  * hhz
                T0_t     = T0_t     + vv(6) * sz * hh_t  * hhz
                T0_ss    = T0_ss    + vv(6) * sz * hh_ss * hhz
                T0_tt    = T0_tt    + vv(6) * sz * hh_tt * hhz
                T0_st    = T0_st    + vv(6) * sz * hh_st * hhz
                T0_p     = T0_p     + vv(6) * sz * hh    * hhz_p
                T0_pp    = T0_pp    + vv(6) * sz * hh    * hhz_pp
#endif
                ! --- Parallel Velocity
#if JOREK_MODEL >= 300
                Vpar0    = Vpar0    + vv(7) * sz * hh    * hhz
                Vpar0_s  = Vpar0_s  + vv(7) * sz * hh_s  * hhz
                Vpar0_t  = Vpar0_t  + vv(7) * sz * hh_t  * hhz
                Vpar0_ss = Vpar0_ss + vv(7) * sz * hh_ss * hhz
                Vpar0_tt = Vpar0_tt + vv(7) * sz * hh_tt * hhz
                Vpar0_st = Vpar0_st + vv(7) * sz * hh_st * hhz
                Vpar0_p  = Vpar0_p  + vv(7) * sz * hh    * hhz_p
                Vpar0_pp = Vpar0_pp + vv(7) * sz * hh    * hhz_pp
#endif
#if JOREK_MODEL == 500
                rn0       = rn0       + vv(8) * sz * hh    * hhz
                rn0_s     = rn0_s     + vv(8) * sz * hh_s  * hhz
                rn0_t     = rn0_t     + vv(8) * sz * hh_t  * hhz
                rn0_ss    = rn0_ss    + vv(8) * sz * hh_ss * hhz
                rn0_tt    = rn0_tt    + vv(8) * sz * hh_tt * hhz
                rn0_st    = rn0_st    + vv(8) * sz * hh_st * hhz
                rn0_p     = rn0_p     + vv(8) * sz * hh    * hhz_p
                rn0_pp    = rn0_pp    + vv(8) * sz * hh    * hhz_pp
#endif

                ! --- Deltas
                do k = 1, n_var
                  delta_g(k) = delta_g(k) + nodes(i)%deltas(i_tor,j,k) * sz * hh    * hhz
                  delta_s(k) = delta_s(k) + nodes(i)%deltas(i_tor,j,k) * sz * hh_s  * hhz
                  delta_t(k) = delta_t(k) + nodes(i)%deltas(i_tor,j,k) * sz * hh_t  * hhz
                end do
                
              end do
            end do
          end do
          
#if JOREK_MODEL == 400
          ! --- Sum up electron and ion temperature for model400 (e.g., to calculate total pressure)
          T0       = Ti0    + Te0   
          T0_s     = Ti0_s  + Te0_s 
          T0_t     = Ti0_t  + Te0_t 
          T0_ss    = Ti0_ss + Te0_ss
          T0_tt    = Ti0_tt + Te0_tt
          T0_st    = Ti0_st + Te0_st
          T0_p     = Ti0_p  + Te0_p 
          T0_pp    = Ti0_pp + Te0_pp
#endif
          
          ! --- Construct Cartesian Derivatives of Variables.
          ps0_R    = (   Z_t * ps0_s - Z_s * ps0_t ) / xjac
          ps0_Z    = ( - R_t * ps0_s + R_s * ps0_t ) / xjac
          ps0_RR   = (ps0_ss * Z_t**2 - 2.d0*ps0_st * Z_s*Z_t + ps0_tt * Z_s**2  &
                      + ps0_s * (Z_st*Z_t - Z_tt*Z_s )                           &
                      + ps0_t * (Z_st*Z_s - Z_ss*Z_t ) )       / xjac**2         &
                    - xjac_R * (ps0_s * Z_t - ps0_t * Z_s)     / xjac**2
          ps0_ZZ   = (ps0_ss * R_t**2 - 2.d0*ps0_st * R_s*R_t + ps0_tt * R_s**2  &
                      + ps0_s * (R_st*R_t - R_tt*R_s )                           &
                      + ps0_t * (R_st*R_s - R_ss*R_t ) )       / xjac**2         &
                    - xjac_Z * (- ps0_s * R_t + ps0_t * R_s )  / xjac**2
          ps0_RZ   = (- ps0_ss * Z_t*R_t - ps0_tt * R_s*Z_s                      &
                      + ps0_st * (Z_s*R_t  + Z_t*R_s  )                          &
                      - ps0_s  * (R_st*Z_t - R_tt*Z_s )                          &
                      - ps0_t  * (R_st*Z_s - R_ss*Z_t )  )     / xjac**2         &
                    - xjac_R * (- ps0_s * R_t + ps0_t * R_s )  / xjac**2
          
          u0_R     = (   Z_t * u0_s - Z_s * u0_t ) / xjac
          u0_Z     = ( - R_t * u0_s + R_s * u0_t ) / xjac
          u0_RR    = (u0_ss * Z_t**2 - 2.d0*u0_st * Z_s*Z_t + u0_tt * Z_s**2  &
                      + u0_s * (Z_st*Z_t - Z_tt*Z_s )                         &
                      + u0_t * (Z_st*Z_s - Z_ss*Z_t ) )      / xjac**2        &
                    - xjac_R * (u0_s * Z_t - u0_t * Z_s)     / xjac**2
          u0_ZZ    = (u0_ss * R_t**2 - 2.d0*u0_st * R_s*R_t + u0_tt * R_s**2  &
                      + u0_s * (R_st*R_t - R_tt*R_s )                         &
                      + u0_t * (R_st*R_s - R_ss*R_t ) )      / xjac**2        &
                    - xjac_Z * (- u0_s * R_t + u0_t * R_s )  / xjac**2
          u0_RZ    = (- u0_ss * Z_t*R_t - u0_tt * R_s*Z_s                     &
                      + u0_st * (Z_s*R_t  + Z_t*R_s  )                        &
                      - u0_s  * (R_st*Z_t - R_tt*Z_s )                        &
                      - u0_t  * (R_st*Z_s - R_ss*Z_t )  )    / xjac**2        &
                    - xjac_R * (- u0_s * R_t + u0_t * R_s )  / xjac**2
          vv2           = BigR**2 *  ( u0_R * u0_R + u0_Z *u0_Z  )
          !----------------- simplified version of 2nd derivatives (for some unknown reason this is more stable!)
          !u0_RR = (  u0_ss * Z_t**2  + u0_tt * Z_s**2  - 2.d0*u0_st * Z_s*Z_t                  ) / xjac**2
          !u0_ZZ = (  u0_ss * R_t**2  + u0_tt * R_s**2  - 2.d0*u0_st * R_s*R_t                  ) / xjac**2
          !u0_RZ = (- u0_ss * Z_t*R_t - u0_tt * R_s*Z_s +      u0_st * (Z_s*R_t + Z_t*R_s) ) / xjac**2
          
          zj0_R    = (   Z_t * zj0_s - Z_s * zj0_t ) / xjac
          zj0_Z    = ( - R_t * zj0_s + R_s * zj0_t ) / xjac
          zj0_RR   = (zj0_ss * Z_t**2 - 2.d0*zj0_st * Z_s*Z_t + zj0_tt * Z_s**2  &
                      + zj0_s * (Z_st*Z_t - Z_tt*Z_s )                           &
                      + zj0_t * (Z_st*Z_s - Z_ss*Z_t ) )       / xjac**2         &
                    - xjac_R * (zj0_s * Z_t - zj0_t * Z_s)     / xjac**2
          zj0_ZZ   = (zj0_ss * R_t**2 - 2.d0*zj0_st * R_s*R_t + zj0_tt * R_s**2  &
                      + zj0_s * (R_st*R_t - R_tt*R_s )                           &
                      + zj0_t * (R_st*R_s - R_ss*R_t ) )       / xjac**2         &
                    - xjac_Z * (- zj0_s * R_t + zj0_t * R_s )  / xjac**2
          zj0_RZ   = (- zj0_ss * Z_t*R_t - zj0_tt * R_s*Z_s                      &
                      + zj0_st * (Z_s*R_t  + Z_t*R_s  )                          &
                      - zj0_s  * (R_st*Z_t - R_tt*Z_s )                          &
                      - zj0_t  * (R_st*Z_s - R_ss*Z_t )  )     / xjac**2         &
                    - xjac_R * (- zj0_s * R_t + zj0_t * R_s )  / xjac**2
          
          w0_R     = (   Z_t * w0_s - Z_s * w0_t ) / xjac
          w0_Z     = ( - R_t * w0_s + R_s * w0_t ) / xjac
          w0_RR    = (w0_ss * Z_t**2 - 2.d0*w0_st * Z_s*Z_t + w0_tt * Z_s**2  &
                      + w0_s * (Z_st*Z_t - Z_tt*Z_s )                         &
                      + w0_t * (Z_st*Z_s - Z_ss*Z_t ) )      / xjac**2        &
                    - xjac_R * (w0_s * Z_t - w0_t * Z_s)     / xjac**2
          w0_ZZ    = (w0_ss * R_t**2 - 2.d0*w0_st * R_s*R_t + w0_tt * R_s**2  &
                      + w0_s * (R_st*R_t - R_tt*R_s )                         &
                      + w0_t * (R_st*R_s - R_ss*R_t ) )      / xjac**2        &
                    - xjac_Z * (- w0_s * R_t + w0_t * R_s )  / xjac**2
          w0_RZ    = (- w0_ss * Z_t*R_t - w0_tt * R_s*Z_s                     &
                      + w0_st * (Z_s*R_t  + Z_t*R_s  )                        &
                      - w0_s  * (R_st*Z_t - R_tt*Z_s )                        &
                      - w0_t  * (R_st*Z_s - R_ss*Z_t )  )    / xjac**2        &
                    - xjac_R * (- w0_s * R_t + w0_t * R_s )  / xjac**2
          !----------------- simplified version of 2nd derivatives (for some unknown reason this is more stable!)
          !w0_RR = (  w0_ss * Z_t**2  + w0_tt * Z_s**2  - 2.d0*w0_st * Z_s*Z_t                  ) / xjac**2
          !w0_ZZ = (  w0_ss * R_t**2  + w0_tt * R_s**2  - 2.d0*w0_st * R_s*R_t                  ) / xjac**2
          !w0_RZ = (- w0_ss * Z_t*R_t - w0_tt * R_s*Z_s +      w0_st * (Z_s*R_t + Z_t*R_s) ) / xjac**2
          
          r0_R     = (   Z_t * r0_s - Z_s * r0_t ) / xjac
          r0_Z     = ( - R_t * r0_s + R_s * r0_t ) / xjac
          r0_RR    = (r0_ss * Z_t**2 - 2.d0*r0_st * Z_s*Z_t + r0_tt * Z_s**2  &
                      + r0_s * (Z_st*Z_t - Z_tt*Z_s )                         &
                      + r0_t * (Z_st*Z_s - Z_ss*Z_t ) )      / xjac**2        &
                    - xjac_R * (r0_s * Z_t - r0_t * Z_s)     / xjac**2
          r0_ZZ    = (r0_ss * R_t**2 - 2.d0*r0_st * R_s*R_t + r0_tt * R_s**2  &
                      + r0_s * (R_st*R_t - R_tt*R_s )                         &
                      + r0_t * (R_st*R_s - R_ss*R_t ) )      / xjac**2        &
                    - xjac_Z * (- r0_s * R_t + r0_t * R_s )  / xjac**2
          r0_RZ    = (- r0_ss * Z_t*R_t - r0_tt * R_s*Z_s                     &
                      + r0_st * (Z_s*R_t  + Z_t*R_s  )                        &
                      - r0_s  * (R_st*Z_t - R_tt*Z_s )                        &
                      - r0_t  * (R_st*Z_s - R_ss*Z_t )  )    / xjac**2        &
                    - xjac_R * (- r0_s * R_t + r0_t * R_s )  / xjac**2
          r0_hat   = BigR**2 * r0
          r0_R_hat = 2.d0 * BigR * BigR_R  * r0 + BigR**2 * r0_R
          r0_Z_hat = BigR**2 * r0_Z
          
          T0_R     = (   Z_t * T0_s  - Z_s * T0_t ) / xjac
          T0_Z     = ( - R_t * T0_s  + R_s * T0_t ) / xjac
#if JOREK_MODEL == 400
          Ti0_R     = (   Z_t * Ti0_s  - Z_s * Ti0_t ) / xjac
          Ti0_Z     = ( - R_t * Ti0_s  + R_s * Ti0_t ) / xjac
          Te0_R     = (   Z_t * Te0_s  - Z_s * Te0_t ) / xjac
          Te0_Z     = ( - R_t * Te0_s  + R_s * Te0_t ) / xjac
#else
          ! --- Set electron and ion temperatures to T/2 for diagnostic purposes
          Te0     = T0     / 2.d0
          Te0_s   = T0_s   / 2.d0
          Te0_t   = T0_t   / 2.d0
          Te0_st  = T0_st  / 2.d0
          Te0_ss  = T0_ss  / 2.d0
          Te0_tt  = T0_tt  / 2.d0
          Te0_p   = T0_p   / 2.d0
          Te0_pp  = T0_pp  / 2.d0
          Te0_R   = T0_R   / 2.d0
          Te0_Z   = T0_Z   / 2.d0

          Ti0     = T0     / 2.d0
          Ti0_s   = T0_s   / 2.d0
          Ti0_t   = T0_t   / 2.d0
          Ti0_st  = T0_st  / 2.d0
          Ti0_ss  = T0_ss  / 2.d0
          Ti0_tt  = T0_tt  / 2.d0
          Ti0_p   = T0_p   / 2.d0
          Ti0_pp  = T0_pp  / 2.d0
          Ti0_R   = T0_R   / 2.d0
          Ti0_Z   = T0_Z   / 2.d0
#endif
          T0_RR    = (T0_ss * Z_t**2 - 2.d0*T0_st * Z_s*Z_t + T0_tt * Z_s**2 &
                      + T0_s * (Z_st*Z_t - Z_tt*Z_s )                                 &
                      + T0_t * (Z_st*Z_s - Z_ss*Z_t ) )        / xjac**2        &
                     - xjac_R * (T0_s * Z_t - T0_t * Z_s)     / xjac**2
          T0_ZZ    = (T0_ss * R_t**2 - 2.d0*T0_st * R_s*R_t + T0_tt * R_s**2 &
                      + T0_s * (R_st*R_t - R_tt*R_s )                                 &
                      + T0_t * (R_st*R_s - R_ss*R_t ) )        / xjac**2        &
                     - xjac_Z * (- T0_s * R_t + T0_t * R_s )  / xjac**2
          T0_RZ    = (- T0_ss * Z_t*R_t - T0_tt * R_s*Z_s                         &
                       + T0_st * (Z_s*R_t  + Z_t*R_s  )                           &
                       - T0_s  * (R_st*Z_t - R_tt*Z_s )                           &
                       - T0_t  * (R_st*Z_s - R_ss*Z_t )  )       / xjac**2      &
                       - xjac_R * (- T0_s * R_t + T0_t * R_s )  / xjac**2
          T0_ps0_R = T0_RR * ps0_Z - T0_RZ * ps0_R + T0_R * ps0_RZ - T0_Z * ps0_RR
          T0_ps0_Z = T0_RZ * ps0_Z - T0_ZZ * ps0_R + T0_R * ps0_ZZ - T0_Z * ps0_RZ
          !----------------- simplified version of 2nd derivatives (for some unknown reason this is more stable!)
          !T0_RR = (  T0_ss * Z_t**2  + T0_tt * Z_s**2  - 2.d0*T0_st * Z_s*Z_t             ) / xjac**2
          !T0_ZZ = (  T0_ss * R_t**2  + T0_tt * R_s**2  - 2.d0*T0_st * R_s*R_t             ) / xjac**2
          !T0_RZ = (- T0_ss * Z_t*R_t - T0_tt * R_s*Z_s +      T0_st * (Z_s*R_t + Z_t*R_s) ) / xjac**2
          
          Vpar0_R  = (   Z_t * Vpar0_s - Z_s * Vpar0_t ) / xjac
          Vpar0_Z  = ( - R_t * Vpar0_s + R_s * Vpar0_t ) / xjac
          Vpar0_RR = (Vpar0_ss * Z_t**2 - 2.d0*Vpar0_st * Z_s*Z_t + Vpar0_tt * Z_s**2  &
                      + Vpar0_s * (Z_st*Z_t - Z_tt*Z_s )                               &
                      + Vpar0_t * (Z_st*Z_s - Z_ss*Z_t ) )         / xjac**2           &
                    - xjac_R * (Vpar0_s * Z_t - Vpar0_t * Z_s)     / xjac**2
          Vpar0_ZZ = (Vpar0_ss * R_t**2 - 2.d0*Vpar0_st * R_s*R_t + Vpar0_tt * R_s**2  &
                      + Vpar0_s * (R_st*R_t - R_tt*R_s )                               &
                      + Vpar0_t * (R_st*R_s - R_ss*R_t ) )         / xjac**2           &
                    - xjac_Z * (- Vpar0_s * R_t + Vpar0_t * R_s )  / xjac**2
          Vpar0_RZ = (- Vpar0_ss * Z_t*R_t - Vpar0_tt * R_s*Z_s                        &
                      + Vpar0_st * (Z_s*R_t  + Z_t*R_s  )                                 &
                      - Vpar0_s  * (R_st*Z_t - R_tt*Z_s )                                 &
                      - Vpar0_t  * (R_st*Z_s - R_ss*Z_t )  )       / xjac**2           &
                    - xjac_R * (- Vpar0_s * R_t + Vpar0_t * R_s )  / xjac**2
          
          !delta_u_R  = (   Z_t * delta_s(2) - Z_s * delta_t(2) ) / xjac
          !delta_u_Z  = ( - R_t * delta_s(2) + R_s * delta_t(2) ) / xjac
          !delta_ps_R = (   Z_t * delta_s(1) - Z_s * delta_t(1) ) / xjac
          !delta_ps_Z = ( - R_t * delta_s(1) + R_s * delta_t(1) ) / xjac
          
          ! --- Pressure
          P0       = r0    * T0
          P0_R     = r0_R  * T0 + r0 * T0_R
          P0_Z     = r0_Z  * T0 + r0 * T0_Z
          P0_s     = r0_s  * T0 + r0 * T0_s
          P0_t     = r0_t  * T0 + r0 * T0_t
          P0_p     = r0_p  * T0 + r0 * T0_p
          P0_pp    = r0_pp * T0 + r0 * T0_pp + 2.d0 * r0_p * T0_p
          P0_RR    = r0_RR * T0 + r0 * T0_RR + 2.d0 * r0_R * T0_R
          P0_ZZ    = r0_ZZ * T0 + r0 * T0_ZZ + 2.d0 * r0_Z * T0_Z
          P0_RZ    = r0_RZ * T0 + r0 * T0_RZ + r0_R * T0_Z + r0_Z * T0_R
          
          ! --- Some things related to the magnetic field
          BB2      = (F0*F0 + ps0_R * ps0_R + ps0_Z * ps0_Z ) / BigR**2
          B_R      = + ps0_Z / BigR
          B_Z      = - ps0_R / BigR
          B_tor    = + F0    / BigR
          psi_norm = get_psi_n(eq, ps0)
          Btheta  = sqrt(ps0_R*ps0_R + ps0_Z * ps0_Z) / BigR
          psi_abs = sqrt(ps0_R*ps0_R + ps0_Z * ps0_Z)

          if (psi_abs > 1.d-6) then
            FFprime_loc = zj0 + (R**2.d0) * (ps0_R*P0_R + ps0_Z*P0_Z)/(psi_abs**2.d0)
          else
            FFprime_loc = zj0 !--- not fully correct, but better than to put 0...
          endif
          Jpol = FFprime_loc * Btheta

          ! --- Some input profiles
          if ( eta_T_dependent ) then
            eta_T     = eta   * (corr_neg_temp(T0)/T_0)**(-1.5d0)
            deta_dT   = - eta   * (1.5d0)  * corr_neg_temp(T0)**(-2.5d0) * T_0**(1.5d0)
            d2eta_d2T =   eta   * (3.75d0) * corr_neg_temp(T0)**(-3.5d0) * T_0**(1.5d0)
          else
            eta_T     = eta
            deta_dT   = 0.d0
            d2eta_d2T = 0.d0
          end if
          
          if ( visco_T_dependent ) then
            visco_T   = visco * (corr_neg_temp(T0)/T_0)**(-1.5d0)
            dvisco_dT = - visco * (1.5d0)  * corr_neg_temp(T0)**(-2.5d0) * T_0**(1.5d0)
          else
            visco_T   = visco
            dvisco_dT = 0.d0
          end if
          
          if ( ZKpar_T_dependent ) then
            ZKpar_T   = ZK_par * (corr_neg_temp(T0)/T_0)**(+2.5d0)
            dZKpar_dT = ZK_par * (2.5d0)  * corr_neg_temp(T0)**(+1.5d0) * T_0**(-2.5d0)
            if (ZKpar_T .gt. ZK_par_max) then
              ZKpar_T   = Zk_par_max
              dZKpar_dT = 0.d0
            end if
          else
            ZKpar_T   = ZK_par
            dZKpar_dT = 0.d0
          end if
          
          D_prof  = get_dperp (psi_norm)
          ZK_prof = get_zkperp(psi_norm)
          
          ! --- Other parameters (combination of the main variables)
          Er       = 0.d0
          Vtheta   = 0.d0
	  mach_par = 0.d0
          mach_pol = 0.d0
          vsound   = 0.d0
          Vneo     = 0.d0
          ki_neo   = 0.d0
          mu_neo   = 0.d0
          Vperp_e  = 0.d0
          Vperp_i  = 0.d0
          V_ExB    = 0.d0 
          Vstar_e  = 0.d0
          Vstar_i  = 0.d0
          
          if ( (psi_abs > 1.d-6) .and. (r0 > 1.d-6) .and. (abs(Btheta) > 1.d-6) ) then
            
            Er       = -(u0_R * ps0_R + u0_Z * ps0_Z) / psi_abs   ! radial electric field
            
            Vsound   = sqrt(GAMMA*T0) / sqrt(BB2)                 ! sound speed
            Mach_par = Vpar0 / Vsound                             ! parallel Mach number
            Mach_pol = Vtheta / Vsound                            ! poloidal Mach number
            
            Vtheta   = -1./Btheta * (  ( u0_R + tauIC/r0 * (T0_R*r0 + r0_R*T0) ) * ps0_R  +        &
              ( u0_Z + tauIC/r0 * (T0_Z*r0 + r0_Z*T0) ) * ps0_Z) + Vpar0 * Btheta
            
            Vperp_i  = -1./Btheta * (  ( u0_R + tauIC/r0 * (T0_R*r0 + r0_R*T0) ) * ps0_R  +        &
              ( u0_Z + tauIC/r0 * (T0_Z*r0 + r0_Z*T0) ) * ps0_Z )
            
            Vperp_e  = -1./Btheta * (  ( u0_R - tauIC/r0 * (T0_R*r0 + r0_R*T0) ) * ps0_R  +        &
              ( u0_Z - tauIC/r0 * (T0_Z*r0 + r0_Z*T0) ) * ps0_Z )
            
            V_ExB    = -1./Btheta* ( u0_R*ps0_R + u0_Z*ps0_Z )
            
            Vstar_i  = -1./Btheta * (  tauIC/r0 * (T0_R*r0 + r0_R*T0) * ps0_R  +                   &
              tauIC/r0 * (T0_Z*r0 + r0_Z*T0) * ps0_Z )
            
            Vstar_e  = +1./Btheta * (  tauIC/r0 * (T0_R*r0 + r0_R*T0) * ps0_R  +                   &
              tauIC/r0 * (T0_Z*r0 + r0_Z*T0) * ps0_Z )
            ! ### Warning : in jorek_model=400, Vstar_i .ne. -Vstar_e since T_i .ne. T_e
          end if
          
          if (NEO) then
            if (num_neo_file) then ! (read neoclassical profiles from ascii file)
              call neo_coef( eq%xpoint, eq%xcase, Z, eq%Z_xpoint, Ps0 ,eq%psi_axis, eq%psi_bnd,    &
                mu_neo, ki_neo)
              Vneo = ki_neo / Btheta * tauIC  * ( ps0_R*T0_R + ps0_Z*T0_Z )
            else ! (use constant neoclassical coefficients from namelist input file)
              mu_neo = amu_neo_const
              ki_neo = aki_neo_const
              Vneo   = aki_neo_const / Btheta * tauIC * ( ps0_R*T0_R + ps0_Z*T0_Z )
            end if
          end if
          
          ln_Lambda = 18. ! approximate value for Coulomb logarithm
          
          E_par = - F0 / R * ( eta_T * zj0 / R**2                                                  &
                             + tauIC / r0 * ( (P0_R * Ps0_Z - P0_Z * Ps0_R) / R )                  &
                             + F0 * P0_p / R**2 )
          
          E_crit = C_LIGHT**2 * EL_CHG**3 * ln_Lambda * MU_ZERO**2.5 ** rho_0**1.5 * r0 / ( 4 * PI * MASS_ELECTRON * MASS_PROTON * central_mass )
          
#if JOREK_MODEL == 303 || JOREK_MODEL == 333 || JOREK_MODEL == 400 || JOREK_MODEL == 500
          call bootstrap_current(R, Z, eq%R_axis, eq%Z_axis, eq%psi_axis, eq%R_xpoint, eq%Z_xpoint, eq%psi_bnd, psi_norm, ps0, ps0_R,    &
            ps0_Z, r0,  r0_R, r0_Z, Ti0, Ti0_R, Ti0_Z, Te0, Te0_R, Te0_Z, J_boot)
#else
          J_boot = 0.d0
#endif

#if JOREK_MODEL == 500

   T_rad = corr_neg_temp(T0)/(2.d0*EL_CHG*MU_ZERO*central_density * 1.d20)
  !write(*,*) 'T_rad = ', T_rad
  if ( units == SI_UNITS ) then

   coef_rad_1 = 1.d0

  else if ( units == JOREK_UNITS ) then

   coef_rad_1 = 2.d0/(3.d0)*MU_ZERO**1.5d0*(central_mass*MASS_PROTON)**0.5d0*(central_density * 1.d20)**2.5d0

  endif

   LradDcont_T = coef_rad_1*5.37d-37*(1.d1)**(-1.5d0)*(1.d0)**2*sqrt(T_rad) ! Only Bremsstrahlung contribution

   LradDrays_T = coef_rad_1*(1.d1)**(-29.44d0*exp(-(log10(T_rad)-4.4283d0)**2.d0/(2.d0*(2.8428d0)**2.d0)) &
                                    -60.947d0*exp(-(log10(T_rad)+2.0835d0)**2.d0/(2.d0*(0.9048d0)**2.d0)) &
                                    -24.067d0*exp(-(log10(T_rad)+0.7363d0)**2.d0/(2.d0*(2.1700d0)**2.d0)))

 !write(*,*) 'Lbrem = ', LradDcont_T
 !write(*,*) 'Lrays = ', LradDrays_T

  !--------------------------------------------------------
  ! --- Radiation from background impurity
  !--------------------------------------------------------

    Arad_bg = 2.4d-31
    Brad_bg = 20.
    Crad_bg = 0.8

  if ( units == SI_UNITS ) then

    frad_bg = nimp_bg*Arad_bg*exp(-((log(T_rad)-log(Brad_bg))**2.)/Crad_bg**2.)

  else if ( units == JOREK_UNITS ) then

    frad_bg = (2./3.)*(1./(central_mass*MASS_PROTON))*((MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)**(1.5d0))*nimp_bg*Arad_bg*exp(-((log(T_rad)-log(Brad_bg))**2.)/Crad_bg**2.)

  endif
  !--------------------------------------------------------

#endif


          ! --- Factors for switching between JOREK normalized and SI units.
          if ( units == SI_UNITS ) then
             rho_norm      = central_density *1.d20 * central_mass * mass_proton   ! rho_0 = central mass density
             fact_time     = sqrt(MU_zero*rho_norm)                                ! time factor
             fact_mu_zero  = MU_zero                                               ! division by mu_zero for P and J
             fact_ne       = central_density * 1.d20                               ! factor for n_e
             fact_rho      = central_density * 1.d20 * central_mass*MASS_PROTON    ! factor for rho
             fact_T        = 1.d0 / ( MU_zero * central_density * 1.d20 * EL_CHG ) ! factor for T
             fact_vpar     = sqrt(BB2) / fact_time                                 ! factor for Vpar
             fact_resistiv = sqrt ( MU_zero / rho_norm )                           ! factor for eta == 1 / (factor for visco)
             fact_Er       = F0 / fact_time
          else if ( units == JOREK_UNITS ) then
             fact_time     = 1.d0
             fact_mu_zero  = 1.d0
             fact_ne       = 1.d0
             fact_rho      = 1.d0
             fact_T        = 1.d0
             fact_vpar     = 1.d0
             fact_resistiv = 1.d0
             fact_Er       = 1.d0
          end if
          
          ! --- Now that everything is prepared, evaluate the requested expressions.
          loop_expr: do iexpr = 1, expr_list%n_expr
            
            select case ( trim(expr_list%expr(iexpr)%name) )
              case ( 'R' )
                res = R
                
              case ( 'Z', 'z' )
                res = Z
                
              case ( 'phi' )
                res = phi
                
              case ( 'theta' )
                res = theta
                
              case ( 'theta_star' )
                res = pol_pos%theta_star
                
              case ( 'length' )
                res = pol_pos%length
                
              case ( 'r_minor' )
                res = pol_pos%r_minor
                
              case ( 'xjac' )
                res = xjac
                
              case ( 'x' )
                res = x_cart
                
              case ( 'y' )
                res = y_cart
                
              case ( 't' )
                res = t_now * fact_time
                
              case ( 'Psi' )
                res = ps0
                
              case ( 'Psi_N' )
                res = psi_norm
                
              case ( 'u' )
                res = u0
                
              case ( 'Phi' )
                res = u0 * F0 !### sign?
                
              case ( 'zj' )
                res = zj0 / fact_mu_zero
                
              case ( 'omega' )
                res = w0
                
              case ( 'rho' )
                res = r0 * fact_rho
                
              case ( 'ne' )
                res = r0 * fact_ne
                
              case ( 'T' )
                res = T0 * fact_T
              
              case ( 'Te' )
                res = T0 * fact_T / 2.d0
              
              case ( 'vpar' )
                res = Vpar0 * fact_vpar
                
              case ( 'eta_T' )
                res = eta_t * fact_resistiv
                
              case ( 'visco_T' )
                res = visco_t / fact_resistiv
                
              case ( 'zkpar_T' )
                res = zkpar_t / fact_time
                
              case ( 'dprof' ) 
                res = d_prof / fact_time
                
              case ( 'zkprof' )
                res = zk_prof / fact_time
                
              case ( 'pres' )
                res = P0 / fact_mu_zero
                
              case ( 'B_abs' )
                res = sqrt(BB2)
                
              case ( 'B_tor' )
                res = B_tor
                
              case ( 'B_R' )
                res = B_R
                
              case ( 'B_Z' )
                res = B_Z
                
              case ( 'B_theta' )
                res = Btheta
                
              case ( 'currdens' )
                res = zj0 / R / fact_mu_zero

              case ( 'FFprime_loc' )
                res = FFprime_loc

              case ( 'Jpol' )
                res = Jpol / fact_mu_zero
                
              case ( 'Er' )
                res = Er * fact_Er
                
              case ( 'Vtheta_i' )
                res = Vtheta / fact_time
                
              case ( 'Mach_par' )
                res = Mach_par
                
              case ( 'Mach_pol' )
                res = Mach_pol
                
              case ( 'V_sound' )
                res = Vsound / fact_time
                
              case ( 'V_neo' )
                res = Vneo / fact_time
                
              case ( 'Vperp_e' )
                res = Vperp_e / fact_time
                
              case ( 'Vperp_i' )
                res = Vperp_i / fact_time
                
              case ( 'V_ExB' )
                res = V_ExB / fact_time
                
              case ( 'Vstar_e' )
                res = Vstar_e / fact_time
                
              case ( 'Vstar_i' )
                res = Vstar_i / fact_time
                
              case ( 'ki_neo' )
                res = ki_neo
                
              case ( 'mu_neo' )
                res = mu_neo / fact_time
                
              case ( 'T_e' )
                res = Te0 * fact_T
                
              case ( 'T_i' )
                res = Ti0 * fact_T
                
              case ( 'E_||' )
                res = E_par / fact_time
                
              case ( 'E_crit' )
                res = E_crit / fact_time
                
#if JOREK_MODEL == 303 || JOREK_MODEL == 333 || JOREK_MODEL == 400 || JOREK_MODEL == 500
              case ( 'J_bootstrap' )
                res = J_boot ! ### check if no normalization needed
#endif

#if JOREK_MODEL == 500
              case ( 'radiation' )

                if (rn0 .lt. 0.d0) then
                  res = r0 * fact_ne * r0 * fact_ne * LradDcont_T &
                       + r0 * fact_ne * frad_bg
                else
                  res = r0 * fact_ne * rn0 * fact_ne * LradDrays_T &
                       + r0 * fact_ne * r0 * fact_ne * LradDcont_T &
                       + r0 * fact_ne * frad_bg
                endif

              case ( 'brem' )
                res = r0 * fact_ne * r0 * fact_ne * LradDcont_T
#endif

              case default
                write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': Illegal expression ("' //      &
                  trim(expr_list%expr(iexpr)%name) // '")'
                ierr = 100
                return
                
            end select
            
            result(itorpos, ipolpos, jpolpos, iexpr) = res
            
          end do loop_expr
        end do loop_tor
      end do loop_pol2
    end do loop_pol1
    
  end subroutine eval_expr
  
  
  
  
  
end module mod_expression
