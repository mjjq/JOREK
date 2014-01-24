!> This module allows to calculate "arbitrary" expressions at any position or list of positions in
!! the computational domain for diagnostic purposes.
!!
!! To obtain good performance, evaluate several expressions at many toroidal positions.
module mod_expression
  
  
  
  
  
  use parameters
  use mod_position
  use phys_module
  use diffusivities
  use corr_neg
  
  
  
  
  
  implicit none
  
  
  
  
  
  public
  
  
  
  
  
  ! --- Constants
  character(len=14), parameter, private :: THIS_MOD_NAME = 'mod_expression'
  
  ! --- Constants specifying available expressions
  !   --- Coordinates etc.
  integer, parameter :: EXPR_R        = 001001
  integer, parameter :: EXPR_Z        = 001002
  integer, parameter :: EXPR_PHI      = 001003
  integer, parameter :: EXPR_XJAC     = 001004
  
  !   --- Basic variables
  integer, parameter :: EXPR_PSI      = 002001
  integer, parameter :: EXPR_PSIN     = 002002
  integer, parameter :: EXPR_U        = 002003
  integer, parameter :: EXPR_ZJ       = 002004
  integer, parameter :: EXPR_W        = 002005
  integer, parameter :: EXPR_RHO      = 002006
  integer, parameter :: EXPR_T        = 002007
  integer, parameter :: EXPR_VPAR     = 002008
  
  !   --- Input profiles
  integer, parameter :: EXPR_ETAT     = 003001
  integer, parameter :: EXPR_VISCOT   = 003002
  integer, parameter :: EXPR_ZKPART   = 003003
  integer, parameter :: EXPR_DPROF    = 003004
  integer, parameter :: EXPR_ZKPROF   = 003005
  
  !   --- Derived variables
  integer, parameter :: EXPR_PRES     = 004001
  integer, parameter :: EXPR_BABS     = 004002
  integer, parameter :: EXPR_BTOR     = 004003
  integer, parameter :: EXPR_BR       = 004004
  integer, parameter :: EXPR_BZ       = 004005
  integer, parameter :: EXPR_CURRDENS = 004006
  
  
  
  
  
  contains
  
  
  
  
  
  !> Evaluate one or several expressions at one or several poloidal and one or several toroidal
  !! positions.
  subroutine eval_expr(eq, si_units, expr, n_expr, pol_pos_list, tor_pos_list, result, ierr)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':eval_expr'
    
    ! --- Routine parameters
    type(t_equil_state),          intent(in)    :: eq
    logical,                      intent(in)    :: si_units !< Output SI or JOREK normalized units?
    integer,                      intent(in)    :: expr(n_expr)
    integer,                      intent(in)    :: n_expr
    type(t_pol_pos_list), target, intent(in)    :: pol_pos_list
    type(t_tor_pos_list), target, intent(in)    :: tor_pos_list
    real*8, allocatable,          intent(inout) :: result(:,:,:)
    integer,                      intent(out)   :: ierr
    
    ! --- Local variables
    type(t_pol_pos), pointer :: pol_pos
    type(t_tor_pos), pointer :: tor_pos
    type(type_element)       :: element
    type(type_node)          :: nodes(n_vertex_max)
    integer :: ipolpos, itorpos, iexpr, ielm, i, j, i_tor
    real*8  :: xjac, xjac_R, xjac_Z, R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, &
      s, t, H(n_vertex_max,n_order+1), H_s(n_vertex_max,n_order+1), H_t(n_vertex_max,n_order+1),   &
      H_st(n_vertex_max,n_order+1), H_ss(n_vertex_max,n_order+1), H_tt(n_vertex_max,n_order+1),    &
      HZ(n_tor), HZ_p(n_tor), HZ_pp(n_tor), phi, res, BigR, BigR_R
    real*8  :: ps0, ps0_s, ps0_t, ps0_ss, ps0_tt, ps0_st, ps0_p, ps0_pp, u0, u0_s, u0_t, u0_ss,    &
      u0_tt, u0_st, u0_p, u0_pp, zj0, zj0_s, zj0_t, zj0_ss, zj0_tt, zj0_st, zj0_p, zj0_pp, w0,     &
      w0_s, w0_t, w0_ss, w0_tt, w0_st, w0_p, w0_pp, r0, r0_s, r0_t, r0_ss, r0_tt, r0_st, r0_p,     &
      r0_pp, T0, T0_s, T0_t, T0_ss, T0_tt, T0_st, T0_p, T0_pp, Vpar0, Vpar0_s, Vpar0_t, Vpar0_ss,  &
      Vpar0_tt, Vpar0_st, Vpar0_p, Vpar0_pp, psi_norm
    real*8  :: ps0_R, ps0_Z, ps0_RR, ps0_ZZ, ps0_RZ, u0_R, u0_Z, u0_RR, u0_ZZ, u0_RZ, vv2, zj0_R,  &
      zj0_Z, zj0_RR, zj0_ZZ, zj0_RZ, w0_R, w0_Z, w0_RR, w0_ZZ, w0_RZ, r0_R, r0_Z, r0_RR, r0_ZZ,    &
      r0_RZ, r0_hat, r0_R_hat, r0_Z_hat, T0_R, T0_Z, T0_RR, T0_ZZ, T0_RZ, T0_ps0_R, T0_ps0_Z,      &
      Vpar0_R, Vpar0_Z, Vpar0_RR, Vpar0_ZZ, Vpar0_RZ, P0, P0_R, P0_Z, P0_s, P0_t, P0_p, P0_pp,     &
      P0_RR, P0_ZZ, P0_RZ, BB2, B_tor, B_R, B_Z
    real*8  :: eta_T, deta_dT, d2eta_d2T, visco_T, dvisco_dT, ZKpar_T, dZKpar_dT, D_prof, ZK_prof
    
    ierr = 0
    
    ! --- Some sanity checks
    if ( n_expr < 1 ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': n_expr < 1 encountered.'
      ierr = -101
      return
    else if ( pol_pos_list%n_pos < 1 ) then
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
    allocate( result(n_expr, pol_pos_list%n_pos, tor_pos_list%n_pos) )
    
    ! --- Loop over different positions in the poloidal plane
    do ipolpos = 1, pol_pos_list%n_pos
      pol_pos => pol_pos_list%pos(ipolpos)
      
      ! --- Coordinates, Jacobian, Basis Functions
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
      xjac   = R_s * Z_t - R_t * Z_s
      xjac_R = ( R_ss * Z_t**2 - 2*R_st * Z_s*Z_t + R_tt * Z_s**2 + R_s * (Z_st*Z_t - Z_tt*Z_s )   &
               + R_t * (Z_st*Z_s - Z_ss*Z_t ) ) / xjac
      xjac_Z = ( Z_ss * R_t**2 - 2*Z_st * R_s*R_t + Z_tt * R_s**2 + Z_s * (R_st*R_t - R_tt*R_s )   &
               + Z_t * (R_st*R_s - R_ss*R_t ) ) / xjac
      call basisfunctions2(s, t, H, H_s, H_t, H_st, H_ss, H_tt)
      BigR   = R
      BigR_R = 1.d0
      
      ! --- Elements and nodes
      element  = pol_pos%element
      nodes(:) = pol_pos%nodes(:)
      
      ! ### Loop over nodes and dofs
      
      ! --- Loop over different toroidal positions
      do itorpos = 1, tor_pos_list%n_pos
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
        
        ps0   = 0.d0; ps0_s   = 0.d0; ps0_t   = 0.d0; ps0_ss   = 0.d0; ps0_tt   = 0.d0; ps0_st   = 0.d0; ps0_p   = 0.d0; ps0_pp   = 0.d0
        u0    = 0.d0; u0_s    = 0.d0; u0_t    = 0.d0; u0_ss    = 0.d0; u0_tt    = 0.d0; u0_st    = 0.d0; u0_p    = 0.d0; u0_pp    = 0.d0
        zj0   = 0.d0; zj0_s   = 0.d0; zj0_t   = 0.d0; zj0_ss   = 0.d0; zj0_tt   = 0.d0; zj0_st   = 0.d0; zj0_p   = 0.d0; zj0_pp   = 0.d0
        w0    = 0.d0; w0_s    = 0.d0; w0_t    = 0.d0; w0_ss    = 0.d0; w0_tt    = 0.d0; w0_st    = 0.d0; w0_p    = 0.d0; w0_pp    = 0.d0
        r0    = 0.d0; r0_s    = 0.d0; r0_t    = 0.d0; r0_ss    = 0.d0; r0_tt    = 0.d0; r0_st    = 0.d0; r0_p    = 0.d0; r0_pp    = 0.d0
        T0    = 0.d0; T0_s    = 0.d0; T0_t    = 0.d0; T0_ss    = 0.d0; T0_tt    = 0.d0; T0_st    = 0.d0; T0_p    = 0.d0; T0_pp    = 0.d0
        Vpar0 = 0.d0; Vpar0_s = 0.d0; Vpar0_t = 0.d0; Vpar0_ss = 0.d0; Vpar0_tt = 0.d0; Vpar0_st = 0.d0; Vpar0_p = 0.d0; Vpar0_pp = 0.d0
        
        ! --- Integrate
        do i = 1, n_vertex_max
          do j = 1, n_order+1
            do i_tor = 1, n_tor
              
              ! --- Poloidal Flux
              ps0      = ps0      + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H   (i,j) * HZ   (i_tor)
              ps0_s    = ps0_s    + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_s (i,j) * HZ   (i_tor)
              ps0_t    = ps0_t    + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_t (i,j) * HZ   (i_tor)
              ps0_ss   = ps0_ss   + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_ss(i,j) * HZ   (i_tor)
              ps0_tt   = ps0_tt   + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_tt(i,j) * HZ   (i_tor)
              ps0_st   = ps0_st   + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_st(i,j) * HZ   (i_tor)
              ps0_p    = ps0_p    + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H   (i,j) * HZ_p (i_tor)
              ps0_pp   = ps0_pp   + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H   (i,j) * HZ_pp(i_tor)
              
              ! --- Stream Function
              u0       = u0       + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H   (i,j) * HZ   (i_tor)
              u0_s     = u0_s     + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_s (i,j) * HZ   (i_tor)
              u0_t     = u0_t     + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_t (i,j) * HZ   (i_tor)
              u0_ss    = u0_ss    + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_ss(i,j) * HZ   (i_tor)
              u0_tt    = u0_tt    + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_tt(i,j) * HZ   (i_tor)
              u0_st    = u0_st    + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_st(i,j) * HZ   (i_tor)
              u0_p     = u0_p     + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H   (i,j) * HZ_p (i_tor)
              u0_pp    = u0_pp    + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H   (i,j) * HZ_pp(i_tor)
              
              ! --- Current
              zj0      = zj0      + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H   (i,j) * HZ   (i_tor)
              zj0_s    = zj0_s    + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_s (i,j) * HZ   (i_tor)
              zj0_t    = zj0_t    + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_t (i,j) * HZ   (i_tor)
              zj0_ss   = zj0_ss   + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_ss(i,j) * HZ   (i_tor)
              zj0_tt   = zj0_tt   + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_tt(i,j) * HZ   (i_tor)
              zj0_st   = zj0_st   + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_st(i,j) * HZ   (i_tor)
              zj0_p    = zj0_p    + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H   (i,j) * HZ_p (i_tor)
              zj0_pp   = zj0_pp   + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H   (i,j) * HZ_pp(i_tor)
              
              ! --- Vorticity
              w0       = w0       + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H   (i,j) * HZ   (i_tor)
              w0_s     = w0_s     + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_s (i,j) * HZ   (i_tor)
              w0_t     = w0_t     + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_t (i,j) * HZ   (i_tor)
              w0_ss    = w0_ss    + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_ss(i,j) * HZ   (i_tor)
              w0_tt    = w0_tt    + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_tt(i,j) * HZ   (i_tor)
              w0_st    = w0_st    + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_st(i,j) * HZ   (i_tor)
              w0_p     = w0_p     + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H   (i,j) * HZ_p (i_tor)
              w0_pp    = w0_pp    + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H   (i,j) * HZ_pp(i_tor)
              
              ! --- Density
              r0       = r0       + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H   (i,j) * HZ   (i_tor)
              r0_s     = r0_s     + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_s (i,j) * HZ   (i_tor)
              r0_t     = r0_t     + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_t (i,j) * HZ   (i_tor)
              r0_ss    = r0_ss    + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_ss(i,j) * HZ   (i_tor)
              r0_tt    = r0_tt    + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_tt(i,j) * HZ   (i_tor)
              r0_st    = r0_st    + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_st(i,j) * HZ   (i_tor)
              r0_p     = r0_p     + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H   (i,j) * HZ_p (i_tor)
              r0_pp    = r0_pp    + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H   (i,j) * HZ_pp(i_tor)
              
              ! --- Temperature
              T0       = T0       + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j) * HZ   (i_tor)
              T0_s     = T0_s     + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_s (i,j) * HZ   (i_tor)
              T0_t     = T0_t     + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_t (i,j) * HZ   (i_tor)
              T0_ss    = T0_ss    + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_ss(i,j) * HZ   (i_tor)
              T0_tt    = T0_tt    + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_tt(i,j) * HZ   (i_tor)
              T0_st    = T0_st    + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_st(i,j) * HZ   (i_tor)
              T0_p     = T0_p     + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j) * HZ_p (i_tor)
              T0_pp    = T0_pp    + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j) * HZ_pp(i_tor)
              
              ! --- Parallel Velocity
              if ( jorek_model >= 300 ) then
                Vpar0    = Vpar0    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j) * HZ   (i_tor)
                Vpar0_s  = Vpar0_s  + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_s (i,j) * HZ   (i_tor)
                Vpar0_t  = Vpar0_t  + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_t (i,j) * HZ   (i_tor)
                Vpar0_ss = Vpar0_ss + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_ss(i,j) * HZ   (i_tor)
                Vpar0_tt = Vpar0_tt + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_tt(i,j) * HZ   (i_tor)
                Vpar0_st = Vpar0_st + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_st(i,j) * HZ   (i_tor)
                Vpar0_p  = Vpar0_p  + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j) * HZ_p (i_tor)
                Vpar0_pp = Vpar0_pp + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j) * HZ_pp(i_tor)
              end if
              
              ! --- Deltas
              !do k=1,n_var
              !  delta_g(k) = delta_g(k) + nodes(i)%deltas(i_tor,j,k) * element%size(i,j) * H   (i,j) * HZ   (i_tor)
              !  delta_s(k) = delta_s(k) + nodes(i)%deltas(i_tor,j,k) * element%size(i,j) * H_s (i,j) * HZ   (i_tor)
              !  delta_t(k) = delta_t(k) + nodes(i)%deltas(i_tor,j,k) * element%size(i,j) * H_t (i,j) * HZ   (i_tor)
              !enddo                        
              
            enddo
          enddo
        enddo
        
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
        
        BB2      = (F0*F0 + ps0_R * ps0_R + ps0_Z * ps0_Z ) / BigR**2
        B_R      = ps0_Z / BigR
        B_Z      = ps0_R / BigR
        B_tor    = F0    / BigR
        
        psi_norm = get_psi_n(eq, ps0)
        
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
          endif
        else
          ZKpar_T   = ZK_par
          dZKpar_dT = 0.d0
        endif
        
        D_prof  = get_dperp (psi_norm)
        ZK_prof = get_zkperp(psi_norm)
        
        do iexpr = 1, n_expr
          
          select case ( expr(iexpr) )
            
            case ( EXPR_R )
              res = R
              
            case ( EXPR_Z )
              res = Z
              
            case ( EXPR_PHI )
              res = phi
              
            case ( EXPR_XJAC )
              res = xjac
              
            case ( EXPR_PSI )
              res = ps0
              
            case ( EXPR_PSIN )
              res = psi_norm
              
            case ( EXPR_U )
              res = u0
              
            case ( EXPR_ZJ )
              res = zj0
              
            case ( EXPR_W )
              res = w0
              
            case ( EXPR_RHO )
              res = r0
              
            case ( EXPR_T )
              res = T0
              
            case ( EXPR_VPAR )
              res = Vpar0
              
            case ( EXPR_ETAT )
              res = eta_t
              
            case ( EXPR_VISCOT )
              res = visco_t
              
            case ( EXPR_ZKPART )
              res = zkpar_t
              
            case ( EXPR_DPROF )
              res = d_prof
              
            case ( EXPR_ZKPROF )
              res = zk_prof
              
            case ( EXPR_PRES )
              res = P0
              
            case ( EXPR_BABS )
              res = sqrt(BB2)
              
            case ( EXPR_BTOR )
              res = B_tor
              
            case ( EXPR_BR )
              res = B_R
              
            case ( EXPR_BZ )
              res = B_Z
              
            case ( EXPR_CURRDENS )
              res = zj0 * R !###
              
            case default
              write(*,*) 'WARNING in '//trim(THIS_ROUTINE_NAME)//': Illegal expression number.'
              ierr = 100
              res  = 0.d0 / 0.d0
              
          end select
          
          result(iexpr, ipolpos, itorpos) = res
          
        end do
      end do
    end do
    
  end subroutine eval_expr
  
  
  
  
  
end module mod_expression
