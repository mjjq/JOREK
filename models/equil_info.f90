!> Defines a data structure to collect all "equilibrium state information" of the simulation.
!!
!! This includes the location of magnetic axis, X-points, limiter point etc.
!!
!! The module also contains routines/functions for tasks like Psi_N calculation etc.
!!
!! @TODO: Better calculation of strike points (just place-holder right now) including i_elm_strike
module equil_info
  
  
  
  use constants,      only: LOWER_XPOINT, UPPER_XPOINT, DOUBLE_NULL
  use data_structure, only: type_node_list, type_element_list, type_bnd_element_list
  
  
  
  implicit none
  
  
  
  public
  
  
  
  !> This data structure contains information about the position of axis, limiter point,
  !! X-point(s), and strike points. It needs to be updated after each time step and after
  !! each iteration of the equilibrium calculation using the routine update_equil_state.
  type t_equil_state
    
    logical          :: initialized = .false.
    
    logical          :: limiter_plasma           !< Is the current state a limiter plasma?
    logical          :: axis_is_psi_minimum      !< Is psi_axis < psi_bnd or > psi_bnd?
    real*8           :: Psi_bnd                  !< Psi of limiting surface or separatrix.
    
    ! --- Magnetic Axis
    real*8           :: R_axis                   !< R coordinate of axis.
    real*8           :: Z_axis                   !< Z coordinate of axis.
    real*8           :: Psi_axis                 !< Poloidal flux value at axis.
    integer          :: i_elm_axis               !< Index of element containing the axis.
    real*8           :: s_axis                   !< s coordinate of axis within element.
    real*8           :: t_axis                   !< t coordinate of axis within element.
    integer          :: ifail_axis               !< Error code for axis determination.
    
    ! --- Limiter Point
    real*8           :: R_lim                    !< R coordinate of limiter point.
    real*8           :: Z_lim                    !< Z coordinate of limiter point.
    real*8           :: Psi_lim                  !< Poloidal flux value at limiter point.
    integer          :: i_elm_lim                !< Index of element containing limiter point.
    real*8           :: s_lim                    !< s coordinate of limiter point within element.
    real*8           :: t_lim                    !< t coordinate of limiter point within element.
    integer          :: ifail_lim                !< Error code for limiter determination.
    
    ! --- X-Point(s)
    logical          :: xpoint                   !< Is this an X-point case? Not necessarily active!
    integer          :: xcase                    !< Upper/lower/double X-point?
    integer          :: active_xpoint            !< Which X-point is active?
    real*8           :: R_xpoint(2)              !< R coordinate of X-point(s).
    real*8           :: Z_xpoint(2)              !< Z coordinate of X-point(s).
    real*8           :: Psi_xpoint(2)            !< Poloidal flux value of X-point(s).
    integer          :: i_elm_xpoint(2)          !< Index of element containing the X-point.
    real*8           :: s_xpoint(2)              !< s coordinate of X-point within element.
    real*8           :: t_xpoint(2)              !< t coordinate of X-point within element.
    integer          :: ifail_xpoint             !< Error code for X-point determination.
    
    ! --- Strike Point(s) derived from axisymmetric field component.
    integer          :: num_strike               !< Number of strike points.
    real*8           :: R_strike(99)             !< R coordinate of strike point(s).
    real*8           :: Z_strike(99)             !< Z coordinate of strike point(s).
    integer          :: i_bndelm_strike(99)      !< Index of boundary element containing strike pt.
    real*8           :: s_strike(99)             !< s coordinate of strike pt within boundary elem.
    
    ! --- Inner/Outer points on the midplane close to the boundary of the computational domain.
    real*8           :: R_midpl(2)               !< R coordinate of "midplane points".
    
  end type t_equil_state
  
  
  
  contains
  
  
  
  !> Re-calculate the equilibrium state.
  subroutine update_equil_state(node_list, element_list, bnd_elm_list, xpoint, xcase, ES)
    
    ! --- Routine parameters.
    type(type_node_list),        intent(in)    :: node_list
    type(type_element_list),     intent(in)    :: element_list
    type(type_bnd_element_list), intent(in)    :: bnd_elm_list
    logical                                    :: xpoint
    integer                                    :: xcase
    type(t_equil_state),         intent(inout) :: ES
    
    ! --- Local variables.
    integer :: my_id, i_out, ifail, i
    real*8  :: R_out, Z_out, s_out, t_out, corr_fact, R1, R2, dR
    
    my_id  = 9999
    
    ! --- Find the magnetic axis.
    call find_axis(my_id, node_list, element_list, ES%psi_axis, ES%R_axis, ES%Z_axis,              &
      ES%i_elm_axis, ES%s_axis, ES%t_axis, ES%ifail_axis)
    
    ! --- Find the X-point(s).
    ES%xpoint       = xpoint
    ES%xcase        = xcase
    ES%ifail_xpoint = 0
    if ( xpoint ) call find_xpoint(my_id, node_list, element_list, ES%psi_xpoint, ES%R_xpoint,     &
      ES%Z_xpoint, ES%i_elm_xpoint, ES%s_xpoint, ES%t_xpoint, ES%xcase, ES%ifail_xpoint)
    
    ! --- Find the limiter point.
    ES%ifail_lim = 0
    call find_limiter(my_id, node_list, element_list, bnd_elm_list, ES%psi_lim, ES%R_lim, ES%Z_lim)
    call find_RZ(node_list, element_list, ES%R_lim, ES%Z_lim, R_out, Z_out, ES%i_elm_lim, ES%s_lim,&
      ES%t_lim, ES%ifail_lim)
    
    ! --- corr_fact such that: psi_axis * corr_fact < psi_bnd * corr_fact
    !    (to account for cases where Psi takes minimum/maximum value at the axis)
    corr_fact = 1.d0
    if ( ES%psi_axis > ES%psi_lim ) corr_fact = -1.d0
    
    if ( xpoint ) then ! (X-point plasma)
      
      if ( (xcase==LOWER_XPOINT) ) then
        
        ES%psi_bnd        = corr_fact * ES%psi_xpoint(1)
        ES%limiter_plasma = .false.
        ES%active_xpoint  = LOWER_XPOINT
        
      else if ( (xcase==UPPER_XPOINT) ) then
        
        ES%psi_bnd        = corr_fact * ES%psi_xpoint(2)
        ES%limiter_plasma = .false.
        ES%active_xpoint  = UPPER_XPOINT
        
      else if ( (xcase==DOUBLE_NULL) ) then
        
        ES%limiter_plasma = .false.
        
        if ( corr_fact*ES%psi_xpoint(1) < corr_fact*ES%psi_xpoint(2) ) then
          ES%psi_bnd       = corr_fact * ES%psi_xpoint(1)
          ES%active_xpoint = LOWER_XPOINT
        else
          ES%psi_bnd       = corr_fact * ES%psi_xpoint(2)
          ES%active_xpoint = UPPER_XPOINT
        end if
        
      else ! This should never happen.
        write(*,*) 'ERROR: ILLEGAL VALUE FOR XCASE:', xcase
        stop
      end if
      
      ! ### We need to take into account the possibility of a limiter plasma here although
      !     we have an X-point, e.g., during a VDE
      
    else ! (limiter plasma)
      
      ES%psi_bnd        = corr_fact * ES%psi_lim
      ES%limiter_plasma = .true.
      ES%active_xpoint  = 0
      
    end if
    
    ES%psi_bnd = ES%psi_bnd * corr_fact ! Undo corr_fact
    
    ! --- psi_axis < psi_bnd or > psi_bnd?
    ES%axis_is_psi_minimum = ( ES%psi_axis < ES%psi_bnd )
    
    ! --- Strike points.
    ES%num_strike          = 0
    ES%R_strike(:)         = -99.d0
    ES%Z_strike(:)         = -99.d0
    ES%i_bndelm_strike(:)  = -9999
    ES%s_strike(:)         = -99.d0
    call find_strike(node_list, bnd_elm_list, ES)
    
    ! --- Inner midplane point.
    R1 = ES%R_axis
    R2 = 0.d0
    i  = 0
    do
      i = i + 1
      call find_RZ(node_list, element_list, (R1+R2)/2.d0, ES%Z_axis, R_out, Z_out, i_out, s_out,   &
        t_out, ifail)
      if ( ifail == 0 ) then
        R1 = (R1 + R2) / 2.d0
      else
        R2 = (R1 + R2) / 2.d0
      end if
      if ( abs(R2-R1) < 1.d-4 ) exit
    end do
    ES%R_midpl(1) = (R2+R1)/2.d0
    
    ! --- Outer midplane point.
    R1 = ES%R_axis
    R2 = 3.d0 * ES%R_axis
    i  = 0
    do
      i = i + 1
      call find_RZ(node_list, element_list, (R1+R2)/2.d0, ES%Z_axis, R_out, Z_out, i_out, s_out,   &
        t_out, ifail)
      if ( ifail == 0 ) then
        R1 = (R1 + R2) / 2.d0
      else
        R2 = (R1 + R2) / 2.d0
      end if
      if ( abs(R2-R1) < 1.d-4 ) exit
    end do
    ES%R_midpl(2) = (R2+R1)/2.d0
    
    ES%initialized = .true.
    
  end subroutine update_equil_state
  
  
  
  !> Readable output of the equilibrium state for the logfile.
  subroutine print_equil_state(ES, verbose)
    
    ! --- Routine parameters.
    type(t_equil_state), intent(in) :: ES
    logical,             intent(in) :: verbose    !< Output much additional information?
    
    ! --- Local variables
    integer :: i
    
    if ( .not. ES%initialized ) then
      write(*,*) 'ERROR in mod_equil_info|print_equil_state: ES used before initialization.'
      write(*,*) 'Call update_equil_state first.'
      stop
    end if
    
    101 format(1x,a)
    102 format(1x,5(a,f10.5))
    103 format(1x,5(a,i10))
    104 format(1x,a,i3,2f10.5)
    105 format(1x,a,i3.3,a,f10.5)
    106 format(1x,a,i3.3,a,i10)
    
    ! --- General description of the plasma.
    write(*,*)
    write(*,*) '=============================================================='
    if ( ES%limiter_plasma ) then
      write(*,101) 'Plasma_type        = Limiter Plasma'
    else
      write(*,101) 'Plasma_type        = X-Point Plasma'
      if ( ES%xcase == LOWER_XPOINT ) then
        write(*,101) 'Xpoint_type        = Lower X-Point'
      else if ( ES%xcase == UPPER_XPOINT ) then
        write(*,101) 'Xpoint_type        = Upper X-Point'
      else
        write(*,101) 'Xpoint_type        = Double X-Point'
        if ( ES%active_xpoint == LOWER_XPOINT ) then
          write(*,101) 'Active_xpoint      = Lower X-Point'
        else
          write(*,101) 'Active_xpoint      = Upper X-Point'
        end if
      end if
    end if
    if (verbose ) then
      if ( ES%axis_is_psi_minimum ) then
        write(*,101) 'Psi_axis           = Minimum'
      else
        write(*,101) 'Psi_axis           = Maximum'
      end if
    end if
    write(*,102) 'Psi_bnd            =', ES%psi_bnd
    
    ! --- Magnetic axis.
    write(*,*) '--- Magnetic Axis --------------------------------------------'
    write(*,102) 'R_axis             =', ES%R_axis
    write(*,102) 'Z_axis             =', ES%Z_axis
    write(*,102) 'Psi_axis           =', ES%Psi_axis
    if ( verbose ) then
      write(*,103) 'i_elm_axis         =', ES%i_elm_axis
      write(*,102) 's_axis             =', ES%s_axis
      write(*,102) 't_axis             =', ES%t_axis
      write(*,103) 'ifail_axis         =', ES%ifail_axis
    end if
    
    ! --- Limiter point.
    if ( ES%limiter_plasma .or. verbose ) then
      write(*,*) '--- Limiter Point --------------------------------------------'
      write(*,102) 'R_lim              =', ES%R_lim
      write(*,102) 'Z_lim              =', ES%Z_lim
      write(*,102) 'Psi_lim            =', ES%Psi_lim
      if ( verbose ) then
        write(*,103) 'i_elm_lim          =', ES%i_elm_lim
        write(*,102) 's_lim              =', ES%s_lim
        write(*,102) 't_lim              =', ES%t_lim
        write(*,103) 'ifail_lim          =', ES%ifail_lim
      end if
    end if
    
    ! --- X-point(s).
    if ( ES%xpoint ) then
      if ( ( ES%xcase == LOWER_XPOINT ) .or. ( ES%xcase == DOUBLE_NULL ) ) then
        write(*,*) '--- Lower X-Point --------------------------------------------'
        write(*,102) 'R_xpoint1          =', ES%R_xpoint(1)
        write(*,102) 'Z_xpoint1          =', ES%Z_xpoint(1)
        write(*,102) 'Psi_xpoint1        =', ES%Psi_xpoint(1)
        if ( verbose ) then
          write(*,103) 'i_elm_xpoint1      =', ES%i_elm_xpoint(1)
          write(*,102) 's_xpoint1          =', ES%s_xpoint(1)
          write(*,102) 't_xpoint1          =', ES%t_xpoint(1)
          write(*,103) 'ifail_xpoint       =', ES%ifail_xpoint
        end if
      end if
      if ( ( ES%xcase == UPPER_XPOINT ) .or. ( ES%xcase == DOUBLE_NULL ) ) then
        write(*,*) '--- Upper X-Point --------------------------------------------'
        write(*,102) 'R_xpoint2          =', ES%R_xpoint(2)
        write(*,102) 'Z_xpoint2          =', ES%Z_xpoint(2)
        write(*,102) 'Psi_xpoint2        =', ES%Psi_xpoint(2)
        if ( verbose ) then
          write(*,103) 'i_elm_xpoint2      =', ES%i_elm_xpoint(2)
          write(*,102) 's_xpoint2          =', ES%s_xpoint(2)
          write(*,102) 't_xpoint2          =', ES%t_xpoint(2)
          write(*,103) 'ifail_xpoint       =', ES%ifail_xpoint
        end if
      end if
    end if
    
    ! --- Strike points.
    if ( verbose ) then
      do i = 1, ES%num_strike
        write(*,'(1x,"--- Strike Point",i3,1x,42("-"))') i
        write(*,105) 'R_strike', i, '        =', ES%R_strike(i)
        write(*,105) 'Z_strike', i, '        =', ES%Z_strike(i)
        write(*,106) 'i_bndelm_strike', i, ' =', ES%i_bndelm_strike(i)
        write(*,105) 's_strike', i, '        =', ES%s_strike(i)
      end do
    end if
    
    ! --- Midplane points.
    if ( verbose ) then
      write(*,*) '--- Midplane Points ------------------------------------------'
      write(*,102) 'R_midpl1           =', ES%R_midpl(1)
      write(*,102) 'R_midpl2           =', ES%R_midpl(2)
    end if
    
    write(*,*) '=============================================================='
    write(*,*)
    
  end subroutine print_equil_state
  
  
  
  !> Nice, readable output of the equilibrium state.
  subroutine save_special_points(ES, filename, append, ioerr)
    
    ! --- Routine parameters.
    type(t_equil_state), intent(in)    :: ES
    character(len=*),    intent(in)    :: filename   !< Output to which file?
    logical,             intent(in)    :: append     !< Append to existing file?
    integer,             intent(inout) :: ioerr      !< I/O Error code.
    
    ! --- Local variables
    integer, parameter :: ifile = 33 !### better solution?
    integer :: i
    
    if ( append ) then
      open(ifile, file=filename, form='formatted', access='append', iostat=ioerr)
      if ( ioerr /= 0 ) then
        write(ifile,*)
        write(ifile,*)
      end if
    else
      open(ifile, file=filename, form='formatted', status='replace', iostat=ioerr)
    end if
    
    if ( ioerr /= 0 ) then
      write(*,*) 'ERROR in mod_equil_info|save_special_points opening file "',trim(filename),'".'
      return
    end if
    
    write(ifile,*) '# Magnetic Axis'
    write(ifile,*) ES%R_axis, ES%Z_axis
    write(ifile,*)
    write(ifile,*)
    
    write(ifile,*) '# Limiter Point'
    write(ifile,*) ES%R_lim, ES%Z_lim
    write(ifile,*)
    write(ifile,*)
    
    if ( ES%xpoint ) then
      if ( (ES%xcase==LOWER_XPOINT) .or. (ES%xcase==DOUBLE_NULL) ) then
        write(ifile,*) '# Lower X-Point'
        write(ifile,*) ES%R_xpoint(1), ES%Z_xpoint(1)
        write(ifile,*)
        write(ifile,*)
      end if
      if ( (ES%xcase==UPPER_XPOINT) .or. (ES%xcase==DOUBLE_NULL) ) then
        write(ifile,*) '# Upper X-Point'
        write(ifile,*) ES%R_xpoint(2), ES%Z_xpoint(2)
        write(ifile,*)
        write(ifile,*)
      end if
    end if
    
    do i = 1, ES%num_strike
      write(ifile,*) '# Strike point', i
      write(ifile,*) ES%R_strike(i), ES%Z_strike(i)
      write(ifile,*)
      write(ifile,*)
    end do
    
    write(ifile,*) '# Inner midplane point'
    write(ifile,*) ES%R_midpl(1), ES%Z_axis
    write(ifile,*)
    write(ifile,*)
    
    write(ifile,*) '# Outer midplane point'
    write(ifile,*) ES%R_midpl(2), ES%Z_axis
    write(ifile,*)
    write(ifile,*)
    
  end subroutine save_special_points
  
  
  
  !> Calculate Psi_N for given Psi.
  real*8 function get_psi_n(ES, psi)
    
    ! --- Routine parameters.
    type(t_equil_state), intent(in) :: ES
    real*8,              intent(in) :: psi    !< Poloidal flux value.
    
    get_psi_n = ( psi - ES%psi_axis ) / ( ES%psi_bnd - ES%psi_axis )
    
  end function get_psi_n
  
  
  
end module equil_info 
