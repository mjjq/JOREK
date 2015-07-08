!> This module makes the whole new_diag package available and adds some easy to use functionality.
module mod_new_diag
  
  
  
  
  
  use parameters
  use mod_position
  use mod_straight_field_line
  use mod_expression
  use mod_four_filter
  use mod_diag_output
  
  
  
  
  
  implicit none
  
  
  
  
  
  public
  
  
  
  
  
  ! --- Constants
  character(len=15), parameter, private :: THIS_MOD_NAME      = 'mod_new_diag'
  
  !    --- Used by routine midplane_profile
  integer,           parameter          :: HIGHFIELD_SIDE     = 0
  integer,           parameter          :: LOWFIELD_SIDE      = 1
  integer,           parameter          :: BOTH_SIDES         = 2
  
  !   --- Used by routine fourier_analysis
  integer,           parameter          :: OUTP_ABS_VALUE     = 0
  integer,           parameter          :: OUTP_REALPART      = 1
  integer,           parameter          :: OUTP_IMAGINARYPART = 2
  integer,           parameter          :: OUTP_PHASE         = 3
  character(len=33), parameter, private :: OUTP_NAMES(0:3) = (/ 'absolute_value', 'real_part     ',&
    'imaginary_part', 'complex_phase ' /)
  
  
  
  
  
  
  contains
  
  
  
  
  
  !> Initialize the new_diag framework
  subroutine init_new_diag(verbose)
    
    ! --- Routine paramters
    logical, intent(in) :: verbose !< Print some information
    
    call init_expr()
    if ( verbose ) call print_exprs(exprs_all)
    
  end subroutine init_new_diag
  
  
  
  
  
  !> Toroidally averaged expressions on the midplane.
  subroutine midplane_profile(node_list, element_list, eq, units, expr_list, res1d, side, n_pts,  &
    ierr, filename, append, comment)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':midplane_profile'
    
    ! --- Routine parameters
    type(type_node_list),           intent(in)    :: node_list    !< List of grid nodes
    type(type_element_list),        intent(in)    :: element_list !< List of grid elements
    type(t_equil_state),            intent(in)    :: eq           !< Plasma equilibrium information
    integer,                        intent(in)    :: units        !< Output in which units?
    type(t_expr_list),              intent(in)    :: expr_list    !< List of expressions to evaluate
    real*8, allocatable,            intent(inout) :: res1d(:,:)   !< Result array
    integer,                        intent(in)    :: side         !< Side of plasma (hfs, lfs, both)
    integer,                        intent(in)    :: n_pts        !< Number of points in profiles
    integer,                        intent(out)   :: ierr         !< Error code
    character(len=*), optional,     intent(in)    :: filename     !< Filename for ascii [optional]
    logical,          optional,     intent(in)    :: append       !< Append or overwrite [optional]
    character(len=*), optional,     intent(in)    :: comment      !< Comment for ascii file [opti.]
    
    ! --- Local variables
    real*8               :: Rstart, Rend
    real*8, allocatable  :: result(:,:,:,:)
    type(t_pol_pos_list) :: pol_pos_list
    type(t_tor_pos_list) :: tor_pos_list
    
    ierr = 0
    
    if ( side == LOWFIELD_SIDE ) then
      Rstart = eq%R_axis     + 1.d-3
      Rend   = eq%R_midpl(2) - 1.d-3
    else if ( side == HIGHFIELD_SIDE ) then
      Rstart = eq%R_midpl(1) + 1.d-3
      Rend   = eq%R_axis     - 1.d-3
    else if ( side == BOTH_SIDES ) then
      Rstart = eq%R_midpl(1) + 1.d-3
      Rend   = eq%R_midpl(2) - 1.d-3
    else
      ierr = 100
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': parameter side has illegal value'
      return
    end if
    pol_pos_list = pol_pos(node_list, element_list, eq, Rstart=Rstart, Rend=Rend, Z=eq%Z_axis,     &
      n=n_pts)
    tor_pos_list = tor_pos(nphi=4*n_plane) !###
    
    call eval_expr(eq, units, expr_list, pol_pos_list, tor_pos_list, result, ierr)
    call apply_four_filter(result, simple_filter(n=0), expr_list%n_coord, ierr)
    call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)
    
    if ( allocated(result) ) deallocate(result)
    call cleanup_pol_pos(pol_pos_list)
    call cleanup_tor_pos(tor_pos_list)
    
    if ( present(filename) ) then
      call write_ascii_1d(ierr, eq, expr_list, res1d, FORM_TABLE, header=.true.,                   &
        filename=filename, append=append, blanks=.true., comment=comment)
    end if
    
  end subroutine midplane_profile
  
  
  
  
  
  !> Construct poloidally and toroidally averaged profiles.
  subroutine average_profiles(node_list, element_list, eq, units, expr_list, res1d, nPsiN, ierr,   &
    filename, append, comment)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':midplane_profile'
    
    ! --- Routine parameters
    type(type_node_list),           intent(in)    :: node_list    !< List of grid nodes
    type(type_element_list),        intent(in)    :: element_list !< List of grid elements
    type(t_equil_state),            intent(in)    :: eq           !< Plasma equilibrium information
    integer,                        intent(in)    :: units        !< Output in which units?
    type(t_expr_list),              intent(in)    :: expr_list    !< List of expressions to evaluate
    real*8, allocatable,            intent(inout) :: res1d(:,:)   !< Result array
    integer,                        intent(in)    :: nPsiN        !< Number of points for profiles
    integer,                        intent(out)   :: ierr         !< Error code
    character(len=*), optional,     intent(in)    :: filename     !< Filename for ascii [optional]
    logical,          optional,     intent(in)    :: append       !< Append or overwrite [optional]
    character(len=*), optional,     intent(in)    :: comment      !< Comment for ascii file [opti.]
    
    ! --- Local variables
    real*8, allocatable  :: result(:,:,:,:)
    type(t_pol_pos_list) :: pol_pos_list
    type(t_tor_pos_list) :: tor_pos_list
    
    ierr = 0
    
    pol_pos_list = pol_pos(node_list, element_list, eq, nPsiN=nPsiN, nTht=6)!###*4*n_plane) !###
    tor_pos_list = tor_pos(nphi=4)!###*n_plane) !###
    
    call eval_expr(eq, units, expr_list, pol_pos_list, tor_pos_list, result, ierr)
    call apply_four_filter(result, simple_filter(m=0,n=0), expr_list%n_coord, ierr)
    call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)
    
    if ( allocated(result) ) deallocate(result)
    call cleanup_pol_pos(pol_pos_list)
    call cleanup_tor_pos(tor_pos_list)
    
    if ( present(filename) ) then
      call write_ascii_1d(ierr, eq, expr_list, res1d, FORM_TABLE, header=.true.,                   &
        filename=filename, append=append, blanks=.true., comment=comment)
    end if
    
  end subroutine average_profiles
  
  
  
  
  
  !> Profiles along a straight line in the poloidal plane.
  subroutine pol_lineout(node_list, element_list, eq, units, expr_list, res1d, phi, Rstart,        &
    Zstart, Rend, Zend, nPts, ierr, filename, append, comment)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':lineout_profiles'
    
    ! --- Routine parameters
    type(type_node_list),           intent(in)    :: node_list    !< List of grid nodes
    type(type_element_list),        intent(in)    :: element_list !< List of grid elements
    type(t_equil_state),            intent(in)    :: eq           !< Plasma equilibrium information
    integer,                        intent(in)    :: units        !< Output in which units?
    type(t_expr_list),              intent(in)    :: expr_list    !< List of expressions to evaluate
    real*8, allocatable,            intent(inout) :: res1d(:,:)   !< Result array
    real*8,                         intent(in)    :: phi          !< Toroidal position
    real*8,                         intent(in)    :: Rstart       !< R-coordinate for start of line
    real*8,                         intent(in)    :: Zstart       !< Z-coordinate for start of line
    real*8,                         intent(in)    :: Rend         !< R-coordinate for end of line
    real*8,                         intent(in)    :: Zend         !< Z-coordinate for end of line
    integer,                        intent(in)    :: nPts         !< Number of points along line
    integer,                        intent(out)   :: ierr         !< Error code
    character(len=*), optional,     intent(in)    :: filename     !< Filename for ascii [optional]
    logical,          optional,     intent(in)    :: append       !< Append or overwrite [optional]
    character(len=*), optional,     intent(in)    :: comment      !< Comment for ascii file [opti.]
    
    ! --- Local variables
    real*8, allocatable  :: result(:,:,:,:)
    type(t_pol_pos_list) :: pol_pos_list
    type(t_tor_pos_list) :: tor_pos_list
    
    ierr = 0
    
    pol_pos_list = pol_pos(node_list, element_list, eq, Rstart=Rstart, Rend=Rend, Zstart=Zstart,   &
      Zend=Zend, n=nPts)
    tor_pos_list = tor_pos(phi=phi)
    
    call eval_expr(eq, units, expr_list, pol_pos_list, tor_pos_list, result, ierr)
    call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)
    
    if ( allocated(result) ) deallocate(result)
    call cleanup_pol_pos(pol_pos_list)
    call cleanup_tor_pos(tor_pos_list)
    
    if ( present(filename) ) then
      call write_ascii_1d(ierr, eq, expr_list, res1d, FORM_TABLE, header=.true.,                   &
        filename=filename, append=append, blanks=.true., comment=comment)
    end if
    
  end subroutine pol_lineout

  !> Integrate expressions along a straight line in the poloidal plane.
  subroutine int_along_pol_lineout(node_list, element_list, eq, units, expr_list, sum, phi, Rstart, &
    Zstart, Rend, Zend, nPts, ierr, filename, append, comment)

    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':lineout_profiles'

    ! --- Routine parameters
    type(type_node_list),           intent(in)    :: node_list    !< List of grid nodes
    type(type_element_list),        intent(in)    :: element_list !< List of grid elements
    type(t_equil_state),            intent(in)    :: eq           !< Plasma equilibrium information
    integer,                        intent(in)    :: units        !< Output in which units?
    type(t_expr_list),              intent(in)    :: expr_list    !< List of expressions to evaluate
    real*8, allocatable,            intent(inout) :: sum(:)       !< Result vector
    real*8,                         intent(in)    :: phi          !< Toroidal position
    real*8,                         intent(in)    :: Rstart       !< R-coordinate for start of line
    real*8,                         intent(in)    :: Zstart       !< Z-coordinate for start of line
    real*8,                         intent(in)    :: Rend         !< R-coordinate for end of line
    real*8,                         intent(in)    :: Zend         !< Z-coordinate for end of line
    integer,                        intent(in)    :: nPts         !< Number of points along line
    integer,                        intent(out)   :: ierr         !< Error code
    character(len=*), optional,     intent(in)    :: filename     !< Filename for ascii [optional]
    logical,          optional,     intent(in)    :: append       !< Append or overwrite [optional]
    character(len=*), optional,     intent(in)    :: comment      !< Comment for ascii file [opti.]

    ! --- Local variables
    real*8, allocatable             :: result(:,:,:,:), res1d(:,:)
    real*8                          :: dl, line_length
    integer                         :: i, n, iexpr
    type(t_pol_pos), pointer        :: pos_i, pos_ip1
    type(t_pol_pos_list), target    :: pol_pos_list
    type(t_tor_pos_list)            :: tor_pos_list

    ierr = 0

    pol_pos_list = pol_pos(node_list, element_list, eq, Rstart=Rstart, Rend=Rend, Zstart=Zstart,   &
      Zend=Zend, n=nPts)
    tor_pos_list = tor_pos(phi=phi)

    call eval_expr(eq, units, expr_list, pol_pos_list, tor_pos_list, result, ierr)
    call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)

    ! --- Integrate along the line
    if ( .not. allocated(sum) ) allocate (sum(expr_list%n_expr))

    do iexpr = 1, expr_list%n_expr

    sum(iexpr) = 0
    line_length = 0

      do i = 1, nPts-1

        pos_i   => pol_pos_list%pos(1,i)
        pos_ip1 => pol_pos_list%pos(1,i+1)
        dl =  ( (pos_ip1%R-pos_i%R)**2 + (pos_ip1%Z-pos_i%Z)**2 )**0.5
        line_length = line_length + dl
        sum(iexpr) = sum(iexpr) + dl*res1d(i,iexpr)

      end do

    ! We do not want to integrate the time along the line
    if ( (expr_list%expr(iexpr)%name).eq.'t' ) sum(iexpr)=res1d(i,iexpr)

    end do

    if ( allocated(result) ) deallocate(result)
    call cleanup_pol_pos(pol_pos_list)
    call cleanup_tor_pos(tor_pos_list)

    if ( present(filename) ) then
      call write_ascii_0d(ierr, eq, expr_list, sum, FORM_TABLE, header=.true., &
        filename=filename, append=append, blanks=.true., comment=comment)
    end if

  end subroutine int_along_pol_lineout
   
  
  
  
  !> Profiles along a toroidal line.
  subroutine tor_lineout(node_list, element_list, eq, units, expr_list, res1d, phi_start, phi_end, &
    R, Z, nPts, ierr, filename, append, comment)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':lineout_profiles'
    
    ! --- Routine parameters
    type(type_node_list),           intent(in)    :: node_list    !< List of grid nodes
    type(type_element_list),        intent(in)    :: element_list !< List of grid elements
    type(t_equil_state),            intent(in)    :: eq           !< Plasma equilibrium information
    integer,                        intent(in)    :: units        !< Output in which units?
    type(t_expr_list),              intent(in)    :: expr_list    !< List of expressions to evaluate
    real*8, allocatable,            intent(inout) :: res1d(:,:)   !< Result array
    real*8,                         intent(in)    :: phi_start    !< Toroidal start position
    real*8,                         intent(in)    :: phi_end      !< Toroidal end position
    real*8,                         intent(in)    :: R            !< R-coordinate of line
    real*8,                         intent(in)    :: Z            !< Z-coordinate of line
    integer,                        intent(in)    :: nPts         !< Number of points along line
    integer,                        intent(out)   :: ierr         !< Error code
    character(len=*), optional,     intent(in)    :: filename     !< Filename for ascii [optional]
    logical,          optional,     intent(in)    :: append       !< Append or overwrite [optional]
    character(len=*), optional,     intent(in)    :: comment      !< Comment for ascii file [opti.]
    
    ! --- Local variables
    real*8, allocatable  :: result(:,:,:,:)
    type(t_pol_pos_list) :: pol_pos_list
    type(t_tor_pos_list) :: tor_pos_list
    
    ierr = 0
    
    pol_pos_list = pol_pos(node_list, element_list, eq, R=R, Z=Z)
    tor_pos_list = tor_pos(phistart=phi_start, phiend=phi_end, nphi=nPts)
    
    call eval_expr(eq, units, expr_list, pol_pos_list, tor_pos_list, result, ierr)
    call reduce_result_to_1d(ierr, result, res1d, i2=1, i3=1)
    
    if ( allocated(result) ) deallocate(result)
    call cleanup_pol_pos(pol_pos_list)
    call cleanup_tor_pos(tor_pos_list)
    
    if ( present(filename) ) then
      call write_ascii_1d(ierr, eq, expr_list, res1d, FORM_TABLE, header=.true.,                   &
        filename=filename, append=append, blanks=.true., comment=comment)
    end if
    
  end subroutine tor_lineout
  
  
  
  
  
  !> Perform a 2D Fourier analysis of the given expressions in straight field line coordinates.
  subroutine fourier_analysis(node_list, element_list, eq, units, expr_list, cp, nPsiN, ierr,      &
    filename_start, output_type2, nsmallsteps)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':fourier_analysis'
    
    type(type_node_list),           intent(in)    :: node_list    !< List of grid nodes
    type(type_element_list),        intent(in)    :: element_list !< List of grid elements
    type(t_equil_state),            intent(in)    :: eq           !< Plasma equilibrium information
    integer,                        intent(in)    :: units        !< Output in which units?
    type(t_expr_list),              intent(in)    :: expr_list    !< List of expressions to evaluate
    complex*16, allocatable,        intent(inout) :: cp(:,:,:,:)  !< Complex Fourier coefficients
    integer,                        intent(in)    :: nPsiN        !< Number of points for profiles
    integer,                        intent(out)   :: ierr         !< Error code
    character(len=*), optional,     intent(in)    :: filename_start !< Start of filename [optional]
    integer,          optional,     intent(in)    :: output_type2 !< Output what? [optional]
    integer,          optional,     intent(in)    :: nsmallsteps  !< Parameter for pol_pos [opt.]
    
    ! --- Local variables
    integer :: nn(4), output_type, m, n, m_max, n_max, i, j
    real*8, allocatable  :: result(:,:,:,:), outp(:,:,:,:), out1d(:,:)
    character(len=256)   :: comment, filename
    type(t_pol_pos_list) :: pol_pos_list
    type(t_tor_pos_list) :: tor_pos_list
    logical :: append
    
    ierr = 0
    
    pol_pos_list = pol_pos(node_list, element_list, eq, nPsiN=nPsiN, nTht=6*4*n_tor,               &
      nsmallsteps=nsmallsteps) !########
    tor_pos_list = tor_pos(nphi=4*n_tor) !########
    
    call eval_expr(eq, units, expr_list, pol_pos_list, tor_pos_list, result, ierr)
    
    call perform_four_trafo(result, cp, POLTOR_TRAFO, ierr)
    
    ! --- Output to files if parameter filename_start is present, otherwise just return cp array.
    if ( present(filename_start) ) then
      
      output_type = OUTP_ABS_VALUE ! preset
      if ( present(output_type2) ) output_type = output_type
      
      nn(:) = (/ size(cp,1), size(cp,2), size(cp,3), size(cp,4) /)
      n_max = nn(1) - 1
      m_max = nn(2) / 2
      
      allocate( outp(nn(1),nn(2),nn(3),nn(4)) )
      allocate( out1d(nn(3),nn(4)) )
      
      ! --- Output absolute value/real part/imaginary part/phase of complex Fourier components?
      if ( output_type == OUTP_ABS_VALUE ) then
        
        outp(:,:,:,:) = abs(cp(:,:,:,:))
        
      else if ( output_type == OUTP_REALPART ) then
        
        outp(:,:,:,:) = real(cp(:,:,:,:))
        
      else if ( output_type == OUTP_IMAGINARYPART ) then
        
        outp(:,:,:,:) = aimag(cp(:,:,:,:))
        
      else if ( output_type == OUTP_PHASE ) then
        
        outp(:,:,:,:) = atan2( real(cp(:,:,:,:)), aimag(cp(:,:,:,:)) ) !#########correct?
        
      else
        
        ierr = 100
        write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': parameter output_type has illegal value'
        return
        
      end if
      
      ! --- Always take the 0/0 component for coordinate expressions!
      do i = 1, expr_list%n_coord
        do j = 1, nn(3)
          outp(:,:,j,i) = real( cp(1,1,j,i) )
        end do
      end do
      
      ! --- Output the Fourier components to ascii files
      do n = 0, (n_tor-1)/2 ! toroidal mode number (needs to be multiplied by n_period!)
        write(filename,'(4a,i3.3,a)') trim(filename_start), '_', trim(OUTP_NAMES(output_type)),    &
          '_n', n*n_period, '.dat'
        do m = -m_max, m_max ! poloidal mode number
          
          if ( m >= 0 ) then
            out1d(:,:) = outp(n+1,m+1,:,:)
          else
            out1d(:,:) = outp(n+1,nn(2)-abs(m)+1,:,:)
          end if
          
          write(comment,'(a,sp,i4.3,a,i4.3)') trim(OUTP_NAMES(output_type))//'s for m/n=', m, '/', &
            n*n_period
          call write_ascii_1d(ierr, eq, expr_list, out1d, FORM_TABLE, .true., filename,            &
            append=(m/=-m_max), comment=trim(comment), blanks=.true.)
          
        end do
      end do
      
    end if
    
    if ( allocated(result) ) deallocate(result)
    if ( allocated(outp )  ) deallocate(outp  )
    if ( allocated(out1d)  ) deallocate(out1d )
    
  end subroutine fourier_analysis
  
  
  
  
  
end module mod_new_diag
