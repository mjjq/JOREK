!> This module contains datastructures describing poloidal and toroidal positions.
!! 
!! Together with mod_expression, this provides a general diagnostic framework for many applications.
module mod_position
  
  
  
  
  
  use constants
  use parameters
  use equil_info
  use data_structure
  
  
  
  
  
  implicit none
  
  
  
  
  
  public
  
  
  
  
  
  ! --- Constants
  character(len=12), parameter, private :: THIS_MOD_NAME = 'mod_position'
  
  
  
  
  
  ! --- Data structures
  
  !> Data structure for a poloidal position
  type t_pol_pos
    real*8             :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, s, t
    integer            :: ielm
    type(type_element) :: element
    type(type_node)    :: nodes(n_vertex_max)
  end type t_pol_pos
  
  !> Data structure for a list of poloidal positions
  type t_pol_pos_list
    integer :: n_pos = 0
    type(t_pol_pos), allocatable :: pos(:)
  end type t_pol_pos_list
  
  !> Data structure for a toroidal position
  type t_tor_pos
    real*8  :: phi
  end type t_tor_pos
  
  !> Data structure for a list of toroidal positions
  type t_tor_pos_list
    integer :: n_pos = 0
    type(t_tor_pos), allocatable :: pos(:)
    logical :: full_period
  end type t_tor_pos_list
  
  
  
  
  
  contains
  
  
  
  
  
  !> Clean up a poloidal position data structure.
  subroutine cleanup_pol_pos(pos_list)
    
    type(t_pol_pos_list), intent(inout) :: pos_list
    
    pos_list%n_pos = 0
    
    if ( allocated(pos_list%pos) ) deallocate( pos_list%pos )
    
  end subroutine cleanup_pol_pos
  
  
  
  
  
  !> Clean up a toroidal position data structure.
  subroutine cleanup_tor_pos(pos_list)
    
    type(t_tor_pos_list), intent(inout) :: pos_list
    
    pos_list%n_pos       = 0
    pos_list%full_period = .false.
    
    if ( allocated(pos_list%pos) ) deallocate( pos_list%pos )
    
  end subroutine cleanup_tor_pos
  
  
  
  
  
  !> Allocate a poloidal position data structure.
  subroutine alloc_pol_pos(pos_list, n_pos)
    
    type(t_pol_pos_list), intent(inout) :: pos_list
    integer,              intent(in)    :: n_pos
    
    call cleanup_pol_pos(pos_list)
    
    pos_list%n_pos = n_pos
    
    allocate( pos_list%pos(n_pos) )
    
  end subroutine alloc_pol_pos
  
  
  
  
  
  !> Allocate a toroidal position data structure.
  subroutine alloc_tor_pos(pos_list, n_pos)
    
    type(t_tor_pos_list), intent(inout) :: pos_list
    integer,              intent(in)    :: n_pos
    
    call cleanup_tor_pos(pos_list)
    
    pos_list%n_pos = n_pos
    
    allocate( pos_list%pos(n_pos) )
    
  end subroutine alloc_tor_pos
  
  
  
  
  
  !> Create a poloidal position data structure (one or several poloidal positions) from:
  !! - (R,Z) or (ielm,s,t)           -> single position
  !! - (Rmin,Rmax,nR,Zmin,Zmax,nZ)   -> 2D array of positions
  !! - (Rstart,Rend,Zstart,Zend,n)   -> Equidistant points along a straight line
  !! - (Psi_N)                       -> flux surface
  !! To be added: All nodes; All nodes with subdivision of elements (for vtk)
  subroutine create_pol_pos(pos_list, ierr, node_list, element_list, eq, R, Z, ielm, s, t, Rmin,  &
    Rmax, nR, Zmin, Zmax, nZ, Rstart, Rend, Zstart, Zend, n, PsiN)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':create_pol_pos'
    
    ! --- Routine parameters
    type(t_pol_pos_list), target, intent(inout) :: pos_list
    integer,                      intent(out)   :: ierr
    type(type_node_list),         intent(in)    :: node_list
    type(type_element_list),      intent(in)    :: element_list
    type(t_equil_state),          intent(in)    :: eq
    real*8,  optional,            intent(in)    :: R, Z, s, t, Rmin, Rmax, Zmin, Zmax, PsiN,       &
      Rstart, Rend, Zstart, Zend
    integer, optional,            intent(in)    :: ielm, nR, nZ, n
    
    ! --- Local variables
    type(t_pol_pos), pointer :: pos
    real*8  :: R_out,Z_out
    integer :: i
    
    ierr = 0
    
    if ( present(R) .and. present(Z) ) then
      
      call alloc_pol_pos(pos_list, 1)
      pos   => pos_list%pos(1)
      pos%R = R
      pos%Z = Z
      call find_RZ(node_list, element_list, R, Z, R_out, Z_out, pos%ielm, pos%s, pos%t, ierr)
      if ( ierr /= 0 ) write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//' when calling find_RZ.'
      
      if ( ierr == 0 ) call fill_pol_pos(pos, node_list, element_list, ierr)
      
    else if ( present(ielm) .and. present(s) .and. present(t) ) then
      
      call alloc_pol_pos(pos_list, 1)
      pos      => pos_list%pos(1)
      pos%ielm = ielm
      pos%s    = s
      pos%t    = t
      
      call fill_pol_pos(pos, node_list, element_list, ierr)
      
    else if ( present(Rmin) .and. present(Rmax) .and. present(nR) .and. present(Zmin) .and.        &
      present(Zmax) .and. present(nZ) ) then
      
      !###
      write(*,*) '### not implemented yet ###'
      stop
      !###
      
    else if ( present(Rstart) .and. present(Rend) .and. present(Zstart) .and. present(Zend) .and.  &
      present(n) ) then
      
      call alloc_pol_pos(pos_list, n)
      do i = 1, n
        pos   => pos_list%pos(i)
        pos%R = Rstart + (Rend-Rstart) * real(i-1)/real(n-1)
        pos%Z = Zstart + (Zend-Zstart) * real(i-1)/real(n-1)
        call find_RZ(node_list, element_list, pos%R, pos%Z, R_out, Z_out, pos%ielm, pos%s, pos%t,  &
          ierr)
        if ( ierr /= 0 ) then
          write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//' when calling find_RZ.'
          exit
        end if
        call fill_pol_pos(pos, node_list, element_list, ierr)
      end do
      
    else
      
      ierr = 99
      write(*,*)
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//':'
      write(*,*) 'No valid representation for poloidal position(s) provided.'
      write(*,*)
      
    end if
    
  end subroutine create_pol_pos
  
  
  
  
  
  !> Fill information (R, R_s, ..., Z_tt, element, nodes) for one poloidal position.
  !!
  !! Requires that ielm, s, t are already set to correct values.
  subroutine fill_pol_pos(pos, node_list, element_list, ierr)
    
    ! --- Routine parameters
    type(t_pol_pos), pointer, intent(inout) :: pos
    type(type_node_list),     intent(in)    :: node_list
    type(type_element_list),  intent(in)    :: element_list
    integer,                  intent(out)   :: ierr
    
    ! --- Local variables
    integer :: i
    
    ierr = 0
    
    call interp_RZ(node_list, element_list, pos%ielm, pos%s, pos%t, pos%R, pos%R_s, pos%R_t,       &
      pos%R_st, pos%R_ss, pos%R_tt, pos%Z, pos%Z_s, pos%Z_t, pos%Z_st, pos%Z_ss, pos%Z_tt)
    
    !### check that ielm in valid range
    
    pos%element = element_list%element(pos%ielm)
    
    do i = 1, n_vertex_max
      pos%nodes(i) = node_list%node(pos%element%vertex(i))
    end do
    
  end subroutine fill_pol_pos
  
  
  
  
  
  !> Create a toroidal position data structure from the specified information:
  !! - phi or iplane            -> single toroidal position
  !! - (phistart, phiend, nphi) -> several equidistant toroidal positions
  !! - nphi                     -> distribute nphi positions over full toroidal turn (or period)
  subroutine create_tor_pos(pos_list, ierr, phi, iplane, phistart, phiend, nphi)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':create_tor_pos'
    
    ! --- Routine parameters
    type(t_tor_pos_list), target, intent(inout) :: pos_list
    real*8,  optional,            intent(in)    :: phi, phistart, phiend
    integer, optional,            intent(in)    :: iplane, nphi
    integer,                      intent(out)   :: ierr
    
    ! --- Local variables
    type(t_tor_pos), pointer :: pos
    integer :: i
    
    ierr = 0
    
    if ( present(phi) ) then
      
      call alloc_tor_pos(pos_list, 1)
      pos     => pos_list%pos(1)
      pos%phi = phi
      
    else if ( present(iplane) ) then
      
      call alloc_tor_pos(pos_list, 1)
      pos      => pos_list%pos(1)
      pos%phi = 2.d0 * PI * real(iplane) / real(n_plane * n_period) !###check###
      
    else if ( present(phistart) .and. present(phiend) .and. present(nphi) ) then
      
      call alloc_tor_pos(pos_list, nphi)
      do i = 1, nphi
        pos     => pos_list%pos(i)
        pos%phi = phistart + (phiend-phistart) * real(i-1)/real(nphi-1)
      end do
      
    else if ( present(nphi) ) then
      
      call alloc_tor_pos(pos_list, nphi)
      do i = 1, nphi
        pos     => pos_list%pos(i)
        pos%phi = (2.d0*PI)/real(n_period) * real(i-1)/real(nphi)
        ! IMPORTANT: Position phi=0 included and position phi=2pi/n_period ommitted!
      end do
      pos_list%full_period = .true.
      
    else
      
      ierr = 99
      write(*,*)
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//':'
      write(*,*) 'No valid representation for toroidal position(s) provided.'
      write(*,*)
      
    end if
    
  end subroutine create_tor_pos
  
  
  
  
  
end module mod_position
