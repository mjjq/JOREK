!> This module contains datastructures describing lists of poloidal and toroidal positions.
!! 
!! Together with mod_expression, this provides a general diagnostic framework for many applications.
module mod_position
  
  
  
  
  
  use constants
  use parameters
  use equil_info
  use data_structure
  use mod_straight_field_line
  
  
  
  
  
  implicit none
  
  
  
  
  
  public
  
  
  
  
  
  ! --- Constants
  character(len=12), parameter, private :: THIS_MOD_NAME = 'mod_position'
  
  
  
  
  
  ! --- Data structures
  
  !> Data structure for a poloidal position
  type t_pol_pos
    logical            :: outside = .false. !< Outside computational domain?
    real*8             :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, s, t
    integer            :: ielm
    type(type_element) :: element
    type(type_node)    :: nodes(n_vertex_max)
    ! --- The following quantities will only be available in certain cases of poloidal positions
    real*8             :: theta_star !< Straight field line angle (for flux surfaces)
    real*8             :: r_minor    !< Minor radius from A = r_minor^2 pi (for flux surfaces)
    real*8             :: length     !< Length along a line (for line)
  end type t_pol_pos
  
  !> Data structure for a list of poloidal positions
  type t_pol_pos_list
    type(t_pol_pos), allocatable :: pos(:,:)
    integer :: n_pos(2) = 0
    logical :: full_turn !< Do the positions cover a full poloidal turn?
  end type t_pol_pos_list
  
  !> Data structure for a toroidal position
  type t_tor_pos
    real*8  :: phi
  end type t_tor_pos
  
  !> Data structure for a list of toroidal positions
  type t_tor_pos_list
    type(t_tor_pos), allocatable :: pos(:)
    integer :: n_pos = 0
    logical :: full_period !< Do the positions cover a full toroidal period?
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
    integer,              intent(in)    :: n_pos(2)
    
    integer :: i, j
    
    call cleanup_pol_pos(pos_list)
    
    pos_list%n_pos(:) = n_pos(:)
    
    allocate( pos_list%pos(n_pos(1),n_pos(2)) )
    
    do j = 1, n_pos(2)
      do i = 1, n_pos(1)
        pos_list%pos(i,j)%theta_star = 0.d0/0.d0
        pos_list%pos(i,j)%length     = 0.d0/0.d0
        pos_list%pos(i,j)%r_minor    = 0.d0/0.d0
      end do
    end do
    
  end subroutine alloc_pol_pos
  
  
  
  
  
  !> Allocate a toroidal position data structure.
  subroutine alloc_tor_pos(pos_list, n_pos)
    
    type(t_tor_pos_list), intent(inout) :: pos_list
    integer,              intent(in)    :: n_pos
    
    call cleanup_tor_pos(pos_list)
    
    pos_list%n_pos = n_pos
    
    allocate( pos_list%pos(n_pos) )
    
  end subroutine alloc_tor_pos
  
  
  
  
  
  !> Simple function to generate poloidal positions (wrapper for routine create_pol_pos).
  function pol_pos(node_list, element_list, eq, R, Z, ielm, s, t, Rmin, Rmax, nR, Zmin, Zmax, nZ,  &
    Rstart, Rend, Zstart, Zend, n, PsiN, nTht, PsiNmin, PsiNmax, nPsiN) result(pos_list)
    type(t_pol_pos_list), target :: pos_list
    
    ! --- Routine parameters
    type(type_node_list),         intent(in)    :: node_list
    type(type_element_list),      intent(in)    :: element_list
    type(t_equil_state),          intent(in)    :: eq
    real*8,  optional,            intent(in)    :: R, Z, s, t, Rmin, Rmax, Zmin, Zmax, PsiN,       &
      Rstart, Rend, Zstart, Zend, PsiNmin, PsiNmax
    integer, optional,            intent(in)    :: ielm, nR, nZ, n, nTht, nPsiN
    
    ! --- Local variables
    integer :: ierr
    
    call create_pol_pos(pos_list, ierr, node_list, element_list, eq, R, Z, ielm, s, t, Rmin, Rmax, &
      nR, Zmin, Zmax, nZ, Rstart, Rend, Zstart, Zend, n, PsiN, nTht, PsiNmin, PsiNmax, nPsiN)
    
  end function pol_pos
  
  
  
  
  
  !> Create a poloidal position data structure (one or several poloidal positions) from:
  !! - (R,Z) or (ielm,s,t)           -> single position
  !! - (Rmin,Rmax,nR,Zmin,Zmax,nZ)   -> 2D array of positions
  !! - (Rstart,Rend,Zstart,Zend,n)   -> Equidistant points along a straight line
  !! - (PsiN,nTht)                   -> flux surface (equidistant points in theta*)
  !! - (PsiNmin,PsiNmax,nPsiN,nTht)  -> flux surfaces (equidistant points in theta*)
  !! To be added: Single node; All nodes; All nodes with subdivision of elements (for vtk)
  recursive subroutine create_pol_pos(pos_list, ierr, node_list, element_list, eq, R, Z, ielm, s,  &
    t, Rmin, Rmax, nR, Zmin, Zmax, nZ, Rstart, Rend, Zstart, Zend, n, PsiN, nTht, PsiNmin, PsiNmax,&
    nPsiN)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':create_pol_pos'
    
    ! --- Routine parameters
    type(t_pol_pos_list), target, intent(inout) :: pos_list
    integer,                      intent(out)   :: ierr
    type(type_node_list),         intent(in)    :: node_list
    type(type_element_list),      intent(in)    :: element_list
    type(t_equil_state),          intent(in)    :: eq
    real*8,  optional,            intent(in)    :: R, Z, s, t, Rmin, Rmax, Zmin, Zmax, PsiN,       &
      Rstart, Rend, Zstart, Zend, PsiNmin, PsiNmax
    integer, optional,            intent(in)    :: ielm, nR, nZ, n, nTht, nPsiN
    
    ! --- Local variables
    type(t_theta_mapping) :: mapping
    type(t_pol_pos), pointer :: pos
    real*8  :: R_out, Z_out, hh, gx, gy, gg, ax, ay, full_length
    integer :: i, j
    real*8, allocatable :: surface(:) !< Poloidal surface inside flux surface (for r_minor)
    
    ierr = 0
    
    ! --- Single position given by R and Z.
    if ( present(R) .and. present(Z) ) then
      
      call alloc_pol_pos(pos_list, (/1,1/))
      pos   => pos_list%pos(1,1)
      pos%R = R
      pos%Z = Z
      call find_RZ(node_list, element_list, R, Z, R_out, Z_out, pos%ielm, pos%s, pos%t, ierr)
      if ( ierr /= 0 ) write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//' when calling find_RZ.'
      
      if ( ierr == 0 ) call fill_pol_pos(pos, node_list, element_list, ierr)
      
    ! --- Single position given by ielm, s, and t.
    else if ( present(ielm) .and. present(s) .and. present(t) ) then
      
      call alloc_pol_pos(pos_list, (/1,1/))
      pos      => pos_list%pos(1,1)
      pos%ielm = ielm
      pos%s    = s
      pos%t    = t
      
      call fill_pol_pos(pos, node_list, element_list, ierr)
      
    ! --- Rectangular array of positions in R and Z.
    else if ( present(Rmin) .and. present(Rmax) .and. present(nR) .and. present(Zmin) .and.        &
      present(Zmax) .and. present(nZ) ) then
      
      call alloc_pol_pos(pos_list, (/nR,nZ/))
      !$OMP parallel do default(none) firstprivate(i,j,pos,R_out,Z_out,ierr)                       &
      !$OMP shared(nR,nZ,node_list,element_list,pos_list,Rmin,Rmax,Zmin,Zmax) schedule(dynamic)
      do i = 1, nR
        do j = 1, nZ
          pos   => pos_list%pos(i,j)
          pos%R = Rmin + (Rmax-Rmin) * real(i-1)/real(nR-1)
          pos%Z = Rmin + (Zmax-Zmin) * real(j-1)/real(nZ-1)
          call find_RZ(node_list, element_list, pos%R, pos%Z, R_out, Z_out, pos%ielm, pos%s, pos%t,&
            ierr)
          if ( ierr /= 0 ) then
            pos%outside = .true.
          else
            call fill_pol_pos(pos, node_list, element_list, ierr)
          end if
        end do
      end do
      !$OMP end parallel do
      
    ! --- Equidistant points along a straight line in R and Z.
    else if ( present(Rstart) .and. present(Rend) .and. present(Zstart) .and. present(Zend) .and.  &
      present(n) ) then
      
      call alloc_pol_pos(pos_list, (/1,n/))
      full_length = sqrt( (Rend-Rstart)**2 + (Zend-Zstart)**2 )
      do i = 1, n
        pos   => pos_list%pos(1,i)
        pos%R = Rstart + (Rend-Rstart) * real(i-1)/real(n-1)
        pos%Z = Zstart + (Zend-Zstart) * real(i-1)/real(n-1)
        pos%length = full_length * real(i-1)/real(n-1)
        call find_RZ(node_list, element_list, pos%R, pos%Z, R_out, Z_out, pos%ielm, pos%s, pos%t,  &
          ierr)
        if ( ierr /= 0 ) then
          write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//' calling find_RZ.'
          exit
        end if
        call fill_pol_pos(pos, node_list, element_list, ierr)
        if ( ierr /= 0 ) exit
      end do
      
    ! --- Equidistant points along a straight line in R at fixed Z (e.g., midplane profiles).
    else if ( present(Rstart) .and. present(Rend) .and. present(Z) .and. present(n) ) then
      
      call create_pol_pos(pos_list, ierr, node_list, element_list, eq, Rstart=Rstart, Rend=Rend,   &
        Zstart=Z, Zend=Z, n=n)
      
    ! --- Single (closed) flux surface.
    else if ( present(PsiN) .and. present(nTht) ) then
      
      call create_pol_pos(pos_list, ierr, node_list, element_list, eq, PsiNmin=PsiN, PsiNmax=PsiN, &
        nPsiN=1)
      
    ! --- Several (closed) flux surfaces.
    else if ( present(PsiNmin) .and. present(PsiNmax) .and. present(nPsiN) .and. present(nTht) )   &
      then
      
      if ( (min(PsiNmin,PsiNmax) < 0.d0) .and. (max(PsiNmin,PsiNmax) > 1.d0) ) then
        write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': PsiNmin and PsiNmax must be between 0 and 1.'
        ierr = 300
        return
      end if 
      
      call determine_theta_mag(mapping, node_list, element_list, eq, (/PsiNmin,PsiNmax/), nPsiN,   &
        nTht, ierr)
      if ( ierr /= 0 ) then
        write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//' calling determine_theta_mag.'
        ierr = 400
        return
      end if 
      
      allocate(surface(nPsiN))
      surface(:) = 0.d0
      call alloc_pol_pos(pos_list, (/nTht,nPsiN/))
      pos_list%full_turn = .true.
      
      do j = 1, nTht
        do i = 1, nPsiN
          pos   => pos_list%pos(j,i)
          pos%R = mapping%rre(i,j-1)
          pos%Z = mapping%zze(i,j-1)
          pos%theta_star = 2.d0 * PI * real(j-1) / real(nTht-1) !######### check if nTht or nTht-1
          call find_RZ(node_list, element_list, pos%R, pos%Z, R_out, Z_out, pos%ielm, pos%s, pos%t,&
            ierr)
          if ( ierr /= 0 ) then
            write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//' calling find_RZ.'
            exit
          end if
          call fill_pol_pos(pos, node_list, element_list, ierr)
          if ( ierr /= 0 ) exit
          
          ! --- Calculate poloidal surface inside flux surface (for r_minor)
          if ( j /= 0 ) then
            gx = pos_list%pos(j,i)%R - pos_list%pos(j-1,i)%R
            gy = pos_list%pos(j,i)%Z - pos_list%pos(j-1,i)%Z
            gg = sqrt( gx**2 + gy**2 )
            ax = pos_list%pos(j,i)%R - eq%R_axis
            ay = pos_list%pos(j,i)%Z - eq%Z_axis
            hh = ( gy * ax - gx * ay ) / gg
            surface(i) = surface(i) + hh * gg / 2.d0
          end if
        end do
      end do
      
      ! --- Fill in r_minor
      do j = 1, nTht
        do i = 1, nPsiN
          pos_list%pos(j,i)%r_minor = sqrt( surface(i) / PI )
        end do
      end do
      
      call cleanup_mapping(mapping)
      deallocate(surface)
      
    ! --- Several (closed) flux surfaces.
    else if ( present(nPsiN) .and. present(nTht) ) then
      
      call create_pol_pos(pos_list, ierr, node_list, element_list, eq, nTht=nTht, nPsiN=nPsiN,     &
        psiNmin=0.005d0, psiNmax=0.990d0)
      
    else
      
      ierr = 99
      write(*,*)
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//':'
      write(*,*) 'No valid representation for poloidal position(s) provided.'
      write(*,*)
      
    end if
    
  end subroutine create_pol_pos
  
  
  
  
  
  !> Auxilliary routine used by create_pol_pos: Fill information (R, R_s, ..., Z_tt, element, nodes)
  !! for a single poloidal position. Requires that ielm, s, t are already set to correct values.
  subroutine fill_pol_pos(pos, node_list, element_list, ierr)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':fill_pol_pos'
    
    ! --- Routine parameters
    type(t_pol_pos), pointer, intent(inout) :: pos
    type(type_node_list),     intent(in)    :: node_list
    type(type_element_list),  intent(in)    :: element_list
    integer,                  intent(out)   :: ierr
    
    ! --- Local variables
    integer :: i
    
    ierr = 0
    
    if ( (pos%ielm < 0) .or. (pos%ielm > element_list%n_elements) ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': pos%ielm has illegal value.'
      ierr = 300
      return
    end if
    
    call interp_RZ(node_list, element_list, pos%ielm, pos%s, pos%t, pos%R, pos%R_s, pos%R_t,       &
      pos%R_st, pos%R_ss, pos%R_tt, pos%Z, pos%Z_s, pos%Z_t, pos%Z_st, pos%Z_ss, pos%Z_tt)
    
    pos%element = element_list%element(pos%ielm)
    
    do i = 1, n_vertex_max
      pos%nodes(i) = node_list%node(pos%element%vertex(i))
    end do
    
  end subroutine fill_pol_pos
  
  
  
  
  
  !> Simple function to generate toroidal positions (wrapper for routine create_tor_pos).
  function tor_pos(phi, iplane, phistart, phiend, nphi) result(pos_list)
    type(t_tor_pos_list), target :: pos_list
    
    ! --- Routine parameters
    real*8,  optional,            intent(in)    :: phi, phistart, phiend
    integer, optional,            intent(in)    :: iplane, nphi
    
    ! --- Local variables
    integer :: ierr
    
    call create_tor_pos(pos_list, ierr, phi, iplane, phistart, phiend, nphi)
    
  end function tor_pos
  
  
  
  
  
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
