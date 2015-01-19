!> Common parameters and variables for the JOREK free boundary extension.
!!
!! @see vacuum_response, vacuum_equilibrium
module vacuum
  
  implicit none
  
  !> @name General parameters
  logical, parameter  :: vacuum_debug          = .false. !< Enable additional output and tests
  logical, parameter  :: vacuum_decouple_modes = .false. !< Option to switch off 3D wall mode coupling
  integer             :: n_dof_bnd                       !< Total number of boundary dofs per harmonic
  integer             :: n_dof_starwall                  !< Total number of boundary dofs in STARWALL response
  integer, parameter  :: ivar_psi = 1                    !< Index of Psi variable
  integer, parameter  :: ivar_j   = 3                    !< Index of j variable
  
  !> @name Resistive wall only
  real*8              :: wall_resistivity                !< Resistivity of the external wall
  real*8              :: wall_thickness        = 1.d0    !< Thickness of the external wall
  logical             :: wall_curr_initialized = .false. !< Have the wall currents been initialized?
  integer             :: n_wall_curr                     !< Number of wall current potentials.
  real*8, allocatable :: wall_curr(:)                    !< Wall current potentials (\f$Y_k\f$).
  real*8, allocatable :: dwall_curr(:)                   !< Change of wall current potentials (\f$\delta Y_k\f$).
  real*8, allocatable :: old_dpsibnd_vec(:)              !< Previous delta Psi values required for time-stepping with zeta/=0

  !> @name JOREK vacuum response matrices
  !! Response matrices derived from STARWALL response (w=wall, p=plasma)
  real*8, allocatable :: response_m_a(:,:)               !< \f$\hat{A}\f$ in the documentation
  real*8, allocatable :: response_d_b(:)                 !< \f$\hat{B}\f$ in the documentation
  real*8, allocatable :: response_d_c(:)                 !< \f$\hat{C}\f$ in the documentation
  real*8, allocatable :: response_m_d(:,:)               !< \f$\hat{D}\f$ in the documentation
  real*8, allocatable :: response_m_e(:,:)               !< \f$\hat{E}\f$ in the documentation
  real*8, allocatable :: response_m_f(:,:)               !< \f$\hat{F}\f$ in the documentation
  real*8, allocatable :: response_m_g(:,:)               !< \f$\hat{G}\f$ in the documentation
  real*8, allocatable :: response_m_h(:,:)               !< \f$\hat{H}\f$ in the documentation
  real*8, allocatable :: response_m_j(:,:)               !< \f$\hat{J}\f$ in the documentation
  real*8, allocatable :: response_m_k(:,:)               !< \f$\hat{K}\f$ in the documentation
  real*8, allocatable :: response_m_l(:,:)               !< \f$\hat{L}\f$ in the documentation
  real*8, allocatable :: response_m_eq(:,:)              !< Response matrix for vacuum_equil

  !> @name Equilibrium coil contributions
  integer             :: n_coils                         !< number of poloidal field coils
  real*8, allocatable :: I_coils(:)                      !< coil currents 
  real*8              :: vertical_FB                     !< a variable for the feedback control of the plasma's vertical position
  real*8, allocatable :: bext_tan(:,:)                   !< external tangential field
  real*8, allocatable :: bext_nor(:,:)                   !< external normal field
  real*8, allocatable :: bext_psi(:,:)                   !< external poloidal flux
  
  ! ### various variables, some need to be removed
  real*8, allocatable :: R_coils(:), Z_coils(:)          ! ### old
  real*8, allocatable :: dR_coils(:), dZ_coils(:)        ! ### old
  real*8, allocatable :: coil_voltages(:)                !< Coil voltages
  
  type :: t_starwall_response
    integer :: n_bnd
    integer :: nd_bez
    integer :: ncoil
    integer :: npot_w
    integer :: n_w
    integer :: ntri_w
    integer :: n_tor
    integer, allocatable :: i_tor(:)
    real*8,  allocatable :: d_yy(:)
    real*8,  allocatable :: a_ye(:,:)
    real*8,  allocatable :: a_ey(:,:)
    real*8,  allocatable :: a_ee(:,:)
    real*8,  allocatable :: a_id(:,:)
    real*8,  allocatable :: a_nw(:,:)
    real*8,  allocatable :: s_ww(:,:)
    real*8,  allocatable :: s_ww_inv(:,:)
    real*8,  allocatable :: xyzpot_w(:,:)
    integer, allocatable :: jpot_w(:,:)
  end type t_starwall_response
  type(t_starwall_response) :: sr
  
  
  contains
  
  
  
  !> Preset freeboundary related input parameters to reasonable default values.
  subroutine vacuum_preset(my_id, freeboundary_equil, freeboundary, resistive_wall)
    
    integer, intent(in)    :: my_id
    logical, intent(out)   :: freeboundary_equil, freeboundary, resistive_wall
    
    ! --- Preset namelist input parameters.
    freeboundary_equil   = .false.
    freeboundary         = .false.
    resistive_wall       = .false.
    wall_resistivity     = 0.d0
    
  end subroutine vacuum_preset
  
  
  
  !> Initialize vacuum, ensure consistency of input parameters.
  subroutine vacuum_init(my_id, freeboundary_equil, freeboundary, resistive_wall)
    
    integer, intent(in)    :: my_id
    logical, intent(inout) :: freeboundary_equil, freeboundary, resistive_wall
    
    ! --- Make input parameters consistent.
    freeboundary   = freeboundary .or. freeboundary_equil
    resistive_wall = freeboundary .and. resistive_wall
        
    ! --- Initialize some variables.
    sr%n_bnd  = 0
    sr%nd_bez = 0
    sr%ncoil  = 0
    sr%npot_w = 0
    sr%n_w    = 0
    sr%ntri_w = 0
    sr%n_tor  = 0
    
  end subroutine vacuum_init
  
  
  
  !> Allows to decide if free- or fixed-boundary conditions apply for a certain combination of
  !! toroidal harmonic i_tor and variable i_var. This is important, as the STARWALL response
  !! must not be provided for all toroidal modes and fixed boundary conditions remain implemented
  !! for some of the JOREK variables.
  logical function is_freebound(i_tor, i_var)
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: i_tor
    integer, intent(in) :: i_var
    
    ! --- Local variables
    integer :: i
    
    ! --- Free boundary conditions only if STARWALL response provided for this toroidal harmonic
    is_freebound = .false.
    do i = 1, sr%n_tor
      is_freebound = is_freebound .or. ( i_tor == sr%i_tor(i) )
    end do
    
    ! --- Free boundary conditions only for certain variables
    is_freebound = is_freebound .and. ( (i_var == ivar_j) .or. (i_var == ivar_psi) )
    
  end function is_freebound
  
  
  
  !> Import some vacuum-related data from the restart file  
  !!
  !! @todo Does not work currently if variable freeboundary is changed between export and import!
  subroutine import_restart_vacuum(file_handle, freeboundary, resistive_wall)
    
    ! --- Routine parameters
    integer, intent(in) :: file_handle
    logical, intent(in) :: freeboundary
    logical, intent(in) :: resistive_wall
    
    ! --- Local variables
    logical :: resistive_wall_rst
    integer :: ierr
    
    if ( freeboundary ) then
      
      read(file_handle, iostat=ierr) resistive_wall_rst
      if ( ierr /= 0 ) then
        write(*,*) 'WARNING: Restarting a simulation with freeboundary=.t. which was run with'
        write(*,*) '  freeboundary=.f. so far.'
        return
      else if ( resistive_wall .neqv. resistive_wall_rst ) then
        write(*,*) 'ERROR: It is currently not possible to restart a JOREK simulation with a'
        write(*,*) '  modified setting for resistive_wall.'
        stop
      end if
      
      if ( resistive_wall ) then
        
        read(file_handle) n_wall_curr, n_dof_starwall
        
        if ( allocated(wall_curr) ) deallocate(wall_curr)
        allocate( wall_curr(n_wall_curr) )
        read(file_handle) wall_curr(:)
        
        if ( allocated(dwall_curr) ) deallocate(dwall_curr)
        allocate( dwall_curr(n_wall_curr) )
        read(file_handle) dwall_curr(:)
        
        if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
        allocate( old_dpsibnd_vec(n_dof_starwall) )
        old_dpsibnd_vec = 0.d0 !###
        
        if ( vacuum_debug .and. resistive_wall ) then
          write(*,*) 'DEBUG: Checksums'
          write(*,*) 'wall_curr', sum(abs(wall_curr))
          write(*,*) 'dwall_curr', sum(abs(dwall_curr))
          write(*,*) 'END: Checksums'
        end if
        
        wall_curr_initialized = .true.
        
      end if
      
    end if
    
    if ( vacuum_debug .and. resistive_wall ) then
      write(*,*) 'DEBUG: Checksums'
      if ( allocated(wall_curr)  ) write(*,*) 'wall_curr ', sum(abs(wall_curr))
      if ( allocated(dwall_curr) ) write(*,*) 'dwall_curr', sum(abs(dwall_curr))
      write(*,*) 'END: Checksums'
    end if
    
  end subroutine import_restart_vacuum


  !> Import some vacuum-related data from the HDF5 restart file  
  !!
  !! @todo Does not work currently if variable freeboundary is changed between export and import!
  subroutine import_HDF5_restart_vacuum(file_id, freeboundary, resistive_wall)
    
#ifdef USE_HDF5
    use hdf5
    use HDF5_io_module

    ! --- Routine parameters
    integer(HID_T), intent(in) :: file_id
    logical,        intent(in) :: freeboundary
    logical,        intent(in) :: resistive_wall
    
    ! --- Local variables
    logical   :: resistive_wall_rst
    character :: t_resistive_wall
    
    if ( freeboundary ) then
       call HDF5_char_reading(file_id,t_resistive_wall,"resistive_wall")
       if (t_resistive_wall == "T") then
          resistive_wall_rst = .TRUE.
       else
          resistive_wall_rst = .FALSE.
       endif
       if ( .not. resistive_wall_rst ) then
          write(*,*) 'WARNING: Restarting a simulation with freeboundary=.t. which was run with'
          write(*,*) '  freeboundary=.f. so far.'
          return
       else if ( resistive_wall .neqv. resistive_wall_rst ) then
          write(*,*) 'ERROR: It is currently not possible to restart a JOREK simulation with a'
          write(*,*) '  modified setting for resistive_wall.'
          stop
       end if
       
       if ( resistive_wall ) then
          call HDF5_integer_saving(file_id,n_wall_curr,"n_wall_curr")
          call HDF5_integer_saving(file_id,n_dof_starwall,"n_dof_starwall")
          
          if ( allocated(wall_curr) ) deallocate(wall_curr)
          allocate( wall_curr(n_wall_curr) )
          call HDF5_array1D_reading(file_id,wall_curr,"wall_curr")
          
          if ( allocated(dwall_curr) ) deallocate(dwall_curr)
          allocate( dwall_curr(n_wall_curr) )
          call HDF5_array1D_reading(file_id,dwall_curr,"dwall_curr")
          
          if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
          allocate( old_dpsibnd_vec(n_dof_starwall) )
          old_dpsibnd_vec = 0.d0 !###
          
          if ( vacuum_debug .and. resistive_wall ) then
             write(*,*) 'DEBUG: Checksums'
             write(*,*) 'wall_curr', sum(abs(wall_curr))
             write(*,*) 'dwall_curr', sum(abs(dwall_curr))
             write(*,*) 'END: Checksums'
          end if
          
          wall_curr_initialized = .true.
          
       end if
       
    end if
     
    if ( vacuum_debug .and. resistive_wall ) then
       write(*,*) 'DEBUG: Checksums'
       if ( allocated(wall_curr)  ) write(*,*) 'wall_curr ', sum(abs(wall_curr))
       if ( allocated(dwall_curr) ) write(*,*) 'dwall_curr', sum(abs(dwall_curr))
       write(*,*) 'END: Checksums'
    end if

#endif    

  end subroutine import_HDF5_restart_vacuum
  
  
  !> Export some vacuum-related data to the restart file  
  subroutine export_restart_vacuum(file_handle, freeboundary, resistive_wall)
    
    ! --- Routine parameters
    integer, intent(in) :: file_handle
    logical, intent(in) :: freeboundary
    logical, intent(in) :: resistive_wall
    
    if ( freeboundary ) then
      
      write(file_handle) resistive_wall
      if ( resistive_wall ) then
        
        if ( (.not. allocated(wall_curr)) .or. (.not. allocated(dwall_curr)) .or.                    &
          (.not. allocated(old_dpsibnd_vec)) ) then
          write(*,*) 'ERROR in mod_vacuum.f90:export_restart_vacuum: Arrays not allocated.'
          stop
        end if
        
        write(file_handle) n_wall_curr, n_dof_starwall
        
        write(file_handle) wall_curr(:)
        
        write(file_handle) dwall_curr(:)
        
      end if
      
    end if
    
    if ( vacuum_debug .and. resistive_wall ) then
      write(*,*) 'DEBUG: Checksums'
      if ( allocated(wall_curr)  ) write(*,*) 'wall_curr ', sum(abs(wall_curr))
      if ( allocated(dwall_curr) ) write(*,*) 'dwall_curr', sum(abs(dwall_curr))
      write(*,*) 'END: Checksums'
    end if
    
  end subroutine export_restart_vacuum
  
  
  !> Export some vacuum-related data to the HDF5 restart file  
  subroutine export_HDF5_restart_vacuum(file_id, freeboundary, resistive_wall)
    
#ifdef USE_HDF5
    use hdf5
    use HDF5_io_module

    ! --- Routine parameters
    integer(HID_T), intent(in) :: file_id
    logical,        intent(in) :: freeboundary
    logical,        intent(in) :: resistive_wall
    
    ! --- Local variables
    character           :: t_resistive_wall

    if ( freeboundary ) then
       if (resistive_wall) then
          t_resistive_wall = "T"
       else
          t_resistive_wall = "F"
       end if
       call HDF5_char_saving(file_id,t_resistive_wall,"resistive_wall"//char(0))

       if ( resistive_wall ) then
          if ( (.not. allocated(wall_curr)) .or. (.not. allocated(dwall_curr)) .or. &
               (.not. allocated(old_dpsibnd_vec)) ) then
             write(*,*) 'ERROR in mod_vacuum.f90:export_restart_vacuum: Arrays not allocated.'
             stop
          end if
          call HDF5_integer_saving(file_id,n_wall_curr,"n_wall_curr"//char(0))
          call HDF5_integer_saving(file_id,n_dof_starwall,"n_dof_starwall"//char(0))
          call HDF5_array1D_saving(file_id,wall_curr,n_wall_curr,"wall_curr"//char(0))
          call HDF5_array1D_saving(file_id,dwall_curr,n_wall_curr,"dwall_curr"//char(0))
       end if
    end if
    
    if ( vacuum_debug .and. resistive_wall ) then
       write(*,*) 'DEBUG: Checksums'
       if ( allocated(wall_curr)  ) write(*,*) 'wall_curr ', sum(abs(wall_curr))
       if ( allocated(dwall_curr) ) write(*,*) 'dwall_curr', sum(abs(dwall_curr))
       write(*,*) 'END: Checksums'
    end if
#endif

  end subroutine export_HDF5_restart_vacuum
  
  
  !> Broadcast vacuum information between MPI processes
  subroutine broadcast_vacuum(my_id, resistive_wall)
    use mpi_mod
    implicit none
    
    
    ! --- Routine parameters
    integer, intent(in) :: my_id
    logical, intent(in) :: resistive_wall
    
    ! --- Local variables
    integer :: ierr
    
    call MPI_BCAST(n_dof_starwall,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
    if ( resistive_wall ) then
      
      call MPI_BCAST(n_wall_curr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      
      if ( n_wall_curr == 0 ) return
      
      if ( my_id /= 0 ) then
        if ( allocated(wall_curr) ) deallocate(wall_curr)
        allocate( wall_curr(n_wall_curr) )
        
        if ( allocated(dwall_curr) ) deallocate(dwall_curr)
        allocate( dwall_curr(n_wall_curr) )
        
        if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
        allocate( old_dpsibnd_vec(n_dof_starwall) )
      end if
      
      call MPI_BCAST(wall_curr,n_wall_curr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dwall_curr,n_wall_curr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(old_dpsibnd_vec,n_dof_starwall,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      
    end if
    
  end subroutine broadcast_vacuum
  
  
  
end module vacuum
