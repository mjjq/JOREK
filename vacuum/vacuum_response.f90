!> Implements the interplay of the plasma with a conducting wall.
!!
!! The plasma-wall interaction is characterized by vacuum response matrices which are calculated by
!! the STARWALL code and imported into JOREK by the routine read_starwall_response(). The matrices
!! allow to express the magnetic field parallel to the interface (boundary of the JOREK domain) in
!! terms of the poloidal flux at the interface and the wall currents. The vacuum
!! response enters into the boundary integral of the current equation which vanishes for fixed
!! boundary conditions. The boundary integral is implemented in the routine
!! vacuum_boundary_integral().
!!
!! @note The variable s in a boundary element may correspond to s or t of the respective
!! 2D element depending on element orientation.
#include "pastix_fortran.h"
module vacuum_response
  
  use vacuum
  
  implicit none
  
  
  
  contains
  
  
  
  !> Is the filehandle associated with a formatted file?
  logical function is_formatted(filehandle)
    implicit none
    ! --- Routine parameters
    integer,          intent(in) :: filehandle
    ! --- Local variables
    character(len=64) :: format_type
    
    inquire(filehandle,form=format_type)
    is_formatted = ( trim(format_type) == 'FORMATTED' )
  end function is_formatted
  
  
  
  !> Get the vacuum response for an ideal or resistive wall.
  subroutine get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list,                    &
    freeboundary_equil, resistive_wall)
  
    use mod_parameters,     only: n_tor
    use data_structure, only: type_node_list, type_bnd_element_list, type_bnd_node_list
    use mpi_mod

    implicit none    
    
    integer,                     intent(in) :: my_id              !< MPI proc ID
    type(type_node_list),        intent(in) :: node_list          !< List of boundary nodes
    type(type_bnd_element_list), intent(in) :: bnd_elm_list       !< List of boundary elements
    type(type_bnd_node_list),    intent(in) :: bnd_node_list      !< List of boundary nodes
    logical,                     intent(in) :: freeboundary_equil !< Use free boundary equilibrium?
    logical,                     intent(in) :: resistive_wall     !< Resistive or ideal wall?

    integer :: i,j, ierr, dim
    logical :: exists
 
    do i=1, bnd_node_list%n_bnd_nodes
      exists = .false.
      do j=1, i-1
        if (bnd_node_list%bnd_node(i)%index_jorek .eq. bnd_node_list%bnd_node(j)%index_jorek) then
          exists  = .true.
          exit
        endif
      enddo
      if (.not. exists) then
        n_dof_bnd = n_dof_bnd + bnd_node_list%bnd_node(i)%n_dof ! Number of boundary degrees of freedom per harmonic
      endif
    enddo

    write(*,'(A,i5)') 'total number of degrees of freedom on the boundary : ',n_dof_bnd
    
    ! --- Write out the boundary information for STARWALL.
    if (my_id == 0) call export_boundary(node_list, bnd_elm_list, bnd_node_list)
    
    if ( vacuum_debug .and. (my_id == 0) ) call log_starwall_response(sr)
    
    if ( my_id == 0 ) call read_starwall_response(sr, 'starwall-response.dat')

    call broadcast_starwall_response(my_id, sr)
    
    if ( my_id == 0 ) call log_starwall_response(sr)
    
  end subroutine get_vacuum_response
  
  
  
  
  
  
  
  !> Read an integer parameter from a the STARWALL response file.
  integer function read_intparam(filehandle, parameter_name)
    
    implicit none
    
    ! --- Routine parameters
    integer,          intent(in) :: filehandle
    character(len=*), intent(in) :: parameter_name
    
    ! --- Local variables
    character(len=12) :: marker
    character(len=24) :: name
    integer           :: ierr
    
    if ( is_formatted(filehandle) ) then
      read(filehandle,'(A12,A24,I12)',iostat=ierr) marker, name, read_intparam
    else
      read(filehandle, iostat=ierr) marker, name, read_intparam
    end if
    
    if ( (ierr /= 0) .or. (trim(adjustl(marker)) /= '#@intparam') .or. (trim(adjustl(name)) /= trim(parameter_name)) ) then
      write(*,*) 'ERROR: Could not read parameter "', trim(parameter_name) ,'" from STARWALL response.'
      stop
    end if
    
    if ( vacuum_debug ) write(*,'(3x,"Read: ",A24,"=",I12)',iostat=ierr) name, read_intparam
    
  end function read_intparam
  
  
  
  
  
  
  !> Read an array from the STARWALL respone file
  subroutine read_array(filehandle, array_name, dim, int1d, int2d, float1d, float2d)
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: filehandle
    character(len=*), intent(in) :: array_name
    integer, intent(in) :: dim(2)
    integer, allocatable, optional, intent(inout)   :: int1d(:)
    integer, allocatable, optional, intent(inout)   :: int2d(:,:)
    real*8,  allocatable, optional, intent(inout)   :: float1d(:)
    real*8,  allocatable, optional, intent(inout)   :: float2d(:,:)
    
    ! --- Local variables
    character(len=12) :: marker
    integer           :: nd, d(2), ierr
    character(len=24) :: name, datatype, requested_type
    logical           :: error
    
    if ( present(int1d) .or. present(int2d) ) then
      requested_type = 'int'
    else
      requested_type = 'float'
    end if
    
    if ( is_formatted(filehandle) ) then
      read(filehandle,'(A12,A24,I12,A24,2I12)',iostat=ierr) marker, name, nd, datatype, d
    else
      read(filehandle,iostat=ierr) marker, name, nd, datatype, d
    end if
    marker   = adjustl(marker)
    name     = adjustl(name)
    datatype = adjustl(datatype)
    
    error = ( ierr /= 0 ) .or. ( trim(marker) /= '#@array' ) .or. ( trim(name) /= trim(array_name) )                     &
      .or. ( dim(1) /= d(1) ) .or. ( dim(2) /= d(2) ) .or. ( trim(datatype) /= trim(requested_type) )
    
    if ( error ) then
      write(*,*) 'ERROR: Could not read array ', trim(array_name), ' from STARWALL response.'
      stop
    end if
    
    if ( present(int1d) ) then
      
      if ( allocated(int1d) ) deallocate( int1d )
      allocate( int1d(dim(1)) )
      if ( is_formatted(filehandle) ) then
        read(filehandle,*) int1d(:)
      else
        read(filehandle) int1d(:)
      end if
      
    else if ( present(int2d) ) then
      
      if ( allocated(int2d) ) deallocate( int2d )
      allocate( int2d(dim(1),dim(2)) )
      if ( is_formatted(filehandle) ) then
        read(filehandle,*) int2d(:,:)
      else
        read(filehandle) int2d(:,:)
      end if
      
    else if ( present(float1d) ) then
      
      if ( allocated(float1d) ) deallocate( float1d )
      allocate( float1d(dim(1)) )
      if ( is_formatted(filehandle) ) then
        read(filehandle,*) float1d(:)
      else
        read(filehandle) float1d(:)
      end if
      
    else if ( present(float2d) ) then
      
      if ( allocated(float2d) ) deallocate( float2d )
      allocate( float2d(dim(1),dim(2)) )
      if ( is_formatted(filehandle) ) then
        read(filehandle,'(4ES24.16)') float2d(:,:)
      else
        read(filehandle) float2d(:,:)
      end if
      
    end if
    
    if ( vacuum_debug ) write(*,'(3x,"Read: ",A24,"> type ",A," size ",2I7)') name, trim(datatype), d(1:nd)

    
  end subroutine read_array
  
  
  
  
  
  
  !> Read the STARWALL response matrices from a single file.
  subroutine read_starwall_response(sr, filename)
    
    use constants
    use mod_parameters, only: n_tor, n_period
    
    implicit none
    
    ! --- Routine parameters
    type(t_starwall_response), intent(inout) :: sr
    character(len=*),          intent(in)    :: filename
    
    ! --- Local variables
    integer, parameter :: filehandle = 60
    character(len=512) :: comment
    integer            :: file_version, i, j, i_starw, n, is_sin, err, i_tmp
    real*8             :: r_tmp
    real*8, allocatable :: tmp(:)
    
    ! --- Open file
    !   --- Try to open as unformatted file
    write(*,*) 'Trying to open response as unformatted file...'
    open(filehandle, file=trim(filename), form='unformatted', status='old', action='read', &
      iostat=err)
    if ( err == 0 ) then
      read(filehandle,iostat=err) i_tmp, r_tmp
      if ( (i_tmp/=42) .or. (r_tmp/=42.d0) ) then
        err=-42
        close(filehandle)
      end if
    end if
    
    !   --- Try to open as formatted file
    if ( err /= 0 ) then
      write(*,*) '  ... failed.'
      write(*,*) 'Trying to open response as formatted file...'
      open(filehandle, file=trim(filename), form='formatted', status='old', action='read', &
        iostat=err)
    end if
    
    if ( err /= 0 ) then
      write(*,*) '  ... failed.'
      write(*,*) 'ERROR: STARWALL response file (',trim(filename),') could not be opened.'
      stop
    end if
    write(*,*) '  ... succeeded.'
    
    ! --- Read data from STARWALL response file
    if ( is_formatted(filehandle) ) then
      read(filehandle,'(A)') comment
    else
      read(filehandle) comment
    end if
    
    file_version = read_intparam(filehandle, 'file_version')
    if ( file_version > 1 ) then
      write(*,*) 'ERROR: STARWALL response file version ', file_version, ' is not supported.'
      stop
    end if

    sr%n_bnd  = read_intparam(filehandle, 'n_bnd')
    sr%nd_bez = read_intparam(filehandle, 'nd_bez')
    sr%ncoil  = read_intparam(filehandle, 'ncoil')
    sr%npot_w = read_intparam(filehandle, 'npot_w')
    sr%n_w    = read_intparam(filehandle, 'n_w')
    sr%ntri_w = read_intparam(filehandle, 'ntri_w')
    sr%n_tor  = read_intparam(filehandle, 'n_tor')

    call read_array(filehandle, 'i_tor',    (/sr%n_tor,0/),          int1d=sr%i_tor)
    call read_array(filehandle, 'yy',       (/sr%n_w,0/),            float1d=sr%d_yy)
    call read_array(filehandle, 'ye',       (/sr%n_w,sr%nd_bez/),    float2d=sr%a_ye)
    call read_array(filehandle, 'ey',       (/sr%nd_bez,sr%n_w/),    float2d=sr%a_ey)
    call read_array(filehandle, 'ee',       (/sr%nd_bez,sr%nd_bez/), float2d=sr%a_ee)
    call read_array(filehandle, 's_ww',     (/sr%n_w,sr%n_w/),       float2d=sr%s_ww)
    call read_array(filehandle, 's_ww_inv', (/sr%n_w,sr%n_w/),       float2d=sr%s_ww_inv)
    call read_array(filehandle, 'xyzpot_w', (/sr%npot_w,3/),         float2d=sr%xyzpot_w)
    call read_array(filehandle, 'jpot_w',   (/sr%ntri_w,3/),         int2d=sr%jpot_w)

    close(filehandle)
    if ( vacuum_debug) write(*,*) 'Finished reading import of vacuum response.'

    ! --- Import normalization
    sr%a_ee(:,:) = sr%a_ee(:,:) * 2.d0*PI
    sr%a_ye(:,:) = sr%a_ye(:,:) * 2.d0*PI
    if ( vacuum_debug) write(*,*) 'Applied import normalization.'

    ! --- STARWALL Cartesian coordinates -> JOREK Cartesian coordinates (replace y <-> z)
    allocate( tmp(sr%npot_w) )
    tmp(:)           = sr%xyzpot_w(:,2)
    sr%xyzpot_w(:,2) = sr%xyzpot_w(:,3)
    sr%xyzpot_w(:,3) = tmp(:)
    deallocate( tmp )
    
    ! --- Compute ideal-wall and no-wall response matrices.
    if ( allocated(sr%a_id) ) deallocate(sr%a_id)
    if ( allocated(sr%a_nw) ) deallocate(sr%a_nw)
    allocate( sr%a_id(sr%nd_bez,sr%nd_bez), sr%a_nw(sr%nd_bez,sr%nd_bez) )
    sr%a_nw(:,:) = sr%a_ee(:,:)
    sr%a_id(:,:) = sr%a_ee(:,:) - matmul( sr%a_ey(:,:), sr%a_ye(:,:) )
    
    ! --- Transform STARWALL harmonics to account for periodicity
    j = 0
    do i = 1, sr%n_tor
      i_starw = sr%i_tor(i)
      n       = i_starw / 2
      is_sin  = i_starw - 2 * n
      i_starw = 2 * n/n_period + is_sin
      if ( (mod(n, n_period) /= 0) .or. (i_starw < 1) .or. (i_starw > n_tor) ) then
        write(*,*) 'WARNING: STARWALL harmonic has no JOREK equivalent!'
        write(*,*) 'i_starw    =', sr%i_tor(i)
        write(*,*) 'n_period   =', n_period 
        write(*,*) 'n_tor      =', n_tor
      else
        j = j + 1
        sr%i_tor(j) = i_starw
      end if
    end do
    sr%n_tor = j
    if ( vacuum_debug) write(*,*) 'End of routine read_starwall_response.'

  end subroutine read_starwall_response
  
  
  
  
  
  
  !> Broadcast the STARWALL response matrices to the other MPI procs.
  subroutine broadcast_starwall_response(my_id, sr)
    
  use mpi_mod
    implicit none
    
    
    ! --- Routine parameters
    integer,                   intent(in)    :: my_id
    type(t_starwall_response), intent(inout) :: sr
    
    ! --- Local parameters
    integer :: ierr
    real*8  :: checksum
    
    if ( vacuum_debug ) write(*,*) my_id, 'Entering broadcast_starwall_response.'
    
    ! --- Broadcast parameters.
    call MPI_bcast(sr%n_bnd,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%nd_bez,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%ncoil,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%npot_w,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_w,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%ntri_w,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_tor,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    n_dof_starwall = sr%nd_bez
    n_wall_curr    = sr%n_w
    
    ! --- Allocate matrics.
    if ( my_id /= 0 ) then
      if (allocated(sr%i_tor)   ) deallocate(sr%i_tor);    allocate(sr%i_tor(sr%n_tor))
      if (allocated(sr%d_yy)    ) deallocate(sr%d_yy);     allocate(sr%d_yy(sr%n_w))
      if (allocated(sr%a_ye)    ) deallocate(sr%a_ye);     allocate(sr%a_ye(sr%n_w,sr%nd_bez))
      if (allocated(sr%a_ey)    ) deallocate(sr%a_ey);     allocate(sr%a_ey(sr%nd_bez,sr%n_w))
      if (allocated(sr%a_ee)    ) deallocate(sr%a_ee);     allocate(sr%a_ee(sr%nd_bez,sr%nd_bez))
      if (allocated(sr%a_id)    ) deallocate(sr%a_id);     allocate(sr%a_id(sr%nd_bez,sr%nd_bez))
      if (allocated(sr%a_nw)    ) deallocate(sr%a_nw);     allocate(sr%a_nw(sr%nd_bez,sr%nd_bez))
      if (allocated(sr%s_ww)    ) deallocate(sr%s_ww);     allocate(sr%s_ww(sr%n_w,sr%n_w))
      if (allocated(sr%s_ww_inv)) deallocate(sr%s_ww_inv); allocate(sr%s_ww_inv(sr%n_w,sr%n_w))
      if (allocated(sr%xyzpot_w)) deallocate(sr%xyzpot_w); allocate(sr%xyzpot_w(sr%npot_w,3))
      if (allocated(sr%jpot_w)  ) deallocate(sr%jpot_w);   allocate(sr%jpot_w(sr%ntri_w,3))
    end if
    
    ! --- Broadcast matrices.
    call MPI_bcast(sr%i_tor,    sr%n_tor,            MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%d_yy,     sr%n_w,              MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%a_ye,     sr%n_w*sr%nd_bez,    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%a_ey,     sr%nd_bez*sr%n_w,    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%a_ee,     sr%nd_bez*sr%nd_bez, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%a_id,     sr%nd_bez*sr%nd_bez, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%a_nw,     sr%nd_bez*sr%nd_bez, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%s_ww,     sr%n_w*sr%n_w,       MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%s_ww_inv, sr%n_w*sr%n_w,       MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%xyzpot_w, sr%npot_w*3,         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%jpot_w,   sr%ntri_w*3,         MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
    
    if ( vacuum_debug ) write(*,'("Checksum",I4,ES24.16)') my_id, sum(sr%i_tor) + sum(sr%d_yy)     &
       + sum(sr%a_ye) + sum(sr%a_ey) + sum(sr%a_ee) + sum(sr%a_id) + sum(sr%a_nw) + sum(sr%s_ww)   &
       + sum(sr%s_ww_inv ) + sum(sr%xyzpot_w) + sum(sr%jpot_w) + sr%n_bnd + sr%nd_bez + sr%ncoil   &
       + sr%npot_w + sr%n_w + sr%ntri_w + sr%n_tor
    
    if ( vacuum_debug ) write(*,*) my_id, 'Exiting broadcast_starwall_response.'
    
  end subroutine broadcast_starwall_response
  
  
  
  
  
  
  !> Write out information about the STARWALL response matrices.
  subroutine log_starwall_response(sr)
    
    use mod_parameters, only: n_period
    
    implicit none
    
    type(t_starwall_response), intent(in) :: sr
    
    32 format(3x,77('-'))
    33 format(3x,a,i8)
    34 format(3x,'sum(',a,')=',es24.16)
    35 format(3x,'sum(',a,')=',i24)
    36 format(3x,'sum(',a,')= ---not allocated---')
    write(*,*)
    write(*,32)
    write(*,33) 'STARWALL RESPONSE INFORMATION:'
    write(*,32)
    write(*,33) 'n_bnd =', sr%n_bnd
    write(*,33) 'nd_bez=', sr%nd_bez
    write(*,33) 'ncoil =', sr%ncoil
    write(*,33) 'npot_w=', sr%npot_w
    write(*,33) 'n_w   =', sr%n_w
    write(*,33) 'ntri_w=', sr%ntri_w
    write(*,33) 'n_tor =', sr%n_tor
    if (allocated(sr%i_tor)) write(*,33) 'i_tor ='//trim(modes_to_str(sr%i_tor,sr%n_tor,n_period))
    if ( vacuum_debug ) then
      write(*,32)
      if (allocated(sr%i_tor   )) then; write(*,35) 'i_tor   ', sum(sr%i_tor   ); else; write(*,36) 'i_tor   '; end if
      if (allocated(sr%d_yy    )) then; write(*,34) 'd_yy    ', sum(sr%d_yy    ); else; write(*,36) 'd_yy    '; end if
      if (allocated(sr%a_ye    )) then; write(*,34) 'a_ye    ', sum(sr%a_ye    ); else; write(*,36) 'a_ye    '; end if
      if (allocated(sr%a_ey    )) then; write(*,34) 'a_ey    ', sum(sr%a_ey    ); else; write(*,36) 'a_ey    '; end if
      if (allocated(sr%a_ee    )) then; write(*,34) 'a_ee    ', sum(sr%a_ee    ); else; write(*,36) 'a_ee    '; end if
      if (allocated(sr%a_id    )) then; write(*,34) 'a_id    ', sum(sr%a_id    ); else; write(*,36) 'a_id    '; end if
      if (allocated(sr%a_nw    )) then; write(*,34) 'a_nw    ', sum(sr%a_nw    ); else; write(*,36) 'a_nw    '; end if
      if (allocated(sr%s_ww    )) then; write(*,34) 's_ww    ', sum(sr%s_ww    ); else; write(*,36) 's_ww    '; end if
      if (allocated(sr%s_ww_inv)) then; write(*,34) 's_ww_inv', sum(sr%s_ww_inv); else; write(*,36) 's_ww_inv'; end if
      if (allocated(sr%xyzpot_w)) then; write(*,34) 'xyzpot_w', sum(sr%xyzpot_w); else; write(*,36) 'xyzpot_w'; end if
      if (allocated(sr%jpot_w  )) then; write(*,35) 'jpot_w  ', sum(sr%jpot_w  ); else; write(*,36) 'jpot_w  '; end if
    end if
    write(*,32)
    write(*,*)
    
  end subroutine log_starwall_response
  
  
  
  
  
  
  !> Write out resistive-wall data as a VTK-file.
  !!
  !! Scalar quantities (at wall triangle nodes):
  !! * pot_w: Wall current potentials
  !! * dpot_w: Change of wall current potentials in previous time-step
  !!
  !! Vector quantities (at triangles):
  !! * jsurf_w: Surface currents on the wall
  subroutine write_wall_vtk(index, resistive_wall)
    
    use phys_module, only: nout
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: index !< Time step index
    logical, intent(in) :: resistive_wall
    
    ! --- Local variables
    real*8              :: phi1, phi2, phi3, r1(3), r2(3), r3(3), r21(3), r32(3), r21_cross_r32(3)
    integer             :: filehandle = 60, i
    character(len=18)   :: filename
    real*8, allocatable :: tripot_w(:)
    
    if ( mod(index,nout) /= 0 ) return
    
    ! --- VTK file header
    write(filename,'(A,I5.5,A)') 'wallcurr.',index,'.vtk'
    open(filehandle, file=filename, status='replace', action='write')
    140 format(a)
    141 format(a,i8,a)
    142 format(3es16.8)
    143 format(a,2i8)
    144 format(4i8)
    write(filehandle,140) '# vtk DataFile Version 2.0'
    write(filehandle,140) 'testdata'
    write(filehandle,140) 'ASCII'
    write(filehandle,140) 'DATASET POLYDATA'
    
    ! --- Triangle node positions
    write(filehandle,141) 'POINTS', sr%npot_w, ' float'
    do i = 1, sr%npot_w
      write(filehandle,142) sr%xyzpot_w(i,:)
    end do
    
    ! --- Node indices corresponding to triangles
    write(filehandle,143) 'POLYGONS', sr%ntri_w, sr%ntri_w * 4
    do i = 1, sr%ntri_w
      write(filehandle,144) 3, sr%jpot_w(i,:) - 1
    end do
    
    ! --- Wall current potentials
    write(filehandle,141) 'POINT_DATA', sr%npot_w
    write(filehandle,140) 'SCALARS pot_w float'
    write(filehandle,140) 'LOOKUP_TABLE default'
    call reconstruct_triangle_potentials(tripot_w, wall_curr)
    do i = 1, sr%npot_w
      write(filehandle,142) tripot_w(i)
    end do
    
    ! --- Change of wall current potentials in previous time-step
    write(filehandle,140) 'SCALARS dpot_w float'
    write(filehandle,140) 'LOOKUP_TABLE default'
    call reconstruct_triangle_potentials(tripot_w, dwall_curr)
    do i = 1, sr%npot_w
      write(filehandle,142) tripot_w(i)
    end do
    
    ! --- Wall current vectors
    write(filehandle,141) 'CELL_DATA', sr%ntri_w
    write(filehandle,140) 'VECTORS jsurf_w float'
    call reconstruct_triangle_potentials(tripot_w, wall_curr)
    do i = 1, sr%ntri_w
      ! --- Wall potential at triangle nodes
      phi1   = tripot_w(sr%jpot_w(i,1))
      phi2   = tripot_w(sr%jpot_w(i,2))
      phi3   = tripot_w(sr%jpot_w(i,3))
      ! --- Position of triangle nodes
      r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
      r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
      r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
      r21(:) = r1(:)-r2(:)
      r32(:) = r2(:)-r3(:)
      r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
        r21(1)*r32(2) - r21(2)*r32(1) /)
      write(filehandle,142) ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sum(r21_cross_r32**2)
    end do
    
    ! --- Close file, clean up
    if ( allocated(tripot_w) ) deallocate( tripot_w )
    close(filehandle)
    
  end subroutine write_wall_vtk
  
  
  
  
  
  
  !> Implements the boundary integral in the current equation which vanishes for fixed boundary
  !! conditions.
  !!
  !! The magnetic field parallel to the interface (boundary of JOREK computational domain) is
  !! expressed by the STARWALL vacuum response in terms of the poloidal magnetic field at the
  !! interface (ideal and resistive wall) and the wall currents (resistive wall).
  subroutine vacuum_boundary_integral(my_id, bnd_node_list, node_list, bnd_elm_list,               &
    freeboundary_equil, resistive_wall, index_min, index_max, rhs_loc, tstep, index_now)
      
    use data_structure, only: type_node_list, type_bnd_node_list, type_bnd_element_list, type_bnd_element
    use mod_parameters,     only: n_plane, n_var, n_tor
    use gauss,          only: n_gauss, xgauss, wgauss
    use global_distributed_matrix, only: irn_glob, jcn_glob, a_glob, ndof_glob, det_row_col, det_sparse_pos
    use basis_at_gaussian, only: H1, H1_s, HZ
    use phys_module, only: t_now, t_start
    
    implicit none
    
    ! --- Routine parameters
    integer,                     intent(in)    :: my_id                !< MPI process ID
    type(type_node_list),        intent(in)    :: node_list            !< List of grid nodes
    type(type_bnd_node_list),    intent(in)    :: bnd_node_list        !< List of boundary grid nodes
    type(type_bnd_element_list), intent(in)    :: bnd_elm_list         !< List of boundary elements
    logical,                     intent(in)    :: freeboundary_equil   !< Use free boundary equilibrium?
    logical,                     intent(in)    :: resistive_wall       !< Resistive or ideal wall?
    integer,                     intent(in)    :: index_min, index_max !< Responsibility of MPI proc
    real*8,                      intent(inout) :: rhs_loc(ndof_glob)   !< Part of RHS of MPI proc
    real*8,                      intent(in)    :: tstep                !< delta t, timestep
    integer,                     intent(in)    :: index_now            !< Current timestep index
    
    ! --- Local variables
    real*8, allocatable :: psibnd_vec(:)    ! Vector of the values of Psi at the boundary
    real*8, allocatable :: dpsibnd_vec(:)   ! Vector of the values of deltaPsi at the boundary
    real*8, allocatable :: psibnd_coils(:)  ! Vector of the values of Psi_coil at the boundary
    real*8   :: amat_contrib, rhs_contrib   ! Vacuum response contribution to lhs and rhs
    real*8   :: testfunc_l                  ! j^*_l in documentation
    real*8   :: basfunc_i                   ! b_i in documentation
    real*8   :: dA                          ! factor from definition of dA
    real*8   :: x_s(n_gauss), y_s(n_gauss)  ! values of dR/ds and dZ/ds at Gaussian points
    real*8   :: common_prefactor
    integer  :: m_bndelem                   ! Boundary element index
    type(type_bnd_element) :: bndelem_m     ! Boundary element corresponding to index m_bndelem
    integer  :: ms                          ! Gauss point index
    integer  :: m_plane                     ! Toroidal plane index
    integer  :: sparsepos_jp, sparsepos_pp  ! Position of lhs contribution in the sparse matrix
    !   --- Test function related quantities
    integer  :: l_vertex, l_dof, l_dir, l_node, l_node_bnd, l_index, l_tor, l_row_j, l_row_psi
    real*8   :: l_size
    !   --- Quantities related to the boundary dof at which response is calculated
    integer  :: i_vertex, i_dof, i_dir, i_node, i_node_bnd, i_index, i_starwall, i_tor, i_resp, i_resp_0
    real*8   :: i_size
    !   --- Quantities related to the boundary dof contributing to the response
    integer  :: j_dof, j_dir, j_node, j_node_bnd, j_index, j_starwall, j_tor, j_col_psi, j_resp

    integer  :: i_resp_old, j_resp_old
#ifdef __GFORTRAN__
    real*8 :: wgauss_copy(4)
#endif
    !integer :: rate, t0, t1 !### timing ###
    logical, save  :: PF_perturbation = .true. 

    if ( vacuum_debug ) write(*,*) my_id, 'Before:', sum(abs(rhs_loc)), sum(abs(A_glob))
    
    ! --- Determine vectors of the psi and deltapsi boundary values.
    call det_psibnd_vec(bnd_node_list, node_list, psibnd_vec, dpsibnd_vec, psibnd_coils)
    
    ! --- Update the derived response matrices
    call update_response(tstep, freeboundary_equil, resistive_wall)
    
    ! --- Perform the time-stepping for the wall currents.
    if ( resistive_wall .and. (index_now>1) ) call evolve_wall_currents(my_id, psibnd_vec, dpsibnd_vec)
    
    if ( vacuum_debug ) then
      write(*,*) my_id, 'psibnd_vec:  ', sum(abs(psibnd_vec)), sum(psibnd_vec)
      write(*,*) my_id, 'dpsibnd_vec: ', sum(abs(dpsibnd_vec)), sum(dpsibnd_vec)
    end if

    if ( my_id == 0 ) call boundary_check()
    
    !-- Add perturbation in PF coils currents to speed-up VDEs
    if (t_start .gt. PF_pert_start_time)  PF_perturbation = .false.
    if (PF_perturbation .and. (t_now .ge. PF_pert_start_time) ) then
      if ( my_id == 0 ) write(*,*) 'Perturbing PF coil currents'
      I_coils(1:n_coils) = I_coils(1:n_coils) + coils0(1:n_coils)%pert
      PF_perturbation    = .false.
    endif
    
#ifdef __GFORTRAN__
    wgauss_copy(1:4) = wgauss(1:4)
#endif
    ! --- Sum over boundary elements
    !$omp parallel do                                                                              &
    !$omp default(none)                                                                            &
    !$omp shared(a_glob, rhs_loc, bnd_elm_list, bnd_node_list, node_list, index_min, index_max,    &
    !$omp   response_m_e, response_m_f, response_m_g, response_m_h, response_m_j, H1, HZ, sr,      &
    !$omp   bext_tan, I_coils, wall_curr, dwall_curr, psibnd_vec, dpsibnd_vec, psibnd_coils,       &
    !$omp   starwall_equil_coils,                                                                  &
#ifdef __GFORTRAN__
    !$omp   wgauss_copy,                                                                                &
#endif
    !$omp   resistive_wall)     &
    !$omp private(m_bndelem, bndelem_m, x_s, y_s, l_vertex, l_dof, l_node, l_dir, l_node_bnd,      &
    !$omp   l_index, l_size, l_tor, l_row_j, l_row_psi, ms, dA, m_plane, common_prefactor,         &
    !$omp   testfunc_l, i_vertex, i_dof, i_node, i_dir, i_node_bnd, i_index, i_size, i_starwall,   &
    !$omp   i_tor, i_resp, i_resp_old, i_resp_0, basfunc_i, j_node_bnd, j_dof, j_node, j_dir, j_index,         &
    !$omp   j_starwall, j_tor, j_resp, j_resp_old,j_col_psi, sparsepos_jp, sparsepos_pp, amat_contrib,        &
    !$omp   rhs_contrib)
    L_MB: do m_bndelem = 1, bnd_elm_list%n_bnd_elements

      bndelem_m = bnd_elm_list%bnd_element(m_bndelem)

      ! --- Determine the values of R,s and Z,s at the Gaussian points.
      call det_coord_bnd(bndelem_m, node_list, R_S=x_s, Z_S=y_s)

      ! --- Select a test function (the weak form equation must hold for every test function)
      L_LV: do l_vertex = 1, 2 ! (loop over nodes in element m_bndelem)

        L_LD: do l_dof = 1, 2 ! (loop over node dofs)

          l_node      = bndelem_m%vertex(l_vertex)
          l_dir       = bndelem_m%direction(l_vertex,l_dof)
          l_node_bnd  = bndelem_m%bnd_vertex(l_vertex)
          l_index     = node_list%node(l_node)%index(l_dir)
          l_size      = bndelem_m%size(l_vertex,l_dof)

          if ( (l_index < index_min) .or. (l_index > index_max) ) cycle ! This MPI proc responsible?

          L_LS: do l_tor = 1, n_tor ! (loop over toroidal harmonics)

            ! --- Determine the row in the main matrix.
            l_row_psi = det_row_col(l_index, ivar_psi, l_tor)
            l_row_j   = det_row_col(l_index, ivar_j,   l_tor)

            ! --- Loop over Gaussian points -- integration in s-direction
            L_MS: do ms = 1, n_gauss

              ! --- Integration factor from the definition of dA:
              !     int dA = sum_{m_bndelem} int ds int dphi sqrt{(R,s)^2 + (Z,s)^2}
              dA = sqrt(x_s(ms)**2 + y_s(ms)**2)

              ! --- Loop over toroidal planes -- integration in phi-direction
              L_MP: do m_plane = 1, n_plane

                ! --- Evaluate test function at current position
                testfunc_l = H1(l_vertex,l_dof,ms) * l_size * HZ(l_tor,m_plane)

                ! --- Sum over boundary dofs at which response is calculated
                L_IV: do i_vertex = 1, 2 ! (loop over nodes in element m_bndelem)

                  i_node      = bndelem_m%vertex(i_vertex)
                  i_node_bnd  = bndelem_m%bnd_vertex(i_vertex)

                  L_ID: do i_dof = 1, 2 ! (loop over node dofs)

                    i_dir       = bndelem_m%direction(i_vertex,i_dof)
                    i_index     = node_list%node(i_node)%index(i_dir)
                    i_size      = bndelem_m%size(i_vertex,i_dof)

                    L_IS: do i_starwall = 1, sr%n_tor ! (loop over STARWALL harmonics)

                      i_tor    = sr%i_tor(i_starwall)

                      i_resp_old   = response_index(i_node_bnd,i_starwall,i_dof)

                      i_resp   = (bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) - 1)*sr%n_tor &
                               + bnd_node_list%bnd_node(i_node_bnd)%n_dof*(i_starwall-1) &
                               + bnd_node_list%bnd_node(i_node_bnd)%index_starwall(i_dof) - bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) + 1

!                      if (i_resp_old .ne. i_resp) write(*,'(A,8i5)') 'PANIC! : ',i_node, i_starwall, i_dof,bnd_node_list%bnd_node(i_node_bnd)%index_starwall,i_resp_old, i_resp
 
                      i_resp_0 = response_index_eq(i_node_bnd,i_dof)

                      ! --- Determine basis function
                      basfunc_i = H1(i_vertex,i_dof,ms) * i_size * HZ(i_tor,m_plane)

#ifdef __GFORTRAN__
                      common_prefactor = wgauss_copy(ms) * dA * testfunc_l * basfunc_i
#else
                      common_prefactor = wgauss(ms) * dA * testfunc_l * basfunc_i
#endif
                      ! --- Sum over boundary dofs contributing to the response
                      L_JB: do j_node_bnd = 1, bnd_node_list%n_bnd_nodes ! (loop over boundary nodes)

                        j_node      = bnd_node_list%bnd_node(j_node_bnd)%index_jorek

                        L_JD: do j_dof = 1, 2 ! (loop over node dofs)

                          j_dir       = bnd_node_list%bnd_node(j_node_bnd)%direction(j_dof)
                          j_index     = node_list%node(j_node)%index(j_dir)

                          L_JS: do j_starwall = 1, sr%n_tor ! (loop over STARWALL harmonics)

                            j_tor  = sr%i_tor(j_starwall)

                            j_resp_old = response_index(j_node_bnd,j_starwall,j_dof)

                            j_resp   = (bnd_node_list%bnd_node(j_node_bnd)%index_starwall(1) - 1)*sr%n_tor &
                                     +  bnd_node_list%bnd_node(j_node_bnd)%n_dof*(j_starwall-1) &
                                     +  bnd_node_list%bnd_node(j_node_bnd)%index_starwall(j_dof)-bnd_node_list%bnd_node(j_node_bnd)%index_starwall(1) + 1

!                      if (j_resp_old .ne. j_resp) write(*,'(A,8i5)') 'PANIC! : ',j_node, j_starwall, j_dof,bnd_node_list%bnd_node(j_node_bnd)%index_starwall,j_resp_old, j_resp

                            ! --- Option to switch off mode coupling due to a 3D wall
                            if ( vacuum_decouple_modes .and. (j_tor /= i_tor) ) cycle

                            ! --- Determine the column in the main matrix
                            j_col_psi = det_row_col(j_index, ivar_psi, j_tor)

                            ! --- Determine the position in the sparse matrix data structure
                            !     which corresponds to the matrix entry at l_row_j, j_col_psi.
                            sparsepos_jp = det_sparse_pos(l_row_j,   j_col_psi, index_min)
                            sparsepos_pp = det_sparse_pos(l_row_psi, j_col_psi, index_min)

                            ! --- Vacuum response contribution to the lhs of the current equation
                            amat_contrib = - common_prefactor * response_m_e(i_resp, j_resp)
                            !$omp atomic
                            A_glob(sparsepos_jp) = A_glob(sparsepos_jp) + amat_contrib

                          end do L_JS
                        end do L_JD
                      end do L_JB

                      ! --- Contribution of vacuum response to the rhs of the current equation

                      rhs_contrib = sum( response_m_h(i_resp,:) * psibnd_vec(:)   )                 &
                                  + sum( response_m_j(i_resp,:) * dpsibnd_vec(:)  )

                      if ( (l_tor == 1) .and. (sr%i_tor(1) == 1) .and. (.not. starwall_equil_coils)) &
                        rhs_contrib = rhs_contrib - sum( bext_tan(i_resp_0, :) * I_coils(:) )       &
                                    - sum( response_m_h(i_resp,:) * psibnd_coils(:) )               

                      if ( resistive_wall ) &
                        rhs_contrib = rhs_contrib + sum( response_m_f(i_resp, :) *  wall_curr(:) )  &
                                                  + sum( response_m_g(i_resp, :) * dwall_curr(:) )   


                      rhs_contrib = rhs_contrib * common_prefactor
                      !$omp atomic
                      rhs_loc(l_row_j) = rhs_loc(l_row_j) + rhs_contrib

                    end do L_IS
                  end do L_ID
                end do L_IV

              end do L_MP

            end do L_MS

          end do L_LS
        end do L_LD
      end do L_LV

    end do L_MB
    !$omp end parallel do
    
    !### timing ###
    !call system_clock(count=t1)
    !write(*,*) 'vacuum_boundary_integral main loop:', real(t1 - t0 ) / real(rate), 's'
    !write(68+my_id,*) real(t1 - t0 ) / real(rate)
    !###
    
    if ( vacuum_debug ) write(*,*) my_id, 'After:', sum(abs(rhs_loc)), sum(abs(A_glob))
  
    if ( allocated(psibnd_vec ) ) deallocate( psibnd_vec  )
    if ( allocated(dpsibnd_vec) ) deallocate( dpsibnd_vec )
    
  end subroutine vacuum_boundary_integral
  
  
  
  
  
  
  !> Determine the values of \f$R\f$, \f$Z\f$, \f$\partial R/\partial s\f$, and
  !! \f$ \partial Z/\partial s \f$ on the Gaussian points of a given boundary element.
  subroutine det_coord_bnd(bndelem, node_list, R, Z, R_s, Z_s)
    
    use gauss,             only: n_gauss, xgauss, wgauss
    use data_structure,    only: type_node, type_bnd_element, type_node_list
    use basis_at_gaussian, only: H1, H1_s, HZ
    
    implicit none
    
    ! --- Routine parameters
    type(type_bnd_element), intent(in)  :: bndelem       !< Boundary element to be considered.
    type(type_node_list),   intent(in)  :: node_list     !< List of grid nodes
    real*8, optional,       intent(out) :: R(n_gauss)    !< Values of R on Gaussian points
    real*8, optional,       intent(out) :: Z(n_gauss)    !< Values of Z on Gaussian points
    real*8, optional,       intent(out) :: R_s(n_gauss)  !< Values of R,s on Gaussian points
    real*8, optional,       intent(out) :: Z_s(n_gauss)  !< Values of Z,s on Gaussian points
    
    ! --- Local variables
    integer         :: k_vertex, k_dof, k_node, k_dir
    real*8          :: k_size
    type(type_node) :: node_k
    
    if ( present(R  ) ) R   = 0.d0
    if ( present(Z  ) ) Z   = 0.d0
    if ( present(R_s) ) R_s = 0.d0
    if ( present(Z_s) ) Z_s = 0.d0
    
    do k_vertex = 1, 2
      do k_dof = 1, 2
        k_node      = bndelem%vertex(k_vertex)
        k_dir       = bndelem%direction(k_vertex,k_dof)
        k_size      = bndelem%size(k_vertex,k_dof)
        node_k      = node_list%node(k_node)
        if ( present(R  ) ) R  (:)  = R  (:)  + node_k%x(k_dir,1) * k_size * H1  (k_vertex,k_dof,:)
        if ( present(Z  ) ) Z  (:)  = Z  (:)  + node_k%x(k_dir,2) * k_size * H1  (k_vertex,k_dof,:)
        if ( present(R_s) ) R_s(:)  = R_s(:)  + node_k%x(k_dir,1) * k_size * H1_s(k_vertex,k_dof,:)
        if ( present(Z_s) ) Z_s(:)  = Z_s(:)  + node_k%x(k_dir,2) * k_size * H1_s(k_vertex,k_dof,:)
      end do
    end do
    
  end subroutine det_coord_bnd
  
  
  
  
  
  
  !> Determine vectors of the psi and deltapsi values at the boundary.
  subroutine det_psibnd_vec(bnd_node_list, node_list, psibnd_vec, dpsibnd_vec, psibnd_coils)

    use data_structure, only: type_node_list, type_bnd_node_list

    implicit none

    ! --- Routine parameters
    type(type_node_list),     intent(in)  :: node_list      !< List of grid nodes
    type(type_bnd_node_list), intent(in)  :: bnd_node_list  !< List of boundary grid nodes
    real*8, allocatable,      intent(out) :: psibnd_vec(:)  !< Vector of Psi boundary values
    real*8, allocatable,      intent(out) :: dpsibnd_vec(:) !< Vector of deltaPsi boundary values
    real*8, allocatable, optional, intent(out) :: psibnd_coils(:)!< Vector of Psi_coil boundary values

    ! --- Local variables
    integer :: jnode, jnode_glob, j_starwall, jtor, jbas, jdir, j_resp, j_resp_0, j_resp_old

    write(*,*) 'det_psibnd_vec: n_dof_starwall : ',n_dof_starwall

    if ( allocated(psibnd_vec) ) deallocate(psibnd_vec)
    allocate( psibnd_vec(n_dof_starwall) )
    psibnd_vec(:) = 0.d0

    if ( allocated(dpsibnd_vec) ) deallocate(dpsibnd_vec)
    allocate( dpsibnd_vec(n_dof_starwall) )
    dpsibnd_vec(:) = 0.d0

    if ( present(psibnd_coils) ) then
      if ( allocated(psibnd_coils) ) deallocate(psibnd_coils)
      allocate( psibnd_coils(n_dof_starwall) )
      psibnd_coils(:) = 0.d0
    end if

    ! --- Determine vector of (delta)psi boundary values.
    do jnode = 1, bnd_node_list%n_bnd_nodes       ! loop over nodes

      jnode_glob = bnd_node_list%bnd_node(jnode)%index_jorek

      do j_starwall = 1, sr%n_tor     ! loop over STARWALL harmonics

        jtor = sr%i_tor(j_starwall)   ! (mode corresponding to STARWALL harmonic)

        do jbas = 1, 2                ! loop over basis functions

          jdir         = bnd_node_list%bnd_node(jnode)%direction(jbas)
          j_resp_old   = response_index(jnode,j_starwall,jbas)
          j_resp_0     = response_index_eq(jnode,jbas)

          j_resp = (bnd_node_list%bnd_node(jnode)%index_starwall(1) - 1)*sr%n_tor &
                 +  bnd_node_list%bnd_node(jnode)%n_dof*(j_starwall-1) &
                 + (bnd_node_list%bnd_node(jnode)%index_starwall(jbas)-bnd_node_list%bnd_node(jnode)%index_starwall(1)) + 1

!          if (j_resp_old .ne. j_resp) write(*,'(A4i5)') 'PANIC jresp: ',j_resp_old,j_resp, &
!           bnd_node_list%bnd_node(jnode)%index_starwall(1), bnd_node_list%bnd_node(jnode)%index_starwall(jbas)

          psibnd_vec ( j_resp ) = node_list%node(jnode_glob)%values(jtor, jdir, ivar_psi)
          dpsibnd_vec( j_resp ) = node_list%node(jnode_glob)%deltas(jtor, jdir, ivar_psi)

          if ( (present(psibnd_coils)) .and. (allocated(I_coils)) .and. (jtor==1) ) then
            j_resp_0 = 2*(jnode-1) + jbas
            psibnd_coils( j_resp_0 ) = sum( bext_psi(j_resp_0,:) * I_coils(:) )
          end if

        end do
      end do
    end do
    
  end subroutine det_psibnd_vec
  
  
  
  
  
  
  !> Initialize the currents in the resistive wall and the external coil currents.
  subroutine init_wall_currents(my_id, resistive_wall)
    
    use nodes_elements, only: bnd_node_list, node_list
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: my_id
    logical, intent(in) :: resistive_wall
    
    ! --- Local variables
    real*8, allocatable :: psibnd_vec(:)    ! Vector of Psi values at the boundary
    real*8, allocatable :: dpsibnd_vec(:)   ! Vector of deltaPsi values at the boundary
    real*8, allocatable :: wall_and_coil_curr(:)
    integer             :: k
    
    call det_psibnd_vec(bnd_node_list, node_list, psibnd_vec, dpsibnd_vec)
    
    if ( resistive_wall ) then
      
      if ( .not. allocated(wall_curr) ) then
        allocate( wall_curr(n_wall_curr) )
        do k = 1, n_wall_curr
          wall_curr(k) = 0.d0! - sum( sr%a_ye(k,:) * psibnd_vec(:) )
        end do
      end if
      
      if ( .not. allocated(dwall_curr) ) then
        allocate( dwall_curr(n_wall_curr) )
        do k = 1, n_wall_curr
          dwall_curr(:) = 0.d0!- sum( sr%a_ye(k,:) * dpsibnd_vec(:) )
        end do
      end if         
      
    end if
    
    ! --- Also initialize the old_dpsibnd_vec
    if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
    allocate( old_dpsibnd_vec(n_dof_starwall) )
    old_dpsibnd_vec(:) = dpsibnd_vec(:)
    
    if ( my_id == 0 ) call write_wall_vtk(0, resistive_wall)
    deallocate( psibnd_vec, dpsibnd_vec )
    
    if ( vacuum_debug ) write(*,*) 'Wall currents initialized.'
    
  end subroutine init_wall_currents
  
  
  
  
  
  
  !> Perform the time-evolution of the wall currents (resistive wall).
  subroutine evolve_wall_currents(my_id, psibnd_vec, dpsibnd_vec)
    
    use phys_module, only: index_now, resistive_wall
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: my_id
    real*8,  intent(in) :: psibnd_vec (n_dof_starwall) !< Vector of Psi boundary values
    real*8,  intent(in) :: dpsibnd_vec(n_dof_starwall) !< Vector of deltaPsi boundary values
    
    ! --- Local variables
    integer :: k
    
    if ( vacuum_debug ) write(*,*) 'wall_curr(before)', sum(abs(wall_curr)), sum(wall_curr)
    
    if ( (.not. allocated(old_dpsibnd_vec)) .or. (size(old_dpsibnd_vec,1)/= n_dof_starwall) ) then
      if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
      allocate( old_dpsibnd_vec(n_dof_starwall) )
    end if
    
    do k = 1, n_wall_curr
      dwall_curr(k) = sum( response_m_a(k,:) * dpsibnd_vec(:) ) &
        + response_d_b(k) * wall_curr(k) &
        + response_d_c(k) * dwall_curr(k) &
        + sum( response_m_d(k,:) * old_dpsibnd_vec(:) ) !&
        !+ sum( response_m_k(k,:) * coil_voltages(:) )
    end do
    wall_curr(:) = wall_curr(:) + dwall_curr(:)
    
    if ( my_id == 0 ) then
      call write_wall_vtk(index_now, resistive_wall)
      if ( vacuum_debug .and. resistive_wall ) then
        call log_wall_curr()
        !call log_coil_curr()
      end if
    end if
    
    if ( vacuum_debug ) write(*,*) 'wall_curr(after)', sum(abs(wall_curr)), sum(wall_curr)
    
  end subroutine evolve_wall_currents
  
  
  
  
  
  
  !> Reconstruct the potential values at the wall triangle nodes.
  subroutine reconstruct_triangle_potentials(tripot_w, wall_curr)
    
    implicit none
    
    ! --- Routine parameters
    real*8, allocatable, intent(inout) :: tripot_w(:)
    real*8, allocatable, intent(in)    :: wall_curr(:)
    
    ! --- Local variables
    integer :: i, j
    
    if ( allocated(tripot_w) ) deallocate(tripot_w); allocate( tripot_w(sr%npot_w) )
    
    if ( allocated(wall_curr) ) then
      do i = 1, sr%npot_w
        j = i + sr%ncoil
        tripot_w(i) = sum(sr%s_ww(j,:) * wall_curr(:))
      end do
      tripot_w(1) = 0.d0
    else
      !###
    end if
    
  end subroutine reconstruct_triangle_potentials
  
  
  
  
  
  
  !> Write wall current potentials to logfile.
  subroutine log_wall_curr()
    
    implicit none
    
    ! --- Local variables
    real*8, allocatable :: tripot_w(:)
    
    call reconstruct_triangle_potentials(tripot_w, wall_curr)
    write(*,'(" Wall current potentials (min, max): ",ES16.8,"...",ES16.8)') minval(tripot_w), maxval(tripot_w)
    deallocate( tripot_w )
    
  end subroutine log_wall_curr
  
  
  
  
  
  
  !> Update the derived response matrices.
  !!
  !! This is necessary
  !! - right after the start or restart of the code
  !! - when wall resistivity, tstep, or some other parameters have changed.
  subroutine update_response(tstep, freeboundary_equil, resistive_wall)
    
    use phys_module, only: time_evol_theta, time_evol_zeta
    
    implicit none
    
    ! --- Routine parameters
    real*8,                      intent(in) :: tstep              !< delta t, timestep
    logical,                     intent(in) :: freeboundary_equil !< Use free boundary equilibrium?
    logical,                     intent(in) :: resistive_wall     !< Resistive or ideal wall?
    
    ! --- Local variables
    integer :: i, j, k, j2, k2
    real*8  :: a, b
    real*8, allocatable :: tmp_d_s(:)
    real*8  :: theta, zeta
    logical :: update_required
    
    ! --- Local variables to store the previous values of some parameters.
    real*8,  save :: old_thick
    real*8,  save :: old_res
    real*8,  save :: old_tstep
    real*8,  save :: old_theta
    real*8,  save :: old_zeta
    logical, save :: old_reswall
    
    theta = time_evol_theta
    zeta  = time_evol_zeta
    
    ! --- Update response matrices only, if parameter values changed or matrices not allocated
    update_required = ( old_thick   /= wall_thickness      ) &
                 .or. ( old_res     /= wall_resistivity    ) &
                 .or. ( old_tstep   /= tstep               ) &
                 .or. ( old_theta   /= theta               ) &
                 .or. ( old_zeta    /= zeta                ) &
                 .or. ( old_reswall .neqv. resistive_wall  ) &
                 .or. ( .not. allocated(response_m_a)      ) &
                 .or. ( .not. allocated(response_d_b)      ) &
                 .or. ( .not. allocated(response_d_c)      ) &
                 .or. ( .not. allocated(response_m_d)      ) &
                 .or. ( .not. allocated(response_m_e)      ) &
                 .or. ( .not. allocated(response_m_f)      ) &
                 .or. ( .not. allocated(response_m_g)      ) &
                 .or. ( .not. allocated(response_m_h)      ) &
                 .or. ( .not. allocated(response_m_j)      ) &
                 .or. ( .not. allocated(response_m_k)      ) &
                 .or. ( .not. allocated(response_m_l)      ) &
                 .or. ( .not. allocated(response_m_eq)     )
    
    if ( update_required ) then
      ! --- Remember parameter values.
      old_thick   = wall_thickness
      old_res     = wall_resistivity
      old_tstep   = tstep
      old_theta   = theta
      old_zeta    = zeta
      old_reswall = resistive_wall
      
      ! --- Allocate matrices if required
      if ( .not. allocated(response_m_eq) ) &
        allocate( response_m_eq(sr%nd_bez/sr%n_tor, sr%nd_bez/sr%n_tor) )
      if ( .not. allocated(response_m_a) ) &
        allocate( response_m_a(n_wall_curr, n_dof_starwall) )
      if ( .not. allocated(response_d_b) ) &
        allocate( response_d_b(n_wall_curr) )
      if ( .not. allocated(response_d_c) ) &
        allocate( response_d_c(n_wall_curr) )
      if ( .not. allocated(response_m_d) ) &
        allocate( response_m_d(n_wall_curr, n_dof_starwall) )
      if ( .not. allocated(response_m_e) ) &
        allocate( response_m_e(n_dof_starwall, n_dof_starwall) )
      if ( .not. allocated(response_m_f) ) &
        allocate( response_m_f(n_dof_starwall, n_wall_curr) )
      if ( .not. allocated(response_m_g) ) &
        allocate( response_m_g(n_dof_starwall, n_wall_curr) )
      if ( .not. allocated(response_m_h) ) &
        allocate( response_m_h(n_dof_starwall, n_dof_starwall) )
      if ( .not. allocated(response_m_j) ) &
        allocate( response_m_j(n_dof_starwall, n_dof_starwall) )
      if ( .not. allocated(response_m_k) ) &
        allocate( response_m_k(n_wall_curr, sr%ncoil) )
      if ( .not. allocated(response_m_l) ) &
        allocate( response_m_l(n_dof_starwall, sr%ncoil) )
      
      ! --- Derived response matrix for equilibrium (extract n=0 part from STARWALL EE matrix)
      
      response_m_eq = 0.d0

      do j = 1, sr%nd_bez, 2
        j2 = (j-1)*sr%n_tor+1
        do k = 1, sr%nd_bez, 2
          k2 = (k-1)*sr%n_tor+1
          response_m_eq(j:j+1,k:k+1) = sr%a_ee(j2:j2+1,k2:k2+1)
        end do
      end do
      
      ! --- Derived response matrices for time-evolution
      if ( resistive_wall ) then
        
        allocate( tmp_d_s(n_wall_curr) )
        
        tmp_d_s(:) = 1.d0 + zeta + tstep * theta * wall_resistivity / wall_thickness * sr%d_yy(:)
        
        do j = 1, n_dof_starwall
          response_m_a(:,j) = -(1.d0+zeta) * sr%a_ye(:,j) / tmp_d_s(:)
        end do
        
        response_d_b(:) = - tstep * wall_resistivity / wall_thickness * sr%d_yy(:) / tmp_d_s(:)
        
        response_d_c(:) = zeta / tmp_d_s(:)
        
        do j = 1, n_dof_starwall
          response_m_d(:,j) = zeta * sr%a_ye(:,j) / tmp_d_s(:)
        end do
        
        response_m_e(:,:) = sr%a_ee(:,:) + matmul( sr%a_ey(:,:), response_m_a(:,:) )
        
        do k = 1, n_wall_curr
          response_m_f(:,k) = sr%a_ey(:,k) * ( 1.d0 + response_d_b(k) )
        end do
        
        do k = 1, n_wall_curr
          response_m_g(:,k) = sr%a_ey(:,k) * response_d_c(k)
        end do
        
        response_m_h(:,:) = sr%a_ee(:,:)
        
        response_m_j(:,:) = matmul( sr%a_ey(:,:), response_m_d(:,:) )
        
        do k = 1, n_wall_curr
          response_m_k(k,:) = -tstep * sr%d_yy(k) * sr%s_ww(:,k)
        end do
        
        response_m_l(:,:) = matmul( sr%a_ey(:,:), response_m_k(:,:) )
        
        deallocate( tmp_d_s )
        
      else ! (Ideal wall)
        
        response_m_a(:,:) = 0.d0
        response_d_b(:)   = 0.d0
        response_d_c(:)   = 0.d0
        response_m_d(:,:) = 0.d0
        response_m_e(:,:) = sr%a_id(:,:)
        response_m_f(:,:) = 0.d0
        response_m_g(:,:) = 0.d0
        response_m_h(:,:) = sr%a_id(:,:)
        response_m_j(:,:) = 0.d0
        response_m_k(:,:) = 0.d0 !####
        response_m_l(:,:) = 0.d0 !####
        
      end if
      
      if ( vacuum_debug ) then
        write(*,*) 'DEBUG: Checksums'
        write(*,*) 'm_a:', sum(abs(response_m_a)), sum(response_m_a)
        write(*,*) 'd_b:', sum(abs(response_d_b)), sum(response_d_b)
        write(*,*) '1+d_b:', sum(abs(1.d0+response_d_b)), sum(1.d0+response_d_b)
        write(*,*) 'd_c:', sum(abs(response_d_c)), sum(response_d_c)
        write(*,*) 'm_d:', sum(abs(response_m_d)), sum(response_m_d)
        write(*,*) 'm_e:', sum(abs(response_m_e)), sum(response_m_e)
        write(*,*) 'm_f:', sum(abs(response_m_f)), sum(response_m_f)
        write(*,*) 'm_g:', sum(abs(response_m_g)), sum(response_m_g)
        write(*,*) 'm_h:', sum(abs(response_m_h)), sum(response_m_h)
        write(*,*) 'm_j:', sum(abs(response_m_j)), sum(response_m_j)
        write(*,*) 'm_k:', sum(abs(response_m_k)), sum(response_m_k)
        write(*,*) 'm_l:', sum(abs(response_m_l)), sum(response_m_l)
        write(*,*) 'm_eq:', sum(abs(response_m_eq(1:128,1:128))), sum(response_m_eq(1:128,1:128))
        write(*,*) 'END: Checksums'
      end if
      
    end if
    
  end subroutine update_response
  
  
  
  
  
  
  !> Read a STARWALL response matrix from a file.
  subroutine read_response_matrix( matrix, dim_expected, filename, err )
    
    ! --- Routine parameters
    real*8, allocatable, intent(inout) :: matrix(:,:)     !< Matrix to be read
    integer,             intent(in)    :: dim_expected(2) !< Matrix dimension expected
    character(len=*),    intent(in)    :: filename        !< Filename to read from
    integer,             intent(out)   :: err             !< Error code
    
    ! --- Local variables
    integer :: dim(2), i, j, i2, j2
    
    err = 0
    
    open(42, FILE=trim(filename), status='old', action='read', iostat=err)
    if ( err /= 0 ) return

    read(42,*) dim

    if ( (dim(1) /= dim_expected(1)) .or. (dim(2) /= dim_expected(2)) ) then
      write(*,*) 'FATAL ERROR: Matrix dimension not as expected. Different resolutions?'
      write(*,'(1x,A,2I7)') 'dim_expected=',dim_expected
      stop
    end if
    
    if ( allocated(matrix) ) deallocate(matrix)
    allocate( matrix(dim(1),dim(2)) )
    matrix = 0.d0
    
    do i = 1, dim(1)
      do j = 1, dim(2)
        
        read (42,*) i2, j2, matrix(i,j)
        
        if ( ( i2 /= i ) .or. ( j2 /= j ) ) then
          write(*,*) 'FATAL ERROR: Matrix indices not as expected. Different resolutions?'
          stop
        end if
      
      end do
    end do
      
    close(42)

  end subroutine read_response_matrix
  
  
  
  
  
  
  !> Read a diagonal STARWALL response matrix from a file.
  subroutine read_response_diagonal( diagonal, dim_expected, filename, err )
    
    ! --- Routine parameters
    real*8, allocatable, intent(inout) :: diagonal(:)     !< Matrix to be read
    integer,             intent(in)    :: dim_expected    !< Matrix dimension expected
    character(len=*),    intent(in)    :: filename        !< Filename to read from
    integer,             intent(out)   :: err             !< Error code
    
    ! --- Local variables
    integer :: dim, i, i2
    
    open(42, FILE=trim(filename), status='old', action='read', iostat=err)
    if ( err /= 0 ) return

    read(42,*) dim

    if ( dim /= dim_expected ) then
      write(*,*) 'FATAL ERROR: Matrix dimension not as expected. Different resolutions?'
      stop
    end if
    
    if ( allocated(diagonal) ) deallocate(diagonal)
    allocate( diagonal(dim) )
    diagonal = 0.d0
    
    do i = 1, dim
        
      read (42,*) i2, diagonal(i)
      
      if ( i2 /= i ) then
        write(*,*) 'FATAL ERROR: Matrix indices not as expected. Different resolutions?'
        stop
      end if
        
    end do
      
    close(42)
    
  end subroutine read_response_diagonal
  
  
  
  
  
  
  !> Determine the index in the response matrix for a certain boundary degree of freedom.
  integer recursive function response_index(inode, i_starwall, ibas)
    
    ! --- Routine parameters
    integer, intent(in)    :: inode      !< Boundary index of the node
    integer, intent(in)    :: i_starwall !< STARWALL harmonic
    integer, intent(in)    :: ibas       !< Basis function (1 or 2)
    
    if ( (i_starwall < 0) .or. (i_starwall > sr%n_tor) ) then
      write(*,*) 'response_index: illegal value i_starwall=', i_starwall
      stop
    end if
    response_index = 2*sr%n_tor*(inode-1) + 2*(i_starwall-1) + ibas
    
    if ( response_index < 1 ) then
      write(*,*) 'FATAL: RESPONSE_INDEX < 1 DETECTED'
      stop
    end if
    
  end function response_index
  
  
  
  
  
  
  !> Determine the index in the response matrix for a certain boundary degree of freedom.
  integer recursive function response_index_eq(inode, ibas)
    
    ! --- Routine parameters
    integer, intent(in)    :: inode      !< Boundary index of the node
    integer, intent(in)    :: ibas       !< Basis function (1 or 2)
    
    response_index_eq = 2*(inode-1) + ibas
    
    if ( response_index_eq < 1 ) then
      write(*,*) 'FATAL: RESPONSE_INDEX_EQ < 1 DETECTED'
      stop
    end if
    
  end function response_index_eq
  
  
  
  
  
  
  !> Return a description for a given toroidal mode index (mode number, cos/sin).
  character(len=12) function mode_to_str(i_tor, n_period)
    
    use phys_module, only: mode, mode_type
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: i_tor     !< Toroidal mode index
    integer, intent(in) :: n_period  !< Periodicity
    
    ! --- Local variables
    integer           :: n                ! Toroidal mode number
    character(len=3)  :: typ              ! sin or cos
    character(len=30) :: i_tor_str, n_str ! Character string representations for i_tor and n
    

    write(i_tor_str,'(I10)') i_tor   ! toroidal mode index
    n   = mode(i_tor)                ! toroidal mode number
    write(n_str,'(I10)') n           !  -"-
    typ = mode_type(i_tor)           ! sin or cos
    
    mode_to_str = trim(adjustl(i_tor_str))//' (n='//trim(adjustl(n_str))//' '//trim(typ)//')'
    
  end function mode_to_str
  
  
  
  
  
  
  !> Return a description of several toroidal modes.
  character(len=1400) function modes_to_str(i_tors, n_tor, n_period)
  
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: n_tor         !< Dimension of i_tors
    integer, intent(in) :: i_tors(n_tor) !< Toroidal mode numbers
    integer, intent(in) :: n_period      !< Periodicity
    
    ! --- Local variables
    integer :: i
    
    modes_to_str = ''
    do i = 1, n_tor
      modes_to_str = trim(modes_to_str)//' '//mode_to_str(i_tors(i), n_period)
      if ( i == n_tor ) exit
      modes_to_str = trim(modes_to_str)//', '
    end do
    
    modes_to_str = adjustl(modes_to_str)
    
  end function modes_to_str
  
  
  
  
  
  
end module vacuum_response
