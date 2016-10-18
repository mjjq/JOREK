!> Module for linearly interpolating in time between values and deltas
!> in JOREK restart files. Contains an action to read the fields.
module mod_jorek_fields_interp_linear
use data_structure
use mod_action
use mod_particle_sim
implicit none
private
public read_jorek_fields_interp_linear, EM_fields_interp_linear

!> Action to read jorek fields and store enough data for linear interpolation
!> in time.
type, extends(action) :: read_jorek_fields_interp_linear
  type(type_node_list), allocatable    :: node_list
  type(type_element_list), allocatable :: element_list
  character(len=80) :: basename = 'jorek' !< Comes before the file number or extension
  integer :: i = 0 !< Number of the restart file to read. Set to -1 to not include
  integer :: rst_format = 0 !< Format of restart file if .rst type
  contains
    procedure :: do => read_jorek_fields_impl
end type read_jorek_fields_interp_linear
interface read_jorek_fields_interp_linear
  module procedure new_read_jorek_fields_interp_linear
end interface read_jorek_fields_interp_linear
contains

!> Constructor to allow for optional and default variables
function new_read_jorek_fields_interp_linear(basename, i, rst_format) result(new)
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional :: i
  integer, intent(in), optional :: rst_format
  type(read_jorek_fields_interp_linear) :: new
  allocate(new%node_list)
  allocate(new%element_list)
  if (present(basename)) new%basename = basename
  if (present(i)) new%i = i
  if (present(rst_format)) new%rst_format = rst_format
  new%name = "ReadJorekFieldsInterpLinear"
  new%log = .true.
end function new_read_jorek_fields_interp_linear



!> Read jorek fields from a restart file
subroutine read_jorek_fields_impl(this, sim)
  use import_restart
  class(read_jorek_fields_interp_linear), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  character(len=80) :: restart_file
  integer :: i, ierr
  logical :: file_exists
  integer :: DUMMY_INT
  real*8  :: DUMMY_REAL

  logical, save :: neighbours_updated = .false.

  ! Read only one file
  if (this%i .eq. -1) then
    write(restart_file,'(A,A)') trim(this%basename), '.h5'
    inquire(file=trim(restart_file), exist=file_exists)
    if (file_exists) then
      call import_hdf5_restart(this%node_list,this%element_list,restart_file,this%rst_format,ierr)
    else
      write(*,*) "ERROR: file ", trim(restart_file), " does not exist"
      call exit(1)
    end if
  else ! Linearly interpolating case
    write(*,*) "ERROR: reading with interpolation not implemented yet"
    call exit(1)
    ! TODO: set the time to run the event at next
    ! TODO: recalculate simulation time and timestep
    ! If not, keep looping (up to 10) to find one, and use the merge import
    ! This assumes that the current node_list contains the values
    ! at time istep (but does not need to contain the deltas, these are calculated)
    write(restart_file,'(A,i5.5,A)') trim(this%basename), this%i+1, '.h5'
    inquire(file=trim(restart_file), exist=file_exists)
    if (file_exists) then
      ! If so, import it and we're done
      call import_hdf5_restart(this%node_list,this%element_list,trim(restart_file),this%rst_format,ierr)
      this%i = this%i+1
    else ! loop over the next few files of this name format
      do i=this%i+2,this%i+10
        write(restart_file,'(A,i5.5,A)') trim(this%basename), i, '.h5'
        inquire(file=trim(restart_file), exist=file_exists)
        if (file_exists) then
          call import_hdf5_restart(this%node_list,this%element_list,trim(restart_file),this%rst_format,ierr)
          this%i=i
          exit
        endif
      enddo
    end if
  end if

  ! After reading we need to update_neighbours, but only do it the first time
  ! This will not work for simulations with refinement!
  if (.not. neighbours_updated) then
    call update_neighbours(this%element_list, this%node_list)
    neighbours_updated = .true.
  end if

  ! Call find_RZ once to initialise elements_minmax
  call find_RZ(this%node_list,this%element_list,0.d0,0.d0,DUMMY_REAL,DUMMY_REAL,DUMMY_INT,DUMMY_REAL,DUMMY_REAL,DUMMY_INT)
end subroutine read_jorek_fields_impl

!> Calculates the electric and magnetic fields at a specific position
!> in the jorek element `i_elm` at `st`.
!> Linear interpolation with element%deltas is performed according to
!> `delta_fraction`, which starts at 1 and goes to 0 for no mixing.
!> If it is 1 we get the fields of the previous timesteps.
pure subroutine EM_fields_interp_linear(fields, i_elm, st, phi, E, B, psi, U, delta_fraction)
use data_structure
use parameters
use constants
use phys_module, only : F0, tstep

interface
  pure subroutine interp_PRZ(node_list, element_list, i_elm, i_v, n_v, &
          s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
    import :: type_node_list, type_element_list
    type (type_node_list),    intent(in)  :: node_list
    type (type_element_list), intent(in)  :: element_list
    integer,                  intent(in)  :: i_elm
    integer,                  intent(in)  :: n_v, i_v(n_v)
    real*8,                   intent(in)  :: s, t, phi
    real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v)
    real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
    real*8,                   intent(out) :: P_phi(n_v)
  end subroutine interp_PRZ
  pure subroutine interp_PRZ_delta(node_list, element_list, i_elm, i_v, n_v, &
          s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
    import :: type_node_list, type_element_list
    type (type_node_list),    intent(in)  :: node_list
    type (type_element_list), intent(in)  :: element_list
    integer,                  intent(in)  :: i_elm
    integer,                  intent(in)  :: n_v, i_v(n_v)
    real*8,                   intent(in)  :: s, t, phi
    real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v)
    real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
    real*8,                   intent(out) :: P_phi(n_v)
  end subroutine interp_PRZ_delta
end interface

! Routine parameters
type (read_jorek_fields_interp_linear), intent(in) :: fields
integer, intent(in) :: i_elm !< JOREK element index
real*8, intent(in)  :: st(2) !< element-local coordinates
real*8, intent(in)  :: phi !< toroidal angle
real*8, intent(in), optional :: delta_fraction !< linear interpolation factor. goes from 1 to 0 in time

real*8, intent(out) :: E(3), B(3), psi, U !< Fields and potentials

! Internal parameters
integer :: i_var(2)
real*8                    :: P(2), P_s(2), P_t(2), P_phi(2) ! Placeholder for evaluating variables and derivatives locally
! Values
real*8                    :: R, R_s, R_t, Z, Z_s, Z_t
! Deltas
real*8                    :: Pd(2), Pd_s(2), Pd_t(2), Pd_phi(2) ! Placeholder for evaluating variables and derivatives locally
! Others
real*8                    :: inv_st_jac, R_inv
real*8                    :: psi_R, psi_Z, U_R, U_Z, U_phi
real*8                    :: psi_time !< Derivative of psi to time (linearly interpolated from values and deltas)

! Select psi and U
i_var = (/1,2/)

! Interpolate the fields to get psi and U at the current position (and the
! changes u_n - u(n-1))
call       interp_PRZ(fields%node_list,fields%element_list,i_elm,i_var,2,st(1),st(2),phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
call interp_PRZ_delta(fields%node_list,fields%element_list,i_elm,i_var,2,st(1),st(2),phi,Pd,Pd_s,Pd_t,Pd_phi,R,R_s,R_t,Z,Z_s,Z_t)

R_inv = 1.d0/R
inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)

if (present(delta_fraction)) then
  psi_R    = -(  Pd_s(1) * Z_t - Pd_t(1) * Z_s ) * inv_st_jac * delta_fraction
  psi_Z    = -(- Pd_s(1) * R_t + Pd_t(1) * R_s ) * inv_st_jac * delta_fraction
  U_R      = -(  Pd_s(2) * Z_t - Pd_t(2) * Z_s ) * inv_st_jac * delta_fraction
  U_Z      = -(- Pd_s(2) * R_t + Pd_t(2) * R_s ) * inv_st_jac * delta_fraction
  U_phi    = -Pd_phi(2) * delta_fraction
  psi = -Pd(1) * delta_fraction
  U   = -Pd(2) * delta_fraction
else
  psi_R = 0.d0
  psi_Z = 0.d0
  U_R   = 0.d0
  U_Z   = 0.d0
  U_phi = 0.d0
  psi = 0.d0
  U   = 0.d0
endif

! Calculate the derivatives to R and Z (at delta_fraction if presen)
psi_R    = psi_R + (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
psi_Z    = psi_Z + (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
U_R      = U_R   + (  P_s(2) * Z_t - P_t(2) * Z_s ) * inv_st_jac
U_Z      = U_Z   + (- P_s(2) * R_t + P_t(2) * R_s ) * inv_st_jac

! Update psi and U
psi = psi + P(1)
U   = U   + P(2)

! Use the current timestep size (very rough)
psi_time = Pd(1)/tstep

! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
B     = (/ + psi_Z, - psi_R, F0 /) * R_inv
! The local electric field, obtained from E=-Grad (u F0)-\partial_t A
! See http://jorek.eu/wiki/doku.php?id=u_phi
E     = (/ - F0 * U_R, - F0 * U_Z, - F0 * U_phi * R_inv - R * psi_time /)
end subroutine EM_fields_interp_linear



!> Import a binary restart file and merges it with the values currently known
!> This can then be used to interpolate linearly between any two restart files
subroutine import_merge_restart(node_list,element_list, restart_file, format_rst, ierr)
  use data_structure
  use phys_module
  use import_restart
  implicit none

  ! --- Routine parameters
  type(type_node_list),    intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  character(len=*),        intent(in)    :: restart_file !< Filename of new restart file to import
  integer,                 intent(out)   :: ierr
  integer,                 intent(in)    :: format_rst !< Restart file format

  ! --- Internal variables
  real*8, allocatable, dimension(:,:,:,:) :: values
  integer :: inode
  real*8 :: tstart_old

  ! Save the old values to calculate the new deltas
  allocate(values(n_tor,n_order+1,n_var,node_list%n_nodes))
  do inode=1,node_list%n_nodes
    values(:,:,:,inode) = node_list%node(inode)%values(:,:,:)
  enddo
  tstart_old = t_start

  ! Import new values
  call import_binary_restart(node_list,element_list, restart_file, format_rst, ierr)

  ! Calculate deltas as values_new - values_old
  do inode=1,node_list%n_nodes
    node_list%node(inode)%deltas = node_list%node(inode)%values - values(:,:,:,inode)
  enddo

  ! Set timestep to time between restart files
  tstep = t_start - tstart_old

  deallocate(values)
end subroutine import_merge_restart
end module mod_jorek_fields_interp_linear
