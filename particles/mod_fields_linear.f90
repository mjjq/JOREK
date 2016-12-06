!> Module for linearly interpolating in time between values and deltas
!> in JOREK restart files. Contains an action to read the fields.
module mod_fields_linear
use data_structure
use mod_particle_sim
use mod_fields
use mod_interp4
use mod_action
use mod_interp_PRZ
implicit none
private
public jorek_fields_interp_linear, read_jorek_fields_interp_linear

!> Action to read in the fields into sim%fields
type, extends(action) :: read_jorek_fields_interp_linear
  character(len=80) :: basename = 'jorek' !< Comes before the file number or extension
  integer :: i = 0 !< Number of the restart file to read. Set to -1 to not include
  integer :: rst_format = 0 !< Format of restart file if .rst type
  contains
    procedure :: do => do_read
end type
interface read_jorek_fields_interp_linear
  module procedure new_read_jorek_fields_interp_linear
end interface read_jorek_fields_interp_linear

!> store enough data for linear interpolation
!> in time.
!> This does not really work for changing grids in time!
!> Use at your own peril in that case. (e.g. i_elm and s,t for a specific spatial position might depend on time)
!>
!> The reason behind not using deltas is that we do not have to alter much code
!> and can import two restarts which are not consecutive and still interpolate.
type, extends(fields_base) :: jorek_fields_interp_linear
  type(type_node_list), allocatable    :: node_list
  type(type_element_list), allocatable :: element_list
  real*8 :: time_now !< Time of current restart file (SI)
  real*8 :: time_prev !< Time of previous restart file (SI)
  contains
    procedure :: interp_PRZ => do_interp_PRZ
end type jorek_fields_interp_linear
contains

!> Interpolate a variable at a specific position (with phi), with first derivatives only
pure subroutine do_interp_PRZ(this, time, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
  use mod_interp_PRZ
  class(jorek_fields_interp_linear),  intent(in)  :: this
  real*8,                   intent(in)  :: time !< Time at which to calculate this variable
  integer,                  intent(in)  :: i_elm
  integer,                  intent(in)  :: n_v, i_v(n_v)
  real*8,                   intent(in)  :: s, t, phi
  real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v), P_phi(n_v), P_time(n_v)
  real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t

  real*8 :: f, df !< result = (1-df)*values_now - df*deltas, df = 1-f
  real*8, dimension(n_v) :: Pd, Pd_s, Pd_t, Pd_phi
  df = (time - this%time_prev)/(this%time_now - this%time_prev)
  f  = 1.d0 - df

  call       interp_PRZ(this%node_list,this%element_list,i_elm,i_v,n_v,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
  call interp_PRZ_delta(this%node_list,this%element_list,i_elm,i_v,n_v,s,t,phi,Pd,Pd_s,Pd_t,Pd_phi,R,R_s,R_t,Z,Z_s,Z_t)

  P      = f*P     - df*Pd
  P_s    = f*P_s   - df*Pd_s
  P_t    = f*P_t   - df*Pd_t
  P_phi  = f*P_phi - df*Pd_phi
  P_time = Pd/(this%time_now - this%time_prev) ! linearisation
end subroutine do_interp_PRZ

!> Constructor to allow for optional and default variables
function new_read_jorek_fields_interp_linear(basename, i, rst_format) result(new)
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional :: i
  integer, intent(in), optional :: rst_format
  type(read_jorek_fields_interp_linear) :: new
  if (present(basename)) new%basename = basename
  if (present(i)) new%i = i
  if (present(rst_format)) new%rst_format = rst_format
  new%name = "ReadJorekFieldsInterpLinear"
  new%log = .true.
end function new_read_jorek_fields_interp_linear



!> Read jorek fields from the next restart file (regardless of when the file number)
subroutine do_read(this, sim)
  use mod_import_restart
  class(read_jorek_fields_interp_linear), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  character(len=80) :: restart_file
  integer :: i, ierr
  logical :: file_exists

  logical, save :: neighbours_updated = .false.

  ! Check that the right fields are allocated in sim and allocate if needed
  if (allocated(sim%fields)) then
    select type (f => sim%fields)
    type is (jorek_fields_interp_linear) ! do nothing
    class default
      write(*,*) "WARNING: wrong type of fields in particle%sim, reallocating"
      deallocate(sim%fields)
      allocate(jorek_fields_interp_linear::sim%fields)
    end select
  else
    allocate(jorek_fields_interp_linear::sim%fields)
  end if

  ! Continue for jorek_fields_interp_linear
  select type (f => sim%fields)
  type is (jorek_fields_interp_linear)
    ! Allocate node and element_list if needed
    if (.not. allocated(f%node_list))    allocate(f%node_list)
    if (.not. allocated(f%element_list)) allocate(f%element_list)

    ! Read only one file
    if (this%i .eq. -1) then
      write(restart_file,'(A,A)') trim(this%basename), '.h5'
      inquire(file=trim(restart_file), exist=file_exists)
      if (file_exists) then
        call import_hdf5_restart(f%node_list,f%element_list,restart_file,this%rst_format,ierr)
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
        call import_hdf5_restart(f%node_list,f%element_list,trim(restart_file),this%rst_format,ierr)
        this%i = this%i+1
      else ! loop over the next few files of this name format
        do i=this%i+2,this%i+10
          write(restart_file,'(A,i5.5,A)') trim(this%basename), i, '.h5'
          inquire(file=trim(restart_file), exist=file_exists)
          if (file_exists) then
            call import_hdf5_restart(f%node_list,f%element_list,trim(restart_file),this%rst_format,ierr)
            this%i=i
            exit
          endif
        enddo
      end if
    end if

    ! After reading we need to update_neighbours, but only do it the first time
    ! This will not work for simulations with refinement!
    if (.not. neighbours_updated) then
      call update_neighbours(f%element_list, f%node_list)
      neighbours_updated = .true.
    end if
  class default
    write(*,*) "ERROR, do_read called with wrong sim%fields"
  end select
end subroutine do_read



!> Import a binary restart file and merges it with the values currently known
!> This can then be used to interpolate linearly between any two restart files
subroutine merge_restart(node_list,element_list, restart_file, format_rst, ierr)
  use data_structure
  use phys_module
  use mod_import_restart
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
end subroutine merge_restart
end module mod_fields_linear
