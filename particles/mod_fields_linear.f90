!> Module for linearly interpolating in time between values and deltas
!> in JOREK restart files. Contains an action to read the fields.
module mod_fields_linear
use data_structure
use mod_particle_sim
use mod_event
use mod_fields
use mod_interp4
use mod_interp_PRZ
implicit none
private
public jorek_fields_interp_linear, read_jorek_fields_interp_linear

!> Action to read in the fields into sim%fields
type, extends(action) :: read_jorek_fields_interp_linear
  character(len=80) :: basename = 'jorek' !< Comes before the file number or extension
  integer :: i = 0 !< Number of the restart file to read. Set to -1 to not include. Corresponds to the index of NOW
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
  real*8 :: time_now  = 0.d0 !< Time of current restart file (SI units)
  real*8 :: time_prev = 0.d0!< Time of previous restart file (SI units)
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

  real*8 :: f, df, dt !< result = (1-df)*values_now - df*deltas, df = 1-f
  real*8, dimension(n_v) :: Pd, Pd_s, Pd_t, Pd_phi

  call       interp_PRZ(this%node_list,this%element_list,i_elm,i_v,n_v,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
  if (abs(this%time_prev-this%time_now) .gt. 1d-10) then
    dt = 1.d0/(this%time_now - this%time_prev)
    df = (this%time_now - time)*dt
    call interp_PRZ_delta(this%node_list,this%element_list,i_elm,i_v,n_v,s,t,phi,Pd,Pd_s,Pd_t,Pd_phi,R,R_s,R_t,Z,Z_s,Z_t)
    P      = P     - df*Pd
    P_s    = P_s   - df*Pd_s
    P_t    = P_t   - df*Pd_t
    P_phi  = P_phi - df*Pd_phi
    P_time = Pd*dt
  else
    P_time = 0.d0
  end if

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
subroutine do_read(this, sim, ev)
  use mod_import_restart
  use phys_module
  class(read_jorek_fields_interp_linear), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), optional :: ev
  character(len=80) :: restart_file
  integer :: i, ierr
  logical :: file_exists

  real*8 :: t_norm
  t_norm = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

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

    ! Read only one file
    if (this%i .eq. -1) then
      write(restart_file,'(A,A)') trim(this%basename), '_restart.h5'
      inquire(file=trim(restart_file), exist=file_exists)
      if (file_exists) then
        call import_hdf5_restart(f%node_list,f%element_list,restart_file,this%rst_format,ierr)
        f%time_now  = 0.d0
        f%time_prev = 0.d0
      else
        write(*,*) "ERROR: file ", trim(restart_file), " does not exist"
        call exit(1)
      end if
    else ! Linearly interpolating case
      ! If nothing has been loaded (i.e. fields%time_prev = 0.d0) load the initial file
      if (abs(f%time_prev) .lt. 1.d-50) then
        write(restart_file,'(A,i5.5,A)') trim(this%basename), this%i, '.h5'
        inquire(file=trim(restart_file), exist=file_exists)
        if (file_exists) then
          call import_hdf5_restart(f%node_list,f%element_list,trim(restart_file),this%rst_format,ierr)
          if (ierr .ne. 0) then
            write(*,*) "ERROR: cannot open restart file"
            call exit(1)
          else
            f%time_prev = t_start*t_norm ! set by import_hdf5_restart
            ! Set sim%time to this time also, to start at the right point
            sim%time = f%time_prev
            write(*,"(A,f9.8,A)") "Read initial restart file, values at t=", f%time_prev, " [s]"
          end if
        else
          write(*,*) "ERROR: cannot read initial file ", trim(restart_file)
          call exit(1)
        end if
      end if
      
      ! Find the following file (next timestep number)
      do i=this%i+1,this%i+10
        write(restart_file,'(A,i5.5,A)') trim(this%basename), i, '.h5'
        inquire(file=trim(restart_file), exist=file_exists)
        if (file_exists) then
          call merge_restart(f%node_list, f%element_list, trim(restart_file), this%rst_format, ierr)
          if (ierr .ne. 0) then
            write(*,*) "ERROR: cannot open restart file"
            call exit(1)
          else
            f%time_now = t_start*t_norm ! Set by import_merge_restart
            this%i=i
            ! set the time to run this event at next
            write(*,"(A,f9.8,A)") " Read next restart file, values until t=", f%time_now, " [s]"
            if (present(ev)) then
              ev%start = f%time_now
              write(*,*) "Set time for next restart file read to ", ev%start
            end if
            exit ! the file-finding loop
          end if
        endif
      enddo
      if (i .gt. this%i+10) then
        write(*,*) "ERROR: cannot find any next restart files"
        call exit(1)
      end if
    end if

    call update_neighbours(f%element_list, f%node_list)
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
  !$omp parallel do default(shared) private(inode)
  do inode=1,node_list%n_nodes
    values(:,:,:,inode) = node_list%node(inode)%values(:,:,:)
  enddo
  !$omp end parallel do
  tstart_old = t_start

  ! Import new values
  call import_hdf5_restart(node_list,element_list, restart_file, format_rst, ierr)

  ! Calculate deltas as values_new - values_old
  !$omp parallel do default(shared) private(inode)
  do inode=1,node_list%n_nodes
    node_list%node(inode)%deltas = node_list%node(inode)%values - values(:,:,:,inode)
  enddo
  !$omp end parallel do

  ! Set timestep to time between restart files
  tstep = t_start - tstart_old

  deallocate(values)
end subroutine merge_restart
end module mod_fields_linear
