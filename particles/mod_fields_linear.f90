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
public jorek_fields_interp_linear, read_jorek_fields_interp_linear, last_file_before_time

!> Action to read in the fields into sim%fields
type, extends(action) :: read_jorek_fields_interp_linear
  character(len=80) :: basename = 'jorek' !< Comes before the file number or extension
  integer :: i = 0 !< Number of the restart file to read. Set to -1 to not include. Corresponds to the index of NOW
  integer :: rst_format = 0 !< Format of restart file if .rst type
  logical :: stop_at_end = .true. !< Whether to stop the simulation at the end of the file list
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
  logical :: static = .false. !< If true do not do time interpolation
  contains
    procedure :: interp_PRZ => do_interp_PRZ
end type jorek_fields_interp_linear
contains

!> Interpolate a variable at a specific position (with phi), with first derivatives only
pure subroutine do_interp_PRZ(this, time, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
  use mod_interp_PRZ
  use constants, only: mu_zero, mass_proton
  use phys_module, only: tstep, central_mass, central_density
  class(jorek_fields_interp_linear),  intent(in)  :: this
  real*8,                   intent(in)  :: time !< Time at which to calculate this variable
  integer,                  intent(in)  :: i_elm
  integer,                  intent(in)  :: n_v, i_v(n_v)
  real*8,                   intent(in)  :: s, t, phi
  real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v), P_phi(n_v), P_time(n_v)
  real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t

  real*8 :: f, df, dt !< result = (1-df)*values_now - df*deltas, df = 1-f
  real*8, dimension(n_v) :: Pd, Pd_s, Pd_t, Pd_phi
  real*8 :: t_norm
  t_norm = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

  call       interp_PRZ(this%node_list,this%element_list,i_elm,i_v,n_v,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
  if (abs(this%time_prev-this%time_now) .gt. 1d-10 .and. .not. this%static) then
    dt = 1.d0/(this%time_now - this%time_prev)
    df = (this%time_now - time)*dt
    call interp_PRZ_delta(this%node_list,this%element_list,i_elm,i_v,n_v,s,t,phi,Pd,Pd_s,Pd_t,Pd_phi,R,R_s,R_t,Z,Z_s,Z_t)
    P      = P     - df*Pd
    P_s    = P_s   - df*Pd_s
    P_t    = P_t   - df*Pd_t
    P_phi  = P_phi - df*Pd_phi
    P_time = Pd*dt
  else
    call interp_PRZ_delta(this%node_list,this%element_list,i_elm,i_v,n_v,s,t,phi,Pd,Pd_s,Pd_t,Pd_phi,R,R_s,R_t,Z,Z_s,Z_t)
    P_time = Pd/(tstep*t_norm)
  end if

end subroutine do_interp_PRZ


!> Constructor to allow for optional and default variables
function new_read_jorek_fields_interp_linear(basename, i, rst_format, stop_at_end) result(new)
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional :: i
  integer, intent(in), optional :: rst_format
  logical, intent(in), optional :: stop_at_end
  type(read_jorek_fields_interp_linear) :: new
  if (present(basename)) new%basename = basename
  if (present(i)) new%i = i
  if (present(rst_format)) new%rst_format = rst_format
  if (present(stop_at_end)) new%stop_at_end = stop_at_end
  new%name = "ReadJorekFieldsInterpLinear"
  new%log = .true.
end function new_read_jorek_fields_interp_linear


!> Find the number of the latest restart file < time (SI units)
!> Perform a bisection method of all the jorek$num.h5 files in the directory
!> Perhaps better to use xtime if this is always present, combined with a filter
!> for all step numbers that are in the current directory
function last_file_before_time(time) result(file_number)
  use phys_module
  use mpi
  real*8, intent(in) :: time
  integer :: file_number
  integer :: i, ierr, my_id, u
  character(len=5) :: my_id_s, num_s
  integer, dimension(:), allocatable :: filenums_tmp, filenums
  integer :: n, i_lower, i_guess, i_upper, io
  real*8 :: t_norm, t_lower, t_guess, t_upper

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  t_norm = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds
  if (my_id .eq. 0) write(*,*) "Looking for jorek restart file just before time ", time

  ! Get list of filenumbers
  write(my_id_s,"(i0.5)") my_id
  call execute_command_line("ls jorek[0-9]*.h5 | grep -o '[0-9]\{5\}' > .jorek_filenums."//my_id_s)
  open(newunit=u,file=".jorek_filenums."//my_id_s)
  allocate(filenums_tmp(100000)) ! assumes 5-digit numbers
  n=0
  do i=1,100000
    read(u,*,iostat=io) filenums_tmp(i)
    if (io/=0) exit
    n = n + 1
  end do
  allocate(filenums(n))
  filenums(:) = filenums_tmp(1:n)
  deallocate(filenums_tmp)
  close(u, status='delete')

  if (n .le. 0) then
    write(*,*) "No files found!"
    file_number = 0
    return
  end if

  ! Calculate upper and lower bounds
  i_lower = 1 ! index into filenumber array
  write(num_s,'(i0.5)') filenums(i_lower)
  t_lower = get_jorek_hdf5_time('jorek'//num_s//'.h5')*t_norm
  i_upper = n
  write(num_s,'(i0.5)') filenums(i_upper)
  t_upper = get_jorek_hdf5_time('jorek'//num_s//'.h5')*t_norm
  i_guess = nint((time-t_lower)/(t_upper-t_lower)*real(i_upper - i_lower)) + i_lower

  do i=1,20
    if (i_guess .le. 1 .or. i_guess .gt. n) then
      if (my_id .eq. 0) write(*,*) "ERROR: requested time out of range"
      exit
    end if
    if (i_guess .eq. i_lower .or. i_guess .eq. i_upper) then
      file_number = filenums(i_lower)
      return
    end if

    write(num_s,'(i0.5)') filenums(i_guess)
    t_guess = get_jorek_hdf5_time('jorek'//num_s//'.h5')*t_norm
    if (my_id .eq. 0) write(*,"(i5,A,g11.4,A,i5,A,g11.4,A,i5,A,g11.4,A)") i_lower, " (", t_lower, &
      ")    ", i_guess, " (", t_guess, &
      ")    ", i_upper, " (", t_upper, ")    "
    ! Based on the value of t_guess, replace either the lower or upper bound
    if (t_guess .le. time) then
      t_lower = t_guess
      i_lower = i_guess
    else
      t_upper = t_guess
      i_upper = i_guess
    end if
    i_guess = i_lower + (i_upper-i_lower)/2
  end do
end function last_file_before_time

!> Get '/t_now' from a file. Does not alter the units in any way
function get_jorek_hdf5_time(filename) result(time)
  use hdf5
  use hdf5_io_module
  character*(*)      , intent(in)  :: filename
  real*8 :: time
  integer(HID_T) :: file
  integer :: hdferr
  call h5open_f(hdferr)
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr)
  call HDF5_real_reading(file,time,'/t_now')
  call h5fclose_f(file,hdferr)
  call h5close_f(hdferr)
end function get_jorek_hdf5_time



!> Read jorek fields from the next restart file
subroutine do_read(this, sim, ev)
  use mod_import_restart
  use phys_module
  use mpi
  class(read_jorek_fields_interp_linear), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), optional :: ev
  character(len=80) :: restart_file
  integer :: i, ierr, my_id
  logical :: file_exists, next_file_found

  real*8 :: t_norm
  t_norm = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

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
  if (.not. allocated(sim%fields%node_list)) allocate(sim%fields%node_list)
  if (.not. allocated(sim%fields%element_list)) allocate(sim%fields%element_list)
  
  ! Continue for jorek_fields_interp_linear
  select type (f => sim%fields)
  type is (jorek_fields_interp_linear)

    ! Read only one file
    if (this%i .eq. -1 .or. f%static) then
      if (this%i .eq. -1) then
        write(restart_file,'(A,A)') trim(this%basename), '_restart.h5'
      else
        write(restart_file,'(A,i5.5,A)') trim(this%basename), this%i, '.h5'
      end if
      inquire(file=trim(restart_file), exist=file_exists)
      if (file_exists) then
        call import_hdf5_restart(f%node_list,f%element_list,restart_file,this%rst_format,my_id,ierr)
        f%static = .true.
      else
        if (my_id .eq. 0) write(*,*) "ERROR: file ", trim(restart_file), " does not exist"
        call exit(1)
      end if
    else ! Linearly interpolating case
      ! If nothing has been loaded (i.e. fields%time_prev = 0.d0) load the initial file
      if (abs(f%time_prev) .lt. 1.d-50) then
        write(restart_file,'(A,i5.5,A)') trim(this%basename), this%i, '.h5'
        inquire(file=trim(restart_file), exist=file_exists)
        if (file_exists) then
          call import_hdf5_restart(f%node_list,f%element_list,trim(restart_file),this%rst_format,my_id,ierr)
          if (ierr .ne. 0) then
            if (my_id .eq. 0) write(*,*) "ERROR: cannot open restart file"
            call exit(1)
          else
            f%time_now = t_start*t_norm ! set by import_hdf5_restart
            ! Set sim%time to this time also, to start at the right point
            if (sim%time .gt. 1d-16) then ! check if this is the right file if we have already set a time
              if (sim%time .le. f%time_now) then
                if (my_id .eq. 0) write(*,*) "ERROR: restart file read that is too far in the future"
              end if
            else ! otherwise set the time to the time of this file
              sim%time = f%time_now
            end if
            if (my_id .eq. 0) write(*,"(A,f9.8,A)") "Read initial restart file, set t=", f%time_now, " [s]"
          end if
        else
          if (my_id .eq. 0) write(*,*) "ERROR: cannot read initial file ", trim(restart_file)
          call exit(1)
        end if
      end if
      
      ! Find the following file (next timestep number)
      next_file_found=.false.
      do i=this%i+1,this%i+20 ! check 20 files ahead
        write(restart_file,'(A,i5.5,A)') trim(this%basename), i, '.h5'
        inquire(file=trim(restart_file), exist=file_exists)
        if (file_exists) then
          next_file_found=.true.
          call merge_restart(f%node_list, f%element_list, trim(restart_file), this%rst_format,my_id, ierr)
          if (ierr .ne. 0) then
            if (my_id .eq. 0) write(*,*) "ERROR: cannot open restart file"
            call exit(1)
          else
            f%time_prev = f%time_now
            f%time_now = t_start*t_norm ! Set by import_merge_restart
            if (my_id .eq. 0 .and. i-this%i .gt. 1) write(*,"(i2,A)") i-this%i, " JOREK steps between restarts"
            this%i=i
            ! set the time to run this event at next
            if (my_id .eq. 0) write(*,"(A,f9.8,A)") " Read next restart file, values until t=", f%time_now, " [s]"
            if (present(ev)) then
              ev%start = f%time_now
              if (my_id .eq. 0) write(*,*) "Set time for next restart file read to ", ev%start
            end if
            exit ! the file-finding loop
          end if
        endif
      enddo
      if (.not. next_file_found .and. this%stop_at_end) then
        if (my_id .eq. 0) write(*,*) "WARNING: cannot find any next restart files. Stopping."
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if
      if (.not. next_file_found .and. .not. this%stop_at_end) then
        if (my_id .eq. 0) write(*,*) "WARNING: cannot find any next restart files. Continuing with &
        the last values without time-dependence"
        f%static = .true.
      end if
    end if

    call update_neighbours(f%element_list, f%node_list)
  class default
    if (my_id .eq. 0) write(*,*) "ERROR, do_read called with wrong sim%fields"
  end select
end subroutine do_read



!> Import a binary restart file and merges it with the values currently known
!> This can then be used to interpolate linearly between any two restart files
subroutine merge_restart(node_list,element_list, restart_file, format_rst,my_id, ierr)
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
  integer,                 intent(in)    :: my_id

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
  call import_hdf5_restart(node_list,element_list, restart_file, format_rst, my_id, ierr)

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
