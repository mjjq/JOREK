!> Module for Hermite-Birkhoff spline interpolation in time
!> in JOREK restart files. Contains an action to read the fields.
!>
!> This module calculates a finite-difference approximation to the time-derivative
!> and stores this in the deltas. We need 3 restart files to contain all the info.
!>
!> We have 3 different node lists (the one in fields_base and 2 others) and cycle
!> through these. The reader does some preprocessing to calculate the derivatives
!> and store these in deltas. If the number of steps between restart files is 1
!> we only need to read every other file and can keep the same algorithm.
!>
!> WARNING: There is no guarantee that fields%node_list contains the values at the
!> current timestep. The only use for this is to get the structure of the grid,
!> which is assumed not to change. For all other uses please call the interp_PRZ
!> function.
module mod_fields_hermite_birkhoff
use data_structure
use mod_particle_sim
use mod_event
use mod_fields
use mod_interp4
use mod_interp_PRZ
use mod_fields_linear
use mod_hermite_birkhoff
implicit none
private
public jorek_fields_interp_hermite_birkhoff, read_jorek_fields_interp_hermite_birkhoff

!> Action to read in the fields into sim%fields
type, extends(action) :: read_jorek_fields_interp_hermite_birkhoff
  character(len=80) :: basename = 'jorek' !< Comes before the file number or extension
  integer :: i = 0 !< Number of the restart file to read. Set to -1 to not include. Corresponds to the index of NOW
  integer :: rst_format = 0 !< Format of restart file if .rst type
  logical :: stop_at_end = .true. !< Whether to stop the simulation at the end of the file list
  contains
    procedure :: do => do_read
end type
interface read_jorek_fields_interp_hermite_birkhoff
  module procedure new_read_jorek_fields_interp_hermite_birkhoff
end interface read_jorek_fields_interp_hermite_birkhoff

!> Store enough data for interpolation in time.
!> This does not work for changing grids in time!
!> Use at your own peril in that case. (e.g. i_elm and s,t for a specific spatial position might depend on time)
!> {node,element}_list_{B,C} contain the other lists to interpolate with.
!> The order is always cyclical, and can be described with the number describing
!> the oldest file. (1=A, 2=B, 3=C)
type, extends(fields_base) :: jorek_fields_interp_hermite_birkhoff
  real*8 :: t_A, t_B, t_C !< Time of restart files (SI).
  real*8 :: dt_AB, dt_BC !< Time between restart files (SI)
  integer*1 :: oldest = 2 !< B is the oldest file in the beginning, A is the newest

  type(type_node_list), allocatable    :: node_list_B, node_list_C
  type(type_element_list), allocatable :: element_list_B, node_list_C
  contains
    procedure :: interp_PRZ => do_interp_PRZ
end type jorek_fields_interp_hermite_birkhoff
contains

!> Interpolate a variable at a specific position (with phi), with first derivatives only
pure subroutine do_interp_PRZ(this, time, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
  use mod_interp_PRZ
  use constants, only: mu_zero, mass_proton
  use phys_module, only: tstep, central_mass, central_density
  class(jorek_fields_interp_hermite_birkhoff),  intent(in)  :: this
  real*8,                   intent(in)  :: time !< Time at which to calculate this variable
  integer,                  intent(in)  :: i_elm
  integer,                  intent(in)  :: n_v, i_v(n_v)
  real*8,                   intent(in)  :: s, t, phi
  real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v), P_phi(n_v), P_time(n_v)
  real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t

  real*8, dimension(n_v,4,4) :: V !< n_var, (P, P_s, P_t, P_phi), (interp, interp, delta, delta)
  real*8 :: t_norm
  integer :: i, j

  t_norm = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

  if (this%other_latest) then
    i = 1
    j = 2
  else
    i = 2
    j = 1
  end if
  call       interp_PRZ(this%node_list,this%element_list,            i_elm,i_v,n_v,s,t,phi, &
    V(:,1,i),   V(:,2,i),   V(:,3,i),   V(:,4,i),   R,R_s,R_t,Z,Z_s,Z_t)
  call       interp_PRZ(this%node_list_other,this%element_list_other,i_elm,i_v,n_v,s,t,phi, &
    V(:,1,j),   V(:,2,j),   V(:,3,j),   V(:,4,j),   R,R_s,R_t,Z,Z_s,Z_t)
  call interp_PRZ_delta(this%node_list,this%element_list,            i_elm,i_v,n_v,s,t,phi, &
    V(:,1,i+2), V(:,2,i+2), V(:,3,i+2), V(:,4,i+2), R,R_s,R_t,Z,Z_s,Z_t)
  call interp_PRZ_delta(this%node_list_other,this%element_list_other,i_elm,i_v,n_v,s,t,phi, &
    V(:,1,j+2), V(:,2,j+2), V(:,3,j+2), V(:,4,j+2), R,R_s,R_t,Z,Z_s,Z_t)
  ! Transform to time-derivatives:
  V(:,:,3) = V(:,:,3)/this%dt_prev
  V(:,:,4) = V(:,:,4)/this%dt_next

  call HB_interp(this%time_prev, this%time_next, n_v, V(:,1,1), V(:,1,2), V(:,1,3), V(:,1,4), time, P)
  call HB_interp(this%time_prev, this%time_next, n_v, V(:,2,1), V(:,2,2), V(:,2,3), V(:,2,4), time, P_s)
  call HB_interp(this%time_prev, this%time_next, n_v, V(:,3,1), V(:,3,2), V(:,3,3), V(:,3,4), time, P_t)
  call HB_interp(this%time_prev, this%time_next, n_v, V(:,4,1), V(:,4,2), V(:,4,3), V(:,4,4), time, P_phi)
  call HB_interp_dt(this%time_prev, this%time_next, n_v, V(:,1,1), V(:,1,2), V(:,1,3), V(:,1,4), time, P_time)
end subroutine do_interp_PRZ


!> Constructor to allow for optional and default variables
function new_read_jorek_fields_interp_hermite_birkhoff(basename, i, rst_format, stop_at_end) result(new)
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional :: i
  integer, intent(in), optional :: rst_format
  logical, intent(in), optional :: stop_at_end
  type(read_jorek_fields_interp_hermite_birkhoff) :: new
  if (present(basename)) new%basename = basename
  if (present(i)) new%i = i
  if (present(rst_format)) new%rst_format = rst_format
  if (present(stop_at_end)) new%stop_at_end = stop_at_end
  new%name = "ReadJorekFieldsInterpHermiteBirkhoff"
  new%log = .true.
end function new_read_jorek_fields_interp_hermite_birkhoff



!> Read jorek fields from the next restart file
subroutine do_read(this, sim, ev)
  use mod_import_restart
  use phys_module
  use mpi
  class(read_jorek_fields_interp_hermite_birkhoff), intent(inout) :: this
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
    type is (jorek_fields_interp_hermite_birkhoff) ! do nothing
    class default
      write(*,*) "WARNING: wrong type of fields in particle%sim, reallocating"
      deallocate(sim%fields)
      allocate(jorek_fields_interp_hermite_birkhoff::sim%fields)
    end select
  else
    allocate(jorek_fields_interp_hermite_birkhoff::sim%fields)
  end if
  if (.not. allocated(sim%fields%node_list)) allocate(sim%fields%node_list)
  if (.not. allocated(sim%fields%element_list)) allocate(sim%fields%element_list)

  ! Continue for jorek_fields_interp_hermite_birkhoff
  select type (f => sim%fields)
  type is (jorek_fields_interp_hermite_birkhoff)
    if (.not. allocated(f%node_list_B)) allocate(f%node_list_B)
    if (.not. allocated(f%element_list_B)) allocate(f%element_B)
    if (.not. allocated(f%node_list_C)) allocate(f%node_list_C)
    if (.not. allocated(f%element_list_C)) allocate(f%element_C)

    ! If nothing has been loaded (i.e. fields%time_prev = 0.d0) load the initial file
    if (abs(f%time_prev) .lt. 1.d-50) then
      write(restart_file,'(A,i5.5,A)') trim(this%basename), this%i, '.h5'
      inquire(file=trim(restart_file), exist=file_exists)
      if (file_exists) then
        call import_hdf5_restart(f%node_list,f%element_list,trim(restart_file),this%rst_format,my_id,ierr)
        f%dt_next = tstep*t_norm
        call update_neighbours(f%element_list, f%node_list)
        if (ierr .ne. 0) then
          if (my_id .eq. 0) write(*,*) "ERROR: cannot open restart file"
          call exit(1)
        else
          f%time_prev = t_start*t_norm ! set by import_hdf5_restart
          ! Set sim%time to this time also, to start at the right point
          if (sim%time .gt. 1d-16) then ! check if this is the right file if we have already set a time
            if (sim%time .le. f%time_prev) then
              if (my_id .eq. 0) write(*,*) "ERROR: restart file read that is too far in the future"
            end if
          else ! otherwise set the time to the time of this file
            sim%time = f%time_prev
          end if
          if (my_id .eq. 0) write(*,"(A,f9.8,A)") "Read initial restart file, set t=", f%time_prev, " [s]"
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
        if (.not. f%other_latest) then ! if f%node_list contained the latest values
          call import_hdf5_restart(f%node_list_other,f%element_list_other,trim(restart_file),this%rst_format,my_id,ierr)
          call update_neighbours(f%element_list_other, f%node_list_other)
        else
          call import_hdf5_restart(f%node_list,f%element_list,trim(restart_file),this%rst_format,my_id,ierr)
          call update_neighbours(f%element_list, f%node_list)
        end if
        f%dt_prev = f%dt_next
        f%dt_next = tstep*t_norm
        f%other_latest = .not. f%other_latest
        if (ierr .ne. 0) then
          if (my_id .eq. 0) write(*,*) "ERROR: cannot open restart file"
          call exit(1)
        else
          f%time_next = t_start*t_norm ! Set by import_merge_restart
          if (my_id .eq. 0 .and. i-this%i .gt. 1) write(*,"(i2,A)") i-this%i, " JOREK steps between restarts"
          if (my_id .eq. 0 .and. i-this%i .le. 2) write(*,*) "WARNING: not enough JOREK steps between restart files for HB interpolation."
          this%i=i
          ! set the time to run this event at next
          if (my_id .eq. 0) write(*,"(A,f9.8,A)") " Read next restart file, values until t=", f%time_next, " [s]"
          if (present(ev)) then
            ev%start = f%time_next
            if (my_id .eq. 0) write(*,*) "Set time for next restart file read to ", ev%start
          end if
          exit ! the file-finding loop
        end if
      endif
    enddo
    if (.not. next_file_found) then
      if (my_id .eq. 0) write(*,*) "WARNING: cannot find any next restart files. Stopping."
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

  class default
    if (my_id .eq. 0) write(*,*) "ERROR, do_read called with wrong sim%fields"
  end select
end subroutine do_read
end module mod_fields_hermite_birkhoff
