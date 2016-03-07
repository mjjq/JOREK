!> jorek2_particles is a post_processing tool for test particles
!! It uses the fields in the file jorek_restart.rst
!! Set t_step_particles and n_step_particles, and nout for output control
program jorek2_particles

use data_structure
use phys_module
use basis_at_gaussian
use nodes_elements
use mod_particles
use mod_import_export_particles
use mod_initialise_particles
use openadas
use clock_module
use mpi_mod

implicit none


type (type_particle_list) :: particle_list
type (type_particle_list) :: particle_list_GC

integer    :: i_tor, my_id, n_cpu, ierr, i_step, i_begin, i_end
integer*4  :: rank, comm_size
integer    :: required, provided, StatInfo
character*17 :: particle_file, restart_file

real*8, dimension(:), allocatable :: energy_list, momentum_list

interface
  subroutine update_particles(my_id, particle_list, t_step, n_step, energy_list, momentum_list, toroidal_field_factor, field_interp_time)
    use mod_particles
    ! -- Routine parameters
    type (type_particle_list) :: particle_list      !< The particles we will march forward in time
    real*8,  intent(in)       :: t_step             !< The size of each timestep
    integer, intent(in)       :: n_step             !< The number of timesteps we will perform
    integer, intent(in)       :: my_id              !< Id of the current process
    real*8,  intent(out), dimension(:), optional :: energy_list !< Energy of the particles at the next-to(!) final timestep
    real*8,  intent(out), dimension(:), optional :: momentum_list !< Generalized toroidal momentum of the particles at the next-to(!) final timestep
    real*8,  intent(in),  optional :: toroidal_field_factor !< Multiply B_phi with this WARNING: use only for testing!
    logical, intent(in),  optional :: field_interp_time !< Interpolate the fields linearly in time as if the first step was in the previous fields (almost) and the last in the current
  end subroutine update_particles
end interface


required = MPI_THREAD_MULTIPLE

call MPI_Init_thread(required, provided, StatInfo)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
my_id = rank
call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
n_cpu = comm_size

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '* JOREK2_particles                    *'
  write(*,*) '***************************************'
endif

call initialise_parameters(my_id, "__NO_FILENAME__")
call read_adas                                     ! read openadas data for ionisation, recombination and radiation rates
call initialise_basis                              ! define the basis functions at the Gaussian points
call coronal                                       ! calculate the coronal equilibria from the adas data

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
  if (my_id .eq. 0) write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
enddo

if (my_id .eq. 0) then
  if (t_particles_begin .eq. -1) then ! special value for old behaviour, default
    !value of t_particles_begin
    restart_file = 'jorek_restart.rst'
  else
    write(restart_file,'(A,i5.5,A)') 'jorek', t_particles_begin, '.rst'
  endif
  call import_binary_restart(node_list,element_list, restart_file, rst_format, ierr)
  if (ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,ierr)
endif
call broadcast_elements(my_id, element_list)       ! elements
call broadcast_nodes(my_id, node_list)             ! nodes
call broadcast_phys(my_id)                         ! physics parameters
call update_neighbours(element_list,node_list)     ! update neighbour information in the element_list
call MPI_Barrier(MPI_COMM_WORLD,ierr) ! for output niceness

if (len_trim(particle_restart_file) .eq. 0) then
  call initialise_particles(my_id, n_cpu, particle_list, particle_list_GC)
else
  call import_particles(particle_restart_file, particle_list)
endif
allocate(energy_list(particle_list%n_particles), momentum_list(particle_list%n_particles))
call MPI_Barrier(MPI_COMM_WORLD,ierr)

! Output particles at start
write(particle_file,'(A3,i9.9,A4)') 'pos',max(t_particles_begin,0),'.vtk'
call particles_vtk(particle_list,particle_file)
call MPI_Barrier(MPI_COMM_WORLD,ierr)

! TODO add full support for tstep_n (also in calc_EB.f90) or get the time from
! the jorek restart files
! If t_particles_begin is set ignore nout_particles and n_step_particles
if (t_particles_begin .gt. -1) then
  i_begin = t_particles_begin + 1 ! Nota bene! we will start at the second restart file as this contains the fields of the first too
  i_end   = t_particles_end
  ! Set nout_particles to the number of steps required to go t_step forward (floored)
  nout_particles = int(tstep_n(1)/t_step_particles,4)
  ! Set t_step_particles to the closest integer divisor of t_step so we don't
  ! miss a substep
  t_step_particles = tstep_n(1)/nout_particles
  write(*,*) "Using actual particle timestep", t_step_particles

  if (tstep_n(2) .gt. 0) write(*,*) "WARNING: No full support for tstep_n"
else
  i_begin=1
  i_end=n_step_particles/nout_particles
endif

! Loop n_step_particles/nout_particles in old mode, t_particles_end-t_particles_begin in new mode
do i_step=i_begin,i_end
  if (t_particles_begin .gt. -1) then
    if (my_id .eq. 0) then
      write(restart_file,'(A,i5.5,A)') 'jorek', i_step, '.rst'
      call import_binary_restart(node_list,element_list, restart_file, rst_format, ierr)
      if (ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,ierr)
    endif
    call broadcast_nodes(my_id, node_list)
    call broadcast_elements(my_id, element_list)
    call update_neighbours(element_list,node_list)
  endif

  ! Do substepping
  call update_particles(my_id,particle_list,t_step_particles,nout_particles,energy_list,momentum_list,field_interp_time=(t_particles_begin .gt. -1))

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  write(particle_file,'(A4,i0.9,A4)') 'part',i_step,'.rst'
  call export_particles(particle_list,particle_file)

  ! Output by each processor
  if (write_energies) call write_list("energy",i_step,my_id,energy_list)
  if (write_momenta) call write_list("momentum",i_step,my_id,momentum_list)
enddo


call MPI_FINALIZE(IERR)
contains
subroutine write_list(ftype,fnum,my_id,list)
  implicit none
  ! Input parameters
  integer, intent(in) :: fnum, my_id
  character(len=*), intent(in) :: ftype
  real*8, dimension(:) :: list

  character*40 :: filename

  write(filename,"(A,I0.9,A,I0.4,A)") trim(ftype), fnum, "_", my_id, ".dat"
  open(file=filename,status="replace",unit=21,access="stream",form='unformatted')
  write(21) list
  close(21)
end subroutine write_list
end program jorek2_particles
