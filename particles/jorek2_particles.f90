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
use mod_redistribute_particles
use mod_import_restart_linear
use openadas
use mod_coronal
use clock_module
use mpi_mod
use mod_random_seed
!$ use omp_lib

implicit none


type (type_particle_list) :: particle_list
type (type_particle_list) :: particle_list_GC

integer    :: i_tor, my_id, n_cpu, ierr, i_step, i_step_out, i_begin, i_end, i, index_prev
integer*4  :: rank, comm_size
integer    :: required, provided, StatInfo
character*17 :: particle_file, restart_file

real*8, dimension(:), allocatable :: energy_list, momentum_list
real*8 :: wstart, wend, wtime ! wall time on this cpu

type(type_ADF11_all) :: adf11(1:N_species)
type(type_coronal)   :: coronal(1:N_species)

interface
  subroutine update_particles(my_id,particle_list,t_step,n_step,adf11,energy_list,momentum_list,toroidal_field_factor,field_interp_time)
    use mod_particles
    use openadas
    ! -- Routine parameters
    integer, intent(in)       :: my_id              !< Id of the current process
    type (type_particle_list) :: particle_list      !< The particles we will march forward in time
    real*8,  intent(in)       :: t_step             !< The size of each timestep
    integer, intent(in)       :: n_step             !< The number of timesteps we will perform
    type (type_adf11_all), intent(in) :: adf11
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

! Filename hardcoded here because we can only read one file from stdin easily
call initialise_parameters(my_id, "in_jorek")
call initialise_particle_parameters(my_id, "__NO_FILENAME__")
call initialise_basis

! Seed random numbers for particle initialisation
if (my_id .eq. 0) then
  do i=1,n_species
    if (particle_seed(i) .eq. 0) call gen_random_seed(particle_seed(i))
  enddo
endif

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
call broadcast_particle_parameters(my_id)          ! particle parameters
call update_neighbours(element_list,node_list)     ! update neighbour information in the element_list
do i=1,n_species ! For each particle read adas files and calculate coronal equilibrium)
  if (len_trim(adas_suffix(i)) .ne. 0) then
    adf11(i)   = read_adf11(adas_suffix(i))    ! read openadas data for ionisation, recombination and radiation rates
    coronal(i) = coronal_equilibrium(adf11(i)) ! calculate the coronal equilibria from the adas data
  endif
enddo

call MPI_Barrier(MPI_COMM_WORLD,ierr) ! for output niceness

if (len_trim(particle_restart_file) .eq. 0) then
  call initialise_particles(my_id, n_cpu, coronal, particle_list, particle_list_GC)
else
  call import_particles(particle_restart_file, particle_list)
endif
allocate(energy_list(particle_list%n_particles), momentum_list(particle_list%n_particles))
call MPI_Barrier(MPI_COMM_WORLD,ierr)

! Output particles at start
write(particle_file,'(A3,i9.9,A4)') 'pos',max(t_particles_begin,0),'.vtk'
call particles_vtk(particle_list,particle_file)
call MPI_Barrier(MPI_COMM_WORLD,ierr)

! If t_particles_begin is set ignore nout_particles and n_step_particles
if (t_particles_begin .gt. -1) then
  i_begin = t_particles_begin ! subtract this again later, bit hacky
  i_end   = t_particles_end
  i_step_out = 0
else
  i_begin=1
  i_end=n_step_particles/nout_particles
endif

! Loop n_step_particles/nout_particles in old mode, t_particles_end-t_particles_begin in new mode
do i_step=i_begin,i_end
  if (t_particles_begin .gt. -1) then
    if (i_step_out .gt. i_step) cycle
    if (my_id .eq. 0) then
      call import_next_restart(node_list,element_list, i_step, i_step_out, rst_format) !  returns i_step of new file
    endif
    call MPI_Bcast(i_step_out, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) ! broadcast this number to every node
    call broadcast_phys(my_id)                         ! physics parameters, because tstep might have changed
    call broadcast_nodes(my_id, node_list)
    call broadcast_elements(my_id, element_list)
    call update_neighbours(element_list,node_list)
    ! Set nout_particles to the number of steps required to go t_step forward (floored)
    nout_particles = int(tstep/t_step_particles,4)
    ! Set t_step_particles to the closest integer divisor of t_step so we don't miss a substep
    t_step_particles = tstep/nout_particles
  else
    i_step_out = i_step
  endif

  ! Do substepping
  call cpu_time(wstart) ! correct if no OMP
  !$ wstart = omp_get_wtime()
  call update_particles(my_id,particle_list,t_step_particles,nout_particles,adf11(1),& ! XXX hardcoded first adf11 index
      energy_list,momentum_list,field_interp_time=(t_particles_begin .gt. -1))
  call cpu_time(wend) ! correct if no OMP
  !$ wend = omp_get_wtime()
  wtime = wend - wstart
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  call cpu_time(wstart) ! correct if no OMP
  !$ wstart = omp_get_wtime()
  ! Redistribute particles over processors based on the walltime for the updating
  !call redistribute_particles(particle_list,wtime)
  call cpu_time(wend) ! correct if no OMP
  !$ wend = omp_get_wtime()
  wtime = wend - wstart
  !if (my_id .eq. 0) write(*,*) "Particle load-balancing took ", wtime, " seconds"
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  write(particle_file,'(A4,i0.9,A4)') 'part',i_step_out,'.rst'
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
