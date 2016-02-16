!> Project particles onto the elements in JOREK and export to a vtk file
!! It uses the elements of file jorek_restart.rst
!! Run as project_particles_vtk < jorek_in
program project_particles_vtk
use data_structure
use phys_module
use basis_at_gaussian
use nodes_elements
use mod_particles
use mod_project_particles
use mod_import_export_particles
implicit none

type (type_particle_list) :: particle_list

integer    :: i, j, i_tor, my_id, n_cpu, ierr, i_step, i_begin, i_end
integer*4  :: rank, comm_size
integer    :: required, provided, StatInfo
real*8     :: boxwidth(3), boxcenter(3) !< size and center of box in RZphi space
real*8     :: substep
character*17 :: particle_file, filenum, restart_file


required = MPI_THREAD_MULTIPLE

call MPI_Init_thread(required, provided, StatInfo)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
my_id = rank
call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
n_cpu = comm_size

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '* JOREK project particles to vtk      *'
  write(*,*) '***************************************'
endif

call initialise_parameters(my_id, "__NO_FILENAME__")
call initialise_basis                              ! define the basis functions at the Gaussian points

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
  if (my_id .eq. 0) write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
enddo

restart_file = 'jorek_restart.rst'
call import_binary_restart(node_list,element_list, restart_file, rst_format, ierr)
if (ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,ierr)

call broadcast_elements(my_id, element_list)       ! elements
call broadcast_nodes(my_id, node_list)             ! nodes
call broadcast_phys(my_id)                         ! physics parameters

! TODO add full support for tstep_n (also in calc_EB.f90) or get the time from
! the jorek restart files
! If t_particles_begin is set ignore nout_particles and n_step_particles
if (t_particles_begin .gt. -1) then
  i_begin = t_particles_begin + 1 ! Nota bene! we will start at the second restart file as this contains the fields of the first too
  i_end   = t_particles_end
  ! Set nout_particles to the number of steps required to go t_step forward (floored)
  nout_particles = tstep_n(1)/t_step_particles
  ! Set t_step_particles to the closest integer divisor of t_step so we don't
  ! miss a substep
  t_step_particles = tstep_n(1)/nout_particles

  if (tstep_n(2) .gt. 0) write(*,*) "WARNING: No full support for tstep_n"
else
  write(*,*) "nout_particles mode not supported"
endif

! Loop n_step_particles/nout_particles in old mode, t_particles_end-t_particles_begin in new mode
do i_step=i_begin,i_end
  ! Read the particle file
  write(particle_file,'(A4,i0.9,A4)') 'part',i_step,'.rst'
  call import_particles(particle_file, particle_list)

  ! Project onto element_list (saves into the first value in node_list)
  call project_particles(particle_list)

  ! Write the first value of node_list to a vtk file
  call 
enddo


call MPI_FINALIZE(IERR)
end program project_particles_vtk
