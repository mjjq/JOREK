!> Project particles onto the elements in JOREK and export to a vtk file
!> It uses the elements of file jorek_restart.rst
!> Run as project_particles_vtk particle_file.rst < jorek_in
!> It reads some parameters from the vtk.nml namelist in the current directory
program project_particles_vtk
use data_structure
use phys_module
use basis_at_gaussian
use nodes_elements
use mod_particles
use mod_project_particles
use mod_particle_io
use mpi
use mod_vtk
implicit none

type (type_particle_list) :: particle_list

integer    :: my_id, n_cpu, ierr
integer*4  :: rank, comm_size
integer    :: required, provided, StatInfo
character*17 :: particle_file, restart_file

!> Parameters
! Not all of these are relevant, but they are here to prevent an error if using
! the same namelist as for regular jorek2vtk
integer :: nsub, i_tor, i_plane, i_t
logical :: without_n0_mode, SI_units
logical :: include_fluxes, include_neo, include_magnetic_field, include_velocity_field,&
           include_bootstrap, include_psi_norm
namelist /vtk_params/ nsub, i_tor, i_plane, without_n0_mode, SI_units, &
                      include_fluxes, include_neo, include_magnetic_field, include_velocity_field,&
                      include_bootstrap, include_psi_norm

required = MPI_THREAD_SINGLE

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



! --- Preset parameters (only these are used!)
nsub                   = 3	 ! Number of subdivisions of the cubic finite elements into linear pieces

! --- Read parameters from namelist file 'vtk.nml' if it exists
open(42, file='vtk.nml', action='read', status='old', iostat=ierr)
if ( ierr == 0 ) then
  write(*,*) 'Reading parameters from vtk.nml namelist.'
  read(42,vtk_params)
  close(42)
end if

write(*,*)
write(*,*) 'Parameters:'
write(*,*) '-----------'
write(*,*) 'nsub            =', nsub
write(*,*) '-----------'
write(*,*) 'n_tor           =', n_tor
write(*,*) 'n_period        =', n_period
write(*,*)

do i_t=1, n_tor
  mode(i_t) = + int(i_t / 2) * n_period
  if (my_id .eq. 0) write(*,*) ' toroidal mode numbers : ',i_t,mode(i_t)
enddo

restart_file = 'jorek_restart.rst'
call import_binary_restart(node_list,element_list, restart_file, rst_format, ierr)
if (ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,ierr)

call broadcast_elements(my_id, element_list)       ! elements
call broadcast_nodes(my_id, node_list)             ! nodes
call broadcast_phys(my_id)                         ! physics parameters

! Get the filename as the first cli argument
if (command_argument_count() < 1) then
  write(*,*) "Expected a filename argument"
  call exit(1)
endif

if (n_tor > n_var) then
  write(*,*) "Too few variables in node_list to save (n_tor > n_var)"
  call exit(2)
endif

! Get particle filename from commandline
call get_command_argument(1, particle_file)
call import_particles(particle_file, particle_list)

! Project onto element_list (saves into the first n_tor values in node_list!)
call project_particles(node_list, element_list, particle_list)

particle_file = particle_file(1:index(particle_file,'.h5',.true.))//'vtk' !  .true. searches backwards
write(*,*) "Done projecting, writing output to ", particle_file
call write_particle_distribution_to_vtk(node_list,element_list,particle_file,nsub)


call MPI_FINALIZE(ierr)
contains

end program project_particles_vtk
