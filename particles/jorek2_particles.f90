!> jorek2_particles is a post_processing tool for test particles
!! It uses the fields in the file jorek_restart.rst
!! Set t_step_particles and n_step_particles, and nout for output control
!! Initialization of particles is set using the boxcenter and boxwidth variables
!! (hardcoded now)
program jorek2_particles

use data_structure
use phys_module
use basis_at_gaussian
use nodes_elements
use mod_particles
use openadas
use clock_module
use mpi_mod

implicit none


type (type_particle_list):: particle_list

integer    :: i, j, i_tor, my_id, n_cpu, ierr, i_step
integer*4  :: rank, comm_size
integer    :: required, provided, StatInfo
real*8     :: boxwidth(3), boxcenter(3) !< size and center of box in RZphi space
character*17 :: particle_file, filenum

real*8, dimension(:), allocatable :: energy_list, momentum_list

interface
  subroutine update_particles(my_id, particle_list, t_step, n_step, energy_list, momentum_list, toroidal_field_factor)
    use mod_particles
    ! -- Routine parameters
    type (type_particle_list) :: particle_list      !< The particles we will march forward in time
    real*8,  intent(in)       :: t_step             !< The size of each timestep
    integer, intent(in)       :: n_step             !< The number of timesteps we will perform
    integer, intent(in)       :: my_id              !< Id of the current process
    real*8,  intent(out), dimension(:), optional :: energy_list !< Energy of the particles at the next-to(!) final timestep
    real*8,  intent(out), dimension(:), optional :: momentum_list !< Generalized toroidal momentum of the particles at the next-to(!) final timestep
    real*8,  intent(in),  optional :: toroidal_field_factor !< Multiply B_phi with this WARNING: use only for testing!
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

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
  write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
enddo

if (my_id .eq. 0) then
  call import_binary_restart(node_list,element_list, 'jorek_restart.rst', rst_format, ierr)
endif

call broadcast_elements(my_id, element_list)       ! elements
call broadcast_nodes(my_id, node_list)             ! nodes
call broadcast_phys(my_id)                         ! physics parameters
call initialise_basis                              ! define the basis functions at the Gaussian points
call update_neighbours(element_list,node_list)     ! update neighbour information in the element_list
call read_adas                                     ! read openadas data for ionisation, recombination and radiation rates
call coronal                                       ! calculate the coronal equilibria from the adas data

! Boxsize is R,Z location and dPhi extent from phi=0
boxcenter = (/2.83d0, 0.d0, 0.d0/)
boxwidth = (/1.03d0, 1.9d0, TWOPI/)


call initialise_particles(my_id, n_cpu, node_list, element_list, particle_list, boxcenter, boxwidth, n_particles)
allocate(energy_list(particle_list%n_particles), momentum_list(particle_list%n_particles))
call MPI_Barrier(MPI_COMM_WORLD,ierr)


write(particle_file,'(A4,i7.7,A4)') 'part',0,'.vtk'
call particles_vtk(particle_list,particle_file)

! TODO be helpful when nout_particles > n_step_particles
do i_step=1,n_step_particles/nout_particles
  call update_particles(my_id,particle_list,t_step_particles,nout_particles,energy_list,momentum_list)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  write(particle_file,'(A4,i7.7,A4)') 'part',i_step*nout_particles,'.vtk'
  call particles_vtk(particle_list,particle_file)

enddo


call MPI_FINALIZE(IERR)
end program jorek2_particles
