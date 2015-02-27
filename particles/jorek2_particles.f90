program particle_transport

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

integer    :: i, in, i_tor, my_id, n_cpu, ierr
integer*4  :: rank, comm_size
integer    :: required, provided, StatInfo
real*8     :: boxsize(3)

write(*,*) '***************************************'
write(*,*) '* JOREK2 : Particles                  *'
write(*,*) '***************************************'

call begplt('part.ps')
call clck_init()

my_id=0

#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif

call MPI_Init_thread(required, provided, StatInfo)

call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
my_id = rank

call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
n_cpu = comm_size

call initialise_parameters(my_id, "__NO_FILENAME__")

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
  write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
enddo

call import_restart(node_list,element_list, 'jorek_restart.rst', rst_format, ierr)

call initialise_basis                              ! define the basis functions at the Gaussian points

call update_neighbours(element_list,node_list)     ! update neighbour information in the element_list

call read_adas                                     ! read openadas data for ionisation, recombination and radiation rates

call coronal                                       ! calculate the coronal equilibria from the adas data

boxsize = 0.
do i=1,node_list%n_nodes
  boxsize(1) = max(boxsize(1),node_list%node(i)%x(1,1))
  boxsize(3) = max(boxsize(3),abs(node_list%node(i)%x(1,2)))
enddo
boxsize(2) = boxsize(1)

write(*,'(A,3e14.6)') ' boxsize : ',boxsize

call initialise_particles(my_id, n_cpu, node_list, element_list, particle_list, boxsize)

t_step_particles = 0.01
n_step_particles = 10000

call export_particles(particle_list,'part000.rst')
call particles_vtk(particle_list,'part000.vtk')

call update_particles(my_id,particle_list,t_step_particles,n_step_particles)

call export_particles(particle_list,'part100.rst')

call particles_vtk(particle_list,'part100.vtk')

call finplt

end
