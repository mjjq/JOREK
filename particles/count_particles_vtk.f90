!> Count particles per element and export to VTK file
!! It uses the elements of file jorek_restart.rst
!! Run as count_particles_vtk < jorek_in
!! It reads some parameters from the vtk.nml namelist in the current directory
!! The resutl is quite inaccurate, but good for quick tests
program project_particles_vtk
use phys_module
use nodes_elements
use mod_particles
use mod_import_export_particles
use mod_vtk
use mpi_mod
implicit none

type (type_particle_list) :: particle_list
integer    :: ierr, provided
character*17 :: particle_file, restart_file
integer    :: i_t, i



call MPI_Init_thread(MPI_THREAD_SINGLE, provided, ierr)

write(*,*) '***************************************'
write(*,*) '* JOREK project particles to vtk      *'
write(*,*) '***************************************'

call initialise_parameters(0, "__NO_FILENAME__")

do i_t=1, n_tor
  mode(i_t) = + int(i_t / 2) * n_period
  write(*,*) ' toroidal mode numbers : ',i_t,mode(i_t)
enddo

restart_file = 'jorek_restart.rst'
call import_binary_restart(node_list,element_list, restart_file, rst_format, ierr)
if (ierr .ne. 0) call exit(1)

! Get the filename as the first cli argument
if (command_argument_count() < 1) then
  write(*,*) "Expected a filename argument"
  call exit(1)
endif

do i=1,command_argument_count()
  ! Get particle filename from commandline
  call get_command_argument(i, particle_file)
  call import_particles(particle_file, particle_list)

  ! Write the first value of node_list to a vtk file
  particle_file = particle_file(1:index(particle_file,'.rst',.true.))//'vtk' !  .true. searches backwards
  call write_particle_counts_to_vtk(node_list,element_list,particle_list,particle_file)
  write(*,*) "Done counting, wrote output to ", particle_file
enddo

contains
subroutine write_particle_counts_to_vtk(node_list,element_list,particle_list,filename)
use data_structure
use basis_at_gaussian ! for HZ (initialise_basis must be called before use)
use mod_vtk
use mod_particles
use basis_at_gaussian
use constants
implicit none

!> Input parameters
type(type_node_list), intent(in)      :: node_list
type(type_element_list), intent(in)   :: element_list
type (type_particle_list), intent(in) :: particle_list
character*(*), intent(in)             :: filename

integer :: nnos, nnoel, nel, i, j, ielm, inode, k
real*4,allocatable    :: xyz (:,:), scalars(:,:)
real*8 :: s, t, R, R_s, R_t, Z, Z_s, Z_t
real*8 :: P, P_s, P_t, P_st, P_ss, P_tt
character*12, allocatable :: scalar_names(:)
integer :: i_elm, weight

real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)

real*8  :: total_volume, total_area, volume, area, xjac, wst
integer :: ms, mt

type(type_node) :: node
type(type_element) :: element

call initialise_basis                              ! define the basis functions at the Gaussian points

nnos = element_list%n_elements
allocate(xyz(3,nnos),scalars(nnos,1))
allocate(scalar_names(1))
scalar_names(1) = "num_density"

scalars = 0.0
xyz     = 0
total_area = 0.0
total_volume = 0.0

! Count number of particles in each element
!$omp parallel do default(none) shared(particle_list) private(i_elm,weight) &
!$omp   reduction(+:scalars)
do i=1,particle_list%n_particles
  if (particle_list%particle(i)%lost) cycle ! skip this iteration
  i_elm  = particle_list%particle(i)%i_elm
  if (i_elm .lt. 1) cycle
  weight = particle_list%particle(i)%weight
  scalars(i_elm,1) = scalars(i_elm,1) + real(weight,4)
enddo
!$omp end parallel do
write(*,*) particle_list%n_particles, int(sum(scalars))

! Create points for each element
!$omp parallel do default(none) shared(xyz,node_list,element_list, wgauss, &
!$omp                                  H, H_s, H_t, scalars) &
!$omp   private(R, Z, area, volume, wst, ms, mt, xjac, i, j, R_s, R_t, Z_s, Z_t, node, element, &
!$omp           x_g, x_s, x_t, y_g, y_s, y_t) &
!$omp   reduction(+:total_area,total_volume)
do i_elm=1,element_list%n_elements
  call interp_RZ2(node_list,element_list,i_elm,0.5d0,0.5d0,R,R_s,R_t,Z,Z_s,Z_t)
  xyz(1:3,i_elm) = (/R, Z, 0.d0/)

  ! Calculate volume to calculate density
  volume = 0.
  area   = 0.
  x_g = 0.d0; x_s = 0.d0; x_t = 0.d0; y_g = 0.d0; y_s = 0.d0; y_t = 0.d0

  element = element_list%element(i_elm)
  do i=1,n_vertex_max
    node = node_list%node(element%vertex(i))
    do j=1,n_order+1
      do ms=1, n_gauss
        do mt=1, n_gauss
          x_g(ms,mt) = x_g(ms,mt) + node%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
          y_g(ms,mt) = y_g(ms,mt) + node%x(j,2) * element%size(i,j) * H(i,j,ms,mt)

          x_s(ms,mt) = x_s(ms,mt) + node%x(j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + node%x(j,1) * element%size(i,j) * H_t(i,j,ms,mt)
          y_s(ms,mt) = y_s(ms,mt) + node%x(j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + node%x(j,2) * element%size(i,j) * H_t(i,j,ms,mt)
        enddo
      enddo
    enddo
  enddo

  do ms=1, n_gauss
    do mt=1, n_gauss
      wst = wgauss(ms)*wgauss(mt)
      xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

      area   = area   + xjac * wst
      volume = volume + TWOPI * x_g(ms,mt) * xjac * wst
    enddo
  enddo

  ! Use volume to convert scalar counts to number density
  scalars(i_elm,1) = scalars(i_elm,1) / volume

  total_area = total_area + area
  total_volume = total_volume + volume
enddo  ! n_elements
!$omp end parallel do

write(*,*) "Area: ", total_area
write(*,*) "Volume: ", total_volume

call write_vtk(filename,xyz,scalar_names=scalar_names,scalars=scalars)

call MPI_FINALIZE(IERR)
end subroutine write_particle_counts_to_vtk
end program project_particles_vtk
