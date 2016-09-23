!> Count particles per element and export to VTK file
!! It uses the elements of file jorek_restart.rst
!! Run as count_particles_vtk (filename)+ < jorek_in
!! It reads some parameters from the vtk.nml namelist in the current directory
!! The result is quite inaccurate, but good for quick tests
program count_particles_vtk
use mpi
implicit none

type (type_particle_list) :: particle_list
integer    :: ierr, provided
character*20 :: particle_file, restart_file, output_file
integer    :: i_t, i, i_p, i_start, i_end
character*2 :: str1

!> Parameters
! Not all of these are relevant, but they are here to prevent an error if using
! the same namelist as for regular jorek2vtk
integer :: nsub, i_tor, i_plane
logical :: without_n0_mode, SI_units
logical :: include_fluxes, include_neo, include_magnetic_field, include_velocity_field,&
           include_bootstrap, include_psi_norm
namelist /vtk_params/ nsub, i_tor, i_plane, without_n0_mode, SI_units, &
                      include_fluxes, include_neo, include_magnetic_field, include_velocity_field,&
                      include_bootstrap, include_psi_norm



call MPI_Init_thread(MPI_THREAD_SINGLE, provided, ierr)

write(*,*) '***************************************'
write(*,*) '* JOREK project particles to vtk      *'
write(*,*) '***************************************'

call initialise_parameters(0, "__NO_FILENAME__")

! --- Preset parameters (only these are used!)
i_plane = 0 ! if 0, count over the entire torus
! if >= 1, count only in immediate surrounding of this plane, and produce output
! in part[0-9]+_{i_plane}.vtk
! if -1, do this for i_plane = 1,n_plane

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
write(*,*) 'i_plane         =', i_plane
write(*,*) '-----------'
write(*,*) 'n_tor           =', n_tor
write(*,*) 'n_period        =', n_period
write(*,*)

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

  if (i_plane .eq. -1) then
    i_start = 1
    i_end = n_plane
  else
    i_start = i_plane
    i_end = i_plane
  endif
  do i_p = i_start, i_end
    if (i_p .eq. 0) then
      output_file = 'count'//particle_file(5:index(particle_file,'.rst',.true.))//'vtk' !  .true. searches backwards
    else
      write(str1,'(i0.2)') i_p
      output_file = 'count'//particle_file(5:index(particle_file,'.rst',.true.)-1)//'_'//str1//'.vtk' !  .true. searches backwards
    endif
    call write_particle_counts_to_vtk(node_list,element_list,particle_list,output_file,i_p)
    write(*,*) "Done counting, wrote output to ", output_file
  enddo
enddo

call MPI_Finalize(ierr)
contains

!> This routine counts particles in each element and writes the output to filename in VTK format
subroutine write_particle_counts_to_vtk(node_list,element_list,particle_list,filename,i_plane)
use phys_module
use data_structure
use basis_at_gaussian ! for HZ (initialise_basis must be called before use)
use mod_vtk
use mod_particles
use constants
implicit none

!> Input parameters
type(type_node_list), intent(in)      :: node_list
type(type_element_list), intent(in)   :: element_list
type (type_particle_list), intent(in) :: particle_list
character*(*), intent(in)             :: filename
integer, intent(in)                   :: i_plane !< Sum up particles near this plane, or in all if 0

integer :: nnos, i, j
real*4,allocatable    :: xyz (:,:), scalars(:,:)
real*8 :: R, R_s, R_t, Z, Z_s, Z_t
character*12, allocatable :: scalar_names(:)
integer :: i_elm

real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)

real*8  :: total_volume, total_area, volume, area, xjac, wst, phif, weight
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
do i=1,particle_list%n_particles
  if (particle_list%particle(i)%lost) cycle ! skip this iteration
  i_elm  = particle_list%particle(i)%i_elm
  if (i_elm .lt. 1) cycle
  if (i_plane .gt. 0) then
    ! project phi onto [0,2pi/n_period)*n_period/TWOPI = [0,1)
    phif = modulo(particle_list%particle(i)%x(3),TWOPI/n_period)*n_period/TWOPI
    if (i_plane .eq. 1) then
      if (phif .gt. 0.5d0/n_plane          .and. phif .lt. (1.d0-0.5d0/n_plane))    cycle
    else
      if (phif .lt. (i_plane-0.5d0)/n_plane .or. phif .gt. (i_plane+0.5d0)/n_plane) cycle
    endif
  endif
  weight = particle_list%particle(i)%weight
  scalars(i_elm,1) = scalars(i_elm,1) + real(weight,4)
enddo

! Create points for each element
do i_elm=1,element_list%n_elements
  call interp_RZ2(node_list,element_list,i_elm,0.5d0,0.5d0,R,R_s,R_t,Z,Z_s,Z_t)
  xyz(1:3,i_elm) = (/real(R,4), real(Z,4), 0.0_4/)

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
  scalars(i_elm,1) = scalars(i_elm,1) / real(volume,4)

  total_area = total_area + area
  total_volume = total_volume + volume
enddo  ! n_elements

write(*,*) "Area: ", total_area
write(*,*) "Volume: ", total_volume

call write_vtk(filename,xyz,scalar_names=scalar_names,scalars=scalars)

end subroutine write_particle_counts_to_vtk
end program count_particles_vtk
