!> Project particles onto the elements in JOREK and export to a vtk file
!! It uses the elements of file jorek_restart.rst
!! Run as project_particles_vtk particle_file.rst < jorek_in
!! It reads some parameters from the vtk.nml namelist in the current directory
program project_particles_vtk
use data_structure
use phys_module
use basis_at_gaussian
use nodes_elements
use mod_particles
use mod_project_particles
use mod_import_export_particles
use mpi_mod
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
nsub                   = 2	 ! Number of subdivisions of the cubic finite elements into linear pieces
i_tor                  = -1	 ! If i_tor > 0, only this mode will be included in the vtk file...
i_plane                = 1	 ! ... otherwise, all modes will be summed up at the toroidal plane i_plane
without_n0_mode        = .false. ! If true, do not include the n=0 mode (i_tor=1)

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
write(*,*) 'i_tor           =', i_tor
write(*,*) 'i_plane         =', i_plane
write(*,*) 'without_n0_mode =', without_n0_mode
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

! Get particle filename from commandline
call get_command_argument(1, particle_file)
call import_particles(particle_file, particle_list)

! Project onto element_list (saves into the first value in node_list!)
call project_particles(node_list, element_list, particle_list)

! Write the first value of node_list to a vtk file
particle_file = particle_file(1:index(particle_file,'.rst',.true.))//'vtk' !  .true. searches backwards
write(*,*) "Done projecting, writing output to ", particle_file
call write_particle_distribution_to_vtk(node_list,element_list,particle_file,nsub,i_tor,i_plane,without_n0_mode)


call MPI_FINALIZE(IERR)
contains
subroutine write_particle_distribution_to_vtk(node_list,element_list,filename,nsub,i_tor,i_plane,without_n0_mode)
use data_structure
use basis_at_gaussian ! for HZ (initialise_basis must be called before)
use mod_vtk
implicit none

!> Input parameters
type(type_node_list), intent(in)    :: node_list
type(type_element_list), intent(in) :: element_list
character*(*), intent(in)           :: filename
integer, intent(in) :: nsub, i_tor, i_plane
logical, intent(in) :: without_n0_mode

integer :: nnos, nnoel, nel, i, j, ielm, inode, k
real*4,allocatable    :: xyz (:,:), scalars(:,:), vectors(:,:,:)
real*8 :: s, t, R, R_s, R_t, Z, Z_s, Z_t
real*8 :: P, P_s, P_t, P_st, P_ss, P_tt
integer,allocatable   :: ien (:,:)
integer, parameter :: n_scalars = 1, n_vectors = 0
character*12, allocatable :: vector_names(:), scalar_names(:)

integer :: i_min, i_max, i_t
integer, parameter :: etype = 9 ! for vtk_quad

nnos = nsub*nsub*node_list%n_nodes
allocate(xyz(3,nnos),scalars(nnos,n_scalars),vectors(nnos,3,n_vectors))
allocate(scalar_names(n_scalars),vector_names(n_vectors))
scalar_names(1) = "num_density"

nnoel = 4
nel   = (nsub-1)*(nsub-1)*element_list%n_elements
allocate(ien(nnoel,nel))

inode   = 0
ielm    = 0
scalars = 0.d0
vectors = 0.d0
xyz     = 0
ien     = 0

! Create points for each element
do i=1,element_list%n_elements
  do j=1,nsub
    s = float(j-1)/float(nsub-1)
    ! Create nsub^2 points per element at regularly spaced intervals
    do k=1,nsub
      t = float(k-1)/float(nsub-1)
      call interp_RZ2(node_list,element_list,i,s,t,R,R_s,R_t,Z,Z_s,Z_t)
      inode = inode+1
      xyz(1:3,inode) = (/ R, Z, 0.d0/)

      if (i_tor .gt. -1) then 
        i_min = i_tor
        i_max = i_tor
      else
        i_min = 1
        i_max = n_tor
      endif

      do i_t = i_min, i_max
        if ( ( i_t == 1 ) .and. ( without_n0_mode ) ) cycle ! Do not include the n=0 mode
        call interp(node_list,element_list,i,1,i_t,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
        scalars(inode,1) = scalars(inode,1) + P * HZ(i_t,i_plane)
      enddo
    enddo
  enddo

  do j=1,nsub-1
     do k=1,nsub-1
        ielm	  = ielm+1
        ien(1,ielm) = inode - nsub*nsub + nsub*(j-1) + k-1       ! 0 based indices for VTK
        ien(2,ielm) = inode - nsub*nsub + nsub*(j  ) + k-1
        ien(3,ielm) = inode - nsub*nsub + nsub*(j  ) + k
        ien(4,ielm) = inode - nsub*nsub + nsub*(j-1) + k
     enddo
  enddo
enddo  ! n_elements

! ------------- Write to VTK
call write_vtk(filename,xyz,&
  ien, etype,&
  scalar_names,scalars,&
  vector_names,vectors)
end subroutine write_particle_distribution_to_vtk
end program project_particles_vtk
