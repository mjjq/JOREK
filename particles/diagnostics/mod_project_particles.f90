!> Module containing routines to project a function of particles onto the JOREK
!> finite elements.
!>
!> The general routine allows the user to specify a (pure) transformation function
!> to calculate the quantity to be projected.
!> Included transformation function (and default) is 1.
!>
!> Samples are taken every t_s.
!> Projections are performed every t_p = n t_s. (integer n >= 1)
!> If multiple samples are taken they are averaged over the integration time.
!> 
!> Projected results can be written to vtk or hdf5 by specifying the appropriate
!> options at creation time of project_particles_base.
!> Using project_to_vtk or project_to_h5 is deprecated but still supported.
!> 
module mod_project_particles
use mod_io_actions
use data_structure
use mod_particle_sim
implicit none
include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure
private
public project_particles
public write_particle_distribution_to_vtk, write_particle_distribution_to_h5 !< public for testing reasons, please don't use directly
public prepare_mumps_par !< public for testing reasons
public project_to_vtk, project_to_h5 !< DEPRECATED
public DMUMPS_STRUC

type :: vtk_grid
  integer :: nsub !< Number of subdivisions to make
  real*4,allocatable, private  :: xyz (:,:) !< positions of points in vtk file
  integer,allocatable, private :: ien (:,:) !< connectivity in vtk file
end type

!> Action to project all particle distributions and save them to vtk
type, extends(io_action) :: projection
  type(type_node_list), allocatable    :: node_list !< node lists to save particle projections in
  type(type_element_list), allocatable :: element_list
  real*8 :: smoothing !< Smoothing factor used for this projection
  type (DMUMPS_STRUC), private :: mumps_par !< matrix is factored by mumps and stored here
  integer :: n = 1 !< Number of analysis steps per projection step
  !> Output storage (optional)
  type(vtk_grid), allocatable, private :: vtk_grid !< if allocated output to vtk
  logical, public :: to_h5 = .false. !< Output to hdf5 file
contains
  procedure :: do => project
end type projection
interface projection
  module procedure new_projection
end interface projection

! Here for legacy reasons, please don't use anymore
type, extends(projection) :: project_to_vtk
end type
interface project_to_vtk
  module procedure new_project_to_vtk
end interface project_to_vtk
type, extends(projection) :: project_to_h5
end type project_to_h5
interface project_to_h5
  module procedure new_project_to_h5
end interface project_to_h5

contains
!> Constructor for project_particles
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_projection(node_list, element_list, smoothing, n, to_h5, to_vtk, &
    nsub, filename, basename, decimal_digits, fractional_digits) result(new)
  use mpi
  type(projection) :: new
  type(type_node_list), intent(in)       :: node_list
  type(type_element_list), intent(in)    :: element_list
  real*8, intent(in)                     :: smoothing
  ! TODO add aggregating function
  integer, intent(in), optional          :: n !< Number of analysis steps to each projection step (1 if omitted)
  logical, intent(in), optional          :: to_h5 !< Write HDF5 output after projecting (false if omitted)
  logical, intent(in), optional          :: to_vtk !< Write vtk output after projecting (false if omitted)
  integer, intent(in), optional          :: nsub !< number of subdivisions of the finite elements
  character(len=*), intent(in), optional :: filename
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional          :: decimal_digits
  integer, intent(in), optional          :: fractional_digits
  integer :: my_id, ierr

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

  allocate(new%node_list,    source=node_list)
  allocate(new%element_list, source=element_list)
  new%smoothing = smoothing
  if (present(to_vtk)) then
    if (to_vtk) then
      allocate(new%vtk_grid)
      new%vtk_grid%nsub = 4
      if (present(nsub)) new%vtk_grid%nsub = nsub
      ! Precalculate the node positions in the vtk file and the connectivity
      if (my_id .eq. 0) call prepare_write_particle_distribution_to_vtk(new%node_list, &
          new%element_list, new%vtk_grid%nsub, new%vtk_grid%xyz, new%vtk_grid%ien)
      ! We don't set the extension here since this is dynamically set in the do
      ! action (to support both vtk and h5 output)
    end if
  end if
  if (present(to_h5)) new%to_h5 = to_h5

  new%basename = "proj"
  new%n = 1
  if (present(n)) new%n = n
  if (present(filename)) new%filename = filename
  if (present(basename)) new%basename = basename
  if (present(decimal_digits)) new%decimal_digits = decimal_digits
  if (present(fractional_digits)) new%fractional_digits = fractional_digits
  new%name = "Project"
  new%log = .true.

  call prepare_mumps_par(new%node_list, new%element_list, new%mumps_par, new%smoothing)
end function new_projection


!> Constructor for project_to_vtk
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_project_to_vtk(node_list, element_list, smoothing, nsub, filename, basename, decimal_digits, fractional_digits) result(new)
  use mpi
  type(project_to_vtk) :: new
  type(projection) :: new_base
  type(type_node_list), intent(in)       :: node_list
  type(type_element_list), intent(in)    :: element_list
  real*8, intent(in)                     :: smoothing
  integer, intent(in), optional          :: nsub !< number of subdivisions of the finite elements
  character(len=*), intent(in), optional :: filename
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional          :: decimal_digits
  integer, intent(in), optional          :: fractional_digits
  integer :: my_id, ierr
  new_base = projection(node_list, element_list, smoothing, &
      to_h5=.false., to_vtk=.true., nsub=nsub, filename=filename, &
      basename=basename, decimal_digits=decimal_digits, &
      fractional_digits=fractional_digits)
  write(*,*) "DEPRECATION WARNING: project_to_vtk has been replaced with type(projection). Please consider updating your code"
  call copy_project_h5_vtk(new, new_base)
end function new_project_to_vtk


!> Subroutine needed for backwards compatibility with new_project_to_{vtk,h5}
!> Copy over all variables one by one...
subroutine copy_project_h5_vtk(out, in)
  class(projection), intent(out) :: out
  type(projection), intent(in) :: in
  ! Projection
  out%node_list = in%node_list
  out%element_list = in%element_list
  out%MUMPS_PAR = in%MUMPS_PAR
  out%n = in%n
  out%vtk_grid = in%vtk_grid
  out%to_h5 = in%to_h5
  ! IO_action
  out%filename = in%filename
  out%basename = in%basename
  out%decimal_digits = in%decimal_digits
  out%fractional_digits = in%fractional_digits
  out%extension = in%extension
end subroutine copy_project_h5_vtk



!> Constructor for project_to_h5
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_project_to_h5(node_list, element_list, smoothing, filename, basename, decimal_digits, fractional_digits) result(new)
  use mpi
  type(project_to_h5) :: new
  type(projection) :: new_base
  type(type_node_list), intent(in)       :: node_list
  type(type_element_list), intent(in)    :: element_list
  real*8, intent(in)                     :: smoothing
  character(len=*), intent(in), optional :: filename
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional          :: decimal_digits
  integer, intent(in), optional          :: fractional_digits
  integer :: my_id, ierr
  new_base = projection(node_list, element_list, smoothing, &
      to_h5=.true., to_vtk=.false., filename=filename, &
      basename=basename, decimal_digits=decimal_digits, &
      fractional_digits=fractional_digits)
  write(*,*) "DEPRECATION WARNING: project_to_h5 has been replaced with type(projection). Please consider updating your code"
  call copy_project_h5_vtk(new, new_base)
end function new_project_to_h5


subroutine project(this, sim, ev)
  use mpi
  use mod_event
  !$ use omp_lib
  class(projection), intent(inout) :: this
  type(particle_sim), intent(inout)    :: sim
  type(event), intent(inout), optional :: ev

  call project_only(this, sim) ! TODO support n != 1

  if (this%to_h5) call save_to_h5(this, sim)
  if (allocated(this%vtk_grid)) call save_to_vtk(this, sim)
end subroutine project


subroutine project_only(this, sim)
  use mpi
  use mod_event
  !$ use omp_lib
  class(projection), intent(inout) :: this
  type(particle_sim), intent(inout)    :: sim
  integer :: i, my_id, ierr
  real*8 :: t0, t1, ostart, oend, mmm(3), mmm2(3)

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  call cpu_time(t0)
  !$ ostart = omp_get_wtime()
  ! Safety checks
  if (.not. allocated(sim%groups)) return

  do i=1,min(size(sim%groups),n_var) ! only project the first n_var groups
    if (.not. allocated(sim%groups(i)%particles)) cycle
    call project_particles(this%node_list, this%element_list, this%mumps_par, &
        sim%groups(i)%particles, i)
    ! results are saved only on mpi process 0
  end do
  ! Communicate these to all processors
  call broadcast_elements(my_id, this%element_list)
  call broadcast_nodes(my_id, this%node_list)

  call cpu_time(t1)
  !$ oend = omp_get_wtime()
  !$ mmm = mpi_minmeanmax(t1-t0)
  !$ mmm2 = mpi_minmeanmax(oend-ostart)
  if (my_id .eq. 0) then
    write(*,"(A,3g12.5)") "projection cpu time", mmm
    !$ write(*,"(A,3g12.5)") "projection wall time", mmm2
  end if
end subroutine project_only


!> Save an already-projected set to a vtk file with current parameters
subroutine save_to_vtk(this, sim)
  use mpi
  use mod_event
  !$ use omp_lib
  class(projection), intent(inout) :: this
  type(particle_sim), intent(inout)    :: sim
  integer :: i, my_id, ierr
  real*8 :: t0, t1, ostart, oend
  character(len=120) :: filename

  this%extension = 'vtk'
  if (.not. allocated(this%vtk_grid)) then
    write(*,*) "ERROR: Trying to write unprepared vtk file"
    return
  end if

  if (len_trim(this%filename) .eq. 0) then
    filename = this%get_filename(sim%time)
  else
    filename = this%filename
  end if

  call cpu_time(t0)
  !$ ostart = omp_get_wtime()
  if (my_id .eq. 0) then
    ! write only on the host
    call write_particle_distribution_to_vtk(this%node_list, this%element_list, &
      trim(filename), this%vtk_grid%nsub, min(size(sim%groups),n_var), this%vtk_grid%xyz, this%vtk_grid%ien)

    write(*,*) "Written projection to ", trim(filename)
  end if
  call cpu_time(t1)
  !$ oend = omp_get_wtime()
  if (my_id .eq. 0) then
    write(*,"(A,2g12.5)") "writing cpu time", t1-t0
    !$ write(*,"(A,2g12.5)") "writing wall time", oend-ostart
  end if
end subroutine save_to_vtk

!> Action for projecting all particles and writing output to vtk
!> DEPRECATED: use new projection type instead
subroutine project_and_save_to_vtk(this, sim, ev)
  use mod_event
  class(project_to_vtk), intent(inout) :: this
  type(particle_sim), intent(inout)    :: sim
  type(event), intent(inout), optional :: ev
  call project_only(this, sim)
  call save_to_vtk(this, sim)
end subroutine project_and_save_to_vtk


!> Action for projecting all particles and writing output to vtk
!> DEPRECATED: use new projection type instead
subroutine project_and_save_to_h5(this, sim, ev)
  use mod_event
  class(project_to_h5), intent(inout)  :: this
  type(particle_sim), intent(inout)    :: sim
  type(event), intent(inout), optional :: ev
  call project_only(this, sim)
  call save_to_h5(this, sim)
end subroutine project_and_save_to_h5


!> Action for projecting all particles and writing output to vtk
subroutine save_to_h5(this, sim)
  use mpi
  use mod_event
  !$ use omp_lib
  class(projection), intent(inout)  :: this
  type(particle_sim), intent(inout)    :: sim
  integer :: i, my_id, ierr
  character(len=120) :: filename
  real*8 :: t0, t1, ostart, oend

  this%extension = 'h5'

  if (len_trim(this%filename) .eq. 0) then
    filename = this%get_filename(sim%time)
  else
    filename = this%filename
  end if

  call cpu_time(t0)
  !$ ostart = omp_get_wtime()
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  if (my_id .eq. 0) then
    ! write only on the host
    call write_particle_distribution_to_h5(this%node_list, this%element_list, &
      trim(filename), min(size(sim%groups),n_var), sim%time)

    write(*,*) "Written projection to ", trim(filename)
  end if
  call cpu_time(t1)
  !$ oend = omp_get_wtime()
  if (my_id .eq. 0) then
    write(*,"(A,2g12.5)") "writing cpu time", t1-t0
    !$ write(*,"(A,2g12.5)") "writing wall time", oend-ostart
  end if
end subroutine save_to_h5



!> Project particles by weight onto the elements
!> Saves output in node_list%values(1)
!> The projection is done by solving
!> $$\int [p(x) u(x) + \lambda p'(x) u'(x)] dV =
!> \int \sum \delta(x_i-x) w_i u(x) dV $$
!> where p(x) is in the bernstein representation, $u(x)$ are the test functions
!> composed of two basis functions and $w_i$ is the particle weight.
!> A smoothing factor lambda is included
!> x is a vector (R,Z,phi) and dV is r dr dphi
!> divide by 1 or 2pi on both sides (LHS gets 2pi for n=0 mode, 1pi for other modes)
!>
!> See also [project_particles]
subroutine prepare_mumps_par(node_list, element_list, mumps_par, smoothing, skip_factorisation)
use phys_module
use data_structure
use basis_at_gaussian
use mod_basisfunctions
use mpi
implicit none

type (type_node_list), intent(in)    :: node_list !< A copy of the node list which will be used to save variables
type (type_element_list), intent(in) :: element_list
type (DMUMPS_STRUC), intent(inout)   :: mumps_par
real*8, intent(in)                   :: smoothing
logical, intent(in), optional        :: skip_factorisation

type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8     :: ELM(n_vertex_max*(n_order+1),n_vertex_max*(n_order+1)), RHS(n_vertex_max*(n_order+1),element_list%n_elements)
real*8     :: wgauss2(n_gauss)
real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)
real*8     :: R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, x(3)
real*8     :: v, v_x, v_y, psi, psi_x, psi_y, wst, area, volume
integer    :: i, j, k, l, m, i_tor, ilarge, index_large_i, index_large_k, inode, knode
integer    :: nz_AA, n_AA, i_elm, index_ij, index_kl
integer    :: ms, mt, n_p, my_id, ierr

! Initialise MUMPS
mumps_par%COMM = MPI_COMM_WORLD
mumps_par%JOB  = -1
mumps_par%SYM  = 0
mumps_par%PAR  = 1
call DMUMPS(mumps_par)
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

nz_AA = element_list%n_elements * (n_vertex_max * (n_order+1))**2
n_AA = maxval(node_list%node(1:node_list%n_nodes)%index(4))

! Only perform the construction of the matrix on the host
if (my_id .eq. 0) then
write(*,*) ' number of unknowns      : ',n_AA, node_list%n_nodes * (n_order+1)
write(*,*) ' nz_AA                   : ',nz_AA

! Allocate space for elements
allocate(mumps_par%A(nz_AA),mumps_par%irn(nz_AA),mumps_par%jcn(nz_AA))
allocate(mumps_par%rhs(n_AA))
mumps_par%irn = 0
mumps_par%jcn = 0
mumps_par%A   = 0.d0
mumps_par%RHS = 0.d0

! Copy wgauss into wgauss2 to get around gfortran not recognizing it as a shared
! thing https://groups.google.com/forum/#!topic/comp.lang.fortran/VKhoAm8m9KE
wgauss2 = wgauss

area = 0.
volume = 0.

write(*,*) '*******************************************'
write(*,*) '* constructing particle projection matrix *'
write(*,*) '*******************************************'
!$omp parallel do default(none) & ! instead of none, bugfix for gfortran: 
!$omp shared(element_list, node_list, H, H_s, H_t, mumps_par, wgauss2, smoothing) &
!$omp private(ELM, i_elm, element, nodes, i, j, k, l, ms, mt, &
!$omp         x_g, y_g, x_s, x_t, y_s, y_t, wst, xjac, v, v_x, v_y, &
!$omp         index_ij, index_kl, psi, psi_x, psi_y, ilarge, &
!$omp         inode, index_large_i, knode, index_large_k) &
!$omp reduction(+:area,volume)
do i_elm=1,element_list%n_elements
  ELM = 0.d0

  element = element_list%element(i_elm)
  do m=1,n_vertex_max
    nodes(m) = node_list%node(element%vertex(m))
  enddo

  ! Set up gauss points in this element
  x_g = 0.d0; x_s = 0.d0; x_t = 0.d0; y_g = 0.d0; y_s = 0.d0; y_t = 0.d0
  do i=1,n_vertex_max
    do j=1,n_order+1
      do ms=1, n_gauss
        do mt=1, n_gauss
          x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
          y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H(i,j,ms,mt)

          x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_t(i,j,ms,mt)

          y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_t(i,j,ms,mt)
        enddo
      enddo
    enddo
  enddo

  ! Perform gauss integration of LHS
  do ms=1, n_gauss
    do mt=1, n_gauss
      wst = wgauss2(ms)*wgauss2(mt)
      xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

      area   = area   + xjac * wst
      volume = volume + TWOPI * x_g(ms,mt) * xjac * wst

      do i=1,n_vertex_max
        do j=1,n_order+1
          index_ij = (i-1)*(n_order+1) + j

          v   = h(i,j,ms,mt)  * element%size(i,j)
          v_x = (  y_t(ms,mt) * h_s(i,j,ms,mt) - y_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac
          v_y = (- x_t(ms,mt) * h_s(i,j,ms,mt) + x_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac
          do k=1,n_vertex_max
            do l=1,n_order+1
              index_kl = (k-1)*(n_order+1) + l

              psi   = h(k,l,ms,mt) * element%size(k,l)
              psi_x = (   y_t(ms,mt) * h_s(k,l,ms,mt) - y_s(ms,mt) * h_t(k,l,ms,mt) ) * element%size(k,l) / xjac
              psi_y = ( - x_t(ms,mt) * h_s(k,l,ms,mt) + x_s(ms,mt) * h_t(k,l,ms,mt) ) * element%size(k,l) / xjac

              ELM(index_ij,index_kl) = ELM(index_ij,index_kl) + psi * v * xjac * x_g(ms,mt) * wst &
                                     + smoothing * (psi_x * v_x + psi_y * v_y) * xjac * x_g(ms,mt) * wst
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  ! Save contribution of this element in MUMPS format
  do i=1,n_vertex_max
    inode = element_list%element(i_elm)%vertex(i)
    do j=1,n_order+1
      index_ij = (i-1)*(n_order+1) + j
      index_large_i = node_list%node(inode)%index(j)  ! base index in the main matrix

      do k=1,n_vertex_max
        knode = element_list%element(i_elm)%vertex(k)
        do l=1,n_order+1
          index_kl = (k-1)*(n_order+1) + l

          index_large_k = node_list%node(knode)%index(l)  ! base index in the main matrix

          ! Explicitly calculate the index
          ilarge = l + (k-1)*(n_order+1) + (j-1)*(n_vertex_max*(n_order+1)) + (i-1)*(n_vertex_max*(n_order+1)**2) &
                     + (i_elm-1)*(n_vertex_max*(n_order+1))**2

          mumps_par%irn(ilarge) = index_large_i
          mumps_par%jcn(ilarge) = index_large_k
          mumps_par%A(ilarge)   = ELM(index_ij,index_kl)
        enddo
      enddo
    enddo
  enddo
enddo
!$omp end parallel do

write(*,'(A,e14.6)') ' Area        : ',area
write(*,'(A,e14.6)') ' Volume      : ',volume
end if

! Perform the analysis and factorisation with all nodes
mumps_par%JOB       = 4
mumps_par%n         = n_AA
mumps_par%nz        = nz_AA
mumps_par%icntl(5)  = 0 ! assembled form
mumps_par%icntl(18) = 0 ! centralized input matrix (i.e. only on cpu 0)
mumps_par%icntl(7)  = 7 ! compute symmetric permutation (PORD or SCOTCH autoselect)
mumps_par%icntl(8)  = 7 ! scaling
mumps_par%icntl(14) = 80 ! memory relaxation parameter
mumps_par%icntl(4)  = 1 ! 2=Print errors, warnings and main statistics
if (present(skip_factorisation) .and. skip_factorisation) then
else
  call DMUMPS(mumps_par)
endif
end subroutine prepare_mumps_par


!> Perform the actual projection of a set of particles on variable ivar_out in node_list.
subroutine project_particles(node_list, element_list, mumps_par, particles, ivar_out, skip_proj)
use phys_module
use data_structure
use mod_basisfunctions
use mpi
use mod_interp_PRZ
!$ use omp_lib
implicit none

type (type_node_list), intent(inout) :: node_list !< A copy of the node list which will be used to save variables
type (type_element_list), intent(in) :: element_list
type (DMUMPS_STRUC), intent(inout)   :: mumps_par
class (particle_base), intent(in), dimension(:)    :: particles
integer, intent(in) :: ivar_out
logical, optional :: skip_proj !< Do not project but write the RHS to the nodes

real*8, allocatable :: RHS(:,:), sum_rhs(:), my_rhs(:)
real*8     :: v, R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, x(3), HH(4,4), HH_s(4,4), HH_t(4,4)
integer    :: i, j, k, m, i_tor, index_large_i, inode, n_AA
integer    :: i_elm, index, index_ij, my_id, ierr
real*8, dimension(0) :: P, P_s, P_t, P_phi

call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

allocate(RHS(n_vertex_max*(n_order+1),element_list%n_elements),my_rhs(mumps_par%n), sum_rhs(mumps_par%n))
! Create the RHS and calculate the projection
! faster if we use mumps with multiple right-hand sides? TODO
do i_tor=1, n_tor
  RHS = 0.d0
  my_rhs = 0.d0
  !$omp parallel do default(none) &
  !$omp shared(particles, element_list, node_list, mode, i_tor, RHS) &
  !$omp private(x, xjac, HH, HH_s, HH_t, R_g, R_s, R_t, Z_g, Z_s, Z_t, &
  !$omp         i, j, index_ij, v, m, i_elm, index_large_i, inode)
  do m=1,size(particles,1)
    if (particles(m)%i_elm .eq. 0) cycle

    associate(particle => particles(m))
    x(1:2) = particle%st
    x(3)   = particle%x(3)

    call interp_RZ(node_list,element_list,particle%i_elm,x(1),x(2),R_g,R_s,R_t,Z_g,Z_s,Z_t)
    xjac = R_s*Z_t - R_t*Z_s
    call basisfunctions3(x(1), x(2), HH, HH_s, HH_t)
    do i=1,n_vertex_max
      do j=1,n_order+1
        index_ij = (i-1)*(n_order+1) + j

        v   = HH(i,j) * element_list%element(particle%i_elm)%size(i,j)
        if (mode(i_tor) .gt. 1) then ! mode(1) = 0, mode(2) -> cos, mode(3) -> sin
          if (mod(i_tor,2) .eq. 0) then
            v = v * cos(mode(i_tor)*x(3))
          else
            v = v * sin(mode(i_tor)*x(3))
          endif
          v = v / PI ! int cos^2(nx) from 0 to 2pi = pi for n > 0
        else
          v = v / TWOPI ! int 1 from 0 to 2pi = 2pi
        endif

        !$omp atomic
        RHS(index_ij,particle%i_elm) = RHS(index_ij,particle%i_elm) + v * particle%weight
      enddo
    enddo
    end associate
  enddo
  !$omp end parallel do

  ! Fill RHS of Projection matrix
  do i_elm=1,element_list%n_elements
    do i=1,n_vertex_max
      inode = element_list%element(i_elm)%vertex(i)
      do j=1,n_order+1
        index_ij = (i-1)*(n_order+1) + j
        index_large_i = node_list%node(inode)%index(j)  ! base index in the main matrix

        my_rhs(index_large_i) = my_rhs(index_large_i) + RHS(index_ij, i_elm)
      enddo
    enddo
  enddo

  ! Gather the RHS's to the root process
  ! cannot do it directly into mumps_par%rhs because this is not allocated in every process
  call MPI_Reduce(my_rhs,sum_rhs,mumps_par%n, &
      MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  ! Compute the solution of Ax=b (b = RHS)
  mumps_par%JOB = 3
  mumps_par%icntl(21) = 0 ! solution is available only on host
  mumps_par%icntl(4)  = 1 ! print only errors
  if (my_id .eq. 0) then
    mumps_par%rhs = sum_rhs
  end if
  if (present(skip_proj) .and. skip_proj) then
  else
    call DMUMPS(mumps_par)
  end if

  if (my_id .eq. 0) then
    do i=1,node_list%n_nodes
      do k=1,n_order+1
        index = node_list%node(i)%index(k)
        node_list%node(i)%values(i_tor,k,ivar_out) = mumps_par%rhs(index)
      enddo    ! order
    enddo      ! nodes
  end if
enddo ! i_tor
deallocate(RHS,sum_rhs,my_rhs)
end subroutine project_particles


!> Calculate the structure of the vtk file without putting in any scalars
subroutine prepare_write_particle_distribution_to_vtk(node_list,element_list,nsub,xyz,ien)
use data_structure
implicit none

!> Input parameters
type(type_node_list), intent(in)    :: node_list
type(type_element_list), intent(in) :: element_list
integer, intent(in)                 :: nsub !< Number of subdivisions of each element
real*4,allocatable, intent(out)     :: xyz (:,:)
integer,allocatable, intent(out)    :: ien (:,:)

integer :: nnos, nnoel, nel, i, j, ielm, inode, k
real*8 :: s, t, R, R_s, R_t, Z, Z_s, Z_t

nnos = nsub*nsub*node_list%n_nodes
allocate(xyz(3,nnos))

nnoel = 4
nel   = (nsub-1)*(nsub-1)*element_list%n_elements
allocate(ien(nnoel,nel))

inode   = 0
ielm    = 0
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
      xyz(1:3,inode) = real([R, Z, 0.d0], 4)
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
end subroutine prepare_write_particle_distribution_to_vtk



!> Save a particle distribution in a sort of minimized JOREK restart format.
!> This contains only: n_tor, n_period, n_fields->n_var, n_vertex_max, vertex, x, size, n_elements, values
!> neighbours, RCS_version
!> model=-1 to signify that the variables have not the expected meaning
subroutine write_particle_distribution_to_h5(node_list,element_list,filename,n_fields,time)
use data_structure
use hdf5
use hdf5_io_module
use tr_module
!$ use omp_lib
implicit none

!> Input parameters
type(type_node_list), intent(in)    :: node_list
type(type_element_list), intent(in) :: element_list
character*(*), intent(in)           :: filename
integer, intent(in) :: n_fields !< number of different particle groups to output
real*8, intent(in) :: time

 
#include "version.h"
! --- Local variables
integer :: i
character*50             :: version_control

integer(HID_T)     :: file_id
integer            :: ind, ierr

! type_node, node_list%n_nodes
real(RKIND), allocatable :: t_x(:,:,:)                   ! n_order+1, n_dim
real(RKIND), allocatable :: t_values(:,:,:,:)            ! n_tor, n_order+1, n_fields

! element, element_list%n_elements
integer,     allocatable :: t_vertex(:,:)                ! n_vertex_max
integer,     allocatable :: t_neighbours(:,:)            ! n_vertex_max
real(RKIND), allocatable :: t_size(:,:,:)                ! n_vertex_max,n_order+1

! type_node, node_list%n_nodes
call tr_allocate(t_x,1,node_list%n_nodes,1,n_order+1,1,n_dim, &
     "node_list%x",CAT_UNKNOWN)
call tr_allocate(t_values,1,node_list%n_nodes,1,n_tor,1,n_order+1,1,n_fields, &
     "node_list%values",CAT_UNKNOWN)

! element_list%n_elements
call tr_allocate(t_vertex,1,element_list%n_elements,1,n_vertex_max,"vertex",CAT_UNKNOWN)
call tr_allocate(t_neighbours,1,element_list%n_elements,1,n_vertex_max,"neighbours",CAT_UNKNOWN)
call tr_allocate(t_size,1,element_list%n_elements,1,n_vertex_max,1,n_order+1,"size",CAT_UNKNOWN)

do i=1,node_list%n_nodes
   t_x(i,:,:)        = node_list%node(i)%x
   t_values(i,:,:,:) = node_list%node(i)%values(:,:,1:n_fields)
end do

do i=1,element_list%n_elements
  t_vertex(i,:)       = element_list%element(i)%vertex
  t_neighbours(i,:)   = element_list%element(i)%neighbours
  t_size(i,:,:)       = element_list%element(i)%size
end do

! -> Create and open HDF5 file
call HDF5_create(trim(filename),file_id,ierr)
if (ierr.ne.0) then
  write(*,*) ' ==> error for opening of HDF5 file',filename
end if
  
! -> Save version of revision control system
write(version_control,'(A)') trim(adjustl(RCS_VERSION))
version_control = trim(adjustl(version_control))
call HDF5_char_saving(file_id,version_control,"RCS_version"//char(0))

! -> Save parameters
call HDF5_integer_saving(file_id,-1,'jorek_model'//char(0)) ! Indicate that this is not a normal JOREK model
call HDF5_integer_saving(file_id,n_fields,'n_var'//char(0))
call HDF5_integer_saving(file_id,n_dim,'n_dim'//char(0))
call HDF5_integer_saving(file_id,n_order,'n_order'//char(0))
call HDF5_integer_saving(file_id,n_tor,'n_tor'//char(0))
call HDF5_integer_saving(file_id,n_period,'n_period'//char(0))
call HDF5_integer_saving(file_id,n_vertex_max,'n_vertex_max'//char(0))
call HDF5_integer_saving(file_id,n_nodes_max,'n_nodes_max'//char(0))
call HDF5_integer_saving(file_id,n_elements_max,'n_elements_max'//char(0))

! -> 
call HDF5_integer_saving(file_id,node_list%n_nodes,'n_nodes'//char(0))
call HDF5_integer_saving(file_id,element_list%n_elements,'n_elements'//char(0))
call HDF5_integer_saving(file_id,node_list%n_dof,'n_dof'//char(0))

call HDF5_array3D_saving(file_id,t_x, &
     node_list%n_nodes,n_order+1,n_dim,'x'//char(0))
call HDF5_array4D_saving(file_id,t_values, &
     node_list%n_nodes,n_tor,n_order+1,n_fields,'values'//char(0))

call HDF5_array2D_saving_int(file_id,t_vertex, &
     element_list%n_elements,n_vertex_max,'vertex'//char(0))
call HDF5_array2D_saving_int(file_id,t_neighbours, &
     element_list%n_elements,n_vertex_max,'neighbours'//char(0))
call HDF5_array3D_saving(file_id,t_size, &
     element_list%n_elements,n_vertex_max,n_order+1,'size'//char(0))
call HDF5_real_saving(file_id,time,'t_now'//char(0))

! -> close file
call HDF5_close(file_id)
end subroutine write_particle_distribution_to_h5


!> Helper subroutine to write a particle distribution, saved in variables 1:n,
!> to `filename`. `nsub` is the number of subdivisions to make per element.
!> Should be called only by process 1
subroutine write_particle_distribution_to_vtk(node_list,element_list,filename,nsub,n_fields,xyz,ien)
use data_structure
use phys_module, only: mode
use mod_vtk
use mod_interp4
use mod_basisfunctions
!$ use omp_lib
implicit none

!> Input parameters
type(type_node_list), intent(in)    :: node_list
type(type_element_list), intent(in) :: element_list
character*(*), intent(in)           :: filename
integer, intent(in) :: nsub !< Number of subdivisions of each element
integer, intent(in) :: n_fields !< number of different particle groups to output
real*4, intent(in)  :: xyz(:,:)
integer, intent(in) :: ien(:,:)

integer :: nnos, i, j, k, l, m, inode, ivar
real*4, allocatable :: scalars(:,:), vectors(:,:,:)
integer :: n_scalars, n_vectors = 0
character*12, allocatable :: vector_names(:), scalar_names(:)
real*8 :: s, t

integer :: i_t, kv, iv, kf
integer, parameter :: etype = 9 ! for vtk_quad

n_scalars = n_tor * n_fields
nnos = nsub*nsub*node_list%n_nodes
allocate(scalars(nnos,n_scalars),vectors(nnos,3,n_vectors))
allocate(scalar_names(n_scalars),vector_names(n_vectors))
do i=1,n_fields
  write(scalar_names(n_tor*(i-1)+1),'(A,i0.2)') "rho_", i
  do j=1,(n_tor-1)/2
    write(scalar_names(n_tor*(i-1)+j+1),"(A4,i0.2,A4,i0.2)") "rho_", i, "_cos", mode(2*j)
    write(scalar_names(n_tor*(i-1)+j+2),"(A4,i0.2,A4,i0.2)") "rho_", i, "_sin", mode(2*j+1)
  end do
end do

scalars = 0.d0
vectors = 0.d0

! Create points for each element
!$omp parallel do default(none) shared(element_list,nsub,node_list,n_fields,scalars) &
!$omp private(i,j,k,l,m,inode,ivar,s,t)
do i=1,element_list%n_elements
  do j=1,n_fields
    do k=1,n_tor
      ivar = (j-1)*n_tor + k
      do l=1,nsub
        s = float(l-1)/float(nsub-1)
        do m=1,nsub
          t = float(m-1)/float(nsub-1)
          inode = (i-1)*nsub*nsub+(l-1)*nsub+m
          scalars(inode,ivar) = real(interp5(node_list,element_list,i,j,k,s,t),4)
        end do
      end do
    end do
  end do
end do ! n_elements
!$omp end parallel do

! ------------- Write to VTK
call write_vtk(filename,xyz,&
  ien, etype,&
  scalar_names,scalars,&
  vector_names,vectors)
end subroutine write_particle_distribution_to_vtk
end module mod_project_particles
