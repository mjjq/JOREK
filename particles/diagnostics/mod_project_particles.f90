!> Module containing routines to project a function of particles onto the JOREK
!> finite elements.
!>
!> The general routine allows the user to specify a transformation object
!> to calculate the quantity to be projected.
!> This quantity is multiplied by the particle weight behind the scenes.
!> Included transformation types
!>   - 1 (proj_one)
!>
!> These functions must subclass proj_transform in this module.
!> Other modules defining transformations:
!>   - tools/mod_radiation.f90
!>
!>
!> Samples are taken every t_s.
!> Projections are performed every t_p = n t_s. (integer n >= 1)
!> If multiple samples are taken they are averaged over the integration time.
!> 
!> Projected results can be written to vtk or hdf5 by specifying the appropriate
!> options at creation time of project_particles_base.
!> Using project_to_vtk or project_to_h5 is deprecated but still supported.
!> They support only transformation function proj_one.
module mod_project_particles
use mod_io_actions
use data_structure
use mod_particle_sim
use mod_particle_types
use mod_fields
implicit none
include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure
private
!public project_particles
public projection
public proj_f_interface, proj_one
public write_particle_distribution_to_vtk, write_particle_distribution_to_h5 !< public for testing reasons, please don't use directly
public prepare_mumps_par, sample_rhs !< public for testing reasons
public DMUMPS_STRUC

interface
  pure function proj_f_interface(sim, group, particle)
    import particle_sim, particle_base
    type(particle_sim), intent(in) :: sim
    integer, intent(in) :: group
    class(particle_base), intent(in) :: particle !< Input particle
    real*8 :: proj_f_interface !< Value to be projected
  end function proj_f_interface
end interface


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
  real*8 :: smoothing2 !< hyper-smoothing factor used for this projection
  type (DMUMPS_STRUC), private :: mumps_par !< matrix is factored by mumps and stored here
  integer, public  :: n = 1 !< Number of analysis steps per projection step
  integer, private :: i = 0 !< Current index of analysis steps

  procedure(proj_f_interface), nopass, pointer :: proj_f => proj_one !< Transformation routine

  !> Output storage (optional)
  type(vtk_grid), allocatable, private :: vtk_grid !< if allocated output to vtk
  logical, public :: to_h5 = .false. !< Output to hdf5 file

  !> Temporary storage
  real*8, dimension(:,:,:,:), allocatable :: rhs !< dim (n_vertex_max*(n_order+1),n_elements,n_tor,n_groups)
  !< right-hand side for accumulation during sampling
  !< not private because we access it in unit tests in particle_projection_spec.f90
contains
  procedure :: do => project
  procedure :: close_mumps => close_mumps
end type projection
interface projection
  module procedure new_projection
end interface projection

contains

!> Project the particle density by using transformation function 1
pure function proj_one(sim, group, particle)
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_one
  proj_one = 1.d0
end function proj_one


!> Project the particle charge by using transformation function q
!> (only valid for particles of type kinetic(_leapfrog) or gc
!>
!> TODO: normalize projection with density. This function calculates
!> the integral of q instead of the mean value.
pure function proj_q(sim, group, particle)
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_q
  select type (p => particle)
  type is (particle_kinetic)
    proj_q = real(p%q, 8)
  type is (particle_kinetic_leapfrog)
    proj_q = real(p%q, 8)
  type is (particle_gc)
    proj_q = real(p%q, 8)
  class default
    proj_q = 0.d0
  end select
end function proj_q



!> Constructor for project_particles
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_projection(node_list, element_list, smoothing, smoothing2, proj_f, n, to_h5, to_vtk, &
    nsub, filename, basename, decimal_digits, fractional_digits) result(new)
  use mpi
  type(projection) :: new
  type(type_node_list), intent(in)       :: node_list
  type(type_element_list), intent(in)    :: element_list
  real*8, intent(in), optional           :: smoothing, smoothing2 !< normal and hyper-smoothing
  procedure(proj_f_interface), optional  :: proj_f !< Function to map over particles before projection
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
  new%smoothing = 0.d0
  new%smoothing2 = 0.d0
  if (present(smoothing)) new%smoothing = smoothing
  if (present(smoothing2)) new%smoothing2 = smoothing2
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

  call prepare_mumps_par(new%node_list, new%element_list, new%mumps_par, new%smoothing, new%smoothing2)
  
  if (present(proj_f)) new%proj_f => proj_f
end function new_projection

subroutine close_mumps(this)
  class(projection), intent(inout) :: this
  this%mumps_par%JOB = -2
  call DMUMPS(this%mumps_par)
end subroutine close_mumps


subroutine project(this, sim, ev)
  use mpi
  use mod_event
  !$ use omp_lib
  class(projection), intent(inout) :: this
  type(particle_sim), intent(inout)    :: sim
  type(event), intent(inout), optional :: ev
  integer :: i_group, i_tor

  ! Sample
  call sample_rhs(this, sim)
  this%i = this%i + 1

  ! Project and save
  if (this%i .eq. this%n) then
    ! Normalize RHS with number of samples
    !$omp parallel default(none) shared(this) private(i_group, i_tor)
    do i_group=1,size(this%rhs,4)
      !$omp do
      do i_tor=1,size(this%rhs,3)
        this%rhs(:,:,i_tor,i_group) = this%rhs(:,:,i_tor,i_group) * (1.d0/real(this%n,8))
      end do
      !$omp end do
    end do
    !$omp end parallel
    call project_only(this, sim)

    if (this%to_h5) call save_to_h5(this, sim)
    if (allocated(this%vtk_grid)) call save_to_vtk(this, sim)
    this%i = 0
    this%rhs = 0.d0
  end if
end subroutine project


subroutine project_only(this, sim)
  use mpi
  use mod_event
  !$ use omp_lib
  class(projection), intent(inout) :: this
  type(particle_sim), intent(inout)    :: sim
  integer :: my_id, ierr
  integer :: i_var, i_tor, i_elm, i, j, k
  integer :: index_large_i, inode, index_ij, index
  real*8 :: t0, t1, ostart, oend, mmm(3), mmm2(3)
  real*8, dimension(:), allocatable :: my_rhs, sum_rhs

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  call cpu_time(t0)
  !$ ostart = omp_get_wtime()
  ! Safety checks
  if (.not. allocated(sim%groups)) return

  ! Preparation
  allocate(my_rhs(this%mumps_par%n), sum_rhs(this%mumps_par%n))

  do i_var=1,min(size(sim%groups),n_var) ! only project the first n_var groups
    if (.not. allocated(sim%groups(i_var)%particles)) cycle
    do i_tor=1,n_tor
      my_rhs = 0.d0
      sum_rhs = 0.d0

      ! Fill RHS of Projection matrix
      do i_elm=1,this%element_list%n_elements
        do i=1,n_vertex_max
          inode = this%element_list%element(i_elm)%vertex(i)
          do j=1,n_order+1
            index_ij = (i-1)*(n_order+1) + j
            index_large_i = this%node_list%node(inode)%index(j)  ! base index in the main matrix

            my_rhs(index_large_i) = my_rhs(index_large_i) + this%rhs(index_ij, i_elm, i_tor, i_var)
          enddo
        enddo
      enddo

      ! Gather the RHS's to the root process
      ! cannot do it directly into mumps_par%rhs because this is not allocated in every process
      call MPI_Reduce(my_rhs,sum_rhs,this%mumps_par%n, &
          MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      ! Compute the solution of Ax=b (b = RHS)
      this%mumps_par%JOB = 3
      this%mumps_par%icntl(21) = 0 ! solution is available only on host
      this%mumps_par%icntl(4)  = 3 ! print only errors
      this%mumps_par%icntl(10) = -1 ! iterative refinement of one step seems to work ok
      this%mumps_par%icntl(11) = 1 ! error analysis
      if (my_id .eq. 0) then
        this%mumps_par%rhs = sum_rhs
      end if
      call DMUMPS(this%mumps_par)

      if (my_id .eq. 0) then
        do i=1,this%node_list%n_nodes
          do k=1,n_order+1
            index = this%node_list%node(i)%index(k)
            this%node_list%node(i)%values(i_tor,k,i_var) = this%mumps_par%rhs(index)
          enddo    ! order
        enddo      ! nodes
      endif
    enddo ! i_tor
  enddo ! i_var
  deallocate(sum_rhs,my_rhs)
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


!> Add samples to the right-hand side, stored in `this` by calling this%f_proj
!> for every particle and saving the contribution multiplied by each of the basis
!> functions (poloidal and toroidal) in `this%rhs`
subroutine sample_rhs(this, sim)
  use mpi
  use mod_event
  use mod_interp, only: mode_moivre, interp_RZ
  use constants, only: PI
  use mod_basisfunctions
  !$ use omp_lib
  class(projection), intent(inout)  :: this
  type(particle_sim), intent(inout) :: sim
  integer :: n_sample !< number of groups to sample
  real*8 :: HZ(n_tor)
  real*8 :: v, R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, x(3), HH(4,4), HH_s(4,4), HH_t(4,4)
  integer :: i_group, m, i, j, i_tor
  integer :: i_elm, index_ij

  ! Safety checks
  if (.not. allocated(sim%groups)) return
  n_sample = min(size(sim%groups),n_var)  ! because we have only n_var storage for now
  if (.not. allocated(this%rhs)) then
    allocate(this%rhs(n_vertex_max*(n_order+1),this%element_list%n_elements,n_tor,n_sample))
    this%rhs = 0.d0
  end if

  !$omp parallel default(none) shared(this, sim, n_sample) &
  !$omp private(x, xjac, HH, HH_s, HH_t, R_g, R_s, R_t, Z_g, Z_s, Z_t, &
  !$omp         i_group, m, i, j, i_tor, index_ij, v, HZ)
  do i_group=1,n_sample
    if (.not. allocated(sim%groups(i_group)%particles)) cycle
    !$omp do
    do m=1,size(sim%groups(i_group)%particles,1)
      associate(particle => sim%groups(i_group)%particles(m))
      if (particle%i_elm .eq. 0) cycle

      x(1:2) = particle%st
      x(3)   = particle%x(3)

      call interp_RZ(this%node_list,this%element_list,particle%i_elm,x(1),x(2),R_g,R_s,R_t,Z_g,Z_s,Z_t)
      xjac = R_s*Z_t - R_t*Z_s
      call basisfunctions(x(1), x(2), HH, HH_s, HH_t)
      call mode_moivre(x(3), HZ)
      HZ(1) = HZ(1)*0.5d0 ! int cos^2(nx) from 0 to 2pi = pi for n > 0
      HZ(:) = HZ(:)/PI ! int 1 from 0 to 2pi = 2pi
      do i=1,n_vertex_max
        do j=1,n_order+1
          index_ij = (i-1)*(n_order+1) + j

          v = HH(i,j) * this%element_list%element(particle%i_elm)%size(i,j)
          v = v * this%proj_f(sim, i_group, particle) * particle%weight

          do i_tor=1,n_tor
            !$omp atomic
            this%rhs(index_ij,particle%i_elm,i_tor,i_group) = &
            this%rhs(index_ij,particle%i_elm,i_tor,i_group) + HZ(i_tor) * v
          enddo
        enddo
      enddo
      end associate
    enddo
    !$omp end do
  enddo
  !$omp end parallel
end subroutine sample_rhs


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

  this%extension = '.vtk'
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







!> Action for projecting all particles and writing output to a hdf5 file
subroutine save_to_h5(this, sim)
  use mpi
  use mod_event
  !$ use omp_lib
  class(projection), intent(inout)  :: this
  type(particle_sim), intent(inout)    :: sim
  integer :: i, my_id, ierr
  character(len=120) :: filename
  real*8 :: t0, t1, ostart, oend

  this%extension = '.h5'

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
!> See also [project_particles] and [mod_elt_matrix] for reference of the integration method
subroutine prepare_mumps_par(node_list, element_list, mumps_par, smoothing, smoothing2, skip_factorisation)
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
real*8, intent(in)                   :: smoothing2
logical, intent(in), optional        :: skip_factorisation

type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8     :: ELM(n_vertex_max*(n_order+1),n_vertex_max*(n_order+1))
real*8     :: wgauss2(n_gauss)
real*8, dimension(n_gauss,n_gauss) :: x_g, x_s, x_t, x_ss, x_tt, x_st, &
                                      y_g, y_s, y_t, y_ss, y_tt, y_st
real*8     :: v, v_s, v_t, v_ss, v_st, v_tt, v_x, v_y, v_xx, v_xy, v_yy
real*8     :: p, p_s, p_t, p_ss, p_st, p_tt, p_x, p_y, p_xx, p_xy, p_yy
real*8     :: wst, area, volume, xjac, xjac_x, xjac_y
integer    :: i, j, k, l, m, ilarge, index_large_i, index_large_k, inode, knode
integer    :: nz_AA, n_AA, i_elm, index_ij, index_kl
integer    :: ms, mt, my_id, ierr

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
!$omp parallel do default(none) &
!$omp shared(element_list, node_list, H, H_s, H_t, H_ss, H_st, H_tt, mumps_par, wgauss2, smoothing, smoothing2) &
!$omp private(ELM, i_elm, element, nodes, i, j, k, l, ms, mt, &
!$omp         x_g, x_s, x_t, x_ss, x_st, x_tt, &
!$omp         y_g, y_s, y_t, y_ss, y_st, y_tt, &
!$omp         v, v_s, v_t, v_ss, v_st, v_tt, v_x, v_y, v_xx, v_xy, v_yy, &
!$omp         p, p_s, p_t, p_ss, p_st, p_tt, p_x, p_y, p_xx, p_xy, p_yy, &
!$omp         wst, xjac, xjac_x, xjac_y, &
!$omp         index_ij, index_kl, ilarge, &
!$omp         inode, index_large_i, knode, index_large_k) &
!$omp reduction(+:area,volume)
do i_elm=1,element_list%n_elements
  ELM = 0.d0

  element = element_list%element(i_elm)
  do m=1,n_vertex_max
    nodes(m) = node_list%node(element%vertex(m))
  enddo

  ! Set up gauss points in this element
  x_g = 0.d0; x_s = 0.d0; x_t = 0.d0; x_ss = 0.d0; x_st = 0.d0; x_tt = 0.d0
  y_g = 0.d0; y_s = 0.d0; y_t = 0.d0; y_ss = 0.d0; y_st = 0.d0; y_tt = 0.d0
  do i=1,n_vertex_max
    do j=1,n_order+1
      do ms=1, n_gauss
        do mt=1, n_gauss
          x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
          x_s(ms,mt)  = x_s(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt)  = x_t(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H_t(i,j,ms,mt)

          x_ss(ms,mt) = x_ss(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
          x_st(ms,mt) = x_st(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_st(i,j,ms,mt)
          x_tt(ms,mt) = x_tt(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

          y_g(ms,mt)  = y_g(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H(i,j,ms,mt)
          y_s(ms,mt)  = y_s(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt)  = y_t(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H_t(i,j,ms,mt)

          y_ss(ms,mt) = y_ss(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
          y_st(ms,mt) = y_st(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_st(i,j,ms,mt)
          y_tt(ms,mt) = y_tt(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_tt(i,j,ms,mt)
        enddo
      enddo
    enddo
  enddo

  ! Perform gauss integration of LHS
  do ms=1, n_gauss
    do mt=1, n_gauss
      wst = wgauss2(ms)*wgauss2(mt)
      xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

      xjac_x  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt) - 2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)   &
              + y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt))                                               &
              + x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac

      xjac_y  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt) - 2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)   &
              + x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt))                                               &
              + y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac

      area   = area   + xjac * wst
      volume = volume + TWOPI * x_g(ms,mt) * xjac * wst

      do i=1,n_vertex_max
        do j=1,n_order+1
          index_ij = (i-1)*(n_order+1) + j

          v   = h(i,j,ms,mt)   * element%size(i,j)
          v_s = h_s(i,j,ms,mt) * element%size(i,j)
          v_t = h_t(i,j,ms,mt) * element%size(i,j)

          v_x = (  y_t(ms,mt) * v_s - y_s(ms,mt) * v_t) / xjac
          v_y = (- x_t(ms,mt) * v_s + x_s(ms,mt) * v_t) / xjac

          v_ss = h_ss(i,j,ms,mt) * element%size(i,j)
          v_tt = h_tt(i,j,ms,mt) * element%size(i,j)
          v_st = h_st(i,j,ms,mt) * element%size(i,j)

          v_xx = (v_ss * y_t(ms,mt)**2 - 2.d0*v_st * y_s(ms,mt)*y_t(ms,mt) + v_tt * y_s(ms,mt)**2  &
               + v_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &
               + v_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &
               - xjac_x * (v_s * y_t(ms,mt) - v_t * y_s(ms,mt)) / xjac**2

          v_yy = (v_ss * x_t(ms,mt)**2 - 2.d0*v_st * x_s(ms,mt)*x_t(ms,mt) + v_tt * x_s(ms,mt)**2  &
               + v_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &
               + v_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &
               - xjac_y * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) ) / xjac**2

          v_xy = (- v_ss * y_t(ms,mt)*x_t(ms,mt) - v_tt * x_s(ms,mt)*y_s(ms,mt)                    &
               + v_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                         &
               - v_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                         &
               - v_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2           &
               - xjac_x * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) )   / xjac**2

          do k=1,n_vertex_max
            do l=1,n_order+1
              index_kl = (k-1)*(n_order+1) + l

              p   = h(k,l,ms,mt)   * element%size(k,l)
              p_s = h_s(k,l,ms,mt) * element%size(k,l)
              p_t = h_t(k,l,ms,mt) * element%size(k,l)

              p_ss = h_ss(k,l,ms,mt) * element%size(k,l)
              p_tt = h_tt(k,l,ms,mt) * element%size(k,l)
              p_st = h_st(k,l,ms,mt) * element%size(k,l)

              p_x = (  y_t(ms,mt) * p_s - y_s(ms,mt) * p_t) / xjac
              p_y = (- x_t(ms,mt) * p_s + x_s(ms,mt) * p_t) / xjac

              p_xx = (p_ss * y_t(ms,mt)**2 - 2.d0*p_st * y_s(ms,mt)*y_t(ms,mt) + p_tt * y_s(ms,mt)**2  &
                   + p_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &
                   + p_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &
                   - xjac_x * (p_s * y_t(ms,mt) - p_t * y_s(ms,mt)) / xjac**2

              p_yy = (p_ss * x_t(ms,mt)**2 - 2.d0*p_st * x_s(ms,mt)*x_t(ms,mt) + p_tt * x_s(ms,mt)**2  &
                   + p_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &
                   + p_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &
                   - xjac_y * (- p_s * x_t(ms,mt) + p_t * x_s(ms,mt) ) / xjac**2

              p_xy = (- p_ss * y_t(ms,mt)*x_t(ms,mt) - p_tt * x_s(ms,mt)*y_s(ms,mt)                    &
                   + p_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                         &
                   - p_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                         &
                   - p_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2           &
                   - xjac_x * (- p_s * x_t(ms,mt) + p_t * x_s(ms,mt) )   / xjac**2


              ELM(index_ij,index_kl) = ELM(index_ij,index_kl) &
                                     + p * v * xjac * x_g(ms,mt) * wst &
                                     + smoothing * (p_x * v_x + p_y * v_y) * xjac * x_g(ms,mt) * wst &
                                     + smoothing2 * (v_xx + v_x/x_g(ms,mt) + v_yy)*(p_xx + p_x/x_g(ms,mt) + p_yy) * xjac * wst
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
mumps_par%icntl(2)  = 6 ! print diagnostics, statistics and warnings to stderr
mumps_par%icntl(4)  = 3 ! print many things
mumps_par%icntl(5)  = 0 ! assembled form
mumps_par%icntl(18) = 0 ! centralized input matrix (i.e. only on cpu 0)
mumps_par%icntl(7)  = 7 ! compute symmetric permutation (PORD or SCOTCH autoselect)
mumps_par%icntl(8)  = 8 ! scaling
mumps_par%icntl(14) = 140 ! memory relaxation parameter
mumps_par%icntl(10) = -2 ! iterative refinement
if (present(skip_factorisation) .and. skip_factorisation) then
else
  call DMUMPS(mumps_par)
endif
end subroutine prepare_mumps_par



!> Calculate the structure of the vtk file without putting in any scalars
subroutine prepare_write_particle_distribution_to_vtk(node_list,element_list,nsub,xyz,ien)
use data_structure
use mod_interp, only: interp_RZ
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
      call interp_RZ(node_list,element_list,i,s,t,R,R_s,R_t,Z,Z_s,Z_t)
      inode = inode+1
      xyz(1:3,inode) = real([R, Z, 0.d0], 4)
    enddo
  enddo

  do j=1,nsub-1
     do k=1,nsub-1
        ielm        = ielm+1
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
use mod_interp
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
real*8 :: P, P_s, P_t, P_st, P_ss, P_tt

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
!$omp private(i,j,k,l,m,inode,ivar,s,t,P, P_s, P_t, P_st, P_ss, P_tt)
do i=1,element_list%n_elements
  do j=1,n_fields
    do k=1,n_tor
      ivar = (j-1)*n_tor + k
      do l=1,nsub
        s = float(l-1)/float(nsub-1)
        do m=1,nsub
          t = float(m-1)/float(nsub-1)
          inode = (i-1)*nsub*nsub+(l-1)*nsub+m
          call interp(node_list, element_list, i, ivar, k, s, t, P, P_s, P_t, P_st, P_ss, P_tt)
          scalars(inode,ivar) = real(P,4)
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
