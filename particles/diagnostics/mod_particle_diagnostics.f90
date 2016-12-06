!> This module contains some routines for calculating diagnostics on particles
module mod_particle_diagnostics
use mod_io_actions
use data_structure
use mod_particle_sim
use mod_fields
implicit none

!> Action to calculate pphi_H_mu and write this to an HDF5 file
!> in an extensible (in the time-dimension) dataset
type, extends(io_action) :: write_constants_of_motion
  class(read_fields_base), pointer :: fields
  ! TODO hdf5 handle
contains
  procedure :: do => do_write_constants_of_motion
end type write_constants_of_motion
interface write_constants_of_motion
  module procedure new_write_constants_of_motion
end interface write_constants_of_motion


contains
!> Constructor. Must use this or open the HDF5 file manually.
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_write_constants_of_motion(filename, fields) result(new)
  type(write_constants_of_motion) :: new
  character(len=*), intent(in)    :: filename
  class(read_fields_base), target :: fields
  new%filename = filename
  new%name = "WriteConstantsOfMotion"
  new%log = .true.
  new%fields => fields
end function new_write_constants_of_motion

!> Action to calculate all of these values and write them to an HDF5 file
subroutine do_write_constants_of_motion(this, sim)
  use mpi
  !$ use omp_lib
  class(write_constants_of_motion), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  integer :: i, my_id, ierr
  real*8 :: t0, t1, ostart, oend

  ! Safety checks
  if (.not. allocated(sim%groups)) return
  ! TODO test if number of groups and dimension of hdf5 is same
  ! TODO test the number of particles dimension
  ! TODO test if handle is open, give error about using constructor
  ! TODO warn if there is already data for this timestep

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  ! Calculate stuff
  ! Write it out in parallel
end subroutine do_write_constants_of_motion

subroutine calculate_pphi_H_mu(node_list, element_list, particles, mass, out, mask)
  use data_structure
  use mod_particle_types
  use phys_module, only: F0
  use constants
  type(type_node_list), intent(in)                    :: node_list
  type(type_element_list), intent(in)                 :: element_list
  class(particle_base), intent(in), dimension(:)      :: particles
  real*8, intent(in)                                  :: mass
  real*8, dimension(3,size(particles,1)), intent(out) :: out !< List of values
  logical, dimension(size(particles,1)), intent(out)  :: mask !< Mask containing .f. if particle is lost
  real*8, dimension(1) :: P, P_s, P_t, P_phi, inv_st_jac, psi_R, psi_Z, B(3), B_hat(3), B_norm
  real*8               :: R, R_s, R_t, Z, Z_s, Z_t

  integer :: i

  mask = .true.
  out  = 0.d0
  !$omp parallel do default(none) &
  !$omp shared(particles, node_list, element_list, out, mask) &
  !$omp private(P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t, &
  !$omp inv_st_jac, psi_R, psi_Z, B, B_hat)
  do i=1,size(particles,1)
    if (particles(i)%i_elm .lt. 1) then
      mask(i) = .false.
    else
      select type (pa => particles(i))
      type is (particle_kinetic_leapfrog)
        ! Calculate psi and B
        call interp_PRZ(node_list, element_list, pa%i_elm, &
                      [1], 1, pa%st(1),pa%st(2), &
                      pa%x(3), P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
        inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
        psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
        psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
        ! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
        B        = [+psi_Z, -psi_R, F0] / R
        B_hat = B/norm2(B)

        ! Calculate output variables
        out(1,i) = real(pa%q,8) * EL_CHG * P(1) + mass * ATOMIC_MASS_UNIT * R * pa%v(3)
        out(2,i) = mass * ATOMIC_MASS_UNIT * 0.5d0 * dot_product(pa%v,pa%v)
        out(3,i) = mass * ATOMIC_MASS_UNIT * 0.5d0 * dot_product(&
            pa%v - dot_product(pa%v,B_hat)*B_hat, &
            pa%v - dot_product(pa%v,B_hat)*B_hat &
            )/norm2(B)
      class default
        write(*,*) "ERROR: calculate_pphi_H_mu not implemented for this particle type"
      end select
    endif
  enddo
  !$omp end parallel do
end subroutine calculate_pphi_H_mu



!> Calculate particles present in specific regions on all particles
!> Performs MPI communications to sum values, returns the value
!> corresponding to all particles on node 0, and the value for each node on this node
!> Regions are: DOMAIN_PLASMA, DOMAIN_SOL, DOMAIN_OUTER_SOL,
!> DOMAIN_UPPER_PRIVATE, DOMAIN_LOWER_PRIVATE
function particles_in_regions(node_list, element_list, particles)
  use data_structure
  use phys_module, only: DOMAIN_PLASMA, DOMAIN_SOL, DOMAIN_OUTER_SOL, DOMAIN_UPPER_PRIVATE,        &
      DOMAIN_LOWER_PRIVATE, xpoint, xcase
  use mod_particle_types
  use domains
  use mpi
  implicit none

  type(type_node_list), intent(in)     :: node_list
  type(type_element_list), intent(in)  :: element_list
  class(particle_base), intent(in), dimension(:) :: particles

  integer, dimension(DOMAIN_PLASMA:DOMAIN_LOWER_PRIVATE) :: particles_in_regions, tmp
  integer :: i, ifail, my_id
  integer :: domain, i_elm_axis, i_elm_xpoint(2)
  real*8  :: psi, psi_s, psi_t, psi_st, psi_ss, psi_tt
  real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis, psi_limit
  real*8, dimension(2) :: psi_xpoint, R_xpoint, Z_xpoint, s_xpoint, t_xpoint


  !! Preparation (force my_id to 1 to suppress message)
  call find_axis(1,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

  if (xpoint) then
    call find_xpoint(1,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
    psi_limit  = psi_xpoint(1)
    if( (xcase .eq. 2) .or. ((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
      psi_limit = psi_xpoint(2)
    endif
  else
    psi_limit = 0.d0
  endif

  ! Call which_domain once to setup saved values
  domain = which_domain(node_list, element_list, &
      0.d0, 0.d0, &
      0.d0, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
      R_axis, Z_axis, psi_axis)

  tmp = 0
  !$omp parallel do default(none) &
  !$omp shared(node_list, element_list, particles, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
  !$omp     R_axis, Z_axis, psi_axis) &
  !$omp private(domain, psi, psi_s, psi_t, psi_st, psi_ss, psi_tt, p) &
  !$omp reduction(+:tmp)
  do i=1,size(particles,1)
    associate (p => particles(i))
    if (p%i_elm .le. 0 .or. p%i_elm .gt. element_list%n_elements) cycle
    call interp(node_list, element_list, p%i_elm, 1, 1, & ! force i_harm to 1
        p%st(1), p%st(2), psi, psi_s, psi_t, psi_st, psi_ss, psi_tt)

    domain = which_domain(node_list, element_list, &
        p%x(1), p%x(2), &
        psi, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
        R_axis, Z_axis, psi_axis)
    end associate

    tmp(domain) = tmp(domain) + 1
  end do
  !$omp end parallel do

  ! Save values on nodes
  particles_in_regions = tmp
  ! Mpi communication to get the total answer on node 0
  call MPI_Reduce(tmp, particles_in_regions, DOMAIN_LOWER_PRIVATE-DOMAIN_PLASMA, &
    MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ifail)
  call MPI_Comm_Rank(MPI_COMM_WORLD, my_id, ifail)
end function particles_in_regions


!> Calculate flux coordinates of particles
subroutine get_particle_flux_coordinates(node_list,element_list,particles,fluxcoord,mask)
  use data_structure
  use mod_particle_types
  implicit none
  type(type_node_list), intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  class(particle_base), intent(in), dimension(:) :: particles
  real*8, dimension(size(particles,1)), intent(out)  :: fluxcoord !< List of values of psi for each particle
  logical, dimension(size(particles,1)), intent(out) :: mask !< Mask containing .f. if particle is lost
  real*8, dimension(1) :: P_s, P_t, P_phi
  real*8               :: R, R_s, R_t, Z, Z_s, Z_t

  integer :: i

  mask = .true.
  !$omp parallel do default(none) &
  !$omp shared(particles, node_list, element_list, fluxcoord, mask) &
  !$omp private(P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
  do i=1,size(particles,1)
    associate(p => particles(i))
    if (p%i_elm .lt.  1) then
      mask(i) = .false.
    else
      call interp_PRZ(node_list, element_list, p%i_elm, &
                    (/1/), 1, p%st(1),p%st(2), &
                    p%x(3), fluxcoord(i), P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
    endif
    end associate
  enddo
  !$omp end parallel do
end subroutine get_particle_flux_coordinates
end module mod_particle_diagnostics
