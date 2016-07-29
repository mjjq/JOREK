!> This module contains some routines for calculating diagnostics on particles
module mod_particle_diagnostics
contains
!> Calculate particles present in specific regions on all particles
!> Performs MPI communications to sum values, returns the value
!> corresponding to all particles on node 0, and the value for each node on this node
!> Regions are: DOMAIN_PLASMA, DOMAIN_SOL, DOMAIN_OUTER_SOL,
!> DOMAIN_UPPER_PRIVATE, DOMAIN_LOWER_PRIVATE
function particles_in_regions(node_list, element_list, particle_list)
  use data_structure
  use phys_module, only: DOMAIN_PLASMA, DOMAIN_SOL, DOMAIN_OUTER_SOL, DOMAIN_UPPER_PRIVATE,        &
      DOMAIN_LOWER_PRIVATE, xpoint, xcase
  use mod_particles
  use domains
  use_mpi
  implicit none

  type(type_node_list), intent(in)     :: node_list
  type(type_element_list), intent(in)  :: element_list
  type(type_particle_list), intent(in) :: particle_list

  integer, dimension(DOMAIN_PLASMA:DOMAIN_LOWER_PRIVATE) :: particles_in_regions, tmp
  integer :: i, ifail, my_id
  integer :: domain, i_elm_axis, i_elm_xpoint(2)
  real*8  :: psi, psi_s, psi_t, psi_st, psi_ss, psi_tt
  real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis, psi_limit
  real*8, dimension(2) :: psi_xpoint, R_xpoint, Z_xpoint, s_xpoint, t_xpoint
  type (type_particle) :: p


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
  !$omp shared(node_list, element_list, particle_list, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
  !$omp     R_axis, Z_axis, psi_axis) &
  !$omp private(domain, psi, psi_s, psi_t, psi_st, psi_ss, psi_tt, p) &
  !$omp reduction(+:tmp)
  do i=1,particle_list%n_particles
    p = particle_list%particle(i)
    if (p%i_elm .le. 0 .or. p%i_elm .gt. element_list%n_elements) cycle
    call interp(node_list, element_list, p%i_elm, 1, 1, & ! force i_harm to 1
        p%st(1), p%st(2), psi, psi_s, psi_t, psi_st, psi_ss, psi_tt)

    domain = which_domain(node_list, element_list, &
        p%x(1), p%x(2), &
        psi, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
        R_axis, Z_axis, psi_axis)

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
subroutine get_particle_flux_coordinates(node_list,element_list,particle_list,fluxcoord,mask)
  use data_structure
  use mod_particles
  implicit none
  type(type_node_list), intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  type(type_particle_list), intent(in) :: particle_list
  real*8, dimension(particle_list%n_particles), intent(out) :: fluxcoord !< List of values of psi for each particle
  logical, dimension(particle_list%n_particles), intent(out) :: mask !< Mask containing .f. if particle is lost
  real*8, dimension(1) :: P_s, P_t, P_phi
  real*8               :: R, R_s, R_t, Z, Z_s, Z_t

  integer :: i

  mask = .true.
  !$omp parallel do default(none) &
  !$omp shared(particle_list, node_list, element_list, fluxcoord, mask) &
  !$omp private(P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
  do i=1,particle_list%n_particles
    if (particle_list%particle(i)%lost .or. (particle_list%particle(i)%i_elm .lt.  1)) then
      mask(i) = .false.
    else
      call interp_PRZ(node_list, element_list, particle_list%particle(i)%i_elm, &
                    (/1/), 1, particle_list%particle(i)%st(1),particle_list%particle(i)%st(2), &
                    particle_list%particle(i)%x(3), fluxcoord(i), P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
    endif
  enddo
  !$omp end parallel do
end subroutine get_particle_flux_coordinates
end module mod_particle_diagnostics
