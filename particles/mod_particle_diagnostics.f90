!> This module contains some routines for calculating diagnostics on particles
module mod_particle_diagnostics
contains
!> Calculate particles present in specific regions on all particles
!! Performs MPI communications to sum values, returns the value
!! corresponding to all particles on node 0, and the value for each node on this node
!! Regions are: DOMAIN_PLASMA, DOMAIN_SOL, DOMAIN_OUTER_SOL,
!! DOMAIN_UPPER_PRIVATE, DOMAIN_LOWER_PRIVATE
function particles_in_regions(node_list, element_list, particle_list)
  use data_structure
  use phys_module, only: DOMAIN_PLASMA, DOMAIN_SOL, DOMAIN_OUTER_SOL, DOMAIN_UPPER_PRIVATE,        &
      DOMAIN_LOWER_PRIVATE, xpoint, xcase
  use mod_particles
  use domains
  use mpi_mod
  implicit none

  type(type_node_list), intent(in)     :: node_list
  type(type_element_list), intent(in)  :: element_list
  type(type_particle_list), intent(in) :: particle_list

  integer, dimension(DOMAIN_PLASMA:DOMAIN_LOWER_PRIVATE) :: particles_in_regions
  integer :: i, ifail
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


  particles_in_regions = 0
  !$omp parallel do default(none) &
  !$omp shared(node_list, element_list, particle_list, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
  !$omp     R_axis, Z_axis, psi_axis) &
  !$omp private(domain, psi, psi_s, psi_t, psi_st, psi_ss, psi_tt, p) &
  !$omp reduction(+:particles_in_regions)
  do i=1,particle_list%n_particles
    p = particle_list%particle(i)
    call interp(node_list, element_list, p%i_elm, 1, 1, & ! force i_harm to 1
        p%st(1), p%st(2), psi, psi_s, psi_t, psi_st, psi_ss, psi_tt)

    domain = which_domain(node_list, element_list, &
        p%x(1), p%x(2), &
        psi, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
        R_axis, Z_axis, psi_axis)

    particles_in_regions(domain) = particles_in_regions(domain) + 1
  end do
  !$omp end parallel do


  !! Mpi communication to get the total answer to node 0
  call MPI_Gather(particles_in_regions, 5, MPI_INTEGER, particles_in_regions, 5, MPI_INTEGER, 0, MPI_COMM_WORLD, ifail)
end function particles_in_regions
end module mod_particle_diagnostics
