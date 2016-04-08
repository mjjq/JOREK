module mod_project_particles
contains
!> Project particles by weight onto the elements
!! Saves output in node_list%values(1)
!! The projection is done by integrating (2d equivalent of)
!! $$P_{i,j} = (\int_V \rho(x) B_i^3(s)B_j^3(t) dV)/(\int_V B_i^3(s)B_j^3(t) dV)$$
!! Where the dV transforms into $r dr d\phi$ in cylindrical coordinates.
!! The top integration is done by summing over all the particle positions, which
!! represent delta functions. (i.e. \rho(x) = \sum \delta(x_i - x) where x_i are
!! the particle positions)
!! The normalization is done by integration on the gaussian points.
subroutine project_particles(node_list, element_list, particle_list)
use phys_module
use data_structure
use basis_at_gaussian ! for HZ (initialise_basis must be called before use)
use mod_particles
use constants
implicit none

!> Input parameters
type(type_node_list), intent(inout)   :: node_list
type(type_element_list), intent(in)   :: element_list
type (type_particle_list), intent(in) :: particle_list

!> Preset parameters
integer, parameter :: into = 1 !< Project into this value index

!> Local variables
real*8  :: H2(4,4), H2_s(4,4), H2_t(4,4), ss, s, t
integer :: kv, iv, kf, i_harm, i_tor, i, j
integer :: i_elm

real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)
real*8  :: wst, volume, total_volume, xjac
integer :: ms, mt

type(type_node) :: node
type(type_element) :: element

real*8, dimension(1) :: P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t

call initialise_basis ! define the basis functions H at the Gaussian points

! Reset the density to zero
do i=1,node_list%n_nodes
  node_list%node(i)%values(:,:,into) = 0.d0
enddo

! Loop over particles to project each (contributions are normalized later)
do i=1,particle_list%n_particles
  if (particle_list%particle(i)%lost) cycle ! skip this particle
  i_elm  = particle_list%particle(i)%i_elm
  if (i_elm .lt. 1) cycle

  s = particle_list%particle(i)%st(1)
  t = particle_list%particle(i)%st(2)

  call basisfunctions3(s,t,H2,H2_s,H2_t)

  ! Get the local particle r coordinate for weighting
  do kv = 1,n_vertex_max  ! 4 vertices
    iv = element_list%element(i_elm)%vertex(kv)  ! the node number

    do kf = 1, n_order+1       ! 4 basis functions
      ss  = element_list%element(i_elm)%size(kv,kf)

      ! The contribution to the integral of this basis function (1,kf,into) in the n=0 mode:
      ! B3_i,j(s,t) * weight * R
      node_list%node(iv)%values(1,kf,into) = &
        node_list%node(iv)%values(1,kf,into) + ss * H2(kv,kf) * particle_list%particle(i)%weight * particle_list%particle(i)%x(1)

      ! The contribution of the toroidal fourier components
      do i_tor = 1, (n_tor-1)/2
        i_harm = 2*i_tor
        node_list%node(iv)%values(i_harm,kf,into) = &
        node_list%node(iv)%values(i_harm,kf,into) + &
          ss * H2(kv,kf) * particle_list%particle(i)%weight * &
          particle_list%particle(i)%x(1) * cos(mode(i_harm)*particle_list%particle(i)%x(3))
        node_list%node(iv)%values(i_harm,kf,into) = &
        node_list%node(iv)%values(i_harm,kf,into) + &
          ss * H2(kv,kf) * particle_list%particle(i)%weight * &
          particle_list%particle(i)%x(1) * sin(mode(i_harm+1)*particle_list%particle(i)%x(3))
      enddo
    enddo
  enddo

  ! Now print the local value for debugging
  !call interp_PRZ(node_list, element_list, i_elm, (/1/), 1, s, t, particle_list%particle(i)%x(3), &
  !P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
  !write(*,*) P
enddo


total_volume = 0.d0
! Normalize values by integral of basis function * R
do i_elm=1,element_list%n_elements
  ! Calculate volume integral of basis function * R to normalize
  x_g = 0.d0; x_s = 0.d0; x_t = 0.d0; y_g = 0.d0; y_s = 0.d0; y_t = 0.d0

  ! For each element
  element = element_list%element(i_elm)
  do i=1,n_vertex_max ! 4 vertices
    node = node_list%node(element%vertex(i))
    do j=1,n_order+1 ! 4 basis functions
      ! For each basis function B_i,j calculate the value at the gaussian points
      ! (without TWOPI factor as we did not include it above either!)
      volume = 0.d0
      do ms=1, n_gauss
        do mt=1, n_gauss
          wst = wgauss(ms)*wgauss(mt)
          x_g(ms,mt) = node%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
          y_g(ms,mt) = node%x(j,2) * element%size(i,j) * H(i,j,ms,mt)

          x_s(ms,mt) = node%x(j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = node%x(j,1) * element%size(i,j) * H_t(i,j,ms,mt)
          y_s(ms,mt) = node%x(j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = node%x(j,2) * element%size(i,j) * H_t(i,j,ms,mt)

          xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
          volume = volume + x_g(ms,mt) * wst ! x_g is R at the gaussian point
        enddo
      enddo

      ! Now that we have the volume, normalize the values corresponding to this
      ! basis function
      !node_list%node(element%vertex(i))%values(:,j,into) = &
      !node_list%node(element%vertex(i))%values(:,j,into) / volume
      !write(*,*) volume
      !write(*,*) node_list%node(element%vertex(i))%values(:,j,into)

      total_volume = total_volume + volume
    enddo
  enddo
enddo  ! n_elements

write(*,*) "Projection done"
write(*,*) "Volume: ", total_volume * TWOPI

end subroutine project_particles
end module mod_project_particles
