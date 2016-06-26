module mod_project_particles
contains
!> Project particles by weight onto the elements
!! Saves output in node_list%values(1)
!! The projection is done by solving
!! $$\int [p(x) u(x) + \lambda p'(x) u'(x)] dV =
!! \int \sum \delta(x_i-x) w_i u(x) dV $$
!! where p(x) is in the bernstein representation, $u(x)$ are the test functions
!! composed of two basis functions and $w_i$ is the particle weight.
!! A smoothing factor lambda is included
!! x is a vector (R,Z,phi) and dV is r dr dphi
!! divide by 1 or 2pi on both sides (LHS gets 2pi for n=0 mode, 1pi for other modes)
subroutine project_particles(node_list, element_list, particle_list)
use phys_module
use data_structure
use basis_at_gaussian
use mod_particles
use mpi_mod
implicit none

type (type_particle_list), intent(in):: particle_list
type (type_node_list), intent(inout) :: node_list !< A copy of the node list which will be used to save variables
type (type_element_list), intent(in) :: element_list

type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)
type (type_particle)     :: particle

integer, parameter :: ivar_out = 1

include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure

type (DMUMPS_STRUC) :: projection_matrix

real*8     :: ELM(n_vertex_max*(n_order+1),n_vertex_max*(n_order+1)), RHS(n_vertex_max*(n_order+1),element_list%n_elements)
real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)
real*8     :: R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, x(3), HH(4,4), HH_s(4,4), HH_t(4,4)
real*8     :: v, v_x, v_y, psi, psi_x, psi_y, wst, area, volume
integer    :: i, j, k, l, m, i_tor, ilarge, index_large_i, index_large_k, inode, knode
integer    :: nz_AA, n_AA, n_border, i_elm, index, index_ij, index_kl
integer    :: ms, mt, n_p

real*8, parameter :: smoothing = 1d-3

nz_AA = element_list%n_elements * (n_vertex_max * (n_order+1))**2

n_AA = 0
do inode = 1, node_list%n_nodes
  n_AA = max(n_AA,node_list%node(inode)%index(4))
enddo
write(*,*) ' number of unknowns      : ',n_AA, node_list%n_nodes * (n_order+1)
write(*,*) ' nz_AA                   : ',nz_AA

! Initialise MUMPS
projection_matrix%COMM = MPI_COMM_WORLD
projection_matrix%JOB      = -1
projection_matrix%SYM      = 0
projection_matrix%PAR      = 1
call DMUMPS(projection_matrix)

! Allocate space for elements
allocate(projection_matrix%A(nz_AA),projection_matrix%irn(nz_AA),projection_matrix%jcn(nz_AA))
allocate(projection_matrix%rhs(n_AA))
projection_matrix%irn = 0
projection_matrix%jcn = 0
projection_matrix%A   = 0.d0
projection_matrix%RHS = 0.d0

area = 0.
volume = 0.

write(*,*) 'constructing projection matrix'
!$omp parallel do default(none) &
!$omp shared(element_list, node_list, H, H_s, H_t, projection_matrix) &
!$omp private(ELM, i_elm, element, nodes, i, j, ms, mt, &
!$omp         x_g, y_g, x_s, x_t, y_s, y_t, wst, xjac, &
!$omp         index_ij, index_kl, v, v_x, v_y, psi, psi_x, psi_y, ilarge, &
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
      wst = wgauss(ms)*wgauss(mt)
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

          projection_matrix%irn(ilarge) = index_large_i
          projection_matrix%jcn(ilarge) = index_large_k
          projection_matrix%A(ilarge)   = ELM(index_ij,index_kl)
        enddo
      enddo
    enddo
  enddo
enddo
!$omp end parallel do

write(*,'(A,e14.6)') ' Area        : ',area
write(*,'(A,e14.6)') ' Volume      : ',volume

!
! Perform the analysis and factorisation
!
projection_matrix%JOB       = 4
projection_matrix%n         = n_AA
projection_matrix%nz        = nz_AA
projection_matrix%icntl(5)  = 0
projection_matrix%icntl(18) = 0
projection_matrix%icntl(7)  = 4
projection_matrix%icntl(8)  = 7
projection_matrix%icntl(14) = 80
call DMUMPS(projection_matrix)

!
! Create the RHS and calculate the projection
!
do i_tor=1, n_tor
  n_p = 0
  RHS = 0.d0
  projection_matrix%rhs = 0.d0
  !$omp parallel do default(none) &
  !$omp shared(particle_list, element_list, node_list, mode, i_tor) &
  !$omp private(particle, x, xjac, HH, HH_s, HH_t, R_g, R_s, R_t, Z_g, Z_s, Z_t, &
  !$omp         i, j, index_ij, v, m, i_elm, index_large_i, inode) &
  !$omp reduction(+:n_p,RHS)
  do m=1, particle_list%n_particles
    if (particle_list%particle(m)%lost) cycle
    n_p = n_p + 1

    particle = particle_list%particle(m)

    x(1:2) = particle%st
    x(3)   = particle%x(3)

    call interp3_RZ(node_list,element_list,particle%i_elm,x(1),x(2),R_g,R_s,R_t,Z_g,Z_s,Z_t)
    xjac =  R_s*Z_t - R_t*Z_s
    call basisfunctions3(x(1), x(2), HH, HH_s, HH_t)
    do i=1,n_vertex_max
      do j=1,n_order+1
        index_ij = (i-1)*(n_order+1) + j

        v   = HH(i,j)  * element_list%element(particle%i_elm)%size(i,j)
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

        RHS(index_ij,particle%i_elm) = RHS(index_ij,particle%i_elm) + v * particle%weight
      enddo
    enddo
  enddo
  !$omp end parallel do

  ! Fill RHS of Projection matrix
  do i_elm=1,element_list%n_elements
    do i=1,n_vertex_max
      inode = element_list%element(i_elm)%vertex(i)
      do j=1,n_order+1
        index_ij = (i-1)*(n_order+1) + j
        index_large_i = node_list%node(inode)%index(j)  ! base index in the main matrix

        projection_matrix%rhs(index_large_i) = projection_matrix%rhs(index_large_i) + RHS(index_ij, i_elm)
      enddo
    enddo
  enddo

  if (i_tor .eq. 1) then
    write(*,'(A,i14)')   ' Particles   : ',n_p
    write(*,'(A,e14.6)') ' Avg density : ',float(n_p)/volume
  endif

  ! Compute the solution of Ax=b (b = RHS)
  projection_matrix%JOB = 3
  call DMUMPS(projection_matrix)

  write(*,*) 'Projection ', i_tor, ' finished'
  write(*,*) minval(projection_matrix%RHS), maxval(projection_matrix%RHS)

  do i=1,node_list%n_nodes
    do k=1,n_order+1
      index = node_list%node(i)%index(k)
      node_list%node(i)%values(i_tor,k,ivar_out) = projection_matrix%RHS(index)
    enddo    ! order
  enddo      ! nodes
enddo ! i_tor

! Clean up MUMPS
projection_matrix%JOB      = -2
call DMUMPS(projection_matrix)
end subroutine project_particles
end module mod_project_particles
