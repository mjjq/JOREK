module mod_project_particles
contains
!> Project particles by weight onto the elements
!! Saves output in node_list%values(1)
!! The projection is done by solving
!! $$\int \sum \delta(x_i-x)\delta(y_i-y) w_i u(x,y) dxdy = 
!! \int [p(x,y) u(x,y) + \lambda p'(x,y) u'(x,y)] dxdy$$
!! where p(x,y) is in the bernstein representation, $u(x,y)$ are the test functions
!! composed of two basis functions and $w_i$ is the particle weight.
!! A smoothing factor lambda is included
!! (equation is vaguely written!)
subroutine project_particles(node_list, element_list, particle_list)
use phys_module
use data_structure
use basis_at_gaussian
use mod_particles
use mpi_mod
implicit none

type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)
type (type_particle)     :: particle
type (type_particle_list):: particle_list
type (type_node_list), intent(inout) :: node_list !< A copy of the node list which will be used to save variables
type (type_element_list), intent(in) :: element_list

include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure

type (DMUMPS_STRUC) :: projection_matrix

real*8     :: ELM(n_vertex_max*(n_order+1),n_vertex_max*(n_order+1)), RHS(n_vertex_max*(n_order+1))
real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)
real*8     :: R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, x(3), HH(4,4), HH_s(4,4), HH_t(4,4)
real*8     :: v, v_x, v_y, psi, psi_x, psi_y, wst, area, volume, total_area, total_volume
integer    :: i, j, k, l, m, i_harm, ilarge, index_large_i, index_large_k, inode, knode
integer    :: nz_AA, n_AA, n_border, i_elm, index, ivar_out, index_ij, index_kl
integer    :: ms, mt, n_p, total_particles

nz_AA = element_list%n_elements * (n_vertex_max * (n_order+1))**2

n_border = 0
!do i=1,node_list%n_nodes
!  if (node_list%node(i)%boundary .eq. 1) n_border = n_border+2  ! INCLUDE OTHER BOUNDARY OPTIONS!
!  if (node_list%node(i)%boundary .eq. 2) n_border = n_border+2
!  if (node_list%node(i)%boundary .eq. 3) n_border = n_border+3
!enddo
!nz_AA = nz_AA + n_border

n_AA = 0
do inode = 1, node_list%n_nodes
  n_AA = max(n_AA,node_list%node(inode)%index(4))
enddo

write(*,*) ' number of unknowns      : ',n_AA, node_list%n_nodes * (n_order+1)
write(*,*) ' number of boundary nodes: ',n_border
write(*,*) ' nz_AA                   : ',nz_AA

projection_matrix%COMM = MPI_COMM_WORLD

projection_matrix%JOB      = -1
projection_matrix%SYM      = 0
projection_matrix%PAR      = 1

call DMUMPS(projection_matrix)

allocate(projection_matrix%A(nz_AA),projection_matrix%irn(nz_AA),projection_matrix%jcn(nz_AA))
allocate(projection_matrix%rhs(n_AA))

! Only n=0 for now



projection_matrix%irn = 0
projection_matrix%jcn = 0
projection_matrix%A   = 0.d0
projection_matrix%RHS = 0.d0

ilarge = 0
total_area   = 0.d0
total_volume = 0.d0
total_particles = 0

write(*,*) ' constructing projection matrix'

do i_elm  = 1, element_list%n_elements

  ELM = 0.d0
  RHS = 0.d0

  element = element_list%element(i_elm)

  do m=1,n_vertex_max
    nodes(m) = node_list%node(element%vertex(m))
  enddo

  n_p = 0

!--------------------------------------------------- sum over particle positions
  do m=1, particle_list%n_particles

    if (particle_list%particle(m)%i_elm .ne. i_elm) cycle
    if (particle_list%particle(m)%lost)             cycle

    n_p = n_p + 1

    particle = particle_list%particle(m)

    x(1:2) = particle%st
    x(3)   = particle%x(3)

    call interp3_RZ(node_list,element_list,i_elm,x(1),x(2),R_g,R_s,R_t,Z_g,Z_s,Z_t)

    xjac =  R_s*Z_t - R_t*Z_s

    call basisfunctions3(x(1), x(2), HH, HH_s, HH_t)

    do i=1,n_vertex_max

      do j=1,n_order+1

        index_ij = (i-1)*(n_order+1) + j

        v   = HH(i,j)  * element%size(i,j)
!      v_R = (  Z_t * H_s(i,j) - Z_s * H_t(i,j) ) * element%size(i,j) / xjac
!      v_Z = (- R_t * H_s(i,j) + R_s * H_t(i,j) ) * element%size(i,j) / xjac

        RHS(index_ij) = RHS(index_ij) + v * R_g  !* xjac

      enddo
    enddo

  enddo

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

  area = 0.
  volume = 0.

  do ms=1, n_gauss

    do mt=1, n_gauss

      wst = wgauss(ms)*wgauss(mt)

      xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

      area   = area   + xjac * wst
      volume = volume + x_g(ms,mt) * xjac * wst

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

              ELM(index_ij,index_kl) = ELM(index_ij,index_kl) + psi * v  * xjac * x_g(ms,mt) * wst! &

                                     !+ 0.001 * (psi_x * v_x + psi_y * v_y) * xjac * x_g(ms,mt) * wst

            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  total_volume    = total_volume    + volume
  total_area      = total_area      + area
  total_particles = total_particles + n_p

  write(111,'(i6,2f12.6,i8,2e14.6)') i_elm,sum(x_g)/16.,sum(y_g)/16.,n_p,volume,float(n_p)/volume


  do i=1,n_vertex_max

    inode = element%vertex(i)

    do j=1,n_order+1

      index_ij = (i-1)*(n_order+1) + j
      index_large_i = node_list%node(inode)%index(j)  ! base index in the main matrix

      projection_matrix%rhs(index_large_i) = projection_matrix%rhs(index_large_i) + RHS(index_ij)

      do k=1,n_vertex_max

        knode = element%vertex(k)

        do l=1,n_order+1

          index_kl = (k-1)*(n_order+1) + l

          index_large_k = node_list%node(knode)%index(l)  ! base index in the main matrix

          ilarge = ilarge + 1

          projection_matrix%irn(ilarge) = index_large_i
          projection_matrix%jcn(ilarge) = index_large_k
          projection_matrix%A(ilarge)   = ELM(index_ij,index_kl)

        enddo
      enddo
    enddo
  enddo

enddo

write(*,'(A,e14.6)') ' Area        : ',total_area
write(*,'(A,e14.6)') ' Volume      : ',total_volume
write(*,'(A,i14)')   ' Particles   : ',total_particles
write(*,'(A,e14.6)') ' Avg density : ',float(total_particles)/total_volume

write(*,*) ' constructed projection matrix ',ilarge

projection_matrix%n         = n_AA
projection_matrix%nz        = ilarge
projection_matrix%JOB       = 6
projection_matrix%icntl(5)  = 0
projection_matrix%icntl(18) = 0
projection_matrix%icntl(7)  = 4
projection_matrix%icntl(8)  = 7
projection_matrix%icntl(14) = 80

call DMUMPS(projection_matrix)

write(*,*) 'MUMPS finished'

ivar_out = 1
i_harm   = 1

write(*,*) maxval(projection_matrix%RHS),minval(projection_matrix%RHS)

do i=1,node_list%n_nodes

    do k=1,n_order+1

      index = node_list%node(i)%index(k)

      node_list%node(i)%values(i_harm,k,ivar_out) = projection_matrix%RHS(index)

    enddo    ! order

enddo        ! nodes

projection_matrix%JOB      = -2
call DMUMPS(projection_matrix)

return
end subroutine project_particles
end module mod_project_particles
