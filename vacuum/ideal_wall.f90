subroutine ideal_wall(my_id,node_list,boundary_list)
!-------------------------------------------------------------------
! routine calculates the vacuum reponse matrix of a vacuum region
! up to an ideally conducting wall
!-------------------------------------------------------------------
use data_structure
use vacuum_response_module
implicit none

type(type_node_list)     :: node_list
type(type_boundary_list) :: boundary_list
type(type_node_list)     :: vacuum_node_list
type(type_element_list)  :: vacuum_element_list

real*8  :: r_wall
integer :: my_id, imode, i
logical :: xpoint

if (my_id .eq. 0) then

  write(*,*) '******************************'
  write(*,*) '* ideal wall vacuum response *'
  write(*,*) '******************************'

  r_wall = 2.d0

  call vacuum_grid(node_list,boundary_list,r_wall,vacuum_node_list,vacuum_element_list)

  call plot_grid(vacuum_node_list,vacuum_element_list,boundary_list,.true.,.true.)    ! plot the grid

endif


imode = 2
xpoint = .false.

call vacuum_Poisson(my_id,vacuum_node_list,vacuum_element_list,boundary_list,imode)

call plot_solution(vacuum_node_list,vacuum_element_list,1,-1,1,'vacuum')

return
end

subroutine vacuum_poisson(my_id,node_list,element_list,boundary_list,imode)
!---------------------------------------------------------------
!---------------------------------------------------------------
use data_structure
use mumps_module
use basis_at_gaussian
use gauss
use phys_module
use vacuum_response_module

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_boundary_list):: boundary_list
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8   :: angle, pi, xx, xxs, yy, yys, pp, ps, ss, HH(2,2), HH_s(2,2), HH_ss(2,2)

real*8   :: ELM(n_vertex_max*(n_order+1),n_vertex_max*(n_order+1))
real*8   :: zbig, dl, ws, psi, psi_s, v, tht, xs, ys
real*8   :: x_g(n_gauss), x_s(n_gauss), y_g(n_gauss), y_s(n_gauss), eq_g(n_gauss), eq_s(n_gauss)
real*8   :: am0, am1, am2, am3, am4, adrive
integer  :: ierr, my_id, ife, iv, imode, itor, ibasis, ms
integer  :: n_border, ilarge, n_AA, nz_AA, nz_AA_old, i,j, k,l, ibnd, jbnd, inode1, inode2, jdir
integer  :: n_elements, index_large_i, inode, knode, index_large_k, index_ij, index_kl, index, index_i
integer  :: vertex(2), dir(2), iv2, sms, index_basis, index_basis2, index_basis_bnd, index_basis2_bnd, kbnd, lbnd

pi = 2.*asin(1.d0)

if (my_id .eq. 0) then

  write(*,*) '**************************************'
  write(*,*) '*           Vacuum Poisson           *'
  write(*,*) '**************************************'
  write(*,*) ' n_elements : ',element_list%n_elements
  write(*,*) ' n_nodes    : ',node_list%n_nodes
  write(*,*) ' n_bnd      : ',boundary_list%n_boundary

  nz_AA = element_list%n_elements * (n_vertex_max * (n_order+1))**2

  n_border = 0
  do i=1,node_list%n_nodes
   if (node_list%node(i)%boundary .eq. 1) n_border = n_border+2
   if (node_list%node(i)%boundary .eq. 2) n_border = n_border+2
   if (node_list%node(i)%boundary .eq. 3) n_border = n_border+3
   if (node_list%node(i)%boundary .eq. 4) n_border = n_border+2
  enddo

  nz_AA = nz_AA + n_border

  n_AA  = node_list%n_nodes * (n_order+1)

  n_AA = 0
  do inode = 1, node_list%n_nodes
    n_AA = max(n_AA,node_list%node(inode)%index(4))
  enddo

  write(*,*) ' number of unknowns      : ',n_AA, node_list%n_nodes * (n_order+1)

  write(*,*) ' RHS  : ',n_AA*n_dof_bnd

  if (.not. associated(mumps_par%A))     allocate(mumps_par%A(nz_AA))
  if (.not. associated(mumps_par%rhs))   allocate(mumps_par%rhs(n_AA*n_dof_bnd))
  if (.not. associated(mumps_par%irn))   allocate(mumps_par%irn(nz_AA))
  if (.not. associated(mumps_par%jcn))   allocate(mumps_par%jcn(nz_AA))

  mumps_par%irn = 0
  mumps_par%jcn = 0
  mumps_par%A   = 0.d0
  mumps_par%RHS = 0.d0

  mumps_par%lrhs = n_AA
  mumps_par%nrhs = n_dof_bnd

  n_elements = element_list%n_elements

  ilarge=0

  do ife =1, n_elements

    element = element_list%element(ife)

    do iv = 1, n_vertex_max

      inode     = element%vertex(iv)
      nodes(iv) = node_list%node(inode)

    enddo

    itor = mode(imode)

    call element_matrix_vacuum(element,nodes,itor,ELM)

    do i=1,n_vertex_max

      inode         = element%vertex(i)

      do j=1,n_order+1

        index_ij = (i-1)*(n_order+1) + j     ! index in the ELM matrix

        index_large_i = node_list%node(inode)%index(j)  ! base index in the main matrix

        do k=1,n_vertex_max

          knode         = element%vertex(k)

          do l=1,n_order+1

            index_kl = (k-1)*(n_order+1) + l

            index_large_k = node_list%node(knode)%index(l)  ! base index in the main matrix

            ilarge = ilarge +1

            mumps_par%irn(ilarge) = index_large_i
            mumps_par%jcn(ilarge) = index_large_k
            mumps_par%A(ilarge)   = ELM(index_ij,index_kl)

          enddo
        enddo
      enddo
    enddo

  enddo

  nz_AA_old = nz_AA
  nz_AA = ilarge

  zbig = 1.d10

  mumps_par%nrhs = n_dof_bnd

  do ibnd = 1, boundary_list%n_boundary
  do jbnd = 1, 2

    index_basis = (ibnd-1) * 2 + jbnd - 1

    inode  = ibnd
    ibasis = 1 + 2*(jbnd-1)      ! only ibasis = 1 or 3 for now

    do i=1, node_list%n_nodes
      node_list%node(i)%values = 0.d0
    enddo

    node_list%node(inode)%values(1,ibasis,1) = 1.d0
!----------------------------------------------------- test unit perturbation Fourier harmonic
!    sms= 3.d0
!    do inode=1, boundary_list%n_boundary
!      tht = atan2( node_list%node(inode)%x(1,2), node_list%node(inode)%x(1,1)-10.d0)
!      node_list%node(inode)%values(1,1,1) = cos(sms*tht)
!      node_list%node(inode)%values(1,3,1) = -1.d0 / 3.d0 * sms*sin(sms*tht) * 6.283185307 / float( boundary_list%n_boundary)
!    enddo


!integrate over the boundary elements


    do ife = 1, element_list%n_elements                                                                   ! boundary integrals

      do iv = 1, n_vertex_max                                                                     ! boundary integrals

        iv2  = mod(iv, n_vertex_max) + 1

        inode1 = element_list%element(ife)%vertex(iv)
        inode2 = element_list%element(ife)%vertex(iv2)

        if (     ((node_list%node(inode1)%boundary .eq. 4) .or.(node_list%node(inode1)%boundary .eq. 3)) &
           .and. ((node_list%node(inode2)%boundary .eq. 4) .or.(node_list%node(inode2)%boundary .eq. 3)) ) then

          nodes(1) = node_list%node(inode1)
	  nodes(2) = node_list%node(inode2)
	  vertex   = (/ iv, iv2 /)

          dir      = (/ 1, 3 /) ! not correct, depends on node not on element side  -->  take it from boundary_list

          x_g = 0.d0; x_s = 0.d0; y_g = 0.d0; y_s = 0.d0; eq_g = 0.d0; eq_s = 0.d0

          do i=1,2

            do j=1,2

              do ms=1, n_gauss

                x_g(ms)  = x_g(ms)  + nodes(i)%x(dir(j),1) * element_list%element(ife)%size(vertex(i),dir(j)) * H1(i,j,ms)
                x_s(ms)  = x_s(ms)  + nodes(i)%x(dir(j),1) * element_list%element(ife)%size(vertex(i),dir(j)) * H1_s(i,j,ms)

                y_g(ms)  = y_g(ms)  + nodes(i)%x(dir(j),2) * element_list%element(ife)%size(vertex(i),dir(j)) * H1(i,j,ms)
                y_s(ms)  = y_s(ms)  + nodes(i)%x(dir(j),2) * element_list%element(ife)%size(vertex(i),dir(j)) * H1_s(i,j,ms)

                eq_g(ms)  = eq_g(ms) + nodes(i)%values(1,dir(j),1) * element_list%element(ife)%size(vertex(i),dir(j)) * H1(i,j,ms)
                eq_s(ms)  = eq_s(ms) + nodes(i)%values(1,dir(j),1) * element_list%element(ife)%size(vertex(i),dir(j)) * H1_s(i,j,ms)

              enddo
            enddo
          enddo

          do ms=1, n_gauss

            ws = wgauss(ms)

            psi   = eq_g(ms)
            psi_s = eq_s(ms)

	    dl = sqrt(x_s(ms)**2 + y_s(ms)**2)

            do i=1,2                                                                       ! loop over nodes

              do j=1,2                                                                     ! loop over basis functions

                index_ij = nodes(i)%index(dir(j))                                          ! index in the ELM matrix

                v   =  H1(i,j,ms) * element_list%element(ife)%size(vertex(i),dir(j))       ! test function

                mumps_par%RHS(index_ij+n_AA*index_basis) = mumps_par%RHS(index_ij+n_AA*index_basis) + v * psi_s * ws              ! add to element RHS

              enddo
            enddo

          enddo

        endif

      enddo
    enddo

  enddo
  enddo
  
!----------------------- boundary conditions

  do i=1,node_list%n_nodes

    if ((node_list%node(i)%boundary .eq. 1) .or. (node_list%node(i)%boundary .eq. 3)) then

      index_i = node_list%node(i)%index(1)  ! base index in the main matrix

      mumps_par%irn(ilarge+1) = index_i
      mumps_par%jcn(ilarge+1) = index_i
      mumps_par%A(ilarge+1)   = zbig
      ilarge = ilarge + 1

      index_i = node_list%node(i)%index(2)  ! base index in the main matrix

      mumps_par%irn(ilarge+1) = index_i
      mumps_par%jcn(ilarge+1) = index_i
      mumps_par%A(ilarge+1)   = zbig
      ilarge = ilarge + 1
    endif

    if ((node_list%node(i)%boundary .eq. 2) .or. (node_list%node(i)%boundary .eq. 3)) then

      index_i = node_list%node(i)%index(1)  ! base index in the main matrix

      mumps_par%irn(ilarge+1) = index_i
      mumps_par%jcn(ilarge+1) = index_i
      mumps_par%A(ilarge+1)   = zbig
      ilarge = ilarge + 1

      index_i = node_list%node(i)%index(3)  ! base index in the main matrix

      mumps_par%irn(ilarge+1) = index_i
      mumps_par%jcn(ilarge+1) = index_i
      mumps_par%A(ilarge+1)   = zbig
      ilarge = ilarge + 1
    endif

  enddo

  nz_AA_old = nz_AA
  nz_AA     = ilarge

  mumps_par%n  = n_AA
  mumps_par%nz = nz_AA

endif

mumps_par%JOB = 6
mumps_par%SYM = 0
mumps_par%icntl(7) = 4

call DMUMPS(mumps_par)

vacuum_response = 0.d0

if (my_id .eq. 0) then

  do ibnd = 1, boundary_list%n_boundary

    do jbnd = 1, 2

      index_basis     = (n_order+1)*(ibnd-1) + 2*(jbnd-1) + 1          ! select index 1 and 3 from the 4 dof at each node
      index_basis_bnd = 2*(ibnd-1) +   (jbnd-1) + 1                    ! the index in the vacuum_response matrix

      do kbnd=1,boundary_list%n_boundary

        do lbnd = 1, 2

          index_basis2     = (n_order+1)*(kbnd-1) + 2*(lbnd-1) + 1
          index_basis2_bnd = 2*(kbnd-1) +   (lbnd-1) + 1               ! the index in the vacuum_response matrix

          vacuum_response(index_basis_bnd,index_basis2_bnd,2:n_tor) = mumps_par%rhs(index_basis2 + n_AA * (index_basis_bnd-1))

        enddo

      enddo

      write(*,'(16e11.3)') vacuum_response(index_basis_bnd,:,2)

    enddo
  enddo

  do i=1,node_list%n_nodes

    inode  = i

    do j=1, n_order+1

      index  = node_list%node(inode)%index(j)

      node_list%node(inode)%values(1,j,1) = mumps_par%RHS(index+(boundary_list%n_boundary+0)*n_AA)

    enddo
  enddo

  deallocate(mumps_par%irn,mumps_par%jcn,mumps_par%A,mumps_par%rhs)

!### BEGIN DEBUNG OUTPUT ###

    do ife = 1, element_list%n_elements                                                                   ! boundary integrals

      do iv = 1, n_vertex_max                                                                     ! boundary integrals

        iv2  = mod(iv, n_vertex_max) + 1

        inode1 = element_list%element(ife)%vertex(iv)
        inode2 = element_list%element(ife)%vertex(iv2)

        if (     ((node_list%node(inode1)%boundary .eq. 4) .or.(node_list%node(inode1)%boundary .eq. 3)) &
           .and. ((node_list%node(inode2)%boundary .eq. 4) .or.(node_list%node(inode2)%boundary .eq. 3)) ) then

          nodes(1) = node_list%node(inode1)
    nodes(2) = node_list%node(inode2)
    vertex   = (/ iv, iv2 /)

          dir      = (/ 1, 3 /) ! not correct, depends on node not on element side  -->  take it from boundary_list

          do ms=1, 10
          
            xx=0.; xxs=0.; yy=0.; yys=0.; pp=0.; ps=0.

            ss = (ms-1.)/9.
            
            call basisfunctions1(ss,HH,HH_s,HH_ss)

            do i=1,2 ! vertex number

              do j=1,2 ! basis function

                xx = xx + nodes(i)%x(dir(j),1) * element_list%element(ife)%size(vertex(i),dir(j)) * HH(i,j)
                xxs = xxs + nodes(i)%x(dir(j),1) * element_list%element(ife)%size(vertex(i),dir(j)) * HH_s(i,j)

                yy = yy + nodes(i)%x(dir(j),2) * element_list%element(ife)%size(vertex(i),dir(j)) * HH(i,j)
                yys = yys + nodes(i)%x(dir(j),2) * element_list%element(ife)%size(vertex(i),dir(j)) * HH_s(i,j)

                pp = pp + nodes(i)%values(1,dir(j),1) * element_list%element(ife)%size(vertex(i),dir(j)) * HH(i,j)
                ps = ps + nodes(i)%values(1,dir(j),1) * element_list%element(ife)%size(vertex(i),dir(j)) * HH_s(i,j)
                
              enddo
            enddo
            
            888 format(ES12.5,1x,ES12.5)
            angle = MOD( atan2(yy, xx-R_geo) + 2*PI, 2*PI )
            write(67, 888) angle, pp
            write(68, 888) angle, ps
            write(69, 888) angle, ps / SQRT(xxs**2 + yys**2)
            write(70, 888) angle, SQRT(xxs**2 + yys**2)
            write(71, 888) SQRT(xxs**2 + yys**2)

          enddo
          
        endif

      enddo
    enddo
    

!### END DEBUG OUTPUT ###



endif


return
end

subroutine element_matrix_vacuum(element,nodes,itor,ELM)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!---------------------------------------------------------------
use parameters
use data_structure
use gauss
use basis_at_gaussian

implicit none

type (type_element) :: element
type (type_node)    :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)
real*8     :: ELM(n_vertex_max*(n_order+1),n_vertex_max*(n_order+1))

real*8     :: xjac, wst
real*8     :: v, v_x, v_y, psi, psi_x, psi_y, rhs_ij
integer    :: ms, mt, i, j, k, l, index_ij, index_kl, itor

ELM=0.d0

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g(:,:)   = 0.d0; x_s(:,:)   = 0.d0; x_t(:,:)   = 0.d0;
y_g(:,:)   = 0.d0; y_s(:,:)   = 0.d0; y_t(:,:)   = 0.d0;

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


!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

 do mt=1, n_gauss

   wst = wgauss(ms)*wgauss(mt)

   xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

   do i=1,n_vertex_max

     do j=1,n_order+1

       index_ij = (i-1)*(n_order+1) + j

       v   =  H(i,j,ms,mt) * element%size(i,j)
       v_x = (  y_t(ms,mt) * h_s(i,j,ms,mt) - y_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac
       v_y = (- x_t(ms,mt) * h_s(i,j,ms,mt) + x_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac

       do k=1,n_vertex_max

         do l=1,n_order+1

           psi   = H(k,l,ms,mt) * element%size(k,l)
           psi_x = (   y_t(ms,mt) * h_s(k,l,ms,mt) - y_s(ms,mt) * h_t(k,l,ms,mt) ) * element%size(k,l) / xjac
           psi_y = ( - x_t(ms,mt) * h_s(k,l,ms,mt) + x_s(ms,mt) * h_t(k,l,ms,mt) ) * element%size(k,l) / xjac

           index_kl = (k-1)*(n_order+1) + l

           ELM(index_ij,index_kl) =  ELM(index_ij,index_kl) &
                            + (psi_x * v_x + psi_y * v_y + 0.5d0*(float(itor)/x_g(ms,mt))**2 * v * psi) * x_g(ms,mt)*xjac*wst
                           !+ (psi_x * v_x + psi_y * v_y + 0.5d0*(float(itor)/x_g(ms,mt))**2 * v * psi) * xjac * wst

         enddo
       enddo

     enddo
   enddo

 enddo
enddo

return
end

subroutine vacuum_grid(node_list,boundary_list,r_wall,vacuum_node_list,vacuum_element_list)
!-------------------------------------------------------------------
! routine calculates the vacuum reponse matrix of a vacuum region
! up to an ideally conducting wall
!-------------------------------------------------------------------
use data_structure

implicit none

type(type_boundary_list) :: boundary_list
type(type_node_list)     :: node_list, vacuum_node_list
type(type_element_list)  :: vacuum_element_list

integer :: nr, np, i, j, k, jv, index, iv, ivp, ivm, node_iv, node_ivp, node_ivm
real*8  :: r_wall, radius, ds, delta_Rp, delta_Zp, delta_Rm, delta_Zm, dir_2, dir_3, RZcentre(2), xlength

write(*,*) '******************************'
write(*,*) '* vacuum grid                *'
write(*,*) '******************************'

np = boundary_list%n_boundary
nr = 21

write(*,*) ' nr, np : ',nr,np

do j=1,np
  jv = boundary_list%boundary(j)%vertex(1)
  RZcentre = RZcentre + node_list%node(jv)%x(1,:)
enddo
RZcentre = RZcentre / float(np)

do i=1,nr

  radius = 1.d0 + (r_wall - 1.d0) * float(i-1)/float(nr-1)

  do j=1,np

    index = (i-1)*np + j

    jv = boundary_list%boundary(j)%vertex(1)

    vacuum_node_list%node(index)%x(1,:) = radius * (node_list%node(jv)%x(1,:) - RZcentre) + RZcentre
    vacuum_node_list%node(index)%x(3,:) = radius *  node_list%node(jv)%x(3,:)

    xlength = sqrt((node_list%node(jv)%x(1,1) - RZcentre(1))**2 + ((node_list%node(jv)%x(1,2) - RZcentre(2))**2))

    vacuum_node_list%node(index)%x(2,:) = 1.d0/3.d0 * (r_wall - 1.d0) / float(nr-1) * (node_list%node(jv)%x(1,:) - RZcentre) / xlength

    vacuum_node_list%node(index)%x(4,:) = 1.d0/3.d0 * (r_wall - 1.d0) / float(nr-1) * node_list%node(jv)%x(3,:) / xlength

    vacuum_node_list%node(index)%boundary = 0

    if (i .eq. nr) vacuum_node_list%node(index)%boundary = 2
    if (i .eq. 1)  vacuum_node_list%node(index)%boundary = 4

    do k=1,n_order+1
      vacuum_node_list%node(index)%index(k) =  (n_order+1)*(index-1)+k
    enddo

  enddo

enddo

vacuum_node_list%n_nodes        = nr * np
vacuum_element_list%n_elements  = (nr-1)*np

do i=1,nr-1

  do j=1,np-1

    index = np*(i-1) + j

    vacuum_element_list%element(index)%vertex(1) = (i-1)*np + j
    vacuum_element_list%element(index)%vertex(4) = (i-1)*np + j + 1
    vacuum_element_list%element(index)%vertex(3) = (i  )*np + j + 1
    vacuum_element_list%element(index)%vertex(2) = (i  )*np + j

  enddo

  index =  np*(i-1) + np

  vacuum_element_list%element(index)%vertex(1)  = (i  )*np
  vacuum_element_list%element(index)%vertex(4)  = (i  )*np - np + 1
  vacuum_element_list%element(index)%vertex(3)  = (i  )*np + 1
  vacuum_element_list%element(index)%vertex(2)  = (i  )*np + np

enddo

do k=1 , vacuum_element_list%n_elements   ! fill in the size of the elements

 do iv = 1, 4                    ! over 4 corners of an element

   ivp = mod(iv,4)   + 1         ! vertex with index one higher
   ivm = mod(iv+2,4) + 1         ! vertex with index one below

   node_iv  = vacuum_element_list%element(k)%vertex(iv)
   node_ivp = vacuum_element_list%element(k)%vertex(ivp)
   node_ivm = vacuum_element_list%element(k)%vertex(ivm)

   if ((iv .eq. 1) .or. (iv .eq.3)) then

     delta_Rp = vacuum_node_list%node(node_ivp)%X(1,1) - vacuum_node_list%node(node_iv)%X(1,1)
     delta_Zp = vacuum_node_list%node(node_ivp)%X(1,2) - vacuum_node_list%node(node_iv)%X(1,2)
     dir_2    = delta_Rp * vacuum_node_list%node(node_iv)%X(2,1) + delta_Zp * vacuum_node_list%node(node_iv)%X(2,2)

     delta_Rm = vacuum_node_list%node(node_ivm)%X(1,1) - vacuum_node_list%node(node_iv)%X(1,1)
     delta_Zm = vacuum_node_list%node(node_ivm)%X(1,2) - vacuum_node_list%node(node_iv)%X(1,2)
     dir_3    = delta_Rm * vacuum_node_list%node(node_iv)%X(3,1) + delta_Zm * vacuum_node_list%node(node_iv)%X(3,2)

   else

     delta_Rp = vacuum_node_list%node(node_ivp)%X(1,1) - vacuum_node_list%node(node_iv)%X(1,1)
     delta_Zp = vacuum_node_list%node(node_ivp)%X(1,2) - vacuum_node_list%node(node_iv)%X(1,2)
     dir_3    = delta_Rp * vacuum_node_list%node(node_iv)%X(3,1) + delta_Zp * vacuum_node_list%node(node_iv)%X(3,2)

     delta_Rm = vacuum_node_list%node(node_ivm)%X(1,1) - vacuum_node_list%node(node_iv)%X(1,1)
     delta_Zm = vacuum_node_list%node(node_ivm)%X(1,2) - vacuum_node_list%node(node_iv)%X(1,2)
     dir_2    = delta_Rm * vacuum_node_list%node(node_iv)%X(2,1) + delta_Zm * vacuum_node_list%node(node_iv)%X(2,2)

   endif

   if (dir_2 .ne. 0.d0) then
     dir_2 = dir_2 / abs(dir_2)
   else
     dir_2 = 1.d0
   endif
   if (dir_3 .ne. 0.d0) then
     dir_3 = dir_3 / abs(dir_3)
   else
     dir_3 = -1.d0
     if (iv.eq.1) dir_3 = 1.d0              ! admittedly not very elegant
   endif

   vacuum_element_list%element(k)%size(iv,1) = 1.d0
   vacuum_element_list%element(k)%size(iv,2) = dir_2
   vacuum_element_list%element(k)%size(iv,3) = dir_3
   vacuum_element_list%element(k)%size(iv,4) = vacuum_element_list%element(k)%size(iv,2) * vacuum_element_list%element(k)%size(iv,3)

 enddo

enddo

return
end





subroutine ideal_wall_starwall(my_id,node_list,boundary_list)

  use data_structure
  use vacuum_response_module
  use phys_module
  implicit none

  integer,                  intent(in) :: my_id            ! MPI thread number of current thread
  type(type_node_list),     intent(in) :: node_list        ! List of boundary nodes
  type(type_boundary_list), intent(in) :: boundary_list    ! List of boundary elements

  
  character(len=128) :: file_response_starwall = 'vacuum_response_with_wall_R10_a1.0_w2a_npol16.starwall'
  
  integer :: dim
  integer :: i_n, j_n
  integer :: itor, jtor
  integer :: inode, jnode
  integer :: ibas, jbas
  integer :: icossin, jcossin ! 0: cos, 1: sin
  integer :: iindex, jindex
  integer :: rn_response, cn_response
  integer :: iindex2, jindex2 !###
  real*8  :: response
  
  OPEN(42, FILE=file_response_starwall)
  
  READ(42,*) dim
  WRITE(*,*) 'dim=',dim
  
  vacuum_response = 0.
  
  ! Outer loops (response index)
  do inode = 1, boundary_list%n_boundary ! loop over nodes
    do itor = 2,3 !###                     ! loop over toroidal modes
      do ibas = 1, 2                         ! loop over basis functions
        
        iindex = 2*(inode-1) +   (ibas-1) + 1  ! first index in response matrix
        
        jindex = 0
        
        ! Inner loops (perturbation index)
        do jnode = 1, boundary_list%n_boundary ! loop over nodes
          do jtor = 2, 3 !###                    ! loop over toroidal harmonics
            do jbas = 1, 2                         ! loop over basis functions
              
              jindex = 2*(jnode-1) +   (jbas-1) + 1  ! second index in response matrix
                
              READ(42,*) iindex2, jindex2, response
              
              if ( itor == jtor ) vacuum_response(iindex,jindex,itor) = response
                  
            end do
          end do
        end do
 
      end do
    end do
  end do

  CLOSE(42)

end subroutine ideal_wall_starwall
