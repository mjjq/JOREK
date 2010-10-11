subroutine ideal_wall(my_id, node_list, boundary_list, bnd_node_list)
!-------------------------------------------------------------------
! routine calculates the vacuum reponse matrix of a vacuum region
! up to an ideally conducting wall
!-------------------------------------------------------------------
  
  use data_structure
  use vacuum_response_module
  use mumps_module
  
  implicit none
  
  include 'mpif.h'
  
  integer,                     intent(in) :: my_id            ! MPI thread number of current thread
  type(type_node_list),        intent(in) :: node_list        ! List of boundary nodes
  type(type_bnd_element_list), intent(in) :: boundary_list    ! List of boundary elements
  type(type_bnd_node_list),    intent(in) :: bnd_node_list    ! List of boundary nodes
  
  
  type(type_node_list)        :: vacuum_node_list
  type(type_element_list)     :: vacuum_element_list
  real*8  :: r_wall
  integer :: imode, i
  logical :: xpoint
  

  call initialise_mumps(MPI_COMM_WORLD)                  ! start MUMPS sparse matrix solver
    
  if (my_id .eq. 0) then
  
    write(*,*) '******************************'
    write(*,*) '*       ideal_wall           *'
    write(*,*) '******************************'
  
    r_wall = 2.d0
  
    call vacuum_grid(node_list,boundary_list,r_wall,vacuum_node_list,vacuum_element_list)
  
    call plot_grid(vacuum_node_list,vacuum_element_list,boundary_list,bnd_node_list,.true.,.true.)    ! plot the grid
  
  endif
  
  
  xpoint = .false.
  
  do i=1,(n_tor+1)/2
  
    if (i .eq. 1) then
      imode = 1
    else
      imode = 2*(i-1)
    endif
  
    write(*,*) ' JOREK Vacuum : ',i,imode
  
    call vacuum_Poisson(my_id,vacuum_node_list,vacuum_element_list,boundary_list,node_list,imode)
  
    if (imode .gt. 1) vacuum_response(:,:,imode+1) = vacuum_response(:,:,imode)
  
  enddo
  
  call plot_solution(vacuum_node_list,vacuum_element_list,1,-1,1,'vacuum')
#ifdef USE_MUMPS
  mumps_par%JOB = -2                                     ! clean up this instance of mumps
  call DMUMPS(mumps_par)
#endif
  return
end subroutine ideal_wall



subroutine vacuum_poisson(my_id,node_list,element_list,boundary_list,temp_node_list,imode)
!---------------------------------------------------------------
!---------------------------------------------------------------
  use data_structure
  use mumps_module
  use basis_at_gaussian
  use gauss
  use phys_module
  use vacuum_response_module
  
  implicit none
  type (type_node_list)       :: node_list, temp_node_list
  type (type_element_list) :: element_list
  type (type_bnd_element_list):: boundary_list
  type (type_element)      :: element
  type (type_node)         :: nodes(n_vertex_max)
  
  real*8, allocatable :: xplot(:), yplot(:), yplot2(:), yplot3(:)
  
  real*8   :: ELM(n_vertex_max*(n_order+1),n_vertex_max*(n_order+1))
  real*8   :: zbig, dl, ws, psi, psi_s, v, tht, xs, ys
  real*8   :: x_g(n_gauss), x_s(n_gauss), y_g(n_gauss), y_s(n_gauss), eq_g(n_gauss), eq_s(n_gauss)
  real*8   :: am0, am1, am2, am3, am4, adrive
  real*8   :: angle, angle2, pi, xx, xxs, yy, yys, pp, ps, ps2, ps3, ss, HH(2,2), HH_s(2,2), HH_ss(2,2)
  integer  :: iplot, inode_p, inode_m, inode_p2, inode_m2, inode_p3, inode_m3, inode_p4, inode_m4
  integer  :: inode_p5, inode_m5, inode_p6, inode_m6, inode_p7, inode_m7
  integer  :: ierr, my_id, ife, iv, imode, itor, ibasis, ms
  integer  :: n_border, ilarge, n_AA, nz_AA, nz_AA_old, i,j, k,l, ibnd, jbnd, inode1, inode2, jdir
  integer  :: n_elements, index_large_i, inode, knode, index_large_k, index_ij, index_kl, index, index_i
  integer  :: vertex(2), dir(2), iv2, sms, index_basis, index_basis2, index_basis_bnd, index_basis2_bnd, kbnd, lbnd
  
  integer :: index_pert, index_resp, index_resp2, index_pert2, ibnd2, ipert, jpert
  
  PI = 2.D0 * asin(1.D0)
  
  if (my_id .eq. 0) then
  
    write(*,*) '**************************************'
    write(*,*) '*           Vacuum Poisson           *'
    write(*,*) '**************************************'
    write(*,*) ' n_elements : ',element_list%n_elements
    write(*,*) ' n_nodes    : ',node_list%n_nodes
    write(*,*) ' n_bnd      : ',boundary_list%n_bnd_elements
  
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
  
    do ibnd = 1, boundary_list%n_bnd_elements
    do jbnd = 1, 2
  
      index_basis = (ibnd-1) * 2 + jbnd - 1
  
      inode  = ibnd
      ibasis = 1 + 2*(jbnd-1)      ! only ibasis = 1 or 3 for now
  
      inode_m  = mod(ibnd + boundary_list%n_bnd_elements - 2,boundary_list%n_bnd_elements) + 1
      inode_m2 = mod(ibnd + boundary_list%n_bnd_elements - 3,boundary_list%n_bnd_elements) + 1
      inode_m3 = mod(ibnd + boundary_list%n_bnd_elements - 4,boundary_list%n_bnd_elements) + 1
      inode_m4 = mod(ibnd + boundary_list%n_bnd_elements - 5,boundary_list%n_bnd_elements) + 1
      inode_m5 = mod(ibnd + boundary_list%n_bnd_elements - 6,boundary_list%n_bnd_elements) + 1
      inode_m6 = mod(ibnd + boundary_list%n_bnd_elements - 7,boundary_list%n_bnd_elements) + 1
      inode_m7 = mod(ibnd + boundary_list%n_bnd_elements - 8,boundary_list%n_bnd_elements) + 1
      inode_p  = mod(ibnd,  boundary_list%n_bnd_elements) + 1
      inode_p2 = mod(ibnd+1,boundary_list%n_bnd_elements) + 1
      inode_p3 = mod(ibnd+2,boundary_list%n_bnd_elements) + 1
      inode_p4 = mod(ibnd+3,boundary_list%n_bnd_elements) + 1
      inode_p5 = mod(ibnd+4,boundary_list%n_bnd_elements) + 1
      inode_p6 = mod(ibnd+5,boundary_list%n_bnd_elements) + 1
      inode_p7 = mod(ibnd+6,boundary_list%n_bnd_elements) + 1
      do i=1, node_list%n_nodes
        node_list%node(i)%values = 0.d0
      enddo
  
      node_list%node(inode)%values(1,ibasis,1) = 1.d0
  !    if (jbnd .eq. 1) then
  !      node_list%node(inode)%values(1,ibasis,1)      =  1.d0
  !      node_list%node(inode_m)%values(1,ibasis,1)    =  0.84375d0
  !      node_list%node(inode_m)%values(1,ibasis+2,1)  =  0.1875d0 /2.d0
  !      node_list%node(inode_m2)%values(1,ibasis,1)   =  0.5d0
  !      node_list%node(inode_m2)%values(1,ibasis+2,1) =  0.25d0/2.d0
  !      node_list%node(inode_m3)%values(1,ibasis,1)   =  0.15625d0
  !      node_list%node(inode_m3)%values(1,ibasis+2,1) =  0.1875d0/2.d0
  
  !      node_list%node(inode_p)%values(1,ibasis,1)    =  0.84375d0
  !      node_list%node(inode_p)%values(1,ibasis+2,1)  = -0.1875d0/2.d0
  !      node_list%node(inode_p2)%values(1,ibasis,1)   =  0.5d0
  !      node_list%node(inode_p2)%values(1,ibasis+2,1) = -0.25d0/2.d0
  !      node_list%node(inode_p3)%values(1,ibasis,1)   =  0.15625d0
  !      node_list%node(inode_p3)%values(1,ibasis+2,1) = -0.1875d0/2.d0
  !    endif
  
  !    if (jbnd .eq. 1) then
  !      node_list%node(inode)%values(1,ibasis,1)      =  1.d0
  !      node_list%node(inode_m)%values(1,ibasis,1)    =  0.5d0
  !      node_list%node(inode_m)%values(1,ibasis+2,1)  =  0.25d0
  !      node_list%node(inode_p)%values(1,ibasis,1)    =  0.5d0
  !      node_list%node(inode_p)%values(1,ibasis+2,1)  = -0.25d0
  !    endif
  
  !    if (jbnd .eq. 2) then
  !      node_list%node(inode)%values(1,ibasis-2,1)    =  0.d0
  !      node_list%node(inode)%values(1,ibasis,1)      =  0.333333d0
  !      node_list%node(inode_m)%values(1,ibasis-2,1)  = -0.25d0
  !      node_list%node(inode_m)%values(1,ibasis,1)    = -0.083333d0
  !      node_list%node(inode_p)%values(1,ibasis-2,1)  =  0.25d0
  !      node_list%node(inode_p)%values(1,ibasis,1)    = -0.083333d0
  !    endif
  
  !    if (jbnd .eq. 2) then
  !      node_list%node(inode)%values(1,ibasis-2,1)     =  0.d0
  !      node_list%node(inode)%values(1,ibasis,1)       =  0.333333d0 /2.d0
  
  !      node_list%node(inode_m)%values(1,ibasis-2,1)   = -0.28125d0
  !      node_list%node(inode_m)%values(1,ibasis,1)     =  0.0625d0   /2.d0
  !      node_list%node(inode_m2)%values(1,ibasis-2,1)  = -0.25d0
  !      node_list%node(inode_m2)%values(1,ibasis,1)    = -0.083333d0 /2.d0
  !      node_list%node(inode_m3)%values(1,ibasis-2,1)  = -0.09375d0
  !      node_list%node(inode_m3)%values(1,ibasis,1)    = -0.1041666d0 /2.d0
  
  !      node_list%node(inode_p)%values(1,ibasis-2,1)   =  0.28125d0
  !      node_list%node(inode_p)%values(1,ibasis,1)     =  0.0625d0 /2.d0
  !      node_list%node(inode_p2)%values(1,ibasis-2,1)  =  0.25d0
  !      node_list%node(inode_p2)%values(1,ibasis,1)    = -0.083333d0 /2.d0
  !      node_list%node(inode_p3)%values(1,ibasis-2,1)  =  0.09375d0
  !      node_list%node(inode_p3)%values(1,ibasis,1)    = -0.1041666d0 /2.d0
  !    endif
  
  !    if (jbnd .eq. 1) then
  !      node_list%node(inode)%values(1,ibasis,1)      =  1.d0
  !      node_list%node(inode_m)%values(1,ibasis,1)    =  0.95703125d0
  !      node_list%node(inode_m)%values(1,ibasis+2,1)  =  0.109375d0 /4.d0
  !      node_list%node(inode_m2)%values(1,ibasis,1)   =  0.84375d0
  !      node_list%node(inode_m2)%values(1,ibasis+2,1) =  0.1875d0   /4.d0
  !      node_list%node(inode_m3)%values(1,ibasis,1)   =  0.68359375d0
  !      node_list%node(inode_m3)%values(1,ibasis+2,1) =  0.234375d0 /4.d0
  !      node_list%node(inode_m4)%values(1,ibasis,1)   =  0.5d0
  !      node_list%node(inode_m4)%values(1,ibasis+2,1) =  0.25d0     /4.d0
  !      node_list%node(inode_m5)%values(1,ibasis,1)   =  0.31640625d0
  !      node_list%node(inode_m5)%values(1,ibasis+2,1) =  0.234375d0 /4.d0
  !      node_list%node(inode_m6)%values(1,ibasis,1)   =  0.15625d0
  !      node_list%node(inode_m6)%values(1,ibasis+2,1) =  0.1875d0   /4.d0
  !      node_list%node(inode_m7)%values(1,ibasis,1)   =  0.04296875d0
  !      node_list%node(inode_m7)%values(1,ibasis+2,1) =  0.109375d0 /4.d0
  
  !      node_list%node(inode_p)%values(1,ibasis,1)    =  0.95703125d0
  !      node_list%node(inode_p)%values(1,ibasis+2,1)  =  -0.109375d0 /4.d0
  !      node_list%node(inode_p2)%values(1,ibasis,1)   =  0.84375d0
  !      node_list%node(inode_p2)%values(1,ibasis+2,1) =  -0.1875d0   /4.d0
  !      node_list%node(inode_p3)%values(1,ibasis,1)   =  0.68359375d0
  !      node_list%node(inode_p3)%values(1,ibasis+2,1) =  -0.234375d0 /4.d0
  !      node_list%node(inode_p4)%values(1,ibasis,1)   =  0.5d0
  !      node_list%node(inode_p4)%values(1,ibasis+2,1) =  -0.25d0     /4.d0
  !      node_list%node(inode_p5)%values(1,ibasis,1)   =  0.31640625d0
  !      node_list%node(inode_p5)%values(1,ibasis+2,1) =  -0.234375d0 /4.d0
  !      node_list%node(inode_p6)%values(1,ibasis,1)   =  0.15625d0
  !      node_list%node(inode_p6)%values(1,ibasis+2,1) =  -0.1875d0   /4.d0
  !      node_list%node(inode_p7)%values(1,ibasis,1)   =  0.04296875d0
  !      node_list%node(inode_p7)%values(1,ibasis+2,1) =  -0.109375d0 /4.d0
  !    endif
  
  !    if (jbnd .eq. 2) then
  !      node_list%node(inode)%values(1,ibasis-2,1)     =  0.d0
  !      node_list%node(inode)%values(1,ibasis,1)       =  0.33333333d0  /4.d0
  
  !      node_list%node(inode_m)%values(1,ibasis-2,1)   = -0.19140625d0
  !      node_list%node(inode_m)%values(1,ibasis,1)     =  0.18229166d0  /4.d0
  !      node_list%node(inode_m2)%values(1,ibasis-2,1)  = -0.28125d0
  !      node_list%node(inode_m2)%values(1,ibasis,1)    =  0.0625d0      /4.d0
  !      node_list%node(inode_m3)%values(1,ibasis-2,1)  = -0.29296875d0
  !      node_list%node(inode_m3)%values(1,ibasis,1)    = -0.02604166d0 /4.d0
  !      node_list%node(inode_m4)%values(1,ibasis-2,1)  = -0.25d0
  !      node_list%node(inode_m4)%values(1,ibasis,1)    = -0.083333d0    /4.d0
  !      node_list%node(inode_m5)%values(1,ibasis-2,1)  = -0.17578124d0
  !      node_list%node(inode_m5)%values(1,ibasis,1)    = -0.109375d0    /4.d0
  !      node_list%node(inode_m6)%values(1,ibasis-2,1)  = -0.09375d0
  !      node_list%node(inode_m6)%values(1,ibasis,1)    = -0.1041666d0   /4.d0
  !      node_list%node(inode_m7)%values(1,ibasis-2,1)  = -0.02734375d0
  !      node_list%node(inode_m7)%values(1,ibasis,1)    = -0.067708335d0 /4.d0
  
  !      node_list%node(inode_p)%values(1,ibasis-2,1)   =  0.19140625d0
  !      node_list%node(inode_p)%values(1,ibasis,1)     =  0.18229166d0  /4.d0
  !      node_list%node(inode_p2)%values(1,ibasis-2,1)  =  0.28125d0
  !      node_list%node(inode_p2)%values(1,ibasis,1)    =  0.0625d0      /4.d0
  !      node_list%node(inode_p3)%values(1,ibasis-2,1)  =  0.29296875d0
  !      node_list%node(inode_p3)%values(1,ibasis,1)    = -0.02604166d0 /4.d0
  !      node_list%node(inode_p4)%values(1,ibasis-2,1)  =  0.25d0
  !      node_list%node(inode_p4)%values(1,ibasis,1)    = -0.083333d0    /4.d0
  !      node_list%node(inode_p5)%values(1,ibasis-2,1)  =  0.17578124d0
  !      node_list%node(inode_p5)%values(1,ibasis,1)    = -0.109375d0    /4.d0
  !      node_list%node(inode_p6)%values(1,ibasis-2,1)  =  0.09375d0
  !      node_list%node(inode_p6)%values(1,ibasis,1)    = -0.1041666d0   /4.d0
  !      node_list%node(inode_p7)%values(1,ibasis-2,1)  =  0.02734375d0
  !      node_list%node(inode_p7)%values(1,ibasis,1)    = -0.067708335d0    /4.d0
  !    endif
  !----------------------------------------------------- test unit perturbation Fourier harmonic
  !    sms= 2.d0
  !    do inode=1, boundary_list%n_bnd_elements
  !      tht = atan2( node_list%node(inode)%x(1,2), node_list%node(inode)%x(1,1)-10.d0)
  !      node_list%node(inode)%values(1,1,1) = cos(sms*tht)
  !      node_list%node(inode)%values(1,3,1) = -1.d0 / 3.d0 * sms*sin(sms*tht) * 6.283185307 / float( boundary_list%n_bnd_elements)
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
  
              angle  = atan2(y_g(ms), x_g(ms)-R_geo)
              angle2 = MOD( angle + 2*PI, 2*PI )
              do i=1,2                                                                       ! loop over nodes
  
                do j=1,2                                                                     ! loop over basis functions
  
                  index_ij = nodes(i)%index(dir(j))                                          ! index in the ELM matrix
  
                  v   =  H1(i,j,ms) * element_list%element(ife)%size(vertex(i),dir(j))       ! test function
  
                  mumps_par%RHS(index_ij+n_AA*index_basis) = mumps_par%RHS(index_ij+n_AA*index_basis) + v * psi_s * ws              ! add to element RHS
  
  !testing_only                mumps_par%RHS(index_ij+n_AA*index_basis) = mumps_par%RHS(index_ij+n_AA*index_basis) + v * psi * ws * dl             ! add to element RHS
                enddo
              enddo
  
            enddo
  
          endif
  
        enddo
      enddo
  
    enddo
    enddo
    
  
    nz_AA_old = nz_AA
    nz_AA     = ilarge
  
    mumps_par%n  = n_AA
    mumps_par%nz = nz_AA
  
  endif

#ifdef USE_MUMPS
  mumps_par%JOB = 6
  mumps_par%SYM = 0
  mumps_par%icntl(7) = 4
  
  call DMUMPS(mumps_par)
#endif
  vacuum_response = 0.d0
  
  if (my_id .eq. 0) then
  
    do ibnd = 1, boundary_list%n_bnd_elements
  
      do jbnd = 1, 2
  
        index_basis     = (n_order+1)*(ibnd-1) + 2*(jbnd-1) + 1          ! select index 1 and 3 from the 4 dof at each node
        index_basis_bnd = 2*(ibnd-1) +   (jbnd-1) + 1                    ! the index in the vacuum_response matrix
  
        do kbnd=1,boundary_list%n_bnd_elements
  
          do lbnd = 1, 2
  
            index_basis2     = (n_order+1)*(kbnd-1) + 2*(lbnd-1) + 1
            index_basis2_bnd = 2*(kbnd-1) +   (lbnd-1) + 1               ! the index in the vacuum_response matrix
  
            vacuum_response(index_basis_bnd,index_basis2_bnd,imode) = mumps_par%rhs(index_basis2 + n_AA * (index_basis_bnd-1))
  
            write(24,'(2i5,e16.8)') index_basis_bnd, index_basis2_bnd,vacuum_response(index_basis_bnd,index_basis2_bnd,2)
  
          enddo
  
        enddo
  
  
      enddo
    enddo
  
    write(*,'(A,e14.6)') ' maxval response : ',maxval(abs(vacuum_response))
    do i=1,node_list%n_nodes
  
      inode  = i
  
      do j=1, n_order+1
  
        index  = node_list%node(inode)%index(j)
  
        node_list%node(inode)%values(1,j,1) = mumps_par%RHS(index+(boundary_list%n_bnd_elements+5)*n_AA)
  
      enddo
    enddo
  
    deallocate(mumps_par%irn,mumps_par%jcn,mumps_par%A,mumps_par%rhs)
  
    allocate(xplot(11* boundary_list%n_bnd_elements),yplot(11* boundary_list%n_bnd_elements))
    allocate(yplot2(11* boundary_list%n_bnd_elements),yplot3(11* boundary_list%n_bnd_elements))
  
  !### BEGIN DEBUG OUTPUT ###
  
  
  
    itor = 1
  
    do ipert = 1,boundary_list%n_bnd_elements
    do jpert = 1,2
  
    write(*,*) ' PERT : ',ipert,jpert
    write(22,*) ' PERT : ',ipert,jpert
  
    xplot  = 0.d0
    yplot  = 0.d0
    yplot2 = 0.d0
    yplot3 = 0.d0
    iplot  = 0
  
    index_pert = 2*(ipert-1) + jpert
  
    do ibnd = 1, boundary_list%n_bnd_elements
  
      inode1 = boundary_list%bnd_element(ibnd)%vertex(1)
      inode2 = boundary_list%bnd_element(ibnd)%vertex(2)
  
      index_resp  = 2*(ibnd-1) + 1
      ibnd2       = mod(ibnd,boundary_list%n_bnd_elements) + 1
      index_resp2 = 2*(ibnd2-1) + 1
  
  !    write(*,'(5i5,3e12.4)') inode1,inode2,index_pert,index_resp,index_resp2,vacuum_response(index_pert,index_resp,itor), &
  !                            vacuum_response2(index_pert,index_resp,itor), vacuum_response3(index_pert,index_resp,itor)
  
      do ms=1, 11
  
        xx=0.d0; xxs=0.d0; yy=0.d0; yys=0.d0; pp=0.d0; ps=0.d0; ps2 = 0.d0; ps3=0.d0
  
        ss = float(ms-1)/10.d0
  
        call basisfunctions1(ss,HH,HH_s,HH_ss)
  
        xx  = xx  + temp_node_list%node(inode1)%x(1,1) * boundary_list%bnd_element(ibnd)%size(1,1) * HH(1,1)    &
                  + temp_node_list%node(inode1)%x(3,1) * boundary_list%bnd_element(ibnd)%size(1,2) * HH(1,2)    &
                  + temp_node_list%node(inode2)%x(1,1) * boundary_list%bnd_element(ibnd)%size(2,1) * HH(2,1)    &
                  + temp_node_list%node(inode2)%x(3,1) * boundary_list%bnd_element(ibnd)%size(2,2) * HH(2,2)
        xxs = xxs + temp_node_list%node(inode1)%x(1,1) * boundary_list%bnd_element(ibnd)%size(1,1) * HH_s(1,1)  &
                  + temp_node_list%node(inode1)%x(3,1) * boundary_list%bnd_element(ibnd)%size(1,2) * HH_s(1,2)  &
                  + temp_node_list%node(inode2)%x(1,1) * boundary_list%bnd_element(ibnd)%size(2,1) * HH_s(2,1)  &
                  + temp_node_list%node(inode2)%x(3,1) * boundary_list%bnd_element(ibnd)%size(2,2) * HH_s(2,2)
        yy  = yy  + temp_node_list%node(inode1)%x(1,2) * boundary_list%bnd_element(ibnd)%size(1,1) * HH(1,1)    &
                  + temp_node_list%node(inode1)%x(3,2) * boundary_list%bnd_element(ibnd)%size(1,2) * HH(1,2)    &
                  + temp_node_list%node(inode2)%x(1,2) * boundary_list%bnd_element(ibnd)%size(2,1) * HH(2,1)    &
                  + temp_node_list%node(inode2)%x(3,2) * boundary_list%bnd_element(ibnd)%size(2,2) * HH(2,2)
        yys = yys + temp_node_list%node(inode1)%x(1,2) * boundary_list%bnd_element(ibnd)%size(1,1) * HH_s(1,1)  &
                  + temp_node_list%node(inode1)%x(3,2) * boundary_list%bnd_element(ibnd)%size(1,2) * HH_s(1,2)  &
                  + temp_node_list%node(inode2)%x(1,2) * boundary_list%bnd_element(ibnd)%size(2,1) * HH_s(2,1)  &
                  + temp_node_list%node(inode2)%x(3,2) * boundary_list%bnd_element(ibnd)%size(2,2) * HH_s(2,2)
  
        pp  = pp  + vacuum_response(index_pert,index_resp,itor)    * boundary_list%bnd_element(ibnd)%size(1,1) * HH(1,1)  &
                  + vacuum_response(index_pert,index_resp+1,itor)  * boundary_list%bnd_element(ibnd)%size(1,2) * HH(1,2)  &
                  + vacuum_response(index_pert,index_resp2,itor)   * boundary_list%bnd_element(ibnd)%size(2,1) * HH(2,1)  &
                  + vacuum_response(index_pert,index_resp2+1,itor) * boundary_list%bnd_element(ibnd)%size(2,2) * HH(2,2)
  
        ps = ps   + vacuum_response(index_pert,index_resp,itor)    * boundary_list%bnd_element(ibnd)%size(1,1) * HH_s(1,1)  &
                  + vacuum_response(index_pert,index_resp+1,itor)  * boundary_list%bnd_element(ibnd)%size(1,2) * HH_s(1,2)  &
                  + vacuum_response(index_pert,index_resp2,itor)   * boundary_list%bnd_element(ibnd)%size(2,1) * HH_s(2,1)  &
                  + vacuum_response(index_pert,index_resp2+1,itor) * boundary_list%bnd_element(ibnd)%size(2,2) * HH_s(2,2)
  
  !      ps2 = ps2 + vacuum_response2(index_pert,index_resp,itor)    * boundary_list%bnd_element(ibnd)%size(1,1) * HH(1,1)  &
  !                + vacuum_response2(index_pert,index_resp+1,itor)  * boundary_list%bnd_element(ibnd)%size(1,2) * HH(1,2)  &
  !                + vacuum_response2(index_pert,index_resp2,itor)   * boundary_list%bnd_element(ibnd)%size(2,1) * HH(2,1)  &
  !                + vacuum_response2(index_pert,index_resp2+1,itor) * boundary_list%bnd_element(ibnd)%size(2,2) * HH(2,2)
  !
  !      ps3 = ps3 + vacuum_response3(index_pert,index_resp,itor)    * boundary_list%bnd_element(ibnd)%size(1,1) * HH_s(1,1)  &
  !                + vacuum_response3(index_pert,index_resp+1,itor)  * boundary_list%bnd_element(ibnd)%size(1,2) * HH_s(1,2)  &
  !                + vacuum_response3(index_pert,index_resp2,itor)   * boundary_list%bnd_element(ibnd)%size(2,1) * HH_s(2,1)  &
  !                + vacuum_response3(index_pert,index_resp2+1,itor) * boundary_list%bnd_element(ibnd)%size(2,2) * HH_s(2,2)
  
        angle  = atan2(yy, xx-R_geo)
        angle2 = MOD( angle + 2*PI, 2*PI )
        write(22, '(9ES14.6)') ss,angle2, xx, yy, pp, ps, ps / SQRT(xxs**2 + yys**2), ps2, ps3
  
        xplot(iplot+1) = angle2
        yplot(iplot+1) = pp !/ SQRT(xxs**2 + yys**2) !ps / SQRT(xxs**2 + yys**2)
        yplot2(iplot+1) = ps2
        yplot3(iplot+1) = ps3
        iplot = iplot + 1
      enddo
    enddo
  
  !  call lincol(0)
  !  call lplot6(1,1,xplot,yplot2,iplot,'response')
  !  call lincol(3)
  !  call lplot6(1,1,xplot,yplot,-iplot,'response')
  !  call lincol(2)
  !  call lplot6(1,1,xplot,yplot3,-iplot,'response')
    enddo
    enddo
  
      xplot = 0.d0
      yplot = 0.d0
      iplot = 0
  
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
  
            do ms=1, 11
            
              xx=0.; xxs=0.; yy=0.; yys=0.; pp=0.; ps=0.
  
              ss = 1.d0 - float(ms-1)/10.d0
                    
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
              
              angle  = atan2(yy, xx-R_geo)
              angle2 = MOD( angle + 2*PI, 2*PI )
              write(22, '(8ES14.6)') ss,angle2, xx, yy, pp, ps, ps / SQRT(xxs**2 + yys**2)
  
              xplot(iplot+1) = angle2
              yplot(iplot+1) = - ps / SQRT(xxs**2 + yys**2)
              iplot = iplot + 1
  
            enddo
            
          endif
  
        enddo
      enddo
      
  
  call lincol(1)
  !call lplot6(1,1,xplot,yplot,-iplot,'response')
  !call lincol(0)
  !### END DEBUG OUTPUT ###
  
  endif
  
  return
end subroutine vacuum_poisson



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
                              + (psi_x * v_x + psi_y * v_y + (float(itor)/x_g(ms,mt))**2 * v * psi) * x_g(ms,mt)*xjac*wst
  
           enddo
         enddo
  
       enddo
     enddo
  
   enddo
  enddo
  
  return
end subroutine element_matrix_vacuum



subroutine vacuum_grid(node_list,boundary_list,r_wall,vacuum_node_list,vacuum_element_list)
  !-------------------------------------------------------------------
  ! routine calculates the vacuum reponse matrix of a vacuum region
  ! up to an ideally conducting wall
  !-------------------------------------------------------------------
  use data_structure
  
  implicit none
  
  type(type_bnd_element_list) :: boundary_list
  type(type_node_list)     :: node_list, vacuum_node_list
  type(type_element_list)  :: vacuum_element_list
  
  integer :: nr, np, i, j, k, jv, index, iv, ivp, ivm, node_iv, node_ivp, node_ivm
  real*8  :: r_wall, radius, ds, delta_Rp, delta_Zp, delta_Rm, delta_Zm, dir_2, dir_3, RZcentre(2), xlength
  
  write(*,*) '******************************'
  write(*,*) '* vacuum grid                *'
  write(*,*) '******************************'
  
  np = boundary_list%n_bnd_elements
  nr = 21
  
  write(*,*) ' nr, np : ',nr,np
  
  do j=1,np
    jv = boundary_list%bnd_element(j)%vertex(1)
    RZcentre = RZcentre + node_list%node(jv)%x(1,:)
  enddo
  RZcentre = RZcentre / float(np)
  
  do i=1,nr
  
    radius = 1.d0 + (r_wall - 1.d0) * float(i-1)/float(nr-1)
  
    do j=1,np
  
      index = (i-1)*np + j
  
      jv = boundary_list%bnd_element(j)%vertex(1)
  
      vacuum_node_list%node(index)%x(1,:) = radius * (node_list%node(jv)%x(1,:) - RZcentre) + RZcentre
      vacuum_node_list%node(index)%x(3,:) = radius *  node_list%node(jv)%x(3,:)
  
      vacuum_node_list%node(index)%x(2,:) = 1.d0/3.d0 * (r_wall - 1.d0) / float(nr-1) * (node_list%node(jv)%x(1,:) - RZcentre)
  
      vacuum_node_list%node(index)%x(4,:) = 1.d0/3.d0 * (r_wall - 1.d0) / float(nr-1) * node_list%node(jv)%x(3,:) 
  
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
end subroutine vacuum_grid
