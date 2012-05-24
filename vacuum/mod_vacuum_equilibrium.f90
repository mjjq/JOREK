!> Implements the free boundary equilibrium.
module vacuum_equilibrium
  
  use constants
  use vacuum
  
  implicit none
  
  
  
  contains
  
  
  
  !> Reads the external fields and the poloidal field coils from the STARWALL output
  subroutine import_external_fields(filename)
  
  use vacuum_response
  implicit none
  
  ! --- Routine parameters
  character(len=*), intent(in) :: filename
  ! --- Local variables
  integer, parameter   :: filehandle = 60
  integer              :: file_version, n_bnd_elems, n_bnd_nodes, dim(2)
  character(len=512)   :: comment
  
  if ( sr%i_tor(1) /= 1 ) return ! external fields not necessary in this case
  
  ! --- Read data from STARWALL response file
  open(filehandle, file=trim(filename), form='formatted', status='old', action='read')
  read(filehandle,'(a)') comment
  
  file_version = read_intparam(filehandle, 'file_version')
  if ( file_version > 1 ) then
    write(*,*) 'ERROR: COIL data file version ', file_version, ' is not supported.'
    stop
  end if
  
  n_coils     = read_intparam(filehandle, 'n_coils')
  n_bnd_nodes = read_intparam(filehandle, 'n_bnd_nodes')
  n_bnd_elems = read_intparam(filehandle, 'n_bnd_elems')
  
  dim         = (/ n_coils, 2*n_bnd_nodes /)
  
  call read_array(filehandle, 'B_t', dim, float2d=bext_tan)
  call read_array(filehandle, 'B_n', dim, float2d=bext_nor)
  call read_array(filehandle, 'Psi', dim, float2d=bext_psi)
  
  end subroutine import_external_fields
  
  
  
  
  subroutine vacuum_equil(node_list,bnd_node_list,bnd_elm_list,psi_axis,psi_bnd)
  !---------------------------------------------------------------------
  ! calculates the matrix contribution of the boundary integral of the
  ! induction equation using the vacuum response from STARWALL
  !---------------------------------------------------------------------
  use parameters
  use data_structure
  use gauss
  use basis_at_gaussian
  use phys_module
  use mumps_module
  use vacuum_response
  
  implicit none
  
  type (type_node_list)        :: node_list
  type (type_bnd_node_list)    :: bnd_node_list
  type (type_bnd_element_list) :: bnd_elm_list
  
  real*8     :: x_g(n_gauss), x_s(n_gauss)
  real*8     :: y_g(n_gauss), y_s(n_gauss)
  
  real*8     :: b_tan(n_gauss)
  
  integer    :: my_id, ibnd,i, j, ms,  kp, kbnd, k, l, jdir, kdir, ldir, imode
  integer    :: korder, kv, lv, ilarge_vv, inode, knode, inode_bnd, inode2, inode2_bnd
  integer    :: index_node, index_node2, index_node3, index_node_bnd
  integer    :: index_node2_bnd, index_node3_bnd, ilarge_pp, ijA_position
  integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, kl1, kl2, kl3, kl4, kl5, kl6, kl7
  real*8     :: ws, xjac,  BigR, phi, eps_cyl
  real*8     :: psi_axis, psi_bnd, Z_xpoint(2)
  real*8     :: rhs_glob_1, A_glob_11, A_glob_11_star, RHS_glob_11, ps0
  real*8     :: psi_norm, theta, zeta
  
  real*8     :: s, t, v, v_x, v_y, v_s, v_p, v_ss, v_xx, v_yy, v_xs, v_ys
  real*8     :: ps0_s, psi, psi_s
  real*8     :: br_coils, bz_coils, psi_coils
  
  integer    :: itmp1, itmp2, ilarge
  logical    :: xpoint2
  
  write(*,*) '**************************************************'
  write(*,*) '*     VACUUM boundary integral (equilibrium)     *'
  write(*,*) '**************************************************'
  write(*,*) 'n_coils : ',n_coils
  
  ! circular R=10 testcase
  !R_coils = (/  9.d0, 9.d0, 11.d0, 11.d0 /)
  !Z_coils = (/ -1.d0, +1.d0,  -1.d0, +1.d0 /) 
  !I_coils = (/  1d0, 1.d0, -1.d0, -1.d0 /)
  !I_coils = -0.005 * I_coils   
  
  !TS2 testacse
!   I_coils = (/ 0., 0., 0., 28800., -239040., -239040., -28800., 0., 0., 54000., 54000., 54000., 54000., 66000., 66000., 66000., 66000./)

! I_coils = (/ 0., 0., 0., 0., 0., 0., 0., 0., 0., 800000., 0., 0., 0., 800000., 0., 0., 0./)

!Limiter case
!   I_coils = (/ 0., 0., 0., 0., -440000., -440000., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0./)
! 
! vertical_FB = 100000.*vertical_FB
! 
! I_coils(10) = -vertical_FB
! I_coils(11) = -vertical_FB
! I_coils(12) = -vertical_FB
! I_coils(13) = -vertical_FB
! I_coils(14) = vertical_FB
! I_coils(15) = vertical_FB
! I_coils(16) = vertical_FB
! I_coils(17) = vertical_FB
!End limiter case

!Divertor case
  I_coils = (/ 0., 0., 0., 0., -400000., -400000., 0., 0., 0., 400000., 0., 0., 0., 400000., 0., 0., 0./)
I_coils(10) = I_coils(10)*(1.-vertical_FB)
I_coils(11) = I_coils(11)*(1.-vertical_FB)
I_coils(12) = I_coils(12)*(1.-vertical_FB)
I_coils(13) = I_coils(13)*(1.-vertical_FB)
I_coils(14) = I_coils(14)*(1.+vertical_FB)
I_coils(15) = I_coils(15)*(1.+vertical_FB)
I_coils(16) = I_coils(16)*(1.+vertical_FB)
I_coils(17) = I_coils(17)*(1.+vertical_FB)
!End divertor case

write(*,*) 'I_coils(10)', I_coils(10)

  I_coils = I_coils * 4.d-7 * PI 
  
  ilarge = mumps_par%nz                                      ! the number of items in the mumps_par coordinate storage (before the boundary conditions)
  
  do ibnd = 1, bnd_elm_list%n_bnd_elements                   ! loop over all boundary elements
  
    x_g   = 0.d0; x_s   = 0.d0;                              ! values of (x,y) and derivatives on Gaussian points
    y_g   = 0.d0; y_s   = 0.d0;
  
    b_tan = 0.d0
    
    do i=1,2                                                ! loop over two corners of each boundary element
  
      do j=1,2                                              ! loop over the two basis functions (H, H_s)
  
        do ms=1, n_gauss                                    ! loop over Gaussian points, construct coordinates and values
  
          inode     = bnd_elm_list%bnd_element(ibnd)%vertex(i)
          inode_bnd = bnd_elm_list%bnd_element(ibnd)%bnd_vertex(i)
          jdir      = bnd_elm_list%bnd_element(ibnd)%direction(i,j)
  		
          x_g(ms)  = x_g(ms)  + node_list%node(inode)%x(jdir,1) * bnd_elm_list%bnd_element(ibnd)%size(i,j) * H1(i,j,ms)         
          x_s(ms)  = x_s(ms)  + node_list%node(inode)%x(jdir,1) * bnd_elm_list%bnd_element(ibnd)%size(i,j) * H1_s(i,j,ms)    ! carefull _s is the tangential derivative here!
          y_g(ms)  = y_g(ms)  + node_list%node(inode)%x(jdir,2) * bnd_elm_list%bnd_element(ibnd)%size(i,j) * H1(i,j,ms)
          y_s(ms)  = y_s(ms)  + node_list%node(inode)%x(jdir,2) * bnd_elm_list%bnd_element(ibnd)%size(i,j) * H1_s(i,j,ms)
        
          do k=1,n_coils
  
            index_node_bnd = 2 * (bnd_node_list%bnd_node(inode_bnd)%index_starwall-1) + j
                     
  	  b_tan(ms)  = b_tan(ms) + bext_tan(index_node_bnd,k) * I_coils(k) * bnd_elm_list%bnd_element(ibnd)%size(i,j) * H1(i,j,ms)
          	
  	  enddo

        enddo                                               ! end loop over Gaussian points
      enddo                                                 ! end loop over the two basis function
    enddo                                                   ! end loop over two corners of each boundary element
  
  
    do ms=1, n_gauss                                        ! loop over Gaussian points

      ws = wgauss(ms)
  
      BigR = x_g(ms)
  
  !    call pfcoils(x_g(ms),y_g(ms),br_coils,bz_coils,psi_coils)
  
      do i=1,2                                                                 ! loop over nodes of this piece of boundary
  
        inode     = bnd_elm_list%bnd_element(ibnd)%vertex(i)
        inode_bnd = bnd_elm_list%bnd_element(ibnd)%bnd_vertex(i)
  
        do j=1,2                                                               ! loop over basis functions (test functions)
  
          jdir = bnd_elm_list%bnd_element(ibnd)%direction(i,j)
  	
          index_node = node_list%node(inode)%index(jdir)
  	
          index_node_bnd = 2 * (bnd_node_list%bnd_node(inode_bnd)%index_starwall-1) + j
  
          v =  H1(i,j,ms) * bnd_elm_list%bnd_element(ibnd)%size(i,j)             ! test function
  		
  !        mumps_par%rhs(index_node) = mumps_par%rhs(index_node) + ws * v * (bz_coils * y_s(ms) + br_coils * x_s(ms)) ! boundary integral of external vacuum field (coils)
  
          mumps_par%rhs(index_node) = mumps_par%rhs(index_node) + ws * v * b_tan(ms)* sqrt(x_s(ms)**2 + y_s(ms)**2)  ! boundary integral of external vacuum field (coils)
  
          do kbnd = 1, bnd_node_list%n_bnd_nodes                               ! kbnd really numbers the nodes of the boundary list
  
            do korder = 1, 2                                                   ! loop over basis_functions at node kbnd
                             
              kdir  = bnd_node_list%bnd_node(kbnd)%direction(korder)          ! ibnd marks the boundary element (not the node number) i.e. TO BE CHANGED
              knode = bnd_node_list%bnd_node(kbnd)%index_jorek
  
              index_node3 = node_list%node(knode)%index(kdir)
               
              !###OLD### index_node3_bnd = 2*(bnd_node_list%bnd_node(kbnd)%index_starwall-1) + korder
              index_node3_bnd = bnd_node_list%bnd_node(kbnd)%index_starwall
  
              do k=1,2                                                        ! loop over nodes in element ibnd (resulting from perturbation at node k_bnd)
  
                inode2      = bnd_elm_list%bnd_element(ibnd)%vertex(k)
                inode2_bnd  = bnd_elm_list%bnd_element(ibnd)%bnd_vertex(k)
  
                do l=1,2                                                      ! loop over basis functions
  
                  ldir        = bnd_elm_list%bnd_element(ibnd)%direction(k,l)
  
                  index_node2     = node_list%node(inode2)%index(ldir)
  		
                  !###OLD### index_node2_bnd = 2 * (bnd_node_list%bnd_node(inode2_bnd)%index_starwall-1) + l !#### to be checked
                  index_node2_bnd = bnd_node_list%bnd_node(inode2_bnd)%index_starwall
                            
                  psi   =  H1(k,l,ms)   * bnd_elm_list%bnd_element(ibnd)%size(k,l)          ! test function
                  psi_s =  H1_s(k,l,ms) * bnd_elm_list%bnd_element(ibnd)%size(k,l)          ! test function derivative (along boundary, i.e. t)
  
                  ps0   = node_list%node(bnd_elm_list%bnd_element(ibnd)%vertex(k))%values(1,ldir,1) &
  		      * H1(k,l,ms)   * bnd_elm_list%bnd_element(ibnd)%size(k,l)
                  ps0_s = node_list%node(bnd_elm_list%bnd_element(ibnd)%vertex(k))%values(1,ldir,1) &
  		      * H1_s(k,l,ms) * bnd_elm_list%bnd_element(ibnd)%size(k,l)
  
                  A_glob_11   = v * psi * sqrt(x_s(ms)**2 + y_s(ms)**2)
                  RHS_glob_11 = v * ps0 * sqrt(x_s(ms)**2 + y_s(ms)**2)
  
                  ilarge = ilarge + 1
  		
                  mumps_par%irn(ilarge) = index_node
                  mumps_par%jcn(ilarge) = index_node3
                  mumps_par%A(ilarge)   = (-1.) * (- ws * A_glob_11 * &
                  response_m_eq(response_index(index_node2_bnd,1,l), response_index(index_node3_bnd,1,korder)))
                  !###OLD### mumps_par%A(ilarge)   = - ws * A_glob_11 * vacuum_response(index_node3_bnd,index_node2_bnd,1)
  		
  ! 		  mumps_par%RHS(index_node) = mumps_par%RHS(index_node) - ws * v * RHS_glob_11 * vacuum_response(index_node3_bnd,index_node2_bnd,1) 
  								
                enddo     ! end loop over basis functions (l)
              enddo       ! end of loop elemet nodes (k)
  
            enddo         ! end loop over order (korder)
          enddo           ! end of loop over all boundary elements (kbnd)
  
        enddo             ! end loop over basis functions (j)
      enddo               ! end loop over nodes in this boundary element
  
    enddo                 ! end of loop over Gaussian points
  enddo                   ! end of loop over all boundary elements
  
  mumps_par%nz = ilarge   ! update the size of the MUMPS matrix
  
  return
  end subroutine vacuum_equil
  
   
  !**********************************************************************
  !* routines borrowed from EQUAL (WZ)                                  *
  !**********************************************************************
  
  subroutine pfcoils(R,Z,br,bz,psi)
  
  implicit none
  
  integer :: n_coils, i
  real*8  :: R, Z, br, bz, psi
  real*8,allocatable :: R_coils(:),Z_coils(:), I_coils(:)
  real*8,allocatable :: br_out(:), bz_out(:), psi_out(:), R_out(:), Z_out(:)
  
  !real*8 :: R_coils(4),Z_coils(4), I_coils(4)
  !real*8 :: br_out(4), bz_out(4), psi_out(4), R_out(4), Z_out(4)
  
  stop 'pfcoils programmed stop'
  
  n_coils = 4
  
  allocate(R_coils(n_coils),Z_coils(n_coils),I_coils(n_coils))
  allocate(R_out(n_coils),  Z_out(n_coils),  br_out(n_coils), bz_out(n_coils),psi_out(n_coils))
  
  R_coils = (/  9.d0, 9.d0, 11.d0, 11.d0 /)
  Z_coils = (/ -1.d0, +1.d0,  -1.d0, +1.d0 /) 
  I_coils = (/  1d0, 1.d0, -1.d0, -1.d0 /)
  I_coils = -0.005 * I_coils    
  
  !R_coils = (/  9.d0, 11.d0, 11.d0, 9.d0 /)
  !Z_coils = (/ -1.d0, -1.d0,  1.d0, 1.d0 /) 
  !I_coils = (/  1d0, -0.7d0, -0.7d0, 1.d0 /)
  
  Z_out(1:n_coils) = - Z + Z_coils(1:n_coils) 
  R_out(1:n_coils) = R 
  
  call brbzv(R_out,R_coils,Z_out,br_out,bz_out,n_coils)
  
  call psicalv(R_out,R_coils,Z_out,psi_out,n_coils)
  
  br  = 0.d0
  bz  = 0.d0
  psi = 0.d0
  
  do i=1, n_coils
    br  = br  + br_out(i) * I_coils(i)
    bz  = bz  + bz_out(i) * I_coils(i)
    psi = psi + psi_out(i)* I_coils(i)
  enddo
  
  return
  end subroutine pfcoils
  
  subroutine brbzv(a1,r1,z1,br,bz,n)                                
  !********************************************************************** 
  !**                                                                  ** 
  !**     MAIN PROGRAM:  MHD FITTING CODE                              ** 
  !**                                                                  ** 
  !**                                                                  ** 
  !**     SUBPROGRAM DESCRIPTION:                                      ** 
  !**          psical computes mutual inductance/2/pi between two      ** 
  !**          circular filaments of radii aa1 and r1 and              ** 
  !**          separation of z1, for mks units multiply returned       ** 
  !**          value by 2.0e-07.                                       ** 
  !**                                                                  ** 
  !**     CALLING ARGUMENTS:                                           ** 
  !**       aa1.............first filament radius                      ** 
  !**       r1..............second filament radius                     ** 
  !**       z1..............vertical separation                        ** 
  !**                                                                  ** 
  !**     REFERENCES:                                                  ** 
  !**          (1) f.w. mcclain and b.b. brown, ga technologies        ** 
  !**              report ga-a14490 (1977).                            ** 
  !**                                                                  ** 
  !**     RECORD OF MODIFICATION:                                      ** 
  !**          26/04/83..........first created                         ** 
  !**                                                                  ** 
  !**                                                                  ** 
  !**                                                                  ** 
  !********************************************************************** 
  implicit none
  real*8   ::  a1(1:n),r1(1:n),z1(1:n),br(1:n),bz(1:n)
  integer  ::  n
  
  real*8            :: z2,r2z2,den,sden,xnuma,xnum,x1,xalog,cay,ee,brkt,a2                               
  ! elliptic functions approximation (see Handbook of Mathematical Functions p.591
  real*8, parameter :: ak1=1.38629436112d0, ak2=0.09666344259d0, ak3=0.03590092383d0, ak4=0.03742563713d0, ak5=0.01451196212d0 
  real*8, parameter :: bk1=0.5d0,           bk2=0.12498593597d0, bk3=0.06880248576d0, bk4=0.03328355346d0, bk5=0.00441787012d0
                                                                   
  real*8, parameter :: ae1=0.44325141463d0, ae2=0.06260601220d0, ae3=0.04757383546d0, ae4=0.01736506451d0, &  
                       be1=0.24998368310d0, be2=0.09200180037d0, be3=0.04069697526d0, be4=0.00526449639d0                   
  integer  ::  i                                                                
  	                          
  if (n>0) then
                                                                     
     do i=1,n                        
  
         z2=z1(i)*z1(i)            
         r2z2=z2+r1(i)*r1(i)       
         den=(a1(i)+r1(i))**2+z2   
         sden=sqrt(den)      
         xnuma=(a1(i)-r1(i))**2+z2 
         xnum=max(xnuma,1d-10*den) 
         x1=xnum/den               
                                   
         xalog=-log(x1)            
                                                                       
         cay=((((x1*bk5+bk4)*x1+bk3)*x1+bk2)*x1+bk1)*xalog  &            
             +(((x1*ak5+ak4)*x1+ak3)*x1+ak2)*x1+ak1                     
                                                                       
         ee=(((be4*x1+be3)*x1+be2)*x1+be1)*x1*xalog          &           
            +(((x1*ae4+ae3)*x1+ae2)*x1+ae1)*x1+1.0d0                      
                                                                       
         brkt=ee/xnum                                                   
  !---------------------------------------------------------------------- 
  !--    br  computation                                                -- 
  !---------------------------------------------------------------------- 
         a2=a1(i)*a1(i)                                                 
         br(i)=((a2+r2z2)*brkt-cay)*z1(i)/(r1(i)*sden)                  
  !---------------------------------------------------------------------- 
  !--    bz  computation                                                -- 
  !---------------------------------------------------------------------- 
         bz(i)=((a2-r2z2)*brkt+cay)/sden                                
  
     end do          
  
  endif                                                          
                                                                         
  return                                                            
  end subroutine brbzv
  
  
  subroutine psicalv(a1,r1,z1,psi,n)
  !**********************************************************************
  !**                                                                  **
  !**     MAIN PROGRAM:  MHD FITTING CODE                              **
  !**                                                                  **
  !**                                                                  **
  !**     SUBPROGRAM DESCRIPTION:                                      **
  !**          psical computes mutual inductance/2/pi between two      **
  !**          circular filaments of radii aa1 and r1 and              **
  !**          separation of z1, for mks units multiply returned       **
  !**          value by 2.0e-07.                                       **
  !**                                                                  **
  !**     CALLING ARGUMENTS:                                           **
  !**       aa1.............first filament radius                      **
  !**       r1..............second filament radius                     **
  !**       z1..............vertical separation                        **
  !**                                                                  **
  !**     REFERENCES:                                                  **
  !**          (1) f.w. mcclain and b.b. brown, ga technologies        **
  !**              report ga-a14490 (1977).                            **
  !**                                                                  **
  !**     RECORD OF MODIFICATION:                                      **
  !**          26/04/83..........first created                         **
  !**                                                                  **
  !**                                                                  **
  !**                                                                  **
  !**********************************************************************
  implicit none
  real*8   ::  a1(1:n),r1(1:n),z1(1:n),psi(1:n)
  integer  ::  n
        
  real*8            :: z2,r2z2,xk,den,sden,xnuma,xnum,x1,xalog,cay,ee,a2 
                                
  real*8, parameter :: ak1=1.38629436112d0, ak2=0.09666344259d0, ak3=0.03590092383d0, ak4=0.03742563713d0, ak5=0.01451196212d0 
  real*8, parameter :: bk1=0.5d0,           bk2=0.12498593597d0, bk3=0.06880248576d0, bk4=0.03328355346d0, bk5=0.00441787012d0
                                                                   
  real*8, parameter :: ae1=0.44325141463, ae2=0.06260601220, ae3=0.04757383546, ae4=0.01736506451, &  
                       be1=0.24998368310, be2=0.09200180037, be3=0.04069697526, be4=0.00526449639                   
  integer  :: i
  
  psi=0.d0
  if(n>0) then
     do i=1,n
        z2=z1(i)*z1(i)
        r2z2=z2+r1(i)*r1(i)
        den=(a1(i)+r1(i))**2+z2
        xk=4.d0*a1(i)*r1(i)/den
        sden=sqrt(den)
        xnuma=(a1(i)-r1(i))**2+z2
  !WZ   xnum=max(xnuma,1d-10*den)
        xnum=max(xnuma,1d-20*den)
        x1=xnum/den
  
        xalog=-log(x1)
  
        cay=((((x1*bk5+bk4)*x1+bk3)*x1+bk2)*x1+bk1)*xalog &
            +(((x1*ak5+ak4)*x1+ak3)*x1+ak2)*x1+ak1
  
        ee=(((be4*x1+be3)*x1+be2)*x1+be1)*x1*xalog &
          +(((x1*ae4+ae3)*x1+ae2)*x1+ae1)*x1+1.d0
  !----------------------------------------------------------------------
  !--   psi computation                                                --
  !----------------------------------------------------------------------
        psi(i)= sden*((1.d0-0.5d0*xk)*cay-ee)
     end do
  endif
  
  return
  end subroutine psicalv
  
  
  
    
end module vacuum_equilibrium
