
module mod_plasma_response
  
  use constants
  use mod_parameters
  use data_structure
  use gauss
  use basis_at_gaussian
  use tr_module
  use phys_module
  
  implicit none
  
  contains
  
         
  !------------------------------------------------------------------
  !> Calculates psi_plasma at given (R,Z) points 
  !------------------------------------------------------------------
  subroutine psi_plasma(node_list,element_list,R0, Z0, psi_p)

    implicit none

    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    type (type_element)      :: element
    type (type_node)         :: nodes(n_vertex_max)
    
    real*8, intent(in)    :: R0(:), Z0(:)
    real*8, intent(inout) :: psi_p(:)

    real*8     :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
    real*8     :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
    real*8     :: eq_g(n_plane,n_gauss,n_gauss)
    
    ! --- local variables    
    integer    :: i, j, ms, mt, iv, inode, ife, mp, in
    integer    :: ierr, n_cpu, my_id, ife_delta, ife_min, ife_max, omp_nthreads, omp_tid
    real*8     :: zj0, R, Z, wst, xjac, delta_phi
    real*8     :: kk, greens_funct, Kellip_kk, Eellip_kk
    integer    :: n_points
    
    n_points = size(R0,1)
    
    psi_p     = 0.d0
    delta_phi = 2.d0 * PI / float(n_plane) / float(n_period)
        
    !--- Go through all the elements
    do ife = 1, element_list%n_elements
    
      element = element_list%element(ife)

      do iv = 1, n_vertex_max
        inode     = element%vertex(iv)
        nodes(iv) = node_list%node(inode)
      enddo
      
      x_g(:,:)    = 0.d0; x_s(:,:)    = 0.d0; x_t(:,:)    = 0.d0;
      y_g(:,:)    = 0.d0; y_s(:,:)    = 0.d0; y_t(:,:)    = 0.d0;
      
      !--- Calculate R,Z and derivatives at gausstian points
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
      
      eq_g(:,:,:) = 0.d0
            
      !--- Calculate the current at gaussian points
      do i=1,n_vertex_max
        do j=1,n_order+1

          do mp=1,n_plane
            do ms=1, n_gauss
              do mt=1, n_gauss
                do in=1,n_tor
                  eq_g(mp,ms,mt) = eq_g(mp,ms,mt) + nodes(i)%values(in,j,3) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
                enddo !---ntor
      	      enddo !---gauss
            enddo !---gauss
          enddo !---planes

        enddo !---order
      enddo !---vertex
      
      !---Do gaussian and toroidal planes integration
      do i=1, n_points
        do mp=1,n_plane
          do ms=1, n_gauss
            do mt=1, n_gauss

              wst  = wgauss(ms)*wgauss(mt)
              xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
              R    = x_g(ms,mt)
              Z    = y_g(ms,mt)
              
              !--- Calculate Green's function            
              kk           = sqrt( 4.d0*R*R0(i) / ( (R+R0(i))**2.d0 + (Z-Z0(i))**2.d0 ) ) 
              call comelp( kk, Kellip_kk, Eellip_kk)            
              greens_funct = (0.5d0/PI) * sqrt(R*R0(i))/kk * ( (2.d0-kk**2.d0)*Kellip_kk - 2.d0*Eellip_kk )
            
              zj0  = eq_g(mp,ms,mt)
            
              !--- psi = \int Greens_funct * J_phi * dA       see (4.66 Computational Methods in P.Physics, Jardin)
              psi_p(i) = psi_p(i) - zj0 / R * greens_funct *xjac * wst * delta_phi * n_period/ (2.d0 * PI)
            
            enddo
          enddo
        enddo
      enddo
      
    
    enddo !---elements

  end subroutine psi_plasma
 
 
 
  
  
  
  !---------------------------------------------------------------------
  !> Find the best coil currents for a given fixed boundary equilibrium
  !---------------------------------------------------------------------
  subroutine find_Icoils(node_list,element_list,bnd_node_list,bnd_elm_list)
    
    use vacuum
    
    type (type_node_list),       intent(in) :: node_list
    type (type_element_list),    intent(in) :: element_list
    type(type_bnd_node_list),    intent(in) :: bnd_node_list         !< List of boundary elements
    type(type_bnd_element_list), intent(in) :: bnd_elm_list          !< List of grid elements
    
    type(type_bnd_element) :: bndelem_m
    
    integer              :: n_points, n_points_elm, numb_coils
    integer              :: i, p, m_bndelem, m_pt, count, i_c, j_c
    integer              :: i_v, i_d, i_node_bnd, i_resp, i_dir, i_node, info
    real*8               :: H1(2,2), H1_s(2,2), H1_ss(2,2)
    real*8               :: s, R, Z, psi, psi_c
    real*8,  allocatable :: R_vec(:), Z_vec(:), coeff(:,:), psi_all(:), psi_p(:), psi_ext(:)
    real*8,  allocatable :: A_mat_min(:,:), RHS_min(:)
    integer, allocatable :: ipiv(:)

    write(*,*) ' '
    write(*,*) '************************************************************'
    write(*,*) '******** Calculate best I_coils from fixed bnd *************'
    write(*,*) '************************************************************'
    write(*,*) ' '
        
    !---- Calculate coefficients for each coil and point  **coeff(point,coil) **
    !------------------------------------------------------------------------
    
    numb_coils   = size(bext_tan(1,:),1)
    write(*,*) ' Coils number = ', numb_coils
    
    n_points_elm = 6
    n_points     = n_points_elm * bnd_elm_list%n_bnd_elements  
    
    allocate( R_vec(n_points), Z_vec(n_points), coeff(n_points, numb_coils) )
    allocate( psi_all(n_points), psi_ext(n_points), psi_p(n_points) )    
    R_vec = 0.d0 ; Z_vec = 0.d0; coeff = 0.d0; psi_all=0.d0; psi_ext=0.d0; psi_p=0.d0;
    
    write(*,*) ' n_bnd_elm    = ',  bnd_elm_list%n_bnd_elements
    write(*,*) ' n_points     =  ',  n_points
    
    count = 0
    
    do m_bndelem = 1, bnd_elm_list%n_bnd_elements      !--- create vector R, Z, psi
      
      bndelem_m = bnd_elm_list%bnd_element(m_bndelem)
  
      do m_pt = 1, n_points_elm     
      
        count = count + 1
        
        s = float(m_pt)/float(n_points_elm+1)

        call basisfunctions1(s, H1, H1_s, H1_ss)
  
        do i_c = 1, numb_coils
        
          R = 0.d0;  Z = 0.d0;  psi=0.d0;  psi_c = 0.d0
    
          do i_v = 1, 2
            do i_d = 1, 2 ! degrees of freedom
              i_node     = bndelem_m%vertex(i_v)
              i_node_bnd = bndelem_m%bnd_vertex(i_v)
              i_dir      = bndelem_m%direction(i_v, i_d)
              i_resp     = bnd_node_list%bnd_node(i_node_bnd)%index_starwall(i_d)
                 
              R      = R     + node_list%node(i_node)%x(i_dir,1)        * H1(i_v,i_d)   * bndelem_m%size(i_v,i_d)
              Z      = Z     + node_list%node(i_node)%x(i_dir,2)        * H1(i_v,i_d)   * bndelem_m%size(i_v,i_d)              
              psi    = psi   + node_list%node(i_node)%values(1,i_dir,1) * H1(i_v,i_d)   * bndelem_m%size(i_v,i_d)
              psi_c  = psi_c + bext_psi(i_resp,i_c)                     * H1(i_v,i_d)   * bndelem_m%size(i_v,i_d)
            end do
          end do
    
          R_vec(count) = R;      Z_vec(count) = Z;     psi_all(count) = psi;     coeff(count, i_c) = psi_c;
        enddo    
      
      enddo  
    enddo
    
    !------------------------- Calculate psi external from the fixed boundary
    !-----------------------------------------------------------------
    write(*,*) ' '
    write(*,*) 'Calculating external flux contribution of fixed bnd ...'
   
    call psi_plasma(node_list,element_list, R_vec, Z_vec, psi_p)
    
    psi_ext(:) = psi_all(:) - psi_p(:) 
    
    open(25,file='psi_contributions.txt',status="replace", position="append", action="write")
    do i = 1, n_points
      write(25,'(5ES14.6)') R_vec(i), Z_vec(i), psi_all(i), psi_p(i), psi_ext(i)
    enddo
    close(25)
    
    !------------------------- Solve minimization problem -------------------
    !------------------------------------------------------------------------
    ! Here we minimize  \sum_i ( \psi_ext - \psi_coils  )_i**2
    ! calculate RHS and A matrix
    write(*,*) ' '
    write(*,*) 'Solve minimization problem'
    
    allocate(RHS_min(numb_coils), A_mat_min(numb_coils,numb_coils) )
    allocate(ipiv(numb_coils))
    RHS_min = 0.d0;     A_mat_min = 0.d0
    
    do i_c = 1, numb_coils
      do p = 1, n_points
        RHS_min(i_c) = RHS_min(i_c) + psi_ext(p) * coeff(p, i_c)               
      enddo
      
      do j_c = 1, numb_coils
        do p = 1, n_points  
          A_mat_min(i_c, j_c) = A_mat_min(i_c, j_c) + coeff(p, i_c) * coeff(p, j_c)       
        enddo
      enddo
    enddo

   call dgesv( numb_coils, 1, A_mat_min, numb_coils, ipiv, RHS_min, numb_coils, info )
   write(*,*) 'info = ', info
   
   !----- Compare given and calculated currents
   write(*,*) ' '
   write(*,*) ' Initial coil currents, Calculated currents, Relative differences '
   do i=1, numb_coils
     write(*,'(3ES14.6)') pf_coils(i)%current, RHS_min(i), (pf_coils(i)%current - RHS_min(i)) / (pf_coils(i)%current + 0.1)
   enddo
    
   deallocate(A_mat_min,RHS_min, R_vec, Z_vec, coeff, psi_all, psi_ext, psi_p)
  
  end subroutine find_Icoils  
  
  
  
  
  
  
  ! --- Routines to calculate the elliptic integrals
  subroutine comelp ( hk, ck, ce )
  
    !*****************************************************************************80
    !
    !! COMELP computes complete elliptic integrals K(k) and E(k).
    !
    !  Licensing:
    !
    !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
    !    they give permission to incorporate this routine into a user program 
    !    provided that the copyright is acknowledged.
    !
    !  Modified:
    !
    !    07 July 2012
    !
    !  Author:
    !
    !    Shanjie Zhang, Jianming Jin
    !
    !  Reference:
    !
    !    Shanjie Zhang, Jianming Jin,
    !    Computation of Special Functions,
    !    Wiley, 1996,
    !    ISBN: 0-471-11963-6,
    !    LC: QA351.C45.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) HK, the modulus.  0 <= HK <= 1.
    !
    !    Output, real ( kind = 8 ) CK, CE, the values of K(HK) and E(HK).
    !
      implicit none
    
      real ( kind = 8 ) ae
      real ( kind = 8 ) ak
      real ( kind = 8 ) be
      real ( kind = 8 ) bk
      real ( kind = 8 ) ce
      real ( kind = 8 ) ck
      real ( kind = 8 ) hk
      real ( kind = 8 ) pk
    
      pk = 1.0D+00 - hk * hk
    
      if ( hk == 1.0D+00 ) then
    
        ck = 1.0D+300
        ce = 1.0D+00
    
      else
    
        ak = ((( &
            0.01451196212D+00   * pk &
          + 0.03742563713D+00 ) * pk &
          + 0.03590092383D+00 ) * pk &
          + 0.09666344259D+00 ) * pk &
          + 1.38629436112D+00
    
        bk = ((( &
            0.00441787012D+00   * pk &
          + 0.03328355346D+00 ) * pk &
          + 0.06880248576D+00 ) * pk &
          + 0.12498593597D+00 ) * pk &
          + 0.5D+00
    
        ck = ak - bk * log ( pk )
    
        ae = ((( &
            0.01736506451D+00   * pk &
          + 0.04757383546D+00 ) * pk &
          + 0.0626060122D+00  ) * pk &
          + 0.44325141463D+00 ) * pk &
          + 1.0D+00
    
        be = ((( &
            0.00526449639D+00   * pk &
          + 0.04069697526D+00 ) * pk &
          + 0.09200180037D+00 ) * pk &
          + 0.2499836831D+00  ) * pk
    
        ce = ae - be * log ( pk )
    
      end if
    
      return
  end  subroutine comelp

    
end module mod_plasma_response
