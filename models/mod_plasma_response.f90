!> Contains functions to calculate the fields created by the plasma current alone and by external currents
!!
!!  * The routines B_plasma and psi_plasma calculate the plasma field given a set of (R,Z) points
!!  * Routines to predict the best coil currents for given fixed-bnd equilibrium are also avaialable
!!  * Functions to calculate elliptic integrals are accessible as well
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
  
    
  !!-------------------------------------------------------------------
  !> Calculates Green's fuctions G(r,r') for B_R, B_Z and \psi 
  !!
  !! The Green's functions are used to calculate the fields by
  !! integrating current densities over areas
  !!
  !!      \psi(r) = \int G_psi(r,r') J(r') dA'
  !!-------------------------------------------------------------------
  pure subroutine Greens_functions(R, Z, R_p, Z_p, G_BR, G_BZ, G_psi)
  
    implicit none
    
    ! --- Routine parameters
    real*8,   intent(in)    :: R      !< R
    real*8,   intent(in)    :: Z      !< Z 
    real*8,   intent(in)    :: R_p    !< R'
    real*8,   intent(in)    :: Z_p    !< Z'
    real*8,   intent(inout) :: G_BR   !< Green's function for BR
    real*8,   intent(inout) :: G_BZ   !< Green's function for BZ
    real*8,   intent(inout) :: G_psi  !< Green's function for psi
    
    ! --- Local variables
    real*8   :: rho2, kk, Kellip_kk, Eellip_kk   
    
    ! --- Reference : Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop, NASA
    rho2  =  (R_p+R)**2.d0 + (Z_p-Z)**2.d0      
    kk    =  sqrt( 4.d0*R_p*R / rho2 ) 
 
    call comelp( kk, Kellip_kk, Eellip_kk)  !--- calculate elliptic functions
    
    G_BR  = (0.5d0/PI) / sqrt( rho2 ) * ( Z-Z_p ) / R                      &
          * ( (R_p**2 + R**2 + (Z_p-Z)**2) / ((R_p-R)**2.d0 + (Z_p-Z)**2.d0) * Eellip_kk - Kellip_kk )
    
    G_BZ  = (0.5d0/PI) / sqrt( rho2 )                                      &
          * ( (R_p**2 - R**2 - (Z_p-Z)**2) / ((R_p-R)**2.d0 + (Z_p-Z)**2.d0) * Eellip_kk + Kellip_kk )
    
    G_psi = (0.5d0/PI) * sqrt(R_p*R)/kk * ( (2.d0-kk**2.d0)*Kellip_kk - 2.d0*Eellip_kk )    
  
  end subroutine Greens_functions
  
  
  
  
         
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
    real*8     :: G_BR, G_BZ, G_psi
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
              call Greens_functions(R0(i), Z0(i), R, Z, G_BR, G_BZ, G_psi)
              zj0  = eq_g(mp,ms,mt)
            
              !--- psi = \int Greens_funct * J_phi * dA       see (4.66 Computational Methods in P.Physics, Jardin)
              psi_p(i) = psi_p(i) - zj0 / R * G_psi *xjac * wst * delta_phi * n_period/ (2.d0 * PI)
            
            enddo
          enddo
        enddo
      enddo
      
    
    enddo !---elements

  end subroutine psi_plasma
 
 
 
  
  
  
  
  
  !------------------------------------------------------------------
  !> Calculates B_plasma at given (R,Z) points 
  !------------------------------------------------------------------
  subroutine B_plasma(node_list,element_list,R0, Z0, B_p)

    implicit none

    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    type (type_element)      :: element
    type (type_node)         :: nodes(n_vertex_max)
    
    real*8, intent(inout)    :: R0(:), Z0(:)
    real*8, intent(inout) :: B_p(:,:)

    real*8     :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
    real*8     :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
    real*8     :: eq_g(n_plane,n_gauss,n_gauss)
    
    ! --- local variables    
    integer    :: i, j, ms, mt, iv, inode, ife, mp, in
    integer    :: ierr, n_cpu, my_id, ife_delta, ife_min, ife_max, omp_nthreads, omp_tid
    real*8     :: zj0, R, Z, wst, xjac, delta_phi
    real*8     :: G_BR, G_BZ, G_psi
    integer    :: n_points
    
    n_points = size(R0,1)
    
    B_p       = 0.d0
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
              
              !--- Calculate Green's functions (Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop, NASA)      
              call Greens_functions(R0(i), Z0(i), R, Z, G_BR, G_BZ, G_psi)  
                      
              zj0  = eq_g(mp,ms,mt)
            
              !B       =           !j_phi   !green   !dA
              B_p(i,1) = B_p(i,1) + zj0 / R * G_BR * xjac * wst * delta_phi * n_period/ (2.d0 * PI)
              B_p(i,2) = B_p(i,2) + zj0 / R * G_BZ * xjac * wst * delta_phi * n_period/ (2.d0 * PI)     
            
            enddo
          enddo
        enddo
      enddo
      
    
    enddo !---elements

  end subroutine B_plasma
  
 

  
  
    
  !---------------------------------------------------------------------
  !> Find the best coil currents for a given fixed boundary equilibrium (using psi)
  !---------------------------------------------------------------------
  subroutine find_Icoils(node_list,element_list,bnd_node_list,bnd_elm_list)
    
    use vacuum
    
    implicit none
    
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
  
  
  
  
  
  
  
  !---------------------------------------------------------------------
  !> Find the best coil currents for a given fixed boundary equilibrium (using Btan)
  !---------------------------------------------------------------------
  subroutine find_Icoils2(node_list,element_list,bnd_node_list,bnd_elm_list)
        
    use vacuum
    use vacuum_response
    
    implicit none
    
    type (type_node_list),       intent(in) :: node_list
    type (type_element_list),    intent(in) :: element_list
    type(type_bnd_node_list),    intent(in) :: bnd_node_list         !< List of boundary elements
    type(type_bnd_element_list), intent(in) :: bnd_elm_list          !< List of grid elements
  
    ! --- Local variables

    type(type_bnd_element) :: bndelem_m
    integer  :: l_starwall, l_tor
    integer  :: m_bndelem, m_pt, m_elm, mv1
    integer  :: i_vertex, i_dof, i_node, i_node_bnd, i_resp, i_resp_old, i_resp_0
    real*8   :: i_size, basfunc_i
    real*8   :: H1(2,2), H1_s(2,2), H1_ss(2,2)
    real*8   :: P, P_s, P_t, P_st, P_ss, P_tt
    real*8   :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt
    real*8   :: s_pt, t_pt, s_or_t ! s and t values at current point
    real*8   :: xjac               ! 2D Jacobian
    real*8   :: B_pol(2)           ! Poloidal magnetic field
    real*8   :: e_par(2)           ! Vector tangential to interface
    real*8   :: P_R, P_Z           ! dPsi/dR, dPsi/dZ
    real*8   :: R1, R2, Z1, Z2
    logical  :: s_const            ! Is the bound. elem. an s=const side of the 2D element?
  
    real*8,  allocatable :: R_vec(:), Z_vec(:), coeff(:,:), B_all(:,:), B_p(:,:), B_ext(:,:)
    real*8,  allocatable ::  Btan_ext(:), v_tan(:,:), weights(:)
    real*8,  allocatable :: A_mat_min(:,:), RHS_min(:)
    integer, allocatable :: ipiv(:)
    integer              :: n_points, n_points_elm, numb_coils
    integer              :: i, pt, count, i_c, j_c, info
    real*8               :: Btan_c, Bmax
  
    l_starwall=1;  l_tor=1;
  
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
    allocate( B_all(n_points,2), B_ext(n_points,2), B_p(n_points,2) )
    allocate( v_tan(n_points,2), Btan_ext(n_points), weights(n_points) )    
    R_vec = 0.d0 ; Z_vec = 0.d0; coeff = 0.d0; B_all=0.d0; B_ext=0.d0; B_p=0.d0;
    v_tan = 0.d0 ; weights=0.d0;
    
    write(*,*) ' n_bnd_elm    = ',  bnd_elm_list%n_bnd_elements
    write(*,*) ' n_points     = ',  n_points
    
    count = 0
  
    ! --- For every boundary element, do...
    L_MB: do m_bndelem = 1, bnd_elm_list%n_bnd_elements
    
      bndelem_m = bnd_elm_list%bnd_element(m_bndelem)
      m_elm     = bnd_elm_list%bnd_element(m_bndelem)%element
      mv1       = bnd_elm_list%bnd_element(m_bndelem)%side
    
      R1 = node_list%node(bndelem_m%vertex(1))%x(1,1)
      Z1 = node_list%node(bndelem_m%vertex(1))%x(1,2)
      R2 = node_list%node(bndelem_m%vertex(2))%x(1,1)
      Z2 = node_list%node(bndelem_m%vertex(2))%x(1,2)
    
      ! --- For several points in the boundary element, do...
      L_MP: do m_pt = 1, n_points_elm
    
        count = count + 1
        
        ! --- Determine 1D basis function (and derivatives) at current point
        s_or_t = float(m_pt-1)/float(n_points_elm-1)
    
        call basisfunctions1(s_or_t, H1, H1_s, H1_ss)
    
        ! --- Which s and t values correspond to the current point and is the
        !     boundary element an s=const or t=const side of the 2D element?
        select case (mv1)
        case (1)
          s_pt = s_or_t;  t_pt = 0.d0;    s_const = .false.
        case (2)
          s_pt = 1.d0;    t_pt = s_or_t;  s_const = .true.
        case (3)
          s_pt = s_or_t;  t_pt = 1.d0;    s_const = .false.
        case (4)
          s_pt = 0.d0;    t_pt = s_or_t;  s_const = .true.
        end select
    
        ! --- Determine coordinate values (plus derivatives)
        call interp_RZ(node_list, element_list, m_elm, s_pt, t_pt, R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt)
    
        ! --- 2D Jacobian
        xjac = R_s * Z_t - R_t * Z_s
    
        ! --- Tangential vector to the interface
        if ( s_const ) then
          e_par = (/ R_t, Z_t /) / sqrt( R_t**2 + Z_t**2 ) * (R_t * (R2-R1) + Z_t * (Z2-Z1))/abs(R_t * (R2-R1) + Z_t * (Z2-Z1))
        else
          e_par = (/ R_s, Z_s /) / sqrt( R_s**2 + Z_s**2 ) * (R_s * (R2-R1) + Z_s * (Z2-Z1))/abs(R_s * (R2-R1) + Z_s * (Z2-Z1))
        end if
    
        ! --- Psi value (plus derivatives) at current point (l_tor mode)
        call interp(node_list, element_list, m_elm, 1, l_tor, s_pt, t_pt, P, P_s, P_t, P_st, P_ss, P_tt)
    
        ! --- Poloidal magnetic field at current point
        P_R   = (   P_s * Z_t - P_t * Z_s ) / xjac ! dPsi/dR
        P_Z   = ( - P_s * R_t + P_t * R_s ) / xjac ! dPsi/dZ
        B_pol = (/ P_Z, -P_R /) / R
        
        do i_c = 1, numb_coils
        
          Btan_c = 0.d0
    
          ! --- Sum over boundary dofs at which response is calculated
          L_IV: do i_vertex = 1, 2 ! (loop over nodes in element m_bndelem)
    
            i_node_bnd  = bndelem_m%bnd_vertex(i_vertex)
    
            L_ID: do i_dof = 1, 2 ! (loop over node dofs)
    
              i_size   = bndelem_m%size(i_vertex,i_dof)                      
              i_resp_0 = response_index_eq(i_node_bnd,i_dof)
        
              ! --- Determine basis function
              basfunc_i = H1(i_vertex,i_dof) * i_size
            
              ! --- Btan per coil
              Btan_c = Btan_c + basfunc_i * bext_tan(i_resp_0, i_c)
          
            end do L_ID
          end do L_IV
                      
          coeff(count, i_c) = Btan_c;
        
        enddo  !coils loop
        
        R_vec(count)   = R;     Z_vec(count) = Z;       B_all(count,:) = B_pol(:);   v_tan(count,:) = e_par(:); 
    
      end do L_MP
    
    end do L_MB
    
    
    !------------------------- Calculate B external from the fixed boundary
    !-----------------------------------------------------------------
    write(*,*) ' '
    write(*,*) 'Calculating external field contribution of fixed bnd ...'
   
    call B_plasma(node_list,element_list, R_vec, Z_vec, B_p)
    
    B_ext(:,:)  = B_all(:,:) - B_p(:,:) 
    Btan_ext(:) = B_ext(:,1)*v_tan(:,1) +  B_ext(:,2)*v_tan(:,2)  
    
    !--- Check coils initial guess
    open(25,file='B_ext.txt',status="replace", position="append", action="write")
    do i = 1, n_points
      write(25,'(5ES14.6)') R_vec(i), Z_vec(i), Btan_ext(i), sum(coeff(i,:) * I_coils(:))
    enddo
    close(25)
    
    Bmax = maxval(abs(Btan_ext))
    !--- Calculate weights for each point
    do i = 1, n_points
      if ( abs(Btan_ext(i)) < 0.1*Bmax) then 
        weights(i) = 10.d0
      else 
        weights(i) = Bmax / abs(Btan_ext(i))
      endif
    enddo
   
    !------------------------- Solve minimization problem -------------------
    !------------------------------------------------------------------------
    ! Here we minimize  \sum_i ( \B_ext - \B_coils  )_i**2
    ! calculate RHS and A matrix
    write(*,*) ' '
    write(*,*) 'Solve minimization problem'
    
    allocate(RHS_min(numb_coils), A_mat_min(numb_coils,numb_coils) )
    allocate(ipiv(numb_coils))
    RHS_min = 0.d0;     A_mat_min = 0.d0
    
    do i_c = 1, numb_coils
      do pt = 1, n_points
        RHS_min(i_c) = RHS_min(i_c) + Btan_ext(pt) * coeff(pt, i_c) * weights(pt)              
      enddo
      
      do j_c = 1, numb_coils
        do pt = 1, n_points  
          A_mat_min(i_c, j_c) = A_mat_min(i_c, j_c) + coeff(pt, i_c) * coeff(pt, j_c) * weights(pt)     
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

  end subroutine find_Icoils2  
  
  
  
  
  
  
  !!-------------------------------------------------------------------------------------------
  !> Find the best coil currents for a given fixed boundary equilibrium (using Btan)  for JET
  !! This is a separate routine as the 20 JET coils are constrained to 10 circuits and therefore
  !! 10 degrees of freedom. TO DO: A generalization for circuits and coils dof should be done
  !!-------------------------------------------------------------------------------------------
  subroutine find_Icoils_JET(node_list,element_list,bnd_node_list,bnd_elm_list)
    
    use vacuum
    use vacuum_response
    
    type (type_node_list),       intent(in) :: node_list
    type (type_element_list),    intent(in) :: element_list
    type(type_bnd_node_list),    intent(in) :: bnd_node_list         !< List of boundary elements
    type(type_bnd_element_list), intent(in) :: bnd_elm_list          !< List of grid elements
  
    ! --- Local variables

    type(type_bnd_element) :: bndelem_m
    integer  :: l_starwall, l_tor
    integer  :: m_bndelem, m_pt, m_elm, mv1
    integer  :: i_vertex, i_dof, i_node, i_node_bnd, i_resp, i_resp_old, i_resp_0
    real*8   :: i_size, basfunc_i
    real*8   :: H1(2,2), H1_s(2,2), H1_ss(2,2)
    real*8   :: P, P_s, P_t, P_st, P_ss, P_tt
    real*8   :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt
    real*8   :: s_pt, t_pt, s_or_t ! s and t values at current point
    real*8   :: xjac               ! 2D Jacobian
    real*8   :: B_pol(2)           ! Poloidal magnetic field
    real*8   :: e_par(2)           ! Vector tangential to interface
    real*8   :: P_R, P_Z           ! dPsi/dR, dPsi/dZ
    real*8   :: R1, R2, Z1, Z2
    logical  :: s_const            ! Is the bound. elem. an s=const side of the 2D element?
  
    real*8,  allocatable :: R_vec(:), Z_vec(:), coeff(:,:), B_all(:,:), B_p(:,:), B_ext(:,:)
    real*8,  allocatable ::  Btan_ext(:), v_tan(:,:), weights(:)
    real*8,  allocatable :: A_mat_min(:,:), RHS_min(:)
    integer, allocatable :: ipiv(:), alpha(:,:)
    integer              :: n_points, n_points_elm, numb_coils, n_eq
    integer              :: i, pt, count, i_c, j_c, info, i_eq_i, i_eq_j
    real*8               :: Btan_c, Bmax, cont_ip, cont_jp
  
    l_starwall=1;  l_tor=1;
  
    write(*,*) ' '
    write(*,*) '************************************************************'
    write(*,*) '****** Calculate best JET I_coils from fixed bnd ***********'
    write(*,*) '************************************************************'
    write(*,*) ' '
        
    !---- Calculate coefficients for each coil and point  **coeff(point,coil) **
    !------------------------------------------------------------------------
    
    numb_coils   = size(bext_tan(1,:),1)
    
    if (numb_coils /= 20) then
      write(*,*) 'For using this feature with JET you must use the 20 coils generated by write_jet_pfsystems'
      stop
    endif
    write(*,*) ' Coils number = ', numb_coils
    
    n_points_elm = 6
    n_points     = n_points_elm * bnd_elm_list%n_bnd_elements  
    
    allocate( R_vec(n_points), Z_vec(n_points), coeff(n_points, numb_coils) )
    allocate( B_all(n_points,2), B_ext(n_points,2), B_p(n_points,2) )
    allocate( v_tan(n_points,2), Btan_ext(n_points), weights(n_points) )    
    R_vec = 0.d0 ; Z_vec = 0.d0; coeff = 0.d0; B_all=0.d0; B_ext=0.d0; B_p=0.d0;
    v_tan = 0.d0 ; weights=0.d0;
    
    write(*,*) ' n_bnd_elm    = ',  bnd_elm_list%n_bnd_elements
    write(*,*) ' n_points     = ',  n_points
    
    count = 0
  
      ! --- For every boundary element, do...
    L_MB: do m_bndelem = 1, bnd_elm_list%n_bnd_elements
    
      bndelem_m = bnd_elm_list%bnd_element(m_bndelem)
      m_elm     = bnd_elm_list%bnd_element(m_bndelem)%element
      mv1       = bnd_elm_list%bnd_element(m_bndelem)%side
    
      R1 = node_list%node(bndelem_m%vertex(1))%x(1,1)
      Z1 = node_list%node(bndelem_m%vertex(1))%x(1,2)
      R2 = node_list%node(bndelem_m%vertex(2))%x(1,1)
      Z2 = node_list%node(bndelem_m%vertex(2))%x(1,2)
    
      ! --- For several points in the boundary element, do...
      L_MP: do m_pt = 1, n_points_elm
    
        count = count + 1
        
        ! --- Determine 1D basis function (and derivatives) at current point
        s_or_t = float(m_pt-1)/float(n_points_elm-1)
    
        call basisfunctions1(s_or_t, H1, H1_s, H1_ss)
    
        ! --- Which s and t values correspond to the current point and is the
        !     boundary element an s=const or t=const side of the 2D element?
        select case (mv1)
        case (1)
          s_pt = s_or_t;  t_pt = 0.d0;    s_const = .false.
        case (2)
          s_pt = 1.d0;    t_pt = s_or_t;  s_const = .true.
        case (3)
          s_pt = s_or_t;  t_pt = 1.d0;    s_const = .false.
        case (4)
          s_pt = 0.d0;    t_pt = s_or_t;  s_const = .true.
        end select
    
        ! --- Determine coordinate values (plus derivatives)
        call interp_RZ(node_list, element_list, m_elm, s_pt, t_pt, R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt)
    
        ! --- 2D Jacobian
        xjac = R_s * Z_t - R_t * Z_s
    
        ! --- Tangential vector to the interface
        if ( s_const ) then
          e_par = (/ R_t, Z_t /) / sqrt( R_t**2 + Z_t**2 ) * (R_t * (R2-R1) + Z_t * (Z2-Z1))/abs(R_t * (R2-R1) + Z_t * (Z2-Z1))
        else
          e_par = (/ R_s, Z_s /) / sqrt( R_s**2 + Z_s**2 ) * (R_s * (R2-R1) + Z_s * (Z2-Z1))/abs(R_s * (R2-R1) + Z_s * (Z2-Z1))
        end if
    
        ! --- Psi value (plus derivatives) at current point (l_tor mode)
        call interp(node_list, element_list, m_elm, 1, l_tor, s_pt, t_pt, P, P_s, P_t, P_st, P_ss, P_tt)
    
        ! --- Poloidal magnetic field at current point
        P_R   = (   P_s * Z_t - P_t * Z_s ) / xjac ! dPsi/dR
        P_Z   = ( - P_s * R_t + P_t * R_s ) / xjac ! dPsi/dZ
        B_pol = (/ P_Z, -P_R /) / R
        
        do i_c = 1, numb_coils
        
          Btan_c = 0.d0
    
          ! --- Sum over boundary dofs at which response is calculated
          L_IV: do i_vertex = 1, 2 ! (loop over nodes in element m_bndelem)
    
            i_node_bnd  = bndelem_m%bnd_vertex(i_vertex)
    
            L_ID: do i_dof = 1, 2 ! (loop over node dofs)
    
              i_size   = bndelem_m%size(i_vertex,i_dof)                      
              i_resp_0 = response_index_eq(i_node_bnd,i_dof)
        
              ! --- Determine basis function
              basfunc_i = H1(i_vertex,i_dof) * i_size
            
              ! --- Btan per coil
              Btan_c = Btan_c + basfunc_i * bext_tan(i_resp_0, i_c)
          
            end do L_ID
          end do L_IV
                      
          coeff(count, i_c) = Btan_c;
        
        enddo  !coils loop
        
        R_vec(count)   = R;     Z_vec(count) = Z;       B_all(count,:) = B_pol(:);   v_tan(count,:) = e_par(:); 
    
      end do L_MP
    
    end do L_MB
    
    
    !------------------------- Calculate B external from the fixed boundary
    !-----------------------------------------------------------------
    write(*,*) ' '
    write(*,*) 'Calculating external field contribution of fixed bnd ...'
   
    call B_plasma(node_list,element_list, R_vec, Z_vec, B_p)
    
    B_ext(:,:)  = B_all(:,:) - B_p(:,:) 
    Btan_ext(:) = B_ext(:,1)*v_tan(:,1) +  B_ext(:,2)*v_tan(:,2)  
    
    !--- Check coils initial guess
    open(25,file='B_ext.txt',status="replace", position="append", action="write")
    do i = 1, n_points
      write(25,'(5ES14.6)') R_vec(i), Z_vec(i), Btan_ext(i), sum(coeff(i,:) * I_coils(:))
    enddo
    close(25)
    
    Bmax = maxval(abs(Btan_ext))
    !--- Calculate weights for each point
    do i = 1, n_points
      if ( abs(Btan_ext(i)) < 0.1*Bmax) then 
        weights(i) = 10.d0
      else 
        weights(i) = Bmax / abs(Btan_ext(i))
      endif
    enddo
   
    !------------------------- Solve minimization problem -------------------
    !------------------------------------------------------------------------
    ! Here we minimize  \sum_i ( \B_ext - \B_coils  )_i**2
    ! calculate RHS and A matrix
    write(*,*) ' '
    write(*,*) 'Solve minimization problem'
    
    n_eq = 10
    
    allocate(RHS_min(n_eq), A_mat_min(n_eq,n_eq) )
    allocate(ipiv(n_eq))
    allocate(alpha(n_eq,numb_coils))
    RHS_min = 0.d0;     A_mat_min = 0.d0
    
    !--- Define relations between circuits and coils, multiply by number of turns
    alpha = 0.d0   !first index = circuit,  second = coil
    alpha(1, 1)   = 710.d0
    alpha(2, 2)   = 426.d0
    alpha(3, 3)   =  -8.d0
    alpha(3, 4)   = -20.d0
    alpha(3, 5)   =  -8.d0
    alpha(3, 6)   = -20.d0
    alpha(4, 7)   =  -8.d0
    alpha(4, 8)   =   8.d0
    alpha(3, 9)   =  30.d0
    alpha(3,10)   =  30.d0
    alpha(4,11)   = -20.d0
    alpha(4,12)   =  20.d0
    alpha(1,13)   =   2.d0
    alpha(1,14)   =   2.d0
    alpha(5,15)   =  61.d0
    alpha(5,16)   =  61.d0
    alpha(6,15)   = -61.d0
    alpha(6,16)   =  61.d0
    alpha(7,17)   =  15.99d0
    alpha(8,18)   =  15.d0
    alpha(9,19)   =  15.d0
    alpha(10,20)  =  21.d0
    
    !--- Calculate B_j 
    do i_eq_j = 1, n_eq
      do pt = 1, n_points
        do i_c = 1, numb_coils        
          RHS_min(i_eq_j) = RHS_min(i_eq_j) + Btan_ext(pt) * coeff(pt, i_c) * weights(pt) * alpha(i_eq_j, i_c)             
        enddo
      enddo
    enddo
    
    !--- Calculate A_ij
    do i_eq_i = 1, n_eq          ! loops for A matrix main indeces
      do i_eq_j = 1, n_eq   
        
        do pt = 1, n_points 
          
          cont_ip = 0.d0; cont_jp = 0.d0
          
          do i_c = 1, numb_coils           
            cont_ip = cont_ip + alpha(i_eq_i, i_c) * coeff(pt, i_c) 
            cont_jp = cont_jp + alpha(i_eq_j, i_c) * coeff(pt, i_c)  
          enddo
          
          A_mat_min(i_eq_i, i_eq_j) = A_mat_min(i_eq_i, i_eq_j) + cont_ip * cont_jp * weights(pt)
        
        enddo
        
      enddo
    enddo

   call dgesv( n_eq, 1, A_mat_min, n_eq, ipiv, RHS_min, n_eq, info )
   write(*,*) 'info = ', info
   
   !----- Compare given and calculated currents
   write(*,*) ' '
   write(*,*) ' Found circuit currents '
   do i=1, n_eq
     write(*,'(1ES14.6)') RHS_min(i)
   enddo
   
   
   write(*,*) ' '
   write(*,*) ' Found total coil currents '
   do i=1, numb_coils
     I_coils(i) =  sum( alpha(:,i) * RHS_min(:) )
     write(*,*) I_coils(i)
   enddo
   
   !--- Check coils initial guess
    open(25,file='B_ext_new.txt',status="replace", position="append", action="write")
    do i = 1, n_points
      write(25,'(4ES14.6)') R_vec(i), Z_vec(i), Btan_ext(i),  sum(coeff(i,:) * I_coils(:)) 
    enddo
    close(25)

  end subroutine find_Icoils_JET
  
  
  
  
  
  
  
  ! --- Routine to calculate the elliptic integrals
  pure subroutine comelp ( hk, ck, ce )
  
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
    
      real*8, intent(in)    :: hk
      real*8, intent(inout) :: ck
      real*8, intent(inout) :: ce
      
      real ( kind = 8 ) ae
      real ( kind = 8 ) ak
      real ( kind = 8 ) be
      real ( kind = 8 ) bk
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
