!> Compares \f$B_{tan}\f$ and \f$B_{tan,vacuum}\f$ at the interface.
!!
!! The relative difference is calculated according to:
!! \f$\frac{\int dA |B_{tan,plasma}-B_{tan,vacuum}|}{\int dA |B_{tan,vacuum}|}\f$.
subroutine boundary_check()
  
  use tr_module
  use data_structure,  only: type_bnd_element 
  use phys_module,     only: resistive_wall
  use nodes_elements,  only: node_list, bnd_node_list, element_list, bnd_elm_list
  use vacuum_response, only: NEW_VACUUM, n_dof_starwall, n_starwall_harmonics, starwall_harmonics, &
    response_index, starwall_m_ee, starwall_m_ey, starwall_m_id, wall_curr, det_psibnd_vec
  
  implicit none
  
  ! --- Local variables
  integer, parameter     :: N_POINTS = 11 ! Number of evaluation points per element
  type(type_bnd_element) :: bndelem_m
  real*8, allocatable    :: B_par(:), B_par_v(:)
  real*8, allocatable    :: val_integral(:), err_integral(:)
  real*8, allocatable    :: psibnd_vec(:), dpsibnd_vec(:)
  integer  :: l_starwall, l_tor
  integer  :: m_bndelem, m_pt, m_elm, mv1
  integer  :: i_vertex, i_dof, i_node, i_node_bnd, i_resp
  real*8   :: i_size, basfunc_i
  real*8   :: H1(2,2), H1_s(2,2), H1_ss(2,2)
  real*8   :: P, P_s, P_t, P_st, P_ss, P_tt
  real*8   :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt
  real*8   :: s_pt, t_pt, s_or_t ! s and t values at current point
  real*8   :: xjac               ! 2D Jacobian
  real*8   :: B_pol(2)           ! Poloidal magnetic field
  real*8   :: e_par(2)           ! Vector tangential to interface
  real*8   :: P_R, P_Z           ! dPsi/dR, dPsi/dZ
  logical  :: s_const            ! Is the bound. elem. an s=const side of the 2D element?
  
  if ( .not. NEW_VACUUM ) return ! (boundary check for old vacuum is not implemented anymore)
  
  write(*,*) '************************************'
  write(*,*) '*    check boundary conditions     *'
  write(*,*) '************************************'
  
  call tr_allocate(psibnd_vec,1,n_dof_starwall,"psibnd_vec",CAT_GRID)
  call tr_allocate(dpsibnd_vec,1,n_dof_starwall,"dpsibnd_vec",CAT_GRID)
  call tr_allocate(B_par,1,n_starwall_harmonics,"B_par",CAT_GRID)
  call tr_allocate(B_par_v,1,n_starwall_harmonics,"B_par_v",CAT_GRID)
  call tr_allocate(val_integral,1,n_starwall_harmonics,"val_integral",CAT_GRID)
  call tr_allocate(err_integral,1,n_starwall_harmonics,"err_integral",CAT_GRID)
  
  ! --- Determine vectors with the Psi and deltaPsi values at the boundary.
  call det_psibnd_vec(bnd_node_list, node_list, psibnd_vec, dpsibnd_vec)
  
  val_integral(:) = 0.d0
  err_integral(:) = 0.d0
  
  ! --- For every boundary element, do...
  L_MB: do m_bndelem = 1, bnd_elm_list%n_bnd_elements
    bndelem_m = bnd_elm_list%bnd_element(m_bndelem)
    m_elm     = bnd_elm_list%bnd_element(m_bndelem)%element
    mv1       = bnd_elm_list%bnd_element(m_bndelem)%side
    
    ! --- For several points in the boundary element, do...
    L_MP: do m_pt = 1, N_POINTS + 1
      
      B_par_v(:) = 0.d0
      
      ! --- Determine 1D basis function (and derivatives) at current point
      s_or_t = float(m_pt-1)/float(N_POINTS)
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
      call interp_RZ(node_list, element_list, m_elm, s_pt, t_pt, R, R_s, R_t, R_st, R_ss, R_tt, Z, &
        Z_s, Z_t, Z_st, Z_ss, Z_tt)
      
      ! --- 2D Jacobian
      xjac = R_s * Z_t - R_t * Z_s
      
      ! --- Tangential vector to the interface
      if ( s_const ) then
        e_par = (/ R_t, Z_t /) / sqrt( R_t**2 + Z_t**2 )
      else
        e_par = (/ R_s, Z_s /) / sqrt( R_s**2 + Z_s**2 )
      end if
      
      ! --- Select one STARWALL harmonic
      L_LS: do l_starwall = 1, n_starwall_harmonics
        l_tor = starwall_harmonics(l_starwall)
        
        ! --- Psi value (plus derivatives) at current point (l_tor mode)
        call interp(node_list, element_list, m_elm, 1, l_tor, s_pt, t_pt, P, P_s, P_t, P_st, P_ss, &
          P_tt)
        
        ! --- Poloidal magnetic field at current point
        P_R   = (   P_s * Z_t - P_t * Z_s ) / xjac ! dPsi/dR
        P_Z   = ( - P_s * R_t + P_t * R_s ) / xjac ! dPsi/dZ
        B_pol = (/ P_Z, -P_R /) / R
        
        ! --- Tangential magnetic field B_{||} reconstructed from the plasma
        B_par(l_starwall) = - sum( B_pol * e_par )
        
        ! --- Sum over boundary dofs at which response is calculated
        L_IV: do i_vertex = 1, 2 ! (loop over nodes in element m_bndelem)
          L_ID: do i_dof = 1, 2 ! (loop over node dofs)
            i_node      = bndelem_m%vertex(i_vertex)
            i_node_bnd  = bndelem_m%bnd_vertex(i_vertex)
            i_size      = bndelem_m%size(i_vertex,i_dof)
            i_resp      = response_index(i_node_bnd,l_starwall,i_dof)
            
            ! --- Determine basis function
            basfunc_i = H1(i_vertex,i_dof) * i_size
            
            ! --- Determine B_{||,v} as prescribed by the vacuum.
            if ( resistive_wall ) then
              B_par_v(l_starwall) = B_par_v(l_starwall) + basfunc_i * (     &
                + sum( starwall_m_ee(i_resp, :) * psibnd_vec(:) )           &
                + sum( starwall_m_ey(i_resp, :) * wall_curr(:)  ) )
            else
              B_par_v(l_starwall) = B_par_v(l_starwall) + basfunc_i         &
                * sum( starwall_m_id(i_resp, :) * psibnd_vec(:) )
            end if
            
          end do L_ID
        end do L_IV
        
      end do L_LS
      
      !### DEBUG OUTPUT ###
      !write(88,'(20ES15.5)') (m_bndelem-1 + s_or_t)/REAL(bnd_elm_list%n_bnd_elements), B_par(:)
      !write(89,'(20ES15.5)') (m_bndelem-1 + s_or_t)/REAL(bnd_elm_list%n_bnd_elements), B_par_v(:)
      !###
      
      ! --- Integration of B_par_v values and differences between B_par and B_par_v.
      val_integral(:) = val_integral(:) + abs( B_par_v(:) )
      err_integral(:) = err_integral(:) + abs( B_par(:) - B_par_v(:) )
      
    end do L_MP
    
  end do L_MB
  
  if ( minval(abs(val_integral)) /= 0.d0 ) then
    write(*,'(A,20ES15.5)') 'relative errors in harmonics:', err_integral(:) / val_integral(:)
  end if
  
  !### DEBUG OUTPUT ###
  !write(88,*)
  !write(88,*)
  !write(89,*)
  !write(89,*)
  !if ( minval(abs(val_integral)) /= 0.d0 ) then ! (avoid division by zero in first timestep)
  !  write(87,'(20ES15.5)') err_integral(:) / val_integral(:)
  !end if
  !###
  
  call tr_deallocate(psibnd_vec,"psibnd_vec",CAT_GRID)
  call tr_deallocate(dpsibnd_vec,"dpsibnd_vec",CAT_GRID)
  call tr_deallocate(B_par,"B_par",CAT_GRID)
  call tr_deallocate(B_par_v,"B_par_v",CAT_GRID)
  call tr_deallocate(val_integral,"val_integral",CAT_GRID)
  call tr_deallocate(err_integral,"err_integral",CAT_GRID)
  
end subroutine boundary_check
