!> Implements the interplay of the plasma with a conducting wall.
!!
!! The plasma-wall interaction is characterized by vacuum response matrices which are calculated by
!! the STARWALL code and imported into JOREK by the routine read_starwall_response(). The vacuum
!! response enters into the boundary integral of the current equation which is implemented in the
!! routine vacuum_boundary_integral().
!!
!! IMPORTANT: The variable s in a boundary element does not necessarily correspond to s in the
!! 2D elements depending on element orientation. 
module vacuum_response
  
  use vacuum
  
  implicit none
  
  integer, parameter :: ivar_psi = 1 !< Index of Psi variable  @todo: put into mod_parameters later
  integer, parameter :: ivar_j   = 3 !< Index of j variable    @todo: put into mod_parameters later
  
  character(len=7), parameter :: LOG_BEGIN_MARKER = 'BEGIN: '
  character(len=5), parameter :: LOG_END_MARKER   = 'END: '
  
  
  
  
  
  
  contains
  
  
  
  
  
  
  !> Determines the vacuum response for an ideal or resistive wall
  subroutine get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list,                    &
    freeboundary_equil, use_starwall, resistive_wall)
  
    use parameters,     only: n_tor
    use data_structure, only: type_node_list, type_bnd_element_list, type_bnd_node_list
    
    implicit none
      
    include 'mpif.h'
    
    integer,                     intent(in) :: my_id              !< MPI proc ID
    type(type_node_list),        intent(in) :: node_list          !< List of boundary nodes
    type(type_bnd_element_list), intent(in) :: bnd_elm_list       !< List of boundary elements
    type(type_bnd_node_list),    intent(in) :: bnd_node_list      !< List of boundary nodes
    logical,                     intent(in) :: freeboundary_equil !< Use free boundary equilibrium?
    logical,                     intent(in) :: use_starwall       !< Use STARWALL response?
    logical,                     intent(in) :: resistive_wall     !< Resistive or ideal wall?

    integer :: ierr, dim
    
    
    n_dof_bnd     = bnd_node_list%n_bnd_nodes * 2 ! Number of boundary degrees of freedom per harmonic
    
    write(*,*) LOG_BEGIN_MARKER//' get_vacuum_response'
    
    ! --- Output some information about the boundary.
    230 format(A,' = ',I8)
    231 format(A,' = ',L8)
    write(*,230) 'n_bnd_elements', bnd_elm_list%n_bnd_elements
    write(*,230) 'n_bnd_nodes   ', bnd_node_list%n_bnd_nodes
    write(*,230) 'n_dof_bnd     ', n_dof_bnd
    write(*,231) 'use_starwall  ', use_starwall
    if ( use_starwall ) write(*,231) 'resistive_wall', resistive_wall
    
    ! --- Write out the boundary information for STARWALL.
    if (my_id .eq. 0) call export_boundary(node_list, bnd_elm_list, bnd_node_list)
    
    ! *************************** NEW IMPLEMENTATION ***********************************
    if ( NEW_VACUUM ) then
      if ( my_id == 0 ) then
        if ( use_starwall ) then
          call read_starwall_response(freeboundary_equil, resistive_wall)
        else
          !### JOREK response currently not implemented ###
        end if
      end if
      call broadcast_starwall_response(my_id, resistive_wall)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    else
    ! *************************** OLD IMPLEMENTATION ***********************************
      ! --- Resistive wall
      if ( resistive_wall ) then
        
        ! --- Get the STARWALL response matrices
        call resistive_wall_starwall(my_id,node_list,bnd_elm_list,bnd_node_list)
      
      ! --- Ideal wall 
      else
        
        ! --- Allocate the ideal wall vacuum response matrix
        if ( allocated(vac_response) ) deallocate(vac_response)
        dim = response_index(bnd_node_list%n_bnd_nodes,n_tor,2)
        allocate( vac_response(dim,dim) )
        vac_response = 0.
        
        ! --- Get the STARWALL response matrix
        if ( use_starwall ) then
          call ideal_wall_starwall(my_id,node_list,bnd_elm_list,bnd_node_list)
          
        ! --- Let JOREK determine the response (works only in special cases; for testing)
        else
          call ideal_wall(my_id,node_list,bnd_elm_list,bnd_node_list)
          
        end if
        
        ! --- Send the vacuum response matrix.
        call MPI_bcast(vac_response, response_index(bnd_node_list%n_bnd_nodes,n_tor,2)**2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        
      end if
    end if
    ! *********************************************************************************
      
    write(*,*) LOG_END_MARKER//' get_vacuum_response'
    
  end subroutine get_vacuum_response
  
  
  
  
  
  
  
  !> Read the STARWALL response matrices from files.
  subroutine read_starwall_response(freeboundary_equil, resistive_wall)
  
    use parameters,     only: n_tor, n_period
    
    implicit none
    
    ! --- Routine parameters
    logical,                     intent(in) :: freeboundary_equil !< Free boundary equilibrium?
    logical,                     intent(in) :: resistive_wall     !< Read which kind of response
    
    ! --- Local variables
    integer           :: dim(2)
    real*8            :: TWOPI
    integer           :: i, i_starw, n, is_sin
    
    write(*,*) LOG_BEGIN_MARKER//' Reading STARWALL response matrices'
    TWOPI = 8.d0 * atan(1.d0)
    
    ! --- Determine STARWALL harmonics
    if ( vacuum_debug ) write(*,*) LOG_BEGIN_MARKER//' Determining STARWALL harmonics'
    open(42, file='starwall_harmonics', status='old', action='read')
    read(42,*) n_starwall_harmonics
    if ( allocated(starwall_harmonics) ) deallocate( starwall_harmonics )
    allocate( starwall_harmonics(n_starwall_harmonics) )
    read(42,*) starwall_harmonics
    close(42)
    n_dof_starwall = n_dof_bnd * n_starwall_harmonics
    22 format(1x,a,100i4)
    write(*,22) 'n_starwall_harmonics =', n_starwall_harmonics
    write(*,22) '=> n_dof_starwall    =', n_dof_starwall
    write(*,*) 'starwall_harmonics   = '//trim(modes_to_str(starwall_harmonics,n_starwall_harmonics,1))
    if ( vacuum_debug ) write(*,*) LOG_END_MARKER//' Determining STARWALL harmonics'
    
    ! --- Transform STARWALL harmonics to account for periodicity
    if ( vacuum_debug ) write(*,*) LOG_BEGIN_MARKER//' Transforming STARWALL harmonics (JOREK periodicity)'
    do i = 1, n_starwall_harmonics
      i_starw = starwall_harmonics(i)
      n       = i_starw / 2
      is_sin  = i_starw - 2 * n
      i_starw = 2 * n/n_period + is_sin
      if ( (mod(n, n_period) /= 0) .or. (i_starw < 1) .or. (i_starw > n_tor) ) then
        write(*,*) 'ERROR: STARWALL harmonic has no JOREK equivalent!'
        write(*,*) 'i_starw    =', starwall_harmonics(i)
        write(*,*) 'n_period   =', n_period 
        write(*,*) 'n_tor      =', n_tor
        stop
      end if
      starwall_harmonics(i) = i_starw
    end do
    write(*,*) 'starwall_harmonics   = '//trim(modes_to_str(starwall_harmonics,n_starwall_harmonics,n_period))
    if ( vacuum_debug ) write(*,*) LOG_END_MARKER//' Transforming STARWALL harmonics (JOREK periodicity)'
    
    ! --- Check if i_tor=1 is provided as required for freeboundary equilibrium
    if ( freeboundary_equil .and. (starwall_harmonics(1) /= 1) ) then
      write(*,*) 'ERROR: STARWALL response does not include i_tor=1 mode'
      write(*,*) '  which is required for the freeboundary equilibrium.'
      stop
    end if
    
    ! --- Read the response matrices.
    if ( vacuum_debug ) write(*,*) LOG_BEGIN_MARKER//' Reading matrices'
    if ( resistive_wall ) then
      write(*,*) 'Reading resistive wall matrices.'
      ! (Determine number of wall currents)
      write(*,*) LOG_BEGIN_MARKER//' Determining number of wall currents'
      open(42, FILE='starwall_d_yy', status='old', action='read')
      read(42,*) n_wall_curr
      write(*,'(1x,a,i7)') 'n_wall_curr=', n_wall_curr
      close(42)
      write(*,*) LOG_END_MARKER//' Determining number of wall currents'
      ! (EE matrix)
      dim(:) = n_dof_starwall
      call read_response_matrix( starwall_m_ee, dim, 'starwall_m_ee' )
      ! (EY matrix)
      dim(:) = (/ n_dof_starwall, n_wall_curr /)
      call read_response_matrix( starwall_m_ey, dim, 'starwall_m_ey' )
      ! (YE matrix)
      dim(:) = (/ n_wall_curr, n_dof_starwall /)
      call read_response_matrix( starwall_m_ye, dim, 'starwall_m_ye' )
      ! (YY diagonal matrix)
      call read_response_diagonal( starwall_d_yy, n_wall_curr, 'starwall_d_yy' )
    else
      write(*,*) 'Reading ideal wall matrix.'
      dim(:) = n_dof_starwall
      call read_response_matrix( starwall_m_id, dim, 'starwall_m_id' )
      call read_response_matrix( starwall_m_nw, dim, 'starwall_m_nw' )
    end if
    if ( vacuum_debug ) write(*,*) LOG_END_MARKER//' Reading matrices'
    
    ! --- Perform import normalizations.
    if ( vacuum_debug ) write(*,*) LOG_BEGIN_MARKER//' Performing import normalizations'
    if ( resistive_wall ) then
      starwall_m_ee(:,:) = starwall_m_ee(:,:) * TWOPI
      starwall_m_ey(:,:) = starwall_m_ey(:,:) * TWOPI
    else
      starwall_m_id(:,:) = starwall_m_id(:,:) * TWOPI
      starwall_m_nw(:,:) = starwall_m_nw(:,:) * TWOPI
    end if
    if ( vacuum_debug ) write(*,*) LOG_END_MARKER//' Performing import normalizations'
    
    ! --- DEBUGGING: Output checksums
    if ( vacuum_debug ) then
      write(*,*) LOG_BEGIN_MARKER//' Cecksums'
      if ( resistive_wall ) then
        write(*,*) 'ee:', sum(abs(starwall_m_ee))
        write(*,*) 'ey:', sum(abs(starwall_m_ey))
        write(*,*) 'ye:', sum(abs(starwall_m_ye))
        write(*,*) 'yy:', sum(abs(starwall_d_yy))      
      else
        write(*,*) 'id:', sum(abs(starwall_m_id))
        write(*,*) 'nw:', sum(abs(starwall_m_nw))
      end if
      write(*,*) LOG_END_MARKER//' Cecksums'
    end if
    
    write(*,*) LOG_END_MARKER//' Reading STARWALL response matrices'
    
  end subroutine read_starwall_response
  
  
  
  
  
  
  !> Broadcast the STARWALL response matrices
  subroutine broadcast_starwall_response(my_id, resistive_wall)
    
    implicit none
    
    include 'mpif.h'
    
    ! --- Routine parameters
    integer, intent(in) :: my_id          !< MPI proc ID
    logical, intent(in) :: resistive_wall !< Broadcast which response matrices
    
    ! --- Local parameters
    integer :: ierr
    real*8  :: checksum
    
    if ( vacuum_debug ) write(*,*) LOG_BEGIN_MARKER//' Broadcasting STARWALL response matrices'
    
    ! --- Broadcast matrices and parameters.
    call MPI_bcast(n_dof_starwall, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(n_starwall_harmonics, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    if ( .not. allocated(starwall_harmonics) ) allocate( starwall_harmonics(n_starwall_harmonics) )
    
    call MPI_bcast(starwall_harmonics, n_starwall_harmonics, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    if ( resistive_wall ) then
      call MPI_bcast(n_wall_curr,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      
      if ( .not. allocated(starwall_m_ee) ) allocate( starwall_m_ee(n_dof_starwall,n_dof_starwall) )
      if ( .not. allocated(starwall_m_ey) ) allocate( starwall_m_ey(n_dof_starwall,n_wall_curr) )
      if ( .not. allocated(starwall_m_ye) ) allocate( starwall_m_ye(n_wall_curr,n_dof_starwall) )
      if ( .not. allocated(starwall_d_yy) ) allocate( starwall_d_yy(n_wall_curr) )
      
      call MPI_bcast(starwall_m_ee, n_dof_starwall**2,          MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(starwall_m_ey, n_dof_starwall*n_wall_curr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(starwall_m_ye, n_wall_curr*n_dof_starwall, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(starwall_d_yy, n_wall_curr,                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    else
      if ( .not. allocated(starwall_m_id) ) allocate( starwall_m_id(n_dof_starwall,n_dof_starwall) )
      if ( .not. allocated(starwall_m_nw) ) allocate( starwall_m_nw(n_dof_starwall,n_dof_starwall) )
      
      call MPI_bcast(starwall_m_id, n_dof_starwall**2,          MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(starwall_m_nw, n_dof_starwall**2,          MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    end if
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    
    ! --- Output a debugging checksum (must be the same on all MPI procs)
    if ( vacuum_debug ) then
      checksum = n_dof_starwall + n_starwall_harmonics + sum(starwall_harmonics)
      if ( resistive_wall ) then
        checksum = checksum + n_wall_curr + sum(abs(starwall_m_ee)) + sum(abs(starwall_m_ey)) &
          + sum(abs(starwall_m_ye)) + sum(abs(starwall_d_yy)) 
      else
        checksum = checksum + sum(abs(starwall_m_id)) + sum(abs(starwall_m_nw))
      end if
      write(*,*) my_id, 'debug checksum:', checksum
    end if
    
    if ( vacuum_debug ) write(*,*) LOG_END_MARKER//' Broadcasting STARWALL response matrices'
  
  end subroutine broadcast_starwall_response
  
  
  
  
  
  
  !> Implement the vacuum boundary integral in the current equation (ideal or resistive wall).
  subroutine vacuum_boundary_integral(bnd_node_list, node_list, bnd_elm_list,              &
    freeboundary_equil, resistive_wall, index_min, index_max, rhs_loc, tstep)
      
    use data_structure, only: type_node_list, type_bnd_node_list, type_bnd_element_list, type_bnd_element
    use parameters,     only: n_plane, n_var, n_tor
    use gauss,          only: n_gauss, xgauss, wgauss
    use global_distributed_matrix, only: irn_glob, jcn_glob, a_glob, ndof_glob, det_row_col, det_sparse_pos
    use basis_at_gaussian, only: H1, H1_s, HZ
    
    implicit none
    
    ! --- Routine parameters
    type(type_node_list),        intent(in)    :: node_list            !< List of grid nodes
    type(type_bnd_node_list),    intent(in)    :: bnd_node_list        !< List of boundary grid nodes
    type(type_bnd_element_list), intent(in)    :: bnd_elm_list         !< List of boundary elements
    logical,                     intent(in)    :: freeboundary_equil   !< Use free boundary equilibrium?
    logical,                     intent(in)    :: resistive_wall       !< Resistive or ideal wall?
    integer,                     intent(in)    :: index_min, index_max !< Responsibility of MPI proc
    real*8,                      intent(inout) :: rhs_loc(ndof_glob)   !< Part of RHS of MPI proc
    real*8,                      intent(in)    :: tstep                !< delta t, timestep
    
    ! --- Local variables
    real*8, allocatable :: psibnd_vec(:)    ! Vector of the values of Psi at the boundary
    real*8, allocatable :: dpsibnd_vec(:)   ! Vector of the values of deltaPsi at the boundary
    real*8   :: ws                          ! Gauss weight for integration in s direction
    real*8   :: amat_contrib, rhs_contrib   ! Vacuum response contribution to lhs and rhs
    real*8   :: testfunc_l                  ! j^*_l in documentation
    real*8   :: basfunc_i                   ! b_i in documentation
    real*8   :: sqrt_xs2_ys2                ! factor from definition of dA
    real*8   :: x_s(n_gauss), y_s(n_gauss)  ! values of dR/ds and dZ/ds at Gaussian points
    integer  :: m_bndelem                   ! Boundary element index
    type(type_bnd_element) :: bndelem_m     ! Boundary element corresponding to index m_bndelem
    integer  :: ms                          ! Gauss point index
    integer  :: m_plane                     ! Toroidal plane index
    integer  :: sparsepos_jp, sparsepos_pp  ! Position of lhs contribution in the sparse matrix
    !   --- Test function related quantities
    integer  :: l_vertex, l_dof, l_dir, l_node, l_node_bnd, l_index, l_tor, l_row_j, l_row_psi
    real*8   :: l_size
    !   --- Quantities related to the boundary dof at which response is calculated
    integer  :: i_vertex, i_dof, i_dir, i_node, i_node_bnd, i_index, i_starwall, i_tor, i_resp
    real*8   :: i_size
    !   --- Quantities related to the boundary dof contributing to the response
    integer  :: j_dof, j_dir, j_node, j_node_bnd, j_index, j_starwall, j_tor, j_col_psi, j_resp
    
    !integer :: rate, t0, t1 !### timing ###
    
    write(*,*) LOG_BEGIN_MARKER//' vacuum_boundary_integral'
    
    if ( vacuum_debug ) then
      write(*,*) LOG_BEGIN_MARKER//' Cecksums'
      write(*,*) 'rhs_loc: ', sum(abs(rhs_loc))
      write(*,*) 'A_glob:  ', sum(abs(A_glob))
      write(*,*) 'irn_glob:', sum(abs(irn_glob))
      write(*,*) 'jcn_glob:', sum(abs(jcn_glob))
      write(*,*) LOG_END_MARKER//' Cecksums'
    end if
    
    ! --- Determine vectors of the psi and deltapsi boundary values.
    call det_psibnd_vec(bnd_node_list, node_list, psibnd_vec, dpsibnd_vec)
    
    ! --- Update the derived response matrices
    call update_response(tstep, freeboundary_equil, resistive_wall)
    
    ! --- Perform the time-stepping for the wall currents.
    if ( resistive_wall ) call evolve_wall_currents(psibnd_vec, dpsibnd_vec)
    
    call boundary_check()
    
    !write(35+my_id,'(4ES20.12)') sum(abs(rhs_loc)), sum(abs(A_glob))
    
    !### timing ###
    !call system_clock(count_rate=rate)
    !call system_clock(count=t0)
    !###
    
    ! --- Sum over boundary elements
    !$omp parallel do                                                                              &
    !$omp default(shared)                                                                          &
    !$omp private(m_bndelem, bndelem_m, x_s, y_s, l_vertex, l_dof, l_node, l_dir, l_node_bnd,      &
    !$omp   l_index, l_size, l_tor, l_row_j, l_row_psi, ms, ws, sqrt_xs2_ys2, m_plane,             &
    !$omp   testfunc_l, i_vertex, i_dof, i_node, i_dir, i_node_bnd, i_index, i_size, i_starwall,   &
    !$omp   i_tor, i_resp, basfunc_i, j_node_bnd, j_dof, j_node, j_dir, j_index, j_starwall,       &
    !$omp   j_tor, j_resp, j_col_psi, sparsepos_jp, sparsepos_pp, amat_contrib, rhs_contrib)
    L_MB: do m_bndelem = 1, bnd_elm_list%n_bnd_elements
      bndelem_m = bnd_elm_list%bnd_element(m_bndelem)
      
      ! --- Determine the values of R,s and Z,s at the Gaussian points.
      call det_coord_bnd(bndelem_m, node_list, R_S=x_s, Z_S=y_s)
      
      ! --- Select a test function (the weak form equation must hold for every test function)
      L_LV: do l_vertex = 1, 2 ! (loop over nodes in element m_bndelem)
        L_LD: do l_dof = 1, 2 ! (loop over node dofs)
          l_node      = bndelem_m%vertex(l_vertex)
          l_dir       = bndelem_m%direction(l_vertex,l_dof)
          l_node_bnd  = bndelem_m%bnd_vertex(l_vertex)
          l_index     = node_list%node(l_node)%index(l_dir)
          l_size      = bndelem_m%size(l_vertex,l_dof)
          if ( (l_index < index_min) .or. (l_index > index_max) ) cycle ! This MPI proc responsible?
          L_LS: do l_tor = 1, n_tor ! (loop over toroidal harmonics)
            
            ! --- Determine the row in the main matrix.
            l_row_psi = det_row_col(l_index, ivar_psi, l_tor)
            l_row_j   = det_row_col(l_index, ivar_j,   l_tor)
            
            ! --- Loop over Gaussian points -- integration in s-direction
            L_MS: do ms = 1, n_gauss
              ws = wgauss(ms)
              
              ! --- Integration factor from the definition of dA:
              !     int dA = sum_{m_bndelem} int ds int dphi sqrt{(R,s)^2 + (Z,s)^2}
              sqrt_xs2_ys2 = sqrt(x_s(ms)**2 + y_s(ms)**2)
              
              ! --- Loop over toroidal planes -- integration in phi-direction
              L_MP: do m_plane = 1, n_plane
                
                ! --- Evaluate test function at current position
                testfunc_l = H1(l_vertex,l_dof,ms) * l_size * HZ(l_tor,m_plane)
                
                ! --- Sum over boundary dofs at which response is calculated
                L_IV: do i_vertex = 1, 2 ! (loop over nodes in element m_bndelem)
                  L_ID: do i_dof = 1, 2 ! (loop over node dofs)
                    i_node      = bndelem_m%vertex(i_vertex)
                    i_dir       = bndelem_m%direction(i_vertex,i_dof)
                    i_node_bnd  = bndelem_m%bnd_vertex(i_vertex)
                    i_index     = node_list%node(i_node)%index(i_dir)
                    i_size      = bndelem_m%size(i_vertex,i_dof)
                    L_IS: do i_starwall = 1, n_starwall_harmonics ! (loop over STARWALL harmonics)
                      i_tor  = starwall_harmonics(i_starwall)
                      i_resp = response_index(i_node_bnd,i_starwall,i_dof)
            
                      ! --- Determine basis function
                      basfunc_i = H1(i_vertex,i_dof,ms) *i_size * HZ(i_tor,m_plane)
                      
                      ! --- Sum over boundary dofs contributing to the response
                      L_JB: do j_node_bnd = 1, bnd_node_list%n_bnd_nodes ! (loop over boundary nodes)
                        L_JD: do j_dof = 1, 2 ! (loop over node dofs)
                          j_node      = bnd_node_list%bnd_node(j_node_bnd)%index_jorek
                          j_dir       = bnd_node_list%bnd_node(j_node_bnd)%direction(j_dof)
                          j_index     = node_list%node(j_node)%index(j_dir)
                          L_JS: do j_starwall = 1, n_starwall_harmonics ! (loop over STARWALL harmonics)
                            j_tor  = starwall_harmonics(j_starwall)
                            j_resp = response_index(j_node_bnd,j_starwall,j_dof)
                            
                            ! --- Option to switch off mode coupling due to a 3D wall
                            if ( vacuum_decouple_modes .and. (j_tor /= i_tor) ) cycle
                            
                            ! --- Determine the column in the main matrix
                            j_col_psi = det_row_col(j_index, ivar_psi, j_tor)
                            
                            ! --- Determine the position in the sparse matrix data structure
                            !     which corresponds to the matrix entry at l_row_j, j_col_psi.
                            sparsepos_jp = det_sparse_pos(l_row_j,   j_col_psi, index_min)
                            sparsepos_pp = det_sparse_pos(l_row_psi, j_col_psi, index_min)
              
                            ! --- Vacuum response contribution to the lhs of the current equation
                            amat_contrib = testfunc_l * ws  &
                              * basfunc_i * sqrt_xs2_ys2 * response_m_pp_lhs(i_resp, j_resp)
                            !$omp atomic
                            A_glob(sparsepos_jp) = A_glob(sparsepos_jp) + amat_contrib

                          end do L_JS
                        end do L_JD
                      end do L_JB
                      
                      ! --- Contribution of vacuum response to the rhs of the current equation
                      rhs_contrib = sum( response_m_pp_rhs(i_resp, :) * psibnd_vec(:) )
                      if ( resistive_wall ) rhs_contrib = rhs_contrib                            &
                        + sum( response_m_pw_rhs(i_resp, :) * wall_curr(:) )
                      rhs_contrib = rhs_contrib * testfunc_l * ws * basfunc_i * sqrt_xs2_ys2
                      !$omp atomic
                      rhs_loc(l_row_j) = rhs_loc(l_row_j) + rhs_contrib
                      
                    end do L_IS
                  end do L_ID
                end do L_IV
                
              end do L_MP
              
            end do L_MS
            
          end do L_LS
        end do L_LD
      end do L_LV
      
    end do L_MB
    !$omp end parallel do
    
    !### timing ###
    !call system_clock(count=t1)
    !write(*,*) 'vacuum_boundary_integral main loop:', real(t1 - t0 ) / real(rate), 's'
    !write(68+my_id,*) real(t1 - t0 ) / real(rate)
    !###
    
    if ( vacuum_debug ) then
      write(*,*) LOG_BEGIN_MARKER//' Cecksums'
      write(*,*) 'rhs_loc: ', sum(abs(rhs_loc))
      write(*,*) 'A_glob:  ', sum(abs(A_glob))
      write(*,*) 'irn_glob:', sum(abs(irn_glob))
      write(*,*) 'jcn_glob:', sum(abs(jcn_glob))
      write(*,*) LOG_END_MARKER//' Cecksums'
    end if
  
    if ( allocated(psibnd_vec ) ) deallocate( psibnd_vec  )
    if ( allocated(dpsibnd_vec) ) deallocate( dpsibnd_vec )
    
    write(*,*) LOG_END_MARKER//' vacuum_boundary_integral'
    
  end subroutine vacuum_boundary_integral
  
  
  
  
  
  
  !> Determine coordinate values on Gaussian points in a given boundary element.
  subroutine det_coord_bnd(bndelem, node_list, R, Z, R_s, Z_s)
    
    use gauss,             only: n_gauss, xgauss, wgauss
    use data_structure,    only: type_node, type_bnd_element, type_node_list
    use basis_at_gaussian, only: H1, H1_s, HZ
    
    implicit none
    
    ! --- Routine parameters
    type(type_bnd_element), intent(in)  :: bndelem       !< Boundary element to be considered.
    type(type_node_list),   intent(in)  :: node_list     !< List of grid nodes
    real*8, optional,       intent(out) :: R(n_gauss)    !< Values of R
    real*8, optional,       intent(out) :: Z(n_gauss)    !< Values of Z
    real*8, optional,       intent(out) :: R_s(n_gauss)  !< Values of R,s
    real*8, optional,       intent(out) :: Z_s(n_gauss)  !< Values of Z,s
    
    ! --- Local variables
    integer         :: k_vertex, k_dof, k_node, k_dir
    real*8          :: k_size
    type(type_node) :: node_k
    
    if ( present(R  ) ) R   = 0.d0
    if ( present(Z  ) ) Z   = 0.d0
    if ( present(R_s) ) R_s = 0.d0
    if ( present(Z_s) ) Z_s = 0.d0
    
    do k_vertex = 1, 2
      do k_dof = 1, 2
        k_node      = bndelem%vertex(k_vertex)
        k_dir       = bndelem%direction(k_vertex,k_dof)
        k_size      = bndelem%size(k_vertex,k_dof)
        node_k      = node_list%node(k_node)
        if ( present(R  ) ) R  (:)  = R  (:)  + node_k%x(k_dir,1) * k_size * H1  (k_vertex,k_dof,:)
        if ( present(Z  ) ) Z  (:)  = Z  (:)  + node_k%x(k_dir,2) * k_size * H1  (k_vertex,k_dof,:)
        if ( present(R_s) ) R_s(:)  = R_s(:)  + node_k%x(k_dir,1) * k_size * H1_s(k_vertex,k_dof,:)
        if ( present(Z_s) ) Z_s(:)  = Z_s(:)  + node_k%x(k_dir,2) * k_size * H1_s(k_vertex,k_dof,:)
      end do
    end do
    
  end subroutine det_coord_bnd
  
  
  
  
  
  
  !> Determine vectors of the psi and deltapsi values at the boundary.
  subroutine det_psibnd_vec(bnd_node_list, node_list, psibnd_vec, dpsibnd_vec)
    
    use data_structure, only: type_node_list, type_bnd_node_list
    
    implicit none
    
    ! --- Routine parameters
    type(type_node_list),     intent(in)  :: node_list      !< List of grid nodes
    type(type_bnd_node_list), intent(in)  :: bnd_node_list  !< List of boundary grid nodes
    real*8, allocatable,      intent(out) :: psibnd_vec(:)  !< Vector of psi boundary values
    real*8, allocatable,      intent(out) :: dpsibnd_vec(:) !< Vector of deltapsi boundary values
    
    ! --- Local variables
    integer :: jnode, jnode_glob, j_starwall, jtor, jbas, jdir, j_resp
    
    if ( vacuum_debug ) write(*,*) LOG_BEGIN_MARKER//' det_psibnd_vec'
    
    if ( allocated(psibnd_vec) ) deallocate(psibnd_vec)
    allocate( psibnd_vec(n_dof_starwall) )
    if ( allocated(dpsibnd_vec) ) deallocate(dpsibnd_vec)
    allocate( dpsibnd_vec(n_dof_starwall) )
    
    ! --- Determine vector of (delta)psi boundary values.
    do jnode = 1, bnd_node_list%n_bnd_nodes       ! loop over nodes
      jnode_glob = bnd_node_list%bnd_node(jnode)%index_jorek
      do j_starwall = 1, n_starwall_harmonics     ! loop over STARWALL harmonics
        jtor = starwall_harmonics(j_starwall)     ! (mode corresponding to STARWALL harmonic)
        do jbas = 1, 2                            ! loop over basis functions
          jdir   = bnd_node_list%bnd_node(jnode)%direction(jbas)
          j_resp = response_index(jnode,j_starwall,jbas)
          
          psibnd_vec ( j_resp ) = node_list%node(jnode_glob)%values(jtor, jdir, ivar_psi)
          dpsibnd_vec( j_resp ) = node_list%node(jnode_glob)%deltas(jtor, jdir, ivar_psi)
          
        end do
      end do
    end do
    
    if ( vacuum_debug ) write(*,*) LOG_END_MARKER//' det_psibnd_vec'
    
  end subroutine det_psibnd_vec
  
  
  
  
  
  
  !> Implement the wall current time evolution.
  subroutine evolve_wall_currents(psibnd_vec, dpsibnd_vec)
    
    implicit none
    
    ! --- Routine parameters
    real*8,  intent(in) :: psibnd_vec (n_dof_starwall) !< Vector of psi boundary values
    real*8,  intent(in) :: dpsibnd_vec(n_dof_starwall) !< Vector of deltapsi boundary values
    
    ! --- Local variables
    integer :: i
    
    if ( vacuum_debug ) write(*,*) LOG_BEGIN_MARKER//' evolve_wall_currents'
    
    ! --- Prepare the wall current array on first use.
    if ( .not. allocated(wall_curr) ) then
      allocate( wall_curr(n_wall_curr) )
      wall_curr = 0.d0
    end if
    
    ! --- Influence of wall current on itself (resistive damping)
    wall_curr(:) = response_d_ww(:) * wall_curr(:)
    
    ! --- Influence of plasma on wall currents (induction)
    do i = 1, n_wall_curr
      if ( vacuum_implicit ) then
        wall_curr(i) = wall_curr(i) + sum( response_m_wp(i,:) * dpsibnd_vec(:) )
      else
        wall_curr(i) = wall_curr(i) + sum( response_m_wp(i,:) * psibnd_vec(:)  )
      end if
    end do
    
    if ( vacuum_debug ) write(*,*) LOG_END_MARKER//' evolve_wall_currents'
    
  end subroutine evolve_wall_currents
  
  
  
  
  
  
  !> Update the derived response matrices.
  subroutine update_response(tstep, freeboundary_equil, resistive_wall)
    
    implicit none
    
    ! --- Routine parameters
    real*8,                      intent(in) :: tstep              !< delta t, timestep
    logical,                     intent(in) :: freeboundary_equil !< Use free boundary equilibrium?
    logical,                     intent(in) :: resistive_wall     !< Resistive or ideal wall?
    
    ! --- Local variables
    integer :: i, j
    real*8  :: a, b
    logical :: update_required
    
    ! --- Local variables to store the previous values of some parameters.
    real*8,  save :: old_thick
    real*8,  save :: old_res
    real*8,  save :: old_tstep
    logical, save :: old_reswall
    
    write(*,*) LOG_BEGIN_MARKER//' update_response'
    
    ! --- Update of response matrices is required only, if these parameter values changed (or the matrices are not allocated)
    update_required = ( old_thick   /= wall_thickness      ) &
                 .or. ( old_res     /= wall_resistivity    ) &
                 .or. ( old_tstep   /= tstep               ) &
                 .or. ( old_reswall .neqv. resistive_wall  ) &
                 .or. (                       .not. allocated(response_m_pp_rhs)  ) &
                 .or. (                       .not. allocated(response_m_pp_lhs)  ) &
                 .or. ( resistive_wall .and. (.not. allocated(response_d_ww))     ) &
                 .or. ( resistive_wall .and. (.not. allocated(response_m_pw_rhs)) ) &
                 .or. ( resistive_wall .and. (.not. allocated(response_m_wp))     ) &
                 .or. (                       .not. allocated(response_m_eq)      )
    
    if ( update_required ) then ! (Do nothing, if derived matrices need no update)
    
      ! --- Remember parameter values.
      old_thick   = wall_thickness
      old_res     = wall_resistivity
      old_tstep   = tstep
      old_reswall = resistive_wall
      
      ! --- Allocate matrices if required
      if ( .not. allocated(response_m_eq    ) ) &
        allocate( response_m_eq(n_dof_starwall, n_dof_starwall) )
      if ( .not. allocated(response_m_pp_rhs) ) &
        allocate( response_m_pp_rhs(n_dof_starwall, n_dof_starwall) )
      if ( .not. allocated(response_m_pp_lhs) ) &
        allocate( response_m_pp_lhs(n_dof_starwall, n_dof_starwall) )
      if ( resistive_wall ) then
        if ( .not. allocated(response_d_ww)     ) &
          allocate( response_d_ww(n_wall_curr) )
        if ( .not. allocated(response_m_pw_rhs) ) &
          allocate( response_m_pw_rhs(n_dof_starwall, n_wall_curr) )
        if ( .not. allocated(response_m_wp)     ) &
          allocate( response_m_wp(n_wall_curr, n_dof_starwall) )
      end if
      
      ! --- Derived response matrix for equilibrium
      if ( resistive_wall ) then
        response_m_eq(:,:) = starwall_m_ee(:,:)
      else
        response_m_eq(:,:) = starwall_m_nw(:,:)
      end if
      
      ! --- Derived response matrices for time-evolution
      
      !   --- Resistive wall, implicit wall-current time evolution
      if ( resistive_wall .and. vacuum_implicit ) then
        
        response_d_ww(:)       = 1.d0 - 1.d0 / &
          ( 1.d0 + wall_thickness / ( wall_resistivity * tstep * starwall_d_yy(:) ) )
        
        do j = 1, n_dof_starwall
          response_m_wp(:,j)     = - starwall_m_ye(:,j) / &
            ( 1.d0 + tstep * wall_resistivity * starwall_d_yy(:) / wall_thickness )
        end do
        do i = 1, n_dof_starwall
          response_m_pw_rhs(i,:) = starwall_m_ey(i,:) * response_d_ww(:)
        end do
        
        response_m_pp_lhs(:,:) = - starwall_m_ee(:,:) - &
          matmul( starwall_m_ey(:,:), response_m_wp(:,:) )
        response_m_pp_rhs(:,:) = starwall_m_ee(:,:)
        
      !   --- Resistive wall, explicit wall-current time evolution
      else if ( resistive_wall ) then
        
        response_d_ww(:)       = 1.d0 - &
          tstep * wall_resistivity * starwall_d_yy(:) / wall_thickness
        response_m_wp(:,:)     = - starwall_m_ye(:,:)
        response_m_pw_rhs(:,:) = starwall_m_ey(:,:)
        response_m_pp_lhs(:,:) = - starwall_m_ee(:,:)
        response_m_pp_rhs(:,:) = starwall_m_ee(:,:)
        
      !   --- Ideal wall
      else
        
        response_m_pp_lhs(:,:) = - starwall_m_id(:,:)
        response_m_pp_rhs(:,:) = + starwall_m_id(:,:)
        
      end if
      
      ! --- DEBUGGING: Perform some consistency checks and output checksums.
      if ( vacuum_debug ) then
        if ( resistive_wall .and. vacuum_implicit .and. ( wall_resistivity == 0.d0 ) ) then
          write(*,*) LOG_BEGIN_MARKER//' Consistency with ideal wall'
          21 format(1x,a,1x,l,1x,2ES14.5)
          call read_starwall_response(freeboundary_equil, .false.)
          a = sum(abs( response_d_ww - 1.d0 ))
          write(*,21) 'response_d_ww:    ', a < 1.d-5, a
          
          a = sum(abs( response_m_wp + starwall_m_ye ))
          b = sum(abs( starwall_m_ye ))
          write(*,21) 'response_m_wp:    ', a < 1.d-5*b, a/b
          
          a = sum(abs( response_m_pp_rhs - starwall_m_ee ))
          b = sum(abs( starwall_m_ee ))
          write(*,21) 'response_m_pp_rhs:', a < 1.d-5*b, a/b
          
          a = sum(abs( response_m_pp_lhs + starwall_m_id ))
          b = sum(abs( starwall_m_id ))
          write(*,21) 'response_m_pp_lhs:', a < 1.d-5*b, a/b
          
          a = sum(abs( starwall_m_ee-matmul(starwall_m_ey,starwall_m_ye) - starwall_m_id ))
          b = sum(abs( starwall_m_id ))
          write(*,21) 'ee-ey*ye:         ', a < 1.d-5*b, a/b
          write(*,*) LOG_END_MARKER//' Consistency with ideal wall'
        end if
        write(*,*) LOG_BEGIN_MARKER//' Cecksums'
        if ( resistive_wall ) then
          write(*,*) 'd_ww:    ', sum(abs(response_d_ww))
          write(*,*) 'm_wp:    ', sum(abs(response_m_wp))
        end if
        write(*,*) 'm_pp_lhs:', sum(abs(response_m_pp_lhs))
        write(*,*) 'm_pp_rhs:', sum(abs(response_m_pp_rhs))      
        write(*,*) LOG_END_MARKER//' Cecksums'
      end if
      
    end if
    
    write(*,*) LOG_END_MARKER//' update_response'
    
  end subroutine update_response
  
  
  
  
  
  
  !> Read a response matrix from a file.
  subroutine read_response_matrix( matrix, dim_expected, filename )
    
    ! --- Routine parameters
    real*8, allocatable, intent(inout) :: matrix(:,:)     !< Matrix to be read
    integer,             intent(in)    :: dim_expected(2) !< Matrix dimension expected
    character(len=*),    intent(in)    :: filename        !< Filename to read from
    
    ! --- Local variables
    integer :: dim(2), i, j, i2, j2
    
    if ( vacuum_debug ) write(*,*) LOG_BEGIN_MARKER//' read_response_matrix "'//trim(filename)//'"'
    
    open(42, FILE=trim(filename), status='old', action='read')

    read(42,*) dim
    write(*,'(1x,a,2i7)') 'dim=',dim

    if ( (dim(1) /= dim_expected(1)) .or. (dim(2) /= dim_expected(2)) ) then
      write(*,*) 'FATAL ERROR: Matrix dimension not as expected. Different resolutions?'
      stop
    end if
    
    if ( allocated(matrix) ) deallocate(matrix)
    allocate( matrix(dim(1),dim(2)) )
    matrix = 0.d0
    
    do i = 1, dim(1)
      do j = 1, dim(2)
        
        read (42,*) i2, j2, matrix(i,j)
        
        if ( ( i2 /= i ) .or. ( j2 /= j ) ) then
          write(*,*) 'FATAL ERROR: Matrix indices not as expected. Different resolutions?'
          stop
        end if
      
      end do
    end do
      
    close(42)

    if ( vacuum_debug ) write(*,*) LOG_END_MARKER//' read_response_matrix "'//trim(filename)//'"'
    
  end subroutine read_response_matrix
  
  
  
  
  
  
  !> Read a diagonal response matrix from a file.
  subroutine read_response_diagonal( diagonal, dim_expected, filename )
    
    ! --- Routine parameters
    real*8, allocatable, intent(inout) :: diagonal(:)     !< Matrix to be read
    integer,             intent(in)    :: dim_expected    !< Matrix dimension expected
    character(len=*),    intent(in)    :: filename        !< Filename to read from
    
    ! --- Local variables
    integer :: dim, i, i2
    
    if ( vacuum_debug ) write(*,*) LOG_BEGIN_MARKER//' read_response_diagonal "'//trim(filename)//'"'
    
    open(42, FILE=trim(filename), status='old', action='read')

    read(42,*) dim
    write(*,'(1x,a,i7)') 'dim=',dim

    if ( dim /= dim_expected ) then
      write(*,*) 'FATAL ERROR: Matrix dimension not as expected. Different resolutions?'
      stop
    end if
    
    if ( allocated(diagonal) ) deallocate(diagonal)
    allocate( diagonal(dim) )
    diagonal = 0.d0
    
    do i = 1, dim
        
      read (42,*) i2, diagonal(i)
      
      if ( i2 /= i ) then
        write(*,*) 'FATAL ERROR: Matrix indices not as expected. Different resolutions?'
        stop
      end if
        
    end do
      
    close(42)
    
    if ( vacuum_debug ) write(*,*) LOG_END_MARKER//' read_response_diagonal "'//trim(filename)//'"'
    
  end subroutine read_response_diagonal
  
  
  
  
  
  
  !> Determine the response index for a certain boundary degree of freedom.
  integer recursive function response_index(inode, i_starwall, ibas)
    
    ! --- Routine parameters
    integer, intent(in)    :: inode      !< Boundary index of the node
    integer, intent(in)    :: i_starwall !< STARWALL harmonic
    integer, intent(in)    :: ibas       !< Basis function (1 or 2)
    
    if ( NEW_VACUUM ) then
      if ( (i_starwall < 0) .or. (i_starwall > n_starwall_harmonics) ) then
        write(*,*) 'response_index: illegal value i_starwall=', i_starwall
        stop
      end if
      response_index = 2*n_starwall_harmonics*(inode-1) + 2*(i_starwall-1) + ibas
    else
      response_index = 4*(inode-1) + 2*(i_starwall-2) + ibas
    end if
    
    if ( response_index < 1 ) then
      write(*,*) 'FATAL: RESPONSE_INDEX < 1 DETECTED'
      stop
    end if
    
  end function response_index
  
  
  
  
  
  
  !> Character string description of a given toroidal mode. @todo put somewhere else
  character(len=12) function mode_to_str(i_tor, n_period)
  
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: i_tor     !< Toroidal mode index
    integer, intent(in) :: n_period  !< Periodicity
    
    ! --- Local variables
    integer           :: n                ! Toroidal mode number
    character(len=3)  :: typ              ! sin or cos
    character(len=30) :: i_tor_str, n_str ! Character string representations for i_tor and n
    
    write(i_tor_str,'(i)') i_tor

    ! --- Determine toroidal mode number.
    n = int(i_tor / 2) * n_period
    write(n_str,'(i)') n
    
    ! --- Determine mode type (sin or cos).
    if ( mod(i_tor,2) == 0 ) then
      typ = 'sin'
    else
      typ = 'cos'
    end if
    
    mode_to_str = trim(adjustl(i_tor_str))//' (n='//trim(adjustl(n_str))//' '//trim(typ)//')'
    
  end function mode_to_str
  
  
  
  
  
  
  !> Character string description of several toroidal modes. @todo put somewhere else
  character(len=1400) function modes_to_str(i_tors, n_tor, n_period)
  
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: i_tors(n_tor) !< Toroidal mode numbers
    integer, intent(in) :: n_tor         !< Dimension of i_tors
    integer, intent(in) :: n_period      !< Periodicity
    
    ! --- Local variables
    integer :: i
    
    modes_to_str = ''
    do i = 1, n_tor
      modes_to_str = trim(modes_to_str)//' '//mode_to_str(i_tors(i), n_period)
      if ( i == n_tor ) exit
      modes_to_str = trim(modes_to_str)//', '
    end do
    
    modes_to_str = adjustl(modes_to_str)
    
  end function modes_to_str
  
  
  
  
  
  
end module vacuum_response
