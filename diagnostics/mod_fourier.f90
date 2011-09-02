!> Allows to determine the poloidal magnetic coordinate theta_mag from field line tracing
!! in the n=0 component of the magnetic field and to perform a Fourier transformation of
!! physical quantities in (theta_mag, phi)
module fourier
  
  use tr_module 
  use parameters,      only: n_vertex_max, n_order, n_plane, n_tor, n_var, variable_names
  use nodes_elements,  only: node_list, element_list
  use phys_module,     only: F0, xpoint
  
  implicit none
  
  save
  
  private
  public determine_theta_mag, log_mapping, transform_qttys
  public t_theta_mapping
  
  type :: t_theta_mapping
    integer              :: nstpts   !< Number of radial positions
    real*8,  allocatable :: psin(:)  !< Psi_normalized values at the start points
    real*8,  allocatable :: rr(:,:)  !< R coordinates of field line positions at large steps
    real*8,  allocatable :: zz(:,:)  !< Z coordinates of field line positions at large steps
    real*8,  allocatable :: tt(:,:)  !< Geometrical theta at large steps
    real*8,  allocatable :: t2(:,:)  !< Magnetic theta at large steps
    integer, allocatable :: npts(:)  !< Number of large step positions for each start point
    real*8,  allocatable :: rre(:,:) !< R coordinates at equidistant theta_mag(!) positions
    real*8,  allocatable :: zze(:,:) !< Z coordinates at equidistant theta_mag(!) positions
    real*8               :: psi_axis, R_axis, Z_axis, s_axis, t_axis
    integer              :: i_elm_axis
    real*8               :: psi_xpoint, R_xpoint, Z_xpoint, s_xpoint, t_xpoint
    integer              :: i_elm_xpoint
  end type t_theta_mapping
  
  integer :: nequidist_pts !< Number of equidistant points for theta_mag(theta_geo) Fourier trafo
  integer, parameter   :: FFTW_ESTIMATE = 64  !< (constant of the FFTW library)
  real*8,  parameter   :: PI =  3.141592653589793_8
  
  
  
  contains
  
  
  
  !> Write information about the mapping determined by determine_theta_mag to the logfile.
  subroutine log_mapping(mapping)
    
    type(t_theta_mapping), intent(in) :: mapping
    
    write(*,'(A,I12,A)')      'nstpts     =', mapping.nstpts, ' (number of radial positions)'
    write(*,'(A,ES12.4,A)')   'R_axis     =', mapping.R_axis, ' (R position of magnetic axis)'
    write(*,'(A,ES12.4,A)')   'Z_axis     =', mapping.Z_axis, ' (Z position of magnetic axis)'
    write(*,'(A,ES12.4,A)')   'psi_axis   =', mapping.psi_axis, ' (poloidal flux at magnetic axis)'
    write(*,'(A,ES12.4,A)')   'R_xpoint   =', mapping.R_xpoint, ' (R position of x-point)'
    write(*,'(A,ES12.4,A)')   'Z_xpoint   =', mapping.Z_xpoint, ' (Z position of x-point)'
    write(*,'(A,ES12.4,A)')   'psi_xpoint =', mapping.psi_xpoint, ' (poloidal flux at x-point)'
    
    write(*,'(A)', ADVANCE='NO') 'psin       ='
    if ( allocated(mapping.psin) ) then
      write(*,'(99F7.3)', ADVANCE='NO') mapping.psin
    else
      write(*,'(A)', ADVANCE='NO') ' not allocated'
    end if
    write(*,'(A)') ' (norm. pol. flux at rad. positions)'
    
    write(*,'(A)', ADVANCE='NO') 'npts       ='
    if ( allocated(mapping.npts) ) then
      write(*,'(99I7)', ADVANCE='NO') mapping.npts
    else
      write(*,'(A)', ADVANCE='NO') ' not allocated'
    end if
    write(*,'(A)') ' (num. of points from field line tracing)'
    
  end subroutine log_mapping
  
  
  
  !> Fourier-transform the physical variables in the magnetic angle theta_mag.
  subroutine transform_qttys(mapping, vfour, m_pol_range)
    
    type(t_theta_mapping), intent(in)    :: mapping
    complex, allocatable,  intent(inout) :: vfour(:,:,:,:) !< Transformed quantities (m,n,irad,ivar)
    integer,               intent(in)    :: m_pol_range(2) !< Range of poloidal mode numbers
    
    real*8, allocatable :: vve(:,:,:,:) ! Variable values (ipol,itor,irad,ivar)
    real*8  :: G(4,4), G_s(4,4), G_t(4,4), G_st(4,4), G_ss(4,4), G_tt(4,4)
    integer :: i, j, k, l, iharm, nn, kv, iv, kf
    real*8  :: R_out, Z_out, s_out, t_out
    integer :: i_elm_out, ifail
    real*8  :: v, basis_function
    integer*8 :: fftw_plan
    integer :: nequidist_tor
    
    nequidist_tor = n_plane-1
    nequidist_pts = 2*(maxval(m_pol_range))
    
    if ( allocated(vfour) ) deallocate(vfour)
    allocate(vve(nequidist_pts,nequidist_tor,mapping.nstpts,n_var))
    allocate(vfour(nequidist_pts/2+1,nequidist_tor,mapping.nstpts,n_var))
    vve   = 0.d0
    vfour = 0.d0
    
    
    do k = 1, mapping.nstpts ! radial positions
      
      write(*,'(1x,a,i4)') 'Transforming variables on surface', k
      
      do j = 1, nequidist_tor  ! toroidal positions
        
        !$omp parallel do                                                                          &
	!$omp   default(shared)                                                                    &
	!$omp   firstprivate(i,iharm,R_out,Z_out,i_elm_out,s_out,t_out,ifail,v,l,nn,kv,iv,kf,      &
	!$omp     basis_function,G,G_s,G_t,G_st,G_ss,G_tt)
        do i = 1, nequidist_pts  ! poloidal positions
          
          call find_RZ(node_list,element_list,mapping.rre(k,i-1),mapping.zze(k,i-1),R_out,Z_out,   &
            i_elm_out,s_out,t_out,ifail)
          
          do iharm = 1, n_tor    ! toroidal harmonics

            nn = iharm / 2 ! toroidal mode number (without periodicity)

            if ( iharm == 1 ) then
	      basis_function = 1.d0
            else if ( MOD(iharm,2) == 0 ) then
	      basis_function = cos(2.*PI*nn*REAL(j-1)/REAL(nequidist_tor))
            else
	      basis_function = sin(2.*PI*nn*REAL(j-1)/REAL(nequidist_tor))
            end if
	    
            call basisfunctions2(s_out,t_out,G(1:4,1:4),G_s(1:4,1:4),G_t(1:4,1:4),G_st(1:4,1:4),   &
              G_ss(1:4,1:4),G_tt(1:4,1:4))
  
	    do l = 1, n_var

              v = 0.d0
  
              do kv = 1, n_vertex_max  ! 4 vertices
    
                iv = element_list%element(i_elm_out)%vertex(kv)  ! the node number
    
                do kf = 1, n_order+1       ! 4 basis functions
    
                  v = v + node_list%node(iv)%values(iharm,kf,l)                                    &
                    * element_list%element(i_elm_out)%size(kv,kf) * G(kv,kf)
    
                end do
    
              end do

              !$omp atomic
              vve(i,j,k,l) = vve(i,j,k,l) + v * basis_function
	    end do
            
          end do
          
        end do
        !$omp end parallel do
        
      end do
      
      do l = 1, n_var
        call dfftw_plan_dft_r2c_2d(fftw_plan, nequidist_pts, nequidist_tor, vve(:,:,k,l),          &
	  vfour(:,:,k,l), FFTW_ESTIMATE)
        call dfftw_execute(fftw_plan)
        call dfftw_destroy_plan(fftw_plan)
      end do
      
    end do

    vfour = vfour / REAL( nequidist_pts * nequidist_tor )
    
    deallocate( vve )
    
  end subroutine transform_qttys
  
  
  
  !> Determine the magnetic poloidal angle, theta_mag, by performing field line tracing
  !! in the axisymmetric magnetic field component.
  subroutine determine_theta_mag(nstpts, nmaxsteps, deltaphi, nsmallsteps, mapping, m_pol_range,   &
    debug)
  
    integer,               intent(in)    :: nstpts      !< Number of radial positions
    integer,               intent(in)    :: nmaxsteps   !< Maximum number of poloidal steps
    real*8,                intent(in)    :: deltaphi    !< Toroidal step width of the large steps
    integer,               intent(in)    :: nsmallsteps !< Number of small between two large steps
    type(t_theta_mapping), intent(inout) :: mapping
    integer,               intent(in)    :: m_pol_range(2) !< Range of poloidal mode numbers
    logical,               intent(in)    :: debug       !< Output debug information?
    
    integer*8            :: fftw_plan
    integer              :: i, j, k
    real*8 :: test1, test2, test3
    
    nequidist_pts = 2*(maxval(m_pol_range))
    
    ! --- Initialisation
    mapping.nstpts = nstpts
    ALLOCATE( mapping.rr(NSTPTS,0:NMAXSTEPS), mapping.zz(NSTPTS,0:NMAXSTEPS),                      &
              mapping.tt(NSTPTS,0:NMAXSTEPS), mapping.t2(NSTPTS,0:NMAXSTEPS),                      &
	      mapping.npts(nstpts), mapping.psin(nstpts), mapping.rre(NSTPTS,0:nequidist_pts-1),   &
	      mapping.zze(NSTPTS,0:nequidist_pts-1) )
    mapping.rr   = 0.
    mapping.zz   = 0.
    mapping.tt   = 0.
    mapping.t2   = 0.
    mapping.npts = 0
    mapping.psin = 0.
    mapping.rre  = 0.
    mapping.zze  = 0.
    
    ! --- Trace field lines in the axisymmetric magnetic field component
    call trace_fieldlines(mapping, nmaxsteps, nsmallsteps, deltaphi)
    
    ! --- Interpolate theta_mag to equidistant theta_geo positions.
    do k = 1, mapping.nstpts
      if ( mapping.npts(k) < nequidist_pts ) then
        write(*,'(1x,a)')       'WARNING: Low number of points from field line tracing.'
        write(*,'(3x,a)')       'You may need to decrease deltaphi in "four_params.nml"'
        write(*,'(3x,a,f5.3)')  'by a factor > ', real(nequidist_pts)/real(mapping.npts(k))
        write(*,'(3x,a,i3)')    'Field line', k
      end if
      call interpolEquidist(mapping.t2(k,0:mapping.npts(k)), mapping.rr(k,0:mapping.npts(k)),      &
        mapping.npts(k)+1, mapping.rre(k,0:nequidist_pts-1), mapping.t2(k,0),                      &
        mapping.t2(k,0)+2.*PI*REAL(nequidist_pts-1)/REAL(nequidist_pts), nequidist_pts)
      call interpolEquidist(mapping.t2(k,0:mapping.npts(k)), mapping.zz(k,0:mapping.npts(k)),      &
        mapping.npts(k)+1, mapping.zze(k,0:nequidist_pts-1), mapping.t2(k,0),                      &
        mapping.t2(k,0)+2.*PI*REAL(nequidist_pts-1)/REAL(nequidist_pts), nequidist_pts)
    end do
    
    ! --- Output some information for debugging.
    if ( debug ) then
      open(42, FILE='determine_theta_mag.rr-zz.dat', ACTION='WRITE', STATUS='REPLACE')
      do k = 1, mapping.nstpts
        do j = 0, mapping.npts(k)
          write(42,*) mapping.rr(k,j), mapping.zz(k,j)
        end do
        write(42,*)
      end do
      close(42)
      
      open(42, FILE='determine_theta_mag.tt-t2.dat', ACTION='WRITE', STATUS='REPLACE')
      do k = 1, mapping.nstpts
        do j = 0, mapping.npts(k)
          write(42,*) mapping.tt(k,j), mapping.t2(k,j)
        end do
        write(42,*)
      end do
      close(42)
      
      open(42, FILE='determine_theta_mag.rre-zze.dat', ACTION='WRITE', STATUS='REPLACE')
      do k = 1, mapping.nstpts
        do j = 0, nequidist_pts
          write(42,*) mapping.rre(k,MOD(j,nequidist_pts)), mapping.zze(k,MOD(j,nequidist_pts))
        end do
        write(42,*)
      end do
      close(42)
  
      open(42, FILE='determine_theta_mag.rre-zze.2.dat', ACTION='WRITE', STATUS='REPLACE')
      do j = 0, nequidist_pts-1
        do k = 1, mapping.nstpts
          write(42,*) mapping.rre(k,j), mapping.zze(k,j)
        end do
        write(42,*)
      end do
      close(42)

      call log_mapping(mapping)
    end if
    
  end subroutine determine_theta_mag
  
  
  
  !> Trace field lines in the n=0 component of Psi.
  subroutine trace_fieldlines(mapping, nmaxsteps, nsmallsteps, deltaphi)
    
    type(t_theta_mapping)   :: mapping
    integer,  intent(in)    :: nmaxsteps    !< Maximum number of poloidal steps
    integer,  intent(in)    :: nsmallsteps  !< Number of small between two large steps
    real*8,   intent(in)    :: deltaphi     !< Toroidal step width of the large steps
  
    integer :: i, j, k, l
    real*8  :: rn, zn ! R and Z at previous small step position
    real*8  :: rp, zp ! R and Z at next small step position
    real*8  :: rh, zh ! R and Z at half small step position
    real*8  :: dpsi_dzn, dpsi_drn ! derivatives of psi at previous small step
    real*8  :: dpsi_dzh, dpsi_drh ! derivatives of psi at half small step
    real*8  :: R_out, Z_out, s_out, t_out
    integer :: i_elm_out, ifail,my_id
    real*8 :: P, P_s, P_t, P_st, P_ss, P_tt
    real*8 :: x, x_s, x_t, y, y_s, y_t
    real*8 :: xjac
    real*8 :: theta_corr
    real   :: smalldeltaphi
    logical :: do_not_continue
    
    smalldeltaphi = deltaphi / nsmallsteps
    my_id = 0
    
    ! --- Determine position of axis and xpoint/boundary point
    !     Field line tracing will start on nstpts positions between axis and x-point
    call find_axis(my_id,node_list,element_list,mapping.psi_axis,mapping.R_axis,mapping.Z_axis,    &
      mapping.i_elm_axis, mapping.s_axis,mapping.t_axis,ifail)
    if (xpoint) then
      call find_xpoint(my_id,node_list,element_list,mapping.psi_xpoint,mapping.R_xpoint,           &
        mapping.Z_xpoint, mapping.i_elm_xpoint,mapping.s_xpoint,mapping.t_xpoint,ifail)
    else
      ! Determine the position of a boundary node for non-xpoint cases.
      do i = 1, node_list%n_nodes
        if ( node_list%node(i)%boundary == 2 ) then
          mapping.R_xpoint = node_list%node(i)%x(1,1)
          mapping.Z_xpoint = node_list%node(i)%x(1,2)
          call find_RZ(node_list,element_list,mapping.R_xpoint,mapping.Z_xpoint,R_out,Z_out,       &
	    mapping.i_elm_xpoint,mapping.s_xpoint,mapping.t_xpoint,ifail)
          exit
        end if
      end do
    end if
    
    do_not_continue = .false.
    
    !$omp parallel do                                                                              &
    !$omp   schedule(dynamic)                                                                      &
    !$omp   default(shared)                                                                        &
    !$omp   firstprivate(k,theta_corr,R_out,Z_out,i_elm_out,s_out,t_out,ifail,P,P_s,P_t,j,rn,zn,i, &
    !$omp     x,x_s,x_t,y,y_s,y_t,xjac,dpsi_dzn,dpsi_drn,rh,zh,rp,zp)
    FL_STPTS: do k = mapping.nstpts, 1, -1
      if ( do_not_continue ) cycle
      
      ! --- Initialize field line position.
      mapping.rr(k,0) = mapping.R_axis + 0.03*( mapping.R_xpoint - mapping.R_axis ) + &
        0.93*( mapping.R_xpoint - mapping.R_axis ) * ( REAL(k-1) / REAL(mapping.nstpts-1) )
      mapping.zz(k,0) = mapping.Z_axis + 0.03*( mapping.Z_xpoint - mapping.Z_axis ) + &
        0.93*( mapping.Z_xpoint - mapping.Z_axis ) * ( REAL(k-1) / REAL(mapping.nstpts-1) )
      mapping.tt(k,0) = atan3( mapping.zz(k,0)-mapping.Z_axis, mapping.rr(k,0)-mapping.R_axis )
      theta_corr = 0.
      
      ! --- Determine Psi_normalized
      call find_RZ(node_list,element_list,mapping.rr(k,0),mapping.zz(k,0),R_out,Z_out,i_elm_out,   &
        s_out,t_out,ifail)
      call interp0(i_elm_out,1,1,s_out,t_out,P,P_s,P_t)
      mapping.psin(k) = (P-mapping.psi_axis)/(mapping.psi_xpoint-mapping.psi_axis)
      
      FL_LARGESTEPS: do j = 1, NMAXSTEPS
        if ( do_not_continue ) cycle
      
        rn = mapping.rr(k,j-1)
        zn = mapping.zz(k,j-1)
        
        FL_SMALLSTEPS: do i = 1, NSMALLSTEPS
          if ( do_not_continue ) cycle
        
          ! --- Predictor step
          ! - Determine element number and s and t coordinates for given (R, Z) position.
          call find_RZ(node_list,element_list,rn,zn,R_out,Z_out,i_elm_out,s_out,t_out,ifail)
          if ( ifail /= 0 ) then
            write(*,*) 'WARNING: Error in fourier:trace_fieldlines calling find_RZ.'
            do_not_continue = .true.
	    cycle
          end if
          ! - Interpolate Psi to current position.
          call interp0(i_elm_out,1,1,s_out,t_out,P,P_s,P_t)
          ! - Determine derivatives of R and Z with respect to s and t at current position.
          call interp_RZ0(node_list,element_list,i_elm_out,s_out,t_out,x,x_s,x_t,y,y_s,y_t)
          ! - Determine derivatives of Psi with respect to R and Z.
          xjac    = x_s*y_t - x_t*y_s
          dpsi_dzn = ( -x_t*P_s + x_s*P_t ) / xjac
          dpsi_drn = ( y_t*P_s - y_s*P_t ) / xjac
          ! - Advance a half step.
          rh    = rn - 1./rn * dpsi_dzn * SMALLDELTAPHI/2. * rn / F0
          zh    = zn + 1./rn * dpsi_drn * SMALLDELTAPHI/2. * rn / F0
          
          ! --- Corrector step
          ! - Determine element number and s and t coordinates for given (R, Z) position.
          call find_RZ(node_list,element_list,rh,zh,R_out,Z_out,i_elm_out,s_out,t_out,ifail)
          if ( ifail /= 0 ) then
            write(*,*) 'WARNING: Error in fourier:trace_fieldlines calling find_RZ.'
            do_not_continue = .true.
	    cycle
          end if
          ! - Interpolate Psi to current position.
          call interp0(i_elm_out,1,1,s_out,t_out,P,P_s,P_t)
          ! - Determine derivatives of R and Z with respect to s and t at current position.
          call interp_RZ0(node_list,element_list,i_elm_out,s_out,t_out,x,x_s,x_t,y,y_s,y_t)
          ! - Determine derivatives of Psi with respect to R and Z.
          xjac    = x_s*y_t - x_t*y_s
          dpsi_dzh = ( -x_t*P_s + x_s*P_t ) / xjac
          dpsi_drh = ( y_t*P_s - y_s*P_t ) / xjac
          ! - Advance a full step.
          rp    = rn - 1./rh * dpsi_dzh * SMALLDELTAPHI * rh / F0
          zp    = zn + 1./rh * dpsi_drh * SMALLDELTAPHI * rh / F0

        end do FL_SMALLSTEPS
        
        ! --- Write result of the small steps into the position arrays.
        mapping.rr(k,j) = rp
        mapping.zz(k,j) = zp
        mapping.tt(k,j) = atan3( zp-mapping.Z_axis, rp-mapping.R_axis )
        
        ! --- Compensate for theta jumps by 2*PI.
        if ( (mapping.tt(k,j-1)-theta_corr > 3.*PI/2.) .AND. (mapping.tt(k,j) < PI/2.) )           &
	  theta_corr = theta_corr + 2.*PI
        if ( (mapping.tt(k,j) > 3.*PI/2.) .AND. (mapping.tt(k,j-1)-theta_corr < PI/2.) )           &
	  theta_corr = theta_corr - 2.*PI
        mapping.tt(k,j) = mapping.tt(k,j) + theta_corr
        
        ! --- Finish following field line k, if one poloidal turn has already been completed.
        if ( ABS(mapping.tt(k,j) - mapping.tt(k,0)) > 2*PI ) then
          mapping.npts(k) = j
          do l = 0, mapping.npts(k)
            mapping.t2(k,l) = REAL(l) / ( mapping.npts(k) - &
              ( ( mapping.tt(k,mapping.npts(k)) - mapping.tt(k,0) - 2.*PI ) / &
              ( mapping.tt(k,mapping.npts(k)) - mapping.tt(k,mapping.npts(k)-1) ) ) ) * 2.*PI
          end do
          !$omp critical
          write(*,'(" Field line",I4,":    psin=",F7.3,"    npts=",I7)') k, mapping.psin(k),       &
	    mapping.npts(k)
          !$omp end critical
          exit
        end if
        
      end do FL_LARGESTEPS
      
    end do FL_STPTS
    !$omp end parallel do
    
    ! --- Make theta positive.
    do k = 1, mapping.nstpts
      mapping.tt(k,0:mapping.npts(k)) = mapping.tt(k,0:mapping.npts(k))                            &
        - MINVAL( mapping.tt(k,0:mapping.npts(k)) )
    end do
    
  end subroutine trace_fieldlines
  
  
  
  !> Same as interp, but less derivatives are calculated
  recursive subroutine interp0(i_elm, i_var, i_harm, s, t, P, P_s, P_t)
  
    real*8 :: G(4,4), G_s(4,4), G_t(4,4), G_st(4,4), G_ss(4,4), G_tt(4,4)
    real*8 :: s, t, P, P_s, P_t
    integer :: kv, iv, kf, i_harm, i_var, i_elm
    
    call basisfunctions2(s,t,G(1:4,1:4),G_s(1:4,1:4),G_t(1:4,1:4),G_st(1:4,1:4),G_ss(1:4,1:4),     &
      G_tt(1:4,1:4))
  
    P = 0.d0; P_s = 0.d0; P_t = 0.d0
  
    do kv = 1,n_vertex_max  ! 4 vertices
    
      iv = element_list%element(i_elm)%vertex(kv)  ! the node number
    
      do kf = 1, n_order+1       ! 4 basis functions
    
        P   = P   + node_list%node(iv)%values(i_harm,kf,i_var)                                     &
	  * element_list%element(i_elm)%size(kv,kf) * G(kv,kf)
        P_s = P_s + node_list%node(iv)%values(i_harm,kf,i_var)                                     &
	  * element_list%element(i_elm)%size(kv,kf) * G_s(kv,kf)
        P_t = P_t + node_list%node(iv)%values(i_harm,kf,i_var)                                     &
	  * element_list%element(i_elm)%size(kv,kf) * G_t(kv,kf)
    
      end do
    
    end do
    
  end subroutine interp0
  
  
  
  !> Same as routine interp_RZ but with less derivatives
  recursive subroutine interp_RZ0(node_list,element_list,i_elm,s,t,R,R_s,R_t,Z,Z_s,Z_t)

  use data_structure
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),    intent(in)  :: node_list
  type (type_element_list), intent(in)  :: element_list
  real*8,                   intent(in)  :: s,t
  real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
  
  ! --- Local variables
  real*8  :: G(4,4), G_s(4,4), G_t(4,4), G_st(4,4), G_ss(4,4), G_tt(4,4)
  real*8  :: xx1, xx2, ss
  integer :: kv, iv, kf, i_elm
  
  call basisfunctions2(s,t,G,G_s,G_t,G_st,G_ss,G_tt)
  
  R = 0.d0; R_s = 0.d0; R_t = 0.d0
  Z = 0.d0; Z_s = 0.d0; Z_t = 0.d0
  
  do kv = 1,n_vertex_max  ! 4 vertices
  
    iv = element_list%element(i_elm)%vertex(kv)  ! the node number
  
    do kf = 1, n_order+1       ! 4 basis functions
      
      xx1 = node_list%node(iv)%x(kf,1)
      xx2 = node_list%node(iv)%x(kf,2)
      ss  = element_list%element(i_elm)%size(kv,kf)
      
      R    = R    + xx1 * ss * G(kv,kf)
      R_s  = R_s  + xx1 * ss * G_s(kv,kf)
      R_t  = R_t  + xx1 * ss * G_t(kv,kf)
  
      Z    = Z    + xx2 * ss * G(kv,kf)
      Z_s  = Z_s  + xx2 * ss * G_s(kv,kf)
      Z_t  = Z_t  + xx2 * ss * G_t(kv,kf)
  
    end do
  
  end do
  
  return
  end subroutine interp_RZ0
  
  
  
  !> Interpolate y(x) data to equidistant x-positions.
  subroutine interpolEquidist(xOld, yOld, nOld, yNew, xMin, xMax, nNew)
    
    integer, intent(in)    :: nOld         !< Number of old positions
    real*8,  intent(in)    :: xOld(nOld)   !< old x values
    real*8,  intent(in)    :: yOld(nOld)   !< old y values
    integer, intent(in)    :: nNew         !< Number of new positions
    real*8,  intent(inout) :: yNew(nNew)   !< New y values
    real*8,  intent(in)    :: xMin         !< Smallest new x value
    real*8,  intent(in)    :: xMax         !< Largest new x value
    
    integer :: k, j
    real*8  :: xNew
    
    j = 1
    LOOP_OLDPTS: do k = 2, nOld
      LOOP_NEWPTS: do
        xNew = xMin + (xMax - xMin) * REAL(j-1)/REAL(nNew-1)
        if ( xNew <= xOld(k) ) then
          yNew(j) = yOld(k-1) + ( xNew - xOld(k-1) ) / ( xOld(k) - xOld(k-1) )                     &
	    * ( yOld(k) - yOld(k-1) )
          j = j + 1
        else
          cycle LOOP_OLDPTS
        end if
        if ( j > nNew ) exit LOOP_OLDPTS
      end do LOOP_NEWPTS
    end do LOOP_OLDPTS
    
  end subroutine interpolEquidist
  
  
  
  !> Same as atan2 (Fortran built-in), but returns only positive values by adding 2*pi
  real*8 function atan3( dy, dx )
    real*8, intent(in) :: dy, dx
    atan3 = ATAN2(dy,dx)
    if ( atan3 < 0 ) atan3 = atan3 + 2.*PI
  end function atan3
  
  
  
end module fourier
