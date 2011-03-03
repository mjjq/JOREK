MODULE FOURIER
! Allows to determine the poloidal magnetic coordinate theta_mag from field line tracing in the n=0 component of the
! magnetic field and to perform a Fourier transformation of physical quantities in (theta_mag, phi).
! For details see ../README.jorek2_four
  
  USE parameters,      ONLY: n_vertex_max, n_order, n_plane, n_tor, variable_names
  USE nodes_elements,  ONLY: node_list, element_list
  USE phys_module,     ONLY: F0, xpoint
  
  IMPLICIT NONE
  
  SAVE
  
  PRIVATE
  PUBLIC determine_theta_mag, log_mapping, transform_qtty
  PUBLIC t_theta_mapping
  PUBLIC NEQUIDIST_PTS
  
  TYPE :: t_theta_mapping
    integer              :: nstpts   ! Number of radial positions, theta_mag(theta_geo) is determined
    real*8,  allocatable :: psin(:)  ! Psi_normalized values at the start points
    real*8,  allocatable :: rr(:,:)  ! R coordinates of field line positions at large steps
    real*8,  allocatable :: zz(:,:)  ! Z coordinates of field line positions at large steps
    real*8,  allocatable :: tt(:,:)  ! Geometrical theta at large steps
    real*8,  allocatable :: t2(:,:)  ! Magnetic theta at large steps
    integer, allocatable :: npts(:)  ! Number of large step positions for each start point
    real*8,  allocatable :: rre(:,:) ! R coordinates at equidistant theta_mag(!) positions
    real*8,  allocatable :: zze(:,:) ! Z coordinates at equidistant theta_mag(!) positions
    real*8               :: psi_axis, R_axis, Z_axis, s_axis, t_axis
    integer              :: i_elm_axis
    real*8               :: psi_xpoint, R_xpoint, Z_xpoint, s_xpoint, t_xpoint
    integer              :: i_elm_xpoint
  END TYPE t_theta_mapping
  
  integer, parameter   :: NEQUIDIST_PTS = 32  ! Number of equidistant points for theta_mag(theta_geo) Fourier trafo
  integer, parameter   :: FFTW_ESTIMATE = 64  ! constant of the FFTW library
  real*8,  parameter   :: PI = 4.*ATAN(1.)
  
  
  
  CONTAINS
  
  
  
  ! --- Write information about the mapping determined by determine_theta_mag to the logfile.
  SUBROUTINE log_mapping(mapping)
    
    TYPE(t_theta_mapping), INTENT(IN) :: mapping
    
!    write(*,*) '@@> LOG_MAPPING'
    
    write(*,'(A,I,A)')        'nstpts     =', mapping.nstpts, ' (number of radial positions)'
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
    
!    write(*,*) '@@< LOG_MAPPING'
    
  END SUBROUTINE log_mapping
  
  
  
  ! --- Fourier-transform a physical variable in the magnetic angle theta_mag.
  SUBROUTINE transform_qtty(mapping, ivar, vve, vfour)
    
    TYPE(t_theta_mapping), intent(in)    :: mapping
    integer,               intent(in)    :: ivar         ! Transform which physical variable
    real*8, allocatable,   intent(inout) :: vve(:,:,:)   ! Variable values at (tte,phi_n,psin)
    complex, allocatable,  intent(inout) :: vfour(:,:,:) ! Fourier transformed quantity (m,n,psin)
    
    integer :: i, j, k, iharm, nn
    real*8  :: R_out, Z_out, s_out, t_out
    integer :: i_elm_out, ifail
    real*8  :: v, v_s, v_t
    integer*8 :: fftw_plan
    real*8 :: tmp1(NEQUIDIST_PTS,n_plane-1)
    complex :: tmp2(NEQUIDIST_PTS/2+1,n_plane-1)
    integer :: nequidist_tor
    
    write(*,'("ivar=",I2," (",A,")")') ivar, TRIM(variable_names(ivar))
    
    nequidist_tor = n_plane-1
    
    if ( ALLOCATED(vve)   ) DEALLOCATE(vve)
    if ( ALLOCATED(vfour) ) DEALLOCATE(vfour)
    ALLOCATE( vve(NEQUIDIST_PTS,nequidist_tor,mapping.nstpts), vfour(NEQUIDIST_PTS/2+1,nequidist_tor,mapping.nstpts) )
    vve   = 0.
    vfour = 0.
    
    do k = 1, mapping.nstpts ! radial positions
      do j = 1, nequidist_tor  ! toroidal positions
        do i = 1, NEQUIDIST_PTS  ! poloidal positions
          do iharm = 1, n_tor      ! toroidal harmonics

            CALL find_RZ(node_list,element_list,mapping.rre(k,i-1),mapping.zze(k,i-1),R_out,Z_out,i_elm_out,s_out,t_out,ifail)
            CALL interp0(i_elm_out,ivar,iharm,s_out,t_out,v,v_s,v_t)
            
            nn = iharm / 2 ! toroidal mode number

            if ( iharm == 1 ) then   ! n=0 cos-mode
              vve(i,j,k) = vve(i,j,k) + v
            else if ( MOD(iharm,2) == 0 ) then ! n/=0 cos-modes
              vve(i,j,k) = vve(i,j,k) + v * cos(2.*PI*nn*REAL(j-1)/REAL(nequidist_tor))
            else ! n/=0 sin-modes
              vve(i,j,k) = vve(i,j,k) + v * sin(2.*PI*nn*REAL(j-1)/REAL(nequidist_tor))
            end if
            
          end do
        end do
      end do
      
      call dfftw_plan_dft_r2c_2d(fftw_plan, NEQUIDIST_PTS, nequidist_tor, vve(:,:,k), vfour(:,:,k), FFTW_ESTIMATE)
      call dfftw_execute(fftw_plan)
      call dfftw_destroy_plan(fftw_plan)
      
    end do

    vfour = vfour / REAL( NEQUIDIST_PTS * nequidist_tor )
    
  END SUBROUTINE transform_qtty
  
  
  
  ! --- Determine the 'magnetic poloidal angle' theta_mag by performing field line tracing in the n=0 component of Psi.
  SUBROUTINE determine_theta_mag(nstpts, nmaxsteps, deltaphi, nsmallsteps, mapping)
  
    integer,               intent(in)    :: nstpts      ! Number of radial positions (start points for field line tracing)
    integer,               intent(in)    :: nmaxsteps   ! Maximum number of points per radial position (per traced field line)
    real*8,                intent(in)    :: deltaphi    ! Toroidal step width of the large field line tracing steps
    integer,               intent(in)    :: nsmallsteps ! Number of small steps performed between two small steps
    TYPE(t_theta_mapping), intent(inout) :: mapping
    
    real*8 :: tmp1(0:NEQUIDIST_PTS-1)
    complex :: tmp2(0:NEQUIDIST_PTS/2)

    integer*8            :: fftw_plan
    integer              :: i, j, k
    real*8 :: test1, test2, test3
    
!    write(*,*) '@@> DETERMINE_THETA_MAG'
    
    ! --- Initialisation
    mapping.nstpts = nstpts
    ALLOCATE( mapping.rr(NSTPTS,0:NMAXSTEPS), mapping.zz(NSTPTS,0:NMAXSTEPS), mapping.tt(NSTPTS,0:NMAXSTEPS),    &
              mapping.t2(NSTPTS,0:NMAXSTEPS), mapping.npts(nstpts), mapping.psin(nstpts),                        &
              mapping.rre(NSTPTS,0:NEQUIDIST_PTS-1), mapping.zze(NSTPTS,0:NEQUIDIST_PTS-1) )
    mapping.rr   = 0.
    mapping.zz   = 0.
    mapping.tt   = 0.
    mapping.t2   = 0.
    mapping.npts = 0
    mapping.psin = 0.
    mapping.rre  = 0.
    mapping.zze  = 0.
    
    ! --- Trace field lines in the axisymmetric magnetic field component to be able to determine theta_mag.
    CALL trace_fieldlines(mapping, nmaxsteps, nsmallsteps, deltaphi)
    
    ! --- Interpolate theta_mag to equidistant theta_geo positions.
    do k = 1, mapping.nstpts
      CALL interpolEquidist(mapping.t2(k,0:mapping.npts(k)), mapping.rr(k,0:mapping.npts(k)), mapping.npts(k)+1, &
                            mapping.rre(k,0:NEQUIDIST_PTS-1), mapping.t2(k,0), &
                            mapping.t2(k,0)+2.*PI*REAL(NEQUIDIST_PTS-1)/REAL(NEQUIDIST_PTS), NEQUIDIST_PTS)
      CALL interpolEquidist(mapping.t2(k,0:mapping.npts(k)), mapping.zz(k,0:mapping.npts(k)), mapping.npts(k)+1, &
                            mapping.zze(k,0:NEQUIDIST_PTS-1), mapping.t2(k,0), &
                            mapping.t2(k,0)+2.*PI*REAL(NEQUIDIST_PTS-1)/REAL(NEQUIDIST_PTS), NEQUIDIST_PTS)
    end do
    
    ! --- Output some information for debugging.
    OPEN(42, FILE='determine_theta_mag.rr-zz.dat', ACTION='WRITE', STATUS='REPLACE')
    do k = 1, mapping.nstpts
      do j = 0, mapping.npts(k)
        write(42,*) mapping.rr(k,j), mapping.zz(k,j)
      end do
      write(42,*)
    end do
    CLOSE(42)
    
    OPEN(42, FILE='determine_theta_mag.tt-t2.dat', ACTION='WRITE', STATUS='REPLACE')
    do k = 1, mapping.nstpts
      do j = 0, mapping.npts(k)
        write(42,*) mapping.tt(k,j), mapping.t2(k,j)
      end do
      write(42,*)
    end do
    CLOSE(42)
    
    OPEN(42, FILE='determine_theta_mag.rre-zze.dat', ACTION='WRITE', STATUS='REPLACE')
    do k = 1, mapping.nstpts
      do j = 0, NEQUIDIST_PTS
        write(42,*) mapping.rre(k,MOD(j,NEQUIDIST_PTS)), mapping.zze(k,MOD(j,NEQUIDIST_PTS))
      end do
      write(42,*)
    end do
    CLOSE(42)

    OPEN(42, FILE='determine_theta_mag.rre-zze.2.dat', ACTION='WRITE', STATUS='REPLACE')
    do j = 0, NEQUIDIST_PTS-1
      do k = 1, mapping.nstpts
        write(42,*) mapping.rre(k,j), mapping.zze(k,j)
      end do
      write(42,*)
    end do
    CLOSE(42)

    CALL log_mapping(mapping)
    
!    write(*,*) '@@< DETERMINE_THETA_MAG'
    
  END SUBROUTINE determine_theta_mag
  
  
  
  ! --- Trace field lines in the n=0 component of Psi.
  SUBROUTINE trace_fieldlines(mapping, nmaxsteps, nsmallsteps, deltaphi)
    
    TYPE(t_theta_mapping)   :: mapping
    integer,  intent(in)    :: nmaxsteps
    integer,  intent(in)    :: nsmallsteps
    real*8,   intent(in)    :: deltaphi
  
    integer :: i, j, k, l
    real*8  :: rn, zn ! R and Z at previous small step position
    real*8  :: rp, zp ! R and Z at next small step position
    real*8  :: rh, zh ! R and Z at half small step position
    real*8  :: dpsi_dzn, dpsi_drn ! derivatives of psi at previous small step
    real*8  :: dpsi_dzh, dpsi_drh ! derivatives of psi at half small step
    real*8  :: R_out, Z_out, s_out, t_out
    integer :: i_elm_out, ifail
    real*8 :: P, P_s, P_t, P_st, P_ss, P_tt
    real*8 :: x, x_s, x_t, x_st, x_ss, x_tt, y, y_s, y_t, y_st, y_ss, y_tt
    real*8 :: xjac
    real*8 :: theta_corr
    real   :: smalldeltaphi
    
!    write(*,*) '@@> TRACE_FIELDLINES'
    
    smalldeltaphi = deltaphi / nsmallsteps
    
    ! --- Determine position of axis and xpoint/boundary point; start points will be started between both positions.
    call find_axis(node_list,element_list,mapping.psi_axis,mapping.R_axis,mapping.Z_axis,mapping.i_elm_axis, &
      mapping.s_axis,mapping.t_axis,ifail)
    if (xpoint) then
      call find_xpoint(node_list,element_list,mapping.psi_xpoint,mapping.R_xpoint,mapping.Z_xpoint, &
        mapping.i_elm_xpoint,mapping.s_xpoint,mapping.t_xpoint,ifail)
    else
      ! determine the position of one boundary node as a replacement in a case without an xpoint.
      do i = 1, node_list%n_nodes
        if ( node_list%node(i)%boundary == 2 ) then
          mapping.R_xpoint = node_list%node(i)%x(1,1)
          mapping.Z_xpoint = node_list%node(i)%x(1,2)
          CALL find_RZ(node_list,element_list,mapping.R_xpoint,mapping.Z_xpoint,R_out,Z_out,mapping.i_elm_xpoint, &
            mapping.s_xpoint,mapping.t_xpoint,ifail)
          exit
        end if
      end do
    end if
    
    FL_STPTS: do k = 1, mapping.nstpts
      
      ! Initialize field line position.
      mapping.rr(k,0) = mapping.R_axis + 0.03*( mapping.R_xpoint - mapping.R_axis ) + &
        0.93*( mapping.R_xpoint - mapping.R_axis ) * ( REAL(k-1) / REAL(mapping.nstpts-1) )
      mapping.zz(k,0) = mapping.Z_axis + 0.03*( mapping.Z_xpoint - mapping.Z_axis ) + &
        0.93*( mapping.Z_xpoint - mapping.Z_axis ) * ( REAL(k-1) / REAL(mapping.nstpts-1) )
      mapping.tt(k,0) = atan3( mapping.zz(k,0)-mapping.Z_axis, mapping.rr(k,0)-mapping.R_axis )
      theta_corr = 0.
      
      ! Determine Psi_normalized
      CALL find_RZ(node_list,element_list,mapping.rr(k,0),mapping.zz(k,0),R_out,Z_out,i_elm_out,s_out,t_out,ifail)
      CALL interp0(i_elm_out,1,1,s_out,t_out,P,P_s,P_t)
      mapping.psin(k) = (P-mapping.psi_axis)/(mapping.psi_xpoint-mapping.psi_axis)
      
      FL_LARGESTEPS: do j = 1, NMAXSTEPS
      
        rn = mapping.rr(k,j-1)
        zn = mapping.zz(k,j-1)
        
        FL_SMALLSTEPS: do i = 1, NSMALLSTEPS
        
          ! --- Predictor step
          ! - Determine element number and s and t coordinates for given (R, Z) position.
          CALL find_RZ(node_list,element_list,rn,zn,R_out,Z_out,i_elm_out,s_out,t_out,ifail)
          if ( ifail /= 0 ) then
            write(*,*) '@1@{WARNING}@1@: Error in fourier_modes:trace_fieldlines when calling find_RZ. Reducing mapping.nstpts from ',mapping.nstpts,' to ', k-1, '.'
            mapping.nstpts = k - 1
            EXIT FL_STPTS
          end if
          ! - Interpolate Psi to current position.
          CALL interp0(i_elm_out,1,1,s_out,t_out,P,P_s,P_t)
          ! - Determine derivatives of R and Z with respect to s and t at current position.
          CALL interp_RZ(node_list,element_list,i_elm_out,s_out,t_out,x,x_s,x_t,x_st,x_ss,x_tt,y,y_s,y_t,y_st,y_ss,y_tt)
          ! - Determine derivatives of Psi with respect to R and Z.
          xjac    = x_s*y_t - x_t*y_s
          dpsi_dzn = ( -x_t*P_s + x_s*P_t ) / xjac
          dpsi_drn = ( y_t*P_s - y_s*P_t ) / xjac
          ! - Advance a half step.
          rh    = rn - 1./rn * dpsi_dzn * SMALLDELTAPHI/2. * rn / F0
          zh    = zn + 1./rn * dpsi_drn * SMALLDELTAPHI/2. * rn / F0
          
          ! --- Corrector step
          ! - Determine element number and s and t coordinates for given (R, Z) position.
          CALL find_RZ(node_list,element_list,rh,zh,R_out,Z_out,i_elm_out,s_out,t_out,ifail)
          if ( ifail /= 0 ) then
            write(*,*) '@1@{WARNING}@1@: Error in fourier_modes:trace_fieldlines when calling find_RZ. Reducing mapping.nstpts from ',mapping.nstpts,' to ', k-1, '.'
            mapping.nstpts = k - 1
            EXIT FL_STPTS
          end if
          ! - Interpolate Psi to current position.
          CALL interp0(i_elm_out,1,1,s_out,t_out,P,P_s,P_t)
          ! - Determine derivatives of R and Z with respect to s and t at current position.
          CALL interp_RZ(node_list,element_list,i_elm_out,s_out,t_out,x,x_s,x_t,x_st,x_ss,x_tt,y,y_s,y_t,y_st,y_ss,y_tt)
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
        if ( (mapping.tt(k,j-1)-theta_corr > 3.*PI/2.) .AND. (mapping.tt(k,j) < PI/2.) ) theta_corr = theta_corr + 2.*PI
        if ( (mapping.tt(k,j) > 3.*PI/2.) .AND. (mapping.tt(k,j-1)-theta_corr < PI/2.) ) theta_corr = theta_corr - 2.*PI
        mapping.tt(k,j) = mapping.tt(k,j) + theta_corr
        
        ! --- Finish following field line k, if one poloidal turn has already been completed.
        if ( ABS(mapping.tt(k,j) - mapping.tt(k,0)) > 2*PI ) then
          mapping.npts(k) = j
          do l = 0, mapping.npts(k)
            mapping.t2(k,l) = REAL(l) / ( mapping.npts(k) - &
              ( ( mapping.tt(k,mapping.npts(k)) - mapping.tt(k,0) - 2.*PI ) / &
              ( mapping.tt(k,mapping.npts(k)) - mapping.tt(k,mapping.npts(k)-1) ) ) ) * 2.*PI
          end do
          write(*,'("Field line",I4,":    psin=",F7.3,"    npts=",I7)') k, mapping.psin(k), mapping.npts(k)
          exit
        end if
        
      end do FL_LARGESTEPS
      
    end do FL_STPTS
    
    ! --- Make theta positive.
    do k = 1, mapping.nstpts
      mapping.tt(k,0:mapping.npts(k)) = mapping.tt(k,0:mapping.npts(k)) - MINVAL( mapping.tt(k,0:mapping.npts(k)) )
    end do
    
!    write(*,*) '@@< TRACE_FIELDLINES'
    
  END SUBROUTINE trace_fieldlines
  
  
  
  ! --- Same as interp, but less derivatives are calculated
  RECURSIVE SUBROUTINE interp0(i_elm, i_var, i_harm, s, t, P, P_s, P_t)
  
    real*8 :: G(4,4), G_s(4,4), G_t(4,4), G_st(4,4), G_ss(4,4), G_tt(4,4)
    real*8 :: s, t, P, P_s, P_t
    integer :: kv, iv, kf, i_harm, i_var, i_elm
    
    call basisfunctions2(s,t,G(1:4,1:4),G_s(1:4,1:4),G_t(1:4,1:4),G_st(1:4,1:4),G_ss(1:4,1:4),G_tt(1:4,1:4))
  
    P = 0.d0; P_s = 0.d0; P_t = 0.d0
  
    do kv = 1,n_vertex_max  ! 4 vertices
    
      iv = element_list%element(i_elm)%vertex(kv)  ! the node number
    
      do kf = 1, n_order+1       ! 4 basis functions
    
        P    = P    + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G(kv,kf)
        P_s  = P_s  + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_s(kv,kf)
        P_t  = P_t  + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_t(kv,kf)
    
      end do
    
    end do
    
  END SUBROUTINE interp0
  
  
  
  ! --- Interpolate y(x) data to equidistant x-positions.
  SUBROUTINE interpolEquidist(xOld, yOld, nOld, yNew, xMin, xMax, nNew)
    
    integer, INTENT(IN)    :: nOld         ! Number of old positions
    real*8,  INTENT(IN)    :: xOld(nOld)   ! old x values
    real*8,  INTENT(IN)    :: yOld(nOld)   ! old y values
    integer, INTENT(IN)    :: nNew         ! Number of new positions
    real*8,  INTENT(INOUT) :: yNew(nNew)   ! New y values (OUTPUT)
    real*8,  INTENT(IN)    :: xMin         ! Smallest new x value
    real*8,  INTENT(IN)    :: xMax         ! Largest new x value
    
    integer :: k, j
    real*8  :: xNew
    
    j = 1
    LOOP_OLDPTS: do k = 2, nOld
      LOOP_NEWPTS: do
        xNew = xMin + (xMax - xMin) * REAL(j-1)/REAL(nNew-1)
        if ( xNew <= xOld(k) ) then
          yNew(j) = yOld(k-1) + ( xNew - xOld(k-1) ) / ( xOld(k) - xOld(k-1) ) * ( yOld(k) - yOld(k-1) )
          j = j + 1
        else
          cycle LOOP_OLDPTS
        end if
        if ( j > nNew ) exit LOOP_OLDPTS
      end do LOOP_NEWPTS
    end do LOOP_OLDPTS
    
  END SUBROUTINE interpolEquidist
  
  
  
  ! --- Same as atan2 (Fortran built-in), but returns only positive values.
  real*8 FUNCTION atan3( dy, dx )
    real*8, INTENT(IN) :: dy, dx
    atan3 = ATAN2(dy,dx)
    if ( atan3 < 0 ) atan3 = atan3 + 2.*PI
  END FUNCTION atan3
  
  
  
END MODULE FOURIER
