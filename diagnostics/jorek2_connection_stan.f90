program jorek2_poincare
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  use data_structure
  use phys_module
  use basis_at_gaussian
  use elements_nodes_neighbours
  use constants
  use mod_import_restart
  use mod_neighbours
  use mod_interp
  use equil_info
  use mod_element_rtree
  use mpi
  use equil_info
  
  implicit none
  
  ! --- Poincare data
  real*8,allocatable    :: rp(:), zp(:), R_all(:), Z_all(:), C_all(:), C_minus(:), C_plus(:)
  real*4,allocatable    :: C_all_4(:)
  real*4,allocatable    :: R_strike(:),  Z_strike(:), P_strike(:)        ! position of strike points
  real*8,allocatable    :: C_strike(:),  B_strike(:)                   ! connection length, boundary type at strike points
  real*8,allocatable    :: T0_strike(:), T_strike(:)                   ! temperature at start and end of fieldline
  real*8,allocatable    :: ZN0_strike(:), ZN_strike(:)                 ! density at start and end of fieldline
  real*8,allocatable    :: PS0_strike(:)                               ! flux at starting point
  real*8,allocatable    :: R_turn(:,:), Z_turn(:,:), C_turn(:,:), C_turn_tmp(:,:)
  real*8,allocatable    :: T_turn(:,:), PSI_turn(:,:), ZN_turn(:,:)
  
  ! --- Extra data
  integer		:: my_id, ikeep, n_cpu, ierr, nsend, nrecv, ikeep0, inode1, inode2, i_line0
  integer		:: elm_start, elm_end, elm_delta, local_elm_start, local_elm_end
  integer		:: i_elm_half, cpu_fraction
  integer		:: nnos, n_scalars, ivtk, i_var, i_strike, i_strike0
  integer		:: i, j, iside_i, iside_j, ip, i_line, n_lines, i_tor, i_harm, i_var_psi, i_dir, k, m, ns, nt
  integer		:: i_elm, ifail, i_phi, n_phi, i_turn, n_turns, i_elm_out, i_elm_prev, i_elm_tmp,i_steps, n_turn_max(2)
  real*8		:: R_start, Z_start, P_start
  real*8		:: R_line, Z_line, s_line, t_line, p_line
  real*8		:: Zjac
  real*8		:: s_mid, t_mid, p_mid
  real*8		:: s_ini, t_ini
  real*8		:: s_out, t_out
  real*8		:: R, R_s, R_t, R_st, R_ss, R_tt, R_in, R_out, R_keep, Rmid, Rmid_s, Rmid_t
  real*8		:: Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, Z_in, Z_out, Z_keep, Zmid, Zmid_s, Zmid_t
  real*8		:: P, P_s, P_t, P_st, P_ss, P_tt
  real*8		:: tol, psi_s, psi_t
  real*8		:: Rmin, Rmax
  real*8		:: Zmin, Zmax
  real*8		:: delta_phi, delta_phi_local, delta_phi_step, total_phi
  real*8		:: delta_s, small_delta_s
  real*8		:: delta_t, small_delta_t
  real*8		:: small_delta, dl2, total_length, length_max
  real*8		:: zl1, zl2, partial(2)
  real*8		:: psi_bnd, psi_bnd2
  integer		:: bnd_tmp, bnd_tmp_opp
  real*8		:: s_tmp,   s_tmp_opp
  real*8		:: t_tmp,   t_tmp_opp
  real*8		:: value_out
  real*4,allocatable	:: RZkeep(:,:),scalars(:,:)
  real*4		:: ZERO
  integer		:: status(MPI_STATUS_SIZE)
  character		:: buffer*80, lf*1, str1*12, str2*12
  character*12, allocatable :: scalar_names(:)
  integer               :: count_lines, count_lines_tot
  real*8                :: C_average,   C_average_tot
  real*8                :: small_r, theta_pol
  logical, parameter    :: Rtheta_plot = .false.
  integer               :: n_points_max
  real*8                :: R_tmp,Z_tmp
  real*8                :: psi_tmp,P0_s,P0_t,P0_st,P0_ss,P0_tt
  integer               :: ielm_tmp
  integer               :: n_points_start
  real*8                :: Rstart_min, Rstart_max
  real*8                :: phi_start, progress
  real*8                :: xjac, dpsi_dR, dpsi_dZ, phi, R_half, Z_half
  real*8                :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt
  real*8                :: Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
  real*8                :: Fprof,Fprof_s,Fprof_t,Fprof_st,Fprof_ss,Fprof_tt
  real*8                :: P0
  real*8                :: AR0,AR0_s,AR0_t,AR0_st,AR0_ss,AR0_tt,AR0_R,AR0_Z,AR0_p
  real*8                :: AZ0,AZ0_s,AZ0_t,AZ0_st,AZ0_ss,AZ0_tt,AZ0_R,AZ0_Z,AZ0_p
  real*8                :: A30,A30_s,A30_t,A30_st,A30_ss,A30_tt,A30_R,A30_Z,A30_p
  real*8                :: BR, BZ, Bp
  real*8                :: HHZ(n_tor), HHZ_p(n_tor)
  
  write(*,*) '***************************************'
  write(*,*) '* JOREK2_poincare                     *'
  write(*,*) '***************************************'
  
  
  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  ! --- First part: Initialisation
  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  
  ! --- MPI initilisation
  call MPI_INIT(IERR)
  !required=MPI_THREAD_MULTIPLE
  !call MPI_Init_thread(required,provided,StatInfo)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
  
  ! --- Initilise data
  call initialise_parameters(my_id, "__NO_FILENAME__")
  do i_tor=1, n_tor
    mode(i_tor) = + int(i_tor / 2) * n_period
  enddo
  call broadcast_elements(my_id, element_list)                ! elements
  call broadcast_nodes(my_id, node_list)                      ! nodes
  call populate_element_rtree(node_list, element_list)
  call broadcast_phys(my_id)                                  ! physics parameters
  call broadcast_equil_state(my_id)                           ! equil_state
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)
  call initialise_basis                                       ! define the basis functions at the Gaussian points
  
  ! --- Broadcast accross MPIs
  call broadcast_elements(my_id, element_list)                ! elements
  call broadcast_nodes(my_id, node_list)                      ! nodes
  call broadcast_phys(my_id)                                  ! physics parameters
  
  ! --- Define element neighbours
  allocate(element_neighbours(4,element_list%n_elements))
  element_neighbours = 0
  do i=1,element_list%n_elements
    do j=i+1,element_list%n_elements
      if (neighbours(node_list,element_list%element(i),element_list%element(j),iside_i,iside_j)) then
        element_neighbours(iside_i,i) = j
        element_neighbours(iside_j,j) = i
      endif
    enddo
  enddo
  
  ! --- Initialise variables
  n_scalars = 4
  ZERO      = 0.
  i_var_psi = 1                                 ! the index of the magnetic flux variable
  allocate(scalar_names(n_scalars))
  scalar_names  = (/ 'length_tot  ','length_min  ','psi_start   ','T_start     '/)
  
  
  ! --- Initialise toroidal variables
  n_points_start  = 20 ! element_list%n_elements
  n_turns   = 500 !500                         ! number of toroidal turns to follow a fieldline
  n_phi     = 1000 !1000                        ! number of steps per toroidal turn
  phi_start = 2.0*PI/12.0
    
  ! --- How many starting points per element
  ns      = 1!3                                 ! number of (s) starting points within one element
  nt      = 1!3                                 ! number of (t) starting points within one element
  if (n_points_start .ne. element_list%n_elements) ns = 1
  if (n_points_start .ne. element_list%n_elements) nt = 1
  n_lines = n_points_start * ns * nt    ! number of starting points
  
  ! --- Allocate data
  allocate(R_strike(n_lines),Z_strike(n_lines),P_strike(n_lines),C_strike(n_lines),B_strike(n_lines))
  allocate(T0_strike(n_lines),T_strike(n_lines),ZN0_strike(n_lines),ZN_strike(n_lines),PS0_strike(n_lines))
  allocate(R_all(n_lines),Z_all(n_lines),C_all(n_lines),C_minus(n_lines),C_plus(n_lines))
  allocate(C_all_4(n_lines))
  allocate(R_turn(n_turns+1,2),Z_turn(n_turns+1,2),C_turn(n_turns+1,2),C_turn_tmp(n_turns+1,2))
  allocate(T_turn(n_turns+1,2),PSI_turn(n_turns+1,2),ZN_turn(n_turns+1,2))
  n_points_max = 10000000
  allocate(RZkeep(2,n_points_max),scalars(n_points_max,n_scalars))

  ! --- Initialise allocated data
  R_all     = 0.d0; Z_all     = 0.d0; C_all     = 0.d0;  C_all_4    = 0.d0
  C_minus   = 0.d0; C_plus    = 0.d0
  R_strike  = 0.d0; Z_strike  = 0.d0; P_strike  = 0.d0;  C_strike   = 0.d0
  T0_strike = 0.d0; T_strike  = 0.d0; ZN0_strike = 0.d0; ZN_strike  = 0.d0; PS0_strike = 0.d0
  R_turn    = 0.d0; Z_turn    = 0.d0; C_turn    = 0.d0;  C_turn_tmp = 0.d0
  
  ! --- Get domain limits
  Rmin = 1.d20; Rmax = -1.d20; Zmin = 1.d20; Zmax=-1.d20
  do i=1,node_list%n_nodes
    Rmin = min(Rmin,node_list%node(i)%x(1,1,1))
    Rmax = max(Rmax,node_list%node(i)%x(1,1,1))
    Zmin = min(Zmin,node_list%node(i)%x(1,1,2))
    Zmax = max(Zmax,node_list%node(i)%x(1,1,2))
  enddo
  
  ! --- find x-point(s)
  if (xpoint) then
    psi_bnd  = ES%psi_xpoint(1)
    psi_bnd2 = ES%psi_xpoint(2)
    if( ES%active_xpoint .eq. UPPER_XPOINT ) then
      psi_bnd  = ES%psi_xpoint(2)
      psi_bnd2 = ES%psi_xpoint(1)
    endif
    if (xcase .eq. LOWER_XPOINT) psi_bnd2 = psi_bnd
  else
    psi_bnd  = 0.d0
    psi_bnd2 = 0.d0
  endif
   
  if (n_points_start .ne. element_list%n_elements) then
  
    ! --- The elements we are looking at (note element indexing starts at axis going gradually outwards)
    elm_start = 1
    elm_end   = n_points_start
    
    Rstart_min = R_axis + 0.05 !2.5
    Rstart_max = Rstart_min + 1.d0 !  2.1
    
    ! --- The elements our local MPI is looking at
    elm_delta = (elm_end - elm_start) / n_cpu + 1
    local_elm_start = elm_start + my_id*elm_delta
    local_elm_end   = min(elm_end,elm_start+(my_id+1)*elm_delta - 1)
  
  else
    ! --- The elements we are looking at (note element indexing starts at axis going gradually outwards)
    elm_start = 1!element_list%n_elements/4
    elm_end   = element_list%n_elements
    if (Rtheta_plot) elm_start = 1
    
    ! --- The elements our local MPI is looking at
    ! --- We divide the elements in two halves, and put 5/6 of the MPIs on the first half
    ! --- Because elements in the core take much longer to run...
    i_elm_half    = (elm_end - elm_start) / 2
    cpu_fraction = real(5)/real(6) * n_cpu
    if (my_id .lt. cpu_fraction) then
      elm_delta       = (i_elm_half - elm_start) / cpu_fraction
      local_elm_start = elm_start + my_id*elm_delta + 1
      local_elm_end   = min(i_elm_half,elm_start+(my_id+1)*elm_delta)
    else
      elm_delta       = (elm_end - i_elm_half+1) / (n_cpu-cpu_fraction)
      local_elm_start = i_elm_half + (my_id-cpu_fraction)*elm_delta + 1
      local_elm_end   = min(elm_end,i_elm_half+((my_id-cpu_fraction)+1)*elm_delta)
    endif
    !elm_delta = (elm_end - elm_start) / n_cpu
    !local_elm_start = elm_start + my_id*elm_delta + 1
    !local_elm_end   = min(elm_end,elm_start+(my_id+1)*elm_delta)
    
    ! --- Some info print outs
    write(*,*) ' Total number of elements : ',element_list%n_elements
    write(*,*) ' Elements considered (start,end) : ',elm_start, elm_end
    write(*,*) ' Local MPI elements (mpi_id, elm_start, elm_end) : ', my_id, local_elm_start, local_elm_end
  
  endif
  
  ! --- Initialise counters
  i_line      = 0
  i_strike    = 0
  ikeep       = 0
  count_lines = 0
  C_average   = 0.d0
  
  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  ! --- Second part: Loop over starting points and get poincare points
  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  
  ! --- Start going over all local elements
  do i = local_elm_start, local_elm_end
    
    ! --- Print progress
    if (local_elm_end .ne. local_elm_start) then
      progress = real(i-local_elm_start)/real(local_elm_end-local_elm_start)*100.d0
    else
      progress = 100.0
    endif
    if (mod(i,10).eq.0) &
      write(*,'(A,5i6,A)')' Progress: MPI_id, element, start, end :', my_id, i, local_elm_start, local_elm_end,int(progress),'%'
    if (n_points_start .ne. element_list%n_elements) &
      write(*,'(A,5i6,A)')' Progress: MPI_id, element, start, end :', my_id, i, local_elm_start, local_elm_end,int(progress),'%'
    
    i_line = i_line + 1
    
    ! --- Initialise temporary data (for the line)
    R_turn     = 0.d0
    Z_turn     = 0.d0
    C_turn_tmp = 0.d0
    
    ! --- Do both directions of field line
    do i_dir = -1,1,2

      ! --- The toroidal step (with direction)
      delta_phi = 2.d0 * PI * float(i_dir) / float(n_period*n_phi)
      
    
      R_start = Rstart_min + float(i-1)/float(n_points_start-1) * (Rstart_max - Rstart_min)
      Z_start = Z_axis
      P_start = phi_start
      R_turn(1,(i_dir+1)/2+1) = R_start
      Z_turn(1,(i_dir+1)/2+1) = Z_start
      C_turn(1,(i_dir+1)/2+1) = 0.d0
      R_line = R_start
      Z_line = Z_start
      phi = phi_start
 
      i_elm = i   ! added by me
      call var_value(i_elm, var_T,   s_line,t_line,P_start, T_turn  (1,(i_dir+1)/2+1) )
      call var_value(i_elm, var_A3,  s_line,t_line,P_start, PSI_turn(1,(i_dir+1)/2+1) )
      call var_value(i_elm, var_rho, s_line,t_line,P_start, ZN_turn (1,(i_dir+1)/2+1) )
     
      ! --- Loop over toroidal turns
      do i_turn = 1, n_turns
        
        ! --- Record the maximum number of turns
        n_turn_max((i_dir+1)/2+1) = i_turn
      
        ! --- Perform the field line tracing.
        do j = 1, n_phi
        
          ! --- toroidal functions
          HHZ  (1) = 1.d0
          HHZ_p(1) = 0.d0
          do i_tor=1,(n_tor-1)/2
            HHZ  (2*i_tor)   = +cos(mode(2*i_tor)   * phi )
            HHZ  (2*i_tor+1) = +sin(mode(2*i_tor+1) * phi )
            HHZ_p(2*i_tor)   = -sin(mode(2*i_tor)   * phi ) * float(mode(2*i_tor))
            HHZ_p(2*i_tor+1) = +cos(mode(2*i_tor+1) * phi ) * float(mode(2*i_tor+1))
          enddo
          
          ! -----------------
          ! --- half step ---
          ! -----------------
          call find_RZ(node_list,element_list,R_line,Z_line,R_out,Z_out,i_elm,s_out,t_out,ifail)
          if ( ifail .eq. 0 ) then
            ! --- RZ-coords
            call interp_RZ(node_list,element_list,i_elm,s_out,t_out,R,R_s,R_t,Z,Z_s,Z_t)
            xjac = (R_s * Z_t - R_t * Z_s)
            ! --- B-variables
            AR0_p  = 0.d0 ; AR0_R  = 0.d0 ; AR0_Z  = 0.d0
            AZ0_p  = 0.d0 ; AZ0_R  = 0.d0 ; AZ0_Z  = 0.d0
            A30_p  = 0.d0 ; A30_R  = 0.d0 ; A30_Z  = 0.d0
            do i_tor = 1, n_tor
              call interp(node_list,element_list,i_elm,var_AR, i_tor,s_out,t_out,AR0,AR0_s,AR0_t,AR0_st,AR0_ss,AR0_tt)
              call interp(node_list,element_list,i_elm,var_AZ, i_tor,s_out,t_out,AZ0,AZ0_s,AZ0_t,AZ0_st,AZ0_ss,AZ0_tt)
              call interp(node_list,element_list,i_elm,var_A3, i_tor,s_out,t_out,A30,A30_s,A30_t,A30_st,A30_ss,A30_tt)
              AR0_p  = AR0_p  + AR0 * HHZ_p(i_tor)
              AZ0_p  = AZ0_p  + AZ0 * HHZ_p(i_tor)
              A30_p  = A30_p  + A30 * HHZ_p(i_tor)
              if ((xjac .gt. 1.d-6)) then  ! avoid the axis
                AR0_R  = AR0_R  + (   Z_t * AR0_s - Z_s * AR0_t ) / xjac * HHZ(i_tor)
                AR0_Z  = AR0_Z  + ( - R_t * AR0_s + R_s * AR0_t ) / xjac * HHZ(i_tor)
                AZ0_R  = AZ0_R  + (   Z_t * AZ0_s - Z_s * AZ0_t ) / xjac * HHZ(i_tor)
                AZ0_Z  = AZ0_Z  + ( - R_t * AZ0_s + R_s * AZ0_t ) / xjac * HHZ(i_tor)
                A30_R  = A30_R  + (   Z_t * A30_s - Z_s * A30_t ) / xjac * HHZ(i_tor)
                A30_Z  = A30_Z  + ( - R_t * A30_s + R_s * A30_t ) / xjac * HHZ(i_tor)
              endif
            enddo

            ! --- Magnetic field
            call interp(node_list,element_list,i_elm,710,1,s_out,t_out,Fprof,Fprof_s,Fprof_t,Fprof_st,Fprof_ss,Fprof_tt)
            BR = ( A30_Z - AZ0_p )/ R
            BZ = ( AR0_p - A30_R )/ R
            Bp = ( AZ0_R - AR0_Z ) + Fprof / R
            
            R_half = R_line + R_out * delta_phi/2. / Bp * BR
            Z_half = Z_line + R_out * delta_phi/2. / Bp * BZ
          else
            cycle
          endif
          
          ! -----------------
          ! --- full step ---
          ! -----------------
          call find_RZ(node_list,element_list,R_half,Z_half,R_out,Z_out,i_elm,s_out,t_out,ifail)
          if ( ifail .eq. 0 ) then
            ! --- RZ-coords
            call interp_RZ(node_list,element_list,i_elm,s_out,t_out,R,R_s,R_t,Z,Z_s,Z_t)
            xjac = (R_s * Z_t - R_t * Z_s)
            ! --- B-variables
            AR0_p  = 0.d0 ; AR0_R  = 0.d0 ; AR0_Z  = 0.d0
            AZ0_p  = 0.d0 ; AZ0_R  = 0.d0 ; AZ0_Z  = 0.d0
            A30_p  = 0.d0 ; A30_R  = 0.d0 ; A30_Z  = 0.d0
            do i_tor = 1, n_tor
              call interp(node_list,element_list,i_elm,var_AR, i_tor,s_out,t_out,AR0,AR0_s,AR0_t,AR0_st,AR0_ss,AR0_tt)
              call interp(node_list,element_list,i_elm,var_AZ, i_tor,s_out,t_out,AZ0,AZ0_s,AZ0_t,AZ0_st,AZ0_ss,AZ0_tt)
              call interp(node_list,element_list,i_elm,var_A3, i_tor,s_out,t_out,A30,A30_s,A30_t,A30_st,A30_ss,A30_tt)
              AR0_p  = AR0_p  + AR0 * HHZ_p(i_tor)
              AZ0_p  = AZ0_p  + AZ0 * HHZ_p(i_tor)
              A30_p  = A30_p  + A30 * HHZ_p(i_tor)
              if ((xjac .gt. 1.d-6)) then  ! avoid the axis
                AR0_R  = AR0_R  + (   Z_t * AR0_s - Z_s * AR0_t ) / xjac * HHZ(i_tor)
                AR0_Z  = AR0_Z  + ( - R_t * AR0_s + R_s * AR0_t ) / xjac * HHZ(i_tor)
                AZ0_R  = AZ0_R  + (   Z_t * AZ0_s - Z_s * AZ0_t ) / xjac * HHZ(i_tor)
                AZ0_Z  = AZ0_Z  + ( - R_t * AZ0_s + R_s * AZ0_t ) / xjac * HHZ(i_tor)
                A30_R  = A30_R  + (   Z_t * A30_s - Z_s * A30_t ) / xjac * HHZ(i_tor)
                A30_Z  = A30_Z  + ( - R_t * A30_s + R_s * A30_t ) / xjac * HHZ(i_tor)
              endif
            enddo

            ! --- Magnetic field
            call interp(node_list,element_list,i_elm,710,1,s_out,t_out,Fprof,Fprof_s,Fprof_t,Fprof_st,Fprof_ss,Fprof_tt)
            BR = ( A30_Z - AZ0_p )/ R
            BZ = ( AR0_p - A30_R )/ R
            Bp = ( AZ0_R - AR0_Z ) + Fprof / R
            
            R_line = R_line + R_out * delta_phi / Bp * BR
            Z_line = Z_line + R_out * delta_phi / Bp * BZ
          else
            cycle
          endif
          
          ! -----------------
          ! --- record    ---
          ! -----------------
          call find_RZ(node_list,element_list,R_line,Z_line,R_out,Z_out,i_elm,s_out,t_out,ifail)
          if ( ifail /= 0 ) cycle
          R_turn(i_turn+1,(i_dir+1)/2+1) = R_line
          Z_turn(i_turn+1,(i_dir+1)/2+1) = Z_line
          call var_value(i_elm, var_T,   s_out,t_out,P_start, T_turn  (1,(i_dir+1)/2+1) )
          call var_value(i_elm, var_A3,  s_out,t_out,P_start, PSI_turn(1,(i_dir+1)/2+1) )
          call var_value(i_elm, var_rho, s_out,t_out,P_start, ZN_turn (1,(i_dir+1)/2+1) )
          
          phi = phi + delta_phi
          
        end do
 
      end do

    end do

    if (     ( (FF_0 .gt. 0.d0) .and. (PSI_turn(1,1) .lt. psi_bnd2) ) &
        .or. ( (FF_0 .lt. 0.d0) .and. (PSI_turn(1,1) .gt. psi_bnd2) ) )then

      ! --- Compute averaged connection length
      count_lines = count_lines + 1
      C_average = C_average + C_all(i_line)
      
      ! --- Do both directions
      do i_dir=1,2
  
        ! --- Do all turns
        do i_turn=1,n_turn_max(i_dir)+1
  
          if (R_turn(i_turn,i_dir) .gt. 0.d0) then
  
            ! --- Connection lengths in both directions
            zl1 = C_turn(i_turn,i_dir)
            zl2 = C_turn(1,1) + C_turn(1,2) - C_turn(i_turn,i_dir)
  
            ikeep = ikeep + 1
            small_r   = sqrt( (R_turn(i_turn,i_dir)-R_axis)**2 + (Z_turn(i_turn,i_dir)-Z_axis)**2 )
            theta_pol = atan2(Z_turn(i_turn,i_dir)-Z_axis,R_turn(i_turn,i_dir)-R_axis)
            if (theta_pol .lt. 0.d0) theta_pol = theta_pol + 2.d0*PI
            if (ikeep .le. n_points_max) then
              if (Rtheta_plot) then
                call find_RZ(node_list,element_list,R_turn(i_turn,i_dir),Z_turn(i_turn,i_dir),R_tmp,Z_tmp,ielm_tmp,s_tmp,t_tmp,ifail)
                call interp(node_list,element_list,ielm_tmp,1,1,s_tmp,t_tmp,psi_tmp,P0_s,P0_t,P0_st,P0_ss,P0_tt)
                RZkeep(1,ikeep)          = (psi_tmp-psi_axis)/(psi_bnd-psi_axis)!small_r
                RZkeep(2,ikeep)          = theta_pol / (2.d0*PI)
              else
                RZkeep(1,ikeep)          = R_turn(i_turn,i_dir)
                RZkeep(2,ikeep)          = Z_turn(i_turn,i_dir)
              endif
              !scalars(ikeep,1:n_scalars) = (/ min(zl1,zl2),T_turn(1,i_dir) /)
              !scalars(ikeep,1:n_scalars) = (/ C_all_4(i_line),C_minus(i_line),C_plus(i_line),T_turn(1,i_dir) /)
              scalars(ikeep,1:n_scalars) = (/ C_all(i_line),min(zl1,zl2),PSI_turn(1,i_dir),T_turn(1,i_dir) /)
            else
              write(*,*) 'Warning! Exceeded maximum number of points',n_points_max
            endif
  
          endif
        enddo
      
      enddo ! end recording in both directions
  
    endif

  end do

  		      else ! crossing an outer boundary
  
  			R_strike(i_strike) = R_in
  			Z_strike(i_strike) = Z_in
  			P_strike(i_strike) = p_line
  			C_strike(i_strike) = total_length + sqrt(abs(dl2))
  
  			call var_value(i_elm_prev,6,s_line,t_line,p_line,T_strike(i_strike))
  			call var_value(i_elm_prev,5,s_line,t_line,p_line,ZN_strike(i_strike))
  
  			inode1 = element_list%element(i_elm_prev)%vertex(1)
  			inode2 = element_list%element(i_elm_prev)%vertex(2)
  
  			if ((node_list%node(inode1)%boundary .ne. 0) .and. (node_list%node(inode2)%boundary .ne. 0)) then
  			  B_strike(i_strike) = min(node_list%node(inode1)%boundary,node_list%node(inode2)%boundary)
  			else
  			  write(*,*) 'error : leaving domain but not at a boundary (t=0)!',inode1,inode2,R_in,Z_in
  			endif
  
  		      endif
  
  		    endif
  
  		  endif
  
  		else  ! this step remains within the element
  
  		  s_line = s_line + delta_s
  		  t_line = t_line + delta_t
  		  p_line = p_line + delta_phi_step
  
  		  dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
  
  		  small_delta = 1.d0
  
  		  call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)
  
  		endif
  
  		! --- Reset the toroial step
		delta_phi_local = delta_phi_local + small_delta * delta_phi_step
  
  		! --- Record total lengths
  		total_length = total_length + sqrt(abs(dl2))
  		total_phi    = total_phi    + small_delta * delta_phi_step
  
  		! --- Exit if we stepped out of domain
  		if (i_elm .eq. 0) exit
  
  	      enddo  ! end of loop over steps within one element
  
              ! --- Exit if we stepped out of domain
  	      if (i_elm .eq. 0) exit
  
  	    enddo    ! end of a 2Pi turn (or before if end of open field line)
  
            ! --- Exit if we stepped out of domain
  	    if (i_elm .eq. 0) exit
  
  	    ! --- Record Poincare point
	    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)
  	    R_turn(i_turn+1,(i_dir+1)/2+1) = R_in
  	    Z_turn(i_turn+1,(i_dir+1)/2+1) = Z_in
  	    C_turn_tmp(i_turn+1,(i_dir+1)/2+1) = total_length
  	    call var_value(i_elm, 6, s_line,t_line,p_line, T_turn  (i_turn+1,(i_dir+1)/2+1) )
  	    call var_value(i_elm, 1, s_line,t_line,p_line, PSI_turn(i_turn+1,(i_dir+1)/2+1) )
  	    call var_value(i_elm, 5, s_line,t_line,p_line, ZN_turn (i_turn+1,(i_dir+1)/2+1) )
  
  	  enddo  ! end of loop over toroidal turns
  
  	  ! --- Field line still in domain, after n_turn turns
	  if (i_elm .ne. 0) then
  	    R_strike(i_strike) = R_in
  	    Z_strike(i_strike) = Z_in
  	    P_strike(i_strike) = p_line
  	    C_strike(i_strike) = 0.d0	       ! to be done (total_length needs correction)
  	    B_strike(i_strike) = 0
  	    call var_value(i_elm,6,s_line,t_line,p_line,T_strike(i_strike))
  	    call var_value(i_elm,5,s_line,t_line,p_line,ZN_strike(i_strike))
  	  endif
  
  	  ! --- Record connection length
  	  if (i_dir .eq. -1) then  
  	    C_all(i_line)   = total_length
  	    C_minus(i_line) = total_length
  	    C_all_4(i_line) = total_length
  	    partial(1)      = total_length
  	  else
            C_all(i_line)   = C_all(i_line)+total_length!min(C_all(i_line),total_length)
  	    C_plus(i_line)  = total_length
  	    C_all_4(i_line) = C_all_4(i_line)+total_length
  	    partial(2)      = total_length
  	  endif
  
  	enddo  ! end of two directions
  
        ! -----------------------------------
	! --- Record what we want in vtk file
        ! -----------------------------------
	
        ! --- Reverse connection length (not from plasma to position, but from target to position)
	do i_dir=1,2
  	  do i_turn = 1, n_turn_max(i_dir)
  	    C_turn(i_turn,i_dir) = partial(i_dir) - c_turn_tmp(i_turn,i_dir)
  	  enddo
  	enddo
  
        ! --- Keep only field lines starting inside the plasma
  	if ( (PSI_turn(1,1) .lt. psi_bnd2) &
	     .and. (    ((n_turn_max(1) .lt. n_turns) .and. (n_turn_max(2) .lt. n_turns)) &
	            .or. Rtheta_plot ) &
  	     .and. (	( (xcase .eq. LOWER_XPOINT) .and. (Z_turn(1,1) .gt. ES%Z_xpoint(1)) ) &
  	  	    .or.( (xcase .eq. UPPER_XPOINT) .and. (Z_turn(1,1) .lt. ES%Z_xpoint(2)) ) &
  	  	    .or.( (xcase .eq. DOUBLE_NULL ) .and. (Z_turn(1,1) .gt. ES%Z_xpoint(1)) .and. (Z_turn(1,1) .lt. ES%Z_xpoint(2)) ) ) ) then
  
          ! --- Compute averaged connection length
	  count_lines = count_lines + 1
	  C_average = C_average + C_all(i_line)
	  
	  ! --- Do both directions
	  do i_dir=1,2
  
            ! --- Do all turns
  	    do i_turn=1,n_turn_max(i_dir)+1
  
  	      if (R_turn(i_turn,i_dir) .gt. 0.d0) then
  
  	  	! --- Connection lengths in both directions
	  	zl1 = C_turn(i_turn,i_dir)
  	  	zl2 = C_turn(1,1) + C_turn(1,2) - C_turn(i_turn,i_dir)
  
  	  	ikeep = ikeep + 1
		small_r   = sqrt( (R_turn(i_turn,i_dir)-ES%R_axis)**2 + (Z_turn(i_turn,i_dir)-ES%Z_axis)**2 )
		theta_pol = atan2(Z_turn(i_turn,i_dir)-ES%Z_axis,R_turn(i_turn,i_dir)-ES%R_axis)
		if (theta_pol .lt. 0.d0) theta_pol = theta_pol + 2.d0*PI
  	  	if (ikeep .le. n_points_max) then
		  if (Rtheta_plot) then
                    call find_RZ(node_list,element_list,R_turn(i_turn,i_dir),Z_turn(i_turn,i_dir),R_tmp,Z_tmp,ielm_tmp,s_tmp,t_tmp,ifail)
                    call interp(node_list,element_list,ielm_tmp,1,1,s_tmp,t_tmp,psi_tmp,P0_s,P0_t,P0_st,P0_ss,P0_tt)
		    RZkeep(1,ikeep)	     = (psi_tmp-ES%psi_axis)/(psi_bnd-ES%psi_axis)!small_r
  	  	    RZkeep(2,ikeep)	     = theta_pol / (2.d0*PI)
		  else
		    RZkeep(1,ikeep)	     = R_turn(i_turn,i_dir)
  	  	    RZkeep(2,ikeep)	     = Z_turn(i_turn,i_dir)
		  endif
          	  !scalars(ikeep,1:n_scalars) = (/ min(zl1,zl2),T_turn(1,i_dir) /)
  	  	  !scalars(ikeep,1:n_scalars) = (/ C_all_4(i_line),C_minus(i_line),C_plus(i_line),T_turn(1,i_dir) /)
		  scalars(ikeep,1:n_scalars) = (/ C_all(i_line),min(zl1,zl2),PSI_turn(1,i_dir),T_turn(1,i_dir) /)
		else
		  write(*,*) 'Warning! Exceeded maximum number of points',n_points_max
		endif
  
  	      endif
  	    enddo
  	  
	  enddo ! end recording in both directions
  
        endif
  
      enddo  ! end over loop over starting points within one element ( ns)
    enddo    ! end over loop over starting points within one element ( nt)
  
  
  enddo ! end of loop over elements
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  write(*,*)'Finished loop!!!'
  
  ! --- Print averaged connection length
  call MPI_Reduce(count_lines,count_lines_tot,1,MPI_INTEGER,         MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_Reduce(C_average,  C_average_tot,  1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if (my_id .eq. 0) then
    if (count_lines_tot .gt. 0) then
      C_average_tot = C_average_tot / real(count_lines_tot)
      write(*,*)'Averaged connection length:',C_average_tot
    endif
  endif
  

  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  ! --- Third part: write poincare plot to VTK file
  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  
  write(*,*)'Writing connection.vtk'
  
  ! --- An arbitrary unit number for the VTK output file
  ivtk = 22
  
  ! --- The local number of points (for mpi_0)
  ikeep0  = ikeep
  call MPI_Reduce(ikeep,nnos,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if (my_id .eq. 0) write(*,*) ' number of points : ',nnos
  
  ! --- Line feed character
  lf = char(10)
  
  ! --- Open file and write headers
  if (my_id .eq. 0) then
#ifdef IBM_MACHINE
    open(unit=ivtk,file='connection.vtk',form='unformatted',access='stream',status='replace')
#else
    open(unit=ivtk,file='connection.vtk',form='unformatted',access='stream',convert='BIG_ENDIAN',status='replace')
#endif
    buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
    buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
    buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
    buffer = 'DATASET UNSTRUCTURED_GRID'//lf//lf                                          ; write(ivtk) trim(buffer)
    write(str1(1:12),'(i12)') nnos
    buffer = 'POINTS '//str1//'  float'//lf                                               ; write(ivtk) trim(buffer)
  endif
  
  ! --- Write points for local MPI (id=0)
  if (my_id .eq. 0) then
    write(ivtk) ( (/RZkeep(1,i), RZkeep(2,i), ZERO /),i=1,ikeep0)
  endif
  
  ! --- Write points for all other MPIs
  if (my_id .eq. 0) then
    ! --- If this is mpi_0, we receive data from the other MPIs and print it
    do j=1,n_cpu-1
      call mpi_recv(ikeep,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      if (ikeep .gt. 0) then
        nrecv = 2*ikeep
        call mpi_recv(RZkeep,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        write(ivtk) ( (/RZkeep(1,i), RZkeep(2,i), ZERO /),i=1,ikeep)
      endif
    enddo
  else
    ! --- If this is not mpi_0, we send data to the main MPI 0
    call mpi_send(ikeep, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    if (ikeep .gt. 0) then
      nsend = 2*ikeep
      call mpi_send(RZkeep, nsend,MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
    endif
  endif
  
  ! --- Now for the variables on points (headers first)
  if (my_id .eq. 0) then
    write(str1(1:12),'(i12)') nnos
    buffer = lf//lf//'POINT_DATA '//str1//lf                                              ; write(ivtk) trim(buffer)
  endif
  
  ! --- Loop over all variables
  do i_var =1, n_scalars

    ! --- More headers
    if (my_id .eq. 0) then
      buffer = 'SCALARS '//scalar_names(i_var)//' float'//lf                              ; write(ivtk) trim(buffer)
      buffer = 'LOOKUP_TABLE default'//lf                                                 ; write(ivtk) trim(buffer)
    endif
  
    ! --- Write data for local MPI (id=0)
    if (my_id .eq. 0) then
      write(ivtk) (scalars(i,i_var),i=1,ikeep0)
    endif
  
    ! --- Write data for all other MPIs
    if (my_id .eq. 0) then
      ! --- If this is mpi_0, we receive data from the other MPIs and print it
      do j=1,n_cpu-1
        call mpi_recv(ikeep,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
        if (ikeep .gt. 0) then
          nrecv = ikeep
          call mpi_recv(scalars(1:ikeep,i_var),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
          write(ivtk) (scalars(i,i_var),i=1,ikeep)
        endif
      enddo
    else
      ! --- If this is not mpi_0, we send data to the main MPI 0
      call mpi_send(ikeep, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
      if (ikeep .gt. 0) then
        nsend = ikeep
        call mpi_send(scalars(1:ikeep,i_var), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      endif
    endif
  
  enddo
  
  close(ivtk)
  
  deallocate(RZkeep,scalars,scalar_names)
  
  
  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  ! --- Fourth part: write strike points plot to VTK file
  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  
  write(*,*)'Writing strikes.vtk'
  
  ! --- Open text file
  open(23,file='strikes.txt')
  
  ! --- Clean data
  do i=1,i_strike
    if (abs(R_strike(i)) .gt. 10.d0) R_strike(i) = 0.d0
    if (abs(Z_strike(i)) .gt. 10.d0) Z_strike(i) = 0.d0
    ! --- Exclude points started outside the plasma
    if (PS0_strike(i) .gt. psi_bnd2) then
      R_strike(i) = 0.d0
      Z_strike(i) = 0.d0
    endif
  enddo
  
  ! --- The local number of strikes (for mpi_0)
  i_strike0 = i_strike
  call MPI_Reduce(i_strike,nnos,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if (my_id .eq. 0) write(*,*) ' number of points : ',nnos
  
  ! --- Reallocate data for strike points
  n_scalars = 3
  allocate(scalar_names(n_scalars),scalars(100000,n_scalars))
  scalar_names  = (/ 'length   ','psi_start','T_start  '/)
  
  ! --- Copy data
  do i=1,i_strike
    scalars(i,1) = C_strike(i)                 ! needs correction !!! see above
    scalars(i,2) = PS0_strike(i)
    scalars(i,3) = T0_strike(i)
  enddo
  
  ! --- Open file and write headers
  if (my_id .eq. 0) then
#ifdef IBM_MACHINE
    open(unit=ivtk,file='strikes.vtk',form='unformatted',access='stream')
#else
    open(unit=ivtk,file='strikes.vtk',form='unformatted',access='stream',convert='BIG_ENDIAN')
#endif
    buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
    buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
    buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
    buffer = 'DATASET UNSTRUCTURED_GRID'//lf//lf                                          ; write(ivtk) trim(buffer)
    write(str1(1:12),'(i12)') nnos
    buffer = 'POINTS '//str1//'  float'//lf                                               ; write(ivtk) trim(buffer)
  endif
  
  ! --- Write points for local MPI (id=0)
  if (my_id .eq. 0) then
    write(ivtk) ( (/R_strike(i)*cos(P_strike(i)), Z_strike(i), R_strike(i)*sin(P_strike(i)) /),i=1,i_strike0)
  !!!  write(23,'(4e16.8)') ( (/ R_strike(i)*cos(P_strike(i)), Z_strike(i), R_strike(i)*sin(P_strike(i)), C_all_4(i) /),i=1,i_strike0)
    write(23,'(3e16.8)') ( (/ R_strike(i), Z_strike(i), C_all_4(i) /),i=1,i_strike0)
  endif
  
  ! --- Write points for all other MPIs
  if (my_id .eq. 0) then
    ! --- If this is mpi_0, we receive data from the other MPIs and print it
    do j=1,n_cpu-1
      call mpi_recv(i_strike,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      if (i_strike .gt. 0) then
        nrecv = i_strike
        call mpi_recv(R_strike,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        call mpi_recv(Z_strike,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        call mpi_recv(P_strike,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        write(ivtk) ( (/R_strike(i)*cos(P_strike(i)), Z_strike(i), R_strike(i)*sin(P_strike(i)) /),i=1,i_strike)
  !!!      write(23,'(4e16.8)') ( (/R_strike(i)*cos(P_strike(i)), Z_strike(i), R_strike(i)*sin(P_strike(i)), C_all_4(i) /),i=1,i_strike)
        write(23,'(3e16.8)') ( (/R_strike(i), Z_strike(i), C_all_4(i) /),i=1,i_strike)
      endif
    enddo
  else
    ! --- If this is not mpi_0, we send data to the main MPI 0
    call mpi_send(i_strike, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    if (i_strike .gt. 0) then
      nsend = i_strike
      call mpi_send(R_strike, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      call mpi_send(Z_strike, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      call mpi_send(P_strike, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
    endif
  endif
  
  ! --- Now for the variables on points (headers first)
  if (my_id .eq. 0) then
    write(str1(1:12),'(i12)') nnos
    buffer = lf//lf//'POINT_DATA '//str1//lf                                              ; write(ivtk) trim(buffer)
  endif
  
  ! --- Loop over all variables
  do i_var =1, n_scalars
  
    ! --- More headers
    if (my_id .eq. 0) then
      buffer = 'SCALARS '//scalar_names(i_var)//' float'//lf                              ; write(ivtk) trim(buffer)
      buffer = 'LOOKUP_TABLE default'//lf                                                 ; write(ivtk) trim(buffer)
    endif

    ! --- Write data for local MPI (id=0)
    if (my_id .eq. 0) then
      write(ivtk) (scalars(i,i_var),i=1,i_strike0)
    endif

    ! --- Write data for all other MPIs
    if (my_id .eq. 0) then
      ! --- If this is mpi_0, we receive data from the other MPIs and print it
      do j=1,n_cpu-1
        call mpi_recv(i_strike,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
        if (i_strike .gt. 0) then
          nrecv = i_strike
          call mpi_recv(scalars(1:i_strike,i_var),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
          write(ivtk) (scalars(i,i_var),i=1,i_strike)
        endif
      enddo
    else
      ! --- If this is not mpi_0, we send data to the main MPI 0
      call mpi_send(i_strike, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
      if (i_strike .gt. 0) then
        nsend = i_strike
        call mpi_send(scalars(1:i_strike,i_var), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      endif
    endif

  enddo

  ! Clean and exit
  close(ivtk)
  close(23)

  call MPI_FINALIZE(IERR)                                ! clean up MPI

end




































! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------








subroutine var_value(i_elm,i_var,s_in,t_in,p_in,value_out)
  use mod_parameters
  use elements_nodes_neighbours
  use phys_module
  use mod_interp
  
  implicit none
  
  integer :: i_var, i_elm, i_tor, i_harm
  
  real*8 :: s_in, t_in, p_in
  real*8 :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt, Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
  real*8 :: P0,P0_s,P0_t,P0_st,P0_ss,P0_tt
  real*8 :: value_out

  !call interp_RZ(node_list,element_list,i_elm,s_in,t_in,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
  !Zjac = (R_s * Z_t - R_t * Z_s)
  
  call interp(node_list,element_list,i_elm,i_var,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
  
  value_out = P0
  
  do i_tor = 1, (n_tor-1)/2
  
    i_harm = 2*i_tor
  
    call interp(node_list,element_list,i_elm,i_var,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)
  
    value_out = value_out + Pcos * cos(mode(i_harm)*p_in)
  
    call interp(node_list,element_list,i_elm,i_var,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)
  
    value_out = value_out + Psin * sin(mode(i_harm+1)*p_in)
  
  enddo

  return
end
