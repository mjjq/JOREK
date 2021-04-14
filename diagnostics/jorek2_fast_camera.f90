!> Program to convert a JOREK2 restart file into binary VTK format
program jorek2_fast_camera

  use mod_parameters, only: n_var, variable_names, jorek_model
  use mod_import_restart
  use data_structure
  use phys_module
  use basis_at_gaussian
  use diffusivities, only: get_dperp, get_zkperp
  use mod_interp

  implicit none
  !include 'mpif.h'

  type (type_node_list)   , pointer :: node_list
  type (type_element_list), pointer :: element_list
  type (type_surface_list)          :: flux_list

  integer               :: nnoel, nnos, nel, nsub, inode, ielm, n_scalars, n_vectors
  real*4,allocatable    :: xyz (:,:), scalars(:,:), vectors(:,:,:)
  integer,allocatable   :: ien (:,:)
  integer, parameter    :: ivtk = 22 ! an arbitrary unit number for the VTK output file
  integer               :: i, j, k, m, etype, irst, int, i_var, i_tor, i_tor_old, i_plane, index, index_node
  character             :: buffer*80, lf*1, str1*12, str2*12
  character*12, allocatable :: scalar_names(:), vector_names(:)
  real*8                :: s, t
  real*8                :: P,P_s,P_t,P_st,P_ss,P_tt
  real*8                :: R,R_s,R_t,R_st,R_ss,R_tt
  real*8                :: Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
  real*8                :: PPPsi,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt
  real*8                :: ZJ,ZJ_s,ZJ_t,ZJ_st,ZJ_ss,ZJ_tt
  real*8                :: U,U_s,U_t,U_st,U_ss,U_tt
  real*8                :: W,W_s,W_t,W_st,W_ss,W_tt
  real*8                :: RRRHO,RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt
  real*8                :: TTT,TT_s,TT_t,TT_st,TT_ss,TT_tt
  real*8                :: Ti,Ti_s,Ti_t,Ti_st,Ti_ss,Ti_tt
  real*8                :: TTTe,Te_s,Te_t,Te_st,Te_ss,Te_tt
  real*8                :: V, V_s, V_t, V_st, V_ss, V_tt
  real*8                :: psi_00, rho_00, Ti_00, Te_00
  real*8                :: ps_x, ps_y
  real*8                :: u0_x, u0_y
  real*8                :: zj_x, zj_y
  real*8                :: w0_x, w0_y, w0_xx, w0_yy
  real*8                :: RHO_x, RHO_y, RHO_p
  real*8                :: TT_x, TT_y, TT_p
  real*8                :: Ti_x, Ti_y, Ti_p
  real*8                :: Te_x, Te_y, Te_p
  real*8                :: psi_axis,      R_axis,      Z_axis,      s_axis,      t_axis
  real*8                :: psi_xpoint(2), R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2)
  real*8                :: ps0, psi_norm, psi_bnd, grad_psi
  real*8                :: xjac, xjac_x, xjac_y, v_perp, Psi_J, R_p, error, Btot, BigR
  real*8                :: particle_source, D_prof, ZK_prof, source_pellet, ZKpar_T
  integer               :: i_find, i_elm_find(8)
  real*8                :: Router,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
  real*8                :: Zouter,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
  real*8                :: s_find(8), t_find(8)
  real*8                :: Jb
  real*8                :: central_ne
  integer               :: i_elm_axis, i_elm_xpoint(2), k_tor
  logical               :: without_n0_mode
  !====================== --- add the diagnostics Er, Vtheta and [not yet Vneo]
  real*8                :: Er, psi_abs, Vtheta, Btheta, Mach_par,Mach_pol,Vsound
  


  ! --- MPI variables
  integer               :: my_id, n_cpu, ierr, nsend, nrecv
  integer               :: status(MPI_STATUS_SIZE)
  integer               :: pix_start, pix_end, pix_delta
  ! --- Photon Emissivity Coeff (PEC) variables
  integer               :: PEC_size, PEC_index_Ne, PEC_index_Te, PEC_index
  character*44          :: PEC_file
  character*50          :: line
  real*8,allocatable    :: PEC_dens(:), PEC_temp(:), PEC(:)
  ! --- Camera variables
  integer               :: n_pixels_hor, n_pixels_ver, n_pix, ncount
  real*8                :: pixel_dim, focus, vec_size
  real*8                :: X_cam, Y_cam, Z_cam
  real*8,allocatable    :: Xp(:), Yp(:), Zp(:)
  real*8,allocatable    :: Xv(:), Yv(:), Zv(:)
  real*8,allocatable    :: Light(:)
  character*1           :: rgb(3)
  integer               :: itmp, icnt
  ! --- Integration variables
  integer               :: i_elm, ifail, inside
  integer               :: n_step
  real*8                :: step
  real*8                :: X_tmp, Y_tmp, Z_tmp
  real*8                :: RR,    ZZ,    Phi
  real*8                :: R_out, Z_out
  real*8                :: ss,    tt
  real*8,allocatable    :: HZ_tor(:)
  real*8                :: psi, rho, rho_n, Te, eV2Joules, solenoid
  real*8                :: tanh_psi, tanh_zmin, tanh_zpls

  real*8                :: maxlight
  
  
  


  ! ********************************************************************** !
  !                              Start MPI                                 !
  ! ********************************************************************** !  
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs



  if (my_id .eq. 0) write(*,*) '/**************************************/'
  if (my_id .eq. 0) write(*,*) '/*******  jorek2_fast_camera  *********/'
  if (my_id .eq. 0) write(*,*) '/**************************************/'
  allocate(node_list)
  allocate(element_list)

  ! ********************************************************************** !
  !  Read data file with PEC(Ne,Te) grid (Photon Emissivity Coefficients)  !
  ! ********************************************************************** !
  
  !PEC_file = "/scratch/pstanis/MAST/fast_camera/my_pec.dat"
  PEC_file = "./my_pec.dat"
  open(123, file=PEC_file, action='read', iostat=ierr)
  
  if ( ierr .ne. 0 ) then
    if (my_id .eq. 0) write(*,*) 'Failed to open PEC data file ',PEC_file,ierr
    if (my_id .eq. 0) write(*,*) 'You can find the my_pec.dat file in jorek/util/my_pec.dat'
    if (my_id .eq. 0) write(*,*) 'Aborting...'
    call MPI_FINALIZE(ierr)
    stop
  else
    
    ! --- File should start with a comment line
    read(123,'(A)') line
    
    ! --- Second line should be "np"
    read(123,'(A)') line
    
    ! --- Third line should be the value of "np"
    read(123,'(A)') line
    read(line,*) PEC_size
    
    ! --- Allocate vectors
    allocate(PEC_dens(PEC_size),PEC_temp(PEC_size),PEC(PEC_size*PEC_size))
    
    ! --- Then comes the density profile
    read(123,'(A)') line
    do i=1, PEC_size
      read(123,'(A)') line
      read(line,*) PEC_dens(i)
    enddo
    
    ! --- Then comes the temperature profile
    read(123,'(A)') line
    do i=1, PEC_size
      read(123,'(A)') line
      read(line,*) PEC_temp(i)
    enddo
    
    ! --- Then comes the PEC profiles
    read(123,'(A)') line
    do i=1, PEC_size*PEC_size
      read(123,'(A)') line
      read(line,*) PEC(i)
    enddo
    
    close(123)
  endif
    
    
    
    

  ! ********************************************************************** !
  !                     Build lines of sight for camera                    !
  ! ********************************************************************** !
  
  ! --- MAST low-res
  ! --- Number of horizontal and vertical pixels
  n_pixels_hor = 230!640!320!
  n_pixels_ver = 180!368!240!
  ! --- Pixel dimension (assumed to be a square)
  pixel_dim = 6.8e-05!17.d-6!4.6d-5!
  
  ! --- MAST high-res
  ! --- Number of horizontal and vertical pixels
  n_pixels_hor = 640!230!320!
  n_pixels_ver = 500!180!240!
  ! --- Pixel dimension (assumed to be a square)
  pixel_dim = 17.d-6!6.8e-05!4.6d-5!
  
  ! --- Location of camera focus (R,Z)
  X_cam = 2.05
  Y_cam = 0.0 
  Z_cam = 0.0 
  
  ! --- Focal length of camera lens
  focus = 4.8d-3
  
  ! --- Allocate the lines of sight coords (one Point plus one Vector plus Light intensity)
  n_pix = n_pixels_hor*n_pixels_ver
  allocate(Xp(n_pix),Yp(n_pix),Zp(n_pix))
  allocate(Xv(n_pix),Yv(n_pix),Zv(n_pix))
  allocate(Light(n_pix))
  Light  = 0.d0
  ncount = 0
  
  ! --- Calculate each Point/Vector
  ! --- (note we assume the sensor is "in front" of lens, otherwise you need to flip the image afterwards)
  do i=1, n_pixels_hor
    do j=1, n_pixels_ver
      ncount     = ncount+1
      Xp(ncount) = X_cam - focus
      Yp(ncount) = Y_cam + pixel_dim*(n_pixels_ver-1)/2 - pixel_dim*(j-1)
      Zp(ncount) = Z_cam - pixel_dim*(n_pixels_hor-1)/2 + pixel_dim*(i-1)
      Xv(ncount) = Xp(ncount) - X_cam
      Yv(ncount) = Yp(ncount) - Y_cam
      Zv(ncount) = Zp(ncount) - Z_cam
      vec_size   = ( Xv(ncount)**2.d0 + Yv(ncount)**2.d0 + Zv(ncount)**2.d0 )**0.5d0
      Xv(ncount) = Xv(ncount)/vec_size
      Yv(ncount) = Yv(ncount)/vec_size
      Zv(ncount) = Zv(ncount)/vec_size
    enddo
  enddo
  
  
  
  ! ********************************************************************** !
  !       Initialise, import restart, and find axis and xpoint             !
  ! ********************************************************************** !
  
  ! --- Import and initialise
  call initialise_and_broadcast_parameters(my_id, "__NO_FILENAME__")
  do k_tor=1, n_tor
    mode(k_tor) = + int(k_tor / 2) * n_period
  enddo
  call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)
  call initialise_basis
  call broadcast_elements(my_id, element_list)                ! elements
  call broadcast_nodes(my_id, node_list)                      ! nodes
  call broadcast_phys(my_id)                                  ! physics parameters

  
  ! --- Find axis
  call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

  ! --- Find Xpoint
  if (xpoint) then
    call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
    psi_bnd  = psi_xpoint(1)
    if( (xcase .eq. 2) .or. ((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
      psi_bnd = psi_xpoint(2)
    endif
  else
    psi_bnd = 0.d0
  endif
  
  if (central_density .gt. 1.d10) then
    central_ne = central_density
  else
    central_ne = central_density*1.d20
  endif

  
  
  ! ********************************************************************** !
  !                 Integrate radiation on each line of sight              !
  ! ********************************************************************** !
  
  ! --- MPI loops
  pix_delta = n_pix / n_cpu
  pix_start = my_id*pix_delta + 1
  pix_end   = min(n_pix,(my_id+1)*pix_delta)
  pix_delta = pix_end-pix_start+1

  ! --- The step size for integration and number of steps to get to the other side of plasma
  step   = 1.d-2!5.d-3
  n_step = 2.d0*(2.0*X_cam/step)
  
  ! --- Need the harmonic contributions for the toroidal location
  do i_tor=1, n_tor
    mode(i_tor) = int(i_tor / 2) * n_period
  enddo
  allocate(HZ_tor(n_tor))
  
  ! --- For each line of sight...
  do i=pix_start, pix_end
    
    !write(*,'(A,i6,A,i6)')'n_pix : ',i,' out of ',n_pix
    X_tmp  = Xp(i)
    Y_tmp  = Yp(i)
    Z_tmp  = Zp(i)
    ifail  = 1
    inside = 0  
    
    ! --- For each step on the line of sight...
    do j=1, n_step
      
      ! --- If we are entering the plasma core, increase step size
      if ( (ifail .eq. 0) .and. (psi .lt. 0.6) .and. (inside .eq. 0) ) then
        step   = 4.d0*step
        inside = 1      
      endif

      ! --- If we are exiting the plasma core, reduce step size
      if ( (ifail .eq. 0) .and. (psi .ge. 0.6) .and. (inside .eq. 1) ) then
        step   = step/4.d0
        inside = 0      
      endif
      
      ! --- X,Y,Z-coords
      X_tmp = X_tmp + step*Xv(i)
      Y_tmp = Y_tmp + step*Yv(i)
      Z_tmp = Z_tmp + step*Zv(i)
      
      ! --- R,Z,Phi-coords
      RR  = ( X_tmp**2.d0 + Z_tmp**2.d0 )**0.5d0
      ZZ  = Y_tmp
      Phi = atan2(Z_tmp,X_tmp)
      if (Phi .lt. 0.d0) Phi = Phi + 2.d0*PI
      
      ! --- If we hit the solenoid, stop integrating
      solenoid = 0.195
      if (RR .lt. solenoid) exit
      
      ! --- Look for this point in the simulation domain
      call find_RZ(node_list,element_list,RR,ZZ,R_out,Z_out,i_elm,ss,tt,ifail)
      
      ! --- Calculate emissivity at this point and add to integration
      if (ifail .eq. 0) then

        ! --- Need the harmonic contributions for the toroidal location
        HZ_tor(1)   = 1.d0
        do i_tor=1,(n_tor-1)/2
          HZ_tor(2*i_tor)     = cos(mode(2*i_tor)  *Phi)
          HZ_tor(2*i_tor+1)   = sin(mode(2*i_tor+1)*Phi)
        enddo

        ! --- Build variables
        psi   = 0.d0 
        rho   = 0.d0 
        rho_n = 0.d0 
        Te    = 0.d0 
        do i_tor=1,n_tor
          call interp(node_list,element_list,i_elm,1,i_tor,ss,tt,P,P_s,P_t,P_st,P_ss,P_tt)
          psi = psi + P * HZ_tor(i_tor)
          call interp(node_list,element_list,i_elm,5,i_tor,ss,tt,P,P_s,P_t,P_st,P_ss,P_tt)
          rho = rho + P * HZ_tor(i_tor)
          call interp(node_list,element_list,i_elm,6,i_tor,ss,tt,P,P_s,P_t,P_st,P_ss,P_tt)
          Te  = Te  + P * HZ_tor(i_tor)
          if ( ( jorek_model == 500 ) .or. ( jorek_model == 555 ) ) then
            call interp(node_list,element_list,i_elm,8,i_tor,ss,tt,P,P_s,P_t,P_st,P_ss,P_tt)
            rho_n = rho_n + P * HZ_tor(i_tor)
          endif
        enddo
        
        ! --- Normalise psi and denormalise density and temperature
        eV2Joules = 1.602176487d-19
        psi   = (psi-psi_axis)/(psi_bnd-psi_axis)
        rho   = rho*central_ne
        rho_n = rho_n*central_ne
        Te    = Te/(central_ne*MU_ZERO*eV2Joules)

        ! --- Neutral density (assumed 20% of core density outside plasma for MAST), if not using neutrals model
        if ( ( jorek_model .ne. 500 ) .and. ( jorek_model .ne. 555 ) ) then
          !tanh_psi  = 2.5d16 * (0.5 - 0.5* tanh((1.05 - psi)/0.025) ) + 2.5d12
          tanh_psi  = 2.5d16 * (0.5 - 0.5* tanh((0.99 - psi)/0.025) ) + 2.5d12
          tanh_psi  = 2.5d16 * (0.5 - 0.5* tanh((0.99 - psi)/0.025) ) + 2.5d12
          if (ZZ .lt. Z_xpoint(1)) tanh_psi  = 2.5d16
          if (ZZ .gt. Z_xpoint(2)) tanh_psi  = 2.5d16
          tanh_zmin = 3.d0   * (0.5 - 0.5* tanh(-(Z_xpoint(1) - ZZ)/0.1) ) + 1.d0
          tanh_zpls = 3.d0   * (0.5 - 0.5* tanh( (Z_xpoint(2) - ZZ)/0.1) ) + 1.d0
          tanh_zmin = 1.d1   * (0.5 - 0.5* tanh(-(Z_xpoint(1) - ZZ)/0.1) ) + 1.d0
          tanh_zpls = 1.d1   * (0.5 - 0.5* tanh( (Z_xpoint(2) - ZZ)/0.1) ) + 1.d0
          !rho_n = tanh_psi * tanh_zmin * tanh_zpls
          rho_n = 0.4d0*central_ne*( 0.5d0 - 0.5d0*tanh((0.975-psi)/0.05) ) + 0.01d0*central_ne
          rho_n = 0.4d0*central_ne*( 0.5d0 - 0.5d0*tanh((1.0-psi)/0.001) ) + 0.01d0*central_ne
          !rho_n = rho_n * tanh_zmin * tanh_zpls
        endif

        ! --- Calculate PEC(Ne,Te)
        do k=2,PEC_size
          if (rho .lt. PEC_dens(k)) then
            if (abs(PEC_dens(k)-rho) .lt. abs(PEC_dens(k-1)-rho)) then
              PEC_index_Ne = k
            else
              PEC_index_Ne = k-1            
            endif
            exit
          endif
          !if (k .eq. PEC_size) write(*,'(A,3e18.6)') 'Warning! no PEC found for density :',rho,PEC_dens(1),PEC_dens(k)
        enddo
        do k=2,PEC_size
          if (Te .lt. PEC_temp(k)) then
            if (abs(PEC_temp(k)-Te) .lt. abs(PEC_temp(k-1)-Te)) then
              PEC_index_Te = k
            else
              PEC_index_Te = k-1            
            endif
            exit
          endif
          !if (k .eq. PEC_size) write(*,'(A,3e18.6)') 'Warning! no PEC found for temperature:',Te,PEC_temp(1),PEC_temp(k)
        enddo        
        PEC_index = (PEC_index_Ne-1)*PEC_size + PEC_index_Te

        ! --- Integrate Emissivity
        Light(i) = Light(i) + step*rho_n*rho*PEC(PEC_index)
        !if (Light(i) .lt. 1.d-14) write(*,*)'zero light:',rho_n,rho,PEC(PEC_index)
                
      endif
      
    enddo
         
    !write(*,'(A,i6,e)')'Light : ',i,Light(i)
    !write(*,'(A,i6,A,i6,A,e)')'n_pix : ',i,' out of ',n_pix,' has light ',Light(i)
  
  enddo
  
  ! --- MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
  ! --- MPI collect
  if (my_id .eq. 0) then
    do j=1,n_cpu-1
      call mpi_recv(pix_start,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      call mpi_recv(pix_end,  1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      call mpi_recv(pix_delta,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      nrecv = pix_delta
      call mpi_recv(Light(pix_start:pix_end),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
    enddo
  else
    call mpi_send(pix_start, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    call mpi_send(pix_end,   1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    call mpi_send(pix_delta, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    nsend = pix_delta
    call mpi_send(Light(pix_start:pix_end), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
  endif

  ! --- write matlab plot file (temporary!)
  !if (my_id .eq. 0) then
  if (.false.) then
    ! --- Open output files
    open(20, file='plot_fast_camera.py', ACTION = 'write')
    write(20,'(A)')      '#!/usr/bin/env python'
    write(20,'(A)')      'import numpy'
    write(20,'(A)')      'import pylab'
    write(20,'(A)')      'import matplotlib.pyplot as plt'
    write(20,'(A)')      'import matplotlib'
    write(20,'(A)')      'def main():'
    write(20,'(A,i6,A)') ' xx = numpy.zeros(',n_pix,')'
    write(20,'(A,i6,A)') ' yy = numpy.zeros(',n_pix,')'
    write(20,'(A,i6,A)') ' ii = numpy.zeros(',n_pix,')'
        
    maxlight = maxval(Light)
    do i=1, n_pix
      write(20,'(A,i6,A,e18.5)') ' xx[',i-1,'] = ',Zp(i)*1.d5
      write(20,'(A,i6,A,e18.5)') ' yy[',i-1,'] = ',Yp(i)*1.d5
      write(20,'(A,i6,A,e18.5)') ' ii[',i-1,'] = ',Light(i)/maxlight
    enddo
    write(20,'(A)')      ' pylab.figure()'
    write(20,'(A,i6,A)') ' for i in range(',n_pix,'):'
    write(20,'(A)')      '  my_col = str(ii[i]) # **0.4'
    write(20,'(A)')      "  pylab.plot(xx[i],yy[i],'.',c=my_col)"
    write(20,'(A)')      " pylab.axis('equal')"
    write(20,'(A)')      ' pylab.show()'
    write(20,'(A)')      'main()'

    close(20)
  endif

  
  ! --- write image file
  if (my_id .eq. 0) then
    maxlight = maxval(Light)
    ! --- Open image file with PPM P3 format
    open(unit=2,file='fast_camera.ppm',status='unknown')
    write(*,*) 'Now writing PPM (P3) file : ', 'fast_camera.ppm'
    ! --- header
    write(2,'(A)') 'P3'
    write(2,'(2(1x,i4),'' 255 '')')  n_pixels_hor, n_pixels_ver
    ! --- Write to image file
    icnt = 0
    do j=1, n_pixels_ver
      do i=1, n_pixels_hor
        call get_color(0, Light(j + (i-1)*n_pixels_ver)/maxlight, rgb)
        do k = 1, 3
          itmp = ichar(rgb(k))
          icnt = icnt + 4
          if (icnt .LT. 60) then
            write(2,fmt='(1x,i3,$)') itmp     ! "$" is not standard.
          else
            write(2,fmt='(1x,i3)') itmp
            icnt = 0
          endif
        enddo
      enddo
    enddo
    write(2,'(A)') ' '
    close(2)
  endif

  

  ! --- MPI finilise
  call MPI_FINALIZE(ierr)

  if (my_id .eq. 0) write(*,*)'finished...'
 
end program jorek2_fast_camera






! --- Given density, return red/green/blue coefs for colormap
subroutine get_color(colormap, density, rgb)

  implicit none
  
  ! --- Routine parameters
  integer, intent(in)      :: colormap
  real*8,  intent(in)      :: density
  character*1, intent(out) :: rgb(3)
  
  ! --- Internal parameters
  real*8  ::  red,  gre,  blu
  integer :: ired, igre, iblu
  

  ! --- Heat colorbar from white to yellow to red to black
  if (colormap .eq. 1) then
    if (density .lt. 1.d0/3.d0) then
      red = 1.d0     
      gre = 1.d0
      blu = 1.d0 - density * 3.d0
    elseif (density .lt. 2.d0/3.d0) then
      red = 1.d0
      gre = 1.d0 - (density-1.d0/3.d0) * 3.d0
      blu = 0.d0
    else
      red = 1.d0 - (density-2.d0/3.d0) * 3.d0
      gre = 0.d0
      blu = 0.d0
    endif
  elseif (colormap .eq. 2) then
  ! --- Rainbow colorbar from blue to green to white to yellow to red to black
    if (density .lt. 1.d0/5.d0) then
      red = 0.d0   
      gre = density * 5.d0
      blu = 1.d0 - density * 5.d0
    elseif (density .lt. 2.d0/5.d0) then
      red = (density-1.d0/5.d0) * 5.d0
      gre = 1.d0  
      blu = (density-1.d0/5.d0) * 5.d0
    elseif (density .lt. 3.d0/5.d0) then
      red = 1.d0
      gre = 1.d0  
      blu = 1.d0 - (density-2.d0/5.d0) * 5.d0
    elseif (density .lt. 4.d0/5.d0) then
      red = 1.d0
      gre = 1.d0 - (density-3.d0/5.d0) * 5.d0
      blu = 0.d0
    else
      red = 1.d0 - (density-4.d0/5.d0) * 5.d0
      gre = 0.d0
      blu = 0.d0
    endif
  elseif (colormap .eq. 0) then
  ! --- From black to white
    red = 0.d0 + density
    gre = 0.d0 + density
    blu = 0.d0 + density
  else
  ! --- From white to black
    red = 1.d0 - density
    gre = 1.d0 - density
    blu = 1.d0 - density
  endif
  
  ! --- Convert to 255 bit map  
  ired = int(red * 255.0D+00)
  igre = int(gre * 255.0D+00)
  iblu = int(blu * 255.0D+00)
  if (ired .GT. 255) ired = 255
  if (igre .GT. 255) igre = 255
  if (iblu .GT. 255) iblu = 255
  if (ired .LT.   0) ired =   0
  if (igre .LT.   0) igre =   0
  if (iblu .LT.   0) iblu =   0
  rgb(1) = char(ired)
  rgb(2) = char(igre)
  rgb(3) = char(iblu)
      

  return
end subroutine get_color

