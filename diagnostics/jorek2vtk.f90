!> Program to convert a JOREK2 restart file into binary VTK format
program jorek2vtk

use parameters, only: n_var, variable_names
use data_structure
use phys_module
use basis_at_gaussian
use diffusivities, only: get_dperp, get_zkperp

implicit none

type (type_node_list)   , pointer :: node_list
type (type_element_list), pointer :: element_list
type (type_surface_list)          :: flux_list

integer               :: nnoel, nnos, nel, nsub, inode, ielm, n_scalars, n_vectors
real*4,allocatable    :: xyz (:,:), scalars(:,:), vectors(:,:,:)
integer,allocatable   :: ien (:,:)
integer, parameter    :: ivtk = 22 ! an arbitrary unit number for the VTK output file
integer               :: i, j, k, m, etype, irst, int, i_var, i_tor, i_tor_old, i_plane, index, index_node, my_id
character             :: buffer*80, lf*1, str1*12, str2*12
character*12, allocatable :: scalar_names(:), vector_names(:)
real*8                :: s, t
real*8                :: P,P_s,P_t,P_st,P_ss,P_tt
real*8                :: R,R_s,R_t,R_st,R_ss,R_tt
real*8                :: Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
real*8                :: Ps0, Ps0_s, Ps0_t, Ps0_st, Ps0_ss, Ps0_tt, Psi, Ps_s, Ps_t, Ps_st, Ps_ss, Ps_tt
real*8                :: ZJ0, ZJ0_s, ZJ0_t, ZJ0_st, ZJ0_ss, ZJ0_tt, ZJ,  ZJ_s, ZJ_t, ZJ_st, ZJ_ss, ZJ_tt
real*8                :: U0,  U0_s,  U0_t,  U0_st,  U0_ss,  U0_tt,  U,   U_s,  U_t,  U_st,  U_ss,  U_tt
real*8                :: W0,  W0_s,  W0_t,  W0_st,  W0_ss,  W0_tt,  W,   W_s,  W_t,  W_st,  W_ss,  W_tt
real*8                :: ZN0, ZN0_s, ZN0_t, ZN0_st, ZN0_ss, ZN0_tt, RHO, RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt
real*8                :: T0,  T0_s,  T0_t,  T0_st,  T0_ss,  T0_tt,  TT,  TT_s, TT_t, TT_st, TT_ss, TT_tt
real*8                :: Ti0, Ti0_s, Ti0_t, Ti0_st, Ti0_ss, Ti0_tt, Ti,  Ti_s, Ti_t, Ti_st, Ti_ss, Ti_tt
real*8                :: Te0, Te0_s, Te0_t, Te0_st, Te0_ss, Te0_tt, Te, Te_s,  Te_t, Te_st, Te_ss, Te_tt
real*8                :: V0,  V0_s,  V0_t,  V0_st,  V0_ss,  V0_tt,  V, V_s, V_t, V_st, V_ss, V_tt
real*8                :: ps0_x, ps0_y, psi_sum, ps_x, ps_y, ps_p
real*8                :: u0_x,  u0_y,  u_sum,   u_x,  u_y,  u_p
real*8                :: zj0_x, zj0_y, zj_sum,  zj_x, zj_y, zj_p
real*8                :: w0_x,  w0_y,  w_sum,   w0_xx, w0_yy, w_x, w_y, w_p, w_xx, w_yy
real*8                :: zn0_x, zn0_y, zn_sum,  zn_x, zn_y, zn_p
real*8                :: T0_x,  T0_y,  T_sum,   TT_x, TT_y, TT_p
real*8                :: Ti0_x, Ti0_y, Ti_sum,  Ti_x, Ti_y, Ti_p
real*8                :: Te0_x, Te0_y, Te_sum,  Te_x, Te_y, Te_p
real*8                :: AR_Z, AR_p, AZ_R, AZ_p, A3_R, A3_Z, Fprof
real*8                :: psi_axis,      R_axis,      Z_axis,      s_axis,      t_axis
real*8                :: psi_xpoint(2), R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2)
real*8                :: psi_norm, psi_bnd, grad_psi
real*8                :: xjac, xjac_x, xjac_y, v_perp, Psi_J, R_p, error, Btot, BigR
real*8                :: particle_source, D_prof, ZK_prof, source_pellet, ZKpar_T
integer 	      :: i_find, i_elm_find(8)
integer               :: n_fluxes, n_neo, n_bfield, n_vfield, n_pellet
real*8  	      :: Router,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8  	      :: Zouter,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8  	      :: s_find(8), t_find(8)
real*8                :: Jb, minRad,rho_norm,t_norm
integer               :: i_elm_axis, i_elm_xpoint(2), k_tor, ifail, ierr
logical               :: without_n0_mode, SI_units
logical               :: include_fluxes, include_neo, include_magnetic_field, include_velocity_field
real*8                :: toroidal_angle
!====================== --- add the diagnostics Er, Vtheta and Vneo
real*8                :: Er, psi_abs, Vtheta, Btheta, Mach_par,Mach_pol,Vsound, Vneo
real*8                :: amu_neo_node, aki_neo_node
real*8                :: Vperp_e 


namelist /vtk_params/ nsub, i_tor, i_plane, without_n0_mode, SI_units, &
                      include_fluxes, include_neo, include_magnetic_field, include_velocity_field

write(*,*) '***************************************'
write(*,*) '*       jorek2vtk                     *'
write(*,*) '***************************************'
write(*,*) ' if your VTK is smaller than expected,'
write(*,*) ' please consider the new parameters:'
write(*,*) '   include_fluxes, include_neo'
write(*,*) '***************************************'

call flush_it(6)

allocate(node_list)
allocate(element_list)

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Preset parameters
nsub            = 2       ! Number of subdivisions of the cubic finite elements into linear pieces
i_tor           = -1      ! If i_tor > 0, only this mode will be included in the vtk file...
i_plane         = 1       ! ... otherwise, all modes will be summed up at the toroidal plane i_plane
without_n0_mode = .false. ! If true, do not include the n=0 mode (i_tor=1)
SI_units        = .false. ! when true, write variables in SI units
include_fluxes  = .false. ! include energy and density fluxes (or not)
include_neo     = .false. ! include neoclassical and more terms (or not)
include_magnetic_field = .false. ! include vector of magnetic field (or not)
include_velocity_field = .false. ! include vector of velocity field (or not)

! --- Read parameters from namelist file 'vtk.nml' if it exists
open(42, file='vtk.nml', action='read', status='old', iostat=ierr)
if ( ierr == 0 ) then
  write(*,*) 'Reading parameters from vtk.nml namelist.'
  read(42,vtk_params)
  close(42)
end if

write(*,*)
write(*,*) 'Parameters:'
write(*,*) '-----------'
write(*,*) 'nsub            =', nsub
write(*,*) 'i_tor           =', i_tor
write(*,*) 'i_plane         =', i_plane
write(*,*) 'without_n0_mode =', without_n0_mode
write(*,*) 'si_units        =', si_units
write(*,*) 'include_fluxes  =', include_fluxes
write(*,*) 'include_neo     =', include_neo
write(*,*) 'include_magnetic_field =',include_magnetic_field
write(*,*) 'include_velocity_field =',include_velocity_field

write(*,*) '-----------'
write(*,*) 'n_tor           =', n_tor
write(*,*) 'n_period        =', n_period
write(*,*) 'F0              =', F0
write(*,*) 'R_geo,Z_geo     =', R_geo, Z_geo
write(*,*)
call flush_it(6)

! --- Number of scalars to write to the VTK output file
#ifdef fullmhd
n_scalars = n_var + 3
! --- Number of vectors to write to the VTK output file
n_vectors = 0
#else
n_scalars = n_var
n_vectors = 0
n_fluxes  = 0
n_neo     = 0
n_bfield  = 0
n_vfield  = 0
n_pellet  = 0

if (include_fluxes) then
  n_fluxes = 8
  n_scalars = n_scalars + n_fluxes
endif
if (include_neo) then
  n_neo = 12
  n_scalars = n_scalars + n_neo
endif
if (include_magnetic_field) then
  n_bfield  = 1
  n_vectors = n_vectors + n_bfield
endif
if (include_velocity_field) then
  n_vfield  = 1
  n_vectors = n_vectors + n_vfield
endif
if (use_pellet) then
  n_pellet  = 1 
  n_scalars = n_scalars+ n_pellet
endif
#endif

allocate(scalar_names(n_scalars), vector_names(n_vectors))

grad_psi = 0.d0

scalar_names(1:n_var) = variable_names(1:n_var)
if ( SI_units ) then
   scalar_names(3)='j_MA/m2     '
   scalar_names(5)='n_e20m-3    '
   if (jorek_model .eq. 400) then
      scalar_names(6)='Ti_keV      '
      scalar_names(8)='Te_keV      '
   else
      scalar_names(6)='Te+Ti_keV   '
   endif
   scalar_names(7)='Vpar_km/s   '
endif

#ifdef fullmhd
scalar_names(n_var+1:n_var+3) = (/  'B_phi       ', 'B_R         ', 'B_Z         '/)

#else
if ( SI_units ) then

  if (include_fluxes) then
    scalar_names(n_var+1:n_var+n_fluxes) = (/ &
     'P_kPa       ', 'E_flux_Kpar ', 'E_flux_kperp', 'E_flux_Vpar ', &
     'E_flux_Vperp', 'D_flux_Dperp', 'D_flux_Vpar ', 'D_flux_Vperp'/)
  endif     
  if (include_neo) then
    scalar_names(n_var+1+n_fluxes:n_var+n_fluxes+n_neo) = (/ &          
     'Er_kV/m     ', 'Vtheta_km/s ', 'Mach_par    ', 'Mach_pol    ', &
     'Vsound_km/s ', 'Btot_T      ', 'J-bootstrap ', 'Vneo_km/s   ', 'psi_norm    ', 'Vperp_e_km/s', &
     'ki_neo      ', 'mu_neo      '/)
  endif

else

  if (include_fluxes) then  
    scalar_names(n_var+1:n_var+n_fluxes) = (/ &
      'pressure    ', 'E_flux_Kpar ', 'E_flux_kperp', 'E_flux_Vpar ',&
      'E_flux_Vperp', 'D_flux_Dperp', 'D_flux_Vpar ', 'D_flux_Vperp'/)
  endif    
  if (include_neo) then
    scalar_names(n_var+1+n_fluxes:n_var+n_fluxes+n_neo) = (/ &              
      'Er          ', 'Vtheta      ', 'Mach_par    ', 'Mach_pol    ', &
      'Vsound      ', 'Btot        ', 'J-bootstrap ', 'Vneo        ', 'psi_norm    ', 'Vperp_e     ', &
      'ki_neo      ', 'mu_neo      '/)
   endif
   
endif
#endif

if (include_magnetic_field)  vector_names(1:n_bfield)                   = (/ 'B_R  ','B_Z   ','B_phi   '/)
if (include_velocity_field)  vector_names(n_bfield+1:n_bfield+n_vfield) = (/ 'V_R  ','V_Z   ','V_phi   '/)

do k_tor=1, n_tor
  mode(k_tor) = + int(k_tor / 2) * n_period
enddo

call import_restart(node_list,element_list, 'jorek_restart.rst', rst_format, ierr)

call initialise_basis                              ! define the basis functions at the Gaussian points

nnos = nsub*nsub*node_list%n_nodes
allocate(xyz(3,nnos),scalars(nnos,1:n_scalars),vectors(nnos,3,1:n_vectors))

nnoel = 4
nel   = (nsub-1)*(nsub-1)*element_list%n_elements
allocate(ien(nnoel,nel))

inode   = 0
ielm    = 0
scalars = 0.d0
vectors = 0.d0
xyz     = 0
ien     = 0

call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

if (xpoint) then
  call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
  psi_bnd  = psi_xpoint(1)
  if( (xcase .eq. 2) .or. ((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
    psi_bnd = psi_xpoint(2)
  endif
else
  psi_bnd = 0.d0
endif

if (jorek_model .eq. 000) then      
  flux_list%n_psi = 1
  call tr_allocate(flux_list%psi_values,1,flux_list%n_psi,"flux_list%psi_values",CAT_GRID)
  flux_list%psi_values(1) = psi_bnd
  call find_flux_surfaces(xpoint,xcase,node_list,element_list,flux_list)
  call find_theta_surface(node_list, element_list, flux_list, 1, 0.0, R_axis, Z_axis,i_elm_find,s_find,t_find,i_find)
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
		 Router,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,  &
		 Zouter,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
  call tr_deallocate(flux_list%psi_values,"flux_list%psi_values",CAT_GRID)
  minRad = Router - R_axis
  write(*,*)'minor radius : ',minRad
else
  minRad = 0.0
endif
  
! --- You may choose to print your poloidal snapshot at a different toroidal angle
toroidal_angle = 0.d0 ! 2*PI / 6
if (toroidal_angle .ne. 0.d0) then
  do i_tor=1, n_tor
    mode(i_tor) = + int(i_tor / 2) * n_period
  enddo
  HZ(1,i_plane)   = 1.d0
  do i=1,(n_tor-1)/2
    HZ(2*i,i_plane)      = cos(mode(2*i)   * toroidal_angle )
    HZ(2*i+1,i_plane)    = sin(mode(2*i+1) * toroidal_angle )
  enddo
endif

do i=1,element_list%n_elements

   ! if(element_list%element(i)%n_sons.eq.0) then

  do j=1,nsub

    s = float(j-1)/float(nsub-1)

    do k=1,nsub

      t = float(k-1)/float(nsub-1)

      call interp_RZ(node_list,element_list,i,s,t,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)

      xjac  = R_s * Z_t - R_t * Z_s

      if ( xjac == 0.d0 ) xjac = 1.d-8

      BigR  = R

      xjac_x  = (R_ss*Z_t**2 - Z_ss*R_t*Z_t - 2.d0*R_st*Z_s*Z_t   & 
              + Z_st*(R_s*Z_t + R_t*Z_s) + R_tt*Z_s**2 - Z_tt*R_s*Z_s) / xjac

      xjac_y  = (Z_tt*R_s**2 - R_tt*Z_s*R_s - 2.d0*Z_st*R_t*R_s   &
              + R_st*(Z_t*Z_s + Z_s*R_t) + Z_ss*R_t**2 - R_ss*Z_t*R_t) / xjac

      inode = inode+1

      xyz(1:3,inode) = (/ R, Z, 0.d0/)

      !====================== --- specific for axisymmetric quantities
      ! Put here all quantities that are axisymmetric (n=0 mode only) and should not be summed 
      ! over all harmonics: for instance, to compute Vtheta, Er, Vneo, etc.
      ! ===> this corresponds to forcing i_tor = 1 (thus n=0 only)

      ! save old values
      i_tor_old = i_tor
      i_tor     = 1
      ! compute all derivatives, as in loop below
      if ((xjac .gt. 1.d-6)) then

        call interp(node_list,element_list,i,1,i_tor,s,t,Ps0,Ps0_s,Ps0_t,Ps0_st,Ps0_ss,Ps0_tt)
        call interp(node_list,element_list,i,2,i_tor,s,t,U0, U0_s, U0_t, U0_st, U0_ss, U0_tt)
        call interp(node_list,element_list,i,3,i_tor,s,t,ZJ0,ZJ0_s,ZJ0_t,ZJ0_st,ZJ0_ss,ZJ0_tt)
        call interp(node_list,element_list,i,4,i_tor,s,t,W0, W0_s, W0_t, W0_st, W0_ss, W0_tt)
        call interp(node_list,element_list,i,5,i_tor,s,t,ZN0,ZN0_s,ZN0_t,ZN0_st,ZN0_ss,ZN0_tt)
        call interp(node_list,element_list,i,6,i_tor,s,t,T0, T0_s, T0_t, T0_st, T0_ss, T0_tt)

        if ( jorek_model >= 300 ) then
          call interp(node_list,element_list,i,7,i_tor,s,t,V0,V0_s,V0_t,V0_st,V0_ss,V0_tt)
        else
          V0=0; V0_s=0; V0_t=0; V0_st=0; V0_ss=0; V0_tt=0
        end if

        u0_x   = (   Z_t * U0_s  - Z_s * U0_t ) / xjac
        u0_y   = ( - R_t * U0_s  + R_s * U0_t ) / xjac

        ps0_x  = (   Z_t * Ps0_s - Z_s * Ps0_t ) / xjac
        ps0_y  = ( - R_t * Ps0_s + R_s * Ps0_t ) / xjac

        T0_x   = (   Z_t * T0_s  - Z_s * T0_t ) / xjac
        T0_y   = ( - R_t * T0_s  + R_s * T0_t ) / xjac

        zj0_x  = (   Z_t * ZJ0_s - Z_s * ZJ0_t ) / xjac
        zj0_y  = ( - R_t * ZJ0_s + R_s * ZJ0_t ) / xjac

        zn0_x  = (   Z_t * zn0_s - Z_s * zn0_t ) / xjac
        zn0_y  = ( - R_t * zn0_s + R_s * zn0_t ) / xjac

        if (include_neo) then

            !*** compute diagnostics ***
          psi_abs = sqrt(ps0_x*ps0_x + ps0_y * ps0_y)
          Btheta  = (psi_abs/R)
          Vtheta  = 0.d0
          Vperp_e = 0.0
          Vneo    = 0.d0
          Er      = 0.d0
	  mach_par= 0.d0
          mach_pol= 0.d0
          vsound  = 0.d0

          if ((psi_abs .gt. 1.d-6) .and. (ZN0.gt.1.d-6) .and. (abs(Btheta).gt.1.d-6)) then

            Vtheta  = -1./Btheta*((u0_x + tauIC/ZN0*(T0_x*ZN0 + ZN0_x*T0))*ps0_x  + &
                                  (u0_y + tauIC/ZN0*(T0_y*ZN0 + ZN0_y*T0))*ps0_y) + V0*Btheta

            Vperp_e = -1./Btheta*((u0_x - tauIC/ZN0*(T0_x*ZN0 + ZN0_x*T0))*ps0_x  + &
                                  (u0_y - tauIC/ZN0*(T0_y*ZN0 + ZN0_y*T0))*ps0_y)

            if (NEO) then
!                num_neo_file= ( neo_file /= 'none')
              if (num_neo_file) then
!                  write(*,*) 'neo_file=',neo_file
!                  write(*,*) 'using ki and mui profiles from file "'//trim(neo_file)//'"'
                call neo_coef( xpoint, xcase, Z, Z_xpoint, Ps0 ,psi_axis, psi_bnd, &
                               amu_neo_node, aki_neo_node)
                Vneo   = aki_neo_node / Btheta*tauIC  * (ps0_x*T0_x + ps0_y*T0_y)
              else
                Vneo   = aki_neo_const / Btheta*tauIC * (ps0_x*T0_x + ps0_y*T0_y)
              endif

            endif  ! NEO

            Er       = -(U0_x * ps0_x + U0_y * ps0_y)/psi_abs         ! radial electric field
            Btot     = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR       ! total magnetic field (equilibrium)
            Vsound   = sqrt(GAMMA*T0)/Btot                            ! sound speed
            Mach_par = V0/Vsound                                      ! parallel Mach number
            Mach_pol = Vtheta/Vsound                                  ! poloidal Mach number

          endif !psi_abs
 
          ! save those specific values of axisymmetric parameters
          if (grad_psi .ne. 0.d0) then
            scalars(inode,n_var+n_fluxes+1) = Er
            scalars(inode,n_var+n_fluxes+2) = Vtheta
            scalars(inode,n_var+n_fluxes+3) = Mach_par
            scalars(inode,n_var+n_fluxes+4) = Mach_pol
            scalars(inode,n_var+n_fluxes+5) = Vsound
            scalars(inode,n_var+n_fluxes+6) = Btot
            scalars(inode,n_var+n_fluxes+8) = Vneo    ! number n_var+n_fluxes+7 is jbootstrap, see below.
            scalars(inode,n_var+n_fluxes+10)= Vperp_e ! number n_var+n_fluxes+9 is psi_norm, see below.
            if (NEO) then
               if (num_neo_file) then
                  scalars(inode,n_var+n_fluxes+11) = aki_neo_node
                  scalars(inode,n_var+n_fluxes+12) = amu_neo_node
               else
                  scalars(inode,n_var+n_fluxes+11) = aki_neo_const
                  scalars(inode,n_var+n_fluxes+12) = amu_neo_const
               endif
            endif   ! NEO

          endif     ! grad_psi
	
	endif       ! include_neo

      endif         ! xjac

#ifdef fullmhd
      ! Magnetic field components
      call interp(node_list,element_list,i,var_AR,i_tor,s,t,U0,U0_s,U0_t,U0_st,U0_ss,U0_tt)
      call interp(node_list,element_list,i,var_AZ,i_tor,s,t,V0,V0_s,V0_t,V0_st,V0_ss,V0_tt)
      call interp(node_list,element_list,i,var_A3,i_tor,s,t,W0,W0_s,W0_t,W0_st,W0_ss,W0_tt)

       AR_Z = ( - R_t * U0_s + R_s * U0_t ) / xjac
       AZ_R = (   Z_t * V0_s - Z_s * V0_t ) / xjac
       A3_R = (   Z_t * W0_s - Z_s * W0_t ) / xjac
       A3_Z = ( - R_t * W0_s + R_s * W0_t ) / xjac
       AR_p = 0.d0 ; AZ_p = 0.d0

       call interp(node_list,element_list,i,456,i_tor,s,t,W0,W0_s,W0_t,W0_st,W0_ss,W0_tt)
       Fprof = W

       scalars(inode,n_var+1) = ( AZ_R - AR_Z )+ Fprof / R        ! B_phi
       scalars(inode,n_var+2) = ( A3_Z - AZ_p )/ BigR             ! B_R
       scalars(inode,n_var+3) = ( AR_p - A3_R )/ BigR             ! B_Z
#endif

       ! old values back to normal
       i_tor = i_tor_old

       !====================== --- specific for NON-axisymmetric quantities
       ! 2 cases, depending on the value of i_tor chosen
       if ((i_tor .ge. 1) .and. (i_tor .le. n_tor)) then

         do m=1,n_var
           call interp(node_list,element_list,i,m,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
           scalars(inode,m) = P
         enddo

         if ((xjac .gt. 1.d-6)) then

           call interp(node_list,element_list,i,1,i_tor,s,t,Psi,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt)
           call interp(node_list,element_list,i,2,i_tor,s,t,U,U_s,U_t,U_st,U_ss,U_tt)
           call interp(node_list,element_list,i,3,i_tor,s,t,ZJ,ZJ_s,ZJ_t,ZJ_st,ZJ_ss,ZJ_tt)
           call interp(node_list,element_list,i,4,i_tor,s,t,W,W_s,W_t,W_st,W_ss,W_tt)
           call interp(node_list,element_list,i,5,i_tor,s,t,RHO,RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt)
           call interp(node_list,element_list,i,6,i_tor,s,t,TT,TT_s,TT_t,TT_st,TT_ss,TT_tt)
           if ( jorek_model >= 300 ) then
             call interp(node_list,element_list,i,7,i_tor,s,t,V,V_s,V_t,V_st,V_ss,V_tt)
           else
             V=0; V_s=0; V_t=0; V_st=0; V_ss=0; V_tt=0
           end if

           u_x   = (   Z_t * U_s  - Z_s * U_t ) / xjac
           u_y   = ( - R_t * U_s  + R_s * U_t ) / xjac

           ps_x  = (   Z_t * PS_s - Z_s * PS_t ) / xjac
           ps_y  = ( - R_t * PS_s + R_s * PS_t ) / xjac

           TT_x  = (   Z_t * TT_s - Z_s * TT_t ) / xjac
           TT_y  = ( - R_t * TT_s + R_s * TT_t ) / xjac

           zj_x  = (   Z_t * ZJ_s - Z_s * ZJ_t ) / xjac
           zj_y  = ( - R_t * ZJ_s + R_s * ZJ_t ) / xjac

            !*** compute diagnostics ***
           v_perp  = R * sqrt(u_x*u_x + u_y*u_y)

               !	  vectors(inode,:,1) = (/ - R * u0_y ,   + R * u0_x ,   0.d0 /)
               !          vectors(inode,:,2) = (/ + ps_y /R * V, - ps_x /R * V, F0/R * V /)
               !          vectors(inode,:,3) = (/ - R * u0_y + ps_y /R * V, + R * u0_x - ps_x /R * V, F0/R * V /)

           psi_J = (Ps_s * ZJ_t - PS_t * ZJ_s ) / xjac
           R_p   = (2.d0 * R * (R_s * (RHO_t * TT + RHO * TT_t) - R_t * (RHO_s * TT + RHO * TT_s) )) / xjac
           error = psi_J - R_p  ! "error" in Grad_Shafranov equilibrium force balance

         endif  ! xjac check

#ifdef fullmhd
         ! Magnetic field components
         call interp(node_list,element_list,i,var_AR,i_tor,s,t,U,U_s,U_t,U_st,U_ss,U_tt)
         call interp(node_list,element_list,i,var_AZ,i_tor,s,t,V,V_s,V_t,V_st,V_ss,V_tt)
         call interp(node_list,element_list,i,var_A3,i_tor,s,t,W,W_s,W_t,W_st,W_ss,W_tt)

         AR_Z = ( - R_t * U_s + R_s * U_t ) / xjac
         AZ_R = (   Z_t * V_s - Z_s * V_t ) / xjac
         A3_R = (   Z_t * W_s - Z_s * W_t ) / xjac
         A3_Z = ( - R_t * W_s + R_s * W_t ) / xjac

         call interp(node_list,element_list,i,var_AR,i_tor+1,s,t,U,U_s,U_t,U_st,U_ss,U_tt) ! sine
         call interp(node_list,element_list,i,var_AZ,i_tor+1,s,t,V,V_s,V_t,V_st,V_ss,V_tt)
         AR_p = U  * HZ_p(i_tor,i_plane)  
         AZ_p = V  * HZ_p(i_tor,i_plane) 

         if (i_tor == 1) then
           call interp(node_list,element_list,i,456,i_tor,s,t,Fprof,W_s,W_t,W_st,W_ss,W_tt)
           scalars(inode,n_var+1)   = ( AZ_R - AR_Z )  + Fprof / R  ! B_phi
         else
           scalars(inode,n_var+1)   = ( AZ_R - AR_Z )    
         endif
         scalars(inode,n_var+2) = ( A3_Z - AZ_p )/ BigR  ! B_R
         scalars(inode,n_var+3) = ( AR_p - A3_R )/ BigR  ! B_Z
#endif

         else  ! i_tor

           u_sum   = 0.d0; u_x  = 0.d0; u_y  = 0.d0; u_p  = 0.d0
           psi_sum = 0.d0; ps_x = 0.d0; ps_y = 0.d0; ps_p = 0.d0
           zj_sum  = 0.d0; zj_x = 0.d0; zj_y = 0.d0; zj_p = 0.d0
           T_sum   = 0.d0; TT_x = 0.d0; TT_y = 0.d0; TT_p = 0.d0
           zn_sum  = 0.d0; zn_x = 0.d0; zn_y = 0.d0; zn_p = 0.d0;
           Ti_sum  = 0.d0; Ti_x = 0.d0; Ti_y = 0.d0; Ti_p = 0.d0
           Te_sum  = 0.d0; Te_x = 0.d0; Te_y = 0.d0; Te_p = 0.d0
           w_sum   = 0.d0; w_x  = 0.d0; w_y  = 0.d0; w_p  = 0.d0; w_xx = 0.d0; w_yy = 0.d0

           do i_tor = 1, n_tor

             if ( ( i_tor == 1 ) .and. ( without_n0_mode ) ) cycle ! Do not include the n=0 mode

               do m=1,n_var
                  call interp(node_list,element_list,i,m,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
                  scalars(inode,m) = scalars(inode,m) + P * HZ(i_tor,i_plane)
               enddo

               call interp(node_list,element_list,i,1,i_tor,s,t,Psi,Ps_s, Ps_t, Ps_st, Ps_ss, Ps_tt)
               call interp(node_list,element_list,i,2,i_tor,s,t,U  ,U_s,  U_t,  U_st,  U_ss,  U_tt)
               call interp(node_list,element_list,i,3,i_tor,s,t,ZJ ,ZJ_s, ZJ_t, ZJ_st, ZJ_ss, ZJ_tt)
               call interp(node_list,element_list,i,4,i_tor,s,t,W  ,W_s,  W_t,  W_st,  W_ss,  W_tt)
               call interp(node_list,element_list,i,5,i_tor,s,t,RHO,RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt)
               call interp(node_list,element_list,i,6,i_tor,s,t,TT ,TT_s, TT_t, TT_st, TT_ss, TT_tt)
               if ( jorek_model >= 300 ) then
                  call interp(node_list,element_list,i,7,i_tor,s,t,V,V_s,V_t,V_st,V_ss,V_tt)
               else
                  V=0; V_s=0; V_t=0; V_st=0; V_ss=0; V_tt=0
               end if
               if ( jorek_model .eq. 400 ) then
                  call interp(node_list,element_list,i,6,i_tor,s,t,Ti,Ti_s,Ti_t,Ti_st,Ti_ss,Ti_tt)
                  call interp(node_list,element_list,i,8,i_tor,s,t,Te,Te_s,Te_t,Te_st,Te_ss,Te_tt)
               end if

               psi_sum = psi_sum + psi * HZ(i_tor,i_plane)
               zj_sum  = zj_sum  + zj  * HZ(i_tor,i_plane)
               u_sum   = u_sum   + U   * HZ(i_tor,i_plane)
               w_sum   = w_sum   + w   * HZ(i_tor,i_plane)
               zn_sum  = zn_sum  + RHO * HZ(i_tor,i_plane)
               if ( jorek_model .eq. 400 ) then
                 Ti_sum  = Ti_sum + Ti * HZ(i_tor,i_plane)
                 Te_sum  = Te_sum + Te * HZ(i_tor,i_plane)
               end if

#ifdef fullmhd
               ! Magnetic field components
               call interp(node_list,element_list,i,var_AR,i_tor,s,t,U,U_s,U_t,U_st,U_ss,U_tt)
               call interp(node_list,element_list,i,var_AZ,i_tor,s,t,V,V_s,V_t,V_st,V_ss,V_tt)
               call interp(node_list,element_list,i,var_A3,i_tor,s,t,W,W_s,W_t,W_st,W_ss,W_tt)

               AR_Z = ( - R_t * U_s + R_s * U_t ) / xjac
               AZ_R = (   Z_t * V_s - Z_s * V_t ) / xjac
               A3_R = (   Z_t * W_s - Z_s * W_t ) / xjac
               A3_Z = ( - R_t * W_s + R_s * W_t ) / xjac

               call interp(node_list,element_list,i,var_AR,i_tor,s,t,U,U_s,U_t,U_st,U_ss,U_tt)
               call interp(node_list,element_list,i,var_AZ,i_tor,s,t,V,V_s,V_t,V_st,V_ss,V_tt)
               AR_p = U  * HZ_p(i_tor,i_plane)  
               AZ_p = V  * HZ_p(i_tor,i_plane) 

               if (i_tor == 1) then
                 call interp(node_list,element_list,i,456,i_tor,s,t,Fprof,W_s,W_t,W_st,W_ss,W_tt)
                 scalars(inode,n_var+1) = scalars(inode,n_var+1) + ( AZ_R - AR_Z )  + Fprof / R  ! B_phi
               else
                 scalars(inode,n_var+1) = scalars(inode,n_var+1) + ( AZ_R - AR_Z )     * HZ(i_tor,i_plane)
                endif

               scalars(inode,n_var+2) = scalars(inode,n_var+2) + ( A3_Z - AZ_p )/ BigR * HZ(i_tor,i_plane)  ! B_R
               scalars(inode,n_var+3) = scalars(inode,n_var+3) + ( AR_p - A3_R )/ BigR * HZ(i_tor,i_plane)  ! B_Z
#endif 

               if ((xjac .gt. 1.d-6)) then      ! avoid the axis

                  u_x  = u_x   + (   Z_t * U_s - Z_s * U_t )     / xjac * HZ(i_tor,i_plane)
                  u_y  = u_y   + ( - R_t * U_s + R_s * U_t )     / xjac * HZ(i_tor,i_plane)

                  ps_x  = ps_x + (   Z_t * PS_s - Z_s * PS_t )   / xjac * HZ(i_tor,i_plane)
                  ps_y  = ps_y + ( - R_t * PS_s + R_s * PS_t )   / xjac * HZ(i_tor,i_plane)

                  zj_x  = zj_x + (   Z_t * ZJ_s - Z_s * ZJ_t )   / xjac * HZ(i_tor,i_plane)
                  zj_y  = zj_y + ( - R_t * ZJ_s + R_s * ZJ_t )   / xjac * HZ(i_tor,i_plane)

                  TT_x  = TT_x + (   Z_t * TT_s - Z_s * TT_t )   / xjac * HZ(i_tor,i_plane)
                  TT_y  = TT_y + ( - R_t * TT_s + R_s * TT_t )   / xjac * HZ(i_tor,i_plane)
                  TT_p  = TT_p + TT * HZ_p(i_tor,i_plane)

                  if ( jorek_model .eq. 400 ) then
                    Ti_x  = Ti_x + (   Z_t * Ti_s - Z_s * Ti_t )   / xjac * HZ(i_tor,i_plane)
                    Ti_y  = Ti_y + ( - R_t * Ti_s + R_s * Ti_t )   / xjac * HZ(i_tor,i_plane)
                    Ti_p  = Ti_p + Ti * HZ_p(i_tor,i_plane)
                    Te_x  = Te_x + (   Z_t * Te_s - Z_s * Te_t )   / xjac * HZ(i_tor,i_plane)
                    Te_y  = Te_y + ( - R_t * Te_s + R_s * Te_t )   / xjac * HZ(i_tor,i_plane)
                    Te_p  = Te_p + Te * HZ_p(i_tor,i_plane)
                  end if

                  zn_x = zn_x  + (   Z_t * RHO_s - Z_s * RHO_t ) / xjac * HZ(i_tor,i_plane)
                  zn_y = zn_y  + ( - R_t * RHO_s + R_s * RHO_t ) / xjac * HZ(i_tor,i_plane)
                  zn_p = zn_p  + RHO * HZ_p(i_tor,i_plane)

                  w_x  = w_x   + (   Z_t * U_s - Z_s * U_t )     / xjac * HZ(i_tor,i_plane)
                  w_y  = w_y   + ( - R_t * U_s + R_s * U_t )     / xjac * HZ(i_tor,i_plane)

                  w_xx = w_xx  + (w_ss * Z_t**2 - 2.d0*w_st * Z_s*Z_t + w_tt * Z_s**2       &
                       + w_s * (Z_st*Z_t - Z_tt*Z_s )                                        & 
                       + w_t * (Z_st*Z_s - Z_ss*Z_t ) )     / xjac**2                        & 
                       - xjac_x * (w_s* Z_t - w_t * Z_s)  / xjac**2

                  w_yy = w_yy  + (w_ss * R_t**2 - 2.d0*w_st * R_s*R_t + w_tt * R_s**2       &
                       + w_s * (R_st*R_t - R_tt*R_s )                                        &
                       + w_t * (R_st*R_s - R_ss*R_t ) )         / xjac**2                    &
                       - xjac_y * (- w_s * R_t + w_t * R_s )  / xjac**2

               endif ! xjac

            enddo  ! end loop toroidal harmonics

            psi_norm = (scalars(inode,1) - psi_axis)/(psi_bnd - psi_axis)
            if ((psi_norm .lt. 1.d0) .and. (xpoint) .and. (Z .lt. Z_xpoint(1)) .and. (xcase .ne. 2)) then
               psi_norm = 2.d0 - psi_norm
            endif
            if ((psi_norm .lt. 1.d0) .and. (xpoint) .and. (Z .gt. Z_xpoint(2)) .and. (xcase .ne. 1)) then
               psi_norm = 2.d0 - psi_norm
            endif

            if ( jorek_model .eq. 400 ) then
              call bootstrap_current_rhs(BigR, minRad, R_axis, psi_axis, psi_bnd,  &
                                         psi_sum, ps_x, ps_y, zn_sum,  zn_x, zn_y,     &
				         Ti_sum,  Ti_x, Ti_y, Te_sum,  Te_x, Te_y,       Jb)
            else
	      Jb = 0.d0
            endif
 
            v_perp  = R * sqrt(u_x*u_x + u_y * u_y)
            Btot    = sqrt(F0**2 + ps_x**2 + ps_y**2) / BigR
            D_prof  = get_dperp (psi_norm)
            ZK_prof = get_zkperp(psi_norm)

            ZKpar_T = ZK_par * ((max( scalars(inode,6), T_min ))/T_0)**2.5

            grad_psi = sqrt(ps_x*ps_x + ps_y*ps_y)

            !   'E_flux_Kpar ','E_flux_kperp','E_flux_Vpar ','E_flux_Vperp','D_flux_Dperp','D_flux_Vpar ','D_flux_Vperp'/)

            if (include_fluxes) then

              scalars(inode,n_var+1)   = scalars(inode,5) * scalars(inode,6)

              if (grad_psi .ne. 0.d0) then

                scalars(inode,n_var+2)  = ZKpar_T * ( F0 * TT_p / BigR**2  + (TT_x * ps_y - TT_y * ps_x) / BigR ) / Btot

                scalars(inode,n_var+3)  = ZK_prof * (TT_x * ps_x + TT_y * ps_y) / grad_psi

                scalars(inode,n_var+4)  = scalars(inode,5) * scalars(inode,6) * scalars(inode,7) * Btot

                scalars(inode,n_var+5)  = BigR   * (u_x * ps_y - u_y * ps_x) / sqrt(ps_x*ps_x + ps_y*ps_y) * scalars(inode,5) * scalars(inode,6)

                scalars(inode,n_var+6)  = D_prof * (zn_x * ps_x + zn_y * ps_y) / sqrt(ps_x*ps_x + ps_y*ps_y)

                scalars(inode,n_var+7)  = scalars(inode,5) * scalars(inode,7) * Btot
 
                scalars(inode,n_var+8)  = BigR   * (u_x * ps_y - u_y * ps_x) / sqrt(ps_x*ps_x + ps_y*ps_y) * scalars(inode,5)

	      endif ! grad_psi
	    endif   ! include_fluxes

            if (include_neo)   then
	      scalars(inode,n_var+n_fluxes+7) = Jb
              scalars(inode,n_var+n_fluxes+9) = psi_norm
            endif
		
            if (use_pellet) then

!                  call pellet_source(pellet_amplitude, pellet_R, pellet_Z, pellet_psi, pellet_phi, &
!                                     pellet_radius, pellet_delta_psi, pellet_sig, pellet_length,   &
!                                     R,Z,ps0,0.d0,source_pellet)
!                  scalars(inode,n_var+n_fluxes+n_neo+n_pellet) = source_pellet
            endif ! pellet

            !        vectors(inode,:,1) = (/ - R * u0_y ,   + R * u0_x ,   0.d0 /)
            !        vectors(inode,:,2) = (/ + ps_y /R * scalars(inode,7), - ps_x /R * scalars(inode,7), 0.d0 /) * Btot
            !        vectors(inode,:,3) = (/ - R * u0_y + ps_y /R * scalars(inode,7) * Btot, + R * u0_x - ps_x /R * scalars(inode,7) * Btot, 0.d0 /)

         endif

      enddo  ! nsub
   enddo     ! nsub

   do j=1,nsub-1
      do k=1,nsub-1
         ielm        = ielm+1
         ien(1,ielm) = inode - nsub*nsub + nsub*(j-1) + k-1       ! 0 based indices for VTK
         ien(2,ielm) = inode - nsub*nsub + nsub*(j  ) + k-1
         ien(3,ielm) = inode - nsub*nsub + nsub*(j  ) + k
         ien(4,ielm) = inode - nsub*nsub + nsub*(j-1) + k
      enddo
   enddo

enddo  ! n_elements

if (SI_units) then
  
  !===========================================================real values=============
  rho_norm = central_density*1.d20 * central_mass * mass_proton
  t_norm   = sqrt(MU_zero*rho_norm)

  !=================================================real values============
  do i=1,nnos 
    !============================================j_phi in MA/m2
    scalars(i,3) = scalars(i,3)/ MU_zero*1.e-6 
    !============================================density in 1e20m-3
    scalars(i,5) = scalars(i,5) * central_density
    if ( jorek_model .eq. 400 ) then
      !===========================================ion and electron temperatures in keV
      scalars(i,6) = scalars(i,6) / MU_zero / (central_density * 1d20) / EL_CHG /1.e3 !
      scalars(i,8) = scalars(i,8) / MU_zero / (central_density * 1d20) / EL_CHG /1.e3 !
    else
    !===========================================electron temperature in keV
      scalars(i,6) = scalars(i,6) / MU_zero / (central_density * 1d20) / EL_CHG /2./1.e3 !(assumes Te=Ti=T/2)
    endif
    !=====================================Vparal in km/s *Btot!!!
    scalars(i,7) = scalars(i,7)*Btot /t_norm/1.e3
    !=====================Pressure in kPa
    if (include_fluxes) scalars(i,n_var+1) = scalars(i,n_var+1) / MU_zero/1.e3
    if (include_neo) then
      !============================Er in kV/m
      scalars(i,n_var+n_fluxes+1) = F0*scalars(i,n_var+n_fluxes+1) / t_norm/1.e3
      !====================================Vtheta km/s
      scalars(i,n_var+n_fluxes+2) = scalars(i,n_var+n_fluxes+2) / t_norm/1.e3
      !===================================Vsound in km/s
      scalars(i,n_var+n_fluxes+5) = scalars(i,n_var+n_fluxes+5) / t_norm/1.e3
      !===================================Vneo in km/s
      scalars(i,n_var+n_fluxes+8) = scalars(i,n_var+n_fluxes+8) / t_norm/1.e3
      !===================================Vperp_e in km/s
      scalars(i,n_var+n_fluxes+10) = scalars(i,n_var+n_fluxes+10) / t_norm/1.e3
      ! mu_neo in SI units
      scalars(inode,n_var+n_fluxes+12) = scalars(inode,n_var+n_fluxes+12) / sqrt(rho_norm*MU_zero)
    endif
  enddo  ! nnos

endif ! SI_UNITS

!--------------------------------------------------- write the binary VTK file
etype = 9  ! for vtk_quad

lf = char(10) ! line feed character

#ifdef IBM_MACHINE
open(unit=ivtk,file='jorek_tmp.vtk',form='unformatted',access='stream')
#else
!open(unit=ivtk,file='jorek_tmp.vtk',form='unformatted',access='stream',convert='BIG_ENDIAN')
open(unit=ivtk,file='jorek_tmp.vtk',form='binary',convert='BIG_ENDIAN')
#endif

buffer = '# vtk DataFile Version 3.0'//lf    ; write(ivtk) trim(buffer)
buffer = 'vtk output'//lf                    ; write(ivtk) trim(buffer)
buffer = 'BINARY'//lf                        ; write(ivtk) trim(buffer)
buffer = 'DATASET UNSTRUCTURED_GRID'//lf     ; write(ivtk) trim(buffer)

! POINTS SECTION
write(str1(1:12),'(i12)') nnos
buffer = 'POINTS '//str1//'  float'//lf      ; write(ivtk) trim(buffer)
write(ivtk) ((real(xyz(i,j),4),i=1,3),j=1,nnos)

! CELLS SECTION
write(str1(1:12),'(i12)') nel            ! number of elements (cells)
write(str2(1:12),'(i12)') nel*(1+nnoel)  ! size of the following element list (nel*(nnoel+1))
buffer = lf//'CELLS '//str1//' '//str2//lf  ; write(ivtk) trim(buffer)
write(ivtk) (int(nnoel,4),(int(ien(i,j),4),i=1,nnoel),j=1,nel)

! CELL_TYPES SECTION
write(str1(1:12),'(i12)') nel   ! number of elements (cells)
buffer = lf//'CELL_TYPES'//str1//lf         ; write(ivtk) trim(buffer)
write(ivtk) (int(etype,4),i=1,nel)

! POINT_DATA SECTION
write(str1(1:12),'(i12)') nnos
buffer = lf//'POINT_DATA '//str1            ; write(ivtk) trim(buffer)

do i_var =1, n_scalars
  buffer = lf//'SCALARS '//scalar_names(i_var)//' float'//lf ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf
  write(ivtk) trim(buffer)
  write(ivtk) (real(scalars(i,i_var),4),i=1,nnos)
enddo

do i_var =1, n_vectors
  buffer = lf//lf//'VECTORS '//vector_names(i_var)//' float'//lf ; write(ivtk) trim(buffer)
  write(ivtk) ((real(vectors(j,i,i_var),4),i=1,3),j=1,nnos)
enddo

close(ivtk)

write(*,*) 'done.'
 
end program jorek2vtk
