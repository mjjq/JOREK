!> Program to convert a JOREK2 restart file into binary VTK format
program jorek2vtk

use mod_parameters, only: n_var, variable_names
use data_structure
use phys_module
use basis_at_gaussian
use diffusivities, only: get_dperp, get_zkperp
use pellet_module
use mpi_mod
use mod_bootstrap_functions
use corr_neg
use mod_import_restart
use mod_vtk
use mod_interp

implicit none

type (type_node_list)   ,     pointer :: node_list
type (type_element_list),     pointer :: element_list
type (type_bnd_element_list), pointer :: bnd_elm_list    

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
integer               :: n_fluxes, n_neo, n_bfield, n_vfield,n_pellet,n_bootstrap, n_psi_norm
integer               :: s_fluxes, s_neo, s_bfield, s_vfield,s_pellet,s_bootstrap, s_psi_norm
real*8                :: Jb, minRad,rho_norm,t_norm
integer               :: i_elm_axis, i_elm_xpoint(2), k_tor, ifail, ierr
logical               :: without_n0_mode, SI_units
logical               :: include_fluxes, include_neo, include_magnetic_field, include_velocity_field
logical               :: include_bootstrap, include_psi_norm
real*8                :: toroidal_angle
!====================== --- add the diagnostics Er, Vtheta and Vneo
real*8                :: Er, psi_abs, Vtheta, Btheta, Mach_par,Mach_pol,Vsound, Vneo
real*8                :: amu_neo_node, aki_neo_node
real*8                :: Vperp_e, Psi_tot

real*8                :: angle, source_volume, local_density, local_temperature, local_pressure, local_psi, local_source

!====================== --- Variables related to atomic physics terms (model 500 or 555)
logical               :: include_radiation
integer               :: n_radiation,s_radiation
real*8                :: Arad_bg, Brad_bg, Crad_bg, frad_bg, dfrad_bg_dT
real*8                :: T_corr, T_rad, coef_rad_1, Sion_T, eta_Sp, ksiion, Tion, LradDcont_T
real*8                :: LradDrays_T, coef_ion_1, coef_ion_2, coef_ion_3, S_ion_puiss
real*8                :: T_real8

!====================== --- Variables related to neutral density evolution (model 500 or 555)
logical               :: include_neutral_dens
integer               :: n_rn0, s_rn0
real*8                :: r0_corr, rn0_corr
real*8                :: r0_real8, rn0_real8
real*8                :: IonN, RecN, AblN, coef_rec_1, Srec_T

#ifdef fullmhd
!====================== --- Variables related to full mhd 
integer               :: n_fullmhd,s_fullmhd
#endif

integer, parameter :: nplot = 200
integer :: iplot, i_elm
real*8  :: stmp(200)
real*8  :: Rp_start, Zp_start, Rp_end, Zp_end
real*8  :: Rp, Zp, Rmin, Rmax, Zmin, Zmax, s_out, t_out, R_out, Z_out

namelist /vtk_params/ nsub, i_tor, i_plane, without_n0_mode, SI_units, &
                      include_fluxes, include_neo, include_magnetic_field, include_velocity_field,&
                      include_bootstrap, include_psi_norm



write(*,*) '***************************************'
write(*,*) '*       jorek2vtk                     *'
write(*,*) '***************************************'
write(*,*) ' if your VTK is smaller than expected,'
write(*,*) ' please consider the new parameters:'
write(*,*) '   -include_fluxes'
write(*,*) '   -include_neo'
write(*,*) '   -include_magnetic_field'
write(*,*) '   -include_velocity_field'
write(*,*) '   -include_bootstrap'
write(*,*) '   -include_psi_norm'
write(*,*) '***************************************'

call flush_it(6)

allocate(node_list)
allocate(element_list)
allocate(bnd_elm_list)

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Preset parameters
nsub                   = 5       ! Number of subdivisions of the cubic finite elements into linear pieces
i_tor                  = -1      ! If i_tor > 0, only this mode will be included in the vtk file...
i_plane                = 1       ! ... otherwise, all modes will be summed up at the toroidal plane i_plane
without_n0_mode        = .false. ! If true, do not include the n=0 mode (i_tor=1)
SI_units               = .false. ! when true, write variables in SI units
include_fluxes         = .false. ! include energy and density fluxes (or not)
include_neo            = .false. ! include neoclassical and more terms (or not)
include_magnetic_field = .false. ! include vector of magnetic field (or not)
include_velocity_field = .false. ! include vector of velocity field (or not)
include_bootstrap      = .false. ! include bootstrap current and averaged current
include_psi_norm       = .false.  ! include normalized flux

#if (JOREK_MODEL == 500)
include_radiation = .true.
include_neutral_dens = .true.
#endif

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
write(*,*) 'include_bootstrap =',include_bootstrap
write(*,*) 'include_psi_norm =', include_psi_norm

write(*,*) '-----------'
write(*,*) 'n_tor           =', n_tor
write(*,*) 'n_period        =', n_period
write(*,*) 'F0              =', F0
write(*,*) 'R_geo,Z_geo     =', R_geo, Z_geo
write(*,*)
call flush_it(6)

! --- Number of scalars to write to the VTK output file
n_scalars   = n_var
n_vectors   = 0
n_fluxes    = 0
n_neo       = 0
n_bfield    = 0
n_vfield    = 0
n_pellet    = 0
n_bootstrap = 0
n_psi_norm  = 0

if (include_fluxes) then
  n_fluxes  = 8
  s_fluxes  = n_scalars
  n_scalars = n_scalars + n_fluxes
endif
if (include_neo) then
  n_neo     = 10
  s_neo     = n_scalars
  n_scalars = n_scalars + n_neo
endif
if (include_magnetic_field) then
  n_bfield  = 1
  s_bfield  = n_vectors
  n_vectors = n_vectors + n_bfield
endif
if (include_velocity_field) then
  n_vfield  = 1
  s_vfield  = n_vectors
  n_vectors = n_vectors + n_vfield
endif
if (use_pellet) then
  n_pellet  = 2  ! pellet and pressuren
  s_pellet  = n_scalars
  n_scalars = n_scalars + n_pellet
endif
if (include_bootstrap) then
  n_bootstrap = 2
  s_bootstrap = n_scalars
  n_scalars   = n_scalars + n_bootstrap
endif
if (include_psi_norm) then
   n_psi_norm = 1
   s_psi_norm = n_scalars
   n_scalars  = n_scalars + n_psi_norm
endif

#if (JOREK_MODEL == 500)
    n_radiation = 0
 if (include_radiation) then
    n_radiation = 5
    s_radiation = n_scalars
    n_scalars   = n_scalars + n_radiation
 endif

    n_rn0 = 0
 if (include_neutral_dens) then
    n_rn0       = 2
    s_rn0       = n_scalars
    n_scalars   = n_scalars + n_rn0
 endif


#endif

#if fullmhd
 n_fullmhd = 3
 s_fullmhd = n_scalars
 n_scalars = n_scalars + n_fullmhd
#endif /*fullmhd*/

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
      scalar_names(6)='Te_keV      '
   endif
   scalar_names(7)='Vpar_km/s   '

#if (JOREK_MODEL == 500)
   scalar_names(8)='N_dens_1d20  '
#endif

endif

if ( SI_units ) then

  if (include_fluxes) then
     scalar_names(s_fluxes+1:s_fluxes+n_fluxes) = (/ &
      'P_kPa       ', 'E_flux_Kpar ', 'E_flux_kperp', 'E_flux_Vpar ', &
      'E_flux_Vperp', 'D_flux_Dperp', 'D_flux_Vpar ', 'D_flux_Vperp'/)
  endif
  if (include_neo) then
     scalar_names(s_neo+1:s_neo+n_neo) = (/ &
      'Er_kV/m     ', 'Vtheta_km/s ', 'Mach_par    ', 'Mach_pol    ', &
      'Vsound_km/s ', 'Btot_T      ', 'Vneo_km/s   ', 'Vperp_e_km/s', &
      'ki_neo      ', 'mu_neo      '/)
  endif
  if (include_bootstrap)then
     scalar_names(s_bootstrap+1:s_bootstrap+n_bootstrap) = (/'j_b_MA/m2   ', 'j_av_MA/m2  '/)
  endif

else

  if (include_fluxes) then
    scalar_names(s_fluxes+1:s_fluxes+n_fluxes) = (/ &
      'pressure    ', 'E_flux_Kpar ', 'E_flux_kperp', 'E_flux_Vpar ',&
      'E_flux_Vperp', 'D_flux_Dperp', 'D_flux_Vpar ', 'D_flux_Vperp'/)
  endif
  if (include_neo) then
    scalar_names(s_neo+1:s_neo+n_neo) = (/ &
      'Er          ', 'Vtheta      ', 'Mach_par    ', 'Mach_pol    ', &
      'Vsound      ', 'Btot        ', 'Vneo        ', 'Vperp_e     ', &
      'ki_neo      ', 'mu_neo      '/)
   endif
   if (include_bootstrap) then
      if (.not. bootstrap) write(*,*)'VTK WARNING: if you want the bootstrap, please set bootstrap=.t. in your input file!'
      scalar_names(s_bootstrap+1:s_bootstrap+n_bootstrap) = (/ 'j_bootstrap ', 'j_averaged  ' /)
   endif
!======================end SI units
endif

if (use_pellet) then
   scalar_names(s_pellet+1:s_pellet+n_pellet) =          (/ 'Pressure    ', 'Pellet      ' /)
endif
if (include_psi_norm) then
   scalar_names(s_psi_norm+1:s_psi_norm+n_psi_norm) = ('psi_norm    ')
endif

#if (JOREK_MODEL == 500)
 if (include_radiation) then
   scalar_names(s_radiation+1:s_radiation+n_radiation)                                   &
                  = (/ 'Ionis_Wm-3  ', 'Lin_radWm-3 ', 'Brems_Wm-3  ', 'Joule_Wm-3  ', 'Imp_bg_Wm-3 '/)
 endif

 if (include_neutral_dens) then
   scalar_names(s_rn0+1:s_rn0+n_rn0)&
                  = (/ 'IonN_s-1     ', 'RecN_s-1     '/)

 endif

#endif

#ifdef fullmhd
scalar_names(s_fullmhd+1:s_fullmhd+n_fullmhd) = (/  'B_phi       ', 'B_R         ', 'B_Z         '/)
#endif /*fullmhd*/

if (include_magnetic_field)  vector_names(s_bfield+1:s_bfield+n_bfield) = (/ 'B_R     ','B_Z     ','B_phi   '/)
if (include_velocity_field)  vector_names(s_vfield+1:s_vfield+n_vfield) = (/ 'V_R     ','V_Z     ','V_phi   '/)

do k_tor=1, n_tor
  mode(k_tor) = + int(k_tor / 2) * n_period
enddo

call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr)

call initialise_basis                              ! define the basis functions at the Gaussian points

nnos = nsub*nsub*element_list%n_elements
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

minRad = 0.0
if (bootstrap) then
  call bootstrap_find_minRad(node_list, element_list, R_axis, Z_axis, psi_axis, psi_bnd, minRad)
  call bootstrap_get_q_and_ft_splines(node_list, element_list, psi_axis, psi_xpoint, R_xpoint, Z_xpoint)
  call bootstrap_get_averaged_j_spline(node_list, element_list, psi_axis, psi_xpoint, R_xpoint, Z_xpoint)
endif

! --- You may choose to print your poloidal snapshot at a different toroidal angle
toroidal_angle = 0.d0 ! 2*PI / 6
if (toroidal_angle .ne. 0.d0) then
  do k_tor=1, n_tor
    mode(k_tor) = + int(k_tor / 2) * n_period
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
              + R_st*(Z_t*R_s + Z_s*R_t) + Z_ss*R_t**2 - R_ss*Z_t*R_t) / xjac

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
          amu_neo_node = 0.d0
          aki_neo_node = 0.d0

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
            scalars(inode,s_neo+1) = Er
            scalars(inode,s_neo+2) = Vtheta
            scalars(inode,s_neo+3) = Mach_par
            scalars(inode,s_neo+4) = Mach_pol
            scalars(inode,s_neo+5) = Vsound
            scalars(inode,s_neo+6) = Btot
            scalars(inode,s_neo+7) = Vneo
            scalars(inode,s_neo+8) = Vperp_e
            if (NEO) then
               if (num_neo_file) then
                  scalars(inode,s_neo+9) = aki_neo_node
                  scalars(inode,s_neo+10) = amu_neo_node
               else
                  scalars(inode,s_neo+9) = aki_neo_const
                  scalars(inode,s_neo+10) = amu_neo_const
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

      scalars(inode,s_fullmhd+1) = ( AZ_R - AR_Z )+ Fprof / R ! B_phi
      scalars(inode,s_fullmhd+2) = ( A3_Z - AZ_p )/ BigR      ! B_R
      scalars(inode,s_fullmhd+3) = ( AR_p - A3_R )/ BigR      ! B_Z
#endif /*fullmhd*/

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

              !  vectors(inode,:,1) = (/ - R * u0_y ,+ R * u0_x ,   0.d0 /)
              !  vectors(inode,:,2) = (/ + ps_y /R * V, - ps_x /R * V, F0/R * V /)
              !  vectors(inode,:,3) = (/ - R * u0_y + ps_y /R * V, + R * u0_x - ps_x /R * V, F0/R * V /)

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
          scalars(inode,s_fullmhd+1) = ( AZ_R - AR_Z )  + Fprof / R  ! B_phi
        else
          scalars(inode,s_fullmhd+1) = ( AZ_R - AR_Z )
        endif
        scalars(inode,s_fullmhd+2) = ( A3_Z - AZ_p )/ BigR  ! B_R
        scalars(inode,s_fullmhd+3) = ( AR_p - A3_R )/ BigR  ! B_Z
#endif /*fullmhd*/

      else  ! i_tor

        u_sum   = 0.d0; u_x  = 0.d0; u_y  = 0.d0; u_p  = 0.d0
        psi_sum = 0.d0; ps_x = 0.d0; ps_y = 0.d0; ps_p = 0.d0
        zj_sum  = 0.d0; zj_x = 0.d0; zj_y = 0.d0; zj_p = 0.d0
        T_sum   = 0.d0; TT_x = 0.d0; TT_y = 0.d0; TT_p = 0.d0
        zn_sum  = 0.d0; zn_x = 0.d0; zn_y = 0.d0; zn_p = 0.d0;
        Ti_sum  = 0.d0; Ti_x = 0.d0; Ti_y = 0.d0; Ti_p = 0.d0
        Te_sum  = 0.d0; Te_x = 0.d0; Te_y = 0.d0; Te_p = 0.d0
        w_sum   = 0.d0; w_x  = 0.d0; w_y  = 0.d0; w_p  = 0.d0; w_xx = 0.d0; w_yy = 0.d0

#ifdef fullmhd
        !reinitialize Bphi,BR,B_Z for loop over all modes
        scalars(inode,s_fullmhd+1) = 0.  
        scalars(inode,s_fullmhd+2) = 0.  
        scalars(inode,s_fullmhd+3) = 0.  
#endif /*fullmhd*/

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
          endif
          if ( jorek_model .eq. 400 ) then
             call interp(node_list,element_list,i,6,i_tor,s,t,Ti,Ti_s,Ti_t,Ti_st,Ti_ss,Ti_tt)
             call interp(node_list,element_list,i,8,i_tor,s,t,Te,Te_s,Te_t,Te_st,Te_ss,Te_tt)
          endif

          psi_sum = psi_sum + psi * HZ(i_tor,i_plane)
          zj_sum  = zj_sum  + zj  * HZ(i_tor,i_plane)
          u_sum   = u_sum   + U   * HZ(i_tor,i_plane)
          w_sum   = w_sum   + w   * HZ(i_tor,i_plane)
          zn_sum  = zn_sum  + RHO * HZ(i_tor,i_plane)
          if ( jorek_model .eq. 400 ) then
            Ti_sum  = Ti_sum + Ti * HZ(i_tor,i_plane)
            Te_sum  = Te_sum + Te * HZ(i_tor,i_plane)
          else
            Ti_sum  = Ti_sum + 0.5d0*TT * HZ(i_tor,i_plane)
            Te_sum  = Te_sum + 0.5d0*TT * HZ(i_tor,i_plane)
          endif

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
            scalars(inode,s_fullmhd+1) = scalars(inode,n_var+1) + ( AZ_R - AR_Z )  + Fprof / R  ! B_phi
          else
            scalars(inode,s_fullmhd+1) = scalars(inode,n_var+1) + ( AZ_R - AR_Z )     * HZ(i_tor,i_plane)
          endif

          scalars(inode,s_fullmhd+2) = scalars(inode,n_var+2) + ( A3_Z - AZ_p )/ BigR * HZ(i_tor,i_plane)  ! B_R
          scalars(inode,s_fullmhd+3) = scalars(inode,n_var+3) + ( AR_p - A3_R )/ BigR * HZ(i_tor,i_plane)  ! B_Z
#endif /*fullmhd*/

          if ((xjac .gt. 1.d-6)) then  ! avoid the axis

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
             else
               Ti_x  = TT_x / 2.0
               Ti_y  = TT_y / 2.0
               Ti_p  = TT_p / 2.0
               Te_x  = TT_x / 2.0
               Te_y  = TT_y / 2.0
               Te_p  = TT_p / 2.0
             end if

             zn_x = zn_x  + (   Z_t * RHO_s - Z_s * RHO_t ) / xjac * HZ(i_tor,i_plane)
             zn_y = zn_y  + ( - R_t * RHO_s + R_s * RHO_t ) / xjac * HZ(i_tor,i_plane)
             zn_p = zn_p  + RHO * HZ_p(i_tor,i_plane)

             w_x  = w_x   + (   Z_t * w_s - Z_s * w_t )     / xjac * HZ(i_tor,i_plane)
             w_y  = w_y   + ( - R_t * w_s + R_s * w_t )     / xjac * HZ(i_tor,i_plane)

             w_xx = w_xx  + (w_ss * Z_t**2 - 2.d0*w_st * Z_s*Z_t + w_tt * Z_s**2        &
                  + w_s * (Z_st*Z_t - Z_tt*Z_s )                                        &
                  + w_t * (Z_st*Z_s - Z_ss*Z_t ) )     / xjac**2                        &
                  - xjac_x * (w_s* Z_t - w_t * Z_s)  / xjac**2

             w_yy = w_yy  + (w_ss * R_t**2 - 2.d0*w_st * R_s*R_t + w_tt * R_s**2        &
                  + w_s * (R_st*R_t - R_tt*R_s )                                        &
                  + w_t * (R_st*R_s - R_ss*R_t ) )       / xjac**2                      &
                  - xjac_y * (- w_s * R_t + w_t * R_s )  / xjac**2

          endif ! xjac

        enddo  ! end loop toroidal harmonics

        Psi_tot = 0.d0
        do i_tor =1, n_tor
           call interp(node_list,element_list,i,1,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
           Psi_tot = Psi_tot + P * HZ(i_tor,i_plane)
        enddo

        psi_norm = (Psi_tot - psi_axis)/(psi_bnd - psi_axis)
        if ((psi_norm .lt. 1.d0) .and. (xpoint) .and. (Z .lt. Z_xpoint(1)) .and. (xcase .ne. 2)) then
           psi_norm = 2.d0 - psi_norm
        endif
        if ((psi_norm .lt. 1.d0) .and. (xpoint) .and. (Z .gt. Z_xpoint(2)) .and. (xcase .ne. 1)) then
           psi_norm = 2.d0 - psi_norm
        endif

        if (include_bootstrap) then
          call bootstrap_current(minRad, R, Z, R_axis, Z_axis, psi_axis, R_xpoint, Z_xpoint, psi_bnd, psi_norm,&
                                 psi_sum, ps_x, ps_y, zn_sum,  zn_x, zn_y,      &
                                 Ti_sum,  Ti_x, Ti_y, Te_sum,  Te_x, Te_y, Jb   )
          scalars(inode,s_bootstrap+1) = Jb
          scalars(inode,s_bootstrap+2) = bootstrap_spline3_eval(n_spline_vtk-1,psi_knots_vtk,j_knots_vtk,j_spline_vtk,psi_norm)
          ! --- The JOREK bootstrap is not constant on flux surface
          ! --- Because J_jorek = R*J_physical, and the physical bootstrap is constant on a surface
          ! --- Hence, it makes more sense to look at R*j_average to compare with the bootstrap...
          scalars(inode,s_bootstrap+2) = R* scalars(inode,n_var+2+n_fluxes+n_neo+n_pellet) / R_axis
        else
          Jb = 0.d0
        endif

        v_perp  = R * sqrt(u_x*u_x + u_y * u_y)
        Btot    = sqrt(F0**2 + ps_x**2 + ps_y**2) / BigR
        D_prof  = get_dperp (psi_norm)
        ZK_prof = get_zkperp(psi_norm)

        ZKpar_T = ZK_par * ((max( scalars(inode,6), T_min ))/T_0)**2.5

        grad_psi = sqrt(ps_x*ps_x + ps_y*ps_y)

        if ((SI_units) .and. (jorek_model .ge. 300)) then
          scalars(inode,7) = scalars(inode,7) * sign(Btot,F0)   ! with si-units= .f. gives jorek variable, otherwise physical v_par
        endif

        !   'E_flux_Kpar ','E_flux_kperp','E_flux_Vpar ','E_flux_Vperp','D_flux_Dperp','D_flux_Vpar ','D_flux_Vperp'/)

        if (include_fluxes) then

          scalars(inode,s_fluxes+1)   = scalars(inode,5) * scalars(inode,6)

          if (grad_psi .ne. 0.d0) then

            scalars(inode,s_fluxes+2)  = ZKpar_T * ( F0 * TT_p / BigR**2  + (TT_x * ps_y - TT_y * ps_x) / BigR ) / Btot

            scalars(inode,s_fluxes+3)  = ZK_prof * (TT_x * ps_x + TT_y * ps_y) / grad_psi

            scalars(inode,s_fluxes+4)  = scalars(inode,5) * scalars(inode,6) * scalars(inode,7) * Btot

            scalars(inode,s_fluxes+5)  = BigR   * (u_x * ps_y - u_y * ps_x) / sqrt(ps_x*ps_x + ps_y*ps_y) * scalars(inode,5) * scalars(inode,6)

            scalars(inode,s_fluxes+6)  = D_prof * (zn_x * ps_x + zn_y * ps_y) / sqrt(ps_x*ps_x + ps_y*ps_y)

            scalars(inode,s_fluxes+7)  = scalars(inode,5) * scalars(inode,7) * Btot

            scalars(inode,s_fluxes+8)  = BigR   * (u_x * ps_y - u_y * ps_x) / sqrt(ps_x*ps_x + ps_y*ps_y) * scalars(inode,5)

          endif ! grad_psi
        endif ! include_fluxes

        if (include_psi_norm) then
           scalars(inode,s_psi_norm+1) = psi_norm
        endif

        if (use_pellet) then

           local_density     = scalars(inode,5)
           local_temperature = scalars(inode,6)/2.
           local_psi         = scalars(inode,1)
           local_pressure    = local_density * local_temperature
           scalars(inode,s_pellet+1)  = local_pressure

           angle = 0.0
           call pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, &
                pellet_theta, R, Z, local_psi, angle, &
                local_density,local_temperature, &
                central_density, pellet_particles, pellet_density, &
                total_pellet_volume, local_source, source_volume)

           scalars(inode,n_var+n_fluxes+n_neo+n_pellet) = local_source
        endif ! use_pellet

        !	 vectors(inode,:,1) = (/ - R * u0_y ,	+ R * u0_x ,   0.d0 /)
        !	 vectors(inode,:,2) = (/ + ps_y /R * scalars(inode,7), - ps_x /R * scalars(inode,7), 0.d0 /) * Btot
        !	 vectors(inode,:,3) = (/ - R * u0_y + ps_y /R * scalars(inode,7) * Btot, + R * u0_x - ps_x /R * scalars(inode,7) * Btot, 0.d0 /)

     endif ! i_tor from 1 to n_tor

  enddo  ! nsub
enddo     ! nsub

do j=1,nsub-1
   do k=1,nsub-1
      ielm	  = ielm+1
      ien(1,ielm) = inode - nsub*nsub + nsub*(j-1) + k-1       ! 0 based indices for VTK
      ien(2,ielm) = inode - nsub*nsub + nsub*(j  ) + k-1
      ien(3,ielm) = inode - nsub*nsub + nsub*(j  ) + k
      ien(4,ielm) = inode - nsub*nsub + nsub*(j-1) + k
   enddo
enddo

enddo  ! n_elements

#if (JOREK_MODEL == 500)

 if (include_radiation) then

   do i=1,nnos

    coef_ion_3 = 27.2d0*EL_CHG*MU_ZERO*central_density*1.d20
    coef_ion_2 = 0.232d0
    coef_ion_1 = (MU_ZERO*central_mass*MASS_PROTON)**(0.5d0)*0.2917d-13*(central_density*1.d20)**(1.5d0)
    S_ion_puiss = 3.9d-1

    ksiion = ksi_ion * central_density * 1.d20

    T_real8 = scalars(i,6)

    T_corr = corr_neg_temp(T_real8)

    Tion = corr_neg_temp(T_real8,(/1.d-5,0.3/))/(2.d0)

    T_rad = corr_neg_temp(T_real8)/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)

    Sion_T = coef_ion_1*((coef_ion_3/Tion)**S_ion_puiss)*1/(coef_ion_2+coef_ion_3/Tion)*exp(-coef_ion_3/Tion)

    coef_rad_1 = 2.d0/(3.d0)*MU_ZERO**1.5d0*(central_mass*MASS_PROTON)**0.5d0*(central_density*1.d20)**2.5d0

    LradDcont_T = coef_rad_1*5.37d-37*(1.d1)**(-1.5d0)*(1.d0)**2*sqrt(T_rad) ! Only Bremsstrahlung contribution

    LradDrays_T = coef_rad_1*(1.d1)**(-29.44d0*exp(-(log10(T_rad)-4.4283d0)**2.d0/(2.d0*(2.8428d0)**2.d0)) &
                                     -60.947d0*exp(-(log10(T_rad)+2.0835d0)**2.d0/(2.d0*(0.9048d0)**2.d0)) &
                                     -24.067d0*exp(-(log10(T_rad)+0.7363d0)**2.d0/(2.d0*(2.1700d0)**2.d0)))

    eta_Sp = 1.65d-9*17*(1.d-3*T_rad)**(-1.5d0) &
                            *(central_mass*MASS_PROTON*central_density * 1.d20/MU_ZERO)**(0.5d0)

    scalars(i,s_radiation+1) = ksiion * scalars(i,5) * scalars(i,8) * Sion_T
    scalars(i,s_radiation+2) = scalars(i,5) * scalars(i,8) * LradDrays_T
    scalars(i,s_radiation+3) = LradDcont_T * scalars(i,5)**2.d0
    scalars(i,s_radiation+4) = (2/(3 * BigR**2)) * eta_Sp * scalars(i,3)**2.d0

 !--------------------------------------------------------
   ! --- Radiation from background impurity
   !--------------------------------------------------------

    Arad_bg = 2.4d-31
    Brad_bg = 20.
    Crad_bg = 0.8

    frad_bg = (2./3.)*(1./(central_mass*MASS_PROTON))                               &
               *((MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)**(1.5d0)) &
               *nimp_bg*Arad_bg*exp(-((log(T_rad)-log(Brad_bg))**2.)/Crad_bg**2.)

    dfrad_bg_dT = -(1./3.)*((MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)**(0.5d0)) &
                   *(1./EL_CHG)*2.*(nimp_bg*Arad_bg/Crad_bg**2.)*(log(T_rad)-log(Brad_bg))     &
                   *(1./T_rad)*exp(-((log(T_rad)-log(Brad_bg))**2.)/Crad_bg**2.)

    scalars(i,s_radiation+5) = scalars(i,5) * frad_bg

   enddo
  endif

 if (include_neutral_dens) then

   coef_rec_1 = (MU_ZERO*central_mass*MASS_PROTON)**(0.5d0)*(central_density * 1.d20)**(1.5d0)

   coef_ion_3 = 27.2d0*EL_CHG*MU_ZERO*central_density*1.d20
   coef_ion_2 = 0.232d0
   coef_ion_1 = (MU_ZERO*central_mass*MASS_PROTON)**(0.5d0)*0.2917d-13*(central_density*1.d20)**(1.5d0)
   S_ion_puiss = 3.9d-1

   do i=1,nnos

     T_real8 = scalars(i,6)

     T_corr = corr_neg_temp(T_real8)


     Srec_T    = coef_rec_1 * 0.7d-19 * (13.6*(2*EL_CHG*MU_ZERO*central_density*1.d20))**(0.5d0) * (T_corr/(2.d0))**(-0.5d0)

     Tion = corr_neg_temp(T_real8,(/1.d-5,0.3/))/(2.d0)
 
     Sion_T = coef_ion_1*((coef_ion_3/Tion)**S_ion_puiss)*1/(coef_ion_2+coef_ion_3/Tion)*exp(-coef_ion_3/Tion)


     r0_real8  = scalars(i,5)
     rn0_real8 = scalars(i,8)

     r0_corr   = corr_neg_dens(r0_real8)
     rn0_corr  = corr_neg_dens(rn0_real8, (/ 0.d-5, 1.d-5 /))

     IonN      = -(r0_corr) * (rn0_corr) * Sion_T
     RecN      = (r0_corr)**2 * Srec_T 

     scalars(i,s_rn0+1) = IonN
     scalars(i,s_rn0+2) = RecN


   end do
 end if


#endif /*(JOREK_MODEL == 500)*/


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
    scalars(i,7) = scalars(i,7) /t_norm/1.e3
#if (JOREK_MODEL == 500)
    !===================================== Neutral density in 1e20m-3
    scalars(i,8) = scalars(i,8) * central_density
#endif
    !=====================Pressure in kPa
    if (include_fluxes) scalars(i,s_fluxes+1) = scalars(i,s_fluxes+1) / MU_zero/1.e3
    if (include_neo) then
      !============================Er in kV/m
      scalars(i,s_neo+ 1) = F0*scalars(i,s_neo+ 1) / t_norm/1.e3
      !====================================Vtheta km/s
      scalars(i,s_neo+ 2) =    scalars(i,s_neo+ 2) / t_norm/1.e3
      !===================================Vsound in km/s
      scalars(i,s_neo+ 5) =    scalars(i,s_neo+ 5) / t_norm/1.e3
      !===================================Vneo in km/s
      scalars(i,s_neo+ 7) =    scalars(i,s_neo+ 7) / t_norm/1.e3
      !===================================Vperp_e in km/s
      scalars(i,s_neo+ 8) =    scalars(i,s_neo+ 8) / t_norm/1.e3
      ! ===================================mu_neo in SI units
      scalars(i,s_neo+10) =    scalars(i,s_neo+10) / sqrt(rho_norm*MU_zero)
    endif
    !============================================j_bootstrap, javeraged in MA/m2

    if (include_bootstrap) then
    scalars(i,s_bootstrap+1)=scalars(i,s_bootstrap+1)/MU_zero*1.e-6
    scalars(i,s_bootstrap+2)=scalars(i,s_bootstrap+2)/MU_zero*1.e-6
    endif
  !========================================================

#if (JOREK_MODEL == 500)

 if (include_radiation) then

    coef_ion_3 = 27.2d0*EL_CHG*MU_zero*central_density*1.d20
    coef_ion_2 = 0.232d0
    coef_ion_1 = 0.2917d-13 !(MU_ZERO*central_mass*MASS_PROTON)**(0.5d0)*0.2917d-13*(central_density*1.d20)**(1.5d0)
    S_ion_puiss = 3.9d-1

    ksiion = ksi_ion * central_density * 1.d20

    T_real8 = scalars(i,6)*1.e3*2.*EL_CHG*MU_zero*(central_density * 1.d20)
    ! ======= T_real8 in JOREK units

    Tion = corr_neg_temp(T_real8,(/1.d-5,0.3/))/(2.d0)

    T_rad = corr_neg_temp(T_real8)/(2.d0*EL_CHG*MU_zero*central_density*1.d20)

    Sion_T = coef_ion_1*((coef_ion_3/Tion)**S_ion_puiss)*1/(coef_ion_2+coef_ion_3/Tion)*exp(-coef_ion_3/Tion)

    coef_rad_1 = 1.d0 !2.d0/(3.d0)*MU_ZERO**1.5d0*(central_mass*MASS_PROTON)**0.5d0*(central_density*1.d20)**2.5d0

    LradDcont_T = coef_rad_1*5.37d-37*(1.d1)**(-1.5d0)*(1.d0)**2*sqrt(T_rad) ! Only Bremsstrahlung contribution

    LradDrays_T = coef_rad_1*(1.d1)**(-29.44d0*exp(-(log10(T_rad)-4.4283d0)**2.d0/(2.d0*(2.8428d0)**2.d0))  &
                                      -60.947d0*exp(-(log10(T_rad)+2.0835d0)**2.d0/(2.d0*(0.9048d0)**2.d0)) &
                                      -24.067d0*exp(-(log10(T_rad)+0.7363d0)**2.d0/(2.d0*(2.1700d0)**2.d0)))

    eta_Sp = 1.65d-9*17*(1.d-3*T_rad)**(-1.5d0)

       scalars(i,s_radiation+1) = ksiion* (1.5d0)/(MU_zero*central_density*1.d20)      &
                                           * scalars(i,5) * 1.d20 * scalars(i,8) * 1.d20 * Sion_T

       scalars(i,s_radiation+2) = scalars(i,5)* 1.d20 * scalars(i,8) * 1.d20 * LradDrays_T

       scalars(i,s_radiation+3) = LradDcont_T * (scalars(i,5)*1.d20)**2.d0

       scalars(i,s_radiation+4) = eta_Sp * (1.d6*scalars(i,3))**2.d0

 !--------------------------------------------------------
   ! --- Radiation from background impurity
   !--------------------------------------------------------

    Arad_bg = 2.4d-31
    Brad_bg = 20.
    Crad_bg = 0.8

    frad_bg = nimp_bg * Arad_bg*exp(-((log(T_rad)-log(Brad_bg))**2.)/Crad_bg**2.)

    dfrad_bg_dT = -(1./3.)*((MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)**(0.5d0)) &
                   *(1./EL_CHG)*2.*(nimp_bg*Arad_bg/Crad_bg**2.)*(log(T_rad)-log(Brad_bg))     &
                   *(1./T_rad)*exp(-((log(T_rad)-log(Brad_bg))**2.)/Crad_bg**2.)

    scalars(i,s_radiation+5) = scalars(i,5)*1.d20 * frad_bg

 endif

#endif /*(JOREK_MODEL == 500)*/

  enddo  ! nnos

endif ! SI_UNITS

!--------------------------------------------------- write the binary VTK file
etype = 9  ! for vtk_quad

call write_vtk('jorek_tmp.vtk',xyz,ien,etype,scalar_names,scalars,vector_names,vectors)

write(*,*) 'done.'

end program jorek2vtk
