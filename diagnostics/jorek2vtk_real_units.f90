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
real*8                :: Psi,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt
real*8                :: ZJ,ZJ_s,ZJ_t,ZJ_st,ZJ_ss,ZJ_tt
real*8                :: U,U_s,U_t,U_st,U_ss,U_tt
real*8                :: W,W_s,W_t,W_st,W_ss,W_tt
real*8                :: RHO,RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt
real*8                :: TT,TT_s,TT_t,TT_st,TT_ss,TT_tt
real*8                :: Ti,Ti_s,Ti_t,Ti_st,Ti_ss,Ti_tt
real*8                :: Te,Te_s,Te_t,Te_st,Te_ss,Te_tt
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
integer 	      :: i_find, i_elm_find(8)
real*8  	      :: Router,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8  	      :: Zouter,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8  	      :: s_find(8), t_find(8)
real*8                :: Jb, minRad,rho_norm,t_norm
integer               :: i_elm_axis, i_elm_xpoint(2), k_tor, ifail, ierr
logical               :: without_n0_mode
!====================== --- add the diagnostics Er, Vtheta and Vneo
real*8                :: Er, psi_abs, Vtheta, Btheta, Mach_par,Mach_pol,Vsound, Vneo
real*8                :: amu_neo_node, aki_neo_node
real*8                :: Vperp_e !, qprofile, small_r
real*8                :: density_tot,density_in,density_out,pressure,pressure_in,pressure_out

namelist /vtk_params/ nsub, i_tor, i_plane, without_n0_mode
GAMMA=5./3.
write(*,*) 'jorek2vtk'
call flush_it(6)
allocate(node_list)
allocate(element_list)

! --- Initialise input parameters and read the input namelist.
my_id     = 0
mach_par=0.d0
mach_pol=0.d0
vsound=0.d0

call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Preset parameters
nsub      = 5             ! Number of subdivisions of the cubic finite elements into linear pieces
i_tor     = -1            ! If i_tor > 0, only this mode will be included in the vtk file...
!i_tor     = 3 
i_plane   = 1             ! ... otherwise, all modes will be summed up at the toroidal plane i_plane
without_n0_mode = .false. ! If true, do not include the n=0 mode (i_tor=1)
!without_n0_mode = .true.

! --- Read parameters from namelist file 'vtk.nml' if it exists
open(42, file='vtk.nml', action='read', status='old', iostat=ierr)
if ( ierr == 0 ) then
  write(*,*) 'Reading parameters from vtk.txt namelist.'
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
write(*,*) '-----------'
write(*,*) 'n_tor           =', n_tor
write(*,*) 'n_period        =', n_period
write(*,*) 'F0              =', F0
write(*,*) 'R_geo,Z_geo     =', R_geo, Z_geo
write(*,*)
call flush_it(6)

! --- Number of scalars to write to the VTK output file
n_scalars = n_var + 18!19

! --- Number of vectors to write to the VTK output file
n_vectors = 0

allocate(scalar_names(n_scalars), vector_names(n_vectors))

grad_psi = 0.d0

scalar_names(1:n_var) = variable_names(1:n_var)
scalar_names(3)='j_MA/m2     '
scalar_names(5)='n_e20m-3    '
scalar_names(6)='Te_keV      '
scalar_names(7)='Vpar_km/s   '
scalar_names(n_var+1:n_scalars) = (/ &
  'P_kPa       ', 'E_flux_Kpar ', 'E_flux_kperp', 'E_flux_Vpar ', 'E_flux_Vperp', 'D_flux_Dperp', &
  'D_flux_Vpar ', 'D_flux_Vperp', 'Er_kV/m     ', 'Vtheta_km/s ', 'Mach_par    ', 'Mach_pol    ', &
  'Vsound_km/s ', 'Btot_T      ', 'J-bootstrap ', 'Vneo_km/s   ', 'psi_norm    ', 'Vperp_e     '/)!, &
!  'qprofile    '/)

!vector_names = (/ 'v_perp  ','v_par   ','V_tot   '/)

do k_tor=1, n_tor
  mode(k_tor) = + int(k_tor / 2) * n_period
enddo

call import_restart(node_list,element_list, 'jorek_restart.rst',rst_format, ierr)

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

if (jorek_model .eq. 400) then      
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
         ! NB: i_tor = 2 corresponds to the cosine component of the harmonics ?
         !     i_tor = 3 corresponds to the sinus component of the harmonics ?
         
         ! save old values
         i_tor_old = i_tor
         i_tor     = 1
         ! compute all derivatives, as in loop below
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

            u0_x  = (   Z_t * U_s - Z_s * U_t ) / xjac
            u0_y  = ( - R_t * U_s + R_s * U_t ) / xjac

            ps_x  = (   Z_t * PS_s - Z_s * PS_t ) / xjac
            ps_y  = ( - R_t * PS_s + R_s * PS_t ) / xjac

            TT_x  = (   Z_t * TT_s - Z_s * TT_t ) / xjac
            TT_y  = ( - R_t * TT_s + R_s * TT_t ) / xjac

            zj_x  = (   Z_t * ZJ_s - Z_s * ZJ_t ) / xjac
            zj_y  = ( - R_t * ZJ_s + R_s * ZJ_t ) / xjac

            RHO_x = (   Z_t * RHO_s - Z_s * RHO_t ) / xjac
            RHO_y  = ( - R_t * RHO_s + R_s * RHO_t ) / xjac

            !*** compute diagnostics ***
            psi_abs = sqrt(ps_x*ps_x + ps_y * ps_y)
            Btheta  = (psi_abs/R)          
            Vtheta  = 0.0
            Vneo    = 0.0
            Er      = 0.0
            if ((psi_abs .gt. 1.d-6.and. RHO.gt.1.d-6.and.abs(Btheta).gt.1.d-6))then
               Vtheta = -1./Btheta*((u0_x+ tauIC/RHO*(TT_x*RHO+RHO_x*TT))*ps_x + &
                    (u0_y+tauIC/RHO*(TT_y*RHO+RHO_y*TT))*ps_y)+V*Btheta
               Vperp_e = -1./Btheta*((u0_x- tauIC/RHO*(TT_x*RHO+RHO_x*TT))*ps_x + &
                    (u0_y-tauIC/RHO*(TT_y*RHO+RHO_y*TT))*ps_y)

               if (NEO) then
!                  num_neo_file= ( neo_file /= 'none')
                  if (num_neo_file) then
                     write(*,*) 'neo_file=',neo_file
                     write(*,*) 'using ki and mui profiles from file "'//trim(neo_file)//'"'
                     call neo_coef( xpoint, xcase, Z, Z_xpoint, Psi ,psi_axis,psi_bnd, &
                          amu_neo_node, aki_neo_node)
                     Vneo   = aki_neo_node/Btheta*tauIC*(ps_x*TT_x+ps_y*TT_y)
                  else
                     Vneo   = aki_neo_const/Btheta*tauIC*(ps_x*TT_x+ps_y*TT_y)
                  endif
               endif
               Er     = -(u0_x*ps_x + u0_y * ps_y)/psi_abs
!======================Mach number===========================================              
               Btot = sqrt(F0**2 + ps_x**2 + ps_y**2) / BigR 

               Vsound=sqrt(GAMMA*TT)/Btot
               ! why term in boundry conditions ? from Vtheta? BigR**2 * U0_s /ps0_s==================????
               Mach_par=V/Vsound             
               Mach_pol=Vtheta/Vsound
            endif
         endif
         ! save those specific values of axisymmetric parameters
         if (grad_psi .ne. 0.d0) then
            scalars(inode,n_var+9) = Er
            scalars(inode,n_var+10) = Vtheta
            scalars(inode,n_var+11) = Mach_par
            scalars(inode,n_var+12) = Mach_pol
            scalars(inode,n_var+13) = Vsound
            scalars(inode,n_var+14) = Btot
            scalars(inode,n_var+16) = Vneo ! number n_var+15 is jbootstrap, see below.
            scalars(inode,n_var+18) = Vperp_e ! number n_var+17 is psi_norm, see below.
            !scalars(inode,n_var+19) = qprofile 
         endif
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

               u0_x  = (   Z_t * U_s - Z_s * U_t ) / xjac
               u0_y  = ( - R_t * U_s + R_s * U_t ) / xjac

               ps_x  = (   Z_t * PS_s - Z_s * PS_t ) / xjac
               ps_y  = ( - R_t * PS_s + R_s * PS_t ) / xjac

               TT_x  = (   Z_t * TT_s - Z_s * TT_t ) / xjac
               TT_y  = ( - R_t * TT_s + R_s * TT_t ) / xjac

               zj_x  = (   Z_t * ZJ_s - Z_s * ZJ_t ) / xjac
               zj_y  = ( - R_t * ZJ_s + R_s * ZJ_t ) / xjac

               !*** compute diagnostics ***
               v_perp  = R * sqrt(u0_x*u0_x + u0_y * u0_y)

               !	  vectors(inode,:,1) = (/ - R * u0_y ,   + R * u0_x ,   0.d0 /)
               !          vectors(inode,:,2) = (/ + ps_y /R * V, - ps_x /R * V, F0/R * V /)
               !          vectors(inode,:,3) = (/ - R * u0_y + ps_y /R * V, + R * u0_x - ps_x /R * V, F0/R * V /)

               psi_J = (Ps_s * ZJ_t - PS_t * ZJ_s ) / xjac
               R_p   = (2.d0 * R * (R_s * (RHO_t * TT + RHO * TT_t) - R_t * (RHO_s * TT + RHO * TT_s) )) / xjac
               error = psi_J - R_p

            endif

         else

                	u0_x  = 0.d0; u0_y  = 0.d0;
	    psi = 0.d0; ps_x  = 0.d0; ps_y  = 0.d0; 
	        	zj_x  = 0.d0; zj_y  = 0.d0;
                	TT_x  = 0.d0; TT_y  = 0.d0; TT_p  = 0.d0
	    RHO = 0.d0; RHO_x = 0.d0; RHO_y = 0.d0; RHO_p = 0.d0;
            Ti  = 0.d0; Ti_x  = 0.d0; Ti_y  = 0.d0; Ti_p  = 0.d0
	    Te  = 0.d0; Te_x  = 0.d0; Te_y  = 0.d0; Te_p  = 0.d0
                	w0_x  = 0.d0; w0_y  = 0.d0; w0_xx = 0.d0; w0_yy = 0.d0

            do i_tor = 1, n_tor

               if ( ( i_tor == 1 ) .and. ( without_n0_mode ) ) cycle ! Do not include the n=0 mode

               do m=1,n_var
                  call interp(node_list,element_list,i,m,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
                  scalars(inode,m) = scalars(inode,m) + P * HZ(i_tor,i_plane)
               enddo

               call interp(node_list,element_list,i,1,i_tor,s,t,Psi_00,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt)
               call interp(node_list,element_list,i,2,i_tor,s,t,U     ,U_s,U_t,U_st,U_ss,U_tt)
               call interp(node_list,element_list,i,3,i_tor,s,t,ZJ    ,ZJ_s,ZJ_t,ZJ_st,ZJ_ss,ZJ_tt)
               call interp(node_list,element_list,i,4,i_tor,s,t,W     ,W_s,W_t,W_st,W_ss,W_tt)
               call interp(node_list,element_list,i,5,i_tor,s,t,RHO_00,RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt)
               call interp(node_list,element_list,i,6,i_tor,s,t,TT    ,TT_s,TT_t,TT_st,TT_ss,TT_tt)
               if ( jorek_model >= 300 ) then
                  call interp(node_list,element_list,i,7,i_tor,s,t,V,V_s,V_t,V_st,V_ss,V_tt)
               else
                  V=0; V_s=0; V_t=0; V_st=0; V_ss=0; V_tt=0
               end if
               if ( jorek_model .eq. 400 ) then
                  call interp(node_list,element_list,i,6,i_tor,s,t,Ti_00,Ti_s,Ti_t,Ti_st,Ti_ss,Ti_tt)
                  call interp(node_list,element_list,i,8,i_tor,s,t,Te_00,Te_s,Te_t,Te_st,Te_ss,Te_tt)
               end if

               psi   = psi    + psi_00 * HZ(i_tor,i_plane)
               RHO   = RHO    + RHO_00 * HZ(i_tor,i_plane)
               if ( jorek_model .eq. 400 ) then
                 Ti  = Ti     + Ti_00  * HZ(i_tor,i_plane)
                 Te  = Te     + Te_00  * HZ(i_tor,i_plane)
               end if
               if ((xjac .gt. 1.d-6)) then      ! avoid the axis

                  u0_x  = u0_x   + (   Z_t * U_s - Z_s * U_t )     / xjac * HZ(i_tor,i_plane)
                  u0_y  = u0_y   + ( - R_t * U_s + R_s * U_t )     / xjac * HZ(i_tor,i_plane)

                  ps_x  = ps_x   + (   Z_t * PS_s - Z_s * PS_t )   / xjac * HZ(i_tor,i_plane)
                  ps_y  = ps_y   + ( - R_t * PS_s + R_s * PS_t )   / xjac * HZ(i_tor,i_plane)

                  zj_x  = zj_x   + (   Z_t * ZJ_s - Z_s * ZJ_t )   / xjac * HZ(i_tor,i_plane)
                  zj_y  = zj_y   + ( - R_t * ZJ_s + R_s * ZJ_t )   / xjac * HZ(i_tor,i_plane)

                  TT_x  = TT_x   + (   Z_t * TT_s - Z_s * TT_t )   / xjac * HZ(i_tor,i_plane)
                  TT_y  = TT_y   + ( - R_t * TT_s + R_s * TT_t )   / xjac * HZ(i_tor,i_plane)
                  TT_p  = TT_p   + TT * HZ_p(i_tor,1)

                  if ( jorek_model .eq. 400 ) then
                    Ti_x  = Ti_x   + (   Z_t * Ti_s - Z_s * Ti_t )   / xjac * HZ(i_tor,i_plane)
                    Ti_y  = Ti_y   + ( - R_t * Ti_s + R_s * Ti_t )   / xjac * HZ(i_tor,i_plane)
                    Ti_p  = Ti_p   + Ti * HZ_p(i_tor,1)
                    Te_x  = Te_x   + (   Z_t * Te_s - Z_s * Te_t )   / xjac * HZ(i_tor,i_plane)
                    Te_y  = Te_y   + ( - R_t * Te_s + R_s * Te_t )   / xjac * HZ(i_tor,i_plane)
                    Te_p  = Te_p   + Te * HZ_p(i_tor,1)
                  end if

                  RHO_x = RHO_x  + (   Z_t * RHO_s - Z_s * RHO_t ) / xjac * HZ(i_tor,i_plane)
                  RHO_y = RHO_y  + ( - R_t * RHO_s + R_s * RHO_t ) / xjac * HZ(i_tor,i_plane)
                  RHO_p = RHO_p  + RHO * HZ_p(i_tor,1)

                  w0_x  = w0_x   + (   Z_t * U_s - Z_s * U_t )     / xjac * HZ(i_tor,i_plane)
                  w0_y  = w0_y   + ( - R_t * U_s + R_s * U_t )     / xjac * HZ(i_tor,i_plane)

                  w0_xx = w0_xx  + (w_ss * Z_t**2 - 2.d0*w_st * Z_s*Z_t + w_tt * Z_s**2       &
                       + w_s * (Z_st*Z_t - Z_tt*Z_s )                                        & 
                       + w_t * (Z_st*Z_s - Z_ss*Z_t ) )     / xjac**2                        & 
                       - xjac_x * (w_s* Z_t - w_t * Z_s)  / xjac**2

                  w0_yy = w0_yy  + (w_ss * R_t**2 - 2.d0*w_st * R_s*R_t + w_tt * R_s**2       &
                       + w_s * (R_st*R_t - R_tt*R_s )                                        &
                       + w_t * (R_st*R_s - R_ss*R_t ) )         / xjac**2                    &
                       - xjac_y * (- w_s * R_t + w_t * R_s )  / xjac**2

               endif

            enddo  ! end loop toroidal harmonics

            ps0 = scalars(inode,1)
            psi_norm = (ps0 - psi_axis)/(psi_bnd - psi_axis)
            if ((psi_norm .lt. 1.d0) .and. (xpoint) .and. (Z .lt. Z_xpoint(1)) .and. (xcase .ne. 2)) then
               psi_norm = 2.d0 - psi_norm
            endif
            if ((psi_norm .lt. 1.d0) .and. (xpoint) .and. (Z .gt. Z_xpoint(2)) .and. (xcase .ne. 1)) then
               psi_norm = 2.d0 - psi_norm
            endif

            if ( jorek_model .eq. 400 ) then
              call bootstrap_current_rhs(BigR, minRad, R_axis, psi_axis, psi_bnd,  &
                                         psi, ps_x, ps_y, rho,  rho_x,  rho_y,     &
				         Ti,  Ti_x, Ti_y, Te,   Te_x,   Te_y,       Jb)
            else
	      Jb = 0.d0
            endif
 
            v_perp = R * sqrt(u0_x*u0_x + u0_y * u0_y)

            Btot = sqrt(F0**2 + ps_x**2 + ps_y**2) / BigR
            D_prof  = get_dperp (psi_norm)
            ZK_prof = get_zkperp(psi_norm)
            !!!small_r = sqrt((R -R_geo)*(R -R_geo)+(Z -Z_geo)*(Z -Z_geo))
            !small_r= amin*sqrt(psi_norm)
            !!!qprofile = small_r*(F0/BigR)/(R*Btheta)
            !qprofile = F0 / (R*Btheta)
            !scalars(inode,6) = max( scalars(inode,6), T_1 ) ! (workaround to avoid floating invalid error)
            ZKpar_T = ZK_par * ((max( scalars(inode,6), T_1 )+T_1)/T_0)**2.5

            grad_psi = sqrt(ps_x*ps_x + ps_y*ps_y)

            !   'E_flux_Kpar ','E_flux_kperp','E_flux_Vpar ','E_flux_Vperp','D_flux_Dperp','D_flux_Vpar ','D_flux_Vperp'/)

            scalars(inode,n_var+1)   = scalars(inode,5) * scalars(inode,6)

            scalars(inode,n_var+15) = Jb

            if (grad_psi .ne. 0.d0) then

               scalars(inode,n_var+2)  = ZKpar_T * ( F0 * TT_p / BigR**2  + (TT_x * ps_y - TT_y * ps_x) / BigR ) / Btot

               scalars(inode,n_var+3)  = ZK_prof * (TT_x * ps_x + TT_y * ps_y) / grad_psi

               scalars(inode,n_var+4)  = scalars(inode,5) * scalars(inode,6) * scalars(inode,7) * Btot

               scalars(inode,n_var+5)  = BigR   * (u0_x * ps_y - u0_y * ps_x) / sqrt(ps_x*ps_x + ps_y*ps_y) * scalars(inode,5) * scalars(inode,6)

               scalars(inode,n_var+6)  = D_prof * (RHO_x * ps_x + RHO_y * ps_y) / sqrt(ps_x*ps_x + ps_y*ps_y)

               scalars(inode,n_var+7)  = scalars(inode,5) * scalars(inode,7) * Btot

               scalars(inode,n_var+8)  = BigR   * (u0_x * ps_y - u0_y * ps_x) / sqrt(ps_x*ps_x + ps_y*ps_y) * scalars(inode,5)

               scalars(inode,n_var+17) = psi_norm

               !           call pellet_source(pellet_amplitude, pellet_R, pellet_Z, pellet_psi, pellet_phi, &
               !                        pellet_radius, pellet_delta_psi, pellet_sig, pellet_length,   &
               !                        R,Z,ps0,0.d0,source_pellet)
               !          scalars(inode,n_var+8) =source_pellet

            endif

            !        vectors(inode,:,1) = (/ - R * u0_y ,   + R * u0_x ,   0.d0 /)
            !        vectors(inode,:,2) = (/ + ps_y /R * scalars(inode,7), - ps_x /R * scalars(inode,7), 0.d0 /) * Btot
            !        vectors(inode,:,3) = (/ - R * u0_y + ps_y /R * scalars(inode,7) * Btot, + R * u0_x - ps_x /R * scalars(inode,7) * Btot, 0.d0 /)

         endif

      enddo
   enddo

   do j=1,nsub-1
      do k=1,nsub-1
         ielm        = ielm+1
         ien(1,ielm) = inode - nsub*nsub + nsub*(j-1) + k-1       ! 0 based indices for VTK
         ien(2,ielm) = inode - nsub*nsub + nsub*(j  ) + k-1
         ien(3,ielm) = inode - nsub*nsub + nsub*(j  ) + k
         ien(4,ielm) = inode - nsub*nsub + nsub*(j-1) + k
      enddo
   enddo

   ! endif
enddo
!===========================================================real values=============
rho_norm = central_density*1.d20 * 2.d0 * 1.67d-27
t_norm   = sqrt(MU_zero*rho_norm)
!=================================================real values============
 do i=1,nnos
!============================================j_phi in MA/m2
scalars(i,3) = scalars(i,3)/ MU_zero*1.e-6 
!============================================density in 1e20m-3
scalars(i,5) = scalars(i,5) * central_density
!===========================================temperature in keV
scalars(i,6) = scalars(i,6) / MU_zero / (central_density * 1d20) / EL_CHG /2./1.e3 !(assumes Te=Ti=T/2)
!=====================================Vparal in km/s *Btot!!!
scalars(i,7) = scalars(i,7)*scalars(i,n_var+14)/t_norm/1.e3
!=====================Pressure in kPa
scalars(i,8) = scalars(i,8) / MU_zero/1.e3
!============================Er in kV/m
scalars(i,16) = F0*scalars(i,16) / t_norm/1.e3
!====================================Vtheta km/s
scalars(i,17) = scalars(i,17) / t_norm/1.e3
!===================================Vsound in km/s
scalars(i,20) = scalars(i,20) / t_norm/1.e3
!===================================Vneo in km/s
scalars(i,23) = scalars(i,23) / t_norm/1.e3
!===================================Vperp_e in km/s
scalars(i,25) = scalars(i,25) / t_norm/1.e3
enddo

call Integrals_3D(my_id, node_list,element_list,density_tot,density_in,density_out,pressure,pressure_in,pressure_out)

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
