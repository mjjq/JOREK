subroutine grid_double_xpoint(node_list, element_list, n_flux, n_tht, n_open, n_outer, n_inner, n_private, n_up_priv, n_leg, n_up_leg, &
                              SIG_closed, SIG_theta, SIG_open, SIG_outer, SIG_inner, SIG_private, SIG_up_priv,                         &
		              SIG_leg_0, SIG_leg_1, SIG_up_leg_0, SIG_up_leg_1, dPSI_open, dPSI_outer, dPSI_inner, dPSI_private, dPSI_up_priv, xcase)
!-----------------------------------------------------------------------
! subroutine defines a flux surface aligned finite element grid
! inclduing a single x-point
!-----------------------------------------------------------------------

use constants
use tr_module 
use data_structure
use grid_xpoint_data

implicit none

! --- Routine parameters
type (type_node_list),    intent(inout) :: node_list
type (type_element_list), intent(inout) :: element_list
integer,                  intent(in)    :: n_flux, n_open, xcase
integer,                  intent(inout) :: n_tht, n_outer, n_inner, n_private, n_leg, n_up_priv, n_up_leg  
real*8,                   intent(in)    :: SIG_closed, SIG_theta, SIG_open, SIG_outer, SIG_inner, SIG_private, SIG_up_priv
real*8,                   intent(in)    :: SIG_leg_0, SIG_leg_1, SIG_up_leg_0, SIG_up_leg_1
real*8,                   intent(in)    :: dPSI_open, dPSI_outer, dPSI_inner, dPSI_private, dPSI_up_priv

! --- local variables
type (type_surface_list) :: flux_list, sep_list

!type (type_node_list),    pointer :: newnode_list
!type (type_element_list), pointer :: newelement_list
type (type_strategic_points) , pointer     :: stpts
type (type_new_points)       , pointer     :: nwpts

real*8              :: psi_axis, R_axis, Z_axis, s_axis, t_axis, R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2), psi_xpoint(2)
real*8              :: s_find(8), t_find(8), tht_x, theta, delta, tmp1, tmp2
real*8              :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8              :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8              :: PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss
!real*8,allocatable  :: xout(:),xp(:),yp(:)
real*8              :: dR_dt, dZ_dt, RZ_jac, PSI_R, PSI_Z
!integer,allocatable :: keep(:,:,:)
integer             :: i, j, k, l, m, i2, j2
integer             :: n_psi
integer             :: i_surf
integer             :: i_elm_axis, i_elm_xpoint(2), i_elm_find(8), i_sep, i_max, i_find, npl, ifail
integer             :: node, index, node_start, index_xpoint, n_xpoint, j_start, j_end
integer             :: iv, ivp, node_iv, node_ivp, i_elm
integer             :: my_id, ielm_out
real*8              :: Rmid, Zmid, R0,Z0, RP,ZP, dR0, dZ0, dRP, dZP, size_0, size_p, denom
real*8              :: R1, Z1, s_out, t_out, R_out, Z_out
real*8              :: EJAC, RX, RY, SX, SY, CRR, CZZ, CRZ, alpha1, alpha2, alpha_max, alpha_min
real*8              :: rr, ss, drr, dss, tt
real*8              :: rr1, ss1, drr1, dss1
real*8              :: rr2, ss2, drr2, dss2
real*8              :: psi_bnd, psi_bnd2
real*8              :: sigmas(16)
integer             :: n_grids(10)
logical             :: xpoint
!real*8,external     :: root
character*4         :: label

xpoint = .true.
my_id  = 1 ! Just don't want the printout...

write(*,*) ' '
write(*,*) ' '
write(*,*) '*************************************'
write(*,*) '*************************************'
write(*,*) '*************************************'
write(*,*) '*          X-point grid             *'
write(*,*) '*************************************'
write(*,*) '*************************************'
write(*,*) '*************************************'
write(*,*) ' '
write(*,*) ' '


allocate(stpts)
allocate(nwpts)
call tr_register_mem(sizeof(stpts),"stpts",CAT_GRID)
call tr_register_mem(sizeof(nwpts),"nwpts",CAT_GRID)
!-------------------------------------------------------------------------------------------!
!---------------------------- Initialise some internal data --------------------------------!
!-------------------------------------------------------------------------------------------!

!-------------------------------- Reset some parameters if they are inconsistent with XCASE
if (xcase .eq. 1) then
  n_outer   = 0
  n_inner   = 0
  n_up_priv = 0
  n_up_leg  = 0
endif
if (xcase .eq. 2) then
  n_outer   = 0
  n_inner   = 0
  n_private = 0
  n_leg     = 0
endif
if ( (xcase .eq. 3) .and. (mod(n_tht,2) .ne. 0) )  n_tht = n_tht + 1
if ( (xcase .ne. 3) .and. (mod(n_tht,2) .eq. 0) )  n_tht = n_tht + 1

!-------------------------------- Build up some arrays to send as routine parameters (avoid long lists...)
sigmas(1)  = SIG_closed  ; sigmas(2)  = SIG_theta
sigmas(3)  = SIG_open    ; sigmas(4)  = SIG_outer   ; sigmas(5)  = SIG_inner
sigmas(6)  = SIG_private ; sigmas(7)  = SIG_up_priv
sigmas(8)  = SIG_leg_0   ; sigmas(9)  = SIG_leg_1
sigmas(10) = SIG_up_leg_0; sigmas(11) = SIG_up_leg_1
sigmas(12) = dPSI_open   ; sigmas(13) = dPSI_outer  ; sigmas(14) = dPSI_inner
sigmas(15) = dPSI_private; sigmas(16) = dPSI_up_priv

n_grids(1) = n_flux   ; n_grids(2) = n_tht
n_grids(3) = n_open   ; n_grids(4) = n_outer  ; n_grids(5) = n_inner
n_grids(6) = n_private; n_grids(7) = n_up_priv
n_grids(8) = n_leg    ; n_grids(9) = n_up_leg
n_grids(10)= 0 ! keep for n_tht_outer, which determines the angle of second Xpoint

!-------------------------------------------------------------------------------------------!
!----------------------------- Find MagAxis and Xpoint -------------------------------------!
!-------------------------------------------------------------------------------------------!

call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

psi_bnd  = 0.d0
psi_bnd2 = 0.d0
call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
if(xcase .eq. 1) psi_bnd = psi_xpoint(1)
if(xcase .eq. 2) psi_bnd = psi_xpoint(2)
if(xcase .eq. 3) then
  if(psi_xpoint(2) .lt. psi_xpoint(1)) then
    psi_bnd  = psi_xpoint(2)
    psi_bnd2 = psi_xpoint(1)
  else
    psi_bnd  = psi_xpoint(1)
    psi_bnd2 = psi_xpoint(2)  
  endif
  ! if we have a symmetric double-null, force the single separatrix
  if (abs(psi_xpoint(1)-psi_xpoint(2)) .lt. 1.d-4) then
    psi_xpoint(1) = (psi_xpoint(1)+psi_xpoint(2))/2.d0
    psi_xpoint(2) = psi_xpoint(1)
    psi_bnd  = psi_xpoint(1)
    psi_bnd2 = psi_bnd  
  endif
endif




!-------------------------------------------------------------------------------------------!
!--------------- Define the flux values on which grid will be aligned ----------------------!
!-------------------------------------------------------------------------------------------!

!-------------------------------- Write some values
!write(*,*) ' n_flux,   n_open,   n_tht   : ', n_flux,	n_open,   n_tht
!if(xcase .eq. 3) then
!  write(*,*) ' n_outer,   n_inner   : ', n_outer,   n_inner
!endif

!-------------------------------- Define number of psi values and allocate flux_list structure
n_psi	        = n_flux   + n_open   + n_outer   + n_inner   + n_private   + n_up_priv + 1   ! this includes the magnetic axis
flux_list%n_psi = n_psi - 1
call tr_allocate(flux_list%psi_values,1,flux_list%n_psi,"flux_list%psi_values",CAT_GRID)

!-------------------------------- Allocate sep_list structure (that's for plotting only)
sep_list%n_psi =3
if(xcase .eq. 3) sep_list%n_psi =6
if (allocated(sep_list%psi_values)) call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
call tr_allocate(sep_list%psi_values,1,sep_list%n_psi,"sep_list%psi_values",CAT_GRID)

!-------------------------------- Call the routine
call define_flux_values(node_list, element_list, flux_list, sep_list, &
                        xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_axis, n_grids, sigmas)

call plot_flux_surfaces(node_list,element_list,flux_list,.true.,1,psi_xpoint,R_xpoint,Z_xpoint,.true.,xcase)
call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,psi_xpoint,R_xpoint,Z_xpoint,.true.,xcase)

if (allocated(sep_list%flux_surfaces))     deallocate(sep_list%flux_surfaces)





!-------------------------------------------------------------------------------------------!
!-------- Find all strategic points (leg corners, strike points and private middles) -------!
!-------------------------------------------------------------------------------------------!

!-------------------------------- Call the routine
call find_strategic_points(node_list, element_list, flux_list, &
                           xcase, R_xpoint, Z_xpoint, psi_xpoint, R_axis, Z_axis, n_grids, stpts)


!-------------------------------------------------------------------------------------------!
!-------------- Find new grid points by crossing polar and radial coordinates --------------!
!-------------------------------------------------------------------------------------------!

!-------------------------------- Call the routine
call define_new_grid_points(node_list, element_list, flux_list, &
                             xcase, R_xpoint, Z_xpoint, psi_xpoint, n_grids, stpts, sigmas, nwpts)



!-------------------------------------------------------------------------------------------!
!---------------------- Define the final grid (new nodes and new elements) -----------------!
!-------------------------------------------------------------------------------------------!

!-------------------------------- Call the routine
call define_final_grid(node_list, element_list, flux_list, &
		       xcase, n_grids, stpts, nwpts)


deallocate(stpts)
deallocate(nwpts)
call tr_unregister_mem(sizeof(stpts),"stpts",CAT_GRID)
call tr_unregister_mem(sizeof(nwpts),"nwpts",CAT_GRID)

return
end subroutine grid_double_xpoint
