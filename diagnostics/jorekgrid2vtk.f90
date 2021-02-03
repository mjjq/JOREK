!> Program to convert a JOREK2 restart file into binary VTK format
program jorekgrid2vtk

use mod_parameters, only: n_var, variable_names
use data_structure
use phys_module
use basis_at_gaussian
use mpi_mod
use mod_import_restart
use mod_boundary
use mod_vtk
use mod_interp
implicit none

type (type_node_list)   ,     pointer :: node_list
type (type_element_list),     pointer :: element_list
type (type_bnd_element_list), pointer :: bnd_elm_list    
type (type_bnd_node_list),    pointer :: bnd_node_list 

integer               :: nnoel, nnos, nel, nsub, inode, ielm, n_scalars, n_vectors
real*4,allocatable    :: currdens(:), xyz (:,:), scalars(:,:), vectors(:,:,:)
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
real*8                :: dPsi, dPs_s, dPs_t, dPs_st, dPs_ss, dPs_tt
real*8                :: dU,    dU_s,  dU_t,  dU_st,  dU_ss,  dU_tt
real*8                :: ps0_x, ps0_y, psi_sum, ps_x, ps_y, ps_p
real*8                :: u0_x,  u0_y,  u_sum,   u_x,  u_y,  u_p
real*8                :: zj0_x, zj0_y, zj_sum,  zj_x, zj_y, zj_p
real*8                :: w0_x,  w0_y,  w_sum,   w0_xx, w0_yy, w_x, w_y, w_p, w_xx, w_yy
real*8                :: zn0_x, zn0_y, zn_sum,  zn_x, zn_y, zn_p, rho_x, rho_y, rho_p
real*8                :: T0_x,  T0_y,  T_sum,   TT_x, TT_y, TT_p
real*8                :: Ti0_x, Ti0_y, Ti_sum,  Ti_x, Ti_y, Ti_p
real*8                :: Te0_x, Te0_y, Te_sum,  Te_x, Te_y, Te_p
real*8                :: psi_axis,      R_axis,      Z_axis,      s_axis,      t_axis
real*8                :: psi_xpoint(2), R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2)
real*8                :: xjac, xjac_x, xjac_y, v_perp, Psi_J, R_p, error, Btot, BigR
integer               :: i_elm_axis, i_elm_xpoint(2), k_tor, ifail, ierr
real*8                :: toroidal_angle





namelist /vtk_params/ nsub, i_tor, i_plane


write(*,*) '***************************************'
write(*,*) '*        grid2vtk                     *'
write(*,*) '***************************************'

call flush_it(6)

allocate(node_list)
allocate(element_list)
allocate(bnd_elm_list)
allocate(bnd_node_list)

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Preset parameters
nsub                   = 5       ! Number of subdivisions of the cubic finite elements into linear pieces
i_tor                  = -1      ! If i_tor > 0, only this mode will be included in the vtk file...
i_plane                = 1       ! ... otherwise, all modes will be summed up at the toroidal plane i_plane


! --- Read parameters from namelist file 'vtk.nml' if it exists
!open(42, file='vtk.nml', action='read', status='old', iostat=ierr)
!if ( ierr == 0 ) then
!  write(*,*) 'Reading parameters from vtk.nml namelist.'
!  read(42,vtk_params)
!  close(42)
!end if

write(*,*)
write(*,*) 'Parameters:'
write(*,*) '-----------'
write(*,*) 'nsub            =', nsub
write(*,*) 'i_tor           =', i_tor
write(*,*) 'i_plane         =', i_plane


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

allocate(scalar_names(n_scalars), vector_names(n_vectors))


scalar_names(1:n_var) = variable_names(1:n_var)

do k_tor=1, n_tor
  mode(k_tor) = + int(k_tor / 2) * n_period
enddo

call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)

call initialise_basis                              ! define the basis functions at the Gaussian points

nnos = nsub*nsub*element_list%n_elements
allocate(currdens(nnos),xyz(3,nnos),scalars(nnos,1:n_scalars),vectors(nnos,3,1:n_vectors))
currdens = 0.

nnoel = 4
nel   = (nsub-1)*(nsub-1)*element_list%n_elements
allocate(ien(nnoel,nel))

inode   = 0
ielm    = 0
scalars = 0.d0
vectors = 0.d0
xyz     = 0
ien     = 0

call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

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



if (i.eq.133) then
write(*,'(A,2f)')'compare j at 0:',element_list%element(i)%size(2,6)*node_list%node(element_list%element(i)%vertex(2))%x(1,6,1:2)
write(*,'(A,2f)')'compare j at 1:',element_list%element(i)%size(3,6)*node_list%node(element_list%element(i)%vertex(3))%x(1,6,1:2)
s=1 ; t=0
call interp_RZ(node_list,element_list,i,s,t,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
!write(*,'(A,4f)')'checking_t  0  :',R_t,node_list%node(element_list%element(i)%vertex(2))%x(1,3,1),element_list%element(i)%size(2,3),0
!write(*,'(A,4f)')'checking_t  0  :',R_t, Z_t, node_list%node(element_list%element(i)%vertex(2))%x(1,3,1:2)
write(*,'(A,4f)')'checking_tt 0  :',R_tt,Z_tt,node_list%node(element_list%element(i)%vertex(2))%x(1,6,1:2)
s=1 ; t=1
call interp_RZ(node_list,element_list,i,s,t,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
!write(*,'(A,4f)')'checking_t  1  :',R_t, Z_t, node_list%node(element_list%element(i)%vertex(3))%x(1,3,1:2)
write(*,'(A,4f)')'checking_tt 1  :',R_tt,Z_tt,node_list%node(element_list%element(i)%vertex(3))%x(1,6,1:2)
s=1 ; t=0.5
call interp_RZ(node_list,element_list,i,s,t,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
!write(*,'(A,2f)')'checking_t  0.5:',R_t, Z_t
write(*,'(A,2f)')'checking_tt 0.5:',R_tt,Z_tt
endif

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
      
      xyz(1:3,inode) = (/ R,    Z, 0.d0 /)
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

      endif         ! xjac

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


!--------------------------------------------------- write the binary VTK file
etype = 9  ! for vtk_quad

call write_vtk('jorek_tmp.vtk',xyz,ien,etype,scalar_names,scalars,vector_names,vectors)

write(*,*) 'done.'

end program jorekgrid2vtk
