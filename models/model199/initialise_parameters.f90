subroutine initialise_parameters(my_id)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use phys_module

implicit none
integer :: in, my_id, i
real*8  :: psi_plot(1001),zj_plot(1001),dj_plot(1001),dz_plot(1001), z_plot(1001), z, zjz, dj_dpsi, dj_dz, psi_n

namelist /in1/  tstep, nstep, eta, visco, restart,  regrid,         &
                n_R, n_Z, n_radial, n_pol, n_tht, n_flux,           &
                n_open,n_private,n_leg,  nout,                      &
                xr1, sig1, xr2, sig2,                               &
                R_begin, R_end, Z_begin, Z_end,                     &
                R_geo, Z_geo, amin, mf, fbnd, fpsi, mode,           &
                F0,                                                 &
                zjz_0, zjz_1, zj_coef,                              &
                rho_0, rho_1, rho_coef,                             &
                T_0,   T_1,   T_coef,                               &
                FF_0,  FF_1,  FF_coef,                              &
                ZK_par, ZK_perp, D_par, D_perp,                     &
                particlesource, heatsource,                         &
                eta_num, visco_num,                                 &
                ellip,tria_u,tria_l,quad_u,quad_l,                  &
                xampl,xwidth,xsig,xtheta,xshift,xleft, xpoint,      &
		        freeboundary,refinement

if (my_id .eq. 0) then

  tstep    = 1.d0
  tstep_in = 1.d0
  nstep    = 0

  eta   = 1.d-5
  visco = 1.d-5

  restart      = .false.
  import_equil = .false.
  regrid       = .false.
  
  freeboundary = .false.
  wall_type    = 'ideal'

  n_R       = 3
  n_Z       = 3
  n_radial  = 11
  n_pol     = 6

  n_flux    = 11
  n_tht     = 17

  n_open    = 5
  n_leg     = 5
  n_private = 5

  R_geo     = 10.d0
  Z_geo     = 0.d0
  amin      = 1.d0

  F0        = 10.d0

  mf        = 8
  fbnd      = 0.d0; fbnd(1) =2.d0

  ellip  = 1.d0
  tria_u = 0.d0
  tria_l = 0.d0
  quad_l = 0.d0
  quad_u = 0.d0

  xampl  = 0.d0
  xwidth = 0.d0
  xsig   = 0.d0
  xtheta = 0.d0
  xshift = 0.d0
  xleft  = 0.d0
  xpoint = .false.

  xr1  = 9999.d0
  sig1 = 9999.d0
  xr2  = 99999.d0
  sig2 = 99999.d0

  R_begin = -0.1d0
  R_end   =  0.1d0
  Z_begin = -0.1d0
  Z_end   = 0.1d0

  ZK_perp(1) = 1.d-5; ZK_perp(2) = 0.d0; ZK_perp(3)= 0.d0; ZK_perp(4)= 99.d0; ZK_perp(5) = 99.d0
  ZK_par     = 1.d0
  D_perp(1)  = 1.d-5; D_perp(2) = 0.d0; D_perp(3)= 0.d0; D_perp(4)= 99.d0; D_perp(5) = 99.d0
  D_par      = 0.d0

  eta_num   = 0.d0
  visco_num = 0.d0

  heatsource     = 1.e-7
  particlesource = 1.e-5

  zjz_0 =  0.1173d0;   T_0   =  1.d-6  ;   rho_0 =  1.d0   ;   FF_0  =  1.d0
  zjz_1 =  0.0d0   ;   T_1   =  1.d-8  ;   rho_1 =  1.d0   ;   FF_1  =  0.d0

  zj_coef     = 0.d0;  zj_coef(1)  = -1.d0
  T_coef      = 0.d0;  T_coef(1)   = -1.d0
  rho_coef    = 0.d0;  rho_coef(1) =  0.d0
  FF_coef     = 0.d0;  FF_coef(1)  = -1.d0

  t_now       = 0.d0
  t_start     = 0.d0
  index_start = 0

  nout = 9999999

  if (my_id .eq. 0) read(5,in1)

  do in=1, n_tor
    mode(in) = + int(in / 2) * n_period
    write(*,*) ' toroidal mode numbers : ',in,mode(in)
  enddo


  tstep_in = tstep

  if (nstep .gt. 0) allocate(energies(n_tor,2,nstep))
  if (nstep .gt. 0) allocate(xtime(nstep))

endif

call broadcast_phys(my_id)         ! distribute some values

if (my_id .eq. 0) then
  230 FORMAT(A, ' = ' 10ES12.4)
  231 FORMAT(A, ' = ' 10I)
  232 FORMAT(A, ' = ' 10L)
  233 FORMAT(A, ' = ' 4ES12.4, '     ...    ', 4ES12.4)
  234 FORMAT(A, ' = ' 9ES12.4, '     ...')
  235 FORMAT(A, ' = "' A, '"')
  write(*,*)   'USING MODEL 199'
  write(*,*)
  write(*,231) 'n_var         ', n_var
  write(*,231) 'n_dim         ', n_dim
  write(*,231) 'n_order       ', n_order
  write(*,231) 'n_tor         ', n_tor
  write(*,231) 'n_period      ', n_period
  write(*,231) 'n_plane       ', n_plane
  write(*,231) 'n_vertex_max  ', n_vertex_max
  write(*,231) 'n_nodes_max   ', n_nodes_max
  write(*,231) 'n_elements_max', n_elements_max
  write(*,231) 'n_boundary_max', n_boundary_max
  write(*,231) 'n_pieces_max  ', n_pieces_max
  write(*,231) 'n_degrees     ', n_degrees
  write(*,*)
  write(*,230) 'tstep         ', tstep
  write(*,231) 'nstep         ', nstep
  write(*,230) 'eta           ', eta
  write(*,230) 'visco         ', visco
  write(*,232) 'restart       ', restart
  write(*,232) 'regrid        ', regrid
  write(*,231) 'n_R           ', n_R
  write(*,231) 'n_Z           ', n_Z
  write(*,231) 'n_radial      ', n_radial
  write(*,231) 'n_pol         ', n_pol
  write(*,231) 'n_tht         ', n_tht
  write(*,231) 'n_flux        ', n_flux
  write(*,231) 'n_open        ', n_open
  write(*,231) 'n_private     ', n_private
  write(*,231) 'n_leg         ', n_leg
  write(*,231) 'nout          ', nout
  write(*,230) 'xr1           ', xr1
  write(*,230) 'sig1          ', sig1
  write(*,230) 'xr2           ', xr2
  write(*,230) 'sig2          ', sig2
  write(*,230) 'R_begin       ', R_begin
  write(*,230) 'R_end         ', R_end
  write(*,230) 'Z_begin       ', Z_begin
  write(*,230) 'Z_end         ', Z_end
  write(*,230) 'R_geo         ', R_geo
  write(*,230) 'Z_geo         ', Z_geo
  write(*,230) 'amin          ', amin
  write(*,231) 'mf            ', mf
  if ( mf >= 0 ) then
    write(*,234) 'fbnd          ', fbnd(1:MIN(9,mf))
    write(*,234) 'fpsi          ', fpsi(1:MIN(9,mf))
  end if
  write(*,231) 'mode          ', mode
  write(*,230) 'F0            ', F0
  write(*,230) 'zjz_0         ', zjz_0
  write(*,230) 'zjz_1         ', zjz_1
  write(*,230) 'zj_coef       ', zj_coef
  write(*,230) 'rho_0         ', rho_0
  write(*,230) 'rho_1         ', rho_1
  write(*,230) 'rho_coef      ', rho_coef(1:5)
  write(*,230) 'T_0           ', T_0
  write(*,230) 'T_1           ', T_1
  write(*,230) 'T_coef        ', T_coef(1:5)
  write(*,230) 'FF_0          ', FF_0
  write(*,230) 'FF_1          ', FF_1
  write(*,230) 'FF_coef       ', FF_coef(1:8)
  write(*,230) 'ZK_par        ', ZK_par
  write(*,230) 'ZK_perp       ', ZK_perp(1:5)
  write(*,230) 'D_par         ', D_par
  write(*,230) 'D_perp        ', D_perp(1:5)
  write(*,230) 'particlesource', particlesource
  write(*,230) 'heatsource    ', heatsource
  write(*,230) 'eta_num       ', eta_num
  write(*,230) 'visco_num     ', visco_num
  write(*,230) 'ellip         ', ellip
  write(*,230) 'tria_u        ', tria_u
  write(*,230) 'tria_l        ', tria_l
  write(*,230) 'quad_u        ', quad_u
  write(*,230) 'quad_l        ', quad_l
  write(*,230) 'xampl         ', xampl
  write(*,230) 'xwidth        ', xwidth
  write(*,230) 'xsig          ', xsig
  write(*,230) 'xtheta        ', xtheta
  write(*,230) 'xshift        ', xshift
  write(*,230) 'xleft         ', xleft
  write(*,232) 'xpoint        ', xpoint
  write(*,231) 'n_boundary    ', n_boundary
  if ( n_boundary > 0 ) then
    write(*,233) 'r_boundary    ', r_boundary(1:4), r_boundary(n_boundary-3:n_boundary)
    write(*,233) 'z_boundary    ', z_boundary(1:4), z_boundary(n_boundary-3:n_boundary)
    write(*,233) 'psi_boundary  ', psi_boundary(1:4), psi_boundary(n_boundary-3:n_boundary)
  end if
  write(*,232) 'freeboundary  ', freeboundary
  write(*,235) 'wall_type     ', TRIM(wall_type)
  write(*,232) 'refinement    ', refinement
endif


return
end
