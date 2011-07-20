!> Write out all relevant parameters defined in mod_parameters
!! and by the namelist input file.
subroutine log_parameters(my_id)

use phys_module

implicit none

! --- Routine parameters
integer, intent(in) :: my_id !< MPI proc id

! --- Local variables
integer           :: ivar, itor
character(len=10) :: mode_num

if (my_id == 0) then
  
  230 format(1X,A, ' = ', 10ES12.4)
  231 format(1X,A, ' = ', 10I12)
  232 format(1X,A, ' = ', 10L12)
  233 format(1X,A, ' = ', 4ES12.4, '     ...    ', 4ES12.4)
  234 format(1X,A, ' = ', 9ES12.4, '     ...')
  235 format(3x,I3,': ',A)
  236 format(3x,I3,': ',A,'(',A,'*phi)')
  237 format(1X,A, ' = "', A, '"')
  
  write(*,*)
  write(*,*) '*************************************************'
  write(*,*) '*          PREPROCESSOR OPTIONS                 *'
  write(*,*) '*************************************************'
  write(*,'(1x,a)',advance='no') ' USE_MUMPS                : '
#ifdef USE_MUMPS
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif
  
  write(*,'(1x,a)',advance='no') ' USE_PASTIX               : '
#ifdef USE_PASTIX
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif
  
  write(*,'(1x,a)',advance='no') ' USE_MURGE                : '
#ifdef USE_MURGE
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif
  
  write(*,'(1x,a)',advance='no') ' USE_HIPS                 : '
#ifdef USE_HIPS
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif
  
  write(*,'(1x,a)',advance='no') ' USE_BLOCK                : '
#ifdef USE_BLOCK
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif
  
  write(*,'(1x,a)',advance='no') ' MEMTRACE                 : '
#ifdef MEMTRACE 
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,*)
  write(*,*) '*************************************************'
  write(*,*) '*          PARAMETERS OF THE JOREK RUN          *'
  write(*,*) '*************************************************'
  write(*,231) 'jorek_model       ', jorek_model
  write(*,231) 'n_var             ', n_var
  do ivar = 1, n_var
    write(*,235) ivar, trim(variable_names(ivar))
  end do
  write(*,231) 'n_dim             ', n_dim
  write(*,231) 'n_order           ', n_order
  write(*,231) 'n_tor             ', n_tor
  write(*,236) 1, 'cos', '0'
  do itor = 2, n_tor
    write(mode_num,'(I4)') mode(itor)
    write(*,236) itor, mode_type(itor), trim(adjustl(mode_num))
  end do
  write(*,231) 'n_period          ', n_period
  write(*,231) 'n_plane           ', n_plane
  write(*,231) 'n_vertex_max      ', n_vertex_max
  write(*,231) 'n_nodes_max       ', n_nodes_max
  write(*,231) 'n_elements_max    ', n_elements_max
  write(*,231) 'n_boundary_max    ', n_boundary_max
  write(*,231) 'n_pieces_max      ', n_pieces_max
  write(*,231) 'n_degrees         ', n_degrees
  write(*,*)
  write(*,231) 'n_period        ', n_period
  write(*,231) 'n_plane         ', n_plane
  write(*,231) 'n_vertex_max    ', n_vertex_max
  write(*,231) 'n_nodes_max     ', n_nodes_max
  write(*,231) 'n_elements_max  ', n_elements_max
  write(*,231) 'n_boundary_max  ', n_boundary_max
  write(*,231) 'n_pieces_max    ', n_pieces_max
  write(*,231) 'n_degrees       ', n_degrees
  write(*,*)
  write(*,230) 'tstep             ', tstep
  write(*,231) 'nstep             ', nstep
  write(*,230) 'tstep_n           ', tstep_n
  write(*,231) 'nstep_n           ', nstep_n
  write(*,230) 'eta               ', eta
  write(*,230) 'visco             ', visco
  write(*,230) 'visco_par         ', visco_par
  write(*,232) 'restart           ', restart
  write(*,232) 'regrid            ', regrid
  write(*,231) 'n_R               ', n_R
  write(*,231) 'n_Z               ', n_Z
  write(*,231) 'n_radial          ', n_radial
  write(*,231) 'n_pol             ', n_pol
  write(*,231) 'n_tht             ', n_tht
  write(*,231) 'n_flux            ', n_flux
  write(*,232) 'xpoint            ', xpoint
  
  if ( xpoint ) then
    write(*,231) 'n_open            ', n_open
    write(*,231) 'n_private         ', n_private
    write(*,231) 'n_leg             ', n_leg
    write(*,230) 'SIG_closed        ', SIG_closed
    write(*,230) 'SIG_open          ', SIG_open
    write(*,230) 'SIG_private       ', SIG_private
    write(*,230) 'SIG_theta         ', SIG_theta
    write(*,230) 'SIG_leg_0         ', SIG_leg_0
    write(*,230) 'SIG_leg_1         ', SIG_leg_1
    write(*,230) 'dPSI_open         ', dPSI_open
    write(*,230) 'dPSI_private      ', dPSI_private
  end if
  
  write(*,231) 'nout              ', nout
  write(*,230) 'xr1               ', xr1
  write(*,230) 'sig1              ', sig1
  write(*,230) 'xr2               ', xr2
  write(*,230) 'sig2              ', sig2
  write(*,230) 'R_begin           ', R_begin
  write(*,230) 'R_end             ', R_end
  write(*,230) 'Z_begin           ', Z_begin
  write(*,230) 'Z_end             ', Z_end
  write(*,230) 'R_geo             ', R_geo
  write(*,230) 'Z_geo             ', Z_geo
  write(*,230) 'amin              ', amin
  write(*,231) 'mf                ', mf
  
  if ( mf >= 0 ) then
    write(*,234) 'fbnd              ', fbnd(1:MIN(9,mf))
    write(*,234) 'fpsi              ', fpsi(1:MIN(9,mf))
  end if
  
  write(*,230) 'F0                ', F0
  write(*,230) 'zjz_0             ', zjz_0
  write(*,230) 'zjz_1             ', zjz_1
  write(*,230) 'zj_coef           ', zj_coef
  
  if ( .not. num_rho ) then
    write(*,230) 'rho_0             ', rho_0
    write(*,230) 'rho_1             ', rho_1
    write(*,230) 'rho_coef          ', rho_coef(1:5)
  else
    write(*,237) 'rho_file          ', trim(rho_file)
  end if
  
  if ( .not. num_T ) then
    write(*,230) 'T_0               ', T_0
    write(*,230) 'T_1               ', T_1
    write(*,230) 'T_coef            ', T_coef(1:5)
  else
    write(*,237) 'T_file            ', trim(T_file)
  end if
  
  if ( .not. num_ffprime ) then
    write(*,230) 'FF_0              ', FF_0
    write(*,230) 'FF_1              ', FF_1
    write(*,230) 'FF_coef           ', FF_coef(1:8)
  else
    write(*,237) 'ffprime_file      ', trim(ffprime_file)
  end if
  
  write(*,230) 'ZK_par            ', ZK_par
  write(*,230) 'ZK_perp           ', ZK_perp(1:5)
  write(*,230) 'D_par             ', D_par
  write(*,230) 'D_perp            ', D_perp(1:5)
  write(*,230) 'particlesource    ', particlesource
  write(*,230) 'heatsource        ', heatsource
  write(*,230) 'tauIC             ', tauIC
  write(*,230) 'eta_num           ', eta_num
  write(*,230) 'visco_num         ', visco_num
  write(*,230) 'visco_par_num     ', visco_par_num
  write(*,230) 'D_perp_num        ', D_perp_num
  write(*,230) 'ZK_perp_num       ', ZK_perp_num
  write(*,230) 'pellet_amplitude  ', pellet_amplitude
  write(*,230) 'pellet_R          ', pellet_R
  write(*,230) 'pellet_Z          ', pellet_Z
  write(*,230) 'pellet_phi        ', pellet_phi
  write(*,230) 'pellet_radius     ', pellet_radius
  write(*,230) 'pellet_sig        ', pellet_sig
  write(*,230) 'pellet_length     ', pellet_length
  write(*,230) 'pellet_psi        ', pellet_psi
  write(*,230) 'pellet_delta_psi  ', pellet_delta_psi
  write(*,*)
  write(*,230) 'ellip             ', ellip
  write(*,230) 'tria_u            ', tria_u
  write(*,230) 'tria_l            ', tria_l
  write(*,230) 'quad_u            ', quad_u
  write(*,230) 'quad_l            ', quad_l
  write(*,230) 'xampl             ', xampl
  write(*,230) 'xwidth            ', xwidth
  write(*,230) 'xsig              ', xsig
  write(*,230) 'xtheta            ', xtheta
  write(*,230) 'xshift            ', xshift
  write(*,230) 'xleft             ', xleft
  write(*,231) 'n_boundary        ', n_boundary
  
  if ( n_boundary > 0 ) then
    write(*,233) 'r_boundary        ', r_boundary(1:4), r_boundary(n_boundary-3:n_boundary)
    write(*,233) 'z_boundary        ', z_boundary(1:4), z_boundary(n_boundary-3:n_boundary)
    write(*,233) 'psi_boundary      ', psi_boundary(1:4), psi_boundary(n_boundary-3:n_boundary)
  end if
  
  write(*,232) 'freeboundary_equil', freeboundary_equil
  write(*,232) 'freeboundary      ', freeboundary
  if ( freeboundary ) then
    write(*,232) 'use_starwall      ', use_starwall
    write(*,232) 'resistive_wall    ', resistive_wall
  end if
  
  write(*,232) 'refinement        ', refinement
  write(*,*)
  
end if

end subroutine log_parameters
