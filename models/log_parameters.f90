subroutine log_parameters(my_id)
! Write out all relevant parameters from mod_parameters and the input file.

use phys_module

implicit none

integer, INTENT(IN) :: my_id



if (my_id .eq. 0) then
  230 FORMAT(A, ' = ' 10ES12.4)
  231 FORMAT(A, ' = ' 10I)
  232 FORMAT(A, ' = ' 10L)
  233 FORMAT(A, ' = ' 4ES12.4, '     ...    ', 4ES12.4)
  234 FORMAT(A, ' = ' 9ES12.4, '     ...')
  write(*,*)
  write(*,*) 'PARAMETERS OF THE JOREK RUN:'
  write(*,*) '----------------------------'
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
  if ( freeboundary ) then
    write(*,232) 'use_starwall  ', use_starwall
    write(*,232) 'resistive_wall', resistive_wall
  end if
  write(*,232) 'refinement    ', refinement
  write(*,*)
endif

return
end subroutine log_parameters
