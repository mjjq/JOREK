!**********************************************************************
!* program to extract data from a JOREK2 restart file                 *
!**********************************************************************

program jorek2_diagno
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
use basis_at_gaussian
implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
integer :: i, in, i_tor
real*8  :: growth_kin, growth_mag,density,density_in,density_out,pressure,pressure_in,pressure_out
real*8  :: Rplot(2), Zplot(2)
real*8  :: psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis
integer :: ifail, my_id

namelist /in1/  tstep, nstep, tstep_n, nstep_n,                     &
                eta, visco, visco_par,                              &
                restart, regrid,                                    &
                n_R, n_Z, n_radial, n_pol, n_tht, n_flux,           &
                n_open, n_private, n_leg,                           &
                SIG_closed, SIG_open, SIG_private, SIG_theta,       &
                SIG_leg_0, SIG_leg_1, dPSI_open, dPSI_private,      &
                nout, xr1, sig1, xr2, sig2,                         &
                R_begin, R_end, Z_begin, Z_end,                     &
                R_geo, Z_geo, amin, mf, fbnd, fpsi, mode,           &
                R_boundary, Z_boundary, psi_boundary, n_boundary,   &
                F0,                                                 &
                zjz_0, zjz_1, zj_coef,                              &
                rho_0, rho_1, rho_coef,                             &
                T_0,   T_1,   T_coef,                               &
                FF_0,  FF_1,  FF_coef,                              &
                ZK_par, ZK_perp, D_par, D_perp,                     &
                particlesource, heatsource, tauIC,                  &
                eta_num, visco_num, visco_par_num, D_perp_num,      &
                pellet_amplitude, pellet_R, pellet_Z, pellet_phi,   &
                pellet_radius, pellet_sig, pellet_length,           &
                pellet_psi, pellet_delta_psi,                       &
                ellip,tria_u,tria_l,quad_u,quad_l,                  &
                xampl,xwidth,xsig,xtheta,xshift,xleft, xpoint,      &
                rho_file, T_file, ffprime_file, freeboundary_equil, &
                freeboundary, use_starwall, resistive_wall,         &
                refinement,                                         &
                produce_live_data

write(*,*) '***************************************'
write(*,*) '* JOREK2_diagno                       *'
write(*,*) '***************************************'

!read(5,in1)

call import_restart(node_list,element_list)

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
  write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
enddo

call initialise_basis                              ! define the basis functions at the Gaussian points

open(20,file='energies.txt')

write(20,'(A,25(A11,i3.3))') '      i      time',('          M',n_period*((in-1)/2),in=1,n_tor,2), &
                                                 ('          K',n_period*((in-1)/2),in=1,n_tor,2)

do i=2,index_start

 Growth_mag  = 0.5d0*log(abs(energies(n_tor,1,i)/energies(n_tor,1,i-1))) &
             / (xtime(i)-xtime(i-1))
 Growth_kin  = 0.5d0*log(abs(energies(n_tor,2,i)/energies(n_tor,2,i-1))) &
             / (xtime(i)-xtime(i-1))

 write(*,'(i7,f12.3,200e14.6)') i,xtime(i),energies(1:n_tor,:,i),growth_mag,growth_kin

 write(20,'(i7,f12.3,200e14.6)') i,xtime(i),energies(1,1,i),(energies(in,1,i)+energies(in+1,1,i),in=2,n_tor,2), &
                                            energies(1,2,i),(energies(in,2,i)+energies(in+1,2,i),in=2,n_tor,2)

enddo
close(20)

!call Integrals_3D(my_id,node_list,element_list,density,density_in,density_out,pressure,pressure_in,pressure_out)

!------------------lowshape3bis outside
!Rplot(1) = 3.0
!Rplot(2) = 3.676
!Zplot(1) = -2.066
!Zplot(2) = -1.9265

!------------------ lowshape3bis inside
!Rplot(1) = 2.3213
!Rplot(2) = 2.8511
!Zplot(1) = -1.9178
!Zplot(2) = -2.0339

!------------------ lowshape3bis midplane
!Rplot(1) = 1.88
!Rplot(2) = 4.88
!Zplot(1) = 0.09
!Zplot(2) = 0.09

!-----------------lowshape7,8 (midplane)
!Rplot(1) = 1.9
!Rplot(2) = 4.2
!Zplot(1) = 0.07
!Zplot(2) = 0.07

call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

Rplot(1) = 2.
Rplot(2) = 4.2
Zplot(1) = Z_axis
Zplot(2) = Z_axis 

call plot_profiles(node_list,element_list,Rplot,Zplot)


!call export_helena(node_list,element_list)

!----------------------------------------- plot profiles
!call begplt('profiles.ps')

!call plot_velocity_profile(node_list,element_list, 3.d0, 0.d0, 3.d0, 2.d0)

!call finplt

end


