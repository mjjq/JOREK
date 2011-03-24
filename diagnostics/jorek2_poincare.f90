module elements_nodes_neighbours
  
  use data_structure
  
  type (type_node_list)    :: node_list
  type (type_element_list) :: element_list
  integer,allocatable      :: element_neighbours(:,:)
  
end module

program jorek2_poincare
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
use basis_at_gaussian
use elements_nodes_neighbours

implicit none

real*8,allocatable  :: rp(:), zp(:), tp(:), pp(:)
integer :: i, j, iside_i, iside_j, ip, i_lines, n_lines, i_tor, i_harm, i_var_psi, iplot_type
integer :: i_elm, ifail, i_phi, n_phi, i_turn, n_turn, i_elm_out, i_elm_prev, i_elm_tmp,i_steps
real*8  :: R_start, Z_start, P_start, R_line, Z_line, s_line, t_line, p_line, R_mid, Z_mid, s_mid, t_mid, p_mid, s_out, t_out
real*8  :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, P, P_s, P_t, P_st, P_ss, P_tt
real*8  :: tol, delta_phi, Zjac, psi_s, psi_t, R_in, Z_in, R_out, Z_out, Rmin, Rmax, Zmin, Zmax, PI, delta_s, delta_t, R_keep, Z_keep
real*8  :: small_delta, small_delta_s, small_delta_t, delta_phi_local, delta_phi_step
real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis, atmp, cur_pert, psi_xpoint,R_xpoint,Z_xpoint,s_xpoint,t_xpoint, psi_bnd, psi_out
integer :: i_elm_axis, i_elm_xpoint

logical, external :: neighbours

namelist /in1/  tstep, nstep, eta, visco, visco_par, cur_pert,      &
                restart,  regrid,                                   &
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
                eta_num, visco_num, visco_par_num, D_perp_num,      &
                ellip,tria_u,tria_l,quad_u,quad_l,                  &
                xampl,xwidth,xsig,xtheta,xshift,xleft, xpoint

write(*,*) '***************************************'
write(*,*) '* JOREK2_poincare                     *'
write(*,*) '***************************************'

iplot_type = 2 ! 1: Poincare plot in (R,Z) coordinates, 2: in (R,theta) coordinates

read(5,in1)

call import_restart(node_list,element_list)

allocate(element_neighbours(4,element_list%n_elements))

element_neighbours = 0

do i=1,element_list%n_elements
  
  do j=i+1,element_list%n_elements
    
    if (neighbours(element_list%element(i),element_list%element(j),iside_i,iside_j)) then
      element_neighbours(iside_i,i) = j
      element_neighbours(iside_j,j) = i
    endif
    
  enddo
enddo

!-------- possibilities
! - total length of traced field line (connection length)
! - keep all RZ positions of field lines
! - keep only crossing with fixed plane phi=constant
! - (local) q-profile

! step at constant delta_phi

PI = 2.d0 *asin(1.d0)

n_lines = 50
n_turn  = 500
n_phi   = 500
 
delta_phi = 2.d0 * PI / float(n_period*n_phi)
tol       = 1.d-6

i_var_psi = 1

allocate(Rp(n_turn),Zp(n_turn), Tp(n_turn),Pp(n_turn))

Rmin = 1.d20; Rmax = -1.d20; Zmin = 1.d20; Zmax=-1.d20
do i=1,node_list%n_nodes
  Rmin = min(Rmin,node_list%node(i)%x(1,1))
  Rmax = max(Rmax,node_list%node(i)%x(1,1))
  Zmin = min(Zmin,node_list%node(i)%x(1,2))
  Zmax = max(Zmax,node_list%node(i)%x(1,2))
enddo

mode(1) = 0
do i=1,(n_tor-1)/2
  mode(2*i)   = i * n_period
  mode(2*i+1) = i * n_period
enddo
write(*,*) ' modes   : ',mode
write(*,*) ' nperiod : ',n_period

call find_axis(node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

psi_bnd = 0.d0
if ( xpoint ) then
  call find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,ifail)
  if ( ifail == 0 ) psi_bnd = psi_xpoint
end if
  
call begplt('poincare.ps')

open(21,file='poincare.txt')

if (iplot_type .eq. 1) then
  call nframe(21,11,1,Rmin,Rmax,Zmin,Zmax,'Poincare',8,'R [m]',4,'Z [m]',4)
else
  call nframe(1,11,1,0.d0,1.2d0,-PI,PI,'Poincare',8,'r [m]',4,'theta',5)
endif

do i_lines=1,n_lines

  ip = 0

  R_start = R_geo + 0.6 + 0.4*amin*float(i_lines-1)/float(n_lines-1)
!  R_start = R_geo + 0.18 + 0.8*amin*float(i_lines-1)/float(n_lines-1)
  Z_start = Z_geo
  P_start = 0.d0 !PI/2.d0

  write(*,'(2i6,2f8.3)') i_lines,n_lines,R_start,Z_start

  call find_RZ(node_list,element_list,R_start,Z_start,R_out,Z_out,i_elm,s_out,t_out,ifail)
  
  if (ifail .ne. 0) exit

  R_line = R_start
  Z_line = Z_start
  p_line = P_start
  s_line = s_out
  t_line = t_out
  
  do i_turn = 1, n_turn

    do i_phi=1,n_phi
    
      delta_phi_local = 0.d0
      
      i_steps = 0
    
      do while ((delta_phi_local .lt. delta_phi) .and. (i_steps .lt.10) )
      
        i_steps = i_steps + 1
      
        delta_phi_step = delta_phi - delta_phi_local
	
!	write(*,'(5i6,3e16.8)') i_lines,i_turn,i_phi,i_steps,i_elm,delta_phi_step, delta_phi, delta_phi_local
      
        call step(i_elm,s_line,t_line,p_line,delta_phi_step,delta_s,delta_t)

        s_mid = s_line + 0.5d0 * delta_s
        t_mid = t_line + 0.5d0 * delta_t
        p_mid = p_line + 0.5d0 * delta_phi_step
            
        call step(i_elm,s_mid,t_mid,p_mid,delta_phi_step,delta_s,delta_t)
        
        small_delta_s = 1.d0
      
        if  (s_line + delta_s .gt. 1.d0) then
      
          small_delta_s = (1.d0 - s_line)/delta_s
 
        elseif  (s_line + delta_s .lt. 0.d0) then

          small_delta_s = abs(s_line/delta_s)
	
        endif
      
        small_delta_t = 1.d0

        if  (t_line + delta_t .gt. 1.d0)  then
      
          small_delta_t = (1.d0 - t_line)/delta_t
	
        elseif  (t_line + delta_t .lt. 0.d0)  then

          small_delta_t = abs(t_line/delta_t)
      
        endif      

        small_delta = min(small_delta_s, small_delta_t)
	
!	write(*,'(A,5e16.8)') ' small delta : ',small_delta,delta_s,delta_t,s_line,t_line

        if (small_delta .lt. 1.d0)  then 

          s_mid = s_line + 0.5d0 * small_delta * delta_s
          t_mid = t_line + 0.5d0 * small_delta * delta_t
          p_mid = p_line + 0.5d0 * small_delta * delta_phi_step
            
          call step(i_elm,s_mid,t_mid,p_mid,delta_phi_step,delta_s,delta_t)
      
          if (small_delta_s .lt. small_delta_t) then

            if (s_line + delta_s .gt. 1.d0) then
	    
	      s_line = 1.d0
              t_line = t_line + small_delta * delta_t
              p_line = p_line + small_delta * delta_phi_step
	      
              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

	      i_elm_prev = i_elm      
              i_elm      = element_neighbours(2,i_elm_prev)
	      i_elm_tmp  = element_neighbours(4,i_elm)
	      
	      if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (1)'
	
              s_line = 0.d0

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)
	      
	      if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
	        write(*,'(A,2i6,4f8.4)') ' error in element change (1) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
	
            elseif (s_line + delta_s .lt. 0.d0) then
	
              s_line = 0.d0
	      t_line = t_line + small_delta * delta_t
	      p_line = p_line + small_delta * delta_phi_step

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)
      
	      i_elm_prev = i_elm      
              i_elm      = element_neighbours(4,i_elm_prev)
	      i_elm_tmp  = element_neighbours(2,i_elm)
	      
	      if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (2)'

              s_line = 1.d0

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)
	      
	      if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
	        write(*,'(A,2i6,4f8.4)') ' error in element change (2) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
	  
	    endif
	
	  else
	
            if (t_line + delta_t .gt. 1.d0) then
      
              s_line = s_line + small_delta * delta_s
              t_line = 1.d0
              p_line = p_line + small_delta * delta_phi_step

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

	      i_elm_prev = i_elm      
              i_elm      = element_neighbours(3,i_elm_prev)
	      i_elm_tmp  = element_neighbours(1,i_elm)
	      
	      if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (3)'
	
              t_line = 0.d0

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)
	      
	      if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
	        write(*,'(A,2i6,4f8.4)') ' error in element change (3) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
	
            elseif (t_line + delta_t .lt. 0.d0) then

              s_line = s_line + small_delta * delta_s	
	      t_line = 0.d0
	      p_line = p_line + small_delta * delta_phi_step
 
              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

	      i_elm_prev = i_elm      
              i_elm      = element_neighbours(1,i_elm_prev)
	      i_elm_tmp  = element_neighbours(3,i_elm)
	      
	      if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (4)',i_elm_prev,i_elm
 
              t_line = 1.d0

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)
	      
	      if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
	        write(*,'(A,2i6,4f8.4)') ' error in element change (4) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
	
            endif
        
	  endif
         
        else      

          s_line = s_line + delta_s
          t_line = t_line + delta_t
          p_line = p_line + delta_phi_step
	
	  small_delta = 1.d0
      
        endif
      
        delta_phi_local = delta_phi_local + small_delta * delta_phi_step
     
        if (i_elm .eq. 0) exit
      
!        write(*,'(A,5e16.8)') ' s,t : ',s_line,t_line
      
      enddo

      if (i_elm .eq. 0) exit
      
!      if (i_steps .gt. 8) write(*,'(A,5i6)') ' WARNING : isteps ',i_lines,i_turn,i_phi,i_steps,i_elm
            
    enddo ! end of a 2Pi turn

    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)

    R_line = R
    Z_line = Z
            
    ip = ip+1

    call var_value(i_elm,1,s_line,t_line,p_line,psi_out)

    Rp(ip) = R_line
    Zp(ip) = Z_line
    Tp(ip)  = atan2( Z_line - Z_axis, R_line - R_axis)
    Pp(ip)  = (psi_out - psi_axis)/(psi_bnd - psi_axis)

    if (i_elm .eq. 0) exit
     
  enddo
  
  write(*,*) ' points : ',ip

  do i=1,ip
    write(21,'(4e18.8)') Rp(i),Zp(i),Pp(i),Tp(i)
  enddo
  
  call lincol(mod(i_lines,8))
   
  if (iplot_type .eq. 1) then
    call pplot(1,1,Rp,Zp,ip,1)
  else
    call pplot(1,1,Pp,Tp,ip,1)
  endif
  
enddo

close(21)
call finplt

end

subroutine step(i_elm,s_in,t_in,p_in,delta_p,delta_s,delta_t)
use parameters
use elements_nodes_neighbours
use phys_module

implicit none

integer :: i_var_psi, i_elm, i_tor, i_harm

real*8 :: s_in, t_in, p_in, delta_p, delta_s, delta_t
real*8 :: R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
real*8 :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt, Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
real*8 :: P0,P0_s,P0_t,P0_st,P0_ss,P0_tt, psi_s, psi_t, Zjac

i_var_psi = 1

call interp_RZ(node_list,element_list,i_elm,s_in,t_in,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)

Zjac = (R_s * Z_t - R_t * Z_s)

call interp(node_list,element_list,i_elm,i_var_psi,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)

psi_s = P0_s 
psi_t = P0_t 

do i_tor = 1, (n_tor-1)/2

  i_harm = 2*i_tor

  call interp(node_list,element_list,i_elm,i_var_psi,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)

  psi_s = psi_s + Pcos_s * cos(mode(i_harm)*p_in)
  psi_t = psi_t + Pcos_t * cos(mode(i_harm)*p_in)

  call interp(node_list,element_list,i_elm,i_var_psi,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)

  psi_s = psi_s + Psin_s * sin(mode(i_harm+1)*p_in)
  psi_t = psi_t + Psin_t * sin(mode(i_harm+1)*p_in)

enddo

delta_s =   psi_t * R / (Zjac * F0) * delta_p
delta_t = - psi_s * R / (Zjac * F0) * delta_p

return
end

subroutine var_value(i_elm,i_var,s_in,t_in,p_in,value_out)
use parameters
use elements_nodes_neighbours
use phys_module

implicit none

integer :: i_var, i_elm, i_tor, i_harm

real*8 :: s_in, t_in, p_in
real*8 :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt, Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
real*8 :: P0,P0_s,P0_t,P0_st,P0_ss,P0_tt
real*8 :: value_out


!call interp_RZ(node_list,element_list,i_elm,s_in,t_in,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
!Zjac = (R_s * Z_t - R_t * Z_s)

call interp(node_list,element_list,i_elm,i_var,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)

value_out = P0

do i_tor = 1, (n_tor-1)/2

  i_harm = 2*i_tor

  call interp(node_list,element_list,i_elm,i_var,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)

  value_out = value_out + Pcos * cos(mode(i_harm)*p_in)

  call interp(node_list,element_list,i_elm,i_var,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)

  value_out = value_out + Psin * sin(mode(i_harm+1)*p_in)

enddo

return
end
