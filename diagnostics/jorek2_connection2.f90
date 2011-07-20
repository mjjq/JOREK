module elements_nodes_neighbours
  
  use data_structure
  
  type (type_node_list)    :: node_list
  type (type_element_list) :: element_list
  integer,allocatable      :: element_neighbours(:,:)
  
end module



program jorek2_connection2

use data_structure
use phys_module
use basis_at_gaussian
use elements_nodes_neighbours

implicit none

real*8,allocatable  :: rp(:), zp(:), R_all(:), Z_all(:), C_all(:), R_strike(:), Z_strike(:), P_strike(:), C_strike(:)
real*8,allocatable  :: R_turn(:,:), Z_turn(:,:), C_turn(:,:), C_turn_tmp(:,:), T_turn(:,:), PSI_turn(:,:), ZN_turn(:,:)
integer :: i, j, iside_i, iside_j, ip, i_line, n_lines, i_tor, i_harm, i_var_psi, i_dir, k, m, nr, np
integer :: i_elm, ifail, i_phi, n_phi, i_turn, n_turns, i_elm_out, i_elm_prev, i_elm_tmp,i_steps, n_turn_max(2)
real*8  :: R_start, Z_start, P_start, R_line, Z_line, s_line, t_line, p_line, s_mid, t_mid, p_mid, s_out, t_out
real*8  :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, P, P_s, P_t, P_st, P_ss, P_tt
real*8  :: tol, delta_phi, Zjac, psi_s, psi_t, R_in, Z_in, R_out, Z_out, Rmin, Rmax, Zmin, Zmax, PI, delta_s, delta_t, R_keep, Z_keep
real*8  :: small_delta, small_delta_s, small_delta_t, delta_phi_local, delta_phi_step, total_phi
real*8  :: Rmid,Zmid,Rmid_s,Rmid_t,Zmid_s,Zmid_t, dl2, total_length, length_max, s_ini, t_ini, zl1, zl2, partial(2)
real*8  :: psi_xpoint,R_xpoint,Z_xpoint,s_xpoint,t_xpoint, value_out
integer :: i_elm_xpoint, my_id, ierr

logical, external :: neighbours

write(*,*) '***************************************'
write(*,*) '* JOREK2_connection2                  *'
write(*,*) '***************************************'
write(*,*) ' nperiod : ',n_period

my_id=0

call initialise_parameters(my_id)

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

n_turns = 500
n_phi   = 1000

np = 1
nr = 1

delta_phi = 2.d0 * PI / float(n_period*n_phi)
tol       = 1.d-6

i_var_psi = 1

n_lines = element_list%n_elements * np * np

allocate(R_strike(n_lines),Z_strike(n_lines),P_strike(n_lines),C_strike(n_lines))

allocate(R_all(n_lines),Z_all(n_lines),C_all(n_lines))
allocate(R_turn(n_turns+1,2),Z_turn(n_turns+1,2),C_turn(n_turns+1,2),C_turn_tmp(n_turns+1,2))
allocate(T_turn(n_turns+1,2),PSI_turn(n_turns+1,2),ZN_turn(n_turns+1,2))

R_all    = 0.d0; Z_all    = 0.d0; C_all = 0.d0
R_strike = 0.d0; Z_strike = 0.d0; P_strike = 0.d0; C_strike = 0.d0
R_turn   = 0.d0; Z_turn   = 0.d0; C_turn   = 0.d0; C_turn_tmp   = 0.d0;

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

write(*,*) ' modes : ',mode
write(*,*) ' nperiod : ',n_period

call initialise_basis                              ! define the basis functions at the Gaussian points

call find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint)

write(*,*) ' PSI_XPOINT : ',psi_xpoint,i_elm_xpoint
  
call begplt('connection2e.ps')
call nframe(21,11,1,Rmin,Rmax,Zmin,Zmax,'Poincare',8,'R [m]',4,'Z [m]',4)

i_line = 0

write(*,*) ' number of elements : ',element_list%n_elements

open(20,file='connection2e.txt')
open(21,file='connection2e_inside.txt')

do i = 6*element_list%n_elements/8, 8*element_list%n_elements/8

  write(*,*) ' element : ',i

  do k=1, nr

    s_ini = (real(k)-0.5)/real(nr)

    do m=1, np

      t_ini = (real(m)-0.5)/real(np)

      i_line = i_line + 1
      
      R_turn = 0.d0
      Z_turn = 0.d0
      C_turn_tmp = 0.d0
      
      do i_dir = -1,1,2
      
      s_line = s_ini
      t_line = t_ini

      delta_phi = 2.d0 * PI * float(i_dir) / float(n_period*n_phi)

      call interp_RZ(node_list,element_list,i,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)
  
    
      total_length = 0.d0
      total_phi    = 0.d0
    
      i_elm = i
      R_start = R_out
      Z_start = Z_out
      P_start = 0.d0
  
      R_all(i_line) = R_start
      Z_all(i_line) = Z_start
      
      R_turn(1,(i_dir+1)/2+1) = R_start
      Z_turn(1,(i_dir+1)/2+1) = Z_start
      C_turn(1,(i_dir+1)/2+1) = 0.d0

      call var_value(i_elm,6,s_line,t_line,P_start,T_turn(1,(i_dir+1)/2+1))
      call var_value(i_elm,1,s_line,t_line,P_start,PSI_turn(1,(i_dir+1)/2+1))
      call var_value(i_elm,5,s_line,t_line,P_start,ZN_turn(1,(i_dir+1)/2+1))
        
      R_line = R_start
      Z_line = Z_start
      p_line = P_start
  
      do i_turn = 1, n_turns

        n_turn_max((i_dir+1)/2+1) = i_turn

        do i_phi=1,n_phi
    
          delta_phi_local = 0.d0
      
          i_steps = 0
    
          do while ((abs(delta_phi_local) .lt. abs(delta_phi)) .and. (i_steps .lt.10) )
      
            i_steps = i_steps + 1
      
            delta_phi_step = delta_phi - delta_phi_local
	
!	write(*,'(5i6,3e16.8)') i_line,i_turn,i_phi,i_steps,i_elm,delta_phi_step, delta_phi, delta_phi_local
      
            call step(i_elm,s_line,t_line,p_line,delta_phi_step,delta_s,delta_t,R,Z,R_s,R_t,Z_s,Z_t)

            s_mid = s_line + 0.5d0 * delta_s
            t_mid = t_line + 0.5d0 * delta_t
            p_mid = p_line + 0.5d0 * delta_phi_step
            
            call step(i_elm,s_mid,t_mid,p_mid,delta_phi_step,delta_s,delta_t,Rmid,Zmid,Rmid_s,Rmid_t,Zmid_s,Zmid_t)
         
            small_delta_s = 1.d0
      
            if  (s_line + delta_s .gt. 1.d0) then
      
              small_delta_s = (1.d0 - s_line)/delta_s
	      
!	      write(*,*) ' 1 : ',small_delta_s,s_line,delta_s
 
            elseif  (s_line + delta_s .lt. 0.d0) then

              small_delta_s = abs(s_line/delta_s)
	      
!	      write(*,*) ' 2 : ',small_delta_s,s_line,delta_s
	
            endif
      
            small_delta_t = 1.d0

            if  (t_line + delta_t .gt. 1.d0)  then
      
              small_delta_t = (1.d0 - t_line)/delta_t

!	      write(*,*) ' 3 : ',small_delta_t,t_line,delta_t
	
            elseif  (t_line + delta_t .lt. 0.d0)  then

              small_delta_t = abs(t_line/delta_t)

!	      write(*,*) ' 4 : ',small_delta_t,t_line,delta_t
      
            endif      

            small_delta = min(small_delta_s, small_delta_t)
	
!           write(*,'(A,5e16.8)') ' small delta : ',small_delta,delta_s,delta_t,s_line,t_line

            if (small_delta .lt. 1.d0)  then 

              s_mid = s_line + 0.5d0 * small_delta * delta_s
              t_mid = t_line + 0.5d0 * small_delta * delta_t
              p_mid = p_line + 0.5d0 * small_delta * delta_phi_step
            
             call step(i_elm,s_mid,t_mid,p_mid,delta_phi_step,delta_s,delta_t,Rmid,Zmid,Rmid_s,Rmid_t,Zmid_s,Zmid_t)
      
             if (small_delta_s .lt. small_delta_t) then

               if (s_line + delta_s .gt. 1.d0) then
	    
	          s_line = 1.d0
                  t_line = t_line + small_delta * delta_t
                  p_line = p_line + small_delta * delta_phi_step
	      
	          dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
	          dl2 = dl2 * small_delta
	      
                  call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

	          i_elm_prev = i_elm      
                  i_elm      = element_neighbours(2,i_elm_prev)
	           i_elm_tmp  = element_neighbours(4,i_elm)
	      
	          if (i_elm .ne. 0) then
	      
	            if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (1)'
	
                    s_line = 0.d0
 
                    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)
	      
	            if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
	              write(*,'(A,2i6,4f8.4)') ' error in element change (1) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
	      
	          endif
	
                elseif (s_line + delta_s .lt. 0.d0) then
	
                  s_line = 0.d0
	          t_line = t_line + small_delta * delta_t
  	          p_line = p_line + small_delta * delta_phi_step

	          dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
	          dl2 = dl2 * small_delta
 
                  call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)
      
	          i_elm_prev = i_elm      
                  i_elm      = element_neighbours(4,i_elm_prev)
	          i_elm_tmp  = element_neighbours(2,i_elm)

	          if (i_elm .ne. 0) then
	      
                    if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (2)'

                    s_line = 1.d0
 
                    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)
	        
	            if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8))  &
	              write(*,'(A,2i6,4f8.4)') ' error in element change (2) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
	       
	          endif
	    
	        endif
	
	      else
	
                if (t_line + delta_t .gt. 1.d0) then
      
                  s_line = s_line + small_delta * delta_s
                  t_line = 1.d0
                  p_line = p_line + small_delta * delta_phi_step

	          dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
	          dl2 = dl2 * small_delta

                  call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

	          i_elm_prev = i_elm      
                  i_elm      = element_neighbours(3,i_elm_prev)
	          i_elm_tmp  = element_neighbours(1,i_elm)
	      
	          if (i_elm .ne. 0) then
	      
	            if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (3)'
	
                    t_line = 0.d0

                    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)
	      
	            if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8))  &
	              write(*,'(A,2i6,4f8.4)') ' error in element change (3) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
	       
	          endif
	      
                elseif (t_line + delta_t .lt. 0.d0) then

                  s_line = s_line + small_delta * delta_s	
	          t_line = 0.d0
	          p_line = p_line + small_delta * delta_phi_step

	          dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
	          dl2 = dl2 * small_delta
 
                  call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

	          i_elm_prev = i_elm      
                  i_elm      = element_neighbours(1,i_elm_prev)
	          i_elm_tmp  = element_neighbours(3,i_elm)
	      
	          if (i_elm .ne. 0) then
	      
	            if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (4)'
  
                    t_line = 1.d0

                    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)
	      
	            if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8))  &
	              write(*,'(A,2i6,4f8.4)') ' error in element change (4) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
	      
	          endif
	      
                endif
        
	      endif
         
            else      

              s_line = s_line + delta_s
              t_line = t_line + delta_t
              p_line = p_line + delta_phi_step

              dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
	
	      small_delta = 1.d0

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)
      
            endif
      
            delta_phi_local = delta_phi_local + small_delta * delta_phi_step
	
	    total_length = total_length + sqrt(abs(dl2))
	    total_phi    = total_phi    + small_delta * delta_phi_step
     
            if (i_elm .eq. 0) exit
      
!            write(*,'(8e16.8)') total_phi, R_in, Z_in, total_length
      
          enddo

          if (i_elm .eq. 0) exit
      
!          if (i_steps .gt. 8) write(*,'(A,5i6)') ' WARNING : isteps ',i_line,i_turn,i_phi,i_steps,i_elm
            
        enddo ! end of a 2Pi turn (or before if end of open field line)

        if (i_elm .eq. 0) exit

        call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)
	
	R_turn(i_turn+1,(i_dir+1)/2+1) = R_in
	Z_turn(i_turn+1,(i_dir+1)/2+1) = Z_in
	C_turn_tmp(i_turn+1,(i_dir+1)/2+1) = total_length
	
	call var_value(i_elm,6,s_line,t_line,p_line,T_turn(i_turn+1,(i_dir+1)/2+1))
	call var_value(i_elm,1,s_line,t_line,p_line,PSI_turn(i_turn+1,(i_dir+1)/2+1))
	call var_value(i_elm,5,s_line,t_line,p_line,ZN_turn(i_turn+1,(i_dir+1)/2+1))
	
!	write(*,'(2i5,8e16.8)') i_turn,i_dir,R_turn(i_turn+1,(i_dir+1)/2+1), Z_turn(i_turn+1,(i_dir+1)/2+1),PSI_turn(i_turn+1,(i_dir+1)/2+1)
!	write(*,'(A,i5,4e16.8)') ' c_turn_tmp : ',i_turn,R_in,Z_in,C_turn_tmp(i_turn+1,(i_dir+1)/2+1)
               
      enddo  ! end of loop over toroidal turns

      if (i_dir .eq. -1) then  
        C_all(i_line) = total_length
	partial(1)      = total_length
      else
        C_all(i_line) = min(C_all(i_line),total_length)
	partial(2)      = total_length
      endif
       
      enddo  ! end of two directions
      
!------------------------------- correct the connection lengths
!      write(*,*) ' total_length : ',total_length
!      write(*,*) ' n_turn_max   : ',n_turn_max
      
      do i_turn = 1, n_turn_max(1)
        C_turn(i_turn,1) = partial(1) - c_turn_tmp(i_turn,1)
!	write(*,'(i5,8e16.8)') i_turn,C_turn(i_turn,1),C_turn_tmp(i_turn,1),partial(1)
      enddo    
      do i_turn = 1, n_turn_max(2)
        C_turn(i_turn,2) = partial(2) - c_turn_tmp(i_turn,2)
!	write(*,'(i5,8e16.8)') i_turn,C_turn(i_turn,2),C_turn_tmp(i_turn,2),partial(2)
      enddo    
      
!      if (i_elm .eq. 0) then
!        R_strike(i_line) = R_line
!        Z_strike(i_line) = Z_line
!        P_strike(i_line) = p_line
!        C_strike(i_line) = total_length
!      endif

      do i_turn=1,n_turn_max(1)+1
        
	if (R_turn(i_turn,1) .gt. 0.d0) then
  
          zl1 = C_turn(i_turn,1)
	  zl2 = C_turn(1,1) - C_turn(i_turn,1) + C_turn(1,2) 
	  	  
	  if ((PSI_turn(1,1) .lt. psi_xpoint).and. (Z_turn(1,1) .gt. Z_xpoint)) then
	  	  
            if (n_turn_max(1) .lt. n_turns) then
              write(21,'(12e16.8)') R_turn(i_turn,1),Z_turn(i_turn,1),min(zl1,zl2),T_turn(i_turn,1),PSI_turn(i_turn,1),ZN_turn(i_turn,1)
	    else
              write(21,'(12e16.8)') R_turn(i_turn,1),Z_turn(i_turn,1),maxval(partial),T_turn(i_turn,1),PSI_turn(i_turn,1),ZN_turn(i_turn,1)
	    endif 
	     
          else
	  
	    if (n_turn_max(1) .lt. n_turns) then
              write(20,'(12e16.8)') R_turn(i_turn,1),Z_turn(i_turn,1),min(zl1,zl2),T_turn(i_turn,1),PSI_turn(i_turn,1),ZN_turn(i_turn,1)
	    else
              write(20,'(12e16.8)') R_turn(i_turn,1),Z_turn(i_turn,1),maxval(partial),T_turn(i_turn,1),PSI_turn(i_turn,1),ZN_turn(i_turn,1)
	    endif  
	  endif
	  
	endif
      enddo
      
      do i_turn=1,n_turn_max(2)+1
	if (R_turn(i_turn,2) .gt. 0.d0) then
  
          zl1 = C_turn(i_turn,2)
	  zl2 = C_turn(1,2) - C_turn(i_turn,2) + C_turn(1,1) 

	  if ((PSI_turn(1,2) .lt. psi_xpoint) .and. (Z_turn(1,2) .gt. Z_xpoint)) then
	  
	    if (n_turn_max(2) .lt. n_turns) then
              write(21,'(12e16.8)') R_turn(i_turn,2),Z_turn(i_turn,2),min(zl1,zl2),T_turn(i_turn,2),PSI_turn(i_turn,2),ZN_turn(i_turn,2)
	    else
              write(21,'(12e16.8)') R_turn(i_turn,2),Z_turn(i_turn,2),maxval(partial),T_turn(i_turn,2),PSI_turn(i_turn,2),ZN_turn(i_turn,2)
	   endif  
          
	  else
	  
	    if (n_turn_max(2) .lt. n_turns) then
              write(20,'(12e16.8)') R_turn(i_turn,2),Z_turn(i_turn,2),min(zl1,zl2),T_turn(i_turn,2),PSI_turn(i_turn,2),ZN_turn(i_turn,2)
	    else
              write(20,'(12e16.8)') R_turn(i_turn,2),Z_turn(i_turn,2),maxval(partial),T_turn(i_turn,2),PSI_turn(i_turn,2),ZN_turn(i_turn,2)
	    endif  
	  
	  endif
	  
	endif
      enddo

    enddo
  enddo

  
enddo ! end of loop over elements

close(20)
close(21)

call finplt

end program jorek2_connection2



subroutine step(i_elm,s_in,t_in,p_in,delta_p,delta_s,delta_t,R,Z,R_s,R_t,Z_s,Z_t)
use parameters
use elements_nodes_neighbours
use phys_module

implicit none

integer :: i_var_psi, i_elm, i_tor, i_harm

real*8 :: s_in, t_in, p_in, delta_p, delta_s, delta_t
real*8 :: R_out, Z_out, Rs_out, Rt_out, Zs_out, Zt_out
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
end subroutine step



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
end subroutine var_value
