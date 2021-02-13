subroutine equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint2,xcase2, nice_q)
!-----------------------------------------------------------------------
! Solve the Grad-Shafranov equation to determine the plasma equilibrium
!   both freeboundary and fixed boundary solutions
!-----------------------------------------------------------------------
use tr_module 
use mod_parameters
use data_structure
use phys_module
use mod_poiss
use mod_iterate2area
use mod_plasma_response
use equil_info
use vacuum
use mpi_mod
use mod_interp, only: interp
use mod_F_profile
implicit none

          
! --- Routine parameters
integer,                      intent(in)    :: my_id
type (type_node_list),        intent(inout) :: node_list
type (type_element_list),     intent(inout) :: element_list
type (type_bnd_node_list),    intent(inout) :: bnd_node_list
type (type_bnd_element_list), intent(inout) :: bnd_elm_list
logical,                      intent(in)    :: xpoint2
integer,                      intent(in)    :: xcase2
logical,                      intent(in)    :: nice_q

! --- Local variables.
type (type_surface_list) :: surface_list, sep_list
integer    :: ierr, n_iter, iter, i, in, mm, i_elm_axis, i_elm_xpoint(2), i_elm_lim, ifail, i_elm
real*8     :: amplitude, psi, psi_bnd
real*8     :: zn
real*8     ::    dn_dpsi, dn_dz                                                ! 1st order derivatives
real*8     ::    dn_dpsi2, dn_dz2, dn_dpsi_dz                                  ! 2nd order derivatives
real*8     ::    dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,  dn_dz3                   ! 2rd order derivatives
real*8     ::    dn_dpsi4, dn_dpsi_dz3, dn_dpsi2_dz2, dn_dpsi3_dz, dn_dz4      ! 4th order derivatives
real*8     ::    dn_dpsi5, dn_dpsi_dz4, dn_dpsi2_dz3, dn_dpsi3_dz2, dn_dpsi4_dz! 5th order derivatives (z5 not needed)
real*8     :: zT
real*8     ::    dT_dpsi,  dT_dz                                               ! 1st order derivatives
real*8     ::    dT_dpsi2, dT_dz2, dT_dpsi_dz                                  ! 2nd order derivatives
real*8     ::    dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz,  dT_dz3                   ! 2rd order derivatives
real*8     ::    dT_dpsi4, dT_dpsi_dz3, dT_dpsi2_dz2, dT_dpsi3_dz, dT_dz4      ! 4th order derivatives
real*8     ::    dT_dpsi5, dT_dpsi_dz4, dT_dpsi2_dz3, dT_dpsi3_dz2, dT_dpsi4_dz! 5th order derivatives (z5 not needed)
real*8     :: zTi
real*8     ::    dTi_dpsi,  dTi_dz                                                  ! 1st order derivatives
real*8     ::    dTi_dpsi2, dTi_dz2, dTi_dpsi_dz                                    ! 2nd order derivatives
real*8     ::    dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz,  dTi_dz3                    ! 2rd order derivatives
real*8     ::    dTi_dpsi4, dTi_dpsi_dz3, dTi_dpsi2_dz2, dTi_dpsi3_dz, dTi_dz4      ! 4th order derivatives
real*8     ::    dTi_dpsi5, dTi_dpsi_dz4, dTi_dpsi2_dz3, dTi_dpsi3_dz2, dTi_dpsi4_dz! 5th order derivatives (z5 not needed)
real*8     :: zTe
real*8     ::    dTe_dpsi,  dTe_dz                                                  ! 1st order derivatives
real*8     ::    dTe_dpsi2, dTe_dz2, dTe_dpsi_dz                                    ! 2nd order derivatives
real*8     ::    dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz,  dTe_dz3                    ! 2rd order derivatives
real*8     ::    dTe_dpsi4, dTe_dpsi_dz3, dTe_dpsi2_dz2, dTe_dpsi3_dz, dTe_dz4      ! 4th order derivatives
real*8     ::    dTe_dpsi5, dTe_dpsi_dz4, dTe_dpsi2_dz3, dTe_dpsi3_dz2, dTe_dpsi4_dz! 5th order derivatives (z5 not needed)
real*8     :: Ti_prof, Te_prof
real*8     :: zFFprime
real*8     ::    dFF_dpsi, dFF_dz                                                ! 1st order derivatives
real*8     ::    dFF_dpsi2, dFF_dz2, dFF_dpsi_dz                                 ! 2nd order derivatives
real*8     ::    dFF_dpsi3, dFF_dpsi_dz2, dFF_dpsi2_dz,  dFF_dz3                 ! 2rd order derivatives
real*8     ::    dFF_dpsi4, dFF_dpsi_dz3, dFF_dpsi2_dz2, dFF_dpsi3_dz, dFF_dz4   ! 4th order derivatives
real*8     :: F_prof   
real*8     ::   dF_dpsi, dF_dz                                             ! 1st order derivatives
real*8     ::   dF_dpsi2, dF_dz2, dF_dpsi_dz                               ! 2nd order derivatives
real*8     ::   dF_dpsi3, dF_dpsi_dz2, dF_dpsi2_dz,  dF_dz3                ! 2rd order derivatives
real*8     ::   dF_dpsi4, dF_dpsi_dz3, dF_dpsi2_dz2, dF_dpsi3_dz, dF_dz4   ! 4th order derivatives
real*8     :: F_values(0:4,0:4,0:4)        !< variable and its derivatives
real*8     :: xx, x_s, x_t, x_st, x_ss, x_tt, yy, y_s, y_t, y_st, y_ss, y_tt
real*8     :: R_axis, Z_axis, s_axis, t_axis, psi_axis,R, Z, BigR, T0, BigR_s, T0_s
real*8     :: R_lim, Z_lim, s_lim, t_lim, psi_lim, R_out, Z_out, s_out, t_out
real*8     :: R_xpoint(2),Z_xpoint(2),s_xpoint(2),t_xpoint(2), psi_xpoint(2)
real*8     :: zjz, psi_n
real*8     :: GP
real*8     ::    dGP_dpsi,  dGP_dz                                                  ! 1st order derivatives
real*8     ::    dGP_dpsi2, dGP_dz2,      dGP_dpsi_dz                               ! 2nd order derivatives
real*8     ::    dGP_dpsi3, dGP_dpsi_dz2, dGP_dpsi2_dz,  dGP_dz3                    ! 2rd order derivatives
real*8     ::    dGP_dpsi4, dGP_dpsi_dz3, dGP_dpsi2_dz2, dGP_dpsi3_dz,  dGP_dz4     ! 4th order derivatives
real*8     ::    dGP_dpsi5, dGP_dpsi_dz4, dGP_dpsi2_dz3, dGP_dpsi3_dz2, dGP_dpsi4_dz! 5th order derivatives (z5 not needed)
real*8     :: ps0_s, ps0_t, p_s, p_t, p_ss, p_st, p_tt 
real*8     :: zj0_s, zj0_t, equil_error, equil_value, ps0_x, ps0_y, Z_s, Z_t, xjac, direction, Btot
real*8     :: current_tot, current_int, diff, R_xpoint2(2), Z_xpoint2(2)
real*8     :: sigmas(16), dZ_axis, Z_axis_int, Z_axis_old, area_ref
integer    :: n_grids(10)
logical    :: freeboundary_equil2
real*8     :: T_prof, T_0_old, FF_0_old, T_1_old, FF_1_old
real*8, allocatable     :: T_profile(:)
real*8     :: density_prof
real*8, allocatable     :: density_profile(:)

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*           equilibrium               *'
  write(*,*) '***************************************'
  write(*,*) '   freeboundary_equil : ',freeboundary_equil
  write(*,*) '   X-point      : ',xpoint2
  write(*,*) '   Xcase        : ',xcase2
endif

freeboundary_equil2 = freeboundary_equil
freeboundary_equil  = .false.

!------------------------------------ fixed boundary equilibrium
n_iter       = 200
psi_bnd      = 0.d0
Z_xpoint(1)  = -99.d0
Z_xpoint(2)  = +99.d0
vertical_FB  = 0.d0
i_elm_xpoint = 0 

if (my_id == 0) then

  do iter = 1, n_iter
  
    call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
  
    if ((ifail .ne. 0) .and. (iter .le. 5)) then
      call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
      call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
      write(*,'(A,3f10.5)')  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
    endif
    
    if (xpoint2) then
      call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint2,Z_xpoint2,i_elm_xpoint,s_xpoint,t_xpoint,xcase2,ifail)
      if (ifail == 0) then ! (otherwise, keep the values of the previous iteration as a reasonable guess)
        psi_bnd  = psi_xpoint(1)
        if( (xcase2 .eq. 2) .or. ((xcase2 .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
          psi_bnd = psi_xpoint(2)
        endif
        R_xpoint(1) = R_xpoint2(1)
        Z_xpoint(1) = Z_xpoint2(1)
        R_xpoint(2) = R_xpoint2(2)
        Z_xpoint(2) = Z_xpoint2(2)
        if(xcase2 .eq. 1) Z_xpoint(2) = +99.d0
        if(xcase2 .eq. 2) Z_xpoint(1) = -99.d0
      else
        if (freeboundary_equil) then
          Z_xpoint(1) = -99.d0
          Z_xpoint(2) = +99.d0
        endif
      endif
    else
      psi_bnd = 0.d0
    endif
    if (.not. xpoint) then
      call find_limiter(my_id,node_list,element_list,bnd_elm_list,psi_lim,R_lim,Z_lim)
      if ( (Z_lim .gt. Z_xpoint(1)) .and. (Z_lim .lt. Z_xpoint(2)) ) then
        if (n_limiter /= 0) then   ! else n_limiter = 0 and psi_bnd is set to 0
          psi_bnd = psi_lim
          write(*,'(A,3f8.3)') ' LIMITER PLASMA ',psi_lim,R_lim,Z_lim
        endif
      endif
    endif
  
    if(xcase2 .eq. 1) write(*,'(A,3es14.6,i3)') ' PSI_AXIS, PSI_BND : ',psi_axis,psi_bnd,Z_xpoint(1),ifail
    if(xcase2 .eq. 2) write(*,'(A,3es14.6,i3)') ' PSI_AXIS, PSI_BND : ',psi_axis,psi_bnd,Z_xpoint(2),ifail
  
    call poisson(my_id,-1,node_list,element_list,bnd_node_list,bnd_elm_list,3,1,1, &
                 psi_axis,psi_bnd,xpoint2,xcase2,Z_xpoint,freeboundary_equil,refinement,iter)   !----------- for GS use -1
  
    diff = 0.d0
    do i=1, node_list%n_nodes
      diff = diff + abs(node_list%node(i)%deltas(1,1,1))
    enddo  
    diff = diff / float(node_list%n_nodes)
    
    write(*,'(A,I4,A,ES10.3)') ' Iteration ', iter, ': diff=', diff
    
    if ( (iter > 1) .and. (diff < equil_accuracy) ) then
      write(*,'(A,I4,A)') ' Fixed boundary equilibrium converged: after', iter, ' iterations'
      exit
    else if ( iter == n_iter) then
      write(*,'(A,ES10.3)') ' WARNING: Fixed boundary equilibrium not fully converged: diff=', diff
      exit
    end if
  
  enddo

end if ! my_id == 0

!--------------------------------------- freeboundary equilibrium
freeboundary_equil = freeboundary_equil2

current_int = 0.d0; Z_axis_int = 0.d0
 
T_0_old = T_0;  FF_0_old = FF_0;  T_1_old = T_1;  FF_1_old = FF_1

if (freeboundary_equil) then

  if (my_id == 0) then

    write(*,*)
    write(*,*) '------------------------------------------------------'
    write(*,*) '--- Iterative solution of freeboundary equilibrium ---'
    write(*,*) '------------------------------------------------------'
    write(*,*)

    ! Target current and axis
    if (current_ref .gt. 1.d20) then    !choose fix bnd equilibrium final current in case of non specification of target current
      call integral_current(node_list,element_list,psi_axis, psi_bnd, xpoint2, xcase2, Z_xpoint, current_ref)
    endif
   
    if (Z_axis_ref .gt. 1.d20) then     !choose fix bnd equilibrium final Zaxis in case of non specification of target Zaxis
      Z_axis_ref = Z_axis
    endif
    
    ! Target poloidal cross section area for limiter plasmas
    if (freeb_equil_iterate_area .and. (.not. xpoint2)) then
      n_limiter = 0    ! Use the full domain to search psibnd enclosing given area
      call area_inside_flux_contour(node_list,element_list, xpoint2, xcase2, psi_bnd, area_ref, R_lim, Z_lim)
      write(*,*) ' The reference area from fixed boundaray is = ', area_ref
    endif
  
  end if ! my_id == 0

  do iter=1, n_iter_freeb

    if (my_id == 0) then
      
      write(*,*)
      write(*,'(1x,a,i5,a)') '>>> ITERATION', iter, ' <<<'
                     
      ! Update axis annd boundary values
      call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
      
      write(10,'(i6,9e20.12)') iter, current_tot, R_axis, Z_axis, psi_bnd-psi_axis
  
      if ((ifail .ne. 0) .and. (iter .le. 5)) then
        call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
        call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
        write(*,*)  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
      endif
      
      psi_bnd = 0.d0
   
      if (xpoint2) then
        call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase2,ifail)
        if (ifail .ne. 1) then      
          psi_bnd  = psi_xpoint(1)
          if( (xcase2 .eq. 2) .or. ((xcase2 .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
            psi_bnd = psi_xpoint(2)
          endif
          if(xcase2 .eq. 1) Z_xpoint(2) = +99.d0
          if(xcase2 .eq. 2) Z_xpoint(1) = -99.d0
        else
          Z_xpoint(1) = -99.d0 
          Z_xpoint(2) = +99.d0
        endif
      endif
  
      if (.not. xpoint2) then
        call find_limiter(my_id,node_list,element_list,bnd_elm_list,psi_lim,R_lim,Z_lim)
        if ( (Z_lim .gt. Z_xpoint(1)) .and. (Z_lim .lt. Z_xpoint(2)) ) then
          call is_axis_psi_mininum(node_list, element_list, bnd_elm_list)
          if (ES%axis_is_psi_minimum) then
            psi_bnd = min(psi_lim,psi_bnd)
          else
            psi_bnd = max(psi_lim,psi_bnd)
          endif
          write(*,'(A,4f8.3)') ' LIMITER PLASMA ',psi_lim, psi_bnd, R_lim,Z_lim
        endif
      endif

      if (freeb_equil_iterate_area .and. (.not. xpoint2)) then
        call iterate2area(node_list,element_list, psi_axis, psi_lim, xpoint2, xcase2, area_ref, psi_bnd)
      endif
      
      write(*,'(A,1f8.3)') ' Psi_bnd = ', psi_bnd   
      
      ! Calculate current feedback
      call integral_current(node_list,element_list,psi_axis, psi_bnd, xpoint2, xcase2, Z_xpoint, current_tot)
  
      current_int = current_int + (current_tot-current_ref)
      
      if (mod(iter,n_feedback_current) .eq. 0) then
        current_FB_fact  = current_FB_fact * (1. - FB_Ip_position * (current_tot-current_ref)/current_ref &
                                                 - FB_Ip_integral *  current_int/current_ref   )
      endif
      
      !-------------- Multiplying FF' and p' profiles by the same factor to scale total current -------------------------
      FF_0 = FF_0_old * current_FB_fact   
      FF_1 = FF_1_old * current_FB_fact      
        
      T_0  = T_0_old  * current_FB_fact    
      T_1  = T_1_old  * current_FB_fact
      !------------------------------------------------------------------------------------------------------------------
      
      write(*,'(A,1e12.4)') 'Current Feedback factor = ',  current_FB_fact
      
      !Vertical feedback - needed for vertically unstable plasmas        
      Z_axis_int = Z_axis_int + (Z_axis - Z_axis_ref)
      if (iter .eq. 1) then
        dZ_axis = 0.d0
      else
        dZ_axis = Z_axis - Z_axis_old
      end if
      
      if ((mod(iter,n_feedback_vertical) .eq. 0) .and. (iter .ge. start_VFB) ) then
        vertical_FB = FB_Zaxis_position   * (Z_axis-Z_axis_ref) &   ! vertical_FB is used in vacuum_equilibrium.f90 to modify the coils current
                    + FB_Zaxis_integral   * Z_axis_int          &   
                    + FB_Zaxis_derivative * dZ_axis
      endif
        
        Z_axis_old = Z_axis
       
    end if ! my_id == 0
    
    call MPI_bcast(vertical_FB, 1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD,ierr)
  
    ! --- Iterate equation
    call poisson(my_id,-1,node_list,element_list,bnd_node_list,bnd_elm_list,3,1,1, &
                 psi_axis,psi_bnd,xpoint2,xcase2,Z_xpoint,freeboundary_equil,refinement,iter)   !----------- for GS use -1
  
  !  call boundary_check
   
    if (my_id == 0) then
      diff = 0.d0
      do i=1, node_list%n_nodes
        diff = diff + abs(node_list%node(i)%deltas(1,1,1))
      enddo  
      diff = diff / float(node_list%n_nodes)
    
      write(*,'(A,i5,e14.6)') ' iteration, diff : ',iter,diff
    end if ! my_id == 0
  
    call MPI_bcast(diff, 1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD,ierr)
  
    if ( (iter > 1) .and. (diff < equil_accuracy_freeb) ) then
      if (my_id == 0) write(*,'(A,I4,A)') ' Free boundary equilibrium converged: after', iter, ' iterations'
      exit
    else if (iter == n_iter_freeb) then
      if (my_id == 0) write(*,'(A,ES10.3)') ' WARNING: Free boundary equilibrium not fully converged: diff=', diff
      exit
    end if
  
  enddo
  
  if (freeb_equil_iterate_area .and. (.not. xpoint2)) then
    n_limiter = 1  ! set found limiter (defined inside iterate2area)
  endif
  
endif

if (my_id == 0) then
  ! Update psi axis and boundary with new values from the last iteration of equilibrium solvers 
  call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
        
  write(10,'(i6,9e20.12)') iter, current_tot, R_axis, Z_axis, psi_bnd-psi_axis
  
  if ((ifail .ne. 0) .and. (iter .le. 5)) then
    call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
    call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
    write(*,*)  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
  endif
  
  psi_bnd = 0.d0
  
  if (xpoint2) then
    call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase2,ifail)
    if (ifail .ne. 1) then      
      psi_bnd  = psi_xpoint(1)
      if( (xcase2 .eq. 2) .or. ((xcase2 .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
        psi_bnd = psi_xpoint(2)
      endif
      if(xcase2 .eq. 1) Z_xpoint(2) = +99.d0
      if(xcase2 .eq. 2) Z_xpoint(1) = -99.d0
    else
      Z_xpoint(1) = -99.d0 
      Z_xpoint(2) = +99.d0
    endif
  endif
  if (.not. xpoint2) then
    call find_limiter(my_id,node_list,element_list,bnd_elm_list,psi_lim,R_lim,Z_lim)
    if ( (Z_lim .gt. Z_xpoint(1)) .and. (Z_lim .lt. Z_xpoint(2)) ) then
      call is_axis_psi_mininum(node_list, element_list, bnd_elm_list)
      if (ES%axis_is_psi_minimum) then
        psi_bnd = min(psi_lim,psi_bnd)
      else
        psi_bnd = max(psi_lim,psi_bnd)
      endif
      write(*,'(A,4f8.3)') ' LIMITER PLASMA ',psi_lim, psi_bnd, R_lim,Z_lim
    endif
  endif

  if (freeboundary_equil .and. freeb_equil_iterate_area .and. (.not. xpoint2)) then
    call iterate2area(node_list,element_list, psi_axis, psi_lim, xpoint2, xcase2, area_ref, psi_bnd)
  endif
  
  !------------------------------- end of equilibrium, start filling data
  psi_axis = psi_axis - psi_offset_freeb
  psi_bnd  = psi_bnd  - psi_offset_freeb
  
  do i=1,node_list%n_nodes
  
    node_list%node(i)%values(1,1,1) = node_list%node(i)%values(1,1,1) - psi_offset_freeb
    psi = node_list%node(i)%values(1,1,1)
    R   = node_list%node(i)%x(1,1,1)
    Z   = node_list%node(i)%x(1,1,2)
  
    call density(    xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zn,             &
                 dn_dpsi,  dn_dz, &                                             ! 1st order derivatives
                 dn_dpsi2, dn_dz2,      dn_dpsi_dz, &                           ! 2nd order derivatives
                 dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,  dn_dz3, &                 ! 2rd order derivatives
                 dn_dpsi4, dn_dpsi_dz3, dn_dpsi2_dz2, dn_dpsi3_dz,  dn_dz4, &   ! 4th order derivatives
                 dn_dpsi5, dn_dpsi_dz4, dn_dpsi2_dz3, dn_dpsi3_dz2, dn_dpsi4_dz)! 5th order derivatives (z5 not needed)
  
    if ( (jorek_model .eq. 400) .or. (jorek_model .eq. 401) .or. (jorek_model .eq. 711) ) then
      call temperature_i(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, zTi, &
                         dTi_dpsi,  dTi_dz, &                                                ! 1st order derivatives
                         dTi_dpsi2, dTi_dz2,      dTi_dpsi_dz, &                             ! 2nd order derivatives
                         dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz,  dTi_dz3, &                  ! 2rd order derivatives
                         dTi_dpsi4, dTi_dpsi_dz3, dTi_dpsi2_dz2, dTi_dpsi3_dz,  dTi_dz4, &   ! 4th order derivatives
                         dTi_dpsi5, dTi_dpsi_dz4, dTi_dpsi2_dz3, dTi_dpsi3_dz2, dTi_dpsi4_dz)! 5th order derivatives (z5 not needed)
  
      call temperature_e(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, zTe, &
                         dTe_dpsi,  dTe_dz, &                                                ! 1st order derivatives
                         dTe_dpsi2, dTe_dz2,      dTe_dpsi_dz, &                             ! 2nd order derivatives
                         dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz,  dTe_dz3, &                  ! 2rd order derivatives
                         dTe_dpsi4, dTe_dpsi_dz3, dTe_dpsi2_dz2, dTe_dpsi3_dz,  dTe_dz4, &   ! 4th order derivatives
                         dTe_dpsi5, dTe_dpsi_dz4, dTe_dpsi2_dz3, dTe_dpsi3_dz2, dTe_dpsi4_dz)! 5th order derivatives (z5 not needed)
      zT  	= zTi	       + zTe
      dT_dpsi	= dTi_dpsi     + dTe_dpsi
      dT_dpsi2	= dTi_dpsi2    + dTe_dpsi2
      dT_dpsi3	= dTi_dpsi3    + dTe_dpsi3
      dT_dz	= dTi_dz       + dTe_dz
      dT_dz2	= dTi_dz2      + dTe_dz2
      dT_dpsi_dz  = dTi_dpsi_dz  + dTe_dpsi_dz
      dT_dpsi2_dz = dTi_dpsi2_dz + dTe_dpsi2_dz
      dT_dpsi_dz2 = dTi_dpsi_dz2 + dTe_dpsi_dz2 
    else
      call temperature(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, zT, &
                       dT_dpsi,  dT_dz, &                                             ! 1st order derivatives
                       dT_dpsi2, dT_dz2,      dT_dpsi_dz, &                           ! 2nd order derivatives
                       dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz,  dT_dz3, &                 ! 2rd order derivatives
                       dT_dpsi4, dT_dpsi_dz3, dT_dpsi2_dz2, dT_dpsi3_dz,  dT_dz4, &   ! 4th order derivatives
                       dT_dpsi5, dT_dpsi_dz4, dT_dpsi2_dz3, dT_dpsi3_dz2, dT_dpsi4_dz)! 5th order derivatives (z5 not needed)
    endif
  
#ifdef fullmhd
    call F_profile(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
                   F_prof, &
                     dF_dpsi, dF_dz, &                                          ! 1st order derivatives
                     dF_dpsi2, dF_dz2, dF_dpsi_dz, &                            ! 2nd order derivatives
                     dF_dpsi3, dF_dpsi_dz2, dF_dpsi2_dz,  dF_dz3, &             ! 2rd order derivatives
                     dF_dpsi4, dF_dpsi_dz3, dF_dpsi2_dz2, dF_dpsi3_dz, dF_dz4, &! 4th order derivatives
                   zFFprime, &
                     dFF_dpsi, dFF_dz, &                                             ! 1st order derivatives
                     dFF_dpsi2, dFF_dz2, dFF_dpsi_dz, &                              ! 2nd order derivatives
                     dFF_dpsi3, dFF_dpsi_dz2, dFF_dpsi2_dz,  dFF_dz3, &              ! 2rd order derivatives
                     dFF_dpsi4, dFF_dpsi_dz3, dFF_dpsi2_dz2, dFF_dpsi3_dz, dFF_dz4)  ! 4th order derivatives
  
    F_values = 0.d0
    F_values(0,0,0) = F_prof
    F_values(1,0,0) = dF_dpsi  ; F_values(0,0,1) = dF_dZ  
    F_values(2,0,0) = dF_dpsi2 ; F_values(0,0,2) = dF_dZ2 
    F_values(3,0,0) = dF_dpsi3 ; F_values(0,0,3) = dF_dZ3 
    F_values(4,0,0) = dF_dpsi4 ; F_values(0,0,4) = dF_dZ4 
    F_values(1,0,1) = dF_dpsi_dZ  ; F_values(2,0,1) = dF_dpsi2_dZ  
    F_values(1,0,2) = dF_dpsi_dZ2 ; F_values(2,0,2) = dF_dpsi2_dZ2
    F_values(1,0,3) = dF_dpsi_dZ3 ; F_values(3,0,2) = dF_dpsi3_dZ

    call project_var_on_node(node_list, i, 710, F_values)
#else
    call FFprime(    xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zFFprime, &
                 dFF_dpsi, dFF_dz, &                                             ! 1st order derivatives
                 dFF_dpsi2, dFF_dz2, dFF_dpsi_dz, &                              ! 2nd order derivatives
                 dFF_dpsi3, dFF_dpsi_dz2, dFF_dpsi2_dz,  dFF_dz3, &              ! 2rd order derivatives
                 dFF_dpsi4, dFF_dpsi_dz3, dFF_dpsi2_dz2, dFF_dpsi3_dz, dFF_dz4, &! 4th order derivatives
                 .true.)

    ! --- Before dealing with the current, let's do grad(p)  
    GP            = dn_dpsi * zT + zn * dT_dpsi
    ! --- 1st derivatives
    dGP_dpsi      = dn_dpsi2 * zT + 2.0 * dn_dpsi * dT_dpsi + zn * dT_dpsi2
    dGP_dz        = dn_dpsi_dZ * zT + dn_dpsi * dT_dZ + dn_dZ * dT_dpsi + zn * dT_dpsi_dZ
    ! --- 2nd derivatives
    dGP_dpsi2     = dn_dpsi3 * zT + 3.0 * dn_dpsi2 * dT_dpsi + 3.0 * dn_dpsi * dT_dpsi2 + zn * dT_dpsi3
    dGP_dz2       = dn_dpsi_dZ2 * zT + 2.0 * dn_dpsi_dZ * dT_dZ + dn_dpsi * dT_dZ2 + dn_dZ2 * dT_dpsi + 2.0 * dn_dZ * dT_dpsi_dZ + zn * dT_dpsi_dZ2
    dGP_dpsi_dz   = dn_dpsi2_dZ * zT + dn_dpsi2 * dT_dZ + 2.0 * dn_dpsi_dZ * dT_dpsi + 2.0 * dn_dpsi * dT_dpsi_dZ + dn_dZ * dT_dpsi2 + zn * dT_dpsi2_dZ
    ! --- 3rd derivatives
    dGP_dpsi3     = dn_dpsi4 * zT  + 4.0 * dn_dpsi3 * dT_dpsi + 6.0 * dn_dpsi2 * dT_dpsi2 &
                                   + 4.0 * dn_dpsi * dT_dpsi3 + zn * dT_dpsi4
    dGP_dpsi_dz2  = dn_dpsi2_dZ2 * zT  + 2.0 * dn_dpsi_dZ2 * dT_dpsi + 2.0 * dn_dpsi2_dZ * dT_dZ &
                    + 4.0 * dn_dpsi_dZ * dT_dpsi_dZ + dn_dpsi2 * dT_dZ2 + 2.0 * dn_dpsi * dT_dpsi_dZ2 &
                    + dn_dZ2 * dT_dpsi2 + 2.0 * dn_dZ * dT_dpsi2_dZ + zn * dT_dpsi2_dZ2
    dGP_dpsi2_dz  = dn_dpsi3_dZ * zT + dn_dpsi3 * dT_dZ + 3.0 * dn_dpsi2_dZ * dT_dpsi &
                    + 3.0 * dn_dpsi2 * dT_dpsi_dZ + 3.0 * dn_dpsi_dZ * dT_dpsi2 &
                    + 3.0 * dn_dpsi * dT_dpsi2_dZ + dn_dZ * dT_dpsi3 + zn * dT_dpsi3_dZ
    dGP_dz3       = dn_dpsi_dZ3 * zT + 3.0 * dn_dpsi_dZ2 * dT_dZ + 3.0 * dn_dpsi_dZ * dT_dZ2 &
                    + dn_dpsi * dT_dZ3 + dn_dZ3 * dT_dpsi + 3.0 * dn_dZ2 * dT_dpsi_dZ &
                    + 3.0 * dn_dZ * dT_dpsi_dZ2 + zn * dT_dpsi_dZ3
    ! --- 4th derivatives
    dGP_dpsi4     = dn_dpsi5 * zT + 5.0 * dn_dpsi4 * dT_dpsi + 10.0 * dn_dpsi3 * dT_dpsi2 &
                    + 10.0 * dn_dpsi2 * dT_dpsi3 + 5.0 * dn_dpsi * dT_dpsi4 + zn * dT_dpsi5
    dGP_dpsi_dz3  = + dn_dpsi2_dZ3 * zT + 2.0 * dn_dpsi_dZ3 * dT_dpsi + 3.0 * dn_dpsi2_dZ2 * dT_dZ &
                    + 6.0 * dn_dpsi_dZ2 * dT_dpsi_dZ + 3.0 * dn_dpsi2_dZ * dT_dZ2 + 6.0 * dn_dpsi_dZ * dT_dpsi_dZ2 &
                    + dn_dpsi2 * dT_dZ3 + 2.0 * dn_dpsi * dT_dpsi_dZ3  + dn_dZ3 * dT_dpsi2 &
                    + 3.0 * dn_dZ2 * dT_dpsi2_dZ + 3.0 * dn_dZ * dT_dpsi2_dZ2 + zn * dT_dpsi2_dZ3
    dGP_dpsi2_dz2 = + dn_dpsi3_dZ2 * zT + 3.0 * dn_dpsi2_dZ2 * dT_dpsi + 3.0 * dn_dpsi_dZ2 * dT_dpsi2 &
                    + 2.0 * dn_dpsi3_dZ * dT_dZ + 6.0 * dn_dpsi2_dZ * dT_dpsi_dZ + 6.0 * dn_dpsi_dZ * dT_dpsi2_dZ &
                    + dn_dpsi3 * dT_dZ2 + 3.0 * dn_dpsi2 * dT_dpsi_dZ2 + 3.0 * dn_dpsi * dT_dpsi2_dZ2 &
                    + dn_dZ2 * dT_dpsi3 + 2.0 * dn_dZ * dT_dpsi3_dZ + zn * dT_dpsi3_dZ2
    dGP_dpsi3_dz  = + dn_dpsi4_dZ * zT + dn_dpsi4 * dT_dZ + 4.0 * dn_dpsi3_dZ * dT_dpsi &
                    + 6.0 * dn_dpsi2_dZ * dT_dpsi2 + 4.0 * dn_dpsi3 * dT_dpsi_dZ &
                    + 4.0 * dn_dpsi_dZ * dT_dpsi3 + 6.0 * dn_dpsi2 * dT_dpsi2_dZ &
                    + 4.0 * dn_dpsi * dT_dpsi3_dZ + dn_dZ * dT_dpsi4 + zn * dT_dpsi4_dZ
    dGP_dz4       = + dn_dpsi_dZ4 * zT + 4.0 * dn_dpsi_dZ3 * dT_dZ + 6.0 * dn_dpsi_dZ2 * dT_dZ2 &
                    + 4.0 * dn_dpsi_dZ * dT_dZ3 + dn_dpsi * dT_dZ4 + dn_dZ4 * dT_dpsi &
                    + 4.0 * dn_dZ3 * dT_dpsi_dZ + 6.0 * dn_dZ2 * dT_dpsi_dZ2 &
                    + 4.0 * dn_dZ * dT_dpsi_dZ3 + zn * dT_dpsi_dZ4


    ! --- Now we  record**2 the current profile for projection on the node
    ! --- For the record**2, j = - FF' - R^2 Grad(p), but JOREK uses -FFprime...
    F_values = 0.d0
    ! --- Current needs to be derived 4*4*4 times...
    F_values(0,0,0) = zFFprime      -       R**2 *  GP
    ! --- 1st derivatives
    F_values(1,0,0) = dFF_dpsi      -       R**2 * dGP_dpsi
    F_values(0,1,0) =               - 2.0 * R    *  GP
    F_values(0,0,1) = dFF_dZ        -       R**2 * dGP_dZ
    ! --- 2nd derivatives
    F_values(2,0,0) = dFF_dpsi2     -       R**2 * dGP_dpsi2
    F_values(0,2,0) =               - 2.0        *  GP
    F_values(0,0,2) = dFF_dZ2       -       R**2 * dGP_dZ2
    F_values(1,1,0) =               - 2.0 * R    * dGP_dpsi
    F_values(1,0,1) = dFF_dpsi_dZ   -       R**2 * dGP_dpsi_dZ
    F_values(0,1,1) =               - 2.0 * R    * dGP_dZ
    ! --- 3rd derivatives
    F_values(3,0,0) = dFF_dpsi3     -       R**2 * dGP_dpsi3
    F_values(0,0,3) = dFF_dZ3       -       R**2 * dGP_dZ3
    F_values(2,1,0) =               - 2.0 * R    * dGP_dpsi2
    F_values(2,0,1) = dFF_dpsi2_dZ  -       R**2 * dGP_dpsi2_dZ
    F_values(1,2,0) =               - 2.0        * dGP_dpsi
    F_values(1,1,1) =               - 2.0 * R    * dGP_dpsi_dZ
    F_values(1,0,2) = dFF_dpsi_dZ2  -       R**2 * dGP_dpsi_dZ2
    F_values(0,2,1) =               - 2.0        * dGP_dZ
    F_values(0,1,2) =               - 2.0 * R    * dGP_dZ2
    ! --- 4th derivatives
    F_values(4,0,0) = dFF_dpsi4     -       R**2 * dGP_dpsi4
    F_values(0,0,4) = dFF_dZ4       -       R**2 * dGP_dZ4
    F_values(3,0,1) = dFF_dpsi3_dZ  -       R**2 * dGP_dpsi3_dZ
    F_values(3,1,0) =               - 2.0 * R    * dGP_dpsi3
    F_values(2,2,0) =               - 2.0        * dGP_dpsi2
    F_values(2,1,1) =               - 2.0 * R    * dGP_dpsi2_dZ
    F_values(2,0,2) = dFF_dpsi2_dZ2 -       R**2 * dGP_dpsi2_dZ2
    F_values(1,2,1) =               - 2.0        * dGP_dpsi_dZ
    F_values(1,1,2) =               - 2.0 * R    * dGP_dpsi_dZ2
    F_values(1,0,3) = dFF_dpsi_dZ3  -       R**2 * dGP_dpsi_dZ3
    F_values(0,2,2) =               - 2.0        * dGP_dZ2
    F_values(0,1,3) =               - 2.0 * R    * dGP_dZ3

    call project_var_on_node(node_list, i, var_zj, F_values)

#endif /* end of non-fullmhd */
  
  enddo
  
  ! --- Find flux surfaces and plot them; determine the q-profile.  
  if (xpoint2 .and. (n_flux .gt. 1)) then
    
    call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
    call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint2,Z_xpoint2,i_elm_xpoint,s_xpoint,t_xpoint,xcase2,ifail)
    
    n_grids = 0
    sigmas  = 0.d0
    
    ! Build up some arrays to send as routine parameters to define_flux_values
    sigmas(1)  = SIG_closed  ; sigmas(2)  = SIG_theta
    sigmas(3)  = SIG_open    ; sigmas(4)  = SIG_outer   ; sigmas(5)  = SIG_inner
    sigmas(6)  = SIG_private ; sigmas(7)  = SIG_up_priv
    sigmas(8)  = SIG_leg_0   ; sigmas(9)  = SIG_leg_1
    sigmas(10) = SIG_up_leg_0; sigmas(11) = SIG_up_leg_1
    sigmas(12) = dPSI_open   ; sigmas(13) = dPSI_outer  ; sigmas(14) = dPSI_inner
    sigmas(15) = dPSI_private; sigmas(16) = dPSI_up_priv
  
    n_grids(1) = 2*n_flux   ; n_grids(2) = n_tht
    n_grids(3) = 2*n_open   ; n_grids(4) = 2*n_outer  ; n_grids(5) = 2*n_inner
    n_grids(6) = 2*n_private; n_grids(7) = 2*n_up_priv
    n_grids(8) = n_leg      ; n_grids(9) = n_up_leg
    if (xcase .eq. 1) then
      n_grids(4) = 0
      n_grids(5) = 0
      n_grids(7) = 0
      n_grids(9) = 0
    endif
    if (xcase .eq. 2) then
      n_grids(4) = 0
      n_grids(5) = 0
      n_grids(6) = 0
      n_grids(8) = 0
    endif
  
    ! Allocate surface_list structure (that's for plotting only)
    if (xcase2 .eq. 1) surface_list%n_psi = 2*n_flux + 2*n_open + 2*n_private
    if (xcase2 .eq. 2) surface_list%n_psi = 2*n_flux + 2*n_open + 2*n_up_priv
    if (xcase2 .eq. 3) surface_list%n_psi = 2*n_flux + 2*n_open + 2*n_outer + 2*n_inner + 2*n_private + 2*n_up_priv
    if (allocated(surface_list%psi_values)) call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
    call tr_allocate(surface_list%psi_values,1,surface_list%n_psi,"surface_list%psi_values",CAT_GRID)
    
    ! Allocate sep_list structure (that's for plotting only)  
    sep_list%n_psi =3
    if(xcase .eq. 3) sep_list%n_psi =6
    if (allocated(sep_list%psi_values)) call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
    call tr_allocate(sep_list%psi_values,1,sep_list%n_psi,"sep_list%psi_values",CAT_GRID)
    
    ! Define the flux values to be plotted...
    psi_axis = psi_axis+0.01 !Just offset a little, because finding surfaces along the side of an element (on the xpoint grid) can be hard...
    call define_flux_values(node_list, element_list, surface_list, sep_list, &
                            xcase2, R_xpoint, Z_xpoint, psi_xpoint, psi_axis, n_grids, sigmas)
    psi_axis = psi_axis-0.01 !Put it back, it's not used anyway, but just for principle!
    
  else
    surface_list%n_psi = 200  
    if (allocated(surface_list%psi_values)) call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
    call tr_allocate(surface_list%psi_values,1,surface_list%n_psi,"surface_list%psi_values",CAT_GRID)
    
    do i = 1, surface_list%n_psi
      surface_list%psi_values(i) = 1.25d0*(float(i)/float(surface_list%n_psi))**2 * (psi_bnd - psi_axis) + psi_axis
    enddo
    
    call find_flux_surfaces(my_id,xpoint2,xcase2,node_list,element_list,surface_list)
  
    sep_list%n_psi =1
    if (allocated(sep_list%psi_values)) call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
    call tr_allocate(sep_list%psi_values,1,sep_list%n_psi,"sep_list%psi_values",CAT_GRID)
    sep_list%psi_values(1) = psi_bnd
  
    call find_flux_surfaces(my_id,xpoint2,xcase2,node_list,element_list,sep_list)
  endif
  
  if (freeboundary_equil) then
    !call plot_coils(.true.)
    call plot_flux_surfaces(node_list,element_list,surface_list,.false.,4,psi_xpoint,R_xpoint,Z_xpoint,xpoint2,xcase2)
    call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,psi_xpoint,R_xpoint,Z_xpoint,xpoint2,xcase2)
  
    call plot_flux_surfaces(node_list,element_list,surface_list,.true.,4,psi_xpoint,R_xpoint,Z_xpoint,xpoint2,xcase2)
    call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,psi_xpoint,R_xpoint,Z_xpoint,xpoint2,xcase2)
    !call plot_coils(.false.)
  else
    if (xpoint2 .and. (n_flux .gt. 1)) then
      call plot_flux_surfaces(node_list,element_list,surface_list,.true.,1,psi_xpoint,R_xpoint,Z_xpoint,xpoint2,xcase2)
      call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,psi_xpoint,R_xpoint,Z_xpoint,xpoint2,xcase2)
    else
      call plot_flux_surfaces(node_list,element_list,surface_list,.true.,1,psi_xpoint,R_xpoint,Z_xpoint,.false.,0)
      call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,psi_xpoint,R_xpoint,Z_xpoint,xpoint2,xcase2)
    endif
  endif
  
  if (nice_q) then
    if (xpoint2) then
      ES%xpoint = xpoint2
      ES%Z_xpoint = Z_xpoint
    endif
    ES%psi_bnd  = psi_bnd
    ES%psi_axis = psi_axis
    ES%Z_xpoint = Z_xpoint
    ES%xpoint   = xpoint
    ES%xcase    = xcase
    call q_profile(node_list,element_list,surface_list,psi_axis,psi_bnd,psi_xpoint,Z_xpoint)
  endif
  
  !================ Temperature and density profiles =f(psi_norm) similar to q(psi_norm) needed to calculate neoclassical coef===========
  if (allocated(T_profile)) call tr_deallocate(T_profile,"T_profile",CAT_GRID)
  call tr_allocate(T_profile,1,surface_list%n_psi,"T_profile",CAT_GRID)
  if (allocated(density_profile)) call tr_deallocate(density_profile,"density_profile",CAT_GRID)
  call tr_allocate(density_profile,1,surface_list%n_psi,"density_profile",CAT_GRID)
  
  do i=2,surface_list%n_psi
    psi= surface_list%psi_values(i)
    
    if ( (jorek_model .eq. 400) .or. (jorek_model .eq. 401) .or. (jorek_model .eq. 711) ) then
      call temperature_i(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, Ti_prof, &
                         dTi_dpsi,  dTi_dz, &                                                ! 1st order derivatives
                         dTi_dpsi2, dTi_dz2,      dTi_dpsi_dz, &                             ! 2nd order derivatives
                         dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz,  dTi_dz3, &                  ! 2rd order derivatives
                         dTi_dpsi4, dTi_dpsi_dz3, dTi_dpsi2_dz2, dTi_dpsi3_dz,  dTi_dz4, &   ! 4th order derivatives
                         dTi_dpsi5, dTi_dpsi_dz4, dTi_dpsi2_dz3, dTi_dpsi3_dz2, dTi_dpsi4_dz)! 5th order derivatives (z5 not needed)
  
      call temperature_e(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, Te_prof, &
                         dTe_dpsi,  dTe_dz, &                                                ! 1st order derivatives
                         dTe_dpsi2, dTe_dz2,      dTe_dpsi_dz, &                             ! 2nd order derivatives
                         dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz,  dTe_dz3, &                  ! 2rd order derivatives
                         dTe_dpsi4, dTe_dpsi_dz3, dTe_dpsi2_dz2, dTe_dpsi3_dz,  dTe_dz4, &   ! 4th order derivatives
                         dTe_dpsi5, dTe_dpsi_dz4, dTe_dpsi2_dz3, dTe_dpsi3_dz2, dTe_dpsi4_dz)! 5th order derivatives (z5 not needed)
      T_prof      = Ti_prof      + Te_prof
      dT_dpsi     = dTi_dpsi     + dTe_dpsi
      dT_dpsi2    = dTi_dpsi2    + dTe_dpsi2
      dT_dpsi3    = dTi_dpsi3    + dTe_dpsi3
      dT_dz       = dTi_dz       + dTe_dz
      dT_dz2      = dTi_dz2      + dTe_dz2
      dT_dpsi_dz  = dTi_dpsi_dz  + dTe_dpsi_dz
      dT_dpsi2_dz = dTi_dpsi2_dz + dTe_dpsi2_dz
      dT_dpsi_dz2 = dTi_dpsi_dz2 + dTe_dpsi_dz2 
    else
      call temperature(.false.,xcase2,0., Z_xpoint, psi,psi_axis,psi_bnd,T_prof,             &
                       dT_dpsi,  dT_dz, &                                             ! 1st order derivatives
                       dT_dpsi2, dT_dz2,      dT_dpsi_dz, &                           ! 2nd order derivatives
                       dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz,  dT_dz3, &                 ! 2rd order derivatives
                       dT_dpsi4, dT_dpsi_dz3, dT_dpsi2_dz2, dT_dpsi3_dz,  dT_dz4, &   ! 4th order derivatives
                       dT_dpsi5, dT_dpsi_dz4, dT_dpsi2_dz3, dT_dpsi3_dz2, dT_dpsi4_dz)! 5th order derivatives (z5 not needed)
    endif

    call density( .false., xcase2,0., Z_xpoint, psi,psi_axis,psi_bnd,density_prof,             &
                 dn_dpsi,  dn_dz, &                                             ! 1st order derivatives
                 dn_dpsi2, dn_dz2,      dn_dpsi_dz, &                           ! 2nd order derivatives
                 dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,  dn_dz3, &                 ! 2rd order derivatives
                 dn_dpsi4, dn_dpsi_dz3, dn_dpsi2_dz2, dn_dpsi3_dz,  dn_dz4, &   ! 4th order derivatives
                 dn_dpsi5, dn_dpsi_dz4, dn_dpsi2_dz3, dn_dpsi3_dz2, dn_dpsi4_dz)! 5th order derivatives (z5 not needed)
  
    T_profile(i)=T_prof
    density_profile(i)=density_prof
  end do
  
  write(*,*) '***************************************'
  write(*,*) 'output T and rho profiles (in JOREK units) for neoclassical profile calculation'
  ! --- Write out T and rho profiles to "T_rho_profiles.dat".
  open(432, file='T_rho_profiles.dat', action='write', status='replace')
  do i=2, surface_list%n_psi
     write(432,'(3ES13.5)') T_profile(i), density_profile(i)
  end do
  close(432)
  !========================= end modif ===========================================
  
  if (allocated(surface_list%psi_values))    call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
  if (allocated(surface_list%flux_surfaces)) deallocate(surface_list%flux_surfaces)
  if (allocated(sep_list%psi_values))        call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
  if (allocated(sep_list%flux_surfaces))     deallocate(sep_list%flux_surfaces)
  
  if (allocated(T_profile)) call tr_deallocate(T_profile,"T_profile",CAT_GRID)
  if (allocated(density_profile)) call tr_deallocate(density_profile,"density_profile",CAT_GRID)
  
end if ! my_id == 0

if (freeboundary_equil) then
  call broadcast_elements(my_id, element_list)
  call broadcast_nodes(my_id, node_list)  !--- This is required for boundary_check
  call broadcast_boundary(my_id, bnd_elm_list, bnd_node_list)
  call boundary_check(my_id)
endif
return
end subroutine equilibrium
