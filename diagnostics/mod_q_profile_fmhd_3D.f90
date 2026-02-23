module mod_q_profile_fmhd_3D

use constants
use mod_parameters
use data_structure
use equil_info
use phys_module
use mod_interp
use mod_find_rz_nearby

implicit none

!--- Input parameters --------------------!
!-----------------------------------------!
real*8, parameter  :: stepsizePhi = 2 * PI / 8192 !1.d-3
real*8, parameter  :: qmax = 11.d0
real*8, parameter  :: qmin = 0.5d0
integer, parameter :: polturns_stop = 100
integer, parameter :: torturns_stop = 1000
integer, parameter :: reference_n_mode_number = 1
!-----------------------------------------!
!-----------------------------------------!
integer   :: npoints
integer   :: i, j, k, ielm, ifail, ielm_new, ielm_old
integer   :: polturns, torturns
logical   :: stop_tracing

real*8    :: s, t, phi, phinew, phiold, snew, tnew, sold, told
real*8    :: R, R_s, R_t, Z, Z_s, Z_t
real*8    :: BR, BZ, Bp, BB, xjac, Fprof
real*8    :: dum01, dum02, dum03, dum04, dum05
real*8    :: AR, AR_p, AR_s, AR_t, AR_R, AR_Z, AZ, AZ_p, AZ_s, AZ_t, AZ_R, AZ_Z, A3, A3_p, A3_s, A3_t, A3_R, A3_Z, psieq
real*8    :: RR, ZZ, Rnew, Znew, Rold, Zold, A3old
real*8    :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt, Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
integer   :: i_tor, i_harm
integer   :: stepPhiMax, num_poinc_max
real*8    :: psi_start, psi_start_norm, ratio_phi

real*8 :: BR_ref, BZ_ref, Bp_ref, Bp_initial, BB_ref
real*8 :: uR_all, uZ_all, up_all, VXB_para_all


logical, allocatable :: responsible(:)
integer, allocatable :: num_poinc_min(:)

real*8, allocatable :: Rstart(:), Zstart(:)
real*8, allocatable :: phi_arr(:), Raxis3D(:), Zaxis3D(:)
real*8, allocatable :: psi_arr_poinc(:,:,:)
real*8, allocatable :: R_arr_poinc(:,:,:)  
real*8, allocatable :: Z_arr_poinc(:,:,:)  
real*8, allocatable :: q_arr_poinc(:,:,:)  
real*8, allocatable :: ne_arr_poinc(:,:,:) 
real*8, allocatable :: Ti_arr_poinc(:,:,:) 
real*8, allocatable :: Te_arr_poinc(:,:,:) 
real*8, allocatable :: VxBpara_poinc(:,:,:) 

contains

subroutine construct_all(num_points, nphi_out)
  integer, intent(in) :: num_points, nphi_out

  npoints = num_points
  stepPhiMax = max(nint(2*PI*qmax/stepsizePhi) * polturns_stop, nint(2*PI*2/stepsizePhi) * torturns_stop)
  num_poinc_max = max(nint(qmax * polturns_stop), torturns_stop) * 2
  allocate(num_poinc_min(num_points))
  allocate(Rstart(num_points))
  allocate(Zstart(num_points))
  allocate(responsible(num_points))
  allocate(phi_arr(num_points))

  allocate(Raxis3D(stepPhiMax))
  allocate(Zaxis3D(stepPhiMax))

  allocate(psi_arr_poinc(num_poinc_max, num_points, nphi_out))
  allocate(R_arr_poinc(num_poinc_max, num_points, nphi_out))
  allocate(Z_arr_poinc(num_poinc_max, num_points, nphi_out))
  allocate(q_arr_poinc(num_poinc_max, num_points, nphi_out))
  allocate(ne_arr_poinc(num_poinc_max, num_points, nphi_out))
  allocate(Ti_arr_poinc(num_poinc_max, num_points, nphi_out))
  allocate(Te_arr_poinc(num_poinc_max, num_points, nphi_out))
  allocate(VxBpara_poinc(num_poinc_max, num_points, nphi_out))

  responsible = .false.
  Rstart = 0.d0
  Zstart = 0.d0
  Raxis3D = 0.d0
  Zaxis3D = 0.d0
  num_poinc_min = num_poinc_max

  call initia_arrs()

end subroutine construct_all

subroutine deconstruct_all()

  deallocate(num_poinc_min)
  deallocate(Rstart)
  deallocate(Zstart)
  deallocate(responsible)
  deallocate(phi_arr)
  deallocate(Raxis3D)
  deallocate(Zaxis3D)
  deallocate(psi_arr_poinc)
  deallocate(R_arr_poinc)
  deallocate(Z_arr_poinc)
  deallocate(q_arr_poinc)
  deallocate(ne_arr_poinc)
  deallocate(Ti_arr_poinc)
  deallocate(Te_arr_poinc)
  deallocate(VxBpara_poinc)

end subroutine deconstruct_all

subroutine initia_arrs()

  phi_arr = 0.d0
  psi_arr_poinc = 0.d0
  R_arr_poinc = 0.d0
  Z_arr_poinc = 0.d0
  q_arr_poinc = 0.d0
  ne_arr_poinc = 0.d0
  Ti_arr_poinc = 0.d0
  Te_arr_poinc = 0.d0
  VxBpara_poinc = 0.d0
end subroutine initia_arrs

pure function cross_product3D(a, b)
  real*8 :: cross_product3D
  real*8, dimension(2), intent(in) :: a, b

  cross_product3D = a(1) * b(2) - a(2) * b(1)
end function cross_product3D

subroutine trance_field_lines_fmhd_3D (node_list, element_list, phi_start, phi_poinc_record, nphi_out)
  type (type_node_list),    intent(in)    :: node_list
  type (type_element_list), intent(in)    :: element_list
  real*8, intent(in) :: phi_start, phi_poinc_record(:)
  integer, intent(in):: nphi_out
  
  real*8, parameter :: deltaFactor = 1.0
  real*8 :: RZstart_plus(2), RZstart_minus(2)
  real*8 :: RZ_now(2), RZ_old(2)
  real*8 :: vec_center(2), vec_left(2), vec_right(2)
  real*8 :: vec_poloidal(2), vec_above(2), vec_below(2)
  real*8 :: LRZ, LR, LZ
  
  real*8 :: cross_product_cl, cross_product_cr
  real*8 :: cross_product_pa, cross_product_pb
  real*8 :: stepPoloidal, deltaPoloidal
  integer:: minJsteps
  logical:: is_intersected, is_overlapped, is_thesame
  real*8 :: distanceToStartNow, distanceToStartOld, distToAxis
  real*8 :: RaxisStart, ZaxisStart, cosThetaStart, sinThetaStart
  real*8 :: RaxisNow, ZaxisNow, cosThetaNow, sinThetaNow
  real*8 :: RaxisOld, ZaxisOld, cosThetaOld, sinThetaOld

  integer :: jstep_repeat, jphi_start, iphi_out

  call initia_arrs()

  ! --- Trace field lines
  do i = 1, npoints
    if ( .not. responsible(i) ) then
      cycle
    end if
    
    
    call find_RZ(node_list,element_list,Rstart(i),Zstart(i),RR,ZZ,ielm,s,t,ifail)

    call interp(node_list, element_list, ielm, 1, 1, s, t, psi_start, dum01, dum02, dum03, dum04, dum05)
    psi_start_norm = get_psi_n(psi_start, ZZ)

    phi = phi_start
    jphi_start = nint(phi_start / stepsizePhi) + 1
    
    if ( ielm < 1 ) then
      write(*,*) 'Illegal starting point for field line ', i, ', Rstart = ', Rstart(i), ', Zstart = ', Zstart(i), '. Skipping.'
      cycle
    end if
    
    stop_tracing = .false.
    polturns     = 0
    torturns     = 0
    j            = 1
    jstep_repeat = 1

    distanceToStartNow = 1.d5
    distanceToStartOld = 1.d9

    ! Here just assume q0 will not below qmin, it can be adjusted
    minJsteps = nint(2.d0*PI/stepsizePhi*qmin)
    
    do while( .not. stop_tracing )

      call do_step_3D(node_list, element_list)
      j = j + 1
      jstep_repeat = jstep_repeat + 1

      if ( DABS(phi) >= 2.d0*PI ) then
        phi   = phi - 2.d0*PI * sign (1.d0, phi)
        torturns = torturns + 1
      endif
    
      do iphi_out = 1, nphi_out
        if ( (phi - phi_poinc_record(iphi_out)) * (phi - stepsizePhi - phi_poinc_record(iphi_out)) <= 0 ) then

          if(torturns < num_poinc_max) then
            ! record for poincare plots
            ratio_phi = 1.d0 - abs(phi - phi_poinc_record(iphi_out)) / stepsizePhi

            call determine_field_3D(node_list,element_list) ! update A3
            psi_arr_poinc(torturns + 1, i, iphi_out) = ratio_phi * get_psi_n(A3, ZZ) + (1.d0 - ratio_phi) * get_psi_n(A3old, Zold)
            R_arr_poinc(torturns + 1, i, iphi_out) = ratio_phi * RR + (1.d0 - ratio_phi) * Rold
            Z_arr_poinc(torturns + 1, i, iphi_out) = ratio_phi * ZZ + (1.d0 - ratio_phi) * Zold

            call interp_values(node_list, element_list, R_arr_poinc(torturns + 1, i, iphi_out), Z_arr_poinc(torturns + 1, i, iphi_out), phi_poinc_record(iphi_out), var_rho, ne_arr_poinc(torturns + 1, i, iphi_out))
            call interp_values(node_list, element_list, R_arr_poinc(torturns + 1, i, iphi_out), Z_arr_poinc(torturns + 1, i, iphi_out), phi_poinc_record(iphi_out), var_Te, Te_arr_poinc(torturns + 1, i, iphi_out))
            call interp_values(node_list, element_list, R_arr_poinc(torturns + 1, i, iphi_out), Z_arr_poinc(torturns + 1, i, iphi_out), phi_poinc_record(iphi_out), var_Ti, Ti_arr_poinc(torturns + 1, i, iphi_out))

            call interp_values(node_list, element_list, R_arr_poinc(torturns + 1, i, iphi_out), Z_arr_poinc(torturns + 1, i, iphi_out), phi_poinc_record(iphi_out), var_uR, uR_all)
            call interp_values(node_list, element_list, R_arr_poinc(torturns + 1, i, iphi_out), Z_arr_poinc(torturns + 1, i, iphi_out), phi_poinc_record(iphi_out), var_uZ, uZ_all)
            call interp_values(node_list, element_list, R_arr_poinc(torturns + 1, i, iphi_out), Z_arr_poinc(torturns + 1, i, iphi_out), phi_poinc_record(iphi_out), var_up, up_all)

            VXB_para_all = ((uZ_all * ( Bp - Bp_ref ) - up_all * ( BZ - BZ_ref )) * BR_ref + &
                            (up_all * ( BR - BR_ref ) - uR_all * ( Bp - Bp_ref )) * BZ_ref + &
                            (uR_all * ( BZ - BZ_ref ) - uZ_all * ( BR - BR_ref )) * Bp_ref) / BB_ref !Bp_initial

            VxBpara_poinc(torturns + 1, i, iphi_out) = VXB_para_all

          else
            write(*,*) "Warning: torturns is larger than num_poinc_max"
            write(*,*) "torturns = ", torturns
            write(*,*) "num_poinc_max = ", num_poinc_max
          endif

          exit

        end if
      enddo

      if(jstep_repeat > minJsteps) then

        ! For grids far from axis, just use the axis of n = 0 field.
        ! Otherwise near the x-point, the Bp is too weak and error in q tracing will be large
        if(psi_start_norm .ge. 0.8) then
          RaxisStart = ES%R_axis
          ZaxisStart = ES%Z_axis
          RaxisNow = ES%R_axis
          ZaxisNow = ES%Z_axis
          RaxisOld = ES%R_axis
          ZaxisOld = ES%Z_axis
        else
          RaxisStart = Raxis3D(jphi_start)
          ZaxisStart = Zaxis3D(jphi_start)
          RaxisNow = Raxis3D(j + jphi_start - 1)
          ZaxisNow = Zaxis3D(j + jphi_start - 1)
          RaxisOld = Raxis3D(j - 1 + jphi_start - 1)
          ZaxisOld = Zaxis3D(j - 1 + jphi_start - 1)
        endif

        distToAxis = dsqrt((Rstart(i)-RaxisStart)**2+(Zstart(i)-ZaxisStart)**2)
        cosThetaStart = (Rstart(i)-RaxisStart) / distToAxis
        sinThetaStart = (Zstart(i)-ZaxisStart) / distToAxis

        distToAxis = dsqrt((RR-RaxisNow)**2+(ZZ-ZaxisNow)**2)
        cosThetaNow = (RR-RaxisNow) / distToAxis
        sinThetaNow = (ZZ-ZaxisNow) / distToAxis

        distToAxis = dsqrt((Rold-RaxisOld)**2+(Zold-ZaxisOld)**2)
        cosThetaOld = (Rold-RaxisOld) / distToAxis
        sinThetaOld = (Zold-ZaxisOld) / distToAxis

        ! these two segments below should intersect with each other
        ! (0, 0)                     -> (cosThetaStart, sinThetaStart)
        ! (cosThetaOld, sinThetaOld) -> (cosThetaNow, sinThetaNow)

        RZstart_minus(1) = cosThetaStart*(1.d0 - deltaFactor)
        RZstart_minus(2) = sinThetaStart*(1.d0 - deltaFactor)
        RZstart_plus(1) = cosThetaStart*(1.d0 + deltaFactor)
        RZstart_plus(2) = sinThetaStart*(1.d0 + deltaFactor)

        RZ_now(1) = cosThetaNow
        RZ_now(2) = sinThetaNow
        RZ_old(1) = cosThetaOld
        RZ_old(2) = sinThetaOld

        vec_center(:) = RZstart_plus(:) - RZstart_minus(:)
        vec_left(:) = RZ_now(:) - RZstart_minus(:)
        vec_right(:) = RZ_old(:) - RZstart_minus(:)
        cross_product_cl = cross_product3D(vec_center, vec_left)
        cross_product_cr = cross_product3D(vec_center, vec_right)
  
        vec_poloidal(:) = RZ_now(:) - RZ_old(:)
        vec_above(:) = RZstart_plus(:) - RZ_old(:)
        vec_below(:) = RZstart_minus(:) - RZ_old(:)
        cross_product_pa = cross_product3D(vec_poloidal, vec_above)
        cross_product_pb = cross_product3D(vec_poloidal, vec_below)

        is_intersected = .false.
        is_overlapped  = .false.
        is_thesame     = .false.

        distanceToStartNow = dsqrt((cosThetaNow-cosThetaStart)**2+(sinThetaNow-sinThetaStart)**2)
        if(distanceToStartNow.le.1.d-10) is_thesame = .true.

        if(.not.is_thesame) then
          if( ((dabs(cross_product_cl) + dabs(cross_product_cr) + dabs(cross_product_pa) + dabs(cross_product_pb)).le.1.d-10) .and. & 
          (((cosThetaOld-RZstart_minus(1))*(cosThetaOld-RZstart_plus(1)).le.0.d0).or.((cosThetaNow-RZstart_minus(1))*(cosThetaNow-RZstart_plus(1)).le.0.d0))) is_overlapped = .true.
        endif

        if(.not.is_overlapped) then
          if((cross_product_cl * cross_product_cr <= 0.d0) .and. (cross_product_pa * cross_product_pb <= 0.d0)) is_intersected = .true.
        endif

        if (is_intersected .or. is_thesame) then
          polturns = polturns + 1
          jstep_repeat = 1
          if((polturns.ge.polturns_stop).and.(torturns.ge.torturns_stop)) stop_tracing = .true.
          if(mod(i, npoints/10) .eq. 0) then
            write(*,*) 'tracing line ', i, ' of ', npoints
            write(*,*) 'polturns, polturns_stop = ', polturns, polturns_stop
            write(*,*) 'torturns, torturns_stop = ', torturns, torturns_stop
          endif
        elseif(is_overlapped) then
          write(*,*) 'two parallel segments are overlaped when tracing line ', i, ' of ', npoints
          write(*,*) '!!! toroidal & poloidal turns:', torturns, polturns
          write(*,*) '!!! Rstart:', Rstart(i), ', R_axis:', ES%R_axis
          write(*,*) '!!! Zstart:', Zstart(i), ', Z_axis:', ES%Z_axis
          write(*,*) '!!! Rold:', Rold, ', Rnow:', RR
          write(*,*) '!!! Zold:', Zold, ', Znow:', ZZ
          write(*,*) '!!! Phiold:', phiold, ', Phinow:', phi
          polturns = polturns + 1
          jstep_repeat = 1
          if((polturns.ge.polturns_stop).and.(torturns.ge.torturns_stop)) stop_tracing = .true.

        ! elseif(torturns > 10 * max(torturns_stop, polturns_stop)) then
        elseif(j > stepPhiMax) then
          write(*,*) 'more than 10 * max(torturns_stop, polturns) circles when tracing line ', i, ' of ', npoints
          polturns = polturns + 1
          jstep_repeat = 1
          if(polturns.ge.polturns_stop) stop_tracing = .true.
        endif  

      endif
      
    end do

    stepPoloidal = dsqrt((cosThetaNow-cosThetaOld)**2+(sinThetaNow-sinThetaOld)**2)
    deltaPoloidal = dsqrt((cosThetaStart-cosThetaOld)**2+(sinThetaStart-sinThetaOld)**2)

    if ( phi > 0.d0 ) then
      phi_arr(i)      = phiold + (phi - phiold) * (deltaPoloidal/stepPoloidal) + torturns * 2.d0 * PI
    else if ( phi < 0.d0 ) then
      phi_arr(i)      = phiold + (phi - phiold) * (deltaPoloidal/stepPoloidal) - torturns * 2.d0 * PI
    else
      if ( phiold > 0.d0 ) then
        phi_arr(i)      = phiold + (phi - phiold) * (deltaPoloidal/stepPoloidal) + torturns * 2.d0 * PI
      else
        phi_arr(i)      = phiold + (phi - phiold) * (deltaPoloidal/stepPoloidal) - torturns * 2.d0 * PI
      end if
    end if
    phi_arr(i) = phi_arr(i) / polturns
    if(mod(i, npoints/10) .eq. 0) then
      write(*,*) 'Finished tracing line ', i, ' of ', npoints
      write(*,*) '...toroidal & poloidal turns:', torturns, polturns
      write(*,*) '...phi_arr(i)/2pi:', phi_arr(i)/2/PI
      write(*,*) '...Rstart:', Rstart(i), ', R_axis:', ES%R_axis
      write(*,*) '...Zstart:', Zstart(i), ', Z_axis:', ES%Z_axis
      write(*,*) '...Rold:', Rold, ', Rnow:', RR
      write(*,*) '...Zold:', Zold, ', Znow:', ZZ
      write(*,*) '...Phiold:', phiold, ', Phinow:', phi
    endif

    q_arr_poinc(:,i, :) = phi_arr(i) / (2.d0 * PI)

  end do
  end subroutine trance_field_lines_fmhd_3D

  !> A "stupid" stepping routine advancing in R, Z, phi
  subroutine do_step_3D(node_list,element_list)

    type (type_node_list),    intent(in)    :: node_list
    type (type_element_list), intent(in)    :: element_list
    real*8    :: Rnewtmp, Znewtmp, RRtmp, ZZtmp

    call determine_field_3D(node_list,element_list)
    
    A3old    = A3
    Rold     = RR
    Zold     = ZZ
    phiold   = phi
    sold     = s
    told     = t
    ielm_old = ielm
    
    Rnew   = RR  + 0.5d0 * stepsizePhi * R * BR / Bp
    Znew   = ZZ  + 0.5d0 * stepsizePhi * R * BZ / Bp
    phinew = phi + 0.5d0 * stepsizePhi
    
    call find_RZ_nearby(node_list, element_list, Rold, Zold, sold, told, ielm_old, &
                        Rnew, Znew, snew, tnew, ielm_new, ifail)
    if ( ielm_new < 1 ) then
      ! write(*,*) 'find_RZ_nearby failed, using find_RZ instead, (R,Z)-position:', Rnew, Znew
      Rnewtmp = Rnew
      Znewtmp = Znew 
      call find_RZ(node_list, element_list, Rnewtmp, Znewtmp, Rnew, Znew, ielm_new, snew, tnew, ifail)
    end if 
    if ( ielm_new < 1 ) then
      ! write(*,*) 'ATTENTION: both find_RZ_nearby and find_RZ failed!', Rnew, Znew, ielm_new, ifail
      stop_tracing = .true.
      return
    end if
    
    RR   = Rnew
    ZZ   = Znew
    s    = snew
    t    = tnew
    ielm = ielm_new
    phi  = phinew
  
    call determine_field_3D(node_list,element_list)
    
    RR   = Rold   + stepsizePhi * R * BR / Bp
    ZZ   = Zold   + stepsizePhi * R * BZ / Bp
    phi  = phiold + stepsizePhi
    
    call find_RZ_nearby(node_list, element_list, Rold, Zold, sold, told, ielm_old, &
                        RR, ZZ, s, t, ielm, ifail)
    if ( ielm < 1 ) then
      ! write(*,*) 'find_RZ_nearby failed, using find_RZ instead, (R,Z)-position:', RR, ZZ
      RRtmp = RR
      ZZtmp = ZZ 
      call find_RZ(node_list, element_list, RRtmp, ZZtmp, RR, ZZ, ielm, s, t, ifail)
    end if
    if ( ielm < 1 ) then
      ! write(*,*) 'ATTENTION: both find_RZ_nearby and find_RZ failed!', RR, ZZ, ielm, ifail
      stop_tracing = .true.
      return
    end if
  
  end subroutine do_step_3D


  subroutine interp_values(node_list,element_list,R0,Z0,phi0,ivar,result_ivar)

    type (type_node_list),    intent(in)    :: node_list
    type (type_element_list), intent(in)    :: element_list
    real*8, intent(in)  :: R0, Z0, phi0
    integer, intent(in) :: ivar
    real*8, intent(out) :: result_ivar
    real*8              :: Rtmp, Ztmp

    call find_RZ(node_list,element_list,R0,Z0,Rtmp,Ztmp,ielm,s,t,ifail)

    call interp(node_list, element_list, ielm, ivar, 1, s, t, result_ivar, dum01, dum02, dum03, dum04, dum05)

    do i_tor = 1, (n_tor-1)/2
      i_harm = 2*i_tor
      call interp(node_list, element_list, ielm, ivar, i_harm, s, t, Pcos, dum01, dum02, dum03, dum04, dum05)
      call interp(node_list, element_list, ielm, ivar, i_harm + 1, s, t, Psin, dum01, dum02, dum03, dum04, dum05)
      result_ivar = result_ivar + Pcos * cos(mode(i_harm)*phi0) + Psin * sin(mode(i_harm + 1)*phi0)
    enddo
    
  end subroutine interp_values

  subroutine interp_values_n_ne_ref_and_0(node_list,element_list,R0,Z0,phi0,ivar,result_ivar)

    type (type_node_list),    intent(in)    :: node_list
    type (type_element_list), intent(in)    :: element_list
    real*8, intent(in)  :: R0, Z0, phi0
    integer, intent(in) :: ivar
    real*8, intent(out) :: result_ivar
    real*8              :: Rtmp, Ztmp

    call find_RZ(node_list,element_list,R0,Z0,Rtmp,Ztmp,ielm,s,t,ifail)

    ! call interp(node_list, element_list, ielm, ivar, 1, s, t, result_ivar, dum01, dum02, dum03, dum04, dum05)

    result_ivar = 0.d0
    do i_tor = 1, (n_tor-1)/2
      if (i_tor .ne. reference_n_mode_number) then
        i_harm = 2*i_tor
        call interp(node_list, element_list, ielm, ivar, i_harm, s, t, Pcos, dum01, dum02, dum03, dum04, dum05)
        call interp(node_list, element_list, ielm, ivar, i_harm + 1, s, t, Psin, dum01, dum02, dum03, dum04, dum05)
        result_ivar = result_ivar + Pcos * cos(mode(i_harm)*phi0) + Psin * sin(mode(i_harm + 1)*phi0)
      endif
    enddo
    
  end subroutine interp_values_n_ne_ref_and_0


  subroutine interp_values_n_eq_ref_and_0(node_list,element_list,R0,Z0,phi0,ivar,result_ivar)

    type (type_node_list),    intent(in)    :: node_list
    type (type_element_list), intent(in)    :: element_list
    real*8, intent(in)  :: R0, Z0, phi0
    integer, intent(in) :: ivar
    real*8, intent(out) :: result_ivar
    real*8              :: Rtmp, Ztmp

    call find_RZ(node_list,element_list,R0,Z0,Rtmp,Ztmp,ielm,s,t,ifail)

    call interp(node_list, element_list, ielm, ivar, 1, s, t, result_ivar, dum01, dum02, dum03, dum04, dum05)

    do i_tor = 1, (n_tor-1)/2
      if (i_tor .eq. reference_n_mode_number) then
        i_harm = 2*i_tor
        call interp(node_list, element_list, ielm, ivar, i_harm, s, t, Pcos, dum01, dum02, dum03, dum04, dum05)
        call interp(node_list, element_list, ielm, ivar, i_harm + 1, s, t, Psin, dum01, dum02, dum03, dum04, dum05)
        result_ivar = result_ivar + Pcos * cos(mode(i_harm)*phi0) + Psin * sin(mode(i_harm + 1)*phi0)
      endif
    enddo
    
  end subroutine interp_values_n_eq_ref_and_0


  subroutine interp_values_n_eq_0(node_list,element_list,R0,Z0,phi0,ivar,result_ivar)

    type (type_node_list),    intent(in)    :: node_list
    type (type_element_list), intent(in)    :: element_list
    real*8, intent(in)  :: R0, Z0, phi0
    integer, intent(in) :: ivar
    real*8, intent(out) :: result_ivar
    real*8              :: Rtmp, Ztmp

    call find_RZ(node_list,element_list,R0,Z0,Rtmp,Ztmp,ielm,s,t,ifail)

    call interp(node_list, element_list, ielm, ivar, 1, s, t, result_ivar, dum01, dum02, dum03, dum04, dum05)
    
  end subroutine interp_values_n_eq_0
 
 !> Determine magnetic field at given position
  subroutine determine_field_3D(node_list,element_list)

    type (type_node_list),    intent(in)    :: node_list
    type (type_element_list), intent(in)    :: element_list

    real*8 :: A3_s_ref, A3_t_ref, A3_p_ref
    real*8 :: AR_s_ref, AR_t_ref, AR_p_ref
    real*8 :: AZ_s_ref, AZ_t_ref, AZ_p_ref
    real*8 :: AR_R_ref, AR_Z_ref, AZ_R_ref, AZ_Z_ref, A3_R_ref, A3_Z_ref

    A3_p = 0.d0
    AR_p = 0.d0
    AZ_p = 0.d0

    A3_s_ref = 0.d0; A3_t_ref = 0.d0; A3_p_ref = 0.d0
    AR_s_ref = 0.d0; AR_t_ref = 0.d0; AR_p_ref = 0.d0
    AZ_s_ref = 0.d0; AZ_t_ref = 0.d0; AZ_p_ref = 0.d0

    call interp(node_list, element_list, ielm, 1, 1, s, t, A3, A3_s, A3_t, dum01, dum02, dum03)
    call interp(node_list, element_list, ielm, 2, 1, s, t, AR, AR_s, AR_t, dum01, dum02, dum03)
    call interp(node_list, element_list, ielm, 3, 1, s, t, AZ, AZ_s, AZ_t, dum01, dum02, dum03)

    A3_s_ref = A3_s; A3_t_ref = A3_t
    AR_s_ref = AR_s; AR_t_ref = AR_t
    AZ_s_ref = AZ_s; AZ_t_ref = AZ_t

    do i_tor = 1, (n_tor-1)/2
      i_harm = 2*i_tor
      call interp(node_list, element_list, ielm, 1, i_harm, s, t, Pcos, Pcos_s, Pcos_t, Pcos_st, Pcos_ss, Pcos_tt)
      call interp(node_list, element_list, ielm, 1, i_harm + 1, s, t, Psin, Psin_s, Psin_t, Psin_st, Psin_ss, Psin_tt)
      A3_s = A3_s + Pcos_s * cos(mode(i_harm)*phi) + Psin_s * sin(mode(i_harm + 1)*phi)
      A3_t = A3_t + Pcos_t * cos(mode(i_harm)*phi) + Psin_t * sin(mode(i_harm + 1)*phi)
      A3_p = A3_p - Pcos   * sin(mode(i_harm)*phi) * mode(i_harm) + Psin   * cos(mode(i_harm+1)*phi) * mode(i_harm+1)

      if ((reference_n_mode_number > 0) .and. (i_tor == reference_n_mode_number)) then
        A3_s_ref = A3_s_ref + Pcos_s * cos(mode(i_harm)*phi) + Psin_s * sin(mode(i_harm + 1)*phi)
        A3_t_ref = A3_t_ref + Pcos_t * cos(mode(i_harm)*phi) + Psin_t * sin(mode(i_harm + 1)*phi)
        A3_p_ref = A3_p_ref - Pcos   * sin(mode(i_harm)*phi) * mode(i_harm) + Psin   * cos(mode(i_harm+1)*phi) * mode(i_harm+1)
      endif

      call interp(node_list, element_list, ielm, 2, i_harm, s, t, Pcos, Pcos_s, Pcos_t, Pcos_st, Pcos_ss, Pcos_tt)
      call interp(node_list, element_list, ielm, 2, i_harm + 1, s, t, Psin, Psin_s, Psin_t, Psin_st, Psin_ss, Psin_tt)
      AR_s = AR_s + Pcos_s * cos(mode(i_harm)*phi) + Psin_s * sin(mode(i_harm + 1)*phi)
      AR_t = AR_t + Pcos_t * cos(mode(i_harm)*phi) + Psin_t * sin(mode(i_harm + 1)*phi)
      AR_p = AR_p - Pcos   * sin(mode(i_harm)*phi) * mode(i_harm) + Psin   * cos(mode(i_harm+1)*phi) * mode(i_harm+1)

      if ((reference_n_mode_number > 0) .and. (i_tor == reference_n_mode_number)) then
        AR_s_ref = AR_s_ref + Pcos_s * cos(mode(i_harm)*phi) + Psin_s * sin(mode(i_harm + 1)*phi)
        AR_t_ref = AR_t_ref + Pcos_t * cos(mode(i_harm)*phi) + Psin_t * sin(mode(i_harm + 1)*phi)
        AR_p_ref = AR_p_ref - Pcos   * sin(mode(i_harm)*phi) * mode(i_harm) + Psin   * cos(mode(i_harm+1)*phi) * mode(i_harm+1)
      endif

      call interp(node_list, element_list, ielm, 3, i_harm, s, t, Pcos, Pcos_s, Pcos_t, Pcos_st, Pcos_ss, Pcos_tt)
      call interp(node_list, element_list, ielm, 3, i_harm + 1, s, t, Psin, Psin_s, Psin_t, Psin_st, Psin_ss, Psin_tt)
      AZ_s = AZ_s + Pcos_s * cos(mode(i_harm)*phi) + Psin_s * sin(mode(i_harm + 1)*phi)
      AZ_t = AZ_t + Pcos_t * cos(mode(i_harm)*phi) + Psin_t * sin(mode(i_harm + 1)*phi)
      AZ_p = AZ_p - Pcos   * sin(mode(i_harm)*phi) * mode(i_harm) + Psin   * cos(mode(i_harm+1)*phi) * mode(i_harm+1)

      if ((reference_n_mode_number > 0) .and. (i_tor == reference_n_mode_number)) then
        AZ_s_ref = AZ_s_ref + Pcos_s * cos(mode(i_harm)*phi) + Psin_s * sin(mode(i_harm + 1)*phi)
        AZ_t_ref = AZ_t_ref + Pcos_t * cos(mode(i_harm)*phi) + Psin_t * sin(mode(i_harm + 1)*phi)
        AZ_p_ref = AZ_p_ref - Pcos   * sin(mode(i_harm)*phi) * mode(i_harm) + Psin   * cos(mode(i_harm+1)*phi) * mode(i_harm+1)
      endif

    enddo
    call interp_RZ(node_list, element_list, ielm, s, t, R, R_s, R_t, Z, Z_s, Z_t)
    
    xjac = R_s*Z_t - R_t*Z_s
    
    AR_R = (   Z_t * AR_s - Z_s * AR_t ) / xjac
    AR_Z = ( - R_t * AR_s + R_s * AR_t ) / xjac
    AZ_R = (   Z_t * AZ_s - Z_s * AZ_t ) / xjac
    AZ_Z = ( - R_t * AZ_s + R_s * AZ_t ) / xjac
    A3_R = (   Z_t * A3_s - Z_s * A3_t ) / xjac
    A3_Z = ( - R_t * A3_s + R_s * A3_t ) / xjac

    AR_R_ref = (   Z_t * AR_s_ref - Z_s * AR_t_ref ) / xjac
    AR_Z_ref = ( - R_t * AR_s_ref + R_s * AR_t_ref ) / xjac
    AZ_R_ref = (   Z_t * AZ_s_ref - Z_s * AZ_t_ref ) / xjac
    AZ_Z_ref = ( - R_t * AZ_s_ref + R_s * AZ_t_ref ) / xjac
    A3_R_ref = (   Z_t * A3_s_ref - Z_s * A3_t_ref ) / xjac
    A3_Z_ref = ( - R_t * A3_s_ref + R_s * A3_t_ref ) / xjac
   
    !call interp(node_list, element_list, ielm, 1, 1, s, t, psieq, dum01, dum02, dum03, dum04, dum05) 
    !call F_profile(xpoint,xcase,Z,ES%Z_xpoint,psieq,ES%psi_axis,ES%psi_bnd,Fprof,dum01,dum02,dum03,dum04,dum05,&
    !               dum06,dum07,dum08,dum09,dum10,dum11)
   
    call interp(node_list, element_list, ielm, 710, 1, s, t, Fprof, dum01, dum02, dum03, dum04, dum05)
  
    BR = ( A3_Z - AZ_p )/ R
    BZ = ( AR_p - A3_R )/ R
    Bp = ( AZ_R - AR_Z ) + Fprof / R
    BB = sqrt( Bp**2 + BR**2 + BZ**2 )

    BR_ref = ( A3_Z_ref - AZ_p_ref )/ R
    BZ_ref = ( AR_p_ref - A3_R_ref )/ R
    Bp_ref = ( AZ_R_ref - AR_Z_ref ) + Fprof / R
    BB_ref = sqrt( Bp_ref**2 + BR_ref**2 + BZ_ref**2 )

    Bp_initial = Fprof / R


    if(Bp * Fprof <= 0) then
      write(*,*) "Warning, Bp sign inversed. Bp ~ Fprof/R should be satisfied"
      write(*,*) "Bp = ", Bp
      write(*,*) "Fprof/R = ", Fprof/R
      write(*,*) "AZ_R = ", AZ_R
      write(*,*) "AR_Z = ", AR_Z
    endif
    
  end subroutine determine_field_3D


  subroutine find_axis_3D(node_list, element_list, Rguess, Zguess, distMissMatch)
    type (type_node_list),    intent(in)    :: node_list
    type (type_element_list), intent(in)    :: element_list
    real*8, intent(inout) :: Rguess, Zguess
    real*8, intent(out)   :: distMissMatch
    integer :: istep, numstep

    numstep = 2 * PI / stepsizePhi

    call find_RZ(node_list,element_list,Rguess,Zguess,RR,ZZ,ielm,s,t,ifail)
    phi = 0.d0

    Rguess = RR
    Zguess = ZZ

    !--- Start field line tracing from phi = 0 plane to finish a toroidal circle
    do istep = 1, numstep
      call do_step_3D(node_list, element_list)
    enddo

    distMissMatch = dsqrt((RR-Rguess)**2 + (ZZ-Zguess)**2)

  end subroutine find_axis_3D

  subroutine record_axis_3D(node_list, element_list, Raxis3Dstart, Zaxis3Dstart, numStepMax, Raxis3DRec, Zaxis3DRec)
    type (type_node_list),    intent(in)    :: node_list
    type (type_element_list), intent(in)    :: element_list
    real*8, intent(in) :: Raxis3Dstart, Zaxis3Dstart
    integer, intent(in) :: numStepMax
    real*8, intent(out) :: Raxis3DRec(:), Zaxis3DRec(:)
    integer :: istep
    real*8 :: distMissMatch

    call find_RZ(node_list,element_list,Raxis3Dstart,Zaxis3Dstart,RR,ZZ,ielm,s,t,ifail)
    phi = 0.d0

    Raxis3DRec(1) = RR
    Zaxis3DRec(1) = ZZ

    !--- Start field line tracing from phi = 0 plane to finish a toroidal circle
    do istep = 2, numStepMax
      call do_step_3D(node_list, element_list)
      Raxis3DRec(istep) = RR
      Zaxis3DRec(istep) = ZZ
    enddo

    distMissMatch = dsqrt((RR-Raxis3Dstart)**2 + (ZZ-Zaxis3Dstart)**2)
    write(*,*) "Finished recording axis at ..."
    write(*,*) "... Rend = ", RR, ", Zend = ", ZZ, ", PhiEnd/2pi = ", phi/2/PI
    write(*,*) "... distMissMatch = ", distMissMatch
    call save_vector_as_ascii(Raxis3DRec, numStepMax, "Raxis3DRec.txt")
    call save_vector_as_ascii(Zaxis3DRec, numStepMax, "Zaxis3DRec.txt")

  end subroutine record_axis_3D

  subroutine linspace(start, end, n, array)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: start, end
    real(8), intent(out) :: array(n)
    integer :: il
    real(8) :: step

    if (n <= 1) then
        print *, "Error: Number of points must be greater than 1."
        write(*,*) "input n = ", n
        stop
    end if

    step = (end - start) / real(n - 1)
    do il = 1, n
        array(il) = start + step * real(il - 1)
    end do
  end subroutine linspace


  subroutine meshgrid(x, y, nx, ny, XX, YY)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: x(nx), y(ny)
    real(8), intent(out) :: XX(nx, ny), YY(nx, ny)
    integer :: im, jm

    do jm = 1, ny
        do im = 1, nx
            XX(im, jm) = x(im)
            YY(im, jm) = y(jm)
        end do
    end do
  end subroutine meshgrid

  subroutine save_matrix_as_ascii(matrix, nx, ny, filename)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: matrix(nx, ny)
    character(len=*), intent(in) :: filename

    integer :: is, js
    integer :: unit

    ! Open the file for writing
    open(newunit=unit, file=filename, status='replace', action='write', form='formatted')

    ! Write the matrix to the file
    do js = 1, ny
        do is = 1, nx
            write(unit, '(F15.8)', advance='no') matrix(is, js)
            if (is < nx) then
                write(unit, '(" ")', advance='no')
            else
                write(unit, '(" ")')
            end if
        end do
    end do

    ! Close the file
    close(unit)
  end subroutine save_matrix_as_ascii

  subroutine save_vector_as_ascii(vector, nv, filename)
    implicit none
    integer, intent(in) :: nv
    real(8), intent(in) :: vector(nv)
    character(len=*), intent(in) :: filename

    integer :: iv
    integer :: unit

    ! Open the file for writing
    open(newunit=unit, file=filename, status='replace', action='write', form='formatted')

    ! Write the vector to the file
    do iv = 1, nv
        write(unit, '(F15.8)') vector(iv)
    end do

    ! Close the file
    close(unit)
  end subroutine save_vector_as_ascii


  subroutine find_min_index(matrix, nx, ny, min_i, min_j)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: matrix(nx, ny)
    integer, intent(out) :: min_i, min_j

    integer :: imin, jmin
    real(8) :: min_val

    ! Initialize min_val with the first element of the matrix
    min_val = matrix(1, 1)
    min_i = 1
    min_j = 1

    ! Loop through the matrix to find the minimum value and its index
    do jmin = 1, ny
        do imin = 1, nx
            if (matrix(imin, jmin) < min_val) then
                min_val = matrix(imin, jmin)
                min_i = imin
                min_j = jmin
            end if
        end do
    end do
  end subroutine find_min_index

end module mod_q_profile_fmhd_3D