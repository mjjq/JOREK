module qprofile_fmhd_mod

implicit none

contains

subroutine qprofile_fmhd(node_list, element_list, npoints)

use constants
use mod_parameters
use data_structure
use equil_info
use phys_module, only: xpoint, xcase

implicit none

integer, parameter :: npolturns = 10
real*8, parameter :: stepsize = 1.d-4

! --- Routine parameters
type(type_node_list),    intent(in) :: node_list
type(type_element_list), intent(in) :: element_list
integer, intent(in) :: npoints

! --- Local variables
real*8 :: Rstart(npoints)
real*8 :: Zstart(npoints)
integer :: i, j, k, ielm, ifail
real*8  :: Rout, Zout, polturns, torturns
logical :: stop_tracing

real*8  :: s, t, phi, phinew, phiold
real*8  :: AA(3), AA_s(3), AA_t(3)
real*8  :: R, R_s, R_t, Z, Z_s, Z_t
real*8  :: AA_p(3)
real*8  :: BR, BZ, Bp, BB, Bs, Bt, xjac, Fprof
real*8  :: dum01, dum02, dum03, dum04, dum05, dum06, dum07, dum08, dum09, dum10, dum11
real*8  :: AR, AR_p, AR_s, AR_t, AR_R, AR_Z, AZ, AZ_p, AZ_s, AZ_t, AZ_R, AZ_Z, A3, A3_p, A3_s, A3_t, A3_R, A3_Z, psieq
real*8  :: RR, ZZ, Rnew, Znew, Rold, Zold

do i = 1, npoints
  Rstart(i) = ES%R_axis + REAL(i) * 0.85d0/npoints !### TODO: replace by actual width on midplane
  Zstart(i) = ES%Z_axis
end do

do i = 1, npoints
  call find_RZ(node_list,element_list,Rstart(i),Zstart(i),RR,ZZ,ielm,s,t,ifail)
  phi = 0.d0
  
  if ( ielm < 1 ) then
    write(*,*) 'Illegal starting point for field line ', i, '. Skipping.'
    cycle
  end if
  
  stop_tracing = .false.
  polturns     = 0.d0
  torturns     = 0.d0
  j            = 0
  do while( .not. stop_tracing )
    
    call do_step()
    j = j + 1

! write(89,*) RR, ZZ !###
   
    if ( phi > 2.d0*PI ) then
      phi   = phi - 2.d0*PI
      torturns = torturns + 1.d0
write(88,*) RR, ZZ !###
    end if

    if ( (RR - ES%R_axis) > 0.d0 .and. (( ZZ - ES%Z_axis ) * ( Zold - ES%Z_axis )) < 0.d0 .and. j > 1 ) then
      polturns = polturns + 1.d0
    else if ( (RR - ES%R_axis) > 0.d0 .and. (( ZZ - ES%Z_axis ) * ( Zold - ES%Z_axis )) == 0.d0 .and. j > 1 ) then
      polturns = polturns + 0.5d0
    end if
      
!    if ( polturns > REAL(npolturns) ) then
!      stop_tracing = .true.
!    end if
    
!    if ( torturns > (REAL(npolturns) * 5.) ) then
   if ( torturns > 99.d0 ) then
      stop_tracing = .true.
    end if

  end do
write(90,*) Rstart(i), torturns, polturns
write(88,*) !###
write(88,*) !###
!write(89,*) !###
!write(89,*) !###

end do 

contains


  !> A "stupid" stepping routine advancing in R, Z, phi
  subroutine do_step()
  
  call determine_field()
  
  Rnew   = RR  + 0.5d0 * stepsize * BR / BB
  Znew   = ZZ  + 0.5d0 * stepsize * BZ / BB
  phinew = phi + 0.5d0 * stepsize * Bp / (RR * BB)
  
  Rold   = RR
  Zold   = ZZ
  phiold = phi
  
  call find_RZ(node_list,element_list,Rnew,Znew,RR,ZZ,ielm,s,t,ifail)
  if ( ielm < 1 ) then
    stop_tracing = .true.
    return
  end if
  call determine_field()
  
  RR   = Rold   + stepsize * BR / BB
  ZZ   = Zold   + stepsize * BZ / BB
  phi  = phiold + stepsize * Bp / (RR * BB)
  
  call find_RZ(node_list,element_list,Rnew,Znew,RR,ZZ,ielm,s,t,ifail)
  if ( ielm < 1 ) then
    stop_tracing = .true.
    return
  end if

  end subroutine do_step
 

 !> Determine magnetic field at given position
  subroutine determine_field()
  
  call interp_PRZ(node_list, element_list, ielm, (/1,2,3/), 3, s, t, phi, AA, AA_s, AA_t, AA_p, R, R_s, R_t, Z, Z_s, Z_t)
  
  xjac = R_s*Z_t - R_t*Z_s
  
  AR   = AA(2)
  AR_p = AA_p(2)
  AR_s = AA_s(2)
  AR_t = AA_t(2)
  AR_R = (   Z_t * AA_s(2) - Z_s * AA_t(2) ) / xjac
  AR_Z = ( - R_t * AA_s(2) + R_s * AA_t(2) ) / xjac
  AZ   = AA(3)
  AZ_p = AA_p(3)
  AZ_s = AA_s(3)
  AZ_t = AA_t(3)
  AZ_R = (   Z_t * AA_s(3) - Z_s * AA_t(3) ) / xjac
  AZ_Z = ( - R_t * AA_s(3) + R_s * AA_t(3) ) / xjac
  A3   = AA(1)
  A3_p = AA_p(1)
  A3_s = AA_s(1)
  A3_t = AA_t(1)
  A3_R = (   Z_t * AA_s(1) - Z_s * AA_t(1) ) / xjac
  A3_Z = ( - R_t * AA_s(1) + R_s * AA_t(1) ) / xjac

  ! Fprof = 0.d0  
  call interp(node_list, element_list, ielm, 1, 1, s, t, psieq, dum01, dum02, dum03, dum04, dum05) 

  call F_profile(xpoint,xcase,Z,ES%Z_xpoint,psieq,ES%psi_axis,ES%psi_bnd,Fprof,dum01,dum02,dum03,dum04,dum05,&
                 dum06,dum07,dum08,dum09,dum10,dum11)

  !### better to use interp & ivar=456 instead ? 
  !call interp(node_list, element_list, ielm, 456, 1, s, t, Fprof, dum01, dum02, dum03, dum04, dum05)

  BR = ( A3_Z - AZ_p )/ R
  BZ = ( AR_p - A3_R )/ R
  Bp = ( AZ_R - AR_Z ) + Fprof / R
  BB = sqrt( Bp**2 + BR**2 + BZ**2 )
  
  Bs = 0.d0 !###
  Bt = 0.d0 !###
  
  end subroutine determine_field 
   
!  !> Determine magnetic field at given position (n=0 only)
!  subroutine determine_field()
!  
!! call interp_PRZ(node_list, element_list, ielm, (/1,2,3/), 3, s, t, phi, AA, AA_s, AA_t, AA_p, R, R_s, R_t, Z, Z_s, Z_t)
!  call interp(node_list, element_list, ielm, 1, 1, s, t, A3, A3_s, A3_t, dum01, dum02, dum03)
!  call interp(node_list, element_list, ielm, 2, 1, s, t, AR, AR_s, AR_t, dum01, dum02, dum03)
!  call interp(node_list, element_list, ielm, 3, 1, s, t, AZ, AZ_s, AZ_t, dum01, dum02, dum03)
!  call interp_RZ(node_list, element_list, ielm, s, t, R, R_s, R_t, Z, Z_s, Z_t)
!
!  xjac = R_s*Z_t - R_t*Z_s
!  
!  AR_R = (   Z_t * AR_s - Z_s * AR_t ) / xjac
!  AR_Z = ( - R_t * AR_s + R_s * AR_t ) / xjac
!  AZ_R = (   Z_t * AZ_s - Z_s * AZ_t ) / xjac
!  AZ_Z = ( - R_t * AZ_s + R_s * AZ_t ) / xjac
!  A3_R = (   Z_t * A3_s - Z_s * A3_t ) / xjac
!  A3_Z = ( - R_t * A3_s + R_s * A3_t ) / xjac
!  
!  AR_p = 0.d0
!  AZ_p = 0.d0
!
!  ! Fprof = 0.d0  
!  call interp(node_list, element_list, ielm, 1, 1, s, t, psieq, dum01, dum02, dum03, dum04, dum05) 
!
!  call F_profile(xpoint,xcase,Z,ES%Z_xpoint,psieq,ES%psi_axis,ES%psi_bnd,Fprof,dum01,dum02,dum03,dum04,dum05,&
!                 dum06,dum07,dum08,dum09,dum10,dum11)
!
!  !### better to use interp & ivar=456 instead ? 
!  ! interp(node_list, element_list, i_elm, 456, 1, s, t, Fprof, dum01, dum02, dum03, dum04, dum05)
!
!  BR = ( A3_Z - AZ_p )/ R
!  BZ = ( AR_p - A3_R )/ R
!  Bp = ( AZ_R - AR_Z ) + Fprof / R
!  BB = sqrt( Bp**2 + BR**2 + BZ**2 )
!  
!  Bs = 0.d0 !###
!  Bt = 0.d0 !###
!  
!  end subroutine determine_field

end subroutine qprofile_fmhd

end module qprofile_fmhd_mod
