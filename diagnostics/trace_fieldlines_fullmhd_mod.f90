module trace_fieldlines_fullmhd_mod

implicit none

contains

subroutine trace_fieldlines_fullmhd(node_list, element_list, nlines, Rstart, Zstart, phistart, nturns, direction)

use constants
use mod_parameters
use data_structure
use equil_info
use phys_module, only: xpoint, xcase

implicit none

real*8, parameter :: stepsize = 1.d-4 !###

! --- Routine parameters
type(type_node_list),    intent(in) :: node_list
type(type_element_list), intent(in) :: element_list
integer, intent(in) :: nlines
real*8,  intent(in) :: Rstart(nlines)
real*8,  intent(in) :: Zstart(nlines)
real*8,  intent(in) :: phistart(nlines)
integer, intent(in) :: nturns(nlines)
real*8,  intent(in) :: direction

! --- Local variables
integer :: i, j, k, ielm, ifail, turns
real*8  :: Rout, Zout, length
logical :: stop_tracing

real*8  :: s, t, phi, phinew, phiold
real*8  :: AA(3), AA_s(3), AA_t(3)
real*8  :: R, R_s, R_t, Z, Z_s, Z_t
real*8  :: AA_p(3)
real*8  :: BR, BZ, Bp, BB, Bs, Bt, xjac, Fprof
real*8  :: dum01, dum02, dum03, dum04, dum05, dum06, dum07, dum08, dum09, dum10, dum11
real*8  :: AR, AR_p, AR_s, AR_t, AR_R, AR_Z, AZ, AZ_p, AZ_s, AZ_t, AZ_R, AZ_Z, A3, A3_p, A3_s, A3_t, A3_R, A3_Z, psieq
real*8  :: RR, ZZ, Rnew, Znew, Rold, Zold

if ( (direction /= 1.d0) .and. (direction /= -1.d0) ) then
  write(*,*) 'Illegal value for direction'
  stop
end if

do i = 1, nlines
  call find_RZ(node_list,element_list,Rstart(i),Zstart(i),RR,ZZ,ielm,s,t,ifail)
  phi = phistart(i)
  
  if ( ielm < 1 ) then
    write(*,*) 'Illegal starting point for field line ', i, '. Skipping.'
    cycle
  end if
  
  stop_tracing = .false.
  length       = 0.d0
  j            = 0
  turns        = 0
  do while( .not. stop_tracing )
    
    call do_step()
    
    if ( direction*phi > 2.d0*PI ) then
      phi   = phi - direction*2.d0*PI
      turns = turns + 1
      write(88,*) RR, ZZ !###
    end if
    
    if ( turns > nturns(i) ) then
      stop_tracing = .true.
    end if
    
    j = j + 1
  end do
  write(88,*) !###
  write(88,*) !###
  
end do




contains
  
  
  
  !> A "stupid" stepping routine advancing in R, Z, phi
  subroutine do_step()
  
  call determine_field()
  
  Rnew   = RR  + direction * 0.5d0 * stepsize * BR / BB
  Znew   = ZZ  + direction * 0.5d0 * stepsize * BZ / BB
  phinew = phi + direction * 0.5d0 * stepsize * Bp / BB
  
  Rold   = RR
  Zold   = ZZ
  phiold = phi
  
  call find_RZ(node_list,element_list,Rnew,Znew,RR,ZZ,ielm,s,t,ifail)
  if ( ielm < 1 ) then
    stop_tracing = .true.
    return
  end if
  call determine_field()
  
  RR   = Rold   + direction * stepsize * BR / BB
  ZZ   = Zold   + direction * stepsize * BZ / BB
  phi  = phiold + direction * stepsize * Bp / BB
  
  call find_RZ(node_list,element_list,Rnew,Znew,RR,ZZ,ielm,s,t,ifail)
  if ( ielm < 1 ) then
    stop_tracing = .true.
    return
  end if
  
  length = length + stepsize
  
  
  
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
  ! interp(node_list, element_list, i_elm, 456, 1, s, t, Fprof, dum01, dum02, dum03, dum04, dum05)

  BR = ( A3_Z - AZ_p )/ R
  BZ = ( AR_p - A3_R )/ R
  Bp = ( AZ_R - AR_Z ) + Fprof / R
  BB = sqrt( Bp**2 + BR**2 + BZ**2 )
  
  Bs = 0.d0 !###
  Bt = 0.d0 !###
  
  end subroutine determine_field


end subroutine trace_fieldlines_fullmhd


end module trace_fieldlines_fullmhd_mod
