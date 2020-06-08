!> Create a Poincare plot and q-profile for a fMHD JOREK restart file
program q_poincare_fmhd_slow

use constants
use mod_parameters
use data_structure
use equil_info
use phys_module        !  , only: xpoint, xcase, rst_format
use mod_import_restart
use mod_log_params
use mod_interp
use elements_nodes_neighbours

implicit none

integer, parameter :: npolturns = 1
real*8, parameter  :: stepsize = 1.d-4
integer, parameter :: npoints = 20

real*8 :: Rstart(npoints)
real*8 :: Zstart(npoints)
integer :: i, j, k, ielm, ifail, my_id, ierr
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

write(*,*) '***************************************'
write(*,*) '* q_poincare_fmhd_sloooow             *'
write(*,*) '***************************************'

my_id=0

call det_modes()
call initialise_parameters(my_id,  "__NO_FILENAME__")
call log_parameters(my_id)

call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)

write(*,*) '*** start tracing ***'

do i = 1, npoints
  Rstart(i) = ES%R_axis + REAL(i) * 0.8d0/npoints !### TODO: replace by actual width on midplane
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
!write(89,*) (BR / BB), (BZ / BB), (Bp / (RR * BB)) !###
    j = j + 1
   
    if ( phi > 2.d0*PI ) then
      phi   = phi - 2.d0*PI
      torturns = torturns + 1.d0
!write(88,*) RR, ZZ !###
    end if

    if ( (RR - ES%R_axis) > 0.d0 .and. (( ZZ - ES%Z_axis ) * ( Zold - ES%Z_axis )) < 0.d0 .and. j > 1 ) then
      polturns = polturns + 1.d0
    else if ( (RR - ES%R_axis) > 0.d0 .and. (( ZZ - ES%Z_axis ) * ( Zold - ES%Z_axis )) == 0.d0 .and. j > 1 ) then
      polturns = polturns + 0.5d0
    end if
     
    if ( torturns > 399.d0 ) then
      stop_tracing = .true.
    end if
    
  !  if ( j > 2 ) then
  !    stop_tracing = .true.
  !  end if

  end do
write(90,*) RR, get_psi_n(A3, ZZ), torturns, polturns
!write(88,*) !###
!write(88,*) !###

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

end program q_poincare_fmhd_slow
