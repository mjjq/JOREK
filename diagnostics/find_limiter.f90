!> Finds the limiter position.
!!
!! * The routine returns the R, Z, and Psi-value of the limiter point with the smallest Psi-value.
!! * The points checked are all points on the boundary of the JOREK domain and additional limiter
!!   points given in the namelist input file (parameters n_limiter, R_limiter, Z_limiter).
!!
subroutine find_limiter(my_id, node_list, element_list, bnd_elm_list, psi_lim, R_lim, Z_lim)

use phys_module, only: n_limiter, R_limiter, Z_limiter, FF_0
use data_structure
use gauss
use basis_at_gaussian

implicit none

! --- Routine parameters
integer,                      intent(in)  :: my_id
type (type_node_list),        intent(in)  :: node_list
type (type_element_list),     intent(in)  :: element_list
type (type_bnd_element_list), intent(in)  :: bnd_elm_list
real*8,                       intent(out) :: psi_lim
real*8,                       intent(out) :: R_lim
real*8,                       intent(out) :: Z_lim

! --- Local variables
real*8  :: s_lim,t_lim, r_min, r_max
real*8  :: r, psim, psimr, psip, psipr, psma, psmi, psmima, psi_min, psi_max, AA, BB, CC, DD, DET, dummy
real*8  :: RM, RMR, RP, RPR, ZM, ZMR, ZP, ZPR
integer :: i_limiter, i_bnd_lim, ibnd, n1, n2, idir1, idir2, i_min, i_max
real*8  :: psi, psi_s,psi_t,psi_st,psi_ss,psi_tt, s_out, t_out, R_out, Z_out
integer :: ifail, i_elm

real*8, external :: root

if ( my_id == 0 ) then
  write(*,*) '*********************************'
  write(*,*) '*     find_limiter              *'
  write(*,*) '*********************************'
end if

psi_lim =  0.d0
psi_min =  1.d20
psi_max = -1.d20

R_lim = 0.d0
Z_lim = 0.d0
s_lim = 0.d0
t_lim = 0.d0

r_min = 0.d0
i_min = 0
r_max = 0.d0
i_max = 0

do ibnd=1,bnd_elm_list%n_bnd_elements
        
  n1 = bnd_elm_list%bnd_element(ibnd)%vertex(1)
  n2 = bnd_elm_list%bnd_element(ibnd)%vertex(2)
    
  idir1 = bnd_elm_list%bnd_element(ibnd)%direction(1,2) 
  idir2 = bnd_elm_list%bnd_element(ibnd)%direction(2,2) 

  PSIM  =  node_list%node(n1)%values(1,1,1)     * bnd_elm_list%bnd_element(ibnd)%size(1,1)              ! PSI(1,n1)
  PSIMR =  node_list%node(n1)%values(1,idir1,1) * bnd_elm_list%bnd_element(ibnd)%size(1,2) * 3.d0/2.d0  ! PSI(3,n1)
  PSIP  =  node_list%node(n2)%values(1,1,1)     * bnd_elm_list%bnd_element(ibnd)%size(2,1)              ! PSI(1,n2)
  PSIPR =  - node_list%node(n2)%values(1,idir2,1) * bnd_elm_list%bnd_element(ibnd)%size(2,2) * 3.d0/2.d0  ! PSI(3,n2)

  PSMA = MAX(PSIM,PSIP)
  PSMI = MIN(PSIM,PSIP)
  AA =  3.d0 * (PSIM + PSIMR - PSIP + PSIPR ) / 4.d0
  BB =  ( - PSIMR + PSIPR ) / 2.d0
  CC =  ( - 3.d0*PSIM - PSIMR + 3.d0*PSIP - PSIPR) / 4.d0
  DET = BB**2 - 4.d0*AA*CC
  
  if (DET .GE. 0.d0) then
    R = ROOT(AA,BB,CC,DET,1.d0)
    if (ABS(R) .GT. 1.d0) then
      R = ROOT(AA,BB,CC,DET,-1.d0)
    endif
    if (ABS(R) .LE. 1.d0) then
      call CUB1D(PSIM,PSIMR,PSIP,PSIPR,R,PSMIMA,DUMMY)
      psma = max(psma,psmima)
      psmi = min(psmi,psmima)
      
      if (psmi .lt. psi_min) then
        psi_min = psmi
        r_min   = r
        i_min   = ibnd
      endif
      
      if (psma .gt. psi_max) then
        psi_max = psma
        r_max   = r
        i_max   = ibnd
      endif
      
    endif
  endif

enddo

if (FF_0 .gt. 0.d0) then
  if ((i_min .gt. 0) .and. (r_min .le. 1.d0)) then

    n1 = bnd_elm_list%bnd_element(i_min)%vertex(1)
    n2 = bnd_elm_list%bnd_element(i_min)%vertex(2)
    
    idir1 = bnd_elm_list%bnd_element(i_min)%direction(1,2) 
    idir2 = bnd_elm_list%bnd_element(i_min)%direction(2,2) 

    RM  =  node_list%node(n1)%x(1,1)	 * bnd_elm_list%bnd_element(i_min)%size(1,1)		  
    RMR =  node_list%node(n1)%x(idir1,1) * bnd_elm_list%bnd_element(i_min)%size(1,2) * 3.d0/2.d0  
    RP  =  node_list%node(n2)%x(1,1)	 * bnd_elm_list%bnd_element(i_min)%size(2,1)		  
    RPR =  node_list%node(n2)%x(idir2,1) * bnd_elm_list%bnd_element(i_min)%size(2,2) * 3.d0/2.d0 

    call CUB1D(RM,RMR,RP,RPR,r_min,R_lim,DUMMY)

    ZM  =  node_list%node(n1)%x(1,2)	 * bnd_elm_list%bnd_element(i_min)%size(1,1)		  
    ZMR =  node_list%node(n1)%x(idir1,2) * bnd_elm_list%bnd_element(i_min)%size(1,2) * 3.d0/2.d0  
    ZP  =  node_list%node(n2)%x(1,2)	 * bnd_elm_list%bnd_element(i_min)%size(2,1)		  
    ZPR =  node_list%node(n2)%x(idir2,2) * bnd_elm_list%bnd_element(i_min)%size(2,2) * 3.d0/2.d0 
    
    call CUB1D(ZM,ZMR,ZP,ZPR,r_min,Z_lim,DUMMY)
    
    psi_lim = psi_min

  else
    psi_lim = 999.d0
    R_lim   = 0.d0
    Z_lim   = 0.d0
  endif
else
  if ((i_max .gt. 0) .and. (r_max .le. 1.d0)) then

    n1 = bnd_elm_list%bnd_element(i_max)%vertex(1)
    n2 = bnd_elm_list%bnd_element(i_max)%vertex(2)
    
    idir1 = bnd_elm_list%bnd_element(i_max)%direction(1,2) 
    idir2 = bnd_elm_list%bnd_element(i_max)%direction(2,2) 

    RM  =  node_list%node(n1)%x(1,1)	 * bnd_elm_list%bnd_element(i_max)%size(1,1)		  
    RMR =  node_list%node(n1)%x(idir1,1) * bnd_elm_list%bnd_element(i_max)%size(1,2) * 3.d0/2.d0  
    RP  =  node_list%node(n2)%x(1,1)	 * bnd_elm_list%bnd_element(i_max)%size(2,1)		  
    RPR =  node_list%node(n2)%x(idir2,1) * bnd_elm_list%bnd_element(i_max)%size(2,2) * 3.d0/2.d0 

    call CUB1D(RM,RMR,RP,RPR,r_max,R_lim,DUMMY)

    ZM  =  node_list%node(n1)%x(1,2)	 * bnd_elm_list%bnd_element(i_max)%size(1,1)		  
    ZMR =  node_list%node(n1)%x(idir1,2) * bnd_elm_list%bnd_element(i_max)%size(1,2) * 3.d0/2.d0  
    ZP  =  node_list%node(n2)%x(1,2)	 * bnd_elm_list%bnd_element(i_max)%size(2,1)		  
    ZPR =  node_list%node(n2)%x(idir2,2) * bnd_elm_list%bnd_element(i_max)%size(2,2) * 3.d0/2.d0 
    
    call CUB1D(ZM,ZMR,ZP,ZPR,r_max,Z_lim,DUMMY)
    
    psi_lim = psi_max

  else
    psi_lim = 999.d0
    R_lim   = 0.d0
    Z_lim   = 0.d0
  endif
endif

! --- Take into account additional limiter points from the namelist input file
do i_limiter = 1, n_limiter
  Rp = R_limiter(i_limiter)
  Zp = Z_limiter(i_limiter)
  
  call find_RZ(node_list, element_list, Rp, Zp, R_out, Z_out, i_elm, s_out, t_out, ifail)
  call interp(node_list, element_list, i_elm, 1, 1, s_out, t_out, psi, psi_s, psi_t, psi_st,       &
    psi_ss, psi_tt)
  
  if (FF_0 .gt. 0.d0) then
    if (psi .lt. psi_lim) then
      psi_lim = psi
      R_lim   = Rp
      Z_lim   = Zp
    end if
  else
    if (psi .gt. psi_lim) then
      psi_lim = psi
      R_lim   = Rp
      Z_lim   = Zp
    end if
  endif
end do

if ( my_id == 0 ) then
  121 format(1x,a,' =',f12.7)
  write(*,121) 'R_lim  ', R_lim
  write(*,121) 'Z_lim  ', Z_lim
  write(*,121) 'Psi_lim', Psi_lim
end if

return
end subroutine find_limiter
