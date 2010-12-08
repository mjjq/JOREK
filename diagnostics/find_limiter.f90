subroutine find_limiter(node_list,bnd_elm_list,psi_lim,R_lim,Z_lim)
!-----------------------------------------------------------------------
! now only looks for the minimum flux on the boundary. 
! An arbitrary list deifning the limiter surface should be added
!
!-----------------------------------------------------------------------
use data_structure
use gauss
use basis_at_gaussian

implicit none

type (type_node_list)        :: node_list
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

real*8  :: psi_lim,R_lim,Z_lim,s_lim,t_lim, r_min
real*8  :: r, psim, psimr, psip, psipr, psma, psmi, psmima, psi_min, psi_max, AA, BB, CC, DD, DET, dummy
real*8  :: RM, RMR, RP, RPR, ZM, ZMR, ZP, ZPR
integer :: i_bnd_lim, ibnd, n1, n2, idir1, idir2, i_min

real*8, external :: root

write(*,*) '*********************************'
write(*,*) '*     find_limiter              *'
write(*,*) '*********************************'

psi_lim =  0.d0
psi_min =  1.d20
psi_max = -1.d20

R_lim = 0.d0
Z_lim = 0.d0
s_lim = 0.d0
t_lim = 0.d0

r_min = 0.d0
i_min = 0

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
      
    endif
  endif

psi_max = max(psi_max,psma)
  
enddo

if ((i_min .gt. 0) .and. (r_min .le. 1.d0)) then

  n1 = bnd_elm_list%bnd_element(i_min)%vertex(1)
  n2 = bnd_elm_list%bnd_element(i_min)%vertex(2)
  
  idir1 = bnd_elm_list%bnd_element(i_min)%direction(1,2) 
  idir2 = bnd_elm_list%bnd_element(i_min)%direction(2,2) 

  RM  =  node_list%node(n1)%x(1,1)     * bnd_elm_list%bnd_element(i_min)%size(1,1)              
  RMR =  node_list%node(n1)%x(idir1,1) * bnd_elm_list%bnd_element(i_min)%size(1,2) * 3.d0/2.d0  
  RP  =  node_list%node(n2)%x(1,1)     * bnd_elm_list%bnd_element(i_min)%size(2,1)              
  RPR =  node_list%node(n2)%x(idir2,1) * bnd_elm_list%bnd_element(i_min)%size(2,2) * 3.d0/2.d0 

  call CUB1D(RM,RMR,RP,RPR,r_min,R_lim,DUMMY)

  ZM  =  node_list%node(n1)%x(1,2)     * bnd_elm_list%bnd_element(i_min)%size(1,1)              
  ZMR =  node_list%node(n1)%x(idir1,2) * bnd_elm_list%bnd_element(i_min)%size(1,2) * 3.d0/2.d0  
  ZP  =  node_list%node(n2)%x(1,2)     * bnd_elm_list%bnd_element(i_min)%size(2,1)              
  ZPR =  node_list%node(n2)%x(idir2,2) * bnd_elm_list%bnd_element(i_min)%size(2,2) * 3.d0/2.d0 
  
  call CUB1D(ZM,ZMR,ZP,ZPR,r_min,Z_lim,DUMMY)
  
  psi_lim = psi_min

else
  psi_lim = 999.d0
  R_lim   = 0.d0
  Z_lim   = 0.d0
endif

return
END
