subroutine RZ_minmax(node_list,element_list,i_elm,Rmin,Rmax,Zmin,Zmax)

use data_structure

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list

real*8  :: Rmin, Rmax, Zmin, Zmax, psimin, psimax, psma, psmi, psmima, psim, psimr, psip, psipr
real*8  :: aa, bb, cc, det, r, dummy
real*8,external :: root
integer :: i_elm, iv, n, im, n1, n2
real*8  :: s,t,P,P_s,P_t,P_st,P_ss,P_tt, ss1, ss2, ss3
integer :: k

do k=1,2

  psimin = 1d10
  psimax =-1d10

  do iv= 1, n_vertex_max

    im = mod(iv,n_vertex_max) + 1
    n1 = element_list%element(i_elm)%vertex(iv)
    n2 = element_list%element(i_elm)%vertex(im)
    
    ss1 = element_list%element(i_elm)%size(iv,1)

    if ((iv .eq. 1) .or. (iv .eq. 3)) THEN

      ss2 = element_list%element(i_elm)%size(iv,2)

      PSIM  =  node_list%node(n1)%x(1,k) * ss1             ! PSI(1,n1)
      PSIMR =  node_list%node(n1)%x(2,k) * ss2 * 3.d0/2.d0 ! PSI(3,n1)
      PSIP  =  node_list%node(n2)%x(1,k) * ss1             ! PSI(1,n2)
      PSIPR =  node_list%node(n2)%x(2,k) * ss2 * 3.d0/2.d0 ! PSI(3,n2)

    elseif ((iv .eq. 2) .or. (iv .eq. 4)) then

      ss3 = element_list%element(i_elm)%size(iv,3)
      
      PSIM  =   node_list%node(n1)%x(1,k) * ss1             ! PSI(1,n1)
      PSIMR =   node_list%node(n1)%x(3,k) * ss3 * 3.d0/2.d0 ! PSI(2,n1)
      PSIP  =   node_list%node(n2)%x(1,k) * ss1             ! PSI(1,n2)
      PSIPR =   node_list%node(n2)%x(3,k) * ss3 * 3.d0/2.d0 ! PSI(2,n2)

    endif

    if ((PSIM .eq. PSIP) .and. (PSIMR .eq. 0.d0) .and. (PSIPR .eq. 0.d0)) then

      psimin = min(psimin,PSIM)
      psimax = max(psimax,PSIM)

    else

      PSMA = MAX(PSIM,PSIP)
      PSMI = MIN(PSIM,PSIP)
      AA =  3.d0 * (PSIM + PSIMR - PSIP + PSIPR ) / 4.d0
      BB =  ( - PSIMR + PSIPR ) / 2.d0
      CC =  ( - 3.d0*PSIM - PSIMR + 3.d0*PSIP - PSIPR) / 4.d0
      DET = BB**2 - 4.d0*AA*CC
      IF (DET .GE. 0.d0) THEN
        R = ROOT(AA,BB,CC,DET,1.d0)
        IF (ABS(R) .GT. 1.d0) THEN
          R = ROOT(AA,BB,CC,DET,-1.d0)
        ENDIF
        IF (ABS(R) .LE. 1.d0) THEN
          CALL CUB1D(PSIM,PSIMR,PSIP,PSIPR,R,PSMIMA,DUMMY)
          psma = max(psma,psmima)
          psmi = min(psmi,psmima)
        ENDIF
      ENDIF
      psimin = min(psimin,psmi)
      psimax = max(psimax,psma)

    endif

  enddo

  if (k.eq.1) then
    Rmin = psimin
    Rmax = psimax
  else
    Zmin = psimin
    Zmax = psimax
  endif

ENDDO

RETURN
END
