subroutine grid_polar_bezier(Rgeo,Zgeo,amin,acentre,fbnd,fpsi,mf,nr,np,node_list,element_list)
!***********************************************************************
! defines a polar grid using Bezier finite elements (using the HELENA
! cubic Hermite elements formulation)
!
! input :
!          Rgeo,Zgeo : the position of the geometric center
!          amin      : the minor radius
!          acentre   : smallest radius
!          fbnd(1:mf): Fourier series describing the radius as function
!                      of the poloidal angle
!          fpsi(1:mf): Fourier series of flux at the boundary
!          nr        : number of radial points, (nr-1) elements
!          np        : number of poloidal points, np elements
!
! output :
!
!          node_list    : list of nodes with grid information
!          element_list : list of elements with element information
!***********************************************************************
use parameters
use data_structure

implicit none

real*8              :: Rgeo, Zgeo,amin,acentre,fbnd(*),fpsi(*),angle_start
real*8, allocatable :: RR(:,:),ZZ(:,:),PSI(:,:)
real*8              :: r_small, pi, dt, ds, thtj, radius, rm, drm, drmt, drmtr, angle, psi_axis
real*8              :: delta_rm, delta_zm, delta_rp, delta_zp, dir_2, dir_3
integer             :: mf, nr, np, i, j, m, index, index0, node, k, iv, ivp, ivm, node_iv, node_ivp, node_ivm,i_sons
integer             :: n_element_start, n_node_start, n_boundary_start, n_index_start, ielm, iside,iv1, iv2, idir1, idir2
real*8              :: XR_r_0,XR_r_1, SIG_r_0, SIG_r_1, XR_tht_0, XR_tht_1, SIG_tht_0, SIG_tht_1, abltg(3), s_tmp, dr_ds, dtht_dt
real*8, allocatable :: S1(:), S2(:), SP1(:), SP2(:), SP3(:), SP4(:)
real*8, allocatable :: T1(:), T2(:), TP1(:), TP2(:), TP3(:), TP4(:)
real*8, external    :: spwert

type(type_node_list)     :: node_list
type(type_element_list)  :: element_list

allocate(RR(4,nr*np),ZZ(4,nr*np),PSI(4,nr*np))

pi = 2.d0*asin(1.d0)

dt = 2.d0*pi/real(np)
ds = 1.d0/real(nr-1)

angle_start = - 3.d0 * PI /4.d0

n_element_start  = element_list%n_elements
n_node_start     = node_list%n_nodes

n_index_start = 0
do i=1,n_node_start
  n_index_start = max(n_index_start,maxval(node_list%node(i)%index(:)))
enddo

write(*,*) '*************************************'
write(*,*) '*        grid_polar_bezier          *'
write(*,*) '*************************************'
write(*,*) ' existing number of elements          : ',n_element_start
write(*,*) ' existing number of nodes             : ',n_node_start
write(*,*) ' index_start                          : ',n_index_start

psi_axis = -0.1d0

XR_r_0  = 999.d0 !0.80d0         ! mesh accumulation parameters radial position
SIG_r_0 = 999.d0 !0.1d0          ! width of accumulation (Gaussian)
XR_r_1  = 999.d0
SIG_r_1 = 999.d0

XR_tht_0  = 999.d0 ! 0.d0 !0.5d0 !0. ! mesh accumulation parameters poloidal position
SIG_tht_0 = 999.d0 ! 0.03d0          ! width of accumulation (Gaussian)
XR_tht_1  = 999.d0 ! 1.d0 !1.d0
SIG_tht_1 = 999.d0 ! 0.03d0

allocate(S1(nr),S2(nr),SP1(nr),SP2(nr),SP3(nr),SP4(nr))
S2 = 0
allocate(T1(np+1),T2(np+1),TP1(np+1),TP2(np+1),TP3(np+1),TP4(np+1))
T2 = 0
do i=1,nr
  S1(i) = real(i-1)/real(nr-1)
enddo
do j=1,np+1
  T1(j) = real(j-1)/real(np)
enddo

call meshac2(nr,S2,XR_r_0,XR_r_1,SIG_r_0,SIG_r_1,0.6d0,1.0d0)
call spline(nr,S1,S2,0.d0,0.d0,2,SP1,SP2,SP3,SP4)

call meshac2(np+1,T2,XR_tht_0,XR_tht_1,SIG_tht_0,SIG_tht_1,0.6d0,1.0d0)
call spline(np+1,T1,T2,0.d0,0.d0,2,TP1,TP2,TP3,TP4)


do i=1,nr

 radius = spwert(nr,S1(i),SP1,SP2,SP3,SP4,S1,ABLTG)
 dr_ds  = abltg(1)
 
! write(*,'(4e12.4)') S1(i),S2(i),radius,dr_ds
 
 do  j=1,np
  
   node   = np*(i-1) + j

   thtj     = spwert(np+1,T1(j),TP1,TP2,TP3,TP4,T1,ABLTG) * 2.d0 * PI
   dtht_dt  = abltg(1)

   if (i.eq.nr) write(*,'(4e12.4)') T1(j),T2(j),thtj,dtht_dt

   RR(1,node) = Rgeo + amin * radius * fbnd(1) * cos(thtj) / 2.d0
   RR(2,node) =        amin *          fbnd(1) * cos(thtj) / 2.d0
   RR(3,node) =      - amin * radius * fbnd(1) * sin(thtj) / 2.d0
   RR(4,node) =      - amin          * fbnd(1) * sin(thtj) / 2.d0
   ZZ(1,node) = Zgeo + amin * radius * fbnd(1) * sin(thtj) / 2.d0
   ZZ(2,node) =        amin *          fbnd(1) * sin(thtj) / 2.d0
   ZZ(3,node) =        amin * radius * fbnd(1) * cos(thtj) / 2.d0
   ZZ(4,node) =        amin *          fbnd(1) * cos(thtj) / 2.d0

   PSI(1,node) =        radius**8 * fpsi(1) / 2.d0 + psi_axis*(1.d0 - radius**2)
   PSI(2,node) =  8.d0* radius**7 * fpsi(1) / 2.d0 - psi_axis*2.d0* radius
   PSI(3,node) =  0.d0
   PSI(4,node) =  0.d0

!   write(*,'(A,2i6,12e16.8)') ' PSI(1,node) : ', node,1,PSI(1,node),radius,thtj,fpsi(1),fpsi(2)

   do m = 2, mf/2

     if (m .eq. 2) then
       rm   = radius * ( fbnd(2*M-1) * cos((M-1)*THTJ)           + fbnd(2*M) * sin((M-1)*THTJ) )
       drm  =          ( fbnd(2*M-1) * cos((M-1)*THTJ)           + fbnd(2*M) * sin((M-1)*THTJ))
       drmt = radius * (-fbnd(2*M-1) * (M-1)*sin((M-1)*THTJ)     + fbnd(2*M) * (M-1)*cos((M-1)*THTJ))
       drmtr=          (-fbnd(2*M-1) * (M-1)*sin((M-1)*THTJ)     + fbnd(2*M) * (M-1)*cos((M-1)*THTJ))
     else
       rm   =      radius**(M-1) * ( fbnd(2*M-1) * cos((M-1)*THTJ)       + fbnd(2*M) * sin((M-1)*THTJ) )
       drm  =(M-1)*radius**(M-2) * ( fbnd(2*M-1) * cos((M-1)*THTJ)       + fbnd(2*M) * sin((M-1)*THTJ))
       drmt =      radius**(M-1) * (-fbnd(2*M-1) * (M-1)*sin((M-1)*THTJ) + fbnd(2*M) *(M-1)*cos((M-1)*THTJ))
       drmtr=(M-1)*radius**(M-2) * (-fbnd(2*M-1) * (M-1)*sin((M-1)*THTJ) + fbnd(2*M) *(M-1)*cos((M-1)*THTJ))
     endif
     RR(1,node) = RR(1,node) + amin * rm  * cos(thtj)
     ZZ(1,node) = ZZ(1,node) + amin * rm  * sin(thtj)
     RR(2,node) = RR(2,node) + amin * drm * cos(thtj)
     ZZ(2,node) = ZZ(2,node) + amin * drm * sin(thtj)
     RR(3,node) = RR(3,node) - amin * rm  * sin(thtj) + amin * drmt  * cos(thtj)
     ZZ(3,node) = ZZ(3,node) + amin * rm  * cos(thtj) + amin * drmt  * sin(thtj)
     RR(4,node) = RR(4,node) - amin * drm * sin(thtj) + amin * drmtr * cos(thtj)
     ZZ(4,node) = ZZ(4,node) + amin * drm * cos(thtj) + amin * drmtr * sin(thtj)

     PSI(1,node) = PSI(1,node) + radius**8              * (   fpsi(2*m-1) *            cos((m-1)*thtj)       &
                                                            + fpsi(2*m)   *            sin((m-1)*thtj) )
     PSI(2,node) = PSI(2,node) + 8.d0 * radius**7 * drm * (   fpsi(2*m-1) *            cos((m-1)*thtj)       &
                                                            + fpsi(2*m)   *            sin((m-1)*thtj) )
     PSI(3,node) = PSI(3,node) + radius**8              * ( - fpsi(2*m-1) * float(m-1)*sin((m-1)*thtj)       &
                                                            + fpsi(2*m)   * float(m-1)*cos((m-1)*thtj) )
     PSI(4,node) = PSI(4,node) + 8.d0 * radius**7 * drm * ( - fpsi(2*m-1) * float(m-1)*sin((m-1)*thtj)       &
                                                            + fpsi(2*m)   * float(m-1)*cos((m-1)*thtj) )

!     write(*,'(A,2i6,12e16.8)') ' PSI(1,node) : ', node,m,PSI(1,node),radius,thtj,fpsi(2*m-1),fpsi(2*m)

   enddo

   RR(2,node)  = RR(2,node)  * ds/2.d0 * dr_ds
   RR(3,node)  = RR(3,node)  * dt/2.d0 * dtht_dt
   RR(4,node)  = RR(4,node)  * ds/2.d0 * dt/2.d0 * dr_ds * dtht_dt
   ZZ(2,node)  = ZZ(2,node)  * ds/2.d0 * dr_ds
   ZZ(3,node)  = ZZ(3,node)  * dt/2.d0 * dtht_dt
   ZZ(4,node)  = ZZ(4,node)  * ds/2.d0 * dt/2.d0 * dr_ds * dtht_dt
   PSI(2,node) = PSI(2,node) * ds/2.d0 * dr_ds
   PSI(3,node) = PSI(3,node) * dt/2.d0 * dtht_dt
   PSI(4,node) = PSI(4,node) * ds/2.d0 * dt/2.d0 * dr_ds * dtht_dt

!   write(*,*) node, PSI(1,node)

 enddo
enddo

element_list%n_elements  = n_element_start  + (nr-1)*np
node_list%n_nodes        = n_node_start     + nr*np

do i=1,nr-1

 do j=1,np
   node  = np*(i-1) + j
   index = n_element_start + node
   element_list%element(index)%vertex(1) = n_node_start + (i-1)*np + j
   element_list%element(index)%vertex(4) = n_node_start + (i-1)*np + j + 1
   element_list%element(index)%vertex(3) = n_node_start + (i  )*np + j + 1
   element_list%element(index)%vertex(2) = n_node_start + (i  )*np + j

   if (j .eq. np) then          
                 element_list%element(Index)%vertex(4) = (i-1)*np + 1
                 element_list%element(Index)%vertex(3) =  i*np    + 1
   endif    





           !Neighbours of the element (refinement procedure)

	    if(i==1) then	        
	         element_list%element(Index)%neighbours(4) = 0    
              else		 
                 element_list%element(Index)%neighbours(4) = Index - np 
	    end if 	 
	    
	    if(j==np) then
	    	 element_list%element(Index)%neighbours(3) = Index - np + 1  
	      else
	    	 element_list%element(Index)%neighbours(3) = Index + 1       
	    end if	  	    
	    
	    if(i==nr-1) then
	    	 element_list%element(Index)%neighbours(2) = 0   
	      else   
	    	 element_list%element(Index)%neighbours(2) = Index + np 
	    end if 
	        
	    if(j==1) then
	    	 element_list%element(Index)%neighbours(1) = Index + np -1 
	      else   
	    	 element_list%element(Index)%neighbours(1) = Index -1       
	    end if     
	  
            ! Initialization of the genealogy  (refinement procedure)

            element_list%element(Index)%father = 0
	    element_list%element(Index)%n_sons = 0
	    do i_sons = 1, 4
	         element_list%element(Index)%sons(i_sons) = 0
	    end do 
  enddo
enddo

 


!-------------------- translate cubic Hermite to Bezier parameters
!
!  type type_node                                      ! type definition of a node (i.e. a vertex)
!    real*8    :: x(n_order+1,ndim)                      ! x,y coordinates of points and additional nodal geometry
!    integer :: boundary                               ! = 1 for boundary nodes
!  endtype type_node                                   ! x(:,1) : position, x(:,2) : vector u, x(:,3) : vector v, x(4) : vector w
!
!  type type_element
!    integer :: vertex(n_vertex_max)
!    integer :: neighbours(n_vertex_max)
!    real*8    :: size(n_vertex_max,n_order+1)
!  endtype type_element
!-----------------------------------------------------------

do i=1,nr

 do j=1,np

   angle = 2.d0*PI * float(j-1)/float(np)

   index0 = np*(i-1) + j
   index  = n_node_start + np*(i-1) + j

   node_list%node(index)%X(1,1)        = RR(1,index0)
   node_list%node(index)%X(1,2)        = ZZ(1,index0)
   node_list%node(index)%values(1,1,1) = PSI(1,index0)

   node_list%node(index)%X(2,1)        = RR(2,index0)  * 2.d0/3.d0
   node_list%node(index)%X(2,2)        = ZZ(2,index0)  * 2.d0/3.d0
   node_list%node(index)%values(1,2,1) = PSI(2,index0) * 2.d0/3.d0

   node_list%node(index)%X(3,1)        = RR(3,index0)  * 2.d0/3.d0
   node_list%node(index)%X(3,2)        = ZZ(3,index0)  * 2.d0/3.d0
   node_list%node(index)%values(1,3,1) = PSI(3,index0) * 2.d0/3.d0

   node_list%node(index)%X(4,1)        = RR(4,index0)  * 4.d0/9.d0
   node_list%node(index)%X(4,2)        = ZZ(4,index0)  * 4.d0/9.d0
   node_list%node(index)%values(1,4,1) = PSI(4,index0) * 4.d0/9.d0

   if (i .eq. nr) node_list%node(index)%boundary = 2

   do k=1,n_order+1
     node_list%node(index)%index(k) = n_index_start + (n_order+1)*(index0-1)+k
   enddo
  
   node_list%node(index)%constrained=.false.
 enddo

enddo


do k=n_element_start+1 , element_list%n_elements   ! fill in the size of the elements

 do iv = 1, 4                    ! over 4 corners of an element

   ivp = mod(iv,4)   + 1         ! vertex with index one higher
   ivm = mod(iv+2,4) + 1         ! vertex with index one below

   node_iv  = element_list%element(k)%vertex(iv)
   node_ivp = element_list%element(k)%vertex(ivp)
   node_ivm = element_list%element(k)%vertex(ivm)

   if ((iv .eq. 1) .or. (iv .eq.3)) then

     delta_Rp = node_list%node(node_ivp)%X(1,1) - node_list%node(node_iv)%X(1,1)
     delta_Zp = node_list%node(node_ivp)%X(1,2) - node_list%node(node_iv)%X(1,2)
     dir_2    = delta_Rp * node_list%node(node_iv)%X(2,1) + delta_Zp * node_list%node(node_iv)%X(2,2)

     delta_Rm = node_list%node(node_ivm)%X(1,1) - node_list%node(node_iv)%X(1,1)
     delta_Zm = node_list%node(node_ivm)%X(1,2) - node_list%node(node_iv)%X(1,2)
     dir_3    = delta_Rm * node_list%node(node_iv)%X(3,1) + delta_Zm * node_list%node(node_iv)%X(3,2)

   else

     delta_Rp = node_list%node(node_ivp)%X(1,1) - node_list%node(node_iv)%X(1,1)
     delta_Zp = node_list%node(node_ivp)%X(1,2) - node_list%node(node_iv)%X(1,2)
     dir_3    = delta_Rp * node_list%node(node_iv)%X(3,1) + delta_Zp * node_list%node(node_iv)%X(3,2)

     delta_Rm = node_list%node(node_ivm)%X(1,1) - node_list%node(node_iv)%X(1,1)
     delta_Zm = node_list%node(node_ivm)%X(1,2) - node_list%node(node_iv)%X(1,2)
     dir_2    = delta_Rm * node_list%node(node_iv)%X(2,1) + delta_Zm * node_list%node(node_iv)%X(2,2)

   endif

   if (dir_2 .ne. 0.d0) then
     dir_2 = dir_2 / abs(dir_2)
   else
     dir_2 = 1.d0
   endif
   if (dir_3 .ne. 0.d0) then
     dir_3 = dir_3 / abs(dir_3)
   else
     dir_3 = -1.d0
     if (iv.eq.1) dir_3 = 1.d0              ! admittedly not very elegant
   endif

   element_list%element(k)%size(iv,1) = 1.d0
   element_list%element(k)%size(iv,2) = dir_2
   element_list%element(k)%size(iv,3) = dir_3
   element_list%element(k)%size(iv,4) = element_list%element(k)%size(iv,2) * element_list%element(k)%size(iv,3)

!   if ((RR(2,node_iv)**2 + ZZ(2,node_iv)**2) .eq. 0.) element_list%element(k)%size(iv,2) = dir_2
!   if ((RR(3,node_iv)**2 + ZZ(3,node_iv)**2) .eq. 0.) element_list%element(k)%size(iv,3) = dir_3
!   if ((RR(4,node_iv)**2 + ZZ(4,node_iv)**2) .eq. 0.) element_list%element(k)%size(iv,4) = dir_2 * dir_3

!    write(*,'(2i5,12e16.8)') k,iv,element_list%element(k)%size(iv,1:4)

 enddo

enddo

return
end

