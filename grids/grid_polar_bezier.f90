!> Defines a polar grid using Bezier finite elements (using the HELENA
!! cubic Hermite elements formulation)
subroutine grid_polar_bezier(Rgeo,Zgeo,amin,acentre,angle_start,fbnd,fpsi,mf,nr,np,node_list,element_list)

use constants
use tr_module
use mod_parameters
use data_structure
use mod_neighbours, only: update_neighbours
use phys_module, only: psi_axis_init, XR_r, SIG_r, XR_tht, SIG_tht, fix_axis_nodes, force_central_node
use mod_node_indices

implicit none

! --- Routine parameters
real*8,                  intent(in)    :: Rgeo           !< R-position of geometric center
real*8,                  intent(in)    :: Zgeo           !< Z-position of geometric center
real*8,                  intent(in)    :: amin           !< minor radius
real*8,                  intent(in)    :: acentre        !< inner radius of grid
real*8,                  intent(in)    :: angle_start    !< poloidal angle of first element
real*8,                  intent(in)    :: fbnd(*)        !< Fourier series describing the radius as function
                                                         !!   of the poloidal angle
real*8,                  intent(in)    :: fpsi(*)        !< Fourier series of flux at the boundary
integer,                 intent(in)    :: mf             !< Number of Fourier modes in fbnd and fpsi
integer,                 intent(in)    :: nr             !< number of radial points, (nr-1) elements
integer,                 intent(in)    :: np             !< number of poloidal points, np elements
type(type_node_list),    intent(inout) :: node_list      !< list of nodes with grid information
type(type_element_list), intent(inout) :: element_list   !< list of elements with element information

! --- local variables
real*8 :: acentre2, radius2
real*8              :: si
real*8, allocatable :: RR(:,:),ZZ(:,:),PSI(:,:)
real*8              :: dt, ds, thtj, radius
real*8              :: rm, drm_r, drm_rr, FFpsi, dFFpsi_t, dFFpsi_tt, GGpsi, dGGpsi_t, dGGpsi_tt
real*8              :: angle, psi_axis
real*8              :: delta_rm, delta_zm, delta_rp, delta_zp, dir_2, dir_3, size_ratio
integer             :: i, j, m, index, index0, node, k, iv, ivp, ivm, node_iv, node_ivp, node_ivm,i_sons
integer             :: n_element_start, n_node_start, n_index_start, n_max
real*8              :: abltg(3), dr_ds, dr2_ds, dr_ds2, dr2_ds2, dtht_dt, dtht_dt2
real*8, allocatable :: S1(:), S2(:), SP1(:), SP2(:), SP3(:), SP4(:)
real*8, allocatable :: T1(:), T2(:), TP1(:), TP2(:), TP3(:), TP4(:)
real*8, external    :: spwert
logical             :: skip_update_neighbours
logical             :: doing_polar_square

integer :: node_indices( (n_order+1)/2, (n_order+1)/2 ), kk, ll, index_tmp

n_max = n_degrees
if (n_order .gt. 5) n_max = (5+1)**2 / 4 ! we don't care about derivatives >= 3...
call tr_allocate(RR,1,n_max,1,nr*np,"RR",CAT_GRID)
call tr_allocate(ZZ,1,n_max,1,nr*np,"ZZ",CAT_GRID)
call tr_allocate(PSI,1,n_max,1,nr*np,"PSI",CAT_GRID)

dt = 2.d0*pi/real(np)
ds = 1.d0/real(nr-1)

n_element_start  = element_list%n_elements
n_node_start     = node_list%n_nodes
do i=n_node_start+1,n_nodes_max
  node_list%node(i)%x        = 0.d0
  node_list%node(i)%values   = 0.d0
  node_list%node(i)%index    = 0
  node_list%node(i)%boundary = 0
enddo

skip_update_neighbours = .false.
if ( n_element_start /= 0 ) skip_update_neighbours = .true. ! In such a case, the call to update_neighbours needs to be done in another routine like grid_bezier_square_polar

do i=n_element_start+1,n_elements_max
  element_list%element(i)%vertex     = 0
  element_list%element(i)%size       = 0.d0
  element_list%element(i)%neighbours = 0
enddo

n_index_start = 0
do i=1,n_node_start
  n_index_start = max(n_index_start,maxval(node_list%node(i)%index(:)))
enddo
doing_polar_square = .false.
if (n_index_start .gt. 0) doing_polar_square = .true.

write(*,*) '*************************************'
write(*,*) '*        grid_polar_bezier          *'
write(*,*) '*************************************'
write(*,*) ' existing number of elements          : ',n_element_start
write(*,*) ' existing number of nodes             : ',n_node_start
write(*,*) ' index_start                          : ',n_index_start

psi_axis = psi_axis_init ! Initial guess for Psi at the magnetic axis

call tr_allocate(S1,1,nr,"S1",CAT_GRID)
call tr_allocate(S2,1,nr,"S2",CAT_GRID)
call tr_allocate(SP1,1,nr,"SP1",CAT_GRID)
call tr_allocate(SP2,1,nr,"SP2",CAT_GRID)
call tr_allocate(SP3,1,nr,"SP3",CAT_GRID)
call tr_allocate(SP4,1,nr,"SP4",CAT_GRID)
S2 = 0.d0
call tr_allocate(T1,1,np+1,"T1",CAT_GRID)
call tr_allocate(T2,1,np+1,"T2",CAT_GRID)
call tr_allocate(TP1,1,np+1,"TP1",CAT_GRID)
call tr_allocate(TP2,1,np+1,"TP2",CAT_GRID)
call tr_allocate(TP3,1,np+1,"TP3",CAT_GRID)
call tr_allocate(TP4,1,np+1,"TP4",CAT_GRID)
T2 = 0.d0
do i=1,nr
  S1(i) = real(i-1)/real(nr-1)
enddo
do j=1,np+1
  T1(j) = real(j-1)/real(np)
enddo

call meshac2(nr,S2,XR_r(1),XR_r(2),SIG_r(1),SIG_r(2),0.6d0,1.0d0)
call spline(nr,S1,S2,0.d0,0.d0,2,SP1,SP2,SP3,SP4)

call meshac2(np+1,T2,XR_tht(1),XR_tht(2),SIG_tht(1),SIG_tht(2),0.6d0,1.0d0)
call spline(np+1,T1,T2,0.d0,0.d0,2,TP1,TP2,TP3,TP4)

acentre2 = acentre * 2.d0 / fbnd(1) !###

do i=1,nr

  si = spwert(nr,S1(i),SP1,SP2,SP3,SP4,S1,ABLTG)

  radius  = ( acentre2 + (1.d0-acentre2) * si )
  dr_ds   = (1.d0-acentre2) * abltg(1)
  dr_ds2  = (1.d0-acentre2) * abltg(2)
  radius2 = si
  dr2_ds  = abltg(1)
  dr2_ds2 = abltg(2)
  
  do  j=1,np
  
    node   = np*(i-1) + j

    thtj     = angle_start + spwert(np+1,T1(j),TP1,TP2,TP3,TP4,T1,ABLTG) * 2.d0 * PI
    dtht_dt  = abltg(1)
    dtht_dt2 = abltg(2)
    
    rm     = radius
    drm_r  = dr_ds
    drm_rr = dr_ds2

    ! --- R for first harmonic
    FFpsi     =   fbnd(1) * cos(thtj) / 2.d0
    dFFpsi_t  = - fbnd(1) * sin(thtj) / 2.d0 * dtht_dt
    dFFpsi_tt = - fbnd(1) * sin(thtj) / 2.d0 * dtht_dt2 - fbnd(1) * cos(thtj) / 2.d0 * dtht_dt**2

    RR(1,node)   = Rgeo + amin * rm     * FFpsi
    RR(2,node)   =        amin * drm_r  * FFpsi
    RR(3,node)   =        amin * rm     * dFFpsi_t
    RR(4,node)   =        amin * drm_r  * dFFpsi_t
    if (n_order .ge. 5) then
      RR(5,node) =        amin * drm_rr * FFpsi
      RR(6,node) =        amin * rm     * dFFpsi_tt
      RR(7,node) =        amin * drm_rr * dFFpsi_t
      RR(8,node) =        amin * drm_r  * dFFpsi_tt
      RR(9,node) =        amin * drm_rr * dFFpsi_tt
    endif

    FFpsi     = fbnd(1) * sin(thtj) / 2.d0
    dFFpsi_t  = fbnd(1) * cos(thtj) / 2.d0 * dtht_dt
    dFFpsi_tt = fbnd(1) * cos(thtj) / 2.d0 * dtht_dt2 - fbnd(1) * sin(thtj) / 2.d0 * dtht_dt**2

    ! --- Z for first harmonic
    ZZ(1,node)   = Zgeo + amin * rm     * FFpsi
    ZZ(2,node)   =        amin * drm_r  * FFpsi
    ZZ(3,node)   =        amin * rm     * dFFpsi_t
    ZZ(4,node)   =        amin * drm_r  * dFFpsi_t
    if (n_order .ge. 5) then
      ZZ(5,node) =        amin * drm_rr * FFpsi
      ZZ(6,node) =        amin * rm     * dFFpsi_tt
      ZZ(7,node) =        amin * drm_rr * dFFpsi_t
      ZZ(8,node) =        amin * drm_r  * dFFpsi_tt
      ZZ(9,node) =        amin * drm_rr * dFFpsi_tt
    endif

    ! --- PSI for first harmonic
    rm     =             radius**8 * fpsi(1) / 2.d0 + psi_axis*(1.d0 - radius**2)
    drm_r  = (8.d0     * radius**7 * fpsi(1) / 2.d0 - psi_axis* 2.d0 * radius   ) * dr_ds
    drm_rr = (8.d0     * radius**7 * fpsi(1) / 2.d0 - psi_axis* 2.d0 * radius   ) * dr_ds2 &
            +(8.d0*7.d0* radius**6 * fpsi(1) / 2.d0 - psi_axis* 2.d0            ) * dr_ds**2

    PSI(1,node)   = rm
    PSI(2,node)   = drm_r
    PSI(3,node)   = 0.d0
    PSI(4,node)   = 0.d0
    if (n_order .ge. 5) then
      PSI(5,node) = drm_rr
      PSI(6,node) = 0.d0
      PSI(7,node) = 0.d0
      PSI(8,node) = 0.d0
      PSI(9,node) = 0.d0
    endif

!   write(*,'(A,2i6,12e16.8)') ' PSI(1,node) : ', node,1,PSI(1,node),radius,thtj,fpsi(1),fpsi(2)

    do m = 2, mf/2

      if (m .eq. 2) then
        rm     = radius2
        drm_r  = dr2_ds
        drm_rr = dr2_ds2
      else
        rm     =             radius2**(m-1)
        drm_r  = (m-1)      *radius2**(m-2) * dr2_ds
        drm_rr = (m-1)      *radius2**(m-2) * dr2_ds2 &
                +(m-1)*(m-2)*radius2**(m-3) * dr2_ds**2
      endif

      FFpsi     = ( fbnd(2*m-1)           *cos((m-1)*thtj) + fbnd(2*m)          *sin((m-1)*thtj) )
      dFFpsi_t  = (-fbnd(2*m-1) * (m-1)   *sin((m-1)*thtj) + fbnd(2*m) *(m-1)   *cos((m-1)*thtj) ) * dtht_dt
      dFFpsi_tt = (-fbnd(2*m-1) * (m-1)   *sin((m-1)*thtj) + fbnd(2*m) *(m-1)   *cos((m-1)*thtj) ) * dtht_dt2 &
                 +(-fbnd(2*m-1) * (m-1)**2*cos((m-1)*thtj) - fbnd(2*m) *(m-1)**2*sin((m-1)*thtj) ) * dtht_dt**2

      ! --- R for higher harmonics
      GGpsi     = FFpsi     * cos(thtj)
      dGGpsi_t  = dFFpsi_t  * cos(thtj) - FFpsi * sin(thtj) * dtht_dt
      dGGpsi_tt = dFFpsi_tt * cos(thtj) - 2.0 * dFFpsi_t * sin(thtj) * dtht_dt &
                 - FFpsi * cos(thtj) * dtht_dt**2 - FFpsi * sin(thtj) * dtht_dt2

      RR(1,node)   = RR(1,node) + amin * rm     *  GGpsi
      RR(2,node)   = RR(2,node) + amin * drm_r  *  GGpsi
      RR(3,node)   = RR(3,node) + amin * rm     * dGGpsi_t
      RR(4,node)   = RR(4,node) + amin * drm_r  * dGGpsi_t
      if (n_order .ge. 5) then
        RR(5,node) = RR(5,node) + amin * drm_rr *  GGpsi
        RR(6,node) = RR(6,node) + amin * rm     * dGGpsi_tt
        RR(7,node) = RR(7,node) + amin * drm_rr * dGGpsi_t
        RR(8,node) = RR(8,node) + amin * drm_r  * dGGpsi_tt
        RR(9,node) = RR(9,node) + amin * drm_rr * dGGpsi_tt
      endif

      ! --- Z for higher harmonics
      GGpsi     = FFpsi     * sin(thtj)
      dGGpsi_t  = dFFpsi_t  * sin(thtj) + FFpsi * cos(thtj) * dtht_dt
      dGGpsi_tt = dFFpsi_tt * sin(thtj) + 2.0 * dFFpsi_t * cos(thtj) * dtht_dt &
                 - FFpsi * sin(thtj) * dtht_dt**2 + FFpsi * cos(thtj) * dtht_dt2

      ZZ(1,node)   = ZZ(1,node) + amin * rm     *  GGpsi
      ZZ(2,node)   = ZZ(2,node) + amin * drm_r  *  GGpsi
      ZZ(3,node)   = ZZ(3,node) + amin * rm     * dGGpsi_t
      ZZ(4,node)   = ZZ(4,node) + amin * drm_r  * dGGpsi_t
      if (n_order .ge. 5) then                              
        ZZ(5,node) = ZZ(5,node) + amin * drm_rr *  GGpsi
        ZZ(6,node) = ZZ(6,node) + amin * rm     * dGGpsi_tt
        ZZ(7,node) = ZZ(7,node) + amin * drm_rr * dGGpsi_t
        ZZ(8,node) = ZZ(8,node) + amin * drm_r  * dGGpsi_tt
        ZZ(9,node) = ZZ(9,node) + amin * drm_rr * dGGpsi_tt
      endif

      ! --- PSI for higher harmonics
      rm     =               radius**8
      drm_r  = 8.d0        * radius**7 * dr_ds
      drm_rr = 8.d0        * radius**7 * dr_ds2 &
              +8.d0 * 7.d0 * radius**6 * dr_ds**2

      FFpsi     =   fpsi(2*m-1)                *cos((m-1)*thtj) + fpsi(2*m)                *sin((m-1)*thtj)
      dFFpsi_t  = - fpsi(2*m-1) * float(m-1)   *sin((m-1)*thtj) + fpsi(2*m) * float(m-1)   *cos((m-1)*thtj)
      dFFpsi_tt = - fpsi(2*m-1) * float(m-1)**2*cos((m-1)*thtj) - fpsi(2*m) * float(m-1)**2*sin((m-1)*thtj)

      PSI(1,node)   = PSI(1,node) + rm     *  FFpsi
      PSI(2,node)   = PSI(2,node) + drm_r  *  FFpsi
      PSI(3,node)   = PSI(3,node) + rm     * dFFpsi_t
      PSI(4,node)   = PSI(4,node) + drm_r  * dFFpsi_t
      if (n_order .ge. 5) then
        PSI(5,node) = PSI(5,node) + drm_rr *  FFpsi
        PSI(6,node) = PSI(6,node) + rm     * dFFpsi_tt
        PSI(7,node) = PSI(7,node) + drm_rr * dFFpsi_t
        PSI(8,node) = PSI(8,node) + drm_r  * dFFpsi_tt
        PSI(9,node) = PSI(9,node) + drm_rr * dFFpsi_tt
      endif

    enddo

    ! --- Include ds and dt for all derivatives
    RR(2,node)   = RR(2,node) * ds/2.d0
    RR(3,node)   = RR(3,node) * dt/2.d0
    RR(4,node)   = RR(4,node) * ds/2.d0 * dt/2.d0
    if (n_order .ge. 5) then
      RR(5,node) = RR(5,node) * ds/2.d0 * ds/2.d0
      RR(6,node) = RR(6,node) * dt/2.d0 * dt/2.d0
      RR(7,node) = RR(7,node) * ds/2.d0 * ds/2.d0 * dt/2.d0
      RR(8,node) = RR(8,node) * ds/2.d0 * dt/2.d0 * dt/2.d0
      RR(9,node) = RR(9,node) * ds/2.d0 * ds/2.d0 * dt/2.d0 * dt/2.d0
    endif

    ZZ(2,node)   = ZZ(2,node) * ds/2.d0
    ZZ(3,node)   = ZZ(3,node) * dt/2.d0
    ZZ(4,node)   = ZZ(4,node) * ds/2.d0 * dt/2.d0
    if (n_order .ge. 5) then
      ZZ(5,node) = ZZ(5,node) * ds/2.d0 * ds/2.d0
      ZZ(6,node) = ZZ(6,node) * dt/2.d0 * dt/2.d0
      ZZ(7,node) = ZZ(7,node) * ds/2.d0 * ds/2.d0 * dt/2.d0
      ZZ(8,node) = ZZ(8,node) * ds/2.d0 * dt/2.d0 * dt/2.d0
      ZZ(9,node) = ZZ(9,node) * ds/2.d0 * ds/2.d0 * dt/2.d0 * dt/2.d0
    endif

    PSI(2,node)   = PSI(2,node) * ds/2.d0
    PSI(3,node)   = PSI(3,node) * dt/2.d0
    PSI(4,node)   = PSI(4,node) * ds/2.d0 * dt/2.d0
    if (n_order .ge. 5) then
      PSI(5,node) = PSI(5,node) * ds/2.d0 * ds/2.d0
      PSI(6,node) = PSI(6,node) * dt/2.d0 * dt/2.d0
      PSI(7,node) = PSI(7,node) * ds/2.d0 * ds/2.d0 * dt/2.d0
      PSI(8,node) = PSI(8,node) * ds/2.d0 * dt/2.d0 * dt/2.d0
      PSI(9,node) = PSI(9,node) * ds/2.d0 * ds/2.d0 * dt/2.d0 * dt/2.d0
    endif

  enddo
enddo

element_list%n_elements  = n_element_start  + (nr-1)*np
node_list%n_nodes        = n_node_start     + nr*np

if ( node_list%n_nodes > n_nodes_max ) then
  write(*,*) 'ERROR in grid_polar_bezier: hard-coded parameter n_nodes_max is too small'
  stop
else if ( element_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in grid_polar_bezier: hard-coded parameter n_elements_max is too small'
  stop
end if


do i=1,nr-1

 do j=1,np
   node  = np*(i-1) + j
   index = n_element_start + node
   element_list%element(index)%vertex(1) = n_node_start + (i-1)*np + j
   element_list%element(index)%vertex(4) = n_node_start + (i-1)*np + j + 1
   element_list%element(index)%vertex(3) = n_node_start + (i  )*np + j + 1
   element_list%element(index)%vertex(2) = n_node_start + (i  )*np + j

   if (j .eq. np) then          
                 element_list%element(Index)%vertex(4) = n_node_start + (i-1)*np + 1
                 element_list%element(Index)%vertex(3) = n_node_start + i*np    + 1
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
!    real*8    :: x(n_degrees,ndim)                      ! x,y coordinates of points and additional nodal geometry
!    integer :: boundary                               ! = 1 for boundary nodes
!  endtype type_node                                   ! x(:,1) : position, x(:,2) : vector u, x(:,3) : vector v, x(4) : vector w
!
!  type type_element
!    integer :: vertex(n_vertex_max)
!    integer :: neighbours(n_vertex_max)
!    real*8    :: size(n_vertex_max,n_degrees)
!  endtype type_element
!-----------------------------------------------------------

do i=1,nr

 do j=1,np

   angle = 2.d0*PI * float(j-1)/float(np)

   index0 = np*(i-1) + j
   index  = n_node_start + np*(i-1) + j

   ! --- Note: factor 2.0 is because of derivatives ds/2 and dt/2 above
   size_ratio = 2.d0/float(n_order)

   node_list%node(index)%X(1,:,1)        = 0.d0
   node_list%node(index)%X(1,1,1)        = RR(1,index0)
   node_list%node(index)%X(1,2,1)        = RR(2,index0)  * size_ratio
   node_list%node(index)%X(1,3,1)        = RR(3,index0)  * size_ratio
   node_list%node(index)%X(1,4,1)        = RR(4,index0)  * size_ratio**2
   if (n_order .ge. 5) then                              
     node_list%node(index)%X(1,5,1)      = RR(5,index0)  * size_ratio**2
     node_list%node(index)%X(1,6,1)      = RR(6,index0)  * size_ratio**2
     node_list%node(index)%X(1,7,1)      = RR(7,index0)  * size_ratio**3
     node_list%node(index)%X(1,8,1)      = RR(8,index0)  * size_ratio**3
     node_list%node(index)%X(1,9,1)      = RR(9,index0)  * size_ratio**4
   endif                                                 
                                                         
   node_list%node(index)%X(1,:,2)        = 0.d0
   node_list%node(index)%X(1,1,2)        = ZZ(1,index0) 
   node_list%node(index)%X(1,2,2)        = ZZ(2,index0)  * size_ratio
   node_list%node(index)%X(1,3,2)        = ZZ(3,index0)  * size_ratio
   node_list%node(index)%X(1,4,2)        = ZZ(4,index0)  * size_ratio**2
   if (n_order .ge. 5) then                              
     node_list%node(index)%X(1,5,2)      = ZZ(5,index0)  * size_ratio**2
     node_list%node(index)%X(1,6,2)      = ZZ(6,index0)  * size_ratio**2
     node_list%node(index)%X(1,7,2)      = ZZ(7,index0)  * size_ratio**3
     node_list%node(index)%X(1,8,2)      = ZZ(8,index0)  * size_ratio**3
     node_list%node(index)%X(1,9,2)      = ZZ(9,index0)  * size_ratio**4
   endif

   node_list%node(index)%values(1,:,1)   = 0.d0
   node_list%node(index)%values(1,1,1)   = PSI(1,index0)
   node_list%node(index)%values(1,2,1)   = PSI(2,index0) * size_ratio
   node_list%node(index)%values(1,3,1)   = PSI(3,index0) * size_ratio
   node_list%node(index)%values(1,4,1)   = PSI(4,index0) * size_ratio**2
   if (n_order .ge. 5) then
     node_list%node(index)%values(1,5,1) = PSI(5,index0) * size_ratio**2
     node_list%node(index)%values(1,6,1) = PSI(6,index0) * size_ratio**2
     node_list%node(index)%values(1,7,1) = PSI(7,index0) * size_ratio**3
     node_list%node(index)%values(1,8,1) = PSI(8,index0) * size_ratio**3
     node_list%node(index)%values(1,9,1) = PSI(9,index0) * size_ratio**4
   endif

   node_list%node(index)%boundary = 0
   if (i .eq. nr) node_list%node(index)%boundary = 2

   node_list%node(index)%axis_node = .false.
   if ( fix_axis_nodes .and. (.not. doing_polar_square) .and. (i .eq. 1) ) node_list%node(index)%axis_node = .true.

   if (force_central_node .and. (.not. doing_polar_square) .and. (i.eq.1)) then

     node_list%node(index)%index(1) = 1

     if (j.eq.1) n_index_start = n_index_start + 1

     do k=2,n_degrees
       node_list%node(index)%index(k) = n_index_start + k-1
     enddo
     n_index_start = n_index_start +n_degrees-1

   else
     do k=1,n_degrees
       node_list%node(index)%index(k) = n_index_start + k
     enddo
     n_index_start = n_index_start + n_degrees
   endif
  
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

     delta_Rp = node_list%node(node_ivp)%X(1,1,1) - node_list%node(node_iv)%X(1,1,1)
     delta_Zp = node_list%node(node_ivp)%X(1,1,2) - node_list%node(node_iv)%X(1,1,2)
     dir_2    = delta_Rp * node_list%node(node_iv)%X(1,2,1) + delta_Zp * node_list%node(node_iv)%X(1,2,2)

     delta_Rm = node_list%node(node_ivm)%X(1,1,1) - node_list%node(node_iv)%X(1,1,1)
     delta_Zm = node_list%node(node_ivm)%X(1,1,2) - node_list%node(node_iv)%X(1,1,2)
     dir_3    = delta_Rm * node_list%node(node_iv)%X(1,3,1) + delta_Zm * node_list%node(node_iv)%X(1,3,2)

   else

     delta_Rp = node_list%node(node_ivp)%X(1,1,1) - node_list%node(node_iv)%X(1,1,1)
     delta_Zp = node_list%node(node_ivp)%X(1,1,2) - node_list%node(node_iv)%X(1,1,2)
     dir_3    = delta_Rp * node_list%node(node_iv)%X(1,3,1) + delta_Zp * node_list%node(node_iv)%X(1,3,2)

     delta_Rm = node_list%node(node_ivm)%X(1,1,1) - node_list%node(node_iv)%X(1,1,1)
     delta_Zm = node_list%node(node_ivm)%X(1,1,2) - node_list%node(node_iv)%X(1,1,2)
     dir_2    = delta_Rm * node_list%node(node_iv)%X(1,2,1) + delta_Zm * node_list%node(node_iv)%X(1,2,2)

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
   if (n_order .ge. 5) then
     element_list%element(k)%size(iv,5) = 1.d0
     element_list%element(k)%size(iv,6) = 1.d0
     element_list%element(k)%size(iv,7) = dir_3 ! because \vec{m} is similar to \vec{v}
     element_list%element(k)%size(iv,8) = dir_2 ! while   \vec{n} is similar to \vec{u}
     element_list%element(k)%size(iv,9) = element_list%element(k)%size(iv,5) * element_list%element(k)%size(iv,6)
   endif
   if (n_order .gt. 5) then
     do j=10,n_degrees
       ! --- Just use a random size for higher orders, but one we're sure is the same for all elements
       element_list%element(k)%size(iv,j) = element_list%element(k)%size(iv,5) * element_list%element(k)%size(iv,6)
     enddo
     ! --- calculate node_indices
     call calculate_node_indices(node_indices)
     ! --- Loop over all indices
     do kk = 1,(n_order+1)/2
       do ll = 1,(n_order+1)/2
         index_tmp = node_indices(kk,ll)
         if (index_tmp .le. 9) cycle ! we want only derivatives >=3
         if ( (kk .ne. 2) .and. (ll .ne. 2) ) cycle ! we want only the nodes that are next to boundary (the only ones with h_ij = - h_-ij)
         if (kk .eq. 2) then
           element_list%element(k)%size(iv,index_tmp) = dir_2
         endif
         if (ll .eq. 2) then
           element_list%element(k)%size(iv,index_tmp) = dir_3
         endif
       enddo
     enddo
   endif
   if (fix_axis_nodes) then
      j = element_list%element(k)%vertex(iv)
      if (node_list%node(j)%axis_node) then
        element_list%element(k)%size(iv,3) = 0.d0
        element_list%element(k)%size(iv,4) = 0.d0
        if (n_order .ge. 5) then
          element_list%element(k)%size(iv,6) = 0.d0
          element_list%element(k)%size(iv,7) = 0.d0
          element_list%element(k)%size(iv,8) = 0.d0
          element_list%element(k)%size(iv,9) = 0.d0
        endif
        if (n_order .gt. 5) then
          ! --- Loop over all indices (note node_indices have already been calculated above...)
          do kk = 1,(n_order+1)/2
            do ll = 2,(n_order+1)/2 ! start t-index from 2 to keep only the pure s-derivatives
              index_tmp = node_indices(kk,ll)
              if (index_tmp .le. 9) cycle ! we want only derivatives >=3
              element_list%element(k)%size(iv,index_tmp) = 0.d0
            enddo
          enddo
        endif
      endif
   endif

!   if ((RR(2,node_iv)**2 + ZZ(2,node_iv)**2) .eq. 0.) element_list%element(k)%size(iv,2) = dir_2
!   if ((RR(3,node_iv)**2 + ZZ(3,node_iv)**2) .eq. 0.) element_list%element(k)%size(iv,3) = dir_3
!   if ((RR(4,node_iv)**2 + ZZ(4,node_iv)**2) .eq. 0.) element_list%element(k)%size(iv,4) = dir_2 * dir_3

!    write(*,'(2i5,12e16.8)') k,iv,element_list%element(k)%size(iv,1:4)

 enddo

enddo

if ( .not. skip_update_neighbours ) call update_neighbours(node_list,element_list, force_rtree_initialize=.true.)

return
end subroutine grid_polar_bezier


















subroutine transform_to_bi_quintic(node_list,element_list)


use constants
use tr_module
use mod_parameters
use data_structure
use phys_module, only: R_geo, Z_geo

implicit none

! --- Routine variables
type(type_node_list),    intent(inout) :: node_list      !< list of nodes with grid information
type(type_element_list), intent(inout) :: element_list   !< list of elements with element information

! --- Local variables
integer, allocatable          :: n_parents(:)         ! for each node, want the number of parent elements
integer, allocatable          :: node_parents(:,:)    ! for each node, want to know the 4 parent elements
integer, allocatable          :: parent_elm_node(:,:) ! for each node, want to know the corresonding vertex for the 4 parent elements
integer                       :: i_node, i_elm, i_vertex, i_vertex1, i_vertex2, i
integer                       :: i_node_u, i_node_v
integer                       :: index
real*8                        :: size_u_min, size_v_min
real*8                        :: size_tmp
real*8                        :: point1(2), point2(2)
real*8                        :: distance1, distance2
real*8                        :: direction
real*8                        :: scale_uv, scale_ij, scale_wk
real*8                        :: distance_nodes, distance_node_point, distance_points

allocate( n_parents      (  node_list%n_nodes) )
allocate( node_parents   (4,node_list%n_nodes) )
allocate( parent_elm_node(4,node_list%n_nodes) )
n_parents       = 0
node_parents    = 0
parent_elm_node = 0

! --- Find parent elements
do i_node = 1, node_list%n_nodes

  n_parents(i_node) = 0
  do i_elm = 1, element_list%n_elements
  
    do i_vertex = 1, n_vertex_max
    
      if (element_list%element(i_elm)%vertex(i_vertex) .eq. i_node) then
        n_parents(i_node) = n_parents(i_node) + 1
        node_parents   (n_parents(i_node),i_node) = i_elm
        parent_elm_node(n_parents(i_node),i_node) = i_vertex
        exit
      endif
    
    enddo
    
    if ( (node_list%node(i_node)%boundary .eq. 2) .and. (n_parents(i_node) .eq. 2) ) exit
    if ( (node_list%node(i_node)%boundary .ne. 2) .and. (n_parents(i_node) .eq. 4) ) exit
  
  enddo
  
enddo


! --- For cubic elements, we scale by 1/3, for quintic elements, by 1/5
scale_uv = 0.2
! scale by 1/10 of the element side?
scale_ij = 0.001 !1.0
! scale by 1/10 of the element side?
scale_wk = 0.0
  

! --- Redefine all vectors at each node
do i_node = 1, node_list%n_nodes

  ! --- Note: all vectors are always of unit size
  ! --- it's the element size that changes


  ! -------------------
  ! --- Vectors u and v
  ! -------------------

  ! --- We must freeze the u/v vector sizes, which means we 
  ! --- must make sure the vector size is based on the smallest element
  
  ! --- Loop over each parent element and get minimal size
  size_u_min = 1.d15
  size_v_min = 1.d15
  do i=1,n_parents(i_node)
    i_elm = node_parents(i,i_node)
    if (parent_elm_node(i,i_node) .eq. 1) then
      i_node_u = element_list%element(i_elm)%vertex(2)
      i_node_v = element_list%element(i_elm)%vertex(4)
    elseif (parent_elm_node(i,i_node) .eq. 2) then
      i_node_u = element_list%element(i_elm)%vertex(1)
      i_node_v = element_list%element(i_elm)%vertex(3)
    elseif (parent_elm_node(i,i_node) .eq. 3) then
      i_node_u = element_list%element(i_elm)%vertex(4)
      i_node_v = element_list%element(i_elm)%vertex(2)
    elseif (parent_elm_node(i,i_node) .eq. 4) then
      i_node_u = element_list%element(i_elm)%vertex(3)
      i_node_v = element_list%element(i_elm)%vertex(1)
    endif
    distance1 = distance_nodes(node_list,i_node,i_node_u)
    distance2 = distance_nodes(node_list,i_node,i_node_v)
    if (distance1 .ne. 0.d0) size_u_min = min(size_u_min, distance1)
    if (distance2 .ne. 0.d0) size_v_min = min(size_v_min, distance2)
  enddo
  if (size_u_min .eq. 1.d15) size_u_min = 0.d0 
  if (size_v_min .eq. 1.d15) size_v_min = 0.d0 
  
  ! --- Set the direction u and v
  ! --- Positive direction of u is always from node-(1->2) === (4->3)
  ! --- Positive direction of v is always from node-(1->4) === (2->3)
  ! --- First, we set a temporary size of the vectors u,v to make sure they
  ! --- are smaller than the element side
  call normalise_vector(node_list,i_node,2)
  call normalise_vector(node_list,i_node,3)
  node_list%node(i_node)%X(1,2,1) = node_list%node(i_node)%X(1,2,1) * size_u_min * scale_uv
  node_list%node(i_node)%X(1,2,2) = node_list%node(i_node)%X(1,2,2) * size_u_min * scale_uv
  node_list%node(i_node)%X(1,3,1) = node_list%node(i_node)%X(1,3,1) * size_v_min * scale_uv
  node_list%node(i_node)%X(1,3,2) = node_list%node(i_node)%X(1,3,2) * size_v_min * scale_uv
  ! --- Then, we check the direction by locating the tip of the vector relative to the nodes
  do i=1,n_parents(i_node)
    i_elm = node_parents(i,i_node)
    if (parent_elm_node(i,i_node) .eq. 1) then
      i_node_u = element_list%element(i_elm)%vertex(2)
      i_node_v = element_list%element(i_elm)%vertex(4)
    elseif (parent_elm_node(i,i_node) .eq. 2) then
      i_node_u = element_list%element(i_elm)%vertex(1)
      i_node_v = element_list%element(i_elm)%vertex(3)
    elseif (parent_elm_node(i,i_node) .eq. 3) then
      i_node_u = element_list%element(i_elm)%vertex(4)
      i_node_v = element_list%element(i_elm)%vertex(2)
    elseif (parent_elm_node(i,i_node) .eq. 4) then
      i_node_u = element_list%element(i_elm)%vertex(3)
      i_node_v = element_list%element(i_elm)%vertex(1)
    endif
    ! --- compare size of side with distance from tip of vector u
    distance1 = distance_nodes(node_list,i_node,i_node_u)
    point1(1) = node_list%node(i_node)%X(1,1,1) + node_list%node(i_node)%X(1,2,1) ! tip of vector u
    point1(2) = node_list%node(i_node)%X(1,1,2) + node_list%node(i_node)%X(1,2,2) ! tip of vector u
    distance2 = distance_node_point(node_list,i_node_u,point1)
    ! --- reverse vector u if needed
    if (distance2 .lt. distance1) then
      if ( (parent_elm_node(i,i_node) .eq. 2) .or. (parent_elm_node(i,i_node) .eq. 3) ) then
        node_list%node(i_node)%X(1,2,1) = - node_list%node(i_node)%X(1,2,1)
        node_list%node(i_node)%X(1,2,2) = - node_list%node(i_node)%X(1,2,2)
      endif
    endif
    ! --- compare size of side with distance from tip of vector v
    distance1 = distance_nodes(node_list,i_node,i_node_v)
    point1(1) = node_list%node(i_node)%X(1,1,1) + node_list%node(i_node)%X(1,3,1) ! tip of vector v
    point1(2) = node_list%node(i_node)%X(1,1,2) + node_list%node(i_node)%X(1,3,2) ! tip of vector v
    distance2 = distance_node_point(node_list,i_node_v,point1)
    ! --- reverse vector v if needed
    if (distance2 .lt. distance1) then
      if ( (parent_elm_node(i,i_node) .eq. 3) .or. (parent_elm_node(i,i_node) .eq. 4) ) then
        node_list%node(i_node)%X(1,3,1) = - node_list%node(i_node)%X(1,3,1)
        node_list%node(i_node)%X(1,3,2) = - node_list%node(i_node)%X(1,3,2)
      endif
    endif
  enddo
  ! --- Finally, we set u,v back to unit vectors
  call normalise_vector(node_list,i_node,2)
  call normalise_vector(node_list,i_node,3)
  
  ! --- SIZE CONDITION: h_u and h_v must be the same on all parent elements (ie. breaks Bezier to become Hermite)
  ! --- Set the element size for nodes u and v
  do i=1,n_parents(i_node)
    i_elm = node_parents(i,i_node)
    i_vertex = parent_elm_node(i,i_node)
    if (i_vertex .eq. 1) then
      element_list%element(i_elm)%size(i_vertex,2) = + size_u_min * scale_uv
      element_list%element(i_elm)%size(i_vertex,3) = + size_v_min * scale_uv
    elseif (i_vertex .eq. 2) then
      element_list%element(i_elm)%size(i_vertex,2) = - size_u_min * scale_uv
      element_list%element(i_elm)%size(i_vertex,3) = + size_v_min * scale_uv
    elseif (i_vertex .eq. 3) then
      element_list%element(i_elm)%size(i_vertex,2) = - size_u_min * scale_uv
      element_list%element(i_elm)%size(i_vertex,3) = - size_v_min * scale_uv
    elseif (i_vertex .eq. 4) then
      element_list%element(i_elm)%size(i_vertex,2) = + size_u_min * scale_uv
      element_list%element(i_elm)%size(i_vertex,3) = - size_v_min * scale_uv
    endif
  enddo
  
  
  ! ------------
  ! --- Vector w
  ! ------------
  
  ! --- There is no rule for vector w, so we just set it as
  ! --- w = u+v, of length scale_wk*(u+v). But the condition comes on the element size
  
  ! --- Vector definition
  node_list%node(i_node)%X(1,4,1) = scale_wk * ( node_list%node(i_node)%X(1,2,1) + node_list%node(i_node)%X(1,3,1) )
  node_list%node(i_node)%X(1,4,2) = scale_wk * ( node_list%node(i_node)%X(1,2,2) + node_list%node(i_node)%X(1,3,2) )
  ! --- SIZE CONDITION: h_w = h_u*h_v
  do i = 1,n_parents(i_node)
    i_elm    = node_parents(i,i_node)
    i_vertex = parent_elm_node(i,i_node)
    element_list%element(i_elm)%size(i_vertex,4) = element_list%element(i_elm)%size(i_vertex,2) * element_list%element(i_elm)%size(i_vertex,3)
  enddo

enddo


! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!

! --- After definition of u,v,w, we are ready for i and j.
do i_node = 1, node_list%n_nodes
  
  ! ------------
  ! --- Vector i
  ! ------------
  
  ! --- Like w, we define i as i = u+v, of length scale_ij*(u+v)
  ! --- The real condition comes on the element size
  ! --- Which simply that the size should be exactly the same on all elements
  ! --- we set the size to be small, scale_ij*(size_u+size_v) of the first node available,
  ! --- with scale_i = 0.1 ?
  
  ! --- Vector definition
  node_list%node(i_node)%X(1,5,1) = 0.d0 ! node_list%node(i_node)%X(1,2,1) + node_list%node(i_node)%X(1,3,1)
  node_list%node(i_node)%X(1,5,2) = 0.d0 ! node_list%node(i_node)%X(1,2,2) + node_list%node(i_node)%X(1,3,2)
  ! --- SIZE CONDITION: h_i must be the same on all parent nodes
  i_elm = node_parents(1,i_node)
  size_tmp = 0.5 * ( abs(element_list%element(i_elm)%size(1,2)) + abs(element_list%element(i_elm)%size(1,3)) )

  ! --- METHOD: Assume curvature in the middle of edge is 0.5(curv_00+curv_05) 
  ! --- Find the next node
  i = n_parents(1)
  i_elm = node_parents(i,i_node)
  if (parent_elm_node(i,i_node) .eq. 1) then
    i_vertex1 = 1 ; i_vertex2 = 2
    i_node_u = element_list%element(i_elm)%vertex(2)
  elseif (parent_elm_node(i,i_node) .eq. 2) then
    i_vertex1 = 2 ; i_vertex2 = 1
    i_node_u = element_list%element(i_elm)%vertex(1)
  elseif (parent_elm_node(i,i_node) .eq. 3) then
    i_vertex1 = 3 ; i_vertex2 = 4
    i_node_u = element_list%element(i_elm)%vertex(4)
  elseif (parent_elm_node(i,i_node) .eq. 4) then
    i_vertex1 = 4 ; i_vertex2 = 3
    i_node_u = element_list%element(i_elm)%vertex(3)
  endif
  ! --- Formulation: h_i*\vec{i} = -0.25 * ( h_00^u*\vec{u}_00 + h_50^u*\vec{u}_50 )
  point1(1) = -0.25 * ( element_list%element(i_elm)%size(i_vertex1,2)*node_list%node(i_node  )%x(1,2,1) &
                       +element_list%element(i_elm)%size(i_vertex2,2)*node_list%node(i_node_u)%x(1,2,1) )
  point1(2) = -0.25 * ( element_list%element(i_elm)%size(i_vertex1,2)*node_list%node(i_node  )%x(1,2,2) &
                       +element_list%element(i_elm)%size(i_vertex2,2)*node_list%node(i_node_u)%x(1,2,2) )
  point2(1) = 0.d0 ; point2(2) = 0.d0
  distance1 = distance_points(point1,point2)
  node_list%node(i_node)%X(1,5,1) = point1(1)
  node_list%node(i_node)%X(1,5,2) = point1(2)
  call normalise_vector(node_list,i_node,5)
  size_tmp = distance1
  do i = 1,n_parents(i_node)
    i_elm    = node_parents(i,i_node)
    i_vertex = parent_elm_node(i,i_node)
    element_list%element(i_elm)%size(i_vertex,5) = size_tmp * scale_ij
  enddo
  
  
  ! ------------
  ! --- Vector j
  ! ------------
  
  ! --- METHOD: Assume curvature in the middle of edge is 0.5(curv_00+curv_05) 
  ! --- Find the next node
  i = n_parents(1)
  i_elm = node_parents(i,i_node)
  if (parent_elm_node(i,i_node) .eq. 1) then
    i_vertex1 = 1 ; i_vertex2 = 4
    i_node_v = element_list%element(i_elm)%vertex(4)
  elseif (parent_elm_node(i,i_node) .eq. 2) then
    i_vertex1 = 2 ; i_vertex2 = 3
    i_node_v = element_list%element(i_elm)%vertex(3)
  elseif (parent_elm_node(i,i_node) .eq. 3) then
    i_vertex1 = 3 ; i_vertex2 = 2
    i_node_v = element_list%element(i_elm)%vertex(2)
  elseif (parent_elm_node(i,i_node) .eq. 4) then
    i_vertex1 = 4 ; i_vertex2 = 1
    i_node_v = element_list%element(i_elm)%vertex(1)
  endif
  ! --- Formulation: h_j*\vec{j} = -0.25 * ( h_00^v*\vec{v}_00 + h_05^v*\vec{v}_05 )
  point1(1) = -0.25 * ( element_list%element(i_elm)%size(i_vertex1,3)*node_list%node(i_node  )%x(1,3,1) &
                       +element_list%element(i_elm)%size(i_vertex2,3)*node_list%node(i_node_v)%x(1,3,1) )
  point1(2) = -0.25 * ( element_list%element(i_elm)%size(i_vertex1,3)*node_list%node(i_node  )%x(1,3,2) &
                       +element_list%element(i_elm)%size(i_vertex2,3)*node_list%node(i_node_v)%x(1,3,2) )
  point2(1) = 0.d0 ; point2(2) = 0.d0
  distance1 = distance_points(point1,point2)
  node_list%node(i_node)%X(1,6,1) = element_list%element(i_elm)%size(i_vertex1,2) * node_list%node(i_node)%X(1,2,1) !point1(1)
  node_list%node(i_node)%X(1,6,2) = element_list%element(i_elm)%size(i_vertex1,2) * node_list%node(i_node)%X(1,2,2) !point1(2)
  call normalise_vector(node_list,i_node,6)
  ! --- SIZE CONDITION: h_j must be the same on all parent nodes
  size_tmp = distance1
  do i = 1,n_parents(i_node)
    i_elm    = node_parents(i,i_node)
    i_vertex = parent_elm_node(i,i_node)
    element_list%element(i_elm)%size(i_vertex,6) = size_tmp * scale_ij
  enddo
  
  
  ! ------------
  ! --- Vector m
  ! ------------
  
  ! --- m needs to be across from u, so we simply set m = v, of unit size
  ! --- and we set the size the same as that of the v vector
  
  ! --- Vector definition
  node_list%node(i_node)%X(1,7,1) = node_list%node(i_node)%X(1,3,1)
  node_list%node(i_node)%X(1,7,2) = node_list%node(i_node)%X(1,3,2)
  call normalise_vector(node_list,i_node,7)
  ! --- SIZE CONDITION: h_m must opposite on either side of the two parent nodes (like v)
  do i = 1,n_parents(i_node)
    i_elm    = node_parents(i,i_node)
    i_vertex = parent_elm_node(i,i_node)
    element_list%element(i_elm)%size(i_vertex,7) = element_list%element(i_elm)%size(i_vertex,3)
  enddo
  
  
  ! ------------
  ! --- Vector n
  ! ------------
  
  ! --- n needs to be across from v, so we simply set n = u, of unit size
  ! --- and we set the size the same as that of the u vector
  
  ! --- Vector definition
  node_list%node(i_node)%X(1,8,1) = node_list%node(i_node)%X(1,2,1)
  node_list%node(i_node)%X(1,8,2) = node_list%node(i_node)%X(1,2,2)
  call normalise_vector(node_list,i_node,8)
  ! --- SIZE CONDITION: h_n must opposite on either side of the two parent nodes (like u)
  do i = 1,n_parents(i_node)
    i_elm    = node_parents(i,i_node)
    i_vertex = parent_elm_node(i,i_node)
    element_list%element(i_elm)%size(i_vertex,8) = element_list%element(i_elm)%size(i_vertex,2)
  enddo
  
  
  ! ------------
  ! --- Vector k
  ! ------------
  
  ! --- Because it is very similar to vector w, we simply set k = w
  ! --- and we set the size the same as that of the u vector
  
  ! --- Vector definition
  node_list%node(i_node)%X(1,9,1) = node_list%node(i_node)%X(1,4,1)
  node_list%node(i_node)%X(1,9,2) = node_list%node(i_node)%X(1,4,2)
  ! --- SIZE CONDITION: h_k = h_i*h_j       ! OLD: (h_u*h_v)**2 = h_w**2
  do i = 1,n_parents(i_node)
    i_elm    = node_parents(i,i_node)
    i_vertex = parent_elm_node(i,i_node)
    element_list%element(i_elm)%size(i_vertex,9) = element_list%element(i_elm)%size(i_vertex,5) *  element_list%element(i_elm)%size(i_vertex,6)
  enddo

  ! --- IMPORTANT NOTE: We assume that "force_central_node" is used
  ! --- This means we identify axis nodes by checking the %index(1)=1 nodes
  ! --- On these noes, the poloidal vector sizes should be zero, ie. vectors v,w,j,n,k
  ! --- IMPORTANT NOTE: I'm not sure about n and k!
  if (node_list%node(i_node)%index(1) .eq. 1) then
    do i = 1,n_parents(i_node)
      i_elm    = node_parents(i,i_node)
      i_vertex = parent_elm_node(i,i_node)
      element_list%element(i_elm)%size(i_vertex,3) = 0.d0
      element_list%element(i_elm)%size(i_vertex,4) = 0.d0
      element_list%element(i_elm)%size(i_vertex,6) = 0.d0
      !element_list%element(i_elm)%size(i_vertex,8) = 0.d0
      !element_list%element(i_elm)%size(i_vertex,9) = 0.d0
    enddo
  endif


enddo





! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!

! --- Redefine node indexes in the matrix
index = 0
do i_node = 1, node_list%n_nodes
  
  ! --- Careful with force_axis_nodes
  if (node_list%node(i_node)%index(1) .eq. 1) then
    if (index .eq. 0) index = 1
    do i=2,n_degrees
      node_list%node(i_node)%index(i) = index + i-1
    enddo
    index = index + n_degrees-1
  else
    do i=1,n_degrees
      node_list%node(i_node)%index(i) = index + i
    enddo
    index = index + n_degrees
  endif
enddo


return

end subroutine transform_to_bi_quintic








pure real*8 function distance_nodes(node_list,i_node1,i_node2)
use data_structure
implicit none
type(type_node_list), intent(in) :: node_list
integer,              intent(in) :: i_node1,i_node2
distance_nodes = sqrt(  (node_list%node(i_node2)%X(1,1,1) - node_list%node(i_node1)%X(1,1,1))**2 &
                      + (node_list%node(i_node2)%X(1,1,2) - node_list%node(i_node1)%X(1,1,2))**2 )
end function distance_nodes

pure real*8 function distance_node_point(node_list,i_node,point)
use data_structure
implicit none
type(type_node_list), intent(in) :: node_list
integer,              intent(in) :: i_node
real*8,               intent(in) :: point(2)
distance_node_point = sqrt(  (point(1) - node_list%node(i_node)%X(1,1,1))**2 &
                           + (point(2) - node_list%node(i_node)%X(1,1,2))**2 )
end function distance_node_point

pure real*8 function distance_points(point1,point2)
use data_structure
implicit none
real*8, intent(in) :: point1(2),point2(2)
distance_points = sqrt(  (point2(1) - point1(1))**2 + (point2(2) - point1(2))**2 )
end function distance_points

pure subroutine normalise_vector(node_list,i_node,i_vector)
use data_structure
implicit none
type(type_node_list), intent(inout) :: node_list
integer,              intent(in)    :: i_node, i_vector
real*8                              :: vector_size
  vector_size = sqrt( node_list%node(i_node)%X(1,i_vector,1)**2.0 + node_list%node(i_node)%X(1,i_vector,2)**2.0 )
  if (vector_size .eq. 0.d0) return ! yes, because with fix_axis, we set v and w to zero!
  node_list%node(i_node)%X(1,i_vector,1) = node_list%node(i_node)%X(1,i_vector,1) / vector_size
  node_list%node(i_node)%X(1,i_vector,2) = node_list%node(i_node)%X(1,i_vector,2) / vector_size
end subroutine normalise_vector







