subroutine initialise_basis()
!---------------------------------------------------------------
! calculates the basis functions at the Gaussian points
!---------------------------------------------------------------
use constants
use gauss
use basis_at_gaussian
use phys_module

implicit none

! --- local variables
integer :: i,k,l
real*8  :: s,t,phi
real*8 :: G(n_vertex_max,n_order+1,n_gauss,n_gauss)
real*8 :: G_s(n_vertex_max,n_order+1,n_gauss,n_gauss)
real*8 :: G_t(n_vertex_max,n_order+1,n_gauss,n_gauss)
real*8 :: G_st(n_vertex_max,n_order+1,n_gauss,n_gauss)

do k=1,n_gauss

 s = xgauss(k)
 
 call basisfunctions1(s,H1(1:2,1:4,k), H1_s(1:2,1:4,k), H1_ss(1:2,1:4,k)) ! the one-D basis functions

 do l=1,n_gauss
 
   t = xgauss(l)
   
   call basisfunctions2(s,t,H(1:4,1:4,k,l),   H_s(1:4,1:4,k,l), H_t(1:4,1:4,k,l), &
                            H_st(1:4,1:4,k,l),H_ss(1:4,1:4,k,l),H_tt(1:4,1:4,k,l) )

   call basisfunctions(s,t,G(1:4,1:4,k,l),   G_s(1:4,1:4,k,l), G_t(1:4,1:4,k,l), &
                           G_st(1:4,1:4,k,l))

!    write(*,'(2i3,20f8.4)') k,l,s,t,H(1:4,1:4,k,l)
!    write(*,'(2i3,20f8.4)') k,l,s,t,H_s(1:4,1:4,k,l)
!    write(*,'(2i3,20f8.4)') k,l,s,t,H_t(1:4,1:4,k,l)
!    write(*,'(2i3,20f8.4)') k,l,s,t,H_st(1:4,1:4,k,l)

 enddo
enddo

if (abs(maxval(H-G))       .gt. 1d-14) write(*,*) ' error in Basisfunctions H    : ',abs(maxval(H-G))
if (abs(maxval(H_s-G_s))   .gt. 1d-14) write(*,*) ' error in Basisfunctions H_s  : ',abs(maxval(H_s-G_s))
if (abs(maxval(H_t-G_t))   .gt. 1d-14) write(*,*) ' error in Basisfunctions H_t  : ',abs(maxval(H_t-G_t))
if (abs(maxval(H_st-G_st)) .gt. 1d-14) write(*,*) ' error in Basisfunctions H_st : ',abs(maxval(H_st-G_st))

do k=1,n_plane

  phi = 2.d0*PI*float(k-1)/float(n_plane) / float(n_period)

  HZ(1,k)   = 1.d0
  HZ_p(1,k) = 0.d0

  do i=1,(n_tor-1)/2

!    write(*,'(A,12i6)') ' ini basis : ',k,i,2*i,2*i+1,mode(2*i),mode(2*i+1)

    HZ(2*i,k)      =                           cos(mode(2*i)  *phi)
    HZ_p(2*i,k)    = - float(mode(2*i))      * sin(mode(2*i)  *phi)
    HZ_pp(2*i,k)   = - float(mode(2*i))**2   * cos(mode(2*i)  *phi)
    HZ(2*i+1,k)    =                           sin(mode(2*i+1)*phi)
    HZ_p(2*i+1,k)  = + float(mode(2*i+1))    * cos(mode(2*i+1)*phi)
    HZ_pp(2*i+1,k) = - float(mode(2*i+1))**2 * sin(mode(2*i+1)*phi)

  enddo

enddo

return
end subroutine initialise_basis

