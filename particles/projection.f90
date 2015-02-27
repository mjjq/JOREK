module gauss

 integer, parameter :: n_gauss   = 4                  !< Number of Gaussian points
 integer, parameter :: n_gauss_2 = n_gauss * n_gauss  !< Square of n_gauss

 real*8,  parameter :: Xgauss(n_gauss) = (/ 0.0694318442029735d0, 0.3300094782075720d0,            &
   0.6699905217924280d0, 0.9305681557970265d0 /)      !< Positions of Gaussian points

 real*8,  parameter :: Wgauss(n_gauss) = (/ 0.173927422568727d0,  0.326072577431273d0,             &
   0.326072577431273d0,  0.173927422568727d0  /)      !< Weights of Gaussian points

end module gauss

program test_projection

use gauss

implicit none

real*8, allocatable :: part_R(:), R(:), f(:), f4(:,:), fractions(:), g(:)
real*8              :: H(2,2), H_s(2,2), H_ss(2,2), c(4), A(4,4), fg(4), d(4), width, s, R_width
real*8              :: W, W_s, W_ss, weight, error1, error2
integer             :: i, j, k, n_particles, np, info, n_width, nk, n_in

A(1,:) = (/ 52., 22., 18., 13. /)
A(2,:) = (/ 22., 12., 13.,  9. /)
A(3,:) = (/ 18., 13., 52., 22. /)
A(4,:) = (/ 13.,  9., 22., 12. /)
A = A/140.0

call dpotrf('U',4,A,4,info)

call random_seed()

n_particles = 1000
np = 5

allocate(part_R(n_particles))
allocate(R(np),f(np),f4(np,4),g(np),fractions(np-1))

n_width = 11
width   = sqrt(1./float(n_particles))

error1 = 0.
error2 = 0.

nk = 100

do k=1,nk

call random_number(part_R)

c = 0.
n_in = 0

do i=1, n_particles

  weight = 0.
  do j=1,n_width

    s = -1. + 2.*float(j-1)/float(n_width-1)
    call basisfunctions1(abs(s),W,W_s,W_ss)

!    R_width  = min(1.,max(part_R(i) + s * width,0.))
    R_width  = part_R(i) + s * width

    if (n_width .eq. 1) then
      R_width = part_R(i)
      W = 1.
    endif

    weight = weight + W

    if ((R_width .ge. 0.) .and. (R_width .le. 1.)) then

      call basisfunctions1(R_width,H,H_s,H_ss)

      c(1) = c(1) + W * H(1,1)
      c(2) = c(2) + W * H(1,2)
      c(3) = c(3) + W * H(2,1)
      c(4) = c(4) + W * H(2,2)

      n_in = n_in + 1
    endif

  enddo

enddo

c = c / float(n_in) / weight * float(n_width)

!write(*,'(A,4e14.6)') ' c : ',c

call dpotrs('U',4,1,A,4,c,4,info)

!write(*,'(A,4e14.6)') ' c : ',c


fractions = 0
do i=1,n_particles
  j = floor(part_R(i)*float(np-1)) + 1
  fractions(j) = fractions(j) + 1
enddo
fractions = float(np-1) * fractions / float(n_particles)

do i=1,np

  R(i) = float(i-1)/float(np-1)

  call basisfunctions1(R(i),H,H_s,H_ss)

  f(i) = c(1)*H(1,1) + c(2)*H(1,2) + c(3)*H(2,1) + c(4)*H(2,2)

enddo

error1 = error1 + abs((maxval(f)-minval(f))/2.)
error2 = error2 + abs((maxval(fractions)-minval(fractions))/2.)

enddo
write(*,'(4e12.4)') error1/float(nk), error2/float(nk),(maxval(f)-minval(f))/2.,(maxval(fractions)-minval(fractions))/2.

call begplt('plot.ps')
call nframe(11,11,1,0.0,1.0,0.0,2.0,'projection',10,'R',1,'density',7)
call lplot6(1,1,R+0.5/float(np-1),fractions,-np+1,'particles')
call lincol(1)
call lplot6(1,1,R,f,-np,'projection')
call lincol(0)
call finplt

end

subroutine basisfunctions1(s,H,H_s,H_ss)

implicit none

! --- Routine parameters
real*8, intent(in)  :: s          !< s-coordinate in the element
real*8, intent(out) :: H(2,2)     !< Basis functions
real*8, intent(out) :: H_s(2,2)   !< Basis functions derived with respect to s
real*8, intent(out) :: H_ss(2,2)  !< Basis functions derived two times with respect to s

!---------------------------------------------------------- vertex (1)
H(1,1)   =(-1.d0 + s)**2*(1.d0 + 2.d0*s)
H_s(1,1) =6.d0*(-1.d0 + s)*s
H_ss(1,1)=6.d0*(-1.d0 + 2.d0*s)

H(1,2)   =3.d0*(-1.d0 + s)**2*s
H_s(1,2) =3.d0*(-1.d0 + s)*(-1.d0 + 3.d0*s)
H_ss(1,2)=6.d0*(-2.d0 + 3.d0*s)

!---------------------------------------------------------- vertex (2)
H(2,1)   =-s**2*(-3.d0 + 2.d0*s)
H_s(2,1) =-6.d0*(-1.d0 + s)*s
H_ss(2,1)=-6.d0*(-1.d0 + 2.d0*s)

H(2,2)   =-3.d0*(-1.d0 + s)*s**2
H_s(2,2) =-3.d0*s*(-2.d0 + 3.d0*s)
H_ss(2,2)=-6.d0*(-1.d0 + 3.d0*s)

return
end
