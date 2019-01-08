      program eqdsk
!--------------------------------------------------------------------
! little program to construct an input file for jorek out of data
! in a eqdsk file
!                         Guido Huysmans,          date : 14-12-2010
!
! Some documentation can be found here: https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90
!--------------------------------------------------------------------
implicit none

real*8,allocatable :: psi(:),p(:),f(:),q(:),rlim(:),zlim(:),rbnd(:), zbnd(:)
real*8,allocatable :: apsi(:),bpsi(:),cpsi(:),dpsi(:)
real*8,allocatable :: ap(:),bp(:),cp(:),dp(:)
real*8,allocatable :: af(:),bf(:),cf(:),df(:)
real*8,allocatable :: radius(:),theta(:),rad(:)
real*8,allocatable :: ar(:),br(:),cr(:),dr(:)
real*8,allocatable :: dpr(:),df2(:),dg(:),work(:),psirz(:,:)
real*8,allocatable :: xx(:),yy(:),zc(:), r_bnd(:), z_bnd(:), psi_bnd(:)
real*8,allocatable :: df2_ext(:),rho_ext(:),T_ext(:),psi_ext(:),p_ext(:)
real*8,allocatable :: tx(:),ty(:),c(:,:),wrk(:)
integer,allocatable :: iwrk(:)
real*8             :: angle, ellip, tria_u, tria_l, quad_u, quad_l, r0, z0, a0, PI
real*8             :: dummy(3), xdim,zdim,rzero,rgrid1,zmid,rmaxis,zmaxis,ssimag,ssibry,bcentr
real*8             :: xip,xdum1,xdum2,xdum3,xdum4,xdum5
real*8             :: psi_sep, sig_sep, tanh1, zmu0, zn0, zmd
real*8             :: xb ,xe, yb, ye, smth, fp, fout
integer            :: mx,my,kx,ky,nxest,nyest,lwrk,kwrk,ier,iopt,nx,ny, i1, j1
integer            :: nr, nz, n_psi, nbbs, limitr, i,j, nc, n_tht, n_sol, n_ext
character          :: AA*52, tokamak_name*50

!----------------------------- read eqdsk file -----------

write(*,*) ' EQDSK to JOREK2 '

tokamak_name = 'JET' ! 'ITER', 'DIII-D'

write(*,*) 'Tokamak = ', tokamak_name

read(5,'(A52,2i4)') AA,nr,nz

write(*,*) AA
write(*,'(A,2i5)') ' nr, nz : ',nr,nz

read(5,'(5e16.9)') xdim,zdim,rzero,rgrid1,zmid
read(5,'(5e16.9)') rmaxis,zmaxis,ssimag,ssibry,bcentr
read(5,'(5e16.9)') xip,ssimag,xdum1,rmaxis,xdum2
read(5,'(5e16.9)') zmaxis,xdum3,ssibry,xdum4,xdum5

write(*,'(A,2f10.5,A)') ' xdim,  zdim : ',xdim,zdim,' m'
write(*,'(A,2f10.5,A)') ' rzero, zmid : ',rzero,zmid, ' m'
write(*,'(A,f10.5,A)')  ' xip         : ',xip/1e6,' MA'
write(*,'(A,f10.5,A)')  ' zmaxis      : ',zmaxis,' m'
write(*,'(A,f10.5,A)')  ' Bvac        : ',bcentr,' T'

      
write(*,*) ' reading profiles'
  
n_psi=nr
allocate(f(n_psi),p(n_psi),df2(n_psi),dpr(n_psi),psirz(nr,nz),q(n_psi))

read(5,'(5e16.9)') (f(i),i=1,n_psi)
read(5,'(5e16.9)') (p(i),i=1,n_psi)
read(5,'(5e16.9)') (df2(i),i=1,n_psi)
read(5,'(5e16.9)') (dpr(i),i=1,n_psi)
read(5,'(5e16.9)') ((psirz(i,j),i=1,nr),j=1,nz)

read(5,'(5e16.9)') (q(i),i=1,n_psi)

write(*,*) ' reading limiter'

read(5,*)  nbbs,limitr
allocate(rbnd(nbbs),zbnd(nbbs))
read(5,'(5e16.9)') (rbnd(i),zbnd(i),i=1,nbbs)
allocate(rlim(limitr),zlim(limitr))
read(5,'(5e16.9)') (rlim(i),zlim(i),i=1,limitr)

write(*,*) ' done reading'

allocate(xx(nr),yy(nz),psi(n_psi))
do i=1,n_psi
  psi(i) = real(i-1)/real(n_psi-1)
enddo     
do i=1,nr
  xx(i) = rgrid1 + xdim*real(i-1)/real(nr-1)
enddo     
do i=1,nz
  yy(i) = zmid + zdim*(real(i-1)/real(nz-1)-0.5)
enddo     

if (tokamak_name == 'ITER') then

  !--------------------close fit to ITER wall
  ellip  = 2.0
  tria_u = 0.55
  tria_l = 0.65
  quad_u = -0.1
  quad_l = 0.15
  n_tht   = 257
  r0     = 6.2
  z0     = 0.1
  a0     = 2.25 

  !-------------------- contour outside ITER wall
  ellip  = 2.1
  tria_u = 0.58
  tria_l = 0.65
  quad_u = -0.12
  quad_l = -0.
  n_tht   = 257
  r0     = 6.2
  z0     = -0.05
  a0     = 2.34 

else if (tokamak_name == 'JET') then
  
  !-------------------- contour outside JET wall
  ! blue contour in https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90
  ellip  = 1.85
  tria_u = 0.4
  tria_l = 0.4
  quad_u = -0.2
  quad_l = -0.2
  n_tht   = 257
  r0     = 2.9
  z0     = 0.1
  a0     = 1.08

  !-------------------- contour to avoid too long divertor legs
  ! red contour in https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90

  ellip  = 1.7
  tria_u = 0.4
  tria_l = 0.4
  quad_u = -0.4
  quad_l = -0.2
  n_tht   = 257
  r0     = 2.85
  z0     = 0.15
  a0     = 1.1

else if (tokamak_name == 'DIII-D') then
  
  !-------------------- contour outside DIII-D wall
  ellip  = 1.85
  tria_u = 0.4
  tria_l = 0.4
  quad_u = -0.2
  quad_l = -0.2
  n_tht   = 257
  r0     = 1.7
  z0     = 0.
  a0     = 0.7

else

  write(*,*) 'Tokamak name not or wrongly specified, stopping'
  stop

end if  
  
PI = 2.d0 * asin(1.d0)

allocate(r_bnd(n_tht),z_bnd(n_tht),psi_bnd(n_tht))

!--------------------------------- interpolate flux using Dierckx spline routine
iopt= 0 
mx = nr
my = nz
xb = xx(1)
xe = xx(nr)
yb = yy(1)
ye = yy(nz)
kx = 3
ky = 3
!beware, the smoothing parameter smth can strongly affect the interpolation routine and it should always be checked that the resulting interpolation along the JOREK boundary is smooth enough
!refer to the "Hard-coded parameters" section of the wikipage https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90
!smth = 1.d-6 !this value can give a non-smooth interpolation close to the X-point, resulting in artifacts in the flux-aligned mesh computed by JOREK
smth = 3.d-3 ! this is usually the best trade-off for good results both at and away from the Xpoint
nxest = nr-5
nyest = nz-5
lwrk  = 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+my*(ky+1)+my+nxest
kwrk  = 3+mx+my+nxest+nyest

allocate(tx(nxest),ty(nyest),c(nxest,nyest),wrk(lwrk),iwrk(kwrk))

call regrid(iopt,mx,xx,my,yy,transpose(psirz),xb,xe,yb,ye,kx,ky,smth,nxest,nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)

write(*,*) ' Dierckx ier   : ',ier
write(*,*) ' Dierckx fp    : ',fp
write(*,*) ' Dierckx nx,ny : ',nx,ny

lwrk = mx*(kx+1)+my*(ky+1)
kwrk = mx+my
deallocate(wrk,iwrk)
allocate(wrk(lwrk),iwrk(kwrk))

!do i=1,nr
!  call bispev(tx,nx,ty,ny,c,kx,ky,xx(i),1,yy(nz/2),1,fout,wrk,lwrk,iwrk,kwrk,ier)
!  write(*,'(4e16.8,i3)') xx(i),yy(nz/2),fout,psirz(i,nz/2),ier
!enddo

do i=1,n_tht/2
  angle = 2.d0 * PI * float(i-1)/float(n_tht-1)
  r_bnd(i) = r0 + a0 * cos(angle + tria_u*sin(angle) + quad_u*sin(2.d0*angle))
  z_bnd(i) = z0 + a0 * ellip * sin(angle)
  call bispev(tx,nx,ty,ny,c,kx,ky,r_bnd(i),1,z_bnd(i),1,psi_bnd(i),wrk,lwrk,iwrk,kwrk,ier)
enddo
do i=n_tht/2+1,n_tht
  angle = 2.d0 * PI * float(i-1)/float(n_tht-1)
  r_bnd(i) = r0 + a0 * cos(angle + tria_l*sin(angle) + quad_l*sin(2.d0*angle))
  z_bnd(i) = z0 + a0 * ellip * sin(angle)
  call bispev(tx,nx,ty,ny,c,kx,ky,r_bnd(i),1,z_bnd(i),1,psi_bnd(i),wrk,lwrk,iwrk,kwrk,ier)
enddo


write(*,*) ' plotting results'  
nc = 51
allocate(zc(nc))
call begplt('eqdsk.ps')
call lblbot('eqdsk data',10)

call cplot(22,1,0,xx,yy,nr,nz,1,1,psirz,n_psi,zc,-nc,'fluxcontours',12,'R [m]',5,'Z [m]',5)
call lincol(3)
call lplot6(2,1,rlim,zlim,-limitr,'limiter')
call lincol(1)
call lplot6(2,1,rbnd,zbnd,-nbbs,'boundary')
call lincol(2)
call lplot6(2,1,r_bnd,z_bnd,-n_tht,'JOREK boundary')
call lincol(0)
call lplot6(3,2,psi,p,n_psi,'pressure')
call lplot6(3,3,psi,q,n_psi,'q')

call lplot6(2,2,psi,df2,n_psi,'df2')
call lplot6(3,2,psi,p,n_psi,'pressure')
call lplot6(2,3,psi,f,n_psi,'f')
call lplot6(3,3,psi,q,n_psi,'q')


!---------------------------- write JOREK input files
n_sol = (n_psi-1)/2
n_ext = n_psi + n_sol

write(*,*) ' n_psi, n_sol, n_ext : ',n_psi, n_sol, n_ext

allocate(df2_ext(n_ext),rho_ext(n_ext),T_ext(n_ext),psi_ext(n_ext),p_ext(n_ext))

df2_ext(1:n_psi) = df2(1:n_psi)
rho_ext(1:n_psi) = 1.d0
T_ext(1:n_psi)   = p(1:n_psi)

df2_ext(n_psi-1:n_ext) = df2_ext(n_psi-2)
rho_ext(n_psi-1:n_ext) = rho_ext(n_psi-2)
T_ext(n_psi-1:n_ext)   = T_ext(n_psi-2)

psi_sep = 1.d0
sig_sep = 0.02

psi_ext(1:n_psi) = psi(1:n_psi)
do i=n_psi+1,n_ext
  psi_ext(i) = 1.d0 + 0.5 * float(i-n_psi)/float(n_sol)
enddo

zmu0 = 4.d-7 * PI

do i=1,n_ext
  tanh1 = tanh((psi_ext(i) - psi_sep)/sig_sep)
  df2_ext(i) = df2_ext(i) * (0.5d0 - 0.5d0*tanh1)
  rho_ext(i) = rho_ext(i) * (0.5d0 - 0.5d0*tanh1)
  T_ext(i)   = T_ext(i)   * (0.5d0 - 0.5d0*tanh1) * zmu0 
  p_ext(i)   = rho_ext(i) * T_ext(i)
enddo

call lplot6(2,2,psi_ext,df2_ext,n_ext,'df2')
call lplot6(3,2,psi_ext,p_ext,n_ext,'pressure')
call lplot6(2,3,psi_ext,rho_ext,n_ext,'density')
call lplot6(3,3,psi_ext,T_ext,n_ext,'T')

call lincol(1)
call lplot6(2,2,psi,df2,-n_psi,'df2')
call lplot6(3,2,psi,p,-n_psi,'pressure')

open(21,file='jorek_ffprime')
do i=1,n_ext
  write(21,*) psi_ext(i),-df2_ext(i)  ! Minus sign because ff' in JOREK is opposite to the usual ff' for historical reasons
enddo
close(21)
open(21,file='jorek_density')
do i=1,n_ext
  write(21,*) psi_ext(i),rho_ext(i)
enddo
close(21)
open(21,file='jorek_temperature')
do i=1,n_ext
  write(21,*) psi_ext(i),T_ext(i)
enddo
close(21)

open(21,file='jorek_namelist')

write(21,*)             '***************************************'
write(21,'(A)')        '*  namelist produced by eqdsk2jorek   *'
write(21,*)             '***************************************'
write(21,*) AA
write(21,'(A,f8.3,A)') ' magnetic field : ',Bcentr,' T'
write(21,'(A,f8.3,A)') ' current        : ',xip/1d6,' MA'
write(21,*)             '***************************************'
write(21,*)
write(21,*) ' &in1'

write(21,*) ' restart = .f.'
write(21,*) ' regrid  = .f.'
write(21,*) ' tstep   = 5.' 
write(21,*) ' nstep   = 0' 

write(21,*) ' nout = 1'

write(21,*)
write(21,*) ' !_____________________________________boundary definition'
write(21,*) ' mf = 0'
write(21,*) ' n_boundary = ',n_tht
do j=1,n_tht
  write(21,'(A,i3,A,e16.8,A,i3,A,e16.8,A,i3,A,e16.8,A)'), &
            '  R_boundary(',j,') =',r_bnd(j), &
            ', Z_boundary(',j,') =',z_bnd(j), &
            ', psi_boundary(',j,') =',psi_bnd(j),','
enddo

write(21,*) ' ellip  = ',ellip
write(21,*) ' tria_u = ',tria_u
write(21,*) ' tria_l = ',tria_l
write(21,*) ' quad_u = ',quad_u
write(21,*) ' quad_l = ',quad_l

write(21,*) ' xampl  = +0.'
write(21,*) ' xpoint = .t.'

write(21,*) ' freeboundary = .f.'
write(21,*) ' resistive_wall = .f.'

write(21,*) ' R_geo = ',r0
write(21,*) ' Z_geo = ',z0
if (tokamak_name=='JET') then
  write(21,*) ' F0    = ',2.96*bcentr ! By convention, the vacuum toroidal field is given at 2.96m in JET eqdsk files.
                                      ! Note that in principle we should put a minus sign due to opposite conventions
                                      ! for the toroidal angle in JOREK and EFIT, but presently (08/01/19) we reverse
                                      ! the field for reasons explained here: https://www.jorek.eu/wiki/doku.php?id=jet									  
else
  write(21,*) ' F0    = ',r0*bcentr
end if
write(21,*) ' amin  = 1.d0 ! scale factor for plasma size only'

write(21,*)
write(21,*) ' !_____________________________________grid parameters'

write(21,*) ' n_R      = 0'
write(21,*) ' n_Z      = 0'
write(21,*) ' n_radial = 41'
write(21,*) ' n_pol    = 64' 

write(21,*) ' n_flux   = 0'
write(21,*) ' n_tht    = 64'
write(21,*) ' n_open   = 15'
write(21,*) ' n_leg    = 15'
write(21,*) ' n_private = 9'
write(21,*) ' dPSI_open    = 0.04'
write(21,*) ' dPSI_private = 0.02'

write(21,*)
write(21,*) ' !_____________________________________physics parameters'

write(21,*) ' eta   = 1.d-7'
write(21,*) ' visco = 1.d-6'
write(21,*) ' visco_par = 1.d-5'

write(21,*) ' eta_num = 1.d-12'

write(21,*) ' rho_file     = "jorek_density"'
write(21,*) ' T_file       = "jorek_temperature"'
write(21,*) ' ffprime_file = "jorek_ffprime"'

write(21,*) ' D_par     = 0.d0'
write(21,*) ' D_perp(1) = 1.d-5'
write(21,*) ' D_perp(2) = 0.85d0'
write(21,*) ' D_perp(3) = 0.d0'
write(21,*) ' D_perp(4) = 0.01d0'
write(21,*) ' D_perp(5) = 0.92d0'

write(21,*) ' ZK_par     = 1.d0'
write(21,*) ' ZK_perp(1) = 1.d-5'
write(21,*) ' ZK_perp(2) = 0.85d0'
write(21,*) ' ZK_perp(3) = 0.d0'
write(21,*) ' ZK_perp(4) = 0.01d0'
write(21,*) ' ZK_perp(5) = 0.92d0'

write(21,*) ' heatsource     = 1.d-7'
write(21,*) ' particlesource = 5.d-6'

write(21,*) ' &end'

close(21)

call finplt

end
 
 
