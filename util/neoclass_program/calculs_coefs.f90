!attention : tau=Ti/Te défini par tau=1 ligne 306, ft ligne 221

MODULE prec_const
!
!   Precision for real and complex
!
  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(10)
  INTEGER, PARAMETER :: CKIND = RKIND
  INTEGER, PARAMETER :: RKIND4 = SELECTED_REAL_KIND(10)

!	real(rkind), parameter 		:: RMaj=3.
!	real(rkind), parameter 		:: eps=3.
!	real(rkind), parameter		:: BM=2.5
!
!   Some useful constants
!
  REAL(RKIND), PARAMETER :: PI=3.141592653589793238462643383279502884197_rkind
  REAL(RKIND), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_rkind
  REAL(RKIND), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_rkind
  REAL(RKIND), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_rkind
  !                                               
  ! Physical constants     
  !                                                                                            
  real(rkind),parameter :: kb = 1.3806504e-23_rkind
  real(rkind),parameter :: mproton = 1.672623e-27_rkind
  real(rkind),parameter :: charge =  1.60217653e-19_rkind
  real(rkind),parameter :: mu0 = 4.0_rkind*PI*1.e-7_rkind
  real(rkind),parameter :: melectron=9.1096e-31_rkind
! beware: OK with D plasma, not with D-T:
  real(rkind),parameter :: massnumber   = 2._rkind
  real(rkind),parameter :: eps0=8.85e-12_rkind
  real(rkind),parameter :: chargenumber   = 1._rkind

real(rkind),parameter :: mion=massnumber*mproton
real(rkind),parameter :: mratio=melectron/mion

END MODULE prec_const



MODULE neo_subroutines

  implicit none
  contains  

       subroutine neocoeffs(Npoints, eps, Rmaj, Te, Ti, Ne, Ni, qsf, rho, muneoe, muneoi)
!-------------------------------------------------------------------
!       compute neoclassical coefficients

      use prec_const  
    
      implicit none

      integer, intent(in)			  	      :: Npoints
      real(rkind),intent(in)			  	      :: eps
      real(rkind),intent(in)			  	      :: Rmaj
      real(rkind), dimension(:), intent(in)                   :: Te
      real(rkind), dimension(:), intent(in)                   :: Ti
      real(rkind), dimension(:), intent(in)                   :: qsf
      real(rkind), dimension(:), intent(in)                   :: Ne
      real(rkind), dimension(:), intent(in)                   :: Ni
      real(rkind), dimension(:), intent(in)                   :: rho
      real(rkind), dimension(:,:), intent(out)                :: muneoe
      real(rkind), dimension(:,:), intent(out)                :: muneoi

      integer                                                 :: i
      real(rkind), dimension(Npoints)                         :: xie
      real(rkind), dimension(Npoints)                         :: xei
      real(rkind), dimension(Npoints)                         :: epsK
      real(rkind), dimension(Npoints)                         :: ft
      real(rkind), dimension(Npoints)                         :: nuste
      real(rkind), dimension(Npoints)                         :: nusti
      real(rkind), dimension(Npoints)                         :: LneeNRL
      real(rkind), dimension(Npoints)                         :: LneiNRL
      real(rkind), dimension(Npoints)                         :: LniiNRL
      real(rkind), dimension(Npoints)                         :: Tauee
      real(rkind), dimension(Npoints)                         :: Tauei
      real(rkind), dimension(Npoints)                         :: Tauie
      real(rkind), dimension(Npoints)                         :: Tauii
      real(rkind)                                             :: M11aa
      real(rkind)                         		      :: M12aa
      real(rkind)                                             :: M21aa
      real(rkind)                                             :: M22aa
      real(rkind)                                             :: N11aa
      real(rkind)                                             :: N12aa
      real(rkind)                                             :: N21aa
      real(rkind)                                             :: N22aa
      real(rkind), dimension(Npoints)                         :: M11ei
      real(rkind), dimension(Npoints)                         :: M12ei
      real(rkind), dimension(Npoints)                         :: M21ei
      real(rkind), dimension(Npoints)                         :: M22ei
      real(rkind), dimension(Npoints)                         :: N11ei
      real(rkind), dimension(Npoints)                         :: N12ei
      real(rkind), dimension(Npoints)                         :: N21ei
      real(rkind), dimension(Npoints)                         :: N22ei
      real(rkind), dimension(Npoints)                         :: M11ie
      real(rkind), dimension(Npoints)                         :: M12ie
      real(rkind), dimension(Npoints)                         :: M21ie
      real(rkind), dimension(Npoints)                         :: M22ie
      real(rkind), dimension(Npoints)                         :: N11ie
      real(rkind), dimension(Npoints)                         :: N12ie
      real(rkind), dimension(Npoints)                         :: N21ie
      real(rkind), dimension(Npoints)                         :: N22ie
      real(rkind), dimension(Npoints)                         :: L11ee
      real(rkind), dimension(Npoints)                         :: L12ee
      real(rkind), dimension(Npoints)                         :: L21ee
      real(rkind), dimension(Npoints)                         :: L22ee
      real(rkind), dimension(Npoints)                         :: L11ei
      real(rkind), dimension(Npoints)                         :: L12ei
      real(rkind), dimension(Npoints)                         :: L21ei
      real(rkind), dimension(Npoints)                         :: L22ei
      real(rkind), dimension(Npoints)                         :: L11ie
      real(rkind), dimension(Npoints)                         :: L12ie
      real(rkind), dimension(Npoints)                         :: L21ie
      real(rkind), dimension(Npoints)                         :: L22ie
      real(rkind), dimension(Npoints)                         :: L11ii
      real(rkind), dimension(Npoints)                         :: L12ii
      real(rkind), dimension(Npoints)                         :: L21ii
      real(rkind), dimension(Npoints)                         :: L22ii





!
!     Coulomb logarithm (NRL)
!
      do i = 1, Npoints

         LneeNRL(i)=24.-0.5*log((Ne(i)+1e-30)/1e6)+ log(Te(i))
 !        LneeNRL(i)=23.5-0.5*log((Ne(i)+1e-30)/1e6)+ &
 !                   1.25*log(Te(i))-(1.e-5+((log(Te(i))-2.0)**2)/16.0)**(0.5)              

	if (Ti(i)*mratio<Te(i).and.Te(i)<10.*chargenumber**2) then
            LneiNRL(i)=23.-log(chargenumber*(Ne(i)/1e6+1e-30)**(0.5)) &
                      +1.5*log(Te(i))
         elseif (Ti(i)*mratio<10.*chargenumber**2..and.Te(i)>10.*chargenumber**2) then
            LneiNRL(i)=24.-0.5*log(Ne(i)/1e6) &
                      +log(Te(i))
         else
            LneiNRL(i)=30.-log(chargenumber**2/massnumber*(Ni(i)/1e6)**(0.5)) &
                      +1.5*log(Ti(i))
         endif
         LniiNRL(i)=23.-log(chargenumber**2*sqrt(2._rkind*Ni(i)/1e6*chargenumber**2)) &
                   +1.5*log(Ti(i))
      enddo


!
!     Collision time (Braginskii)
!
      Tauee(:)=(3.*eps0**2*sqrt(melectron)*(TWOPI*charge*Te(:))**(1.5))/ &
               (charge**4*(Ne(:)+1e-30)*LneeNRL(:))

      Tauei(:)=(3.*eps0**2*sqrt(melectron)*(TWOPI*charge*Te(:))**(1.5))/ &
               (charge**4*chargenumber**2*(Ni(:)+1e-30)*LneiNRL(:))

      Tauie(:)=(3.*eps0**2*sqrt(mion)*(TWOPI*charge*Ti(:))**(1.5))/ &
               (charge**4*chargenumber**2*(Ne(:)+1e-30)*LneiNRL(:))

      Tauii(:)=(3.*eps0**2*sqrt(mion)*(TWOPI*charge*Ti(:))**(1.5))/ &
               (charge**4*chargenumber**4*(Ni(:)+1e-30)*LniiNRL(:))

!
!     Friction forces
!
      xei(:)=(mratio*Ti(:)/Te(:))**(0.5)
      xie(:)=1./xei(:)

!
!     Friction matrices
!
      M11aa=-2./2.**(1.5)
!     M12aa=1.5*2./2.**(2.5)
      M12aa=-1.5*2./2.**(2.5)
      M21aa=M12aa
      M22aa=-(13./4.+4.+15./2.)/2.**(2.5)
      N11aa=-M11aa
      N12aa=-M12aa
      N21aa=-M21aa
      N22aa=27./4./2.**(2.5)
!
      M11ei=-(1.+mratio)/(1.+xei**2)**(1.5)
!      M12ei=1.5*(1.+mratio)/(1.+xei**2)**(2.5)
      M12ei=-1.5*(1.+mratio)/(1.+xei**2)**(2.5)
      M21ei=M12ei
      M22ei=-(13./4.+4.*xei**2+15./2.*xei**4)/(1.+xei**2)**(2.5)
      N11ei=-M11ei
      N12ei=-xei**2*M12ei
      N21ei=-M21ei
      N22ei=27./4.*sqrt(Te/Ti)*xei**2/(1.+xei**2)**(2.5)
! sqrt(Te/Ti) instead of (Te/Ti) for the symmetry lijab=ljiba (cf.Houlberg 1997)
      M11ie=-(1.+1./mratio)/(1.+xie**2)**(1.5)
      M12ie=-1.5*(1.+1./mratio)/(1.+xie**2)**(2.5)
      M21ie=M12ie
      M22ie=-(13./4.+4.*xie**2+15./2.*xie**4)/(1.+xie**2)**(2.5)
      N11ie=-M11ie
      N12ie=-xie**2*M12ie
      N21ie=-M21ie
      N22ie=27./4.*sqrt(Ti/Te)*xie**2/(1.+xie**2)**(2.5)

!
!     Friction coefficients
!	ajouter les termes Meimp/Taueimp si ajout d'impurité
      L11ee=M11aa/Tauee+M11ei/Tauei+N11aa/Tauee
      L12ee=M12aa/Tauee+M12ei/Tauei+N12aa/Tauee
      L21ee=M21aa/Tauee+M21ei/Tauei+N21aa/Tauee
      L22ee=M22aa/Tauee+M22ei/Tauei+N22aa/Tauee
      L11ii=M11aa/Tauii+M11ie/Tauie+N11aa/Tauii
      L12ii=M12aa/Tauii+M12ie/Tauie+N12aa/Tauii
      L21ii=M21aa/Tauii+M21ie/Tauie+N21aa/Tauii
      L22ii=M22aa/Tauii+M22ie/Tauie+N22aa/Tauii
      L11ei=N11ei/Tauei
      L12ei=N12ei/Tauei
      L21ei=N21ei/Tauei
      L22ei=N22ei/Tauei
      L11ie=N11ie/Tauie
      L12ie=N12ie/Tauie
      L21ie=N21ie/Tauie
      L22ie=N22ie/Tauie

!
!     Viscosity coefficients mu
!
      do i = 1, Npoints
         epsK(i) = (rho(i))/eps
         ft(i) = 1._rkind - (1._rkind - epsK(i))**2/ &
              (sqrt(1._rkind - epsK(i)**2)*(1+1.46_rkind*sqrt(epsK(i))))
!  	 ft(i)=2/PI*sqrt(2*(rho(i))/eps)

         nuste(i)= RMaj*qsf(i)/(sqrt(2._rkind*charge*Te(i)/melectron)*Tauee(i)*(epsK(i))**(1.5))
         nusti(i)= RMaj*qsf(i)/(sqrt(2._rkind*charge*Ti(i)/mion)*Tauii(i)*(epsK(i))**(1.5))
      enddo

      call nui_visc(Npoints, Tauee,Tauii,1,nuste,nusti,epsK,ft,muneoe)
      call nui_visc(Npoints, Tauee,Tauii,2,nuste,nusti,epsK,ft,muneoi)

      end subroutine




      subroutine nui_visc(Npoints, Tauee,Tauii,spec,nuste,nusti,epsK,ft,nuvi)
!-------------------------------------------------------------------
!       compute viscosity coefficients from Kessel94
!
      use prec_const


      implicit none
!
      integer 	     				  :: Npoints
      integer, parameter                  	  :: nb_xint = 201
      integer, parameter                      	  :: nspec = 2
      integer                           	  :: i
      integer                            	  :: j
      integer                            	  :: spec
      real(rkind), dimension(nspec)      	  :: Z_spec
      real(rkind), dimension(:)                   :: Tauee
      real(rkind), dimension(:)        	          :: Tauii
      real(rkind), dimension(nspec) 	          :: rtemp
      real(rkind), dimension(nspec)       	  :: nustar
      real(rkind), dimension(nspec)	          :: vt
      real(rkind), dimension(Npoints)         	  :: fc
      real(rkind), dimension(:)         	  :: ft
      real(rkind)                         	  :: dx
      real(rkind), dimension(:)         	  :: epsK
      real(rkind), dimension(nb_xint)     	  :: xa
      real(rkind)                         	  :: xb
      real(rkind), dimension(Npoints)         	  :: nuaDv
      real(rkind), dimension(Npoints)         	  :: nuasv
      real(rkind), dimension(Npoints)        	  :: nuapv
      real(rkind), dimension(Npoints)        	  :: nuaEv
      real(rkind), dimension(Npoints)       	  :: nuatv
      real(rkind), dimension(Npoints)       	  :: nuatotv
      real(rkind), dimension(Npoints)        	  :: F2
      real(rkind), dimension(Npoints)        	  :: F3
      real(rkind), dimension(nb_xint)    	  :: coef
      real(rkind), dimension(Npoints)        	  :: int11
      real(rkind), dimension(Npoints)         	  :: int12
      real(rkind), dimension(Npoints)         	  :: int22
      real(rkind)                         	  :: fact1
      real(rkind), dimension(Npoints)         	  :: fact2
      real(rkind), dimension(Npoints)         	  :: fact3
      real(rkind), dimension(Npoints)         	  :: Ki11
      real(rkind), dimension(Npoints)         	  :: Ki12
      real(rkind), dimension(Npoints)         	  :: Ki22
      real(rkind)                         	  :: derf
!      real(rkind)                         	  :: chand
      real(rkind), dimension(:,:)       	  :: nuvi
      real(rkind), dimension(:)         	  :: nuste
      real(rkind), dimension(:)         	  :: nusti
      real(rkind), dimension(Npoints)       	  :: aux1
      real(rkind), dimension(Npoints)        	  :: aux2
      real(rkind)		        	  :: tau=1.


!
      Z_spec(1) = -1._rkind
      Z_spec(2) = chargenumber
      rtemp(1) = 1._rkind
      rtemp(2) = tau
      vt(1) = 1._rkind
      vt(2) = sqrt(tau*mratio)
!
      do i = 1, nb_xint
         xa(i) = 0.0001 + real(i-1)*(4.0-0.0001)/real(nb_xint-1)
      enddo
      dx = xa(2)-xa(1)
      fc = 1._rkind - ft
      if (spec==1) then
         call mult(nuste,Tauee,aux1,Npoints)
! mult is defined at the end by : aux1 = nuste*Tauee
      else
         call mult(nusti,Tauii,aux1,Npoints)
      end if
      if (spec==1) then
         fact2 = 3._rkind*sqrt(PI)/4._rkind/Tauee
      else
         fact2 = 3._rkind*sqrt(PI)/4._rkind/Tauii
      end if
!
      coef(1) = 1._rkind/3._rkind*dx*exp(-xa(1)**2)
      do j = 2, nb_xint-1, 2
        coef(j) = 4._rkind/3._rkind*dx*exp(-xa(j)**2)
      enddo
      do j = 3, nb_xint-1, 2
        coef(j) = 2._rkind/3._rkind*dx*exp(-xa(j)**2)
      enddo
      coef(nb_xint) = 1._rkind/3._rkind*dx*exp(-xa(nb_xint)**2)
!
      if (spec==1) then
         fact3 = 8._rkind/3._rkind/sqrt(PI)*ft/(fc*Tauee)
      else
         fact3 = 8._rkind/3._rkind/sqrt(PI)*ft/(fc*Tauii)
      end if

! dset=initialisation:
      call dset(Npoints,0._rkind,Ki11,1)
      call dset(Npoints,0._rkind,Ki12,1)
      call dset(Npoints,0._rkind,Ki22,1)

      do j = 1, nb_xint
	
         call dset(Npoints,0._rkind,nuaDv,1)
         call dset(Npoints,0._rkind,nuasv,1)
         call dset(Npoints,0._rkind,nuapv,1)
         call dset(Npoints,0._rkind,nuaEv,1)
         call dset(Npoints,0._rkind,nuatv,1)
         call dset(Npoints,1._rkind,F2,1)
         call dset(Npoints,1._rkind,F3,1)


         do i = 1, nspec
            xb = xa(j)*vt(spec)/vt(i)
            fact1 = abs(Z_spec(i)/Z_spec(spec))
            call daxpy(Npoints,fact1*(derf(xb)-chand(xb))/xa(j)**3,fact2,1,nuaDv,1)
            call daxpy(Npoints,fact1*2._rkind*(rtemp(spec)/rtemp(i)+ &
			(vt(spec)/vt(i))**2)*chand(xb)/xa(j), fact2,1,nuasv,1)
            call daxpy(Npoints,fact1*2._rkind*chand(xb)/xa(j)**3,fact2,1,nuapv,1)
         enddo
         call daxpy(Npoints,2._rkind,nuasv,1,nuaEv,1)
         call daxpy(Npoints,-2._rkind,nuaDv,1,nuaEv,1)
         call daxpy(Npoints,-1._rkind,nuapv,1,nuaEv,1)
         call daxpy(Npoints,1._rkind,nuaEv,1,nuatv,1)
         call daxpy(Npoints,3._rkind,nuaDv,1,nuatv,1)
         call mult(aux1,nuaDv,aux2,Npoints)
         call daxpy(Npoints,2.48_rkind/xa(j),aux2,1,F2,1)
         call mult(aux1,nuatv,aux2,Npoints)
         call daxpy(Npoints,1.9634_rkind/xa(j),aux2*epsK**(1.5),1,F3,1)
!         call daxpy(Npoints,1.9634_rkind/xa(j),aux2,1,F3,1)

         nuatotv = nuaDv/(F2*F3)


!
         if (spec==1) then
            int11 = coef(j)*nuatotv*Tauee*xa(j)**4
            int12 = int11*xa(j)**2
            int22 = int11*xa(j)**4
         else
            int11 = coef(j)*nuatotv*Tauii*xa(j)**4
            int12 = int11*xa(j)**2
            int22 = int11*xa(j)**4
         end if
         Ki11 = Ki11+fact3*int11
         Ki12 = Ki12+fact3*int12
         Ki22 = Ki22+fact3*int22
      enddo


      nuvi(:,1) = Ki11
      nuvi(:,2) = Ki12 - 5._rkind/2._rkind*Ki11
      nuvi(:,3) = Ki22 - 5._rkind*Ki12 + 25._rkind/4._rkind*Ki11

      end subroutine



      function derf(x)
!-------------------------------------------------------------------
!       compute the error function
!
      use prec_const
      implicit none
!
      real(rkind)                     :: derf
      real(rkind)                     :: x
      real(rkind), dimension(5)       :: p0
      real(rkind), dimension(4)       :: q0
      real(rkind), dimension(9)       :: p1
      real(rkind), dimension(8)       :: q1
      real(rkind), dimension(6)       :: p2
      real(rkind), dimension(5)       :: q2
      real(rkind)                     :: xmin
      real(rkind)                     :: xlarge
      real(rkind)                     :: res
      real(rkind)                     :: xsq
      real(rkind)                     :: xnum
      real(rkind)                     :: xden
      real(rkind)                     :: xi
      integer                         :: isw
      integer                         :: i

      p0(1)=113.8641541510502
      p0(2)=377.4852376853020
      p0(3)=3209.377589138469
      p0(4)=.1857777061846032
      p0(5)=3.161123743870566
      q0(1)=244.0246379344442
      q0(2)=1282.616526077372
      q0(3)=2844.236833439171
      q0(4)=23.60129095234412

      p1(1)=8.883149794388376
      p1(2)=66.11919063714163
      p1(3)=298.6351381974001
      p1(4)=881.9522212417691
      p1(5)=1712.047612634071
      p1(6)=2051.078377826071
      p1(7)=1230.339354797997
      p1(8)=2.153115354744038e-8
      p1(9)=.5641884969886701
      q1(1)=117.6939508913125
      q1(2)=537.1811018620099
      q1(3)=1621.389574566690
      q1(4)=3290.799235733460
      q1(5)=4362.619090143247
      q1(6)=3439.367674143722
      q1(7)=1230.339354803749
      q1(8)=15.74492611070983

      p2(1)=-3.603448999498044e-01
      p2(2)=-1.257817261112292e-01
      p2(3)=-1.608378514874228e-02
      p2(4)=-6.587491615298378e-04
      p2(5)=-1.631538713730210e-02
      p2(6)=-3.053266349612323e-01
      q2(1)=1.872952849923460
      q2(2)=5.279051029514284e-01
      q2(3)=6.051834131244132e-02
      q2(4)=2.335204976268692e-03
      q2(5)=2.568520192289822

      xmin = 1.0e-10
      xlarge = 6.375e0

      isw = 1
      if (x<0._rkind) then
         isw = -1
         x = -x
      endif
      if (x<0.477) then
         if (x>=xmin) then
            xsq = x*x
            xnum = p0(4)*xsq+p0(5)
            xden = xsq+q0(4)
            do i = 1,3
               xnum = xnum*xsq+p0(i)
               xden = xden*xsq+q0(i)
            enddo
            res = x*xnum/xden
         else
            res = x*p0(3)/q0(3)
         end if
      elseif (x<=4.0) then
         xsq = x*x
         xnum = p1(8)*x+p1(9)
         xden = x+q1(8)
         do i = 1,7
            xnum = xnum*x+p1(i)
            xden = xden*x+q1(i)
         enddo
         res = xnum/xden
         res = res*exp(-xsq)
         res = 1._rkind-res
      else
         xsq = x*x
         xi = 1._rkind/xsq
         xnum = p2(5)*xi+p2(6)
         xden = xi+q2(5)
         do i = 1,4
            xnum = xnum*xi+p2(i)
            xden = xden*xi+q2(i)
         enddo
         res = (1._rkind/sqrt(PI)+xi*xnum/xden)/x
         res = res*exp(-xsq)
         res = 1._rkind-res
      endif
      if (isw.eq.-1) res = -res
      derf = res
      return
      end function


      function chand(x)
!-------------------------------------------------------------------
!       compute the Chandrasekar function
!
      use prec_const
      implicit none
!
      real(rkind)                     :: chand
      real(rkind)                     :: derfx
      real(rkind)                     :: x, derf
!
      derfx = 2._rkind/sqrt(PI)*exp(-x**2)
      chand = (derf(x)-x*derfx)/(2._rkind*x**2)
!
      return
      end function

!-------------------------------------

subroutine mult(a,b,c,n)
  !     ------------------------
  ! c = a*b
  !
  use prec_const
  implicit none
  !
  integer                   :: n
  real(rkind), dimension(n) :: a
  real(rkind), dimension(n) :: b
  real(rkind), dimension(n) :: c
  !
  c=a*b     
  !
  return
end subroutine mult




!------------------------

subroutine dset(n,dx,dy,incy)
! -------- dy(:)=dx
	use prec_const
  	implicit none

     	 integer 			:: n, incy, iy, i
      	real(rkind)			:: dx
	real(rkind), dimension(n)	:: dy

      if (n.le.0) return
      iy = 1
      if (incy.lt.0) iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx
        iy = iy + incy
   10 continue
      return
      end subroutine dset

END MODULE neo_subroutines



