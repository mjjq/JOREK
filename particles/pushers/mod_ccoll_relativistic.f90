!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
!> Coulomb collisions for (relativistic) test particles.
!> See Sarkimaki et al, "Adaptive time-stepping Monte Carlo integration of Coulomb collisions",
!> Comp. Phys. Comm.
!<
module mod_ccoll_relativistic
  use constants
  use mod_bessel, only : bessel_k2exp, bessel_k1exp, bessel_k0exp
  use mod_simpson, only : simpson_adaptive, func_real8_1D
  implicit none

  real*8, parameter :: DEFAULT_L0L1_eps    = 1.D-8 !< default tolerance in eval_L0L1
  real*8, parameter :: DEFAULT_L0L1_cutoff = 1.D-7 !< default cutoff in evalL0L1
  
  ! Struct for storing tabulated values of special functions L0 and L1
  type ccoll_tabulatedL0L1
     real*8, allocatable, dimension(:)   :: u
     real*8, allocatable, dimension(:)   :: theta
     real*8, allocatable, dimension(:,:) :: L0
     real*8, allocatable, dimension(:,:) :: L1
  end type ccoll_tabulatedL0L1

  public :: ccoll_tabulatedL0L1, ccoll_compute_L0L1table, ccoll_free_L0L1table, &
       ccoll_read_L0L1table, ccoll_write_L0L1table, &
       ccoll_kinetic_relativistic_push, ccoll_clog, ccoll_coeffs, ccoll_gc_relativistic_push

  private

contains


  !> Evaluates the special functions L0 and L1.
  !> The cutoff parameter (< 1) divides the integral in two
  !> parts at u_cutoff = sqrt((1-theta * ln(cutoff))**2 -1).
  !> This helps to ensure that adaptive integration does not
  !> fail when theta is very small.
  subroutine eval_L0L1(u,theta,L0,L1,eps,cutoff)
    implicit none
    real*8, intent(in) :: u                !< p/mc value where functions are evaluated
    real*8, intent(in) :: theta            !< T/mc^2 value where functions are evaluated
    real*8, intent(in), optional :: eps    !< error tolerance for the adaptive simpsons rule
    real*8, intent(in), optional :: cutoff !< cutoff value

    real*8, intent(out) :: L0 !< evaluated L0 integral (Eq. 17 in referenced paper)
    real*8, intent(out) :: L1 !< evaluated L1 integral (Eq. 18 in referenced paper)

    real*8 :: def_cutoff, u_cutoff
    real*8 :: def_tol, tol
    procedure(func_real8_1d), pointer :: L0_ptr => null(), L1_ptr => null()

    def_tol    = DEFAULT_L0L1_eps
    def_cutoff = DEFAULT_L0L1_cutoff
    if(present(eps))    def_tol    = eps
    if(present(cutoff)) def_cutoff = cutoff

    u_cutoff = sqrt((1.D0 - theta * log(def_cutoff))**2 - 1.D0)

    L0_ptr => L0_integrand
    L1_ptr => L1_integrand

    ! the integral is evaluated in two parts if u > u_cutoff
    if(u_cutoff > u) then
       tol = 0.5D0 * def_tol * u
       L0 = simpson_adaptive(L0_ptr,0.D0,u,tol,10)
       L1 = simpson_adaptive(L1_ptr,0.D0,u,tol,10)
    else
       tol = 0.5D0 * def_tol * u_cutoff
       L0 = simpson_adaptive(L0_ptr,0.D0,u_cutoff,tol,10) + simpson_adaptive(L0_ptr,u_cutoff,u,tol,20)
       L1 = simpson_adaptive(L1_ptr,0.D0,u_cutoff,tol,10) + simpson_adaptive(L1_ptr,u_cutoff,u,tol,20)
    end if

  contains

    function L0_integrand(u) result (val)
      real*8, intent(in) :: u
      real*8 :: val,gamma
      gamma = sqrt( 1.D0 + u**2 )
      val   = exp( ( 1.D0 - gamma ) / theta ) / gamma
    end function L0_integrand

    function L1_integrand(u) result (val)
      real*8, intent(in) :: u
      real*8 :: val,gamma
      gamma = sqrt( 1.D0 + u**2 )
      val   = exp( (1.D0 - gamma ) / theta )
    end function L1_integrand

  end subroutine eval_L0L1

  !> Allocates and initializes tables containing computed  L0 and L1 values for interpolation
  type(ccoll_tabulatedL0L1) function ccoll_compute_L0L1table(uminxp,umaxxp,thminxp,thmaxxp,nu,nth,&
       eps,cutoff)
    implicit none

    real*8, intent(in)  :: uminxp  !> minimum u as umin=10^uminxp 
    real*8, intent(in)  :: umaxxp  !> maximum u as umax=10^umaxxp
    real*8, intent(in)  :: thminxp !> minimum theta as thmin=10^thminxp
    real*8, intent(in)  :: thmaxxp !> maximum theta as thmax=10^thmaxxp
    integer, intent(in) :: nu      !> number of u grid points
    integer, intent(in) :: nth     !> number of theta grid points
    real*8, intent(inout), optional :: eps    !> error tolerance for evaluating L0 and L1
    real*8, intent(inout), optional :: cutoff !> cutoff value for evaluating L0 and L1

    integer :: i, j
    type(ccoll_tabulatedL0L1) :: ccoll_tabulatedL0L1
    
    allocate(ccoll_tabulatedL0L1%u(nu), &
         ccoll_tabulatedL0L1%theta(nth),&
         ccoll_tabulatedL0L1%L0(nu,nth),&
         ccoll_tabulatedL0L1%L1(nu,nth))

    ! Set abscissae
    ccoll_tabulatedL0L1%u     = 10**( uminxp  + (/ ( ( i - 1 ) * ( umaxxp  - uminxp  ) / ( nu  - 1 ), i=1,nu)  /) )
    ccoll_tabulatedL0L1%theta = 10**( thminxp + (/ ( ( i - 1 ) * ( thmaxxp - thminxp ) / ( nth - 1 ), i=1,nth) /) )

    ! Evaluate and store values to the table
    do i=1,nu
       do j=1,nth
          if(present(eps) .and. present(cutoff)) then
             call eval_L0L1(ccoll_tabulatedL0L1%u(i),ccoll_tabulatedL0L1%theta(j),&
                  ccoll_tabulatedL0L1%L0(i,j),ccoll_tabulatedL0L1%L1(i,j),&
                  eps=eps,cutoff=cutoff)
          elseif(present(eps)) then
             call eval_L0L1(ccoll_tabulatedL0L1%u(i),ccoll_tabulatedL0L1%theta(j),&
                  ccoll_tabulatedL0L1%L0(i,j),ccoll_tabulatedL0L1%L1(i,j),&
                  eps=eps)
          elseif(present(cutoff)) then
             call eval_L0L1(ccoll_tabulatedL0L1%u(i),ccoll_tabulatedL0L1%theta(j),&
                  ccoll_tabulatedL0L1%L0(i,j),ccoll_tabulatedL0L1%L1(i,j),&
                  cutoff=cutoff)
          else
             call eval_L0L1(ccoll_tabulatedL0L1%u(i),ccoll_tabulatedL0L1%theta(j),&
                  ccoll_tabulatedL0L1%L0(i,j),ccoll_tabulatedL0L1%L1(i,j))
          end if
       end do
       
       ! This might last some time so keep user updated on progress
       write(*,*) "done: ",i,"/",nu
    end do

    ccoll_compute_L0L1table = ccoll_tabulatedL0L1
    
  end function ccoll_compute_L0L1table


  !> Deinitializes tabulated L0L1 values struct.
  subroutine ccoll_free_L0L1table(dat)
    implicit none
    type(ccoll_tabulatedL0L1), intent(inout) :: dat ! data to be deinitialized
   
    deallocate(dat%u,dat%theta,dat%L0,dat%L1)

  end subroutine ccoll_free_L0L1table


  !> Writes tabulated L0L1 values to a file.
  subroutine ccoll_write_L0L1table(dat,fn)
    implicit none
    type(ccoll_tabulatedL0L1), intent(in) :: dat ! data to be written
    character(len=*), intent(in) :: fn           ! output filename

    integer :: chn=989
    integer :: nu,nth

    nu  = size(dat%u)
    nth = size(dat%theta)

    open(unit=chn,file=fn,action='write')
    write(chn,*) nu
    write(chn,*) nth
    write(chn,*) dat%u
    write(chn,*) dat%theta
    write(chn,*) dat%L0
    write(chn,*) dat%L1
    close(chn)
    
  end subroutine ccoll_write_L0L1table

  !> Reads tabulated L0L1 values from a file.
  type(ccoll_tabulatedL0L1) function ccoll_read_L0L1table(fn)
    character(len=*), intent(in) :: fn !< input filename

    integer :: chn=989
    integer :: nu,nth

    open(unit=chn,file=fn,action='read')
    read(chn,*) nu
    read(chn,*) nth
    allocate(ccoll_read_L0L1table%u(nu), ccoll_read_L0L1table%theta(nth),&
         ccoll_read_L0L1table%L0(nu,nth), ccoll_read_L0L1table%L1(nu,nth))
    read(chn,*) ccoll_read_L0L1table%u
    read(chn,*) ccoll_read_L0L1table%theta
    read(chn,*) ccoll_read_L0L1table%L0
    read(chn,*) ccoll_read_L0L1table%L1
    close(chn)
    
  end function ccoll_read_L0L1table

  !> Computes Coulomb logarithm for a given test particle and plasma species.
  !> The logarithm is estimated as ln{lambda_D/min{bqm,bcl}} where lambda_D is
  !> the Debye length, and bqm and bcl are quantum mechanical and classical
  !> impact parameters, respectively. The Coulomb logarithm for different
  !> plasma species is returned. 
  subroutine ccoll_clog(ma,qa,mb,qb,nb,thb,u,clog)
    implicit none
    real*8, intent(in)  :: ma      !< test particle mass [kg]
    real*8, intent(in)  :: qa      !< test particle charge [C]
    real*8, intent(in)  :: mb(:)   !< list of background species masses [kg]
    real*8, intent(in)  :: qb(:)   !< list of background species charges [C]
    real*8, intent(in)  :: nb(:)   !< list of background densities [1/m^3]
    real*8, intent(in)  :: thb(:)  !< list of normalized background temperatures [T_b/(m_b*c^2)]
    real*8, intent(in)  :: u       !< normalized test particle momentum [p/mc]
    real*8, intent(out) :: clog(:) !< Coulomb logarithm for each species

    real*8  :: debyeLength ! Debye length accounting for all plasma species
    real*8  :: mr          ! Reduced mass 
    real*8  :: bcl         ! Classical impact parameter
    real*8  :: bqm         ! Quantum mechanical impact parameter
    integer :: nspec,i     ! Helper variables
    real*8, dimension(size(clog)) :: ubar ! Mean relative velocity

    debyeLength = sqrt( EPS_ZERO * SPEED_OF_LIGHT**2 / sum( ( nb * qb**2 ) / ( thb * mb ) ) )
    nspec = size(clog)
    ubar = SPEED_OF_LIGHT * sqrt( u**2 / ( 1 + u**2 ) * (/(1,i=1,nspec)/) + 3 * thb )

    do i=1,nspec
       mr  = ma * mb(i) / ( ma + mb(i) )
       bcl = qa * qb(i) / ( 4 * PI * EPS_ZERO * mr * ubar(i)**2 )
       bqm = HBAR / ( 2 * mr * ubar(i) )

       clog(i) = log(debyeLength/max(bcl,bqm))
    end do

  end subroutine ccoll_clog

  !> Computes requested collision coefficients and respective derivatives
  !> for a given test particle and background species.
  subroutine ccoll_coeffs(dat,ma,qa,clog,mb,qb,nb,thb,u,&
       K,dK,Dpar,dDpar,Dperp,dDperp,kappa,dkappa)
    implicit none

    type(ccoll_tabulatedL0L1), intent(in) :: dat !< tabulated L0L1 values
    real*8, intent(in)  :: ma   !< test particle mass [kg]
    real*8, intent(in)  :: qa   !< test particle charge [C]
    real*8, intent(in)  :: clog !< Coulomb logarithm for each species
    real*8, intent(in)  :: mb   !< background species mass [kg]
    real*8, intent(in)  :: qb   !< background species charge [C]
    real*8, intent(in)  :: nb   !< background density [1/m^3]
    real*8, intent(in)  :: thb  !< normalized background temperature [T_b/(m_b*c^2)]
    real*8, intent(in)  :: u    !< normalized test particle momentum [p/mc]

    real*8, intent(out), optional :: K      !< friction coefficient [1/s]
    real*8, intent(out), optional :: Dpar   !< parallel momentum diffusion [1/s]
    real*8, intent(out), optional :: Dperp  !< perpendicular momentum diffusion [1/s]
    real*8, intent(out), optional :: kappa  !< guiding center friction coefficient [1/s]
    real*8, intent(out), optional :: dK     !< derivative of K with respect to u [1/s]
    real*8, intent(out), optional :: dDpar  !< derivative of Dpar with respect to u [1/s]
    real*8, intent(out), optional :: dDperp !< derivative of Dperp with respect to u [1/s]
    real*8, intent(out), optional :: dkappa !< derivative of kappa with respect to u [1/s]

    real*8  :: Gab,mu0,mu1,mu2,gamma,dmu0,dmu1,dmu2 ! Special functions and coefficients
    real*8  :: u2,u3,u4,gamma2,gamma3 ! Helper variables
    
    gamma = sqrt(1.D0+u**2)
    Gab = nb * ( qa * qb )**2 * clog / ( 4 * pi * EPS_ZERO**2 * ma**2 * SPEED_OF_LIGHT**3 )
    u2  = u**2
    u3  = u**3

    if(present(dK) .or. present(dDpar) .or. present(dDperp) .or. present(dkappa)) then
       u4=u**4
       gamma2 = gamma**2
       gamma3 = gamma**3
       call ccoll_mufuncs(dat,u,thb,mu0,mu1,mu2,dmu0,dmu1,dmu2)
    else
       call ccoll_mufuncs(dat,u,thb,mu0,mu1,mu2)
    end if
    
    if(present(K)) then
       K = -Gab * ( mu0 / gamma + ( ma / mb) * mu1 ) / u2
    end if

    if(present(Dpar)) then
       Dpar = Gab * gamma * thb * mu1 / u3
    end if
    if(present(dDpar)) then
       dDpar = ( Gab * thb / ( gamma * u4 ) ) * ( gamma2 * u * dmu1 - ( 1.D0 + 2 * gamma2 ) * mu1 )
    end if

    if(present(Dperp)) then
       Dperp=Gab*(u2*(mu0+gamma*thb*mu2)-thb*mu1)/(2*gamma*u3)
    end if
    if(present(dDperp)) then
       dDperp=(Gab/(2*gamma3*u4))*((4*gamma2-1.D0)*thb*mu1-u2*((2*gamma2-1.D0)*mu0+thb*&
            gamma3*mu2)+gamma2*(u3*(dmu0+thb*gamma*dmu2)-thb*u*dmu1))
    end if

    if(present(kappa)) then
       kappa=-Gab*(ma/mb)*mu1/u2
    end if
    if(present(dkappa)) then
       dkappa=-Gab*(ma/mb)*(dmu1-2*mu1/u)/u2
    end if

  end subroutine ccoll_coeffs
  
  !> Computes the value for particle momentum after collisions with
  !> background species using Euler-Maruyama method with a fixed time step.
  subroutine ccoll_kinetic_relativistic_push(dat,ma,qa,mb,qb,nb,thb,dt,rnd,uin,uout)
    implicit none
    type(ccoll_tabulatedL0L1), intent(in) :: dat !< tabulated L0L1 values
    real*8, intent(in) :: ma      !< test particle mass [kg]
    real*8, intent(in) :: qa      !< test particle charge [C]
    real*8, intent(in) :: mb(:)   !< list of background species masses [kg]
    real*8, intent(in) :: qb(:)   !< list of background species charges [C]
    real*8, intent(in) :: nb(:)   !< list of background densities [1/m^3]
    real*8, intent(in) :: thb(:)  !< list of normalized background temperatures [T_b/(m_b*c^2)]
    real*8, intent(in) :: dt      !< time step length [s]
    real*8, intent(in) :: rnd(3)  !< array with three elements of standard normal random numbers ~ N(0,1)
    real*8, intent(in) :: uin(3)  !< normalized test particle momentum [p/mc]
    
    real*8, intent(out) :: uout(3) !< updated momentum [p/mc]

    real*8, allocatable :: clogab(:)
    real*8 :: K,Dpar,Dperp,Kb,Dparb,Dperpb ! the fokker-planck coefficients
    real*8 :: dW(3)                        ! the change in the Wiener process during dt
    real*8 :: uhat(3)                      ! a unit vector parallel to pin
    real*8 :: u                            ! absolute value of particle momentum normalized to mc
    integer  :: i, nspecies                ! for iterating over plasma species

    ! Wiener process for this step
    dW = sqrt(dt)*rnd

     ! Evaluate and sum Fokker-Planck coefficients
    u = norm2(uin)
    uhat = uin / u
    nspecies = size(mb)
    allocate(clogab(nspecies))
    call ccoll_clog(ma,qa,mb,qb,nb,thb,u,clogab)
    
    K        = 0.D0
    Dpar     = 0.D0
    Dperp    = 0.D0
    do i = 1,nspecies
       call ccoll_coeffs(dat,ma,qa,clogab(i),mb(i),qb(i),nb(i),thb(i),u,K=Kb,Dpar=Dparb,Dperp=Dperpb)
       K     = K     + Kb
       Dpar  = Dpar  + Dparb
       Dperp = Dperp + Dperpb
    end do
    deallocate(clogab)

    ! Use Euler-Maruyama method to get uout
    uout = uin + K * uhat * dt + sqrt( 2 * Dpar ) * dot_product( uhat, dW ) * uhat &
         + sqrt( 2 * Dperp ) * ( dW - dot_product( uhat, dW ) * uhat )

  end subroutine ccoll_kinetic_relativistic_push

  ! Apply Coulomb collisions for guiding center
  subroutine ccoll_gc_relativistic_push(dat,ma,qa,mb,qb,nb,thb,uin,uout,xiin,xiout,dt,rnd,cutoff)
    
    type(ccoll_tabulatedL0L1), intent(in) :: dat !< tabulated L0L1 values
    real*8, intent(in) :: ma     !< test particle mass [kg]
    real*8, intent(in) :: qa     !< test particle charge [C]
    real*8, intent(in) :: mb(:)  !< list of background species masses [kg]
    real*8, intent(in) :: qb(:)  !< list of background species charges [C]
    real*8, intent(in) :: nb(:)  !< list of background densities [1/m^3]
    real*8, intent(in) :: thb(:) !< list of normalized background temperatures [T_b/(m_b*c^2)]
    real*8, intent(in) :: dt     !< time step length [s]
    real*8, intent(in) :: uin    !< test particle momentum  [p/mc]
    real*8, intent(in) :: xiin   !< test particle pitch [ppar/p]
    real*8, intent(in) :: cutoff !< minimum normalized momentum, energies below this are reflected
    real*8, intent(in) :: rnd(2) !< normally ditributed random numbes

    real*8, intent(out) :: uout  !< updated momentum
    real*8, intent(out) :: xiout !< updated pitch

    real*8, allocatable :: clogab(:)
    real*8 :: kappa, Dpar, dDpar, Dperp, nu ! Collision coefficients
    real*8 :: kappab, Dparb, dDparb, Dperpb ! Coll. coefficients species-wise
    integer :: i, nspecies

    nspecies = size(mb)
    allocate(clogab(nspecies))
    call ccoll_clog(ma,qa,mb,qb,nb,thb,uin,clogab)

    kappa  = 0
    Dpar   = 0
    dDpar  = 0
    Dperp  = 0
    do i = 1,nspecies

       call ccoll_coeffs(dat, ma, qa, clogab(i), mb(i), qb(i), nb(i), thb(i), uin, &
            kappa  = kappab,   Dpar = Dparb,   Dperp = Dperpb, dDpar = dDparb)

       kappa  = kappa  + kappab
       Dpar   = Dpar   + Dparb
       dDpar  = dDpar  + dDparb
       Dperp  = Dperp  + Dperpb
    end do

    nu = 2 * Dperp / uin**2
    deallocate(clogab)

    uout  = uin + ( kappa + dDpar + 2 * Dpar / uin ) * dt &
                + sqrt( 2 * Dpar * dt ) * rnd(1)
    xiout = xiin - xiin * nu * dt + sqrt( ( 1.0 - xiin**2 ) * nu * dt ) * rnd(2)

    ! Reflect uout if uout is below the cutoff value
    if(uout .lt. cutoff) then
       uout = 2*cutoff-uout
    end if

    ! Reflect pitch if xiout is outside the interval [-1, 1]
    if(abs(xiout) .gt. 1.D0) then
       ! First make sure xiout is between the interval [-2, 2]. Physics-wise what we do here is not justified,
       ! but neither is having |xiout| > 2 (one should decrease time step if this happens).
       xiout = modulo( xiout, 2.0 )

       ! Reflect (this part is ok physics-wise)
       if(abs(xiout) .gt. 1.D0) then
          xiout = sign(2.D0-abs(xiout), xiout)
       end if
    end if

  end subroutine ccoll_gc_relativistic_push

  !> Evaluates the mu functions (and their derivatives if needed)
  subroutine ccoll_mufuncs(data,u,th,mu0,mu1,mu2,dmu0,dmu1,dmu2)
    implicit none
    type(ccoll_tabulatedL0L1), intent(in) :: data !< initialized L0L1 tables
    real*8, intent(in)  :: u   !< p/mc value
    real*8, intent(in)  :: th  !< T/mc^2
    real*8, intent(out) :: mu0 !< Eq. 14 in the referece paper
    real*8, intent(out) :: mu1 !< Eq. 15 in the referece paper
    real*8, intent(out) :: mu2 !< Eq. 16 in the referece paper
    real*8, intent(out), optional :: dmu0 !< d mu0 / d u
    real*8, intent(out), optional :: dmu1 !< d mu1 / d u
    real*8, intent(out), optional :: dmu2 !< d mu2 / d u

    real*8 :: gamma,gammasq,expBessel2,L0,L1,expgammatheta,tg,th2,u2,tgK

    gammasq       = 1.D0+u**2
    gamma         = sqrt(gammasq)
    expBessel2    = bessel_k2exp(1.D0 / th)
    expgammatheta = exp( ( 1.D0 - gamma ) / th )
    tg  = th * gamma
    th2 = th**2

    call interp_L0L1(data,u,th,L0,L1)

    mu0 = ( gammasq * L0 - th * L1 + ( th - gamma ) * u * expgammatheta ) / expBessel2
    mu1 = ( gammasq * L1 - th * L0 + ( th * gamma - 1 ) * u * expgammatheta ) / expBessel2
    mu2 = ( 2 * tg * L1 + ( 1 + 2 * th2 ) * u * expgammatheta ) / ( th * expBessel2 )

    if(present(dmu0) .or. present(dmu1) .or. present(dmu2)) then
       u2  = u**2
       tgK = tg * expBessel2
       if(present(dmu0)) then
          dmu0 = ( 2 * tg * u * L0 + ( gamma - 2 * th ) * u2 * expgammatheta ) / tgK
       end if
       if(present(dmu1)) then
          !dmu1=(2.D0*tg*u*L1+(2.D0*th2+1.D0)*u2*expgammatheta)/tgK ! This is the explicit form
          dmu1 = mu2 * u / gamma
       end if
       if(present(dmu2)) then
          dmu2 = ( 2 * th2 * u * L1 + ( 2 * th2 * tg + 2 * th2 + tg - u2 ) * expgammatheta ) / ( tgK * th )
       end if
    end if
    
  end subroutine ccoll_mufuncs


  !> Interpolates (bilinear) the special functions L0 and L1
  !> Approximations are used outside the tabulated domain.
  subroutine interp_L0L1(data,u,theta,L0,L1)
    type(ccoll_tabulatedL0L1), intent(in) :: data !< tabulated L0L1 values
    real*8, intent(in)  :: u     !< queried p/mc value
    real*8, intent(in)  :: theta !< queried T/mc^2 value
    real*8, intent(out) :: L0    !< interpolated L0
    real*8, intent(out) :: L1    !< interpolated L1

    real*8  :: th
    integer :: nu, nth

    nu = size(data%u)
    nth = size(data%theta)

    if ( theta.gt.data%theta(nth) ) then
       th=data%theta(nth) 
       print*,'Warning, temperature exceeds tabulated'
    else
       th=theta
    end if
    
    if(u.gt.data%u(nu)) then
       L0 = bessel_k0exp(1.D0/th)
       L1 = bessel_k1exp(1.D0/th)
    else
       if( (u.gt.data%u(1)) .and. (th.gt.data%theta(1)) ) then
          
          L0 = interp2(log10(data%u),log10(data%theta),data%L0,log10(u),log10(th)) 
          L1 = interp2(log10(data%u),log10(data%theta),data%L1,log10(u),log10(th))
       else
          L0 = sqrt(pi*theta/2)*erf(u/sqrt(2*theta))
          L1 = L0
       end if
    end if
    
  end subroutine interp_L0L1

  !> Bilinear interpolation for interpolating tabulated L0 and L1 values
  real(kind=8) function interp2(x, y, f, xq, yq)
    implicit none
    
    real*8, intent(in) :: x(:), y(:) !< abscissae
    real*8, intent(in) :: f(:,:)     !< function values
    real*8, intent(in) :: xq, yq     !< queried point

    real*8  :: dx, dy
    integer :: ix, iy

    dx = x(2) - x(1)
    dy = y(2) - y(1)
    
    ix = floor( ( xq - x(1) ) / dx ) + 1
    iy = floor( ( yq - y(1) ) / dy ) + 1

    interp2 = ( f(ix, iy)         * ( x(ix + 1) - xq ) * ( y(iy + 1) - yq ) &
              + f(ix + 1, iy)     * ( xq - x(ix) )     * ( y(iy + 1) - yq ) &
              + f(ix, iy + 1)     * ( x(ix + 1) - xq ) * ( yq - y(iy) )     &
              + f(ix + 1, iy + 1) * ( xq - x(ix) )     * ( yq - y(iy) )     &
              ) / ( dx * dy )
    
  end function interp2


end module mod_ccoll_relativistic


!< Simple program that generates a file "ccoll.data" that contains tabulated L0 and L1 values
!< that should be applicable for every fusion plasma.
!program ccoll_generate_L0L1
!  use mod_ccoll_relativistic
!  
!  implicit none
!
!  type(ccoll_tabulatedL0L1) :: dat
!  logical :: storage_file_on_disk
!  character(20), parameter :: storage_file='ccolldata'
!
!   real*8, parameter  :: uminxp  = -3.D0
!   real*8, parameter  :: umaxxp  = 2.D0
!   real*8, parameter  :: thminxp = -3.D0
!   real*8, parameter  :: thmaxxp = -1.D0 
!   integer, parameter :: nu      = 401
!   integer, parameter :: nth     = 201
!
!   ! Initialize the look-up tables
!  print*,''
!  write(*,*) 'Initializing look-up tables...'
!  inquire(file=storage_file,exist=storage_file_on_disk)
!  if(storage_file_on_disk) then
!     write(*,*) 'File already exists.'
!  else
!     dat = ccoll_compute_L0L1table(uminxp,umaxxp,thminxp,thmaxxp,nu,nth)
!     call ccoll_write_L0L1table(dat,storage_file)
!  end if
!  print*,'Done!'
!
!end program ccoll_generate_L0L1
