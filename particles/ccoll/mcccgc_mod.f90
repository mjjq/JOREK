!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
!> MCCC operator for guiding center coordinates.
!> These coordinates are location X, momentum p (scalar), and
!> pitch = p_para/p. Implements Euler-Maruyama method.
!<
module mcccgc_mod
  use mccc_mod
  implicit none

  type mcccgc_coefstruct 
     real*8 :: kappa  = 0.D0
     real*8 :: dkappa = 0.D0
     real*8 :: Dpar   = 0.D0
     real*8 :: dDpar  = 0.D0
     real*8 :: Dperp  = 0.D0
     real*8 :: dDperp = 0.D0
     real*8 :: gx     = 0.D0
     real*8 :: nu     = 0.D0
  end type mcccgc_coefstruct

  public :: mcccgc_push, mcccgc_evalCoefs

contains

  !> Computes the collision coefficients and stores them in a struct that is used
  !> by mcccgc_push(). The coefficients are calculated separately like this so that
  !> they are not unnecessarily re-calculated if the time step is rejected.
  !> 
  !> input:
  !> 
  !> mccc_special dat -- struct that is obtained with the mccc_init() call
  !> real*8 ma        -- test particle mass [kg]
  !> real*8 qa        -- test particle charge [C]
  !> real*8 clogab(:) -- list of coulomb logarithms for species a colliding with b [1] 
  !> real*8 mb(:)     -- list of background species masses [kg]
  !> real*8 qb(:)     -- list of background species charges [C]
  !> real*8 nb(:)     -- list of background densities [1/m^3]
  !> real*8 thb(:)    -- list of normalized background temperatures thb=T_b/(m_b*c^2) [1] 
  !> real*8 u         -- normalized guiding center momentum magnitude u=p/mc [1]
  !> real*8 xi        -- guiding center pitch [1]
  !>
  !> output:
  !> 
  !> mcccgc_coefstruct coefs -- struct containing the coefficients needed for evaluating the collisions
  !<
  subroutine mcccgc_evalCoefs(dat,ma,qa,clogab,mb,qb,nb,thb,u,xi,coefs)
    type(mccc_special), intent(in) :: dat
    real*8, intent(in) :: ma,qa,u,xi
    real*8, intent(in) :: clogab(:),mb(:),qb(:),nb(:),thb(:) 

    type(mcccgc_coefstruct), intent(out) :: coefs

    real*8 :: kappab,dkappab,Dparb,dDparb,Dperpb,dDperpb ! Coefficients for species-wise collisions

    integer :: nspecies,i

    nspecies = size(clogab)
    coefs%kappa  = 0
    coefs%dkappa = 0
    coefs%Dpar   = 0
    coefs%dDpar  = 0
    coefs%Dperp  = 0
    coefs%dDperp = 0
    do i = 1,nspecies
       call mccc_coeffs(dat, ma, qa, clogab(i), mb(i), qb(i), nb(i), thb(i), u, &
            kappa  = kappab,   Dpar = Dparb,   Dperp = Dperpb, &
            dkappa = dkappab, dDpar = dDparb, dDperp = dDperpb)
       
       coefs%kappa  = coefs%kappa  + kappab
       coefs%dkappa = coefs%dkappa + dkappab
       coefs%Dpar   = coefs%Dpar   + Dparb
       coefs%dDpar  = coefs%dDpar  + dDparb
       coefs%Dperp  = coefs%Dperp  + Dperpb
       coefs%dDperp = coefs%dDperp + dDperpb
    end do

    coefs%nu = 2 * coefs%Dperp / u**2
  end subroutine mcccgc_evalCoefs
  
  !> Computes the value for guiding center phase space position after collisions with
  !> background species using Euler-Maruyama method with a fixed time step
  !>
  !> input:
  !> 
  !> mcccgc_coefstruct coefs -- struct containing the coefficients from mccca_evalCoefs() call
  !> real*8 cutoff           -- value (>= 0) used to mirror u as "u = 2*cutoff - u" if u < cutoff otherwise [1]
  !> real*8 uin              -- normalized test particle momentum u=p/mc before collisions [1]
  !> real*8 xiin             -- guiding center pitch before collisions [1]
  !> real*8 dtin             -- time step length [s]
  !>
  !> output:
  !> 
  !> real*8 uout  -- normalized guiding center momentum u=p/mc after collisions [1]
  !> real*8 xiout -- guiding center pitch after collisions [1]
  !<
  subroutine mcccgc_push(coefs,cutoff,uin,uout,xiin,xiout,dtin,rnd)
    type(mcccgc_coefstruct), intent(in) :: coefs
    
    real*8, intent(in) :: dtin ! time step length [s]
    real*8, intent(in) :: uin ! p/mc test particle momentum vector (px,py,pz) before collision, normalized to mc 
    real*8, intent(in) :: xiin !
    real*8, intent(in) :: cutoff
    real*8, intent(in) :: rnd(2)
    
    real*8, intent(out) :: uout 
    real*8, intent(out) :: xiout

    real*8 :: time
    real*8 :: bhat(3)
    real*8 :: F, gu, gxi
    real*8 :: kappa_k, kappa_d(2), erru, gx
    real*8 :: dWopt(2),dti,alpha
    integer :: i, windex,tindex, ki,kmax
    logical rejected

    F   = coefs%kappa+coefs%dDpar+2*coefs%Dpar/uin
    
    uout  = uin + ( coefs%kappa + coefs%dDpar + 2 * coefs%Dpar / uin ) * dtin &
                + sqrt( 2 * coefs%Dpar * dtin ) * rnd(1)
    xiout = xiin - xiin * coefs%nu * dtin + sqrt( ( 1.0 - xiin**2 ) * coefs%nu * dtin) * rnd(2)

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


  end subroutine mcccgc_push

end module mcccgc_mod
