module fgauss
contains
FUNCTION FGAUS(ZS,BGF,XR1,XR2,SIG1,SIG2,FACT,DFGAUSS)
!-----------------------------------------------------------------------
!     BGF + (1 - BGF) * (GAUSS1 + FACT * GAUSS2) / FACT
!-----------------------------------------------------------------------
implicit none
real*8 :: ZS, BGF, XR1, XR2, SIG1, SIG2, FACT, DFGAUSS
real*8 :: ZNORM1, ZNORM2, ZEX1, ZEX2, DEX1, DEX2, F1, F2, DF1, DF2
real*8 :: FGAUS

ZNORM1 = 0.39894d0 / SIG1
ZNORM2 = 0.39894d0 / SIG2
ZEX1   = -0.5d0 * (ZS - XR1)**2 / SIG1**2
ZEX2   = -0.5d0 * (ZS - XR2)**2 / SIG2**2
DEX1   = -(ZS-XR1)/SIG1**2
DEX2   = -(ZS-XR2)/SIG2**2

F1     = ZNORM1 * EXP(ZEX1)
F2     = ZNORM2 * EXP(ZEX2)
DF1    = ZNORM1 * DEX1 * EXP(ZEX1)
DF2    = ZNORM2 * DEX2 * EXP(ZEX2)

FGAUS  = BGF + (1.d0 - BGF) * (F1 + FACT * F2) / FACT
DFGAUSS = (1.d0-BGF) * (DF1 + FACT * DF2) / FACT

RETURN
END FUNCTION FGAUS

FUNCTION FGAUS3(ZS,BGF,XR1,XR2,XR3,SIG1,SIG2,SIG3,FACT,DFGAUSS3)
!-----------------------------------------------------------------------
!     BGF + (1 - BGF) * (GAUSS1 + FACT * GAUSS2) / FACT
!-----------------------------------------------------------------------
implicit none
real*8 :: ZS, BGF, XR1, XR2, XR3, SIG1, SIG2, SIG3, FACT, DFGAUSS3
real*8 :: ZNORM1, ZNORM2, ZNORM3, ZEX1, ZEX2, ZEX3
real*8 :: DEX1, DEX2, DEX3, F1, F2, F3, DF1, DF2, DF3
real*8 :: FGAUS3

ZNORM1 = 0.39894d0 / SIG1
ZNORM2 = 0.39894d0 / SIG2
ZNORM3 = 0.39894d0 / SIG3

ZEX1   = -0.5d0 * (ZS - XR1)**2 / SIG1**2
ZEX2   = -0.5d0 * (ZS - XR2)**2 / SIG2**2
ZEX3   = -0.5d0 * (ZS - XR3)**2 / SIG3**2

DEX1   = -(ZS-XR1)/SIG1**2
DEX2   = -(ZS-XR2)/SIG2**2
DEX3   = -(ZS-XR3)/SIG3**2

F1     = ZNORM1 * EXP(ZEX1)
F2     = ZNORM2 * EXP(ZEX2)
F3     = ZNORM3 * EXP(ZEX3)

DF1    = ZNORM1 * DEX1 * EXP(ZEX1)
DF2    = ZNORM2 * DEX2 * EXP(ZEX2)
DF3    = ZNORM3 * DEX3 * EXP(ZEX3)

FGAUS3  = BGF + (1.d0 - BGF) * (F1 + FACT * F2 + F3) / FACT
DFGAUSS3 = (1.d0-BGF) * (DF1 + FACT * DF2 + DF3) / FACT

RETURN
END FUNCTION FGAUS3

end module fgauss
