module corr_neg

implicit none

  interface corr_neg_dens
    module procedure corr_neg_dens1, corr_neg_dens2, corr_neg_dens3
  end interface

  interface dcorr_neg_dens_drho
    module procedure dcorr_neg_dens_drho1, dcorr_neg_dens_drho2, dcorr_neg_dens_drho3
  end interface
      
  interface corr_neg_temp
    module procedure corr_neg_temp1, corr_neg_temp2, corr_neg_temp3
  end interface

  interface dcorr_neg_temp_dT
    module procedure dcorr_neg_temp_dT1, dcorr_neg_temp_dT2, dcorr_neg_temp_dT3
  end interface

  interface d2corr_neg_temp_dT2
    module procedure d2corr_neg_temp_dT21, d2corr_neg_temp_dT22, d2corr_neg_temp_dT23
  end interface


contains

! Source code moved to an include file to be able to inline function calls into the loops
! in element_matrix. This module file is kept for those models which do not use inlining.
!
! For optimized element_matrix construction subroutines, inlining the corr_neg functions improves
! the performance. But the -ipo switch, which would enable inlining across compilation units, makes
! the compilation very slow and certain other subroutines become slower with it. As an alternative
! to -ipo, we separete the function definitions into the corr_neg_include.f90 file, which can be
! included into the source files that define element matrix. This way the compiler will be able to
! inline the corr_neg functions.
!
! The include directive is used here to avoid code duplication.
#include "corr_neg_include.f90"
end module corr_neg
