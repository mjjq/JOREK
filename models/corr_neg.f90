module corr_neg

implicit none

  interface corr_neg_dens
    module procedure corr_neg_dens1, corr_neg_dens2
  end interface
      
  interface corr_neg_temp
    module procedure corr_neg_temp1, corr_neg_temp2
  end interface
contains
#include "corr_neg_include.f90"
end module corr_neg
