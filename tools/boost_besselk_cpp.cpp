/* boost_besselk_cpp is a wrapper used for
 importing the cyl_bessel_k function from boost
 and fixing its template variables.
 WARNING: the boost cyl_bessel_k function
          accepts only scalar variables */
#include "boost_besselk_cpp.h"
/*BOOST modified bessel function of the second kind
  with fractional order
  inputs:
    nu: (double) bessel function fractional order
    x:  (double) bessel function variable
	outputs:
    boost_besselk_cpp: (double) modified bessel
                       function of the 2nd kind at x
                       and of order nu */
double boost_besselk_cpp(double const  &nu, double const &x){
  #ifdef USE_BOOST
    return boost::math::cyl_bessel_k(nu,x);
  #else
    return 0.e0;
  #endif
}
