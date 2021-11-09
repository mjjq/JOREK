/* Header file containing the definitions of 
the boost_besselk_cpp function */
#ifndef BOOST_BESSELK_CPP_H
#define BOOST_BESSELK_CPP_H
  #ifdef USE_BOOST
    #include <boost/math/special_functions/bessel.hpp>
  #endif
  #ifdef __cplusplus
    extern "C" {
  #endif

  /* wrapper of the BOOST modified bessel function of the 
      2nd kind and with fractional order */
  double boost_besselk_cpp(double const &nu, double const &x);
  
  #ifdef __cplusplus
    }
  #endif
#endif

