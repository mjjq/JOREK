/* Header file containing the definitions of 
the boost_besselk_cpp function */
#ifndef BOOST_BESSELK_CPP_H
#define BOOST_BESSELK_CPP_H
  #ifdef USE_BOOST
    #include <boost/math/special_functions/bessel.hpp>
    using boost::math::cyl_bessel_k;
    using boost::math::policies::policy;
    using boost::math::policies::promote_float;
    using boost::math::policies::promote_double;
    using boost::math::policies::evaluation_error;
    using boost::math::policies::errno_on_error;
    /* defining a policy for avoiding a SIGABRT in case the
       cyl_bessel_k function of boost does not converges */
    typedef policy <
      promote_float<true>,  //< promote float to double for increasing precision
      promote_double<true>, //< promote double to long double for increasing precision
      evaluation_error<errno_on_error> // return best guess with EDOM=33
    > no_sigabrt_policy;
  #endif
  #ifdef __cplusplus
    extern "C" {
  #endif
  
  /* wrapper of the BOOST modified bessel function of the 
      2nd kind and with fractional order */
  //double boost_besselk_cpp(double const &nu, double const &x);
  double boost_besselk_cpp(double const &nu, double const &x);
  #ifdef __cplusplus
    }
  #endif
#endif

