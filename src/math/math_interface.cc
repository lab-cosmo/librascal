

#include "math_interface.hh"


namespace rascal {
  namespace math {

    double hyp2f1(double& a,double& b,double& c,double& x ){
      return cephes::hyp2f1( a, b, c, x );
    }
  } // math
} // rascal
