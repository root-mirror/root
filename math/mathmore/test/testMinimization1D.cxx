#include "Math/Polynomial.h"
#include "Math/GSLMinimizer1D.h"
#include "Math/Functor.h"
#include <iostream>


void testMinimization1D() {


  ROOT::Math::Polynomial *polyf = new ROOT::Math::Polynomial(2);

  std::vector<double> p(3);
  p[0] = 1;
  p[1] = -4;
  p[2] = 1;
  polyf->SetParameters(&p[0]);
  //ROOT::Math::Functor1D<ROOT::Math::Base> func(*polyf);
  ROOT::Math::IGenFunction & func = *polyf;


  {
     // default (Brent)
     ROOT::Math::GSLMinimizer1D min;
     min.SetFunction(func,1,-10,10);
     min.Minimize(100,0.01,0.01);
     std::cout << "test Min1D " << min.Name() << "  Return code " << min.Status() << std::endl;

     std::cout.precision(20);

     std::cout << "Found minimum: x = " << min.XMinimum() << "  f(x) = " << min.FValMinimum() << std::endl;

  }
  {
     // Golden Section
     ROOT::Math::GSLMinimizer1D min(ROOT::Math::Minim1D::kGOLDENSECTION);
     min.SetFunction(func,1,-10,10);
     min.Minimize(100,0.01,0.01);
     std::cout << "test Min1D " << min.Name() << "  Return code " << min.Status() << std::endl;

     std::cout.precision(20);

     std::cout << "Found minimum: x = " << min.XMinimum() << "  f(x) = " << min.FValMinimum() << std::endl;

  }


}



int main() {

  testMinimization1D();
  return 0;

}
