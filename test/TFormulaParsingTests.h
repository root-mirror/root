#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TFormula.h"
#include "TGraph.h"
#include "TMath.h"
#include "Math/ChebyshevPol.h"
#include "TError.h"
#include "TFile.h"
#include "TMacro.h"
#include "TSystem.h"

#include <limits>
#include <cstdlib>
#include <stdio.h>
// test of tformula neeeded to be run


class TFormulaParsingTests {


bool verbose;
std::vector<int> failedTests;

// We need a softer way to reason about equality in 32 bits
// Being this a quick test, doing the check at runtime is really no problem.
bool fpEqual(double x, double y, bool epsilon = false)
{
   bool isEqual = epsilon ? std::abs(x-y) <= std::numeric_limits<double>::epsilon() : x == y;
   if (!isEqual) {
       // std::hexfloat not there for older gcc versions
       printf("\nThe numbers differ: %A and %A\n", x, y);
   }
   return isEqual;
}

public:

TFormulaParsingTests(bool _verbose = false) : verbose(_verbose) {}

bool test1() {
   // test composition of functions
   TF1 f1("f1","[0]+[1]*x*x");
   TF1 f2("f2","[0]+[1]*f1");

   f2.SetParameters(1,2,3,4);

   return (f2.Eval(2) == 39.);

}

bool test2() {

   TF1 f1("f1","[0]+[1]*x");
   TF1 f2("f2","[0]+[1]*x*f1");

   TF1 f3("f3",f2.GetExpFormula() );

   f3.SetParameters(1,2,3,4);

   return (f3.Eval(2) == 45.);

}


bool test3() {

   // still tets composition of functions
   TF1 f1("f1","gaus");
   TF1 f2("f2","[0]+[1]*x+f1");


   f2.SetParameters(10,2,5,2,1);

   f1.SetParameters(5,2,1);

   return (f2.Eval(2) == (10. + 2*2 + f1.Eval(2)) );

}

bool test4() {

   // similar but with different name (it contains gaus)
   TF1 f1("fgaus","gaus");
   TF1 f2("f2","[0]+[1]*x+fgaus");


   f2.SetParameters(10,2,5,2,1);

   f1.SetParameters(5,2,1);

   return (f2.Eval(2) == (10. + 2*2 + f1.Eval(2)) );

}
bool test5() {

   // similar but with different name (it contains gaus)
   TF1 f1("gausnfunc","gaus");
   TF1 f2("f2","[0]+[1]*x+gausnfunc");


   f2.SetParameters(10,2,5,2,1);

   f1.SetParameters(5,2,1);

   return (f2.Eval(2) == (10. + 2*2 + f1.Eval(2)) );

}

bool test1a() {
   // this makes infinite loop
   // why re-using same name
   TF1 f1("f1","[0]+[1]*x*x");
   TF1 f2("f1","[0]+[1]*f1");
   return true;
}


bool test6() {
   // test linear function used in fitting
   bool ok = true;
   double x[] = {1,2,3,4,5};
   double y[] = {1,4,7,9,10};
   TGraph g(5,x,y);

   int iret = g.Fit("x++1","Q");
   ok &= (iret == 0);
   iret = g.Fit("1++x","Q");
   return iret == 0;
}

bool test7() {
   // test copying and deleting of linear functions
   TF1 * f1 = new TF1("f1","1++x");
   if (f1->GetNpar() != 2) return false;
   f1->SetParameters(2,3);
   if (f1->Eval(3) != 11) return  false;

   if (verbose) printf("Test7: test linear part1 of function\n");
   TFormula * lin1 = (TFormula*) f1->GetLinearPart(1);
   assert (lin1);
   if (lin1->Eval(3) != 3) return false;

   if (verbose) printf("Test7: test copying linear function\n");

   TF1 * f2 = new TF1(*f1);
   if (f2->Eval(3) != 11) return  false;

   if (verbose) printf("Test7: test linear part1 of copied function\n");
   if (!f2->IsLinear()) return false;
   lin1 = (TFormula*) f2->GetLinearPart(1);
   assert (lin1);
   if (lin1->Eval(3) != 3) return false;

   delete f1;

   if (verbose) printf("Test7: test cloning linear function\n");

   TF1 * f3 = (TF1*) f2->Clone("f3");
   if (f3->Eval(3) != 11) return  false;

   if (verbose) printf("Test7: test deleting the copied function\n");
   delete f2;

   if (verbose) printf("Test7: test linear part1 of cloned function\n");
   if (!f3->IsLinear()) return false;
   lin1 = (TFormula*) f3->GetLinearPart(1);
   assert (lin1);
   if (verbose) printf("Test7: test evaluating linear part1 of cloned function\n");
   if (lin1->Eval(3) != 3) return false;


   if (verbose) printf("Test7: test deleting the cloned function\n");
   delete f3;
   return true;
}

bool test8() {
   // test the operator ^
   bool ok = true;
   TFormula * f = 0;
   f = new TFormula("f","x^y");

   ok &= (f->Eval(2,3) == 8);
   delete f;

   f = new TFormula("f","(x+[0])^y");
   f->SetParameter(0,1);
   ok &= (f->Eval(2,3) == 27);
   delete f;

   f = new TFormula("f","sqrt(x+[0])^y");
   f->SetParameter(0,2);
   ok &= (f->Eval(2,3) == 8);
   delete f;

   f = new TFormula("f","[0]/((x+2)^y)");
   f->SetParameter(0,27);
   ok &= (f->Eval(1,3) == 1);
   delete f;

   f = new TFormula("f","[0]/((x+2)^(y+1))");
   f->SetParameter(0,27);
   ok &= (f->Eval(1,2) == 1);
   delete f;

   // test also nested operators
   f = new TFormula("f","((x+1)^y)^z");
   ok &= (f->Eval(1,3,4) == 4096);
   delete f;

   f = new TFormula("f","x^((y+1)^z)");
   ok &= (f->Eval(2,1,3) == 256);
   delete f;

   return ok;

}

bool test9() {
   // test the exponent notations in numbers
   bool ok = true;

   TFormula * f = 0;
   f = new TFormula("f","x+2.0e1");
   ok &= (f->Eval(1) == 21.);

   f = new TFormula("f","x*2.e-1");
   ok &= (f->Eval(10) == 2.);

   f = new TFormula("f","x*2.e+1");
   ok &= (f->Eval(0.1) == 2.);

   f = new TFormula("f","x*2E2");
   ok &= (f->Eval(0.01) == 2.);

   delete f;
   return ok;


}

bool test10() {
   // test the operator "? : "
   bool ok = true;
   TFormula * f = 0;
   f  = new TFormula("f","(x<0)?-x:x");
   ok &= (f->Eval(-2) == 2);
   ok &= (f->Eval(2) == 2);

   f = new TFormula("f","(x<0)?x:pol2");
   f->SetParameters(1,2,3);
   ok &= (f->Eval(-2) == -2);
   ok &= (f->Eval(2) == 1 + 2*2 + 2*2*3);

   delete f;
   return ok;
}

bool test11() {
   // test with ::
   bool ok = true;

   TFormula f1("f","ROOT::Math::normal_pdf(x,1,2)");
   TFormula f2("f","[0]+TMath::Gaus(x,2,1,true)");
   f2.SetParameter(0,1);
   ok &= ( (f1.Eval(2) +1. ) == f2.Eval(2) );
   return ok;
}
bool test12() {
   // test parameters order
   bool ok = true;
   TFormula * f = 0;
   f = new TFormula("f","[2] + [3]*x + [0]*x^2 + [1]*x^3");
   f->SetParameters(1,2,3,4);
   double result = 3+4*2+1*4+2*8;
   ok &= (f->Eval(2) == result);
   f = new TFormula("f","[b] + [c]*x + [d]*x^2 + [a]*x^3");
   f->SetParameters(1,2,3,4);
   result = 2+3*2+4*4+1*8;
   ok &= (f->Eval(2) == result);
   // change a parameter value
   f->SetParName(2,"par2");
   ok &= (f->Eval(2) == result);
   return ok;
}

bool test13()  {
   // test GetExpFormula
   TFormula f("f","[2] + [0]*x + [1]*x*x");
   f.SetParameters(1,2,3);
   return (f.GetExpFormula() == TString("[p2]+[p0]*x+[p1]*x*x"));
}
bool test14()  {
   // test GetExpFormula
   TFormula f("f","[2] + [0]*x + [1]*x*x");
   f.SetParameters(1,2,3);
   return (f.GetExpFormula("P") == TString("3+1*x+2*x*x"));
}
bool test15()  {
   // test GetExpFormula
   TFormula f("f","[2] + [0]*x + [1]*x*x");
   f.SetParameters(1,2,3);
   return (f.GetExpFormula("CLING") == TString("p[2]+p[0]*x[0]+p[1]*x[0]*x[0] ") ); // need an extra white space
}
bool test16()  {
   // test GetExpFormula
   TFormula f("f","[2] + [0]*x + [1]*x*x");
   f.SetParameters(1,2,3);
   return (f.GetExpFormula("CLING P") == TString("3.000000+1.000000*x[0]+2.000000*x[0]*x[0] ") );
}

bool test17() {
   // test Eval for TF1
   TF1 * f1 = new TF1("f1","[0]*sin([1]*x)");
   f1->SetParameters(2,3);
   TF1 * f0 = new TF1("f0",[](double *x, double *p){ return p[0]*sin(p[1]*x[0]); },0,10,2);
   f0->SetParameters(2,3);
   bool ok = true;
   ok &= fpEqual(f1->Eval(1.5) , f0->Eval(1.5) );
   double xx[1] = {2.5};
   ok &= fpEqual(f1->EvalPar(xx) , f0->Eval(2.5) );
   return ok;
}

bool test18() {
   // test Eval for TF2
   TF2 * f1 = new TF2("f2","[0]*sin([1]*x*y)");
   f1->SetParameters(2,3);
   TF2 * f0 = new TF2("f0",[](double *x, double *p){ return p[0]*sin(p[1]*x[0]*x[1]); },0,10,0,10,2);
   f0->SetParameters(2,3);
   bool ok = true;
   ok &= fpEqual(f1->Eval(1.5,2.5) , f0->Eval(1.5,2.5) );
   double par[2] = {3,4};
   double xx[2] = {0.8,1.6};
   ok &= fpEqual(f1->EvalPar(xx,par) , f0->EvalPar(xx,par) );
   return ok;
}

bool test19() {
   // test Eval for TF3
   TF3 * f1 = new TF3("f3","[0]*sin([1]*x*y*z)");
   f1->SetParameters(2,3);
   TF3 * f0 = new TF3("f0",[](double *x, double *p){ return p[0]*sin(p[1]*x[0]*x[1]*x[2]); },0,10,0,10,0,10,2);
   f0->SetParameters(2,3);
   bool ok = true;
   ok &= fpEqual(f1->Eval(1.5,2.5,3.5) , f0->Eval(1.5,2.5,3.5) );
   double par[2] = {3,4};
   double xx[3] = {0.8,1.6,2.2};
   ok &= fpEqual(f1->EvalPar(xx,par) , f0->EvalPar(xx,par) );
   return ok;
}

bool test20() {
   // test parameter order with more than 10 parameters
   TF2 f2("f2","xygaus+xygaus(5)+xygaus(10)+[offset]");
   double params[16] = {1,0,1,1,1, 2,-1,2,0,2, 2,1,3,-1,2, 10};
   f2.SetParameters(params);
   TF2 f0("f2",[](double *x, double *p){ return p[0]*TMath::Gaus(x[0],p[1],p[2])*TMath::Gaus(x[1],p[3],p[4]) +
            p[5]*TMath::Gaus(x[0],p[6],p[7])*TMath::Gaus(x[1],p[8],p[9]) +
            p[10]*TMath::Gaus(x[0],p[11],p[12])*TMath::Gaus(x[1],p[13],p[14]) + p[15]; },
      -10,10,-10,10,16);
   double xx[2]={1,2};
   //printf(" difference = %f , value %f \n", f2.Eval(1,2) - f0.EvalPar(xx,params), f2.Eval(1,2) );
   return fpEqual( f2.Eval(1,2) , f0.EvalPar(xx,params) );
}

bool test21() {
   // test parsing polynomials (bug ROOT-7312)
   TFormula f("f","pol2+gaus(3)");
   f.SetParameters(1,2,3,1,0,1);
   TF1 f0("f0",[](double *x, double *p){ return p[0]+x[0]*p[1]+x[0]*x[0]*p[2]+p[3]*TMath::Gaus(x[0],p[4],p[5]); },0,1,6);
   f0.SetParameters(f.GetParameters() );
   return fpEqual(f.Eval(2) , f0.Eval(2) );
}

bool test22() {
   // test chebyshev
   TF1 f("f","cheb10+[offset]");
   double p[12] = {1,1,1,1,1,1,1,1,1,1,1,10 };
   f.SetParameters(p);
   return (f.Eval(0.5) == ROOT::Math::ChebyshevN(10, 0.5, p ) + f.GetParameter("offset"));
}

bool test23() {
   // fix function compositions using pre-defined functions
   bool ok = true;
   TF1 f1("f1","gaus");
   TF1 f2("f2","[0]+f1");
   TF1 f0("f0",[](double *x, double *p){ return p[0]+p[1]*TMath::Gaus(x[0],p[2],p[3]); },-3,3,4 );
   f2.SetParameters(10,1,0,1);
   f0.SetParameters(f2.GetParameters() );
   ok &= fpEqual(f2.Eval(1) , f0.Eval(1) );

   TF1 f3("f3","f1+[0]");
   // param order should be the same
   f3.SetParameters( f2.GetParameters() );
   ok &= fpEqual(f3.Eval(1) , f0.Eval(1) );
   return ok;
}

bool test24() {
   // test I/O for parameter ordering
   bool ok = true;
   TF2 f("f","xygaus");
   f.SetParameters(10,0,1,-1,2);
   TF2 * f2 = (TF2*) f.Clone();
   ok &= ( f.Eval(1,1) == f2->Eval(1,1) );
   // test with copy
   TF2 f3(f);
   ok &= ( f.Eval(1,1) == f3.Eval(1,1) );
   return ok;
}

bool test25() {
   // fix parsing of operator^ (ROOT-7349)
   bool ok = true;
   TF1 f1("f1","x^-2.5");
   ok &= (f1.Eval(3.) == TMath::Power(3,-2.5) );
   if (!ok)  std::cout << "Error in test25 - f != x^-2.5 " << f1.Eval(3.) << "  " <<  TMath::Power(3,-2.5) << std::endl;

   TF1 f2("f2","x^+2.5");
   //TF1 f3("f3","std::pow(x,2.5)");  // this needed to be fixed
   TF1 f3("f3","TMath::Power(x,2.5)");
   bool ret =  (f2.Eval(3.) == f3.Eval(3) );
   if (!ret)  std::cout << "Error in test25 - f2 != f3 " << f2.Eval(3.) << "  " <<  f3.Eval(3.) << std::endl;
   ok &= ret;

   //cms test
   TF1 t1("t1","(x<190)?(-18.7813+(((2.49368+(10.3321/(x^0.881126)))*exp(-((x^-1.66603)/0.074916)))-(-17.5757*exp(-((x^-1464.26)/-7.94004e+06))))):(1.09984+(0.394544*exp(-(x/562.407))))");
   double x = 2;
   double y =(x<190)?(-18.7813+(((2.49368+(10.3321/(std::pow(x,0.881126))))*exp(-((std::pow(x,-1.66603))/0.074916)))-(-17.5757*exp(-((std::pow(x,-1464.26))/-7.94004e+06))))):(1.09984+(0.394544*exp(-(x/562.407))));
   // this fails on 32 bits - put a tolerance
   ret = TMath::AreEqualAbs(t1.Eval(2) , y , 1.E-8);
   if (!ret)  std::cout << "Error in test25 - t1 != y " << t1.Eval(2.) << "  " <<  y << std::endl;
   ok &= ret;

   // tests with scientific notations
   auto ff = new TFormula("ff","x+2.e-2^1.2e-1");
   ret = ( ff->Eval(1.) == (1. + std::pow(2.e-2,1.2e-1) ) );
   if (!ret) std::cout << "Error in test25 - ff != expr " << ff->Eval(1.) << "  " <<   (1. + std::pow(2.e-2,1.2e-1) ) << std::endl;
   ok &= ret;

   ff = new TFormula("ff","x^-1.2e1");
   ret = ( ff->Eval(1.5) == std::pow(1.5,-1.2e1) ) ;
   if (!ret) std::cout << "Error in test25 - ff(1.5) != pow " <<  ff->Eval(1.5) << "  " <<  std::pow(1.5,-1.2e1) << std::endl;
   ok &= ret;

   ff = new TFormula("ff","1.5e2^x");
   ret = ( ff->Eval(2) == std::pow(1.5e2,2) );
   if (!ret) std::cout << "Error in test25 - ff(2) != pow " << ff->Eval(2) << "  " <<  std::pow(1.5e2,2) << std::endl;
   ok &= ret;

   ff = new TFormula("ff","1.5e2^x^-1.1e-2");
   ret = ( ff->Eval(2.) == std::pow(1.5e2, std::pow(2,-1.1e-2) ) );
   if (!ret) std::cout << "Error in test25 - ff(2) != pow^pow " << ff->Eval(2.) << "  " <<  std::pow(1.5e2, std::pow(2,-1.1e-2) ) << std::endl;
   ok &= ret;

   // test same prelacements
   ff = new TFormula("ff","pol10(3)+pol2");
   std::vector<double> p = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
   ff->SetParameters(p.data() );
   double sum = 0; for (auto &a : p) { sum+= a;}
   ret = ( ff->Eval(1.) == sum );
   if (!ret) std::cout << "Error in test25 - ff(1) != sum " << ff->Eval(1.) << "  " <<  sum << std::endl;
   ok &= ret;

   return ok;
}

bool test26() {
   // test sign function
   bool ok = true;
   TF1 f("f","x*sign(1.,x+2.)");
   ok &= (f.Eval(2) == 2);
   ok &= (f.Eval(-1) == -1);
   ok &= (f.Eval(-3) == 3);
   TF1 f2("f2","x*TMath::Sign(1,x+2)");
   ok &= (f2.Eval(2) == 2);
   ok &= (f2.Eval(-1) == -1);
   ok &= (f2.Eval(-3) == 3);
   TF1 f3("f3","TMath::SignBit(x-2)");
   ok &= (f3.Eval(1) == 1);
   ok &= (f3.Eval(3) == 0);
   return ok;
}

bool test27() {
   // test ssq function
   bool ok = true;
   TF1 f1("f1","x+sq(x+2)+sq(x+[0])");
   TF1 f2("f2","x+(x+2)^2+(x+[0])^2");
   f1.SetParameter(0,3); f2.SetParameter(0,3);
   ok &= fpEqual(f1.Eval(2) , f2.Eval(2));
   ok &= fpEqual(f1.Eval(-4) , f2.Eval(-4));
   // test nested expressions and conflict with sqrt
   TF1 f3("f3","sqrt(1.+sq(x))");
   ok &= fpEqual(f3.Eval(2) , sqrt(5) );
   TF1 f4("f4","sq(1.+std::sqrt(x))");
   ok &= fpEqual(f4.Eval(2) , TMath::Sq(1.+sqrt(2)) );
   TF1 f5("f5","sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))");
   auto func = [](double *x, double *p){ return TMath::Sqrt(((TMath::Sign(1,p[0])*TMath::Sq(p[0]/x[0]))+(TMath::Sq(p[1])*(TMath::Power(x[0],(p[3]-1)))))+TMath::Sq(p[2])); };
   TF1 f6("f6",func,-10,10,4);
   f5.SetParameters(-1,2,3,4); f6.SetParameters(f5.GetParameters());
   ok &= fpEqual(f5.Eval(2) , f6.Eval(2) );
   return ok;
}

bool test28() {
   bool ok = true;
   // test composition of two functions
   TF1 fsin("fsin", "[0]*sin(x)", 0., 10.);
   fsin.SetParNames( "sin");
   fsin.SetParameter( 0, 2.1);

   TF1  fcos("fcos", "[0]*cos(x)", 0., 10.);
   fcos.SetParNames( "cos");
   fcos.SetParameter( 0, 1.1);

   TF1 fsincos("fsc", "fsin+fcos");

   // keep same order in evaluation
   TF1 f0("f0",[](double *x, double *p){ return p[1]*sin(x[0]) + p[0]*cos(x[0]);},0.,10.,2);
   f0.SetParameters(1.1,2.1);
#ifdef R__B64
   bool epsilon = false;
#else
   bool epsilon = true;
#endif
   ok &= fpEqual(fsincos.Eval(2) , f0.Eval(2), epsilon);
   return ok;

}

bool test29() {
   // test hexadecimal numbers
   bool ok = true;
   TF1 f1("f1","x+[0]*0xaf");
   f1.SetParameter(0,2);
   ok &= (f1.Eval(3) == (3.+2*175.) );

   TF1 f2("f2","0x64^2+x");
   ok &= (f2.Eval(1) == 10001 );

   TF1 f3("f3","x^0x000c+1");
   ok &= (f3.Eval(2) == 4097 );

   return ok;

}

bool test30() {
// handle -- (++ is in linear expressions)
   bool ok = true;
   TF1 f1("f1","x--[0]");
   f1.SetParameter(0,2);
   ok &= (f1.Eval(3) == 5. );

   return ok;

}

bool test31() {
// test whitespaces in par name and cloning
   bool ok = true;
   TF1 f1("f1","x*[0]");
   f1.SetParameter(0,2);
   f1.SetParName(0,"First Param");
   auto f2 = (TF1*) f1.Clone();

   ok &= (f1.Eval(3) == f2->Eval(3) );
   ok &= (TString(f1.GetParName(0) ) == TString(f2->GetParName(0) ) );

   return ok;

}

bool test32() {
// test polynomial are linear and have right number
   bool ok = true;
   TF1 f1("f1","pol2");
   ok &= (f1.GetNumber() == 302);
   ok &= (f1.IsLinear() );

   TF1 f2("f2","gaus(0)+pol1(3)");
   ok &= (f2.GetNumber() == 0);
   ok &= (!f2.IsLinear());
   return ok;
}

bool test33() {
   // test new bigaus pre-defined funcition
   bool ok = true;
   TF2 f1("f1","bigaus",-10,10,-10,10);
   ok &= (f1.GetNumber() == 112);
   ok &= (std::string(f1.GetParName(5)) == "Rho");
   f1.SetParameters(1,0,1,1,2,0.);
   TF2 f2("f2","xygaus",-10,10,-10,10);
   f2.SetParameters(1,0,1,1,2);
   ok &= TMath::AreEqualAbs( f1.Eval(0), f2.Eval(0)/(f2.Integral(-10,10,-20,20) ), 1.E-4 );
   if (!ok) std::cout << "Error in test33 - " << f1.Eval(0) << "  " << f2.Eval(0)/f2.Integral(-10,10,-10,10) << std::endl;
   return ok;
}

bool test34() {
   // test for bug 8105
   bool ok  = true;
   TF1 f1("f1","(1.- gaus)*[3]",-10,10);
   f1.SetParameters(1,0,1,3);
   ok &=  TMath::AreEqualAbs( f1.Eval(1), (1.- TMath::Gaus(1,0,1) )*3., 1.E-10);
   return ok;

}
bool test35() {
   // test for similar pre-defined functions
   bool ok = true;
   TF1 f1("f1","cheb1(0)+cheb10(2)",-1,1);
   std::vector<double> par(13);
   par.assign(13,1.); par[1] = 2; par[2] = 3;
   TF1 g1("g1",[](double *x, double *p){ return ROOT::Math::ChebyshevN(1, x[0], p ) + ROOT::Math::ChebyshevN(10,x[0],p+2 ); }, -1, 1, 13);
   f1.SetParameters(par.data());
   g1.SetParameters(par.data());

   ok &=  TMath::AreEqualRel( f1.Eval(2), g1.Eval(2), 1.E-6);
   if (!ok) std::cout << "Error in test35 - f1 != g1 " << f1.Eval(2) << "  " << g1.Eval(2) << std::endl;

   TF1 f2("f2","cheb10(0)+cheb1(11)",-1,1);
   TF1 g2("g2",[](double *x, double *p){ return ROOT::Math::ChebyshevN(10, x[0], p ) + ROOT::Math::ChebyshevN(1,x[0],p+11 ); }, -1, 1, 13);
   f2.SetParameters(par.data());
   g2.SetParameters(par.data());

   ok &=  TMath::AreEqualRel( f2.Eval(2), g2.Eval(2), 1.E-6);
   if (!ok) std::cout << "Error in test35 - f2 != g2 " << f2.Eval(2.) << "  " << g2.Eval(2.) << std::endl;

   return ok;
}
bool test36() {
   // test for mixed dim functions
   bool ok = true;
   TF2 f1("f1","xygaus(0) + gaus(5)");
   f1.SetParameters(1,0,1,1,2,2,-1,1);
   auto g1 = [](double x, double y){ return TMath::Gaus(x,0,1)*TMath::Gaus(y,1,2)+2.*TMath::Gaus(x,-1,1); };
   ok &=  TMath::AreEqualAbs( f1.Eval(1,1), g1(1,1), 1.E-10);

   TF2 f2("f2","xygaus(0) + gaus[y](5)");
   f2.SetParameters(1,0,1,1,2,2,-1,1);
   auto g2 = [](double x, double y){ return TMath::Gaus(x,0,1)*TMath::Gaus(y,1,2)+2.*TMath::Gaus(y,-1,1); };
   ok &=  TMath::AreEqualAbs( f2.Eval(1,1), g2(1,1), 1.E-10);


   return ok;
}

bool test37() {
  // test for inserting correcting polynomials (bug ROOT-8496)
  bool ok = true;
  TF1 f1("f1","[0]*pol1(1) + pol2(3)*[6]",0,1);
  f1.SetParameters(2,1,2,1,2,3,4);
  auto ref = [](double x) { return 2 * (1 + 2*x ) +  (1 + 2*x + 3*x*x) * 4 ; };

  ok &= TMath::AreEqualAbs( f1.Eval(0.5), ref(0.5), 1.E-10);
  return ok;
}

bool test38() {
  // test for missing parameters  (bug ROOT-8182)
  bool ok = true;
  TF1 f1("f1","[1]",0,1);
  f1.SetParameters(999,2);
  ok &= (f1.Eval(0) == 2.);
  TF1 f2("f2","[A]+[1]*x",0,1);
  f2.SetParameters(999,2,3);
  ok &= (f2.Eval(2) == 7.);
  return ok;
}

bool test39() {
   // test special characters in parameter names (bug ROOT-8303)
   // test with operator ^, @ and predefined functions (pol, gaus, etc..)
   bool ok = true;

   TF1 f1("f1","[s^x]*x+[0]");
   f1.SetParameters(1,2);
   ok &= (f1.Eval(2) == 2*2+1);

   TF1 f2("f2","[0]*x+[s@x]");
   f2.SetParameters(2,1);
   ok &= (f2.Eval(2) == 2*2+1);

   TF1 f3("f2","[0]*x+[pol_par_1]");
   f3.SetParameters(2,1);
   ok &= (f3.Eval(2) == 2*2+1);

   TF1 f4("f2","gaus+[gaus_offset]*x");
   f4.SetParameters(2,2,1,3);
   ok &= (f4.Eval(2) == 2+3*2);

   return ok;
}

bool test40()
{
   // test parsing variables/parameters of user-defined functions

   TF2 f1("f1", "x - y", 0, 5, 0, 5);
   TF2 f2("f2", "f1(y,x)", 0, 5, 0, 5);
   bool ok = (f1.Eval(1, 2) == -1);
   ok &= (f2.Eval(1, 2) == 1);

   TF3 f3("f3", "x + 2*y + 3*z", 0, 5, 0, 5, 0, 5);
   TF1 f4("f4", "f3(x,x,x)", 0, 5);
   ok &= (f3.Eval(2, 2, 2) == 12);
   ok &= (f4.Eval(2) == 12);

   TF1 f5("f5", "[0]*x + [1]", 0, 5);
   TF1 f6("f6", "f5(x,[1],[0])", 0, 5);
   f6.SetParameters(1, 2);
   ok &= (f6.Eval(0) == 1);
   ok &= (f6.Eval(1) == 3);

   // implicit x now
   TF1 f7("f7", "f5([1], [0])", 0, 5);
   f7.SetParameters(1, 2);
   ok &= (f7.Eval(0) == 1);
   ok &= (f7.Eval(1) == 3);

   // now implicit parameters
   TF2 f8("f8", "f5(y)", 0, 5, 0, 5);
   f8.SetParameters(1, 2);
   ok &= (f8.Eval(0, 0) == 2);
   ok &= (f8.Eval(1, 0) == 2);
   ok &= (f8.Eval(0, 1) == 3);
   ok &= (f8.Eval(1, 1) == 3);

   // and test [p0] notation
   TF1 f9("f9", "[p0]*x + [p1]", 0, 5);
   TF1 f10("f10", "f9(x,[p1],[p0])", 0, 5);
   f10.SetParameters(1, 2);
   ok &= (f10.Eval(0) == 1);
   ok &= (f10.Eval(1) == 3);

   // implicit x now
   TF1 f11("f11", "f9([p1], [p0])", 0, 5);
   f11.SetParameters(1, 2);
   ok &= (f11.Eval(0) == 1);
   ok &= (f11.Eval(1) == 3);

   return ok;
}

bool test41()
{
   // Test variable/parameter parsing for parametrized functions

   bool ok = true;

   // old variable-counting method
   TF1 f1("f1", "gaus(0) + gaus(3)", -5, 5);
   f1.SetParameters(1, 0, 1, 1, 1, 1);
   ok &= fpEqual(f1.Eval(0), 1 + TMath::Exp(-.5), true);
   ok &= fpEqual(f1.Eval(1), 1 + TMath::Exp(-.5), true);

   // new param-range method
   TF1 f2("f2", "gaus([0..2]) + gaus([3..5])", -5, 5);
   f2.SetParameters(1, 0, 1, 1, 1, 1);
   ok &= fpEqual(f2.Eval(0), 1 + TMath::Exp(-.5), true);
   ok &= fpEqual(f2.Eval(1), 1 + TMath::Exp(-.5), true);

   TF1 f3("f3", "[0] + gaus([1..3])", -5, 5);
   f3.SetParameters(2, 1, 0, 1);
   ok &= fpEqual(f3.Eval(0), 3, true);
   ok &= fpEqual(f3.Eval(1), 2 + TMath::Exp(-.5), true);

   TF2 f4("f4", "gaus(y)", -5, 5, -5, 5);
   f4.SetParameters(2, 0, 1);
   ok &= fpEqual(f4.Eval(0, 0), 2, true);
   ok &= fpEqual(f4.Eval(1, 0), 2, true);
   ok &= fpEqual(f4.Eval(0, -1), 2 * TMath::Exp(-.5), true);
   ok &= fpEqual(f4.Eval(1, -1), 2 * TMath::Exp(-.5), true);

   TF2 f5("f5", "[0] + gaus(y, [1..3])", -5, 5, -5, 5);
   f5.SetParameters(0, 2, 0, 1);
   ok &= fpEqual(f5.Eval(0, 0), 2, true);
   ok &= fpEqual(f5.Eval(1, 0), 2, true);
   ok &= fpEqual(f5.Eval(0, -1), 2 * TMath::Exp(-.5), true);
   ok &= fpEqual(f5.Eval(1, -1), 2 * TMath::Exp(-.5), true);

   return ok;
}

bool test42()
{
   // Test variable parsing when using form x[N]

   bool ok = true;

   TF2 f1("f1", "x[1] + 1", -5, 5, -5, 5);
   ok &= (f1.Eval(1, 1) == 2);
   ok &= (f1.Eval(0, 1) == 2);
   ok &= (f1.Eval(1, 0) == 1);
   ok &= (f1.Eval(0, 0) == 1);

   TF2 f2("f2", "f1(y,x) + 0*y", -5, 5, -5, 5);
   ok &= (f2.Eval(1, 1) == 2);
   ok &= (f2.Eval(0, 1) == 1);
   ok &= (f2.Eval(1, 0) == 2);
   ok &= (f2.Eval(0, 0) == 1);

   TF2 f3("f3", "f1(x[1], x[0]) + 0*y", -5, 5, -5, 5);
   ok &= (f3.Eval(1, 1) == 2);
   ok &= (f3.Eval(0, 1) == 1);
   ok &= (f3.Eval(1, 0) == 2);
   ok &= (f3.Eval(0, 0) == 1);

   return ok;
}

bool test43()
{
   // test whether value of parameter name carries through

   bool ok = true;

   TF1 f1("f1", "[const] + [linear]*x", -5, 5);
   f1.SetParameters(1, 2);

   TF1 f2("f2", "f1", -5, 5);
   ok &= (f2.Eval(1) == 3);

   TF1 f3("f3", "f1(x, [const], [linear])", -5, 5);
   ok &= (f3.Eval(1) == 3);

   TF1 f4("f4", "f1([const], [linear])", -5, 5);
   ok &= (f4.Eval(1) == 3);

   TF1 f5("f5", "f1(x)", -5, 5);
   ok &= (f5.Eval(1) == 3);

   TF1 f6("f6", "f1([first], [second])");
   // parameters "should" initialize to zero
   ok &= (f6.Eval(1) == 0);

   return ok;
}

bool test44()
{
   // test whether user-defined and parametrized functions can be nested

   bool ok = true;

   TF1 f1("f1", "x**[0]");
   TF1 f2("f2", "x + 1");
   TF2 f3("f3", "f1(f2(x), y)");
   ok &= (f3.Eval(2, 3) == 27);

   TF1 f4("f4", "f2(f2(x))");
   ok &= (f4.Eval(5) == 7);

   TF1 f5("f5", "gaus(f2(x), 1, 0, 1)");
   ok &= fpEqual(f5.Eval(0), TMath::Exp(-.5), true);

   TF1 f6("f6", "gaus(gaus(x, 1, 0, 1), 1, 0, 1)");
   ok &= fpEqual(f6.Eval(0), TMath::Exp(-.5), true);

   return ok;
}

bool test45()
{
   // test dealing with whitespace in parameter names
   // inlcuding cloning tests (see ROOT-8971)
   TF1* func = new TF1("expo","expo");
   func->SetParNames("A", "- 1 / T");
   func->SetParameters(1,1);

   TF1 * func2 = (TF1*) func->Clone("func2");

   bool ok = fpEqual( func2->Eval(2), func->Eval(2), true);
   return ok;
}

bool test46() {
   // test multi-dim formula (like new xyzgaus)
   auto func = new TF3("f3","xyzgaus");
   func->SetParameters(2,1,2,3,4,5,6);
   bool ok = fpEqual( func->Eval(2,4,6), 2.*TMath::Gaus(2,1,2)*TMath::Gaus(4,3,4)*TMath::Gaus(6,5,6) , true);

   auto func2 = new TF3("f3","gaus(x,[0],[1],[2])*gaus(y,1,[3],[4])*gaus(z,1,[5],[6])");
   double x[] = {2,4,6};
   ok &= fpEqual( func->EvalPar(x,nullptr), func2->EvalPar(x, func->GetParameters() ), true );
   return ok;
}

bool test47() {
   // test mod operator
   // one needs to convert always to integer because % works only for int
   auto f1 = new TF1("f1","exp(x)");
   (void)f1; // f1 is used by modf but the compiler doesn't see that.
   auto func = new TF1("modf","int(2*f1(x)) % 3");
   bool ok = func->Eval(1) == 2;
   ok &= func->Eval(3) == 1;
   ok &= func->Eval(1.2) == 0;
   return ok;
}

bool test48() {
   // test creating two identical functions
   // and reading back from a file
   // ROOT-9467
   // The bug woruld need to exit ROOT and when the file already esists
   TString fname = "TFormulaTest48.root";
   int prevErr = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kFatal; 
   TFile* f = TFile::Open(fname);
   gErrorIgnoreLevel = prevErr; 
   if (!f) {
      TFile * fout = TFile::Open(fname,"NEW");
      TF1* f1 = new TF1("f1", "[0] + [1]*x+2.0", 0, 1);
      TF1* f2 = new TF1("f2", "[0] + [1]*x+2.0", 0, 1);
      f1->SetParameters(1,1);
      f2->SetParameters(0,2);

      f1->Write();
      f2->Write();
      fout->Close();
      f = TFile::Open(fname);
   }

   TF1* f1 = dynamic_cast<TF1*>(f->Get("f1"));
   TF1* f2 = dynamic_cast<TF1*>(f->Get("f2"));

   bool ok = f1 != nullptr && f2 != nullptr;
   if (ok) {
      ok &= (f1->Eval(1) == 4. && f1->Eval(1) == f2->Eval(1) );
   }
   return ok;
}

bool test49() {
   // test copy consttructor in case of lazy initialization (i.e. when reading from a file)
   TFile* f = TFile::Open("TFormulaTest49.root","RECREATE");   
   if (!f) {
      Error("test49","Error creating file for test49");
      return false;
   }
   TF1 f1("f1","x*[0]");
   f1.SetParameter(0,2); 
   f1.Write();
   f->Close();
   // read the file 
   f = TFile::Open("TFormulaTest49.root");   
   if (!f) {
      Error("test49","Error reading file for test49");
      return false;
   }
   auto fr = (TF1*) f->Get("f1");
   if (!fr) { 
      Error("test49","Error reading function from file for test49");
      return false;
   }
   // create a copy
   TF1 fr2 = *fr;
   bool ok = (fr->Eval(2.) == 4.);
   ok &= ( fr2.Eval(2.) == fr->Eval(2.) );

   // now read using an indpendent process (ROOT session)
   // this should cause ROOT-9801

   TMacro m;
   m.AddLine("bool TFormulaTest49() { TFile * f = TFile::Open(\"TFormulaTest49.root\");"
             "TF1 *f1 = (TF1*) f->Get(\"f1\"); TF1 f2 = *f1;"
             "bool ok = (f1->Eval(2) == f2.Eval(2.)) && (f2.Eval(2.) == 4.);"
             "if (!ok) Error(\"test49\",\"Error in test49 (lazy initialization)\");" 
             "return ok; }");

   m.SaveSource("TFormulaTest49.C");
   int ret = gSystem->Exec("root.exe -q -l -i TFormulaTest49.C");
   ok |= (ret == 0);
   return ok;
}

bool test50() {
   // test detailed printing of function
   TFormula f1("f1","[A]*sin([B]*x)");
   f1.Print("V");
   bool ok = f1.IsValid(); 

   TF2 f2("f2","[0]*x+[1]*y");
   f2.Print("V");
   ok &= f2.GetFormula()->IsValid();

   // create using lambda expression, need to pass ndim and npar
   TFormula f3("f3","[](double *x, double *p){ return p[0]*x[0] + p[1]; } ",1,2);
   f3.Print("V");
   ok &= f3.IsValid();

   // create again using lambda from TF1, need to pass xmin(0.),xmax(1.), npar (1)
   TF1 f4("f3","[](double *x, double *p){ return p[0]*x[0]; } ",0.,1.,1);  
   f4.Print("V");
   ok &= f3.IsValid();

   return ok;
}

bool test51() {
   //switch off error messages to have test passing
   int prevLevel = gErrorIgnoreLevel; 
   gErrorIgnoreLevel= kFatal; 
   TFormula f("fMissingParenthesis", "exp(x");
   bool ok = !f.IsValid();
   TFormula f2("fEmpty", "");
   ok &= !f2.IsValid();
   TFormula f3("fNonsense", "skmg#$#@!1");
   ok &= !f3.IsValid();
   gErrorIgnoreLevel = prevLevel; 
   return ok;
}

bool test52() {
   // test for bug 10815
   // mixing user previous defined functions (available in gROOT)
   // and pre-defined functions
   bool ok  = true;
   TF1 f1("f1gaus","[0]*gaus(1)",-10,10);
   TF1 f2("f2","[0]*f1gaus",-10,10);
   f1.SetParameters(2,3,1,2);
   f2.SetParameters(3,2,3,1,2);
   ok &=  TMath::AreEqualAbs( f1.Eval(1), 2.*3.*TMath::Gaus(1,1,2), 1.E-10);
   if (!ok) Error("test52","Error testing f1");
   bool ret =  TMath::AreEqualAbs( f2.Eval(1), 3.*2.*3.*TMath::Gaus(1,1,2), 1.E-10);
   if (!ret) Error("test52","Error testing f2");
   ok &= ret; 
   TF1 f3("f3","f1gaus*gaus(4)",-10,10);
   f3.SetParameters(2,3,1,2,3,2,3);
   ret =  TMath::AreEqualAbs( f3.Eval(1), 2.*3.*TMath::Gaus(1,1,2) * 3. * TMath::Gaus(1,2,3), 1.E-10);    
   if (!ret) Error("test52","Error testing f3");
   ok &= ret; 
   // check also after
   TF1 f4("gaus2a","[0]*gaus(1)",-10,10);
   TF1 f5("f2","[0]*gaus2a",-10,10);
   f4.SetParameters(2,3,1,2);
   f5.SetParameters(3,2,3,1,2);
   ret =  TMath::AreEqualAbs( f5.Eval(1), 3.*f4.Eval(1),1.E-10);
   if (!ret) Error("test52","Error testing f4 & f5");
   return ok;

}


///////////////////////////////////////////////////////////////////////////////////////

void PrintError(int itest)  {
   Error("TFormula test","test%d FAILED ",itest);
   failedTests.push_back(itest);
}
void IncrTest(int & itest) {
   if (itest > 0) std::cout << ".\n";
   itest++;
   std::cout << "Test " << itest << " :        ";
}

int runTests(bool debug = false) {

   verbose = debug;

   int itest = 0;


   IncrTest(itest); if (!test1() ) { PrintError(itest); }
   IncrTest(itest); if (!test2() ) { PrintError(itest); }
   IncrTest(itest); if (!test3() ) { PrintError(itest); }
   IncrTest(itest); if (!test4() ) { PrintError(itest); }
   IncrTest(itest); if (!test5() ) { PrintError(itest); }
   IncrTest(itest); if (!test6() ) { PrintError(itest); }
   IncrTest(itest); if (!test7() ) { PrintError(itest); }
   IncrTest(itest); if (!test8() ) { PrintError(itest); }
   IncrTest(itest); if (!test9() ) { PrintError(itest); }
   IncrTest(itest); if (!test10() ) { PrintError(itest); }
   IncrTest(itest); if (!test11() ) { PrintError(itest); }
   IncrTest(itest); if (!test12() ) { PrintError(itest); }
   IncrTest(itest); if (!test13() ) { PrintError(itest); }
   IncrTest(itest); if (!test14() ) { PrintError(itest); }
   IncrTest(itest); if (!test15() ) { PrintError(itest); }
   IncrTest(itest); if (!test16() ) { PrintError(itest); }
   IncrTest(itest); if (!test17() ) { PrintError(itest); }
   IncrTest(itest); if (!test18() ) { PrintError(itest); }
   IncrTest(itest); if (!test19() ) { PrintError(itest); }
   IncrTest(itest); if (!test20() ) { PrintError(itest); }
   IncrTest(itest); if (!test21() ) { PrintError(itest); }
   IncrTest(itest); if (!test22() ) { PrintError(itest); }
   IncrTest(itest); if (!test23() ) { PrintError(itest); }
   IncrTest(itest); if (!test24() ) { PrintError(itest); }
   IncrTest(itest); if (!test25() ) { PrintError(itest); }
   IncrTest(itest); if (!test26() ) { PrintError(itest); }
   IncrTest(itest); if (!test27() ) { PrintError(itest); }
   IncrTest(itest); if (!test28() ) { PrintError(itest); }
   IncrTest(itest); if (!test29() ) { PrintError(itest); }
   IncrTest(itest); if (!test30() ) { PrintError(itest); }
   IncrTest(itest); if (!test31() ) { PrintError(itest); }
   IncrTest(itest); if (!test32() ) { PrintError(itest); }
   IncrTest(itest); if (!test33() ) { PrintError(itest); }
   IncrTest(itest); if (!test34() ) { PrintError(itest); }
   IncrTest(itest); if (!test35() ) { PrintError(itest); }
   IncrTest(itest); if (!test36() ) { PrintError(itest); }
   IncrTest(itest); if (!test37() ) { PrintError(itest); }
   IncrTest(itest); if (!test38() ) { PrintError(itest); }
   IncrTest(itest); if (!test39() ) { PrintError(itest); }
   IncrTest(itest); if (!test40() ) { PrintError(itest); }
   IncrTest(itest); if (!test41() ) { PrintError(itest); }
   IncrTest(itest); if (!test42() ) { PrintError(itest); }
   IncrTest(itest); if (!test43() ) { PrintError(itest); }
   IncrTest(itest); if (!test44() ) { PrintError(itest); }
   IncrTest(itest); if (!test45() ) { PrintError(itest); }
   IncrTest(itest); if (!test46() ) { PrintError(itest); }
   IncrTest(itest); if (!test47() ) { PrintError(itest); }
   IncrTest(itest); if (!test48() ) { PrintError(itest); }
   IncrTest(itest); if (!test49() ) { PrintError(itest); }
   IncrTest(itest); if (!test50() ) { PrintError(itest); }
   IncrTest(itest); if (!test51() ) { PrintError(itest); }
   IncrTest(itest); if (!test52() ) { PrintError(itest); }

   std::cout << ".\n";

   if (failedTests.size() == 0)
      std::cout << "All TFormula Parsing tests PASSED !" << std::endl;
   else {
      Error("TFORMULA Tests","%d tests failed ",int(failedTests.size()) );
      std::cout << "failed tests are : ";
      for (auto & ittest : failedTests) {
         std::cout << ittest << "   ";
      }
      std::cout << std::endl;
   }

   return failedTests.size();

}

};
