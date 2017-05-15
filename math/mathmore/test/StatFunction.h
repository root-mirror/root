// @(#)root/test:$Id$
// Author: Lorenzo Moneta       06/2005
//         Martin Storø Nyfløtt 05/2017

#ifndef STATFUNCTION
#define STATFUNCTION

#include "Math/AllIntegrationTypes.h"
#include "Math/IParamFunction.h"
#include "Math/RootFinderAlgorithms.h"
#include "Math/RootFinder.h"

#include <iomanip>

using namespace ROOT::Math;

#define INF std::numeric_limits<double>::infinity()

// typedef for a free function like gamma(double x, double a, double b)

typedef std::function<double(double, double, double)> AlgoFunc;

// statistical function class
const int N_PAR = 2;

class StatFunction : public ROOT::Math::IParamFunction {
private:
   double DoEvalPar(double x, const double *) const
   {
      // use explicitly cached param values
      return fPdf(x, fParams[0], fParams[1]);
   }

   AlgoFunc fPdf;
   AlgoFunc fCdf;
   AlgoFunc fQuant;
   double fParams[N_PAR];
   double fScaleIg;
   double fScaleDer;
   double fScaleInv;
   int NFuncTest;
   double xmin;
   double xmax;
   double xlow;
   double xup;
   bool fHasLowRange;
   bool fHasUpRange;
   double fStartRoot;

public:
   StatFunction(AlgoFunc pdf, AlgoFunc cdf, AlgoFunc quant, double x1 = -INF, double x2 = INF)
      : fPdf(pdf), fCdf(cdf), fQuant(quant), xmin(0.), xmax(0.), xlow(x1), xup(x2), fHasLowRange(false),
        fHasUpRange(false), fStartRoot(0.)
   {
      fScaleIg = 10;   // scale for integral test
      fScaleDer = 1;   // scale for der test
      fScaleInv = 100; // scale for inverse test
      for (int i = 0; i < N_PAR; ++i) fParams[i] = 0;
      NFuncTest = 100;
      if (xlow > -INF) fHasLowRange = true;
      if (xup < INF) fHasUpRange = true;
   }

   unsigned int NPar() const { return N_PAR; }
   const double *Parameters() const { return fParams; }
   ROOT::Math::IGenFunction *Clone() const { return new StatFunction(fPdf, fCdf, fQuant); }

   void SetParameters(const double *p) { std::copy(p, p + N_PAR, fParams); }

   void SetParameters(double p0, double p1)
   {
      fParams[0] = p0;
      fParams[1] = p1;
   }

   void SetTestRange(double x1, double x2)
   {
      xmin = x1;
      xmax = x2;
   }
   void SetNTest(int n) { NFuncTest = n; }
   void SetStartRoot(double x) { fStartRoot = x; }

   double Pdf(double x) const { return (*this)(x); }

   double Cdf(double x) const { return fCdf(x, fParams[0], fParams[1]); }

   double Quantile(double x) const { return fQuant(x, fParams[0], fParams[1]); }

   // test integral with cdf function
   void TestIntegral(IntegrationOneDim::Type algotype);

   // test derivative from cdf to pdf function
   void TestDerivative();

   // test root finding algorithm for finding inverse of cdf
   void TestInverse1(RootFinder::EType algotype);

   // test root finding algorithm for finding inverse of cdf using derivatives
   void TestInverse2(RootFinder::EType algotype);

   void SetScaleIg(double s) { fScaleIg = s; }
   void SetScaleDer(double s) { fScaleDer = s; }
   void SetScaleInv(double s) { fScaleInv = s; }
};

#endif // STATFUNCTION
