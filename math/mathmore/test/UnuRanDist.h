// class to test unuran dist

#include "TUnuran.h"
#include "Math/Util.h"

#include <iostream>
#include <cassert>
#include <string>

class UnuRanDist {
public:


   UnuRanDist() : first(true) {
      // use default grandom
   }


      //~UnuRanDist() { if (fUnuran) delete fUnuran; }

   void SetSeed(int seed) {
      gRandom->SetSeed(seed);
      first = true;
   }

   double Gaus(double mu, double sigma) {
      if (first) {
         //if (!fUnuran.Init("normal()", "method=cstd;variant=0") ) { // acr method (default)
         if (!fUnuran.Init("normal()", "method=cstd;variant=6") ) {
            assert(0);
         }
         first = false;
      }
      return mu + sigma* fUnuran.Sample();
   }

   int Poisson(double mu) {
      if (first) {
//          std::string smu = ROOT::Math::Util::ToString(mu);
//          std::string dist = "poisson(" + smu + ")";
//          if (!fUnuran.Init(dist, "method=dstd") ) {
//             assert(0);
//          }
          if (!fUnuran.InitPoisson(mu,"method=dstd") ) {
             assert(0);
          }

         first = false;
      }
      else {
         par[0] = mu;
         fUnuran.ReInitDiscrDist(1,par);
      }
      return fUnuran.SampleDiscr();
   }

   int Binomial(int n, double p) {
      if (first) {
//          std::string sn = ROOT::Math::Util::ToString(n);
//          std::string sp = ROOT::Math::Util::ToString(p);
//          std::string dist = "binomial(" + sn + "," + sp + ")";
//          std::cout << dist << std::endl;
         if (!fUnuran.InitBinomial(n,p,"method=dstd") ) {
             assert(0);
          }
         first = false;
      }
      else {
         par[0] = n;
         par[1] = p;
         fUnuran.ReInitDiscrDist(2,par);
      }
      return fUnuran.SampleDiscr();
   }



private:
   TUnuran fUnuran;
   bool  first;
   double par[2];
};
