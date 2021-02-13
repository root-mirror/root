#include "TMath.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TROOT.h"

#include "Fit/UnBinData.h"
#include "Fit/Fitter.h"

#include "Math/IParamFunction.h"
#include "Math/WrappedTF1.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/WrappedParamFunction.h"
#include "Math/MultiDimParamFunctionAdapter.h"

#include "TGraphErrors.h"

#include "TStyle.h"

#include "Math/DistFunc.h"

#include <string>
#include <iostream>

#include "TStopwatch.h"

#define DEBUG

// test fit with many dimension

const int N = 10;
const std::string branchType = "x[10]/D";
const int NPoints = 100000;

// const int N = 50;
// const std::string branchType = "x[50]/D";
// const int NPoints = 10000;

double truePar[2 * N];
double iniPar[2 * N];
// const int nfit = 1;
const int strategy = 0;

double gausnorm(const double *x, const double *p)
{

   double invsig = 1. / p[1];
   double tmp = (x[0] - p[0]) * invsig;
   const double sqrt_2pi = 1. / std::sqrt(2. * 3.14159);
   return std::exp(-0.5 * tmp * tmp) * sqrt_2pi * invsig;
}

double gausnormN(const double *x, const double *p)
{
   double f = 1;
   for (int i = 0; i < N; ++i)
      f *= gausnorm(x + i, p + 2 * i);

   return f;
}

struct MINUIT2 {
   static std::string name() { return "Minuit2"; }
   static std::string name2() { return ""; }
};

// fill fit data structure
ROOT::Fit::UnBinData *FillUnBinData(TTree *tree)
{

   // fill the unbin data set from a TTree
   ROOT::Fit::UnBinData *d = 0;

   // for the large tree access directly the branches
   d = new ROOT::Fit::UnBinData();

   unsigned int n = tree->GetEntries();
#ifdef DEBUG
   std::cout << "number of unbin data is " << n << " of dim " << N << std::endl;
#endif
   d->Initialize(n, N);
   TBranch *bx = tree->GetBranch("x");
   double vx[N];
   bx->SetAddress(vx);
   std::vector<double> m(N);
   for (int unsigned i = 0; i < n; ++i) {
      bx->GetEntry(i);
      d->Add(vx);
      for (int j = 0; j < N; ++j)
         m[j] += vx[j];
   }

#ifdef DEBUG
   std::cout << "average values of means :\n";
   for (int j = 0; j < N; ++j)
      std::cout << m[j] / n << "  ";
   std::cout << "\n";
#endif

   delete tree;
   tree = 0;
   return d;
}

// unbin fit

typedef ROOT::Math::IParamMultiFunction Func;
template <class MinType, class T>
int DoUnBinFit(T *tree, Func &func, bool debug = false)
{

   ROOT::Fit::UnBinData *d = FillUnBinData(tree);
   // need to have done Tree->Draw() before fit
   // FillUnBinData(d,tree);

   // std::cout << "data size type and size  is " << typeid(*d).name() <<  "   " << d->Size() << std::endl;
   std::cout << "Fit data size =  " << d->Size() << " dimension = " << d->NDim() << std::endl;

   // printData(d);

   // create the fitter
   // std::cout << "Fit parameter 2  " << f.Parameters()[2] << std::endl;

   ROOT::Fit::Fitter fitter;
   fitter.Config().SetMinimizer(MinType::name().c_str(), MinType::name2().c_str());

   if (debug)
      fitter.Config().MinimizerOptions().SetPrintLevel(3);
   else
      fitter.Config().MinimizerOptions().SetPrintLevel(0);

   // set tolerance 1 for tree to be same as in TTTreePlayer::UnBinFIt
   fitter.Config().MinimizerOptions().SetTolerance(1);

   // set strategy (0 to avoid MnHesse
   fitter.Config().MinimizerOptions().SetStrategy(strategy);

   // create the function

   fitter.SetFunction(func);
   // need to fix param 0 , normalization in the unbinned fits
   // fitter.Config().ParSettings(0).Fix();

   bool ret = fitter.Fit(*d);
   if (!ret) {
      std::cout << " Fit Failed " << std::endl;
      return -1;
   }
   if (debug)
      fitter.Result().Print(std::cout);

   // check fit result
   double chi2 = 0;
   for (int i = 0; i < N; ++i) {
      double d = (truePar[i] - fitter.Result().Value(i)) / (fitter.Result().Error(i));
      chi2 += d * d;
   }
   double prob = ROOT::Math::chisquared_cdf_c(chi2, N);
   int iret = (prob < 1.0E-6) ? -1 : 0;
   if (iret != 0) {
      std::cout << "Found difference in fitted values - prob = " << prob << std::endl;
      if (!debug)
         fitter.Result().Print(std::cout);
      for (int i = 0; i < N; ++i) {
         double d = (truePar[i] - fitter.Result().Value(i)) / (fitter.Result().Error(i));
         std::cout << "par_" << i << " = " << fitter.Result().Value(i) << " true  = " << truePar[i] << " pull = " << d
                   << std::endl;
      }
   }

   delete d;

   return iret;
}

template <class MinType>
int DoFit(TTree *tree, Func &func, bool debug = false)
{
   return DoUnBinFit<MinType, TTree>(tree, func, debug);
}
// template <class MinType>
// int DoFit(TH1 * h1, Func & func, bool debug = false, bool copyData = false ) {
//    return DoBinFit<MinType, TH1>(h1, func, debug, copyData);
// }
// template <class MinType>
// int DoFit(TGraph * gr, Func & func, bool debug = false, bool copyData = false ) {
//    return DoBinFit<MinType, TGraph>(gr, func, debug, copyData);
// }

template <class MinType, class FitObj>
int FitUsingNewFitter(FitObj *fitobj, Func &func)
{

   std::cout << "\n************************************************************\n";
   std::cout << "\tFit using new Fit::Fitter  " << typeid(*fitobj).name() << std::endl;
   std::cout << "\tMinimizer is " << MinType::name() << "  " << MinType::name2() << " func dim = " << func.NDim()
             << std::endl;

   int iret = 0;
   TStopwatch w;
   w.Start();

#ifdef DEBUG
   std::cout << "initial Parameters " << iniPar << "  " << *iniPar << "   " << *(iniPar + 1) << std::endl;
   func.SetParameters(iniPar);
   iret |= DoFit<MinType>(fitobj, func, true);

#else
   for (int i = 0; i < nfit; ++i) {
      func.SetParameters(iniPar);
      iret = DoFit<MinType>(fitobj, func, false);
      if (iret != 0) {
         std::cout << "Fit failed " << std::endl;
         break;
      }
   }
#endif
   w.Stop();
   std::cout << "\nTime: \t" << w.RealTime() << " , " << w.CpuTime() << std::endl;
   std::cout << "\n************************************************************\n";

   return iret;
}

int testNdimFit()
{

   std::cout << "\n\n************************************************************\n";
   std::cout << "\t UNBINNED TREE (GAUSSIAN MULTI-DIM)  FIT\n";
   std::cout << "************************************************************\n";

   TTree *t1 = new TTree("t2", "a large Tree with gaussian variables");
   double x[N];
   Int_t ev;
   t1->Branch("x", x, branchType.c_str());
   t1->Branch("ev", &ev, "ev/I");

   // generate the true parameters
   for (int j = 0; j < N; ++j) {
      double mu = double(j) / 10.;
      double s = 1.0 + double(j) / 10.;
      truePar[2 * j] = mu;
      truePar[2 * j + 1] = s;
   }

   // fill the tree
   TRandom3 r;
   for (Int_t i = 0; i < NPoints; i++) {
      for (int j = 0; j < N; ++j) {
         double mu = truePar[2 * j];
         double s = truePar[2 * j + 1];
         x[j] = r.Gaus(mu, s);
      }

      ev = i;
      t1->Fill();
   }
   // t1.Draw("x"); // to select fit variable

   for (int i = 0; i < N; ++i) {
      iniPar[2 * i] = 0;
      iniPar[2 * i + 1] = 1;
   }

   // use simply TF1 wrapper
   // ROOT::Math::WrappedMultiTF1 f2(*f1);
   ROOT::Math::WrappedParamFunction<> f2(&gausnormN, N, 2 * N, iniPar);

   int iret = 0;
   iret |= FitUsingNewFitter<MINUIT2>(t1, f2);

   return iret;
}

int main()
{
   return testNdimFit();
}
