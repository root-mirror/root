// test fitting using also RooFit and new Fitter
#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooGaussian.h>

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooGlobalFunc.h"
#include "RooFitResult.h"
#include "RooProdPdf.h"

#include <TF1.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TBranch.h>

#include "TStopwatch.h"

#include "Math/DistFunc.h"

#include "Fit/UnBinData.h"
#include "Fit/BinData.h"
#include "Fit/Fitter.h"



#include "WrapperRooPdf.h"

#include <string>
#include <iostream>
#include <vector>
#include <memory>

#include "MinimizerTypes.h"

#include "Math/WrappedParamFunction.h"
#include <cmath>

const int N = 3; // n must be greater than 1
const int nfit = 1;
const int nEvents = 10000;
double iniPar[2*N];

//#define DEBUG

typedef ROOT::Math::IParamMultiFunction Func;

void fillTree(TTree & t2) {


   double  x[N];
   Int_t ev;
   for (int j = 0; j < N; ++j) {
      std::string xname = "x_" + ROOT::Math::Util::ToString(j);
      std::string xname2 = "x_" + ROOT::Math::Util::ToString(j) + "/D";
      t2.Branch(xname.c_str(),&x[j],xname2.c_str());
   }
   t2.Branch("ev",&ev,"ev/I");
   //fill the tree
   TRandom3 r;
   for (Int_t i=0;i<nEvents;i++) {
      for (int j = 0;  j < N; ++j) {
         double mu = double(j)/10.;
         double s  = 1.0 + double(j)/10.;
         x[j] = r.Gaus(mu,s);
      }

      ev = i;
      t2.Fill();
   }
   // t2.Print();
   // std::cout << "number of branches " << t2.GetNbranches() << std::endl;

   t2.ResetBranchAddresses() ;
}

void FillUnBinData(ROOT::Fit::UnBinData &d, TTree * tree ) {

   // fill the unbin data set from a TTree

   // large tree
   unsigned int n = tree->GetEntries();
   std::cout << "number of unbin data is " << n << " of dim " << N << std::endl;
   d.Append(n,N);

   double vx[N];
   for (int j = 0; j <N; ++j) {
      std::string bname = "x_" + ROOT::Math::Util::ToString(j);
      TBranch * bx = tree->GetBranch(bname.c_str());
      bx->SetAddress(&vx[j]);
   }

   std::vector<double>  m(N);
   for (int unsigned i = 0; i < n; ++i) {
      tree->GetEntry(i);
      d.Add(vx);
      for (int j = 0; j < N; ++j)
         m[j] += vx[j];
   }

#ifdef DEBUG
   std::cout << "average values of means :\n";
   for (int j = 0; j < N; ++j)
      std::cout << m[j]/n << "  ";
   std::cout << "\n";
#endif

   tree->ResetBranchAddresses() ;

   return;

}



// class describing product of gaussian pdf
class MultiGaussRooPdf {

 public:

   // create all pdf
   MultiGaussRooPdf(unsigned int n) :
      x(n), m(n), s(n), g(n), pdf(n)
      {
         //assert(n >= 2);

         // create the gaussians
         for (unsigned int j = 0; j < n; ++j) {

            std::string xname = "x_" + ROOT::Math::Util::ToString(j);
            x[j] = new RooRealVar(xname.c_str(),xname.c_str(),-10000,10000) ;

            std::string mname = "m_" + ROOT::Math::Util::ToString(j);
            std::string sname = "s_" + ROOT::Math::Util::ToString(j);


//             m[j] = new RooRealVar(mname.c_str(),mname.c_str(),iniPar[2*j],-10000,10000) ;
//             s[j] = new RooRealVar(sname.c_str(),sname.c_str(),iniPar[2*j+1],-10000,10000) ;
            m[j] = new RooRealVar(mname.c_str(),mname.c_str(),iniPar[2*j],-10,10) ;
            s[j] = new RooRealVar(sname.c_str(),sname.c_str(),iniPar[2*j+1],-10,10) ;

            std::string gname = "g_" + ROOT::Math::Util::ToString(j);
            g[j] = new RooGaussian(gname.c_str(),"gauss(x,mean,sigma)",*x[j],*m[j],*s[j]);

            std::string pname = "prod_" + ROOT::Math::Util::ToString(j);
            if (j == 0)
               pdf[0] = g[0];
            else if (j == 1) {
               pdf[1] = new RooProdPdf(pname.c_str(),pname.c_str(),RooArgSet(*g[1],*g[0]) );
            }
            else if (j > 1) {
               pdf[j] = new RooProdPdf(pname.c_str(),pname.c_str(),RooArgSet(*g[j],*pdf[j-1]) );
            }
         }


   }

   RooAbsPdf & getPdf() { return *pdf.back(); }

   std::unique_ptr<RooArgSet>  getVars() {
      std::unique_ptr<RooArgSet> vars(new RooArgSet() );
      for (unsigned int i = 0; i < x.size(); ++i)
         vars->add(*x[i]);
      return vars;
   }

   ~MultiGaussRooPdf() {
      // free
      int n = x.size();
      for (int j = 0; j < n; ++j) {
         delete x[j];
         delete m[j];
         delete s[j];
         delete g[j];
         if (j> 0) delete pdf[j]; // no pdf allocated for j = 0
      }
   }

   private:

      std::vector<RooRealVar *> x;
      std::vector<RooRealVar *> m;
      std::vector<RooRealVar *> s;

      std::vector<RooAbsPdf *> g;
      std::vector<RooAbsPdf *> pdf;

};


//unbinned roo fit (large tree)
int  FitUsingRooFit(TTree & tree, RooAbsPdf & pdf, RooArgSet & xvars) {

   int iret = 0;
   std::cout << "\n************************************************************\n";
   std::cout << "\tFit using RooFit (Likelihood Fit) on " << pdf.GetName() << std::endl;



   TStopwatch w;
   w.Start();

   for (int i = 0; i < nfit; ++i) {

      RooDataSet data("unbindata","unbin dataset with x",&tree,xvars) ;



#ifdef DEBUG
      int level = 2;
      std::cout << "num entries = " << data.numEntries() << std::endl;
      bool save = true;
      pdf.getVariables()->Print("v"); // print the parameters
      std::cout << "\n\nDo the fit now \n\n";
#else
      int level = -1;
      bool save = false;
#endif

#ifndef _WIN32 // until a bug 30762 is fixed
      RooFitResult * result = pdf.fitTo(data, RooFit::Minos(0), RooFit::Hesse(0) , RooFit::PrintLevel(level), RooFit::Save(save) );
#else
      RooFitResult * result = pdf.fitTo(data );
#endif

#ifdef DEBUG
      assert(result != 0);
      std::cout << " Roofit status " << result->status() << std::endl;
      result->Print();
#endif
      iret |= (result == 0);

   }

   w.Stop();

   std::cout << "RooFit result " << std::endl;
   RooArgSet * params = pdf.getParameters(xvars);
   params->Print("v");
   delete params;


   std::cout << "\nTime: \t" << w.RealTime() << " , " << w.CpuTime() << std::endl;
   std::cout << "\n************************************************************\n";
   return iret;
}


// unbin fit
template <class MinType>
int DoFit(TTree * tree, Func & func, bool debug = false, bool = false ) {

   ROOT::Fit::UnBinData d;
   // need to have done Tree->Draw() before fit
   FillUnBinData(d,tree);

   //std::cout << "Fit parameter 2  " << f.Parameters()[2] << std::endl;

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(0);
   fitter.Config().SetMinimizer(MinType::name().c_str(),MinType::name2().c_str());
   fitter.Config().MinimizerOptions().SetTolerance(1.); // to be consistent with RooFit

   if (debug)
      fitter.Config().MinimizerOptions().SetPrintLevel(3);


   // create the function

   fitter.SetFunction(func);
   // need to fix param 0 , normalization in the unbinned fits
   //fitter.Config().ParSettings(0).Fix();

   bool ret = fitter.Fit(d);
   if (!ret) {
      std::cout << " Fit Failed " << std::endl;
      return -1;
   }
   fitter.Result().Print(std::cout);
   return 0;

}

template <class MinType, class FitObj>
int FitUsingNewFitter(FitObj * fitobj, Func & func, bool useGrad=false) {

   std::cout << "\n************************************************************\n";
   std::cout << "\tFit using new Fit::Fitter\n";
   std::cout << "\tMinimizer is " << MinType::name() << "  " << MinType::name2() << std::endl;

   int iret = 0;
   TStopwatch w; w.Start();

#ifdef DEBUG
   func.SetParameters(iniPar);
   iret |= DoFit<MinType>(fitobj,func,true, useGrad);

#else
   for (int i = 0; i < nfit; ++i) {
      func.SetParameters(iniPar);
      iret = DoFit<MinType>(fitobj,func, false, useGrad);
      if (iret != 0) break;
   }
#endif
   w.Stop();
   std::cout << "\nTime: \t" << w.RealTime() << " , " << w.CpuTime() << std::endl;
   std::cout << "\n************************************************************\n";

   return iret;
}

template <int N>
struct GausNorm {
   static double F(const double *x, const double *p) {
      return ROOT::Math::normal_pdf(x[N-1], p[2*N-1], p[2*N-2] ) * GausNorm<N-1>::F(x,p);
   }
};
template <>
struct GausNorm<1> {

   static double F(const double *x, const double *p) {
      return ROOT::Math::normal_pdf(x[0], p[1], p[0] );
   }
};


// double gausnorm(
//    return ROOT::Math::normal_pdf(x[0], p[1], p[0] );
//    // //return p[0]*TMath::Gaus(x[0],p[1],p[2]);
//    // double invsig = 1./p[1];
//    // double tmp = (x[0]-p[0]) * invsig;
//    // const double sqrt_2pi = 1./std::sqrt(2.* 3.14159 );
//    // return std::exp(-0.5 * tmp*tmp ) * sqrt_2pi * invsig;
// }



int main() {

   TTree tree("t","a large Tree with many gaussian variables");
   fillTree(tree);

   // set initial parameters
   for (int i = 0; i < N; ++i) {
      iniPar[2*i] = 1;
      iniPar[2*i+1] = 1;
   }

   // start counting the time
   MultiGaussRooPdf multipdf(N);
   RooAbsPdf & pdf = multipdf.getPdf();

   std::unique_ptr<RooArgSet> xvars(multipdf.getVars());

   WrapperRooPdf  wpdf( &pdf, *xvars );


   std::cout << "ndim " << wpdf.NDim() << std::endl;
   std::cout << "npar " << wpdf.NPar() << std::endl;
   for (unsigned int i = 0; i < wpdf.NPar(); ++i)
      std::cout << " par " << i << " is " <<  wpdf.ParameterName(i) << " value " << wpdf.Parameters()[i] << std::endl;


   FitUsingNewFitter<TMINUIT>(&tree,wpdf);

   // reset pdf original values
   wpdf.SetParameters(iniPar);

   FitUsingRooFit(tree,pdf, *xvars);

   // in case of N = 1 do also a simple gauss fit
   // using TF1 gausN
//   if (N == 1) {
   ROOT::Math::WrappedParamFunction<> gausn(&GausNorm<N>::F,2*N,iniPar,iniPar+2*N);
   FitUsingNewFitter<TMINUIT>(&tree,gausn);
   //}

}
