#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2DErrors.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"

#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
#include "HFitInterface.h"
#include "Fit/Fitter.h"

#include "Math/WrappedMultiTF1.h"
#include "Math/WrappedParamFunction.h"
#include "Math/WrappedTF1.h"
//#include "Math/Polynomial.h"
#include "RConfigure.h"

#include <string>
#include <iostream>
#include <cmath>

// print the data
void printData(const ROOT::Fit::BinData & data) {
   std::cout << "Bin data, data dimension is " << data.NDim() << " type is " << data.GetErrorType() << std::endl;
   for (unsigned int i = 0; i < data.Size(); ++i) {
      if (data.GetErrorType() == ROOT::Fit::BinData::kNoError)
         std::cout << data.Coords(i)[0] << "   " << data.Value(i) << std::endl;
      else if (data.GetErrorType() == ROOT::Fit::BinData::kValueError)
         std::cout << data.Coords(i)[0] << "   " << data.Value(i) << " +/-  " << data.Error(i) << std::endl;
      else if (data.GetErrorType() == ROOT::Fit::BinData::kCoordError) {
         double ey = 0;
         const double * ex = data.GetPointError(i,ey);
         std::cout << data.Coords(i)[0] <<  " +/-  " << ex[0] <<  "   " << data.Value(i) << " +/-  " << ey << std::endl;
      }
      else if (data.GetErrorType() == ROOT::Fit::BinData::kAsymError) {
         double eyl = 0; double eyh = 0;
         const double * ex = data.GetPointError(i,eyl,eyh);
         std::cout << data.Coords(i)[0] <<  " +/-  " << ex[0] <<  "   " << data.Value(i) << " -  " << eyl << " + " << eyh << std::endl;
      }

   }
   std::cout << "\ndata size is " << data.Size() << std::endl;
}
void printData(const ROOT::Fit::UnBinData & data) {
   for (unsigned int i = 0; i < data.Size(); ++i) {
      std::cout << data.Coords(i)[0] << "\t";
   }
   std::cout << "\ndata size is " << data.Size() << std::endl;
}

int compareResult(double v1, double v2, std::string s = "", double tol = 0.01) {
   // compare v1 with reference v2
   // give 1% tolerance
   if (std::abs(v1-v2) < tol * std::abs(v2) ) return 0;
   std::cerr << s << " Failed comparison of fit results \t chi2 = " << v1 << "   it should be = " << v2 << std::endl;
   return -1;
}

double chi2FromFit(const TF1 * func )  {
   // return last chi2 obtained from Fit method function
   R__ASSERT(TVirtualFitter::GetFitter() != 0 );
   return (TVirtualFitter::GetFitter()->Chisquare(func->GetNpar(), func->GetParameters() ) );
}

int testHisto1DFit() {


   std::string fname("gaus");
   TF1 * func = (TF1*)gROOT->GetFunction(fname.c_str());
   func->SetParameter(0,10);
   func->SetParameter(1,0);
   func->SetParameter(2,3.0);

   TRandom3 rndm;
   int iret = 0;
   double chi2ref = 0;

   // fill an histogram
   TH1D * h1 = new TH1D("h1","h1",30,-5.,5.);
//      h1->FillRandom(fname.c_str(),100);
   for (int i = 0; i <1000; ++i)
      h1->Fill( rndm.Gaus(0,1) );

   h1->Print();
   //h1->Draw();

//    gSystem->Load("libMinuit2");
//    gSystem->Load("libFit");


   // ROOT::Fit::DataVector<ROOT::Fit::BinPoint> dv;

   ROOT::Fit::BinData d;
   ROOT::Fit::FillData(d,h1,func);


   printData(d);

   // create the function
   ROOT::Math::WrappedMultiTF1 wf(*func);
   // need to do that to avoid gradient calculation
   ROOT::Math::IParamMultiFunction & f = wf;

   double p[3] = {100,0,3.};
   f.SetParameters(p);

   // create the fitter

   ROOT::Fit::Fitter fitter;

   bool ret = fitter.Fit(d, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Chi2 Fit Failed " << std::endl;
      return -1;
   }
   chi2ref = fitter.Result().Chi2();

   // compare with TH1::Fit
   TVirtualFitter::SetDefaultFitter("Minuit2");
   std::cout << "\n******************************\n\t TH1::Fit Result \n" << std::endl;
   func->SetParameters(p);
   auto res = h1->Fit(func,"S");

   iret |= compareResult( res->Chi2(), chi2ref,"1D histogram chi2 fit");


   // test using binned likelihood
   std::cout << "\n\nTest Binned Likelihood Fit" << std::endl;

   ROOT::Fit::DataOptions opt;
   opt.fUseEmpty = true;
   ROOT::Fit::BinData dl(opt);
   ROOT::Fit::FillData(dl,h1,func);

   ret = fitter.LikelihoodFit(dl, f, true);
   f.SetParameters(p);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Binned Likelihood Fit Failed " << std::endl;
      return -1;
   }
   iret |= compareResult(fitter.Result().Chi2(), chi2ref,"1D histogram likelihood fit",0.3);

   // compare with TH1::Fit
   std::cout << "\n******************************\n\t TH1::Fit Result \n" << std::endl;
   func->SetParameters(p);
   h1->Fit(func,"L");

   iret |= compareResult(func->GetChisquare(),fitter.Result().Chi2(),"TH1::Fit likelihood ",0.001);
   //std::cout << "Equivalent Chi2 from TF1::Fit " << func->GetChisquare() << std::endl;

   std::cout << "\n\nTest Chi2 Fit using integral option" << std::endl;

   // need to re-create data
   opt = ROOT::Fit::DataOptions();
   opt.fIntegral = true;
   ROOT::Fit::BinData d2(opt);
   ROOT::Fit::FillData(d2,h1,func);

   f.SetParameters(p);
   ret = fitter.Fit(d2, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Integral Chi2 Fit Failed " << std::endl;
      return -1;
   }
   iret |= compareResult(fitter.Result().Chi2(), chi2ref,"1D histogram integral chi2 fit",0.2);

   // compare with TH1::Fit
   std::cout << "\n******************************\n\t TH1::Fit Result \n" << std::endl;
   func->SetParameters(p);
   h1->Fit(func,"I ");
   iret |= compareResult(func->GetChisquare(),fitter.Result().Chi2(),"TH1::Fit integral ",0.001);

   // test integral likelihood
   std::cout << "\n\nTest Likelihood Fit using integral option" << std::endl;
   opt = ROOT::Fit::DataOptions();
   opt.fIntegral = true;
   opt.fUseEmpty = true;
   ROOT::EExecutionPolicy execPolicy = ROOT::EExecutionPolicy::kSequential;

   // if (ROOT::IsImplicitMTEnabled()) {
   //    execPolicy = ROOT::EExecutionPolicy::kMultiThread;
   // }
   ROOT::Fit::BinData dl2(opt);
   ROOT::Fit::FillData(dl2,h1,func);
   f.SetParameters(p);
   fitter.SetFunction(f);
   ret = fitter.LikelihoodFit(dl2, true, execPolicy);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      fitter.Result().Print(std::cout);
      std::cout << "ERROR: Integral Likelihood Fit Failed " << std::endl;
      //return -1;
   }
   iret |= compareResult(fitter.Result().Chi2(), chi2ref,"1D histogram integral likelihood fit",0.3);

   // compare with TH1::Fit
   std::cout << "\n******************************\n\t TH1::Fit Result \n" << std::endl;
   func->SetParameters(p);
   h1->Fit(func,"IL");
   //std::cout << "Equivalent Chi2 from TF1::Fit " << func->GetChisquare() << std::endl;
   iret |= compareResult(func->GetChisquare(),fitter.Result().Chi2(),"TH1::Fit likelihood integral ",0.001);



   // redo chi2fit
   std::cout << "\n\nRedo Chi2 Hist Fit" << std::endl;
   f.SetParameters(p);
   ret = fitter.Fit(d, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Chi2 Fit Failed " << std::endl;
      return -1;
   }
   iret |= compareResult(fitter.Result().Chi2(), chi2ref,"1D histogram chi2 fit (2)",0.001);



   // test grapherrors fit
   std::cout << "\n\nTest Same Fit from a TGraphErrors - no coord errors" << std::endl;
   TGraphErrors gr(h1);
   ROOT::Fit::BinData dg;
   dg.Opt().fCoordErrors = false;  // do not use coordinate errors (default is using )
   ROOT::Fit::FillData(dg,&gr);

   f.SetParameters(p);
   ret = fitter.Fit(dg, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Chi2 Graph Errors Fit Failed " << std::endl;
      return -1;
   }
   iret |= compareResult(fitter.Result().Chi2(), chi2ref,"TGraphErrors chi2 fit",0.001);


   // fit using error on X
   std::cout << "\n\nTest Same Fit from a TGraphErrors - use coord errors" << std::endl;
   ROOT::Fit::BinData dger;
   // not needed since they are used by default
   //dger.Opt().fCoordErrors = true;  // use coordinate errors
   dger.Opt().fUseEmpty = true;  // this will set error 1 for the empty bins
   ROOT::Fit::FillData(dger,&gr);

   f.SetParameters(p);
   ret = fitter.Fit(dger, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Chi2 Graph Errors Fit Failed " << std::endl;
      return -1;
   }
   iret |= compareResult(fitter.Result().Chi2(), chi2ref,"TGraphErrors effective chi2 fit ",0.7);

   // compare with TGraphErrors::Fit
   std::cout << "\n******************************\n\t TGraphErrors::Fit Result \n" << std::endl;
   func->SetParameters(p);
   // set error = 1 for empty bins
   for (int ip = 0; ip < gr.GetN(); ++ip)
      if (gr.GetErrorY(ip) <= 0) gr.SetPointError(ip, gr.GetErrorX(ip), 1.);

   gr.Fit(func);
   std::cout << "Ndf of TGraphErrors::Fit  = " << func->GetNDF() << std::endl;
   iret |= compareResult(func->GetChisquare(),fitter.Result().Chi2(),"TGraphErrors::Fit ",0.001);


   // test graph fit (errors are 1) do a re-normalization
   std::cout << "\n\nTest Same Fit from a TGraph" << std::endl;
   fitter.Config().SetNormErrors(true);
   TGraph gr2(h1);
   ROOT::Fit::BinData dg2;
   ROOT::Fit::FillData(dg2,&gr2);

   f.SetParameters(p);
   ret = fitter.Fit(dg2, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Chi2 Graph Fit Failed " << std::endl;
      return -1;
   }
   //iret |= compareResult(fitter.Result().Chi2(), chi2ref,"TGraph fit (no errors) ",0.3);


   // compare with TGraph::Fit (no errors)
   std::cout << "\n******************************\n\t TGraph::Fit Result \n" << std::endl;
   func->SetParameters(p);
   gr2.Fit(func);
   std::cout << "Ndf of TGraph::Fit = " << func->GetNDF() << std::endl;

   iret |= compareResult(func->GetChisquare(),fitter.Result().Chi2(),"TGraph::Fit ",0.001);


   // reddo chi2fit using Fumili2
   std::cout << "\n\nRedo Chi2 Hist Fit using FUMILI2" << std::endl;
   f.SetParameters(p);
   fitter.Config().SetMinimizer("Minuit2","Fumili");
   ret = fitter.Fit(d, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Chi2 Fit Failed " << std::endl;
      return -1;
   }
   iret |= compareResult(fitter.Result().Chi2(), chi2ref,"1D Histo Fumili2 fit ");

   // reddo chi2fit using old Fumili
   std::cout << "\n\nRedo Chi2 Hist Fit using FUMILI" << std::endl;
   f.SetParameters(p);
   fitter.Config().SetMinimizer("Fumili");
   ret = fitter.Fit(d, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Chi2 Fit Failed " << std::endl;
      return -1;
   }
   iret |= compareResult(fitter.Result().Chi2(), chi2ref,"1D Histo Fumili fit ");

   // test using GSL multi fit (L.M. method)
   std::cout << "\n\nRedo Chi2 Hist Fit using GSLMultiFit" << std::endl;
   f.SetParameters(p);
   fitter.Config().SetMinimizer("GSLMultiFit");
   ret = fitter.Fit(d, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Chi2 Fit Failed " << std::endl;
      return -1;
   }
   iret |= compareResult(fitter.Result().Chi2(), chi2ref,"1D Histo GSL NLS fit ");

   // test using GSL multi min method
   std::cout << "\n\nRedo Chi2 Hist Fit using GSLMultiMin" << std::endl;
   f.SetParameters(p);
   fitter.Config().SetMinimizer("GSLMultiMin","BFGS2");
   ret = fitter.Fit(d, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Chi2 Fit Failed " << std::endl;
      return -1;
   }
   iret |= compareResult(fitter.Result().Chi2(), chi2ref,"1D Histo GSL Minimizer fit ");

   return iret;
}


class Func1D : public ROOT::Math::IParamFunction {
public:
   void SetParameters(const double *p) { std::copy(p,p+NPar(),fp);}
   const double * Parameters() const { return fp; }
   ROOT::Math::IGenFunction * Clone() const {
      Func1D * f =  new Func1D();
      f->SetParameters(fp);
      return f;
   };
   unsigned int NPar() const { return 3; }
private:
   double DoEvalPar( double x, const double *p) const {
      return p[0]*x*x + p[1]*x + p[2];
   }
   double fp[3];

};

// gradient 2D function
class GradFunc2D : public ROOT::Math::IParamMultiGradFunction {
public:
   void SetParameters(const double *p) { std::copy(p,p+NPar(),fp);}
   const double * Parameters() const { return fp; }
   ROOT::Math::IMultiGenFunction * Clone() const {
      GradFunc2D * f =  new GradFunc2D();
      f->SetParameters(fp);
      return f;
   };
   unsigned int NDim() const { return 2; }
   unsigned int NPar() const { return 5; }

   void ParameterGradient( const double * x, const double * , double * grad) const {
      grad[0] = x[0]*x[0];
      grad[1] = x[0];
      grad[2] = x[1]*x[1];
      grad[3] = x[1];
      grad[4] = 1;
   }

private:

   double DoEvalPar( const double *x, const double * p) const {
      return p[0]*x[0]*x[0] + p[1]*x[0] + p[2]*x[1]*x[1] + p[3]*x[1] + p[4];
   }
//    double DoDerivative(const double *x,  unsigned int icoord = 0) const {
//       assert(icoord <= 1);
//       if (icoord == 0)
//          return 2. * fp[0] * x[0] + fp[1];
//       else
//          return 2. * fp[2] * x[1] + fp[3];
//    }

   double DoParameterDerivative(const double * x, const double * p, unsigned int ipar) const {
      std::vector<double> grad(NPar());
      ParameterGradient(x, p, &grad[0] );
      return grad[ipar];
   }

   double fp[5];

};

int testHisto1DPolFit() {



   std::string fname("pol2");
   TF1 * func = (TF1*)gROOT->GetFunction(fname.c_str());
   func->SetParameter(0,1.);
   func->SetParameter(1,2.);
   func->SetParameter(2,3.0);

   TRandom3 rndm;
   int iret = 0;
   double chi2ref = 0;

   // fill an histogram
   TH1D * h2 = new TH1D("h2","h2",30,-5.,5.);
//      h1->FillRandom(fname.c_str(),100);
   for (int i = 0; i <1000; ++i)
      h2->Fill( func->GetRandom() );

   // fill fit data
   ROOT::Fit::BinData d;
   ROOT::Fit::FillData(d,h2,func);


   printData(d);

   // create the function
   Func1D f;

   double p[3] = {100,0,3.};
   f.SetParameters(p);


   // create the fitter
   //std::cout << "Fit parameter 2  " << f.Parameters()[2] << std::endl;
   std::cout << "\n\nTest histo polynomial fit using Fitter" << std::endl;

   ROOT::Fit::Fitter fitter;
   bool ret = fitter.Fit(d, f);
   if (ret)
      fitter.Result().Print(std::cout,true);
   else {
      std::cout << " Fit Failed " << std::endl;
      return -1;
   }
   chi2ref = fitter.Result().Chi2();

   // compare with TH1::Fit
   std::cout << "\n******************************\n\t TH1::Fit(pol2) Result   \n" << std::endl;
   func->SetParameters(p);
   h2->Fit(func,"F");
   iret |= compareResult(func->GetChisquare(),chi2ref,"TH1::Fit ",0.001);


   std::cout << "\n\nTest histo polynomial linear fit " << std::endl;

   ROOT::Math::WrappedTF1 pf(*func);
   //ROOT::Math::Polynomial pf(2);
   pf.SetParameters(p);

   fitter.Config().SetMinimizer("Linear");
   ret = fitter.Fit(d, pf);
   if (ret)  {
      fitter.Result().Print(std::cout);
      fitter.Result().PrintCovMatrix(std::cout);
   }
   else {
      std::cout << " Fit Failed " << std::endl;
      return -1;
   }
   iret |= compareResult(fitter.Result().Chi2(),chi2ref,"1D histo linear Fit ");

   // compare with TH1::Fit
   std::cout << "\n******************************\n\t TH1::Fit(pol2) Result with TLinearFitter \n" << std::endl;
   func->SetParameters(p);
   h2->Fit(func);
   iret |= compareResult(func->GetChisquare(),fitter.Result().Chi2(),"TH1::Fit linear",0.001);



   return iret;
}

int testHisto2DFit() {

   // fit using a 2d parabola (test also gradient)


   std::string fname("pol2");
   TF2 * func = new TF2("f2d",ROOT::Math::ParamFunctor(GradFunc2D() ), 0.,10.,0,10,5);
   double p0[5] = { 1.,2.,0.5,1.,3. };
   func->SetParameters(p0);
   assert(func->GetNpar() == 5);

   TRandom3 rndm;
   double chi2ref = 0;
   int iret = 0;

   // fill an histogram
   TH2D * h2 = new TH2D("h2d","h2d",30,0,10.,30,0.,10.);
//      h1->FillRandom(fname.c_str(),100);
   for (int i = 0; i <10000; ++i) {
      double x,y = 0;
      func->GetRandom2(x,y);
      h2->Fill(x,y);
   }
   // fill fit data
   ROOT::Fit::BinData d;
   ROOT::Fit::FillData(d,h2,func);


   //printData(d);

   // create the function
   GradFunc2D f;

   double p[5] = { 2.,1.,1,2.,100. };
   f.SetParameters(p);


   // create the fitter
   ROOT::Fit::Fitter fitter;
   //fitter.Config().MinimizerOptions().SetPrintLevel(3);
   fitter.Config().SetMinimizer("Minuit2");

   std::cout <<"\ntest 2D histo fit using gradient" << std::endl;
   bool ret = fitter.Fit(d, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Gradient Fit Failed " << std::endl;
      return -1;
   }
   chi2ref = fitter.Result().Chi2();

   // test without gradient
   std::cout <<"\ntest result without using gradient" << std::endl;
   ROOT::Math::WrappedParamFunction<GradFunc2D *> f2(&f,2,5,p);
   ret = fitter.Fit(d, f2);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << " Chi2 Fit Failed " << std::endl;
      return -1;
   }
   iret |= compareResult(fitter.Result().Chi2(), chi2ref,"2D histogram chi2 fit");

   // test Poisson bin likelihood fit (no gradient)
   std::cout <<"\ntest result without gradient and binned likelihood" << std::endl;
   f.SetParameters(p);
   fitter.SetFunction(static_cast<const ROOT::Math::IParamMultiFunction &>(f) );
   fitter.Config().ParSettings(0).SetLimits(0,100);
   fitter.Config().ParSettings(1).SetLimits(0,100);
   fitter.Config().ParSettings(2).SetLimits(0,100);
   fitter.Config().ParSettings(3).SetLimits(0,100);
   fitter.Config().ParSettings(4).SetLowerLimit(0);
   //fitter.Config().MinimizerOptions().SetPrintLevel(3);
   ret = fitter.LikelihoodFit(d);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Poisson 2D Bin Likelihood  Fit Failed " << std::endl;
      return -1;
   }

   // test binned likelihood gradient
   std::cout <<"\ntest result using gradient and binned likelihood" << std::endl;
   f.SetParameters(p);
   fitter.SetFunction(f);
   //fitter.Config().MinimizerOptions().SetPrintLevel(3);
   fitter.Config().ParSettings(0).SetLimits(0,100);
   fitter.Config().ParSettings(1).SetLimits(0,100);
   fitter.Config().ParSettings(2).SetLimits(0,100);
   fitter.Config().ParSettings(3).SetLimits(0,100);
   fitter.Config().ParSettings(4).SetLowerLimit(0);
   ret = fitter.LikelihoodFit(d);
   if (ret)  {
      // redo fit releasing the parameters
      f.SetParameters(&(fitter.Result().Parameters().front()) );
      ret = fitter.LikelihoodFit(d,f, true);
   }
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Gradient Bin Likelihood  Fit Failed " << std::endl;
      return -1;
   }

   // test with linear fitter
   std::cout <<"\ntest result using linear fitter" << std::endl;
   fitter.Config().SetMinimizer("Linear");
   f.SetParameters(p);
   ret = fitter.Fit(d, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Linear 2D Fit Failed " << std::endl;
      return -1;
   }
   iret |= compareResult(fitter.Result().Chi2(), chi2ref,"2D histogram linear fit");

   // test fitting using TGraph2D ( chi2 will be larger since errors are 1)
   // should test with a TGraph2DErrors
   TGraph2D g2(h2);

   std::cout <<"\ntest using TGraph2D" << std::endl;
   ROOT::Fit::BinData d2;
   ROOT::Fit::FillData(d2,&g2,func);
   //g2.Dump();
   std::cout << "data size from graph " << d2.Size() <<  std::endl;

   f2.SetParameters(p);
   fitter.Config().SetMinimizer("Minuit2");
   ret = fitter.Fit(d2, f2);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << " TGraph2D Fit Failed " << std::endl;
      return -1;
   }
   double chi2ref2 = fitter.Result().Chi2();

   // compare with TGraph2D::Fit
   std::cout << "\n******************************\n\t TGraph::Fit Result \n" << std::endl;
   func->SetParameters(p);
   g2.Fit(func);

   std::cout <<"\ntest using TGraph2D and gradient function" << std::endl;
   f.SetParameters(p);
   ret = fitter.Fit(d2, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << " TGraph2D Grad Fit Failed " << std::endl;
      return -1;
   }
   iret |= compareResult(fitter.Result().Chi2(), chi2ref2,"TGraph2D chi2 fit");

   std::cout <<"\ntest using TGraph2DErrors -  error only in Z" << std::endl;
   TGraph2DErrors g2err(g2.GetN() );
   // need to set error by hand since constructor from TH2 does not exist
   for (int i = 0; i < g2.GetN(); ++i) {
      double x = g2.GetX()[i];
      double y = g2.GetY()[i];
      g2err.SetPoint(i,x,y,g2.GetZ()[i]);
      g2err.SetPointError(i,0,0,h2->GetBinError(h2->FindBin(x,y) ) );
   }
   func->SetParameters(p);
   // g2err.Fit(func);
   f.SetParameters(p);
   ROOT::Fit::BinData d3;
   ROOT::Fit::FillData(d3,&g2err,func);
   ret = fitter.Fit(d3, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << " TGraph2DErrors Fit Failed " << std::endl;
      return -1;
   }

   iret |= compareResult(fitter.Result().Chi2(), chi2ref,"TGraph2DErrors chi2 fit");



   std::cout <<"\ntest using TGraph2DErrors -  with error  in X,Y,Z" << std::endl;
   for (int i = 0; i < g2err.GetN(); ++i) {
      double x = g2.GetX()[i];
      double y = g2.GetY()[i];
      g2err.SetPointError(i,0.5* h2->GetXaxis()->GetBinWidth(1),0.5*h2->GetXaxis()->GetBinWidth(1),h2->GetBinError(h2->FindBin(x,y) ) );
   }
   std::cout << "\n******************************\n\t TGraph2DErrors::Fit Result \n" << std::endl;
   func->SetParameters(p);
   g2err.Fit(func);


   return iret;
}


int testUnBin1DFit() {

   int iret = 0;

   std::string fname("gausn");
   TF1 * func = (TF1*)gROOT->GetFunction(fname.c_str());

   TRandom3 rndm;

   int n = 100;
   ROOT::Fit::UnBinData d(n);

   for (int i = 0; i <n; ++i)
      d.Add( rndm.Gaus(0,1) );


   // printData(d);

   // create the function
   ROOT::Math::WrappedMultiTF1 wf(*func);
   // need to do that to avoid gradient calculation
   ROOT::Math::IParamMultiFunction & f = wf;
   double p[3] = {1,2,10.};
   f.SetParameters(p);

   // create the fitter
   //std::cout << "Fit parameters  " << f.Parameters()[2] << std::endl;

   ROOT::Fit::Fitter fitter;
   fitter.SetFunction(f);
   std::cout << "fix parameter 0 " << " to value " << f.Parameters()[0] << std::endl;
   fitter.Config().ParSettings(0).Fix();
   // set range in sigma sigma > 0
   std::cout << "set lower range to 0 for sigma " << std::endl;
   fitter.Config().ParSettings(2).SetLowerLimit(0);

#ifdef DEBUG
   fitter.Config().MinimizerOptions().SetPrintLevel(3);
#endif

//    double x[1]; x[0] = 0.;
//    std::cout << "fval " << f(x) << std::endl;
//    x[0] = 1.; std::cout << "fval " << f(x) << std::endl;

   //fitter.Config().SetMinimizer("Minuit2");
   // fails with minuit (t.b. investigate)
   fitter.Config().SetMinimizer("Minuit2");


   bool ret = fitter.Fit(d);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Unbinned Likelihood Fit Failed " << std::endl;
      iret |= 1;
   }
   double lref = fitter.Result().MinFcnValue();

   std::cout << "\n\nRedo Fit using FUMILI2" << std::endl;
   f.SetParameters(p);
   fitter.Config().SetMinimizer("Fumili2");
   // need to set function first (need to change this)
   fitter.SetFunction(f);
   fitter.Config().ParSettings(0).Fix(); //need to re-do it
   // set range in sigma sigma > 0
   fitter.Config().ParSettings(2).SetLowerLimit(0);

   ret = fitter.Fit(d);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Unbinned Likelihood Fit using FUMILI2 Failed " << std::endl;
      iret |= 1;
   }

   iret |= compareResult(fitter.Result().MinFcnValue(), lref,"1D unbin FUMILI2 fit");

   std::cout << "\n\nRedo Fit using FUMILI" << std::endl;
   f.SetParameters(p);
   fitter.Config().SetMinimizer("Fumili");
   // fitter.Config().MinimizerOptions().SetPrintLevel(3);
   // need to set function first (need to change this)
   fitter.SetFunction(f);
   fitter.Config().ParSettings(0).Fix(); //need to re-do it
   // set range in sigma sigma > 0
   fitter.Config().ParSettings(2).SetLowerLimit(0);

   ret = fitter.Fit(d);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Unbinned Likelihood Fit using FUMILI Failed " << std::endl;
      iret |= 1;
   }

   iret |= compareResult(fitter.Result().MinFcnValue(), lref,"1D unbin FUMILI fit");


   return iret;
}

int testGraphFit() {

   int iret = 0;

   // simple test of fitting a Tgraph

   double x[5] = {1,2,3,4,5};
   double y[5] = {2.1, 3.5, 6.5, 8.8, 9.5};
   double ex[5] = {.3,.3,.3,.3,.3};
   double ey[5] = {.5,.5,.5,.5,.5};
   double eyl[5] = {.2,.2,.2,.2,.2};
   double eyh[5] = {.8,.8,.8,.8,.8};

   std::cout << "\n********************************************************\n";
   std::cout << "Test simple fit of Tgraph of 5 points" << std::endl;
   std::cout << "\n********************************************************\n";


   double p[2] = {1,1};
   TF1 * func = new TF1("f","pol1",0,10);
   func->SetParameters(p);

   ROOT::Math::WrappedMultiTF1 wf(*func);
   // need to do that to avoid gradient calculation
   ROOT::Math::IParamMultiFunction & f = wf;
   f.SetParameters(p);

   ROOT::Fit::Fitter fitter;
   fitter.SetFunction(f);


   std::cout <<"\ntest TGraph (no errors) " << std::endl;
   TGraph gr(5, x,y);

   ROOT::Fit::BinData dgr;
   ROOT::Fit::FillData(dgr,&gr);

   //printData(dgr);

   f.SetParameters(p);
   bool ret = fitter.Fit(dgr, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Chi2 Graph Fit Failed " << std::endl;
      return -1;
   }
   double chi2ref = fitter.Result().Chi2();


   // compare with TGraph::Fit
   std::cout << "\n******************************\n\t TGraph::Fit Result \n" << std::endl;
   func->SetParameters(p);
   gr.Fit(func,"F"); // use Minuit

   iret |= compareResult(func->GetChisquare(),fitter.Result().Chi2(),"TGraph::Fit ",0.001);


   std::cout <<"\ntest TGraphErrors  " << std::endl;
   TGraphErrors grer(5, x,y,ex,ey);

   ROOT::Fit::BinData dgrer;
   dgrer.Opt().fCoordErrors = true;
   ROOT::Fit::FillData(dgrer,&grer);

   //printData(dgrer);

   f.SetParameters(p);
   ret = fitter.Fit(dgrer, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Chi2 Graph Fit Failed " << std::endl;
      return -1;
   }

   iret |= compareResult(fitter.Result().Chi2(),chi2ref,"TGraphErrors fit with coord errors",0.8);


   // compare with TGraph::Fit
   std::cout << "\n******************************\n\t TGraphErrors::Fit Result \n" << std::endl;
   func->SetParameters(p);
   grer.Fit(func,"F"); // use Minuit
   iret |= compareResult(func->GetChisquare(),fitter.Result().Chi2(),"TGraphErrors::Fit ",0.001);

   chi2ref = fitter.Result().Chi2();

   std::cout <<"\ntest TGraphAsymmErrors  " << std::endl;
   TGraphAsymmErrors graer(5, x,y,ex,ex,eyl, eyh);

   ROOT::Fit::BinData dgraer;
   // option error on coordinate and asymmetric on values
   dgraer.Opt().fCoordErrors = true;
   dgraer.Opt().fAsymErrors = true;
   ROOT::Fit::FillData(dgraer,&graer);
   //printData(dgraer);

   f.SetParameters(p);
   ret = fitter.Fit(dgraer, f);
   if (ret)
      fitter.Result().Print(std::cout);
   else {
      std::cout << "Chi2 Graph Fit Failed " << std::endl;
      return -1;
   }
   //iret |= compareResult(fitter.Result().Chi2(),chi2ref,"TGraphAsymmErrors fit ",0.5);

   // compare with TGraph::Fit
   std::cout << "\n******************************\n\t TGraphAsymmErrors::Fit Result \n" << std::endl;
   func->SetParameters(p);
   graer.Fit(func,"F"); // use Minuit
   iret |= compareResult(func->GetChisquare(),fitter.Result().Chi2(),"TGraphAsymmErrors::Fit ",0.001);



   return iret;
}

int testFitResultScan() {
   TH1D *h1 = new TH1D("h1", "h1", 30, -3., 3.);
   h1->FillRandom("gaus",500);

   auto r = h1->Fit("gaus", "S");
   if (r->Status() != 0)
      return r->Status();
   unsigned int np = 20;
   std::vector<double> x(np);
   std::vector<double> y(np);
   bool ok = r->Scan(1, np, x.data(), y.data(), -1., 1.);
   if (!ok) return -1;
   std::cout << "scan points:  " << np << std::endl;
   for (unsigned int i = 0; i < np; ++i)
      std::cout << "( " << x[i] << " , " << y[i] << " )  ";
   std::cout << std::endl;
   // check that is a parabola
   unsigned int imin = TMath::LocMin(np, y.data());
   double xmin = x[imin];
   std::cout << "Minimum of scan for parameter 1 (Gaussian mean) is " << xmin << " index " << imin << std::endl;
   return !( (std::abs(xmin) < 0.3) && std::abs(double(imin)-np/2.) < 2. );
}

template<typename Test>
int testFit(Test t, std::string name) {
   std::cout << "******************************\n";
   std::cout << name << "\n\n\t\t";
   int iret = t();
   std::cout << "\n" << name << ":\t\t";
   if (iret == 0)
      std::cout << "OK" << std::endl;
   else
      std::cout << "Failed" << std::endl;
   return iret;
}

int runAllTests() {

   int iret = 0;
   iret |= testFit( testHisto1DFit, "Histogram1D Fit");
   iret |= testFit( testHisto1DPolFit, "Histogram1D Polynomial Fit");
   iret |= testFit( testHisto2DFit, "Histogram2D Gradient Fit");
   iret |= testFit( testUnBin1DFit, "Unbin 1D Fit");
   iret |= testFit( testGraphFit, "Graph 1D Fit");
   iret |= testFit( testFitResultScan, "FitResult Scan");

   std::cout << "\n******************************\n";
   if (iret) std::cerr << "\n\t testFit FAILED !!!!!!!!!!!!!!!! \n";
   else std::cout << "\n\t testFit all OK  !\n";
   return iret;
}
int main() {
   runAllTests();

   std::cout << "\n******************************\n";
   std::cout << "Test NOW in Multithreading mode " << std::endl;
   std::cout << "\n******************************\n";
   ROOT::EnableImplicitMT(1);
   runAllTests();
}
