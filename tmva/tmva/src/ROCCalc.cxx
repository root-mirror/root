/*! \class TMVA::ROCCalc
\ingroup TMVA

*/

#include <algorithm>
#include <cstdlib>
#include <cerrno>

#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TList.h"
#include "TSpline.h"
#include "TVector.h"
#include "TMatrixD.h"
#include "TMatrixDSymEigen.h"
#include "TTreeFormula.h"
#include "TXMLEngine.h"
#include "TROOT.h"
#include "TColor.h"
#include "TGraph.h"

#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/ROCCalc.h"
#include "TMVA/Event.h"
#include "TMVA/Version.h"
#include "TMVA/PDF.h"
#include "TMVA/MsgLogger.h"

#include "TMVA/TSpline1.h"
#include "TMVA/TSpline2.h"
#include "TMVA/Types.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////

TMVA::ROCCalc::ROCCalc(TH1* mvaS, TH1* mvaB) :
   fMaxIter(100),
   fAbsTol(0.0),
   fStatus(kTRUE),
   fmvaS(0),
   fmvaB(0),
   fmvaSpdf(0),
   fmvaBpdf(0),
   fSplS(0),
   fSplB(0),
   fSplmvaCumS(0),
   fSplmvaCumB(0),
   fSpleffBvsS(0),
   fnStot(0),
   fnBtot(0),
   fSignificance(0),
   fPurity(0),
   effBvsS(0),
   rejBvsS(0),
   inveffBvsS(0),
   fLogger ( new TMVA::MsgLogger("ROCCalc") )
{
   fUseSplines = kTRUE;
   fNbins      = 100;
   // fmvaS = (TH1*) mvaS->Clone("MVA Signal"); fmvaS->SetTitle("MVA Signal");
   // fmvaB = (TH1*) mvaB->Clone("MVA Backgr"); fmvaB->SetTitle("MVA Backgr");
   fmvaS =  mvaS; fmvaS->SetTitle("MVA Signal");
   fmvaB =  mvaB; fmvaB->SetTitle("MVA Backgr");
   fXmax = fmvaS->GetXaxis()->GetXmax();
   fXmin = fmvaS->GetXaxis()->GetXmin();

   if (TMath::Abs(fXmax-fmvaB->GetXaxis()->GetXmax()) > 0.000001 ||
       TMath::Abs(fXmin-fmvaB->GetXaxis()->GetXmin()) > 0.000001 ||
       fmvaB->GetNbinsX() != fmvaS->GetNbinsX()) {
      Log() << kERROR << "Cannot cal ROC curve etc, as in put mvaS and mvaB have differen #nbins or range "<<Endl;
      fStatus=kFALSE;
   }
   if (!strcmp(fmvaS->GetXaxis()->GetTitle(),"")) fmvaS->SetXTitle("MVA-value");
   if (!strcmp(fmvaB->GetXaxis()->GetTitle(),"")) fmvaB->SetXTitle("MVA-value");
   if (!strcmp(fmvaS->GetYaxis()->GetTitle(),"")) fmvaS->SetYTitle("#entries");
   if (!strcmp(fmvaB->GetYaxis()->GetTitle(),"")) fmvaB->SetYTitle("#entries");
   ApplySignalAndBackgroundStyle(fmvaS, fmvaB);
   //    std::cout<<"mvaS->GetNbinsX()"<<mvaS->GetNbinsX()<<std::endl;
   //    std::cout<<"mvaB->GetNbinsX()"<<mvaB->GetNbinsX()<<std::endl;
   //the output of mvaS->GetNbinsX() is about 40 and if we divide it by 100 the results is 0
   //the I will divide it by 10 anyway doing some tests ROC integral is the same
   fmvaSpdf = mvaS->RebinX(mvaS->GetNbinsX()/10,"MVA Signal PDF");
   fmvaBpdf = mvaB->RebinX(mvaB->GetNbinsX()/10,"MVA Backgr PDF");
   if(fmvaSpdf==0||fmvaBpdf==0)
      {
         Log() << kERROR << "Cannot Rebin Histograms mvaS and mvaB, ROC values will be calculated without Rebin histograms."<<Endl;
         fStatus=kFALSE;
         fmvaSpdf = (TH1*)mvaS->Clone("MVA Signal PDF");
         fmvaBpdf = (TH1*)mvaB->Clone("MVA Backgr PDF");
      }
   fmvaSpdf->SetTitle("MVA Signal PDF");
   fmvaBpdf->SetTitle("MVA Backgr PDF");
   fmvaSpdf->Scale(1./fmvaSpdf->GetSumOfWeights());
   fmvaBpdf->Scale(1./fmvaBpdf->GetSumOfWeights());
   fmvaSpdf->SetMaximum(TMath::Max(fmvaSpdf->GetMaximum(), fmvaBpdf->GetMaximum()));
   fmvaBpdf->SetMaximum(TMath::Max(fmvaSpdf->GetMaximum(), fmvaBpdf->GetMaximum()));
   ApplySignalAndBackgroundStyle(fmvaSpdf, fmvaBpdf);

   fCutOrientation = (fmvaS->GetMean() > fmvaB->GetMean()) ? +1 : -1;

   fNevtS = 0;

}

////////////////////////////////////////////////////////////////////////////////

void TMVA::ROCCalc::ApplySignalAndBackgroundStyle( TH1* sig, TH1* bkg, TH1* any ) {
   Int_t c_SignalLine     = TColor::GetColor( "#0000ee" );
   Int_t c_SignalFill     = TColor::GetColor( "#7d99d1" );
   Int_t c_BackgroundLine = TColor::GetColor( "#ff0000" );
   Int_t c_BackgroundFill = TColor::GetColor( "#ff0000" );
   //  Int_t c_NovelBlue      = TColor::GetColor( "#2244a5" );

   //signal
   // const Int_t FillColor__S = 38 + 150; // change of Color Scheme in ROOT-5.16.
   // convince yourself with gROOT->GetListOfColors()->Print()
   Int_t FillColor__S = c_SignalFill;
   Int_t FillStyle__S = 1001;
   Int_t LineColor__S = c_SignalLine;
   Int_t LineWidth__S = 2;

   // background
   //Int_t icolor = gConfig().fVariablePlotting.fUsePaperStyle ? 2 + 100 : 2;
   Int_t FillColor__B = c_BackgroundFill;
   Int_t FillStyle__B = 3554;
   Int_t LineColor__B = c_BackgroundLine;
   Int_t LineWidth__B = 2;

   if (sig != NULL) {
      sig->SetLineColor( LineColor__S );
      sig->SetLineWidth( LineWidth__S );
      sig->SetFillStyle( FillStyle__S );
      sig->SetFillColor( FillColor__S );
   }

   if (bkg != NULL) {
      bkg->SetLineColor( LineColor__B );
      bkg->SetLineWidth( LineWidth__B );
      bkg->SetFillStyle( FillStyle__B );
      bkg->SetFillColor( FillColor__B );
   }

   if (any != NULL) {
      any->SetLineColor( LineColor__S );
      any->SetLineWidth( LineWidth__S );
      any->SetFillStyle( FillStyle__S );
      any->SetFillColor( FillColor__S );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::ROCCalc::~ROCCalc() {
   // delete Splines and all histograms that were created only for internal use
   if (fSplS)            { delete fSplS; fSplS = 0; }
   if (fSplB)            { delete fSplB; fSplB = 0; }
   if (fSpleffBvsS)      { delete fSpleffBvsS; fSpleffBvsS = 0; }
   if (fSplmvaCumS)      { delete fSplmvaCumS; fSplmvaCumS = 0; }
   if (fSplmvaCumB)      { delete fSplmvaCumB; fSplmvaCumB = 0; }
   if (fmvaScumul)       { delete fmvaScumul; }
   if (fmvaBcumul)       { delete fmvaBcumul; }
   if (effBvsS)          { delete effBvsS; }
   if (rejBvsS)          { delete rejBvsS; }
   if (inveffBvsS)       { delete inveffBvsS; }
   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////
/// get the ROC curve

TH1D* TMVA::ROCCalc::GetROC(){
   // first get the cumulative distributions of the mva distribution
   // --> efficiencies vs cut value
   fNevtS = fmvaS->GetSumOfWeights(); // needed to get the error on the eff.. will only be correct if the histogram is not scaled to "integral == 1" Yet;
   if (fNevtS < 2) {
      Log() << kERROR << "I guess the mva distributions fed into ROCCalc were already normalized, therefore the calculated error on the efficiency will be incorrect !! " << Endl;
      fNevtS = 0;  // reset to zero --> no error will be calculated on the efficiencies
      fStatus=kFALSE;
   }
   fmvaScumul = gTools().GetCumulativeDist(fmvaS);
   fmvaBcumul = gTools().GetCumulativeDist(fmvaB);
   fmvaScumul->Scale( 1.0/TMath::Max(std::numeric_limits<double>::epsilon(),fmvaScumul->GetMaximum()) );
   fmvaBcumul->Scale( 1.0/TMath::Max(std::numeric_limits<double>::epsilon(),fmvaBcumul->GetMaximum()) );
   fmvaScumul->SetMinimum(0);
   fmvaBcumul->SetMinimum(0);
   //   fmvaScumul->Draw("hist");
   //   fmvaBcumul->Draw("histsame");

   // background efficiency versus signal efficiency
   if(effBvsS==0) effBvsS = new TH1D("effBvsS", "ROC-Curve", fNbins, 0, 1 );
   effBvsS->SetXTitle( "Signal eff" );
   effBvsS->SetYTitle( "Backgr eff" );

   // background rejection (=1-eff.) versus signal efficiency
   if(rejBvsS==0) rejBvsS = new TH1D( "rejBvsS", "ROC-Curve", fNbins, 0, 1 );
   rejBvsS->SetXTitle( "Signal eff" );
   rejBvsS->SetYTitle( "Backgr rejection (1-eff)" );

   // inverse background eff (1/eff.) versus signal efficiency
   if(inveffBvsS ==0) inveffBvsS = new TH1D("invBeffvsSeff", "ROC-Curve" , fNbins, 0, 1 );
   inveffBvsS->SetXTitle( "Signal eff" );
   inveffBvsS->SetYTitle( "Inverse backgr. eff (1/eff)" );

   // use root finder
   // spline background efficiency plot
   // note that there is a bin shift when going from a TH1D object to a TGraph :-(
   if (fUseSplines) {
      fSplmvaCumS  = new TSpline1( "spline2_signal",     new TGraph( fmvaScumul ) );
      fSplmvaCumB  = new TSpline1( "spline2_background", new TGraph( fmvaBcumul ) );
      // verify spline sanity
      gTools().CheckSplines( fmvaScumul, fSplmvaCumS );
      gTools().CheckSplines( fmvaBcumul, fSplmvaCumB );
   }

   Double_t effB = 0;
   for (UInt_t bini=1; bini<=fNbins; bini++) {

      // find cut value corresponding to a given signal efficiency
      Double_t effS = effBvsS->GetBinCenter( bini );
      Double_t cut  = Root( effS );

      // retrieve background efficiency for given cut
      if (fUseSplines) effB = fSplmvaCumB->Eval( cut );
      else             effB = fmvaBcumul->GetBinContent( fmvaBcumul->FindBin( cut ) );

      // and fill histograms
      effBvsS->SetBinContent( bini, effB     );
      rejBvsS->SetBinContent( bini, 1.0-effB );
      if (effB>std::numeric_limits<double>::epsilon())
         inveffBvsS->SetBinContent( bini, 1.0/effB );
   }

   // create splines for histogram
   fSpleffBvsS = new TSpline1( "effBvsS", new TGraph( effBvsS ) );

   // search for overlap point where, when cutting on it,
   // one would obtain: eff_S = rej_B = 1 - eff_B

   Double_t effS = 0., rejB = 0., effS_ = 0., rejB_ = 0.;
   Int_t    nbins = 5000;
   for (Int_t bini=1; bini<=nbins; bini++) {

      // get corresponding signal and background efficiencies
      effS = (bini - 0.5)/Float_t(nbins);
      rejB = 1.0 - fSpleffBvsS->Eval( effS );

      // find signal efficiency that corresponds to required background efficiency
      if ((effS - rejB)*(effS_ - rejB_) < 0) break;
      effS_ = effS;
      rejB_ = rejB;
   }
   // find cut that corresponds to signal efficiency and update signal-like criterion
   fSignalCut = Root( 0.5*(effS + effS_) );

   return rejBvsS;
}

////////////////////////////////////////////////////////////////////////////////
/// code to compute the area under the ROC ( rej-vs-eff ) curve

Double_t TMVA::ROCCalc::GetROCIntegral(){
   Double_t effS = 0, effB = 0;
   Int_t    nbins = 1000;
   if (fSpleffBvsS == 0) this->GetROC(); // that will make the ROC calculation if not done yet

   // compute area of rej-vs-eff plot
   Double_t integral = 0;
   for (Int_t bini=1; bini<=nbins; bini++) {

      // get corresponding signal and background efficiencies
      effS = (bini - 0.5)/Float_t(nbins);
      effB = fSpleffBvsS->Eval( effS );
      integral += (1.0 - effB);
   }
   integral /= nbins;

   return integral;
}

////////////////////////////////////////////////////////////////////////////////
/// get the signal efficiency for a particular background efficiency
/// that will be the value of the efficiency retured (does not affect
/// the efficiency-vs-bkg plot which is done anyway.

Double_t TMVA::ROCCalc::GetEffSForEffBof(Double_t effBref, Double_t &effSerr){
   // find precise efficiency value
   Double_t effS=0., effB, effSOld=1., effBOld=0.;
   Int_t    nbins = 1000;
   if (fSpleffBvsS == 0) this->GetROC(); // that will make the ROC calculation if not done yet

   Float_t step=1./nbins;  // stepsize in efficiency binning
   for (Int_t bini=1; bini<=nbins; bini++) {
      // get corresponding signal and background efficiencies
      effS = (bini - 0.5)*step;  // efficiency goes from 0-to-1 in nbins steps of 1/nbins (take middle of the bin)
      effB = fSpleffBvsS->Eval( effS );

      // find signal efficiency that corresponds to required background efficiency
      if ((effB - effBref)*(effBOld - effBref) <= 0) break;
      effSOld = effS;
      effBOld = effB;
   }

   // take mean between bin above and bin below
   effS = 0.5*(effS + effSOld);


   if (fNevtS > 0) effSerr = TMath::Sqrt( effS*(1.0 - effS)/fNevtS );
   else effSerr = 0;

   return effS;
}

////////////////////////////////////////////////////////////////////////////////
/// returns efficiency as function of cut

Double_t TMVA::ROCCalc::GetEffForRoot( Double_t theCut )
{
   Double_t retVal=0;

   // retrieve the class object
   if (fUseSplines) retVal = fSplmvaCumS->Eval( theCut );
   else             retVal = fmvaScumul->GetBinContent( fmvaScumul->FindBin( theCut ) );

   // caution: here we take some "forbidden" action to hide a problem:
   // in some cases, in particular for likelihood, the binned efficiency distributions
   // do not equal 1, at xmin, and 0 at xmax; of course, in principle we have the
   // unbinned information available in the trees, but the unbinned minimization is
   // too slow, and we don't need to do a precision measurement here. Hence, we force
   // this property.
   Double_t eps = 1.0e-5;
   if      (theCut-fXmin < eps) retVal = (fCutOrientation > 0) ? 1.0 : 0.0;
   else if (fXmax-theCut < eps) retVal = (fCutOrientation > 0) ? 0.0 : 1.0;


   return retVal;
}

////////////////////////////////////////////////////////////////////////////////
/// Root finding using Brents algorithm; taken from CERNLIB function RZERO

Double_t TMVA::ROCCalc::Root( Double_t refValue  )
{
   Double_t a  = fXmin, b = fXmax;
   Double_t fa = GetEffForRoot( a ) - refValue;
   Double_t fb = GetEffForRoot( b ) - refValue;
   if (fb*fa > 0) {
      Log() << kWARNING << "<ROCCalc::Root> initial interval w/o root: "
            << "(a=" << a << ", b=" << b << "),"
            << " (Eff_a=" << GetEffForRoot( a )
            << ", Eff_b=" << GetEffForRoot( b ) << "), "
            << "(fa=" << fa << ", fb=" << fb << "), "
            << "refValue = " << refValue << Endl;
      return 1;
   }

   Bool_t   ac_equal(kFALSE);
   Double_t fc = fb;
   Double_t c  = 0, d = 0, e = 0;
   for (Int_t iter= 0; iter <= fMaxIter; iter++) {
      if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) {

         // Rename a,b,c and adjust bounding interval d
         ac_equal = kTRUE;
         c  = a; fc = fa;
         d  = b - a; e  = b - a;
      }

      if (TMath::Abs(fc) < TMath::Abs(fb)) {
         ac_equal = kTRUE;
         a  = b;  b  = c;  c  = a;
         fa = fb; fb = fc; fc = fa;
      }

      Double_t tol = 0.5 * 2.2204460492503131e-16 * TMath::Abs(b);
      Double_t m   = 0.5 * (c - b);
      if (fb == 0 || TMath::Abs(m) <= tol || TMath::Abs(fb) < fAbsTol) return b;

      // Bounds decreasing too slowly: use bisection
      if (TMath::Abs (e) < tol || TMath::Abs (fa) <= TMath::Abs (fb)) { d = m; e = m; }
      else {
         // Attempt inverse cubic interpolation
         Double_t p, q, r;
         Double_t s = fb / fa;

         if (ac_equal) { p = 2 * m * s; q = 1 - s; }
         else {
            q = fa / fc; r = fb / fc;
            p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
            q = (q - 1) * (r - 1) * (s - 1);
         }
         // Check whether we are in bounds
         if (p > 0) q = -q;
         else       p = -p;

         Double_t min1 = 3 * m * q - TMath::Abs (tol * q);
         Double_t min2 = TMath::Abs (e * q);
         if (2 * p < (min1 < min2 ? min1 : min2)) {
            // Accept the interpolation
            e = d;        d = p / q;
         }
         else { d = m; e = m; } // Interpolation failed: use bisection.
      }
      // Move last best guess to a
      a  = b; fa = fb;
      // Evaluate new trial root
      if (TMath::Abs(d) > tol) b += d;
      else                     b += (m > 0 ? +tol : -tol);

      fb = GetEffForRoot( b ) - refValue;

   }

   // Return our best guess if we run out of iterations
   Log() << kWARNING << "<ROCCalc::Root> maximum iterations (" << fMaxIter
         << ") reached before convergence" << Endl;

   return b;
}

////////////////////////////////////////////////////////////////////////////////

TH1* TMVA::ROCCalc::GetPurity( Int_t nStot, Int_t nBtot)
{
   if (fnStot!=nStot || fnBtot!=nBtot || !fSignificance) {
      GetSignificance(nStot, nBtot);
      fnStot=nStot;
      fnBtot=nBtot;
   }
   return fPurity;
}
////////////////////////////////////////////////////////////////////////////////

TH1* TMVA::ROCCalc::GetSignificance( Int_t nStot, Int_t nBtot)
{
   if (fnStot==nStot && fnBtot==nBtot && !fSignificance) return fSignificance;
   fnStot=nStot; fnBtot=nBtot;

   fSignificance = (TH1*) fmvaScumul->Clone("Significance"); fSignificance->SetTitle("Significance");
   fSignificance->Reset(); fSignificance->SetFillStyle(0);
   fSignificance->SetXTitle("mva cut value");
   fSignificance->SetYTitle("Stat. significance S/Sqrt(S+B)");
   fSignificance->SetLineColor(2);
   fSignificance->SetLineWidth(5);

   fPurity = (TH1*) fmvaScumul->Clone("Purity"); fPurity->SetTitle("Purity");
   fPurity->Reset(); fPurity->SetFillStyle(0);
   fPurity->SetXTitle("mva cut value");
   fPurity->SetYTitle("Purity: S/(S+B)");
   fPurity->SetLineColor(3);
   fPurity->SetLineWidth(5);

   Double_t maxSig=0;
   for (Int_t i=1; i<=fSignificance->GetNbinsX(); i++) {
      Double_t S = fmvaScumul->GetBinContent( i ) * nStot;
      Double_t B = fmvaBcumul->GetBinContent( i ) * nBtot;
      Double_t purity;
      Double_t sig;
      if (S+B > 0){
         purity = S/(S+B);
         sig = S/TMath::Sqrt(S+B);
         if (sig > maxSig) {
            maxSig    = sig;
         }
      } else {
         purity=0;
         sig=0;
      }
      cout << "S="<<S<<" B="<<B<< " purity="<<purity<< endl;
      fPurity->SetBinContent( i, purity );
      fSignificance->SetBinContent( i, sig );
   }

   /*
        TLatex*  line1;
        TLatex*  line2;
        TLatex tl;
        tl.SetNDC();
        tl.SetTextSize( 0.033 );
        Int_t maxbin = fSignificance->GetMaximumBin();
        line1 = tl.DrawLatex( 0.15, 0.23, Form("For %1.0f signal and %1.0f background", nStot, nBtot));
        tl.DrawLatex( 0.15, 0.19, "events the maximum S/Sqrt(S+B) is");

        line2 = tl.DrawLatex( 0.15, 0.15, Form("%4.2f when cutting at %5.2f",
        maxSig,
        fSignificance->GetXaxis()->GetBinCenter(maxbin)) );
   */
   return fSignificance;

}
