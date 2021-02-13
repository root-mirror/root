// @(#)root/gpad:$Id$
// Author: Paul Gessinger   25/08/2016

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TRatioPlot.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TF1.h"
#include "TPad.h"
#include "TString.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "THStack.h"

/** \class TRatioPlot
    \ingroup gpad
Class for displaying ratios, differences and fit residuals.

TRatioPlot has two constructors, one which accepts two histograms, and is responsible
for setting up the calculation of ratios and differences. This calculation is in part
delegated to `TEfficiency`. A single option can be given as a parameter, that is
used to determine which procedure is chosen. The remaining option string is then
passed through to the calculation, if applicable. The other constructor uses a
fitted histogram to calculate the fit residual and plot it with the histogram
and the fit function.

## Ratios and differences
The simplest case is passing two histograms without specifying any options. This defaults to using
`TGraphAsymmErrors::Divide`. The `option` variable is passed through, as are the parameters
`c1` and `c2`, that you can set via `TRatioPlot::SetC1` and `TRatioPlot::SetC1`. If you set the
`option` to `divsym` the method `TH1::Divide` will be used instead, also receiving all the parameters.

Using the `option` `diff` or `diffsig`, both histograms will be subtracted, and in the case of diffsig,
the difference will be divided by the  uncertainty. `c1` and `c2` will only be used to
scale the histograms using `TH1::Scale` prior to subtraction.

Available options are for `option`:
| Option     | Description                                                  |
| ---------- | ------------------------------------------------------------ |
| divsym    | uses the histogram `TH1::Divide` method, yields symmetric errors    |
| diff       | subtracts the histograms                                     |
| diffsig    | subtracts the histograms and divides by the uncertainty |

Begin_Macro(source)
../../../tutorials/hist/ratioplot1.C
End_Macro

## Fit residuals
A second constructor only accepts a single histogram, but expects it to have a fitted
function. The function is used to calculate the residual between the fit and the
histogram. Here, it is expected that h1 has a fit function in it's list of functions. The class calculates the
difference between the histogram and the fit function at each point and divides it by the uncertainty. There
are a few option to steer which error is used (as is the case for `diffsig`). The default is to use
the statistical uncertainty from h1 using `TH1::GetBinError`. If the `option` string contains `errasym`, asymmetric
errors will be used. The type of error can be steered by `TH1::SetBinErrorOption`. The corresponding error will be used,
depending on if the function is below or above the bin content. The third option `errfunc` uses the square root of
the function value as the error.


Begin_Macro(source)
../../../tutorials/hist/ratioplot2.C
End_Macro

## Error options for difference divided by uncertainty and fit residual
The uncertainty that is used in the calculation can be steered by providing
options to the `option` argument.

| Option     | Description                                                  |
| ---------- | ------------------------------------------------------------ |
| errasym    | Uses calculated asymmetric errors from `TH1::GetBinErrorUp`/`TH1::GetBinErrorLow`. Note that you need to set `TH1::SetBinErrorOption` first |
| errfunc    | Uses \f$ \sqrt{f(x)} \f$ as the error |

The asymmetric error case uses the upper or lower error depending on the relative size
of the bin contents, or the bin content and the function value.

## Access to internal parts
You can access the internal objects that are used to construct the plot via a series of
methods. `TRatioPlot::GetUpperPad` and `TRatioPlot::GetLowerPad` can be used to draw additional
elements on top of the existing ones.
`TRatioPlot::GetLowerRefGraph` returns a reference to the lower pad's graph that
is responsible for the range, which enables you to modify the range.

\image html gpad_ratioplot.png
*/

////////////////////////////////////////////////////////////////////////////////
/// TRatioPlot default constructor

TRatioPlot::TRatioPlot()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TRatioPlot::~TRatioPlot()
{

   gROOT->GetListOfCleanups()->Remove(this);

   if (fRatioGraph && fRatioGraph->TestBit(kNotDeleted)) delete fRatioGraph;
   if (fConfidenceInterval1 && fConfidenceInterval1->TestBit(kNotDeleted)) delete fConfidenceInterval1;
   if (fConfidenceInterval2 && fConfidenceInterval2->TestBit(kNotDeleted)) delete fConfidenceInterval2;

   for (unsigned int i=0;i<fGridlines.size();++i) {
      delete (fGridlines[i]);
   }

   if (fSharedXAxis && fSharedXAxis->TestBit(kNotDeleted)) delete fSharedXAxis;
   if (fUpperGXaxis && fUpperGXaxis->TestBit(kNotDeleted)) delete fUpperGXaxis;
   if (fLowerGXaxis && fLowerGXaxis->TestBit(kNotDeleted)) delete fLowerGXaxis;
   if (fUpperGYaxis && fUpperGYaxis->TestBit(kNotDeleted)) delete fUpperGYaxis;
   if (fLowerGYaxis && fLowerGYaxis->TestBit(kNotDeleted)) delete fLowerGYaxis;
   if (fUpperGXaxisMirror && fUpperGXaxisMirror->TestBit(kNotDeleted)) delete fUpperGXaxisMirror;
   if (fLowerGXaxisMirror && fLowerGXaxisMirror->TestBit(kNotDeleted)) delete fLowerGXaxisMirror;
   if (fUpperGYaxisMirror && fUpperGYaxisMirror->TestBit(kNotDeleted)) delete fUpperGYaxisMirror;
   if (fLowerGYaxisMirror && fLowerGYaxisMirror->TestBit(kNotDeleted)) delete fLowerGYaxisMirror;

   if (fUpYaxis && fUpYaxis->TestBit(kNotDeleted)) delete fUpYaxis;
   if (fLowYaxis && fLowYaxis->TestBit(kNotDeleted)) delete fLowYaxis;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal method that shares constructor logic

void TRatioPlot::Init(TH1* h1, TH1* h2,Option_t *option)
{

   fH1 = h1;
   fH2 = h2;

   SetupPads();

   TString optionString = TString(option);

   if (optionString.Contains("divsym")) {
      optionString.ReplaceAll("divsym", "");
      fMode = TRatioPlot::CalculationMode::kDivideHist;
   } else if (optionString.Contains("diffsig")) {
      optionString.ReplaceAll("diffsig", "");
      fMode = TRatioPlot::CalculationMode::kDifferenceSign;

      // determine which error style
      if (optionString.Contains("errasym")) {
         fErrorMode = TRatioPlot::ErrorMode::kErrorAsymmetric;
         optionString.ReplaceAll("errasym", "");
      }

      if (optionString.Contains("errfunc")) {
         fErrorMode = TRatioPlot::ErrorMode::kErrorFunc;
         optionString.ReplaceAll("errfunc", "");
      }
   } else if (optionString.Contains("diff")) {
      optionString.ReplaceAll("diff", "");
      fMode = TRatioPlot::CalculationMode::kDifference;
   } else {
      fMode = TRatioPlot::CalculationMode::kDivideGraph; // <- default
   }

   fOption = optionString;


   fH1DrawOpt = "hist";
   fH2DrawOpt = "E";
   fGraphDrawOpt = "AP";


   // build ratio, everything is ready
   if (!BuildLowerPlot()) return;

   // taking x axis information from h1 by cloning it x axis
   fSharedXAxis = (TAxis*)(fH1->GetXaxis()->Clone());
   fUpYaxis = (TAxis*)(fH1->GetYaxis()->Clone());
   fLowYaxis = (TAxis*)(fRatioGraph->GetYaxis()->Clone());
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor for two histograms
///
/// \param h1 First histogram
/// \param h2 Second histogram
/// \param option Steers the error calculation, as well as ratio / difference

TRatioPlot::TRatioPlot(TH1* h1, TH1* h2, Option_t *option)
   : fGridlines()
{
   gROOT->GetListOfCleanups()->Add(this);

   if (!h1 || !h2) {
      Warning("TRatioPlot", "Need two histograms.");
      return;
   }

   Bool_t h1IsTH1=h1->IsA()->InheritsFrom(TH1::Class());
   Bool_t h2IsTH1=h2->IsA()->InheritsFrom(TH1::Class());

   if (!h1IsTH1 && !h2IsTH1) {
      Warning("TRatioPlot", "Need two histograms deriving from TH2 or TH3.");
      return;
   }

   fHistDrawProxy = h1;

   Init(h1, h2, option);

}

////////////////////////////////////////////////////////////////////////////////
/// Constructor which accepts a `THStack` and a histogram. Converts the
/// stack to a regular sum of its containing histograms for processing.
///
/// \param st The THStack object
/// \param h2 The other histogram
/// \param option Steers the calculation of the lower plot

TRatioPlot::TRatioPlot(THStack* st, TH1* h2, Option_t *option)
{
   if (!st || !h2) {
      Warning("TRatioPlot", "Need a histogram and a stack");
      return;
   }

   TList *stackHists = st->GetHists();

   if (stackHists->GetSize() == 0) {
      Warning("TRatioPlot", "Stack does not have histograms");
      return;
   }

   TH1* tmpHist = (TH1*)stackHists->At(0)->Clone();
   tmpHist->Reset();

   for (int i=0;i<stackHists->GetSize();++i) {
      tmpHist->Add((TH1*)stackHists->At(i));
   }

   fHistDrawProxy = st;

   Init(tmpHist, h2, option);

}

////////////////////////////////////////////////////////////////////////////////
/// Constructor for one histogram and a fit.
///
/// \param h1 The histogram
/// \param option Steers the error calculation
/// \param fitres Explicit fit result to be used for calculation. Uses last fit if left empty

TRatioPlot::TRatioPlot(TH1* h1, Option_t *option, TFitResult *fitres)
   : fH1(h1),
     fGridlines()
{
   gROOT->GetListOfCleanups()->Add(this);

   if (!fH1) {
      Warning("TRatioPlot", "Need a histogram.");
      return;
   }

   Bool_t h1IsTH1=fH1->IsA()->InheritsFrom(TH1::Class());

   if (!h1IsTH1) {
      Warning("TRatioPlot", "Need a histogram deriving from TH2 or TH3.");
      return;
   }

   TList *h1Functions = fH1->GetListOfFunctions();

   if (h1Functions->GetSize() < 1) {
      Warning("TRatioPlot", "Histogram given needs to have a (fit) function associated with it");
      return;
   }


   fHistDrawProxy = h1;

   fFitResult = fitres;

   fMode = TRatioPlot::CalculationMode::kFitResidual;

   TString optionString = TString(option);

   // determine which error style
   if (optionString.Contains("errasym")) {
      fErrorMode = TRatioPlot::ErrorMode::kErrorAsymmetric;
      optionString.ReplaceAll("errasym", "");
   }

   if (optionString.Contains("errfunc")) {
      fErrorMode = TRatioPlot::ErrorMode::kErrorFunc;
      optionString.ReplaceAll("errfunc", "");
   }

   fOption = optionString;

   if (!BuildLowerPlot()) return;

   // emulate option behaviour of TH1
   if (fH1->GetSumw2N() > 0) {
      fH1DrawOpt = "E";
   } else {
      fH1DrawOpt = "hist";
   }
   fGraphDrawOpt = "LX"; // <- default

   fSharedXAxis = (TAxis*)(fH1->GetXaxis()->Clone());
   fUpYaxis     = (TAxis*)(fH1->GetYaxis()->Clone());
   fLowYaxis    = (TAxis*)(fRatioGraph->GetYaxis()->Clone());

   //SyncAxesRanges();

   SetupPads();

}

////////////////////////////////////////////////////////////////////////////////
/// Sets the drawing option for h1

void TRatioPlot::SetH1DrawOpt(Option_t *opt)
{
   fH1DrawOpt = opt;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the drawing option for h2

void TRatioPlot::SetH2DrawOpt(Option_t *opt)
{
   TString optString = TString(opt);
   optString.ReplaceAll("same", "");
   optString.ReplaceAll("SAME", "");

   fH2DrawOpt = optString;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the drawing option for the lower graph

void TRatioPlot::SetGraphDrawOpt(Option_t *opt)
{
   fGraphDrawOpt = opt;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the drawing option for the fit in the fit residual case

void TRatioPlot::SetFitDrawOpt(Option_t *opt)
{
   fFitDrawOpt = opt;
}

////////////////////////////////////////////////////////////////////////////////
/// Setup the pads.

void TRatioPlot::SetupPads() {

   // this method will delete all the pads before recreating them

   if (fUpperPad != 0) {
      delete fUpperPad;
      fUpperPad = 0;
   }

   if (fLowerPad != 0) {
      delete fLowerPad;
      fLowerPad = 0;
   }

   if (!gPad) {
      Error("SetupPads", "need to create a canvas first");
      return;
   }

   double pm = fInsetWidth;
   double width = gPad->GetWNDC();
   double height = gPad->GetHNDC();
   double f = height/width;

   fUpperPad = new TPad("upper_pad", "", pm*f, fSplitFraction, 1.-pm*f, 1.-pm);
   fLowerPad = new TPad("lower_pad", "", pm*f, pm, 1.-pm*f, fSplitFraction);

   SetPadMargins();

   // connect to the pads signal
   fUpperPad->Connect("RangeAxisChanged()", "TRatioPlot", this, "RangeAxisChanged()");
   fLowerPad->Connect("RangeAxisChanged()", "TRatioPlot", this, "RangeAxisChanged()");

   fUpperPad->Connect("UnZoomed()", "TRatioPlot", this, "UnZoomed()");
   fLowerPad->Connect("UnZoomed()", "TRatioPlot", this, "UnZoomed()");

   fUpperPad->Connect("Resized()", "TRatioPlot", this, "SubPadResized()");
   fLowerPad->Connect("Resized()", "TRatioPlot", this, "SubPadResized()");

   if (fTopPad != 0) {
      delete fTopPad;
      fTopPad = 0;
   }

   fTopPad = new TPad("top_pad", "", pm*f, pm, 1-pm*f, 1-pm);

   fTopPad->SetBit(kCannotPick);

}

////////////////////////////////////////////////////////////////////////////////
/// Browse.

void TRatioPlot::Browse(TBrowser *b)
{
   Draw(b ? b->GetDrawOption() : "");
   gPad->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the top margin of the upper pad.
///
/// \param margin The new margin

void TRatioPlot::SetUpTopMargin(Float_t margin)
{
   fUpTopMargin = margin;
   SetPadMargins();
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the bottom margin of the upper pad.
///
/// \param margin The new margin

void TRatioPlot::SetUpBottomMargin(Float_t margin)
{
   fUpBottomMargin = margin;
   SetPadMargins();
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the top margin of the lower pad.
///
/// \param margin The new margin

void TRatioPlot::SetLowTopMargin(Float_t margin)
{
   fLowTopMargin = margin;
   SetPadMargins();
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the bottom margin of the lower pad.
///
/// \param margin The new margin

void TRatioPlot::SetLowBottomMargin(Float_t margin)
{
   fLowBottomMargin = margin;
   SetPadMargins();
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the left margin of both pads.
/// \param margin The new margin

void TRatioPlot::SetLeftMargin(Float_t margin)
{
   fLeftMargin = margin;
   SetPadMargins();
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the right margin of both pads.
///
/// \param margin The new margin

void TRatioPlot::SetRightMargin(Float_t margin)
{
   fRightMargin = margin;
   SetPadMargins();
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the margin that separates the two pads. The margin is split according
/// to the relative sizes of the pads
///
/// \param margin The new margin
///
/// Begin_Macro(source)
/// ../../../tutorials/hist/ratioplot6.C
/// End_Macro

void TRatioPlot::SetSeparationMargin(Float_t margin)
{
   Float_t sf = fSplitFraction;
   fUpBottomMargin = margin/2./(1-sf);
   fLowTopMargin = margin/2./sf;
   SetPadMargins();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the separation margin value.

Float_t TRatioPlot::GetSeparationMargin() const
{
   Float_t sf = fSplitFraction;
   Float_t up = fUpBottomMargin * (1-sf);
   Float_t down = fLowTopMargin * sf;
   return up+down;
}

////////////////////////////////////////////////////////////////////////////////
/// Draws the ratio plot to the currently active pad. Therefore it requires that
/// a TCanvas has been created first.
///
/// It takes the following options
///
/// | Option     | Description                                                  |
/// | ---------- | ------------------------------------------------------------ |
/// | grid / nogrid | enable (default) or disable drawing of dashed lines on lower plot |
/// | hideup     | hides the first label of the upper axis if there is not enough space |
/// | fhideup    | always hides the first label of the upper axis |
/// | hidelow (default) | hides the last label of the lower axis if there is not enough space |
/// | fhidelow   | always hides the last label of the lower axis |
/// | nohide     | does not hide a label if there is not enough space |
/// | noconfint  | does not draw the confidence interval bands in the fit residual case |
/// | confint    | draws the confidence interval bands in the fit residual case (default) |

void TRatioPlot::Draw(Option_t *option)
{

   TString drawOpt = option;

   if (drawOpt.Contains("nogrid")) {
      drawOpt.ReplaceAll("nogrid", "");
      fShowGridlines = kFALSE;
   } else if (drawOpt.Contains("grid")) {
      drawOpt.ReplaceAll("grid", "");
      fShowGridlines = kTRUE;
   }

   if (drawOpt.Contains("noconfint")) {
      drawOpt.ReplaceAll("noconfint", "");
      fShowConfidenceIntervals = kFALSE;
   } else if (drawOpt.Contains("confint")) {
      drawOpt.ReplaceAll("confint", "");
      fShowConfidenceIntervals = kTRUE; // <- default
   }

   if (drawOpt.Contains("fhideup")) {
      fHideLabelMode = TRatioPlot::HideLabelMode::kForceHideUp;
   } else if (drawOpt.Contains("fhidelow")) {
      fHideLabelMode = TRatioPlot::HideLabelMode::kForceHideLow;
   } else if (drawOpt.Contains("hideup")) {
      fHideLabelMode = TRatioPlot::HideLabelMode::kHideUp;
   } else if (drawOpt.Contains("hidelow")) {
      fHideLabelMode = TRatioPlot::HideLabelMode::kHideLow;
   } else if (drawOpt.Contains("nohide")) {
      fHideLabelMode = TRatioPlot::HideLabelMode::kNoHide;
   } else {
      fHideLabelMode = TRatioPlot::HideLabelMode::kHideLow; // <- default
   }

   if (!gPad) {
      Error("Draw", "need to create a canvas first");
      return;
   }

   TVirtualPad *padsav = gPad;
   fParentPad = gPad;

   fUpperPad->SetLogy(fParentPad->GetLogy());
   fUpperPad->SetLogx(fParentPad->GetLogx());
   fLowerPad->SetLogx(fParentPad->GetLogx());

   fUpperPad->SetGridx(fParentPad->GetGridx());
   fUpperPad->SetGridy(fParentPad->GetGridy());
   fLowerPad->SetGridx(fParentPad->GetGridx());
   fLowerPad->SetGridy(fParentPad->GetGridy());

   // we are a TPad

   fUpperPad->Draw();
   fLowerPad->Draw();

   fTopPad->SetFillStyle(0);
   fTopPad->Draw();

   fUpperPad->cd();

   fConfidenceInterval2->SetFillColor(fCi1Color);
   fConfidenceInterval1->SetFillColor(fCi2Color);

   if (fMode == TRatioPlot::CalculationMode::kFitResidual) {
      TF1 *func = dynamic_cast<TF1*>(fH1->GetListOfFunctions()->At(0));

      if (func == 0) {
         // this is checked in constructor and should thus not occur
         Error("BuildLowerPlot", "h1 does not have a fit function");
         return;
      }

      fH1->Draw("A"+fH1DrawOpt);
      func->Draw(fFitDrawOpt+"same");

      fLowerPad->cd();

      if (fShowConfidenceIntervals) {
         fConfidenceInterval2->Draw("IA3");
         fConfidenceInterval1->Draw("3");
         fRatioGraph->Draw(fGraphDrawOpt+"SAME");
      } else {
         fRatioGraph->Draw("IA"+fGraphDrawOpt+"SAME");
      }
   } else {

      if (fHistDrawProxy) {
         if (fHistDrawProxy->InheritsFrom(TH1::Class())) {
            ((TH1*)fHistDrawProxy)->Draw("A"+fH1DrawOpt);
         } else if (fHistDrawProxy->InheritsFrom(THStack::Class())) {
            ((THStack*)fHistDrawProxy)->Draw("A"+fH1DrawOpt);
         } else {
            Warning("Draw", "Draw proxy not of type TH1 or THStack, not drawing it");
         }
      }

      fH2->Draw("A"+fH2DrawOpt+"same");

      fLowerPad->cd();

      TString opt = fGraphDrawOpt;
      fRatioGraph->Draw("IA"+fGraphDrawOpt);

   }

   // assign same axis ranges to lower pad as in upper pad
   // the visual axes will be created on paint
   SyncAxesRanges();

   CreateGridline();

   padsav->cd();
   AppendPad();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the reference graph for the lower pad, which means the graph that
/// is responsible for setting the coordinate system. It is the first graph
/// added to the primitive list of the lower pad.
/// This reference can be used to set the minimum and maximum of the lower pad.
/// Note that `TRatioPlot::Draw` needs to have been called first, since the
/// graphs are only created then.
///
/// Begin_Macro(source)
/// ../../../tutorials/hist/ratioplot3.C
/// End_Macro

TGraph* TRatioPlot::GetLowerRefGraph() const
{
   if (fLowerPad == 0) {
      Error("GetLowerRefGraph", "Lower pad has not been defined");
      return 0;
   }

   TList *primlist = fLowerPad->GetListOfPrimitives();
   if (primlist->GetSize() == 0) {
      Error("GetLowerRefGraph", "Lower pad does not have primitives");
      return 0;
   }

   TObjLink *lnk = primlist->FirstLink();

   while (lnk) {
      TObject *obj = lnk->GetObject();

      if (obj->InheritsFrom(TGraph::Class())) {
         return (TGraph*)obj;
      }

      lnk = lnk->Next();
   }

   Error("GetLowerRefGraph", "Did not find graph in list");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the reference object. Its the first TH1 or THStack type object
/// in the upper pads list of primitives.
/// Note that it returns a `TObject`, so you need to test and cast it to use it.

TObject* TRatioPlot::GetUpperRefObject() const
{
   TList *primlist = fUpperPad->GetListOfPrimitives();
   TObject *refobj = 0;
   for (Int_t i=0;i<primlist->GetSize();++i) {
      refobj = primlist->At(i);
      if (refobj->InheritsFrom(TH1::Class()) || refobj->InheritsFrom(THStack::Class())) {
         return refobj;
      }
   }

   Error("GetUpperRefObject", "No upper ref object of TH1 or THStack type found");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Gets the x axis of the object returned by `TRatioPlot::GetUpperRefObject`.

TAxis* TRatioPlot::GetUpperRefXaxis() const
{
   TObject *refobj = GetUpperRefObject();

   if (!refobj) return 0;

   if (refobj->InheritsFrom(TH1::Class())) {
      return ((TH1*)refobj)->GetXaxis();
   } else if (refobj->InheritsFrom(THStack::Class())) {
      return ((THStack*)refobj)->GetXaxis();
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Gets the y axis of the object returned by `TRatioPlot::GetUpperRefObject`.

TAxis* TRatioPlot::GetUpperRefYaxis() const
{
   TObject *refobj = GetUpperRefObject();

   if (!refobj) return 0;

   if (refobj->InheritsFrom(TH1::Class())) {
      return ((TH1*)refobj)->GetYaxis();
   } else if (refobj->InheritsFrom(THStack::Class())) {
      return ((THStack*)refobj)->GetYaxis();
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a grid line

void TRatioPlot::CreateGridline()
{

   if (!fShowGridlines) {
      return; // don't draw them
   }

   TVirtualPad *padsav = gPad;

   fLowerPad->cd();

   unsigned int dest = fGridlinePositions.size();

   Double_t lowYFirst = fLowerPad->GetUymin();
   Double_t lowYLast = fLowerPad->GetUymax();

   double y;
   int outofrange = 0;
   for (unsigned int i=0;i<fGridlinePositions.size();++i) {
      y = fGridlinePositions.at(i);

      if (y < lowYFirst || lowYLast < y) {
         ++outofrange;
      }

   }

   dest = dest - outofrange;

   // clear all
   for (unsigned int i=0;i<fGridlines.size();++i) {
      delete fGridlines.at(i);
   }

   fGridlines.erase(fGridlines.begin(), fGridlines.end());

   for (unsigned int i=0;i<dest;++i) {
      TLine *newline = new TLine(0, 0, 0, 0);
      newline->SetLineStyle(2);
      newline->Draw();
      fGridlines.push_back(newline);
   }

   Double_t first = fSharedXAxis->GetBinLowEdge(fSharedXAxis->GetFirst());
   Double_t last = fSharedXAxis->GetBinUpEdge(fSharedXAxis->GetLast());

   TLine *line;
   unsigned int skipped = 0;
   for (unsigned int i=0;i<fGridlinePositions.size();++i) {
      y = fGridlinePositions[i];

      if (y < lowYFirst || lowYLast < y) {
         // this is one of the ones that was out of range
         ++skipped;
         continue;
      }

      line = fGridlines.at(i-skipped);

      line->SetX1(first);
      line->SetX2(last);
      line->SetY1(y);
      line->SetY2(y);
   }

   padsav->cd();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates the visual axes when painting.

void TRatioPlot::Paint(Option_t * /*opt*/)
{
   // create the visual axes
   CreateVisualAxes();
   CreateGridline();

   if (fIsUpdating) fIsUpdating = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Syncs the axes ranges from the shared ones to the actual ones.

void TRatioPlot::SyncAxesRanges()
{
   // get ranges from the shared axis clone
   Double_t first = fSharedXAxis->GetBinLowEdge(fSharedXAxis->GetFirst());
   Double_t last = fSharedXAxis->GetBinUpEdge(fSharedXAxis->GetLast());

   // set range on computed graph, have to set it twice because
   // TGraph's axis looks strange otherwise
   TAxis *ref = GetLowerRefXaxis();
   ref->SetLimits(first, last);
   ref->SetRangeUser(first, last);

   GetUpperRefXaxis()->SetRangeUser(first, last);

}

////////////////////////////////////////////////////////////////////////////////
/// Build the lower plot according to which constructor was called, and
/// which options were passed.

Int_t TRatioPlot::BuildLowerPlot()
{
   // Clear and delete the graph if not exists
   if (fRatioGraph != 0) {
      fRatioGraph->IsA()->Destructor(fRatioGraph);
      fRatioGraph = 0;
   }

   if (fConfidenceInterval1 == 0) {
      fConfidenceInterval1 = new TGraphErrors();
   }

   if (fConfidenceInterval2 == 0) {
      fConfidenceInterval2 = new TGraphErrors();
   }

   static Double_t divideGridlines[] = {0.7, 1.0, 1.3};
   static Double_t diffGridlines[] = {0.0};
   static Double_t signGridlines[] = {1.0, 0.0, -1.0};

   // Determine the divide mode and create the lower graph accordingly
   // Pass divide options given in constructor
   if (fMode == TRatioPlot::CalculationMode::kDivideGraph) {
      // use TGraphAsymmErrors Divide method to create

      SetGridlines(divideGridlines, 3);

      TH1 *tmpH1 = (TH1*)fH1->Clone();
      TH1 *tmpH2 = (TH1*)fH2->Clone();

      tmpH1->Scale(fC1);
      tmpH2->Scale(fC2);

      TGraphAsymmErrors *ratioGraph = new TGraphAsymmErrors();
      ratioGraph->Divide(tmpH1, tmpH2, fOption.Data());
      fRatioGraph = ratioGraph;

      delete tmpH1;
      delete tmpH2;

   } else if (fMode == TRatioPlot::CalculationMode::kDifference) {
      SetGridlines(diffGridlines, 3);

      TH1 *tmpHist = (TH1*)fH1->Clone();

      tmpHist->Reset();

      tmpHist->Add(fH1, fH2, fC1, -1*fC2);
      fRatioGraph = new TGraphErrors(tmpHist);

      delete tmpHist;
   } else if (fMode == TRatioPlot::CalculationMode::kDifferenceSign) {

      SetGridlines(signGridlines, 3);

      fRatioGraph = new TGraphAsymmErrors();
      Int_t ipoint = 0;
      Double_t res;
      Double_t error;

      Double_t val;
      Double_t val2;

      for (Int_t i=0; i<=fH1->GetNbinsX();++i) {
         val = fH1->GetBinContent(i);
         val2 = fH2->GetBinContent(i);

         if (fErrorMode == TRatioPlot::ErrorMode::kErrorAsymmetric) {

            Double_t errUp = fH1->GetBinErrorUp(i);
            Double_t errLow = fH1->GetBinErrorLow(i);

            if (val - val2 > 0) {
               // h1 > h2
               error = errLow;
            } else {
               // h1 < h2
               error = errUp;
            }

         } else if (fErrorMode == TRatioPlot::ErrorMode::kErrorSymmetric) {
            error = fH1->GetBinError(i);
         } else {
            Warning("BuildLowerPlot", "error mode is invalid");
            error = 0;
         }

         if (error != 0) {

            res = (val - val2) / error;

            ((TGraphAsymmErrors*)fRatioGraph)->SetPoint(ipoint, fH1->GetBinCenter(i), res);
            ((TGraphAsymmErrors*)fRatioGraph)->SetPointError(ipoint,  fH1->GetBinWidth(i)/2., fH1->GetBinWidth(i)/2., 0.5, 0.5);

            ++ipoint;

         }
      }

   } else if (fMode == TRatioPlot::CalculationMode::kFitResidual) {

      SetGridlines(signGridlines, 3);

      TF1 *func = dynamic_cast<TF1*>(fH1->GetListOfFunctions()->At(0));

      if (func == 0) {
         // this is checked in constructor and should thus not occur
         Error("BuildLowerPlot", "h1 does not have a fit function");
         return 0;
      }

      fRatioGraph = new TGraphAsymmErrors();
      Int_t ipoint = 0;

      Double_t res;
      Double_t error;

      std::vector<double> ci1;
      std::vector<double> ci2;

      Double_t *x_arr = new Double_t[fH1->GetNbinsX()];
      std::fill_n(x_arr, fH1->GetNbinsX(), 0);
      Double_t *ci_arr1 = new Double_t[fH1->GetNbinsX()];
      std::fill_n(ci_arr1, fH1->GetNbinsX(), 0);
      Double_t *ci_arr2 = new Double_t[fH1->GetNbinsX()];
      std::fill_n(ci_arr2, fH1->GetNbinsX(), 0);
      for (Int_t i=0; i<fH1->GetNbinsX();++i) {
         x_arr[i] = fH1->GetBinCenter(i+1);
      }

      Double_t cl1 = fCl1;
      Double_t cl2 = fCl2;

      if (fFitResult != 0) {
         // use this to get conf int

         fFitResult->GetConfidenceIntervals(fH1->GetNbinsX(), 1, 1, x_arr, ci_arr1, cl1);
         for (Int_t i=1; i<=fH1->GetNbinsX();++i) {
            ci1.push_back(ci_arr1[i-1]);
         }

         fFitResult->GetConfidenceIntervals(fH1->GetNbinsX(), 1, 1, x_arr, ci_arr2, cl2);
         for (Int_t i=1; i<=fH1->GetNbinsX();++i) {
            ci2.push_back(ci_arr2[i-1]);
         }
      } else {
         (TVirtualFitter::GetFitter())->GetConfidenceIntervals(fH1->GetNbinsX(), 1, x_arr, ci_arr1, cl1);
         for (Int_t i=1; i<=fH1->GetNbinsX();++i) {
            ci1.push_back(ci_arr1[i-1]);
         }
         (TVirtualFitter::GetFitter())->GetConfidenceIntervals(fH1->GetNbinsX(), 1, x_arr, ci_arr2, cl2);
         for (Int_t i=1; i<=fH1->GetNbinsX();++i) {
            ci2.push_back(ci_arr2[i-1]);
         }

      }

      Double_t x;
      Double_t val;

      for (Int_t i=0; i<=fH1->GetNbinsX();++i) {
         val = fH1->GetBinContent(i);
         x = fH1->GetBinCenter(i+1);

         if (fErrorMode == TRatioPlot::ErrorMode::kErrorAsymmetric) {

            Double_t errUp = fH1->GetBinErrorUp(i);
            Double_t errLow = fH1->GetBinErrorLow(i);

            if (val - func->Eval(fH1->GetBinCenter(i)) > 0) {
               // h1 > fit
               error = errLow;
            } else {
               // h1 < fit
               error = errUp;
            }

         } else if (fErrorMode == TRatioPlot::ErrorMode::kErrorSymmetric) {
            error = fH1->GetBinError(i);
         } else if (fErrorMode == TRatioPlot::ErrorMode::kErrorFunc) {

            error = sqrt(func->Eval(x));

         } else {
            Warning("BuildLowerPlot", "error mode is invalid");
            error = 0;
         }

         if (error != 0) {

            res = (fH1->GetBinContent(i)- func->Eval(fH1->GetBinCenter(i) ) ) / error;
            //__("x="<< x << " y=" << res << " err=" << error);

            ((TGraphAsymmErrors*)fRatioGraph)->SetPoint(ipoint, fH1->GetBinCenter(i), res);
            ((TGraphAsymmErrors*)fRatioGraph)->SetPointError(ipoint,  fH1->GetBinWidth(i)/2., fH1->GetBinWidth(i)/2., 0.5, 0.5);

            fConfidenceInterval1->SetPoint(ipoint, x, 0);
            fConfidenceInterval1->SetPointError(ipoint, x, i < (Int_t)ci1.size() ? ci1[i] / error : 0);
            fConfidenceInterval2->SetPoint(ipoint, x, 0);
            fConfidenceInterval2->SetPointError(ipoint, x, i < (Int_t)ci2.size() ? ci2[i] / error : 0);

            ++ipoint;

         }

      }
      delete [] x_arr;
      delete [] ci_arr1;
      delete [] ci_arr2;
   } else if (fMode == TRatioPlot::CalculationMode::kDivideHist){
      SetGridlines(divideGridlines, 3);

      // Use TH1's Divide method
      TH1 *tmpHist = (TH1*)fH1->Clone();
      tmpHist->Reset();

      tmpHist->Divide(fH1, fH2, fC1, fC2, fOption.Data());
      fRatioGraph = new TGraphErrors(tmpHist);

      delete tmpHist;
   } else {
      // this should not occur
      Error("BuildLowerPlot", "Invalid fMode value");
      return 0;
   }

   // need to set back to "" since recreation. we don't ever want
   // title on lower graph

   if (fRatioGraph == 0) {
      Error("BuildLowerPlot", "Error creating lower graph");
      return 0;
   }

   fRatioGraph->SetTitle("");
   fConfidenceInterval1->SetTitle("");
   fConfidenceInterval2->SetTitle("");

   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// (Re-)Creates the TGAxis objects that are used for consistent display of the
/// axes.

void TRatioPlot::CreateVisualAxes()
{
   TVirtualPad *padsav = gPad;
   fTopPad->cd();

   // this is for errors
   TString thisfunc = "CreateVisualAxes";

   // figure out where the axis has to go.
   // Implicit assumption is, that the top pad spans the full other pads
   Double_t upTM = fUpperPad->GetTopMargin();
   Double_t upBM = fUpperPad->GetBottomMargin();
   Double_t upLM = fUpperPad->GetLeftMargin();
   Double_t upRM = fUpperPad->GetRightMargin();

   Double_t lowTM = fLowerPad->GetTopMargin();
   Double_t lowBM = fLowerPad->GetBottomMargin();
   Double_t lowLM = fLowerPad->GetLeftMargin();
   Double_t lowRM = fLowerPad->GetRightMargin();

   Double_t first = fSharedXAxis->GetBinLowEdge(fSharedXAxis->GetFirst());
   Double_t last = fSharedXAxis->GetBinUpEdge(fSharedXAxis->GetLast());

   Double_t upYFirst = fUpperPad->GetUymin();
   Double_t upYLast = fUpperPad->GetUymax();
   Double_t lowYFirst = fLowerPad->GetUymin();
   Double_t lowYLast = fLowerPad->GetUymax();

   Float_t sf = fSplitFraction;

   // check if gPad has the all sides axis set
   Bool_t mirroredAxes = fParentPad->GetFrameFillStyle() == 0;
   Bool_t axistop = fParentPad->GetTickx() == 1 || mirroredAxes;
   Bool_t axisright = fParentPad->GetTicky() == 1 || mirroredAxes;

   Bool_t logx = fUpperPad->GetLogx() || fLowerPad->GetLogx();
   Bool_t uplogy = fUpperPad->GetLogy();
   Bool_t lowlogy = fLowerPad->GetLogy();

   if (uplogy) {

      upYFirst = TMath::Power(10, upYFirst);
      upYLast = TMath::Power(10, upYLast);

      if (upYFirst <= 0 || upYLast <= 0) {
         Error(thisfunc, "Cannot set upper Y axis to log scale");
      }
   }

   if (lowlogy) {
      lowYFirst = TMath::Power(10, lowYFirst);
      lowYLast = TMath::Power(10, lowYLast);

      if (lowYFirst <= 0 || lowYLast <= 0) {
         Error(thisfunc, "Cannot set lower Y axis to log scale");
      }

   }

   // this is different than in y, y already has pad coords converted, x not...
   if (logx) {
      if (first <= 0 || last <= 0) {
         Error(thisfunc, "Cannot set X axis to log scale");
      }
   }

   // determine axes options to create log axes if needed
   TString xopt = "";
   if (logx) xopt.Append("G");
   TString upyopt = "";
   if (uplogy) upyopt.Append("G");
   TString lowyopt = "";
   if (lowlogy) lowyopt.Append("G");

   // only actually create them once, reuse otherwise b/c memory
   if (fUpperGXaxis == 0) {
      fUpperGXaxis = new TGaxis(0, 0, 1, 1, 0, 1, 510, "+U"+xopt);
      fUpperGXaxis->Draw();
   }

   if (fUpperGYaxis == 0) {
      fUpperGYaxis = new TGaxis(0, 0, 1, 1, upYFirst, upYLast, 510, "S"+upyopt);
      fUpperGYaxis->Draw();
   }

   if (fLowerGXaxis == 0) {
      fLowerGXaxis = new TGaxis(0, 0, 1, 1, first, last, 510, "+S"+xopt);
      fLowerGXaxis->Draw();
   }

   if (fLowerGYaxis == 0) {
      fLowerGYaxis = new TGaxis(0, 0, 1, 1, lowYFirst, lowYLast, 510, "-S"+lowyopt);
      fLowerGYaxis->Draw();
   }

   // import infos from TAxis
   ImportAxisAttributes(fUpperGXaxis, GetUpperRefXaxis());
   ImportAxisAttributes(fUpperGYaxis, GetUpperRefYaxis());
   ImportAxisAttributes(fLowerGXaxis, GetLowerRefXaxis());
   ImportAxisAttributes(fLowerGYaxis, GetLowerRefYaxis());

   // lower x axis needs to get title from upper x
   fLowerGXaxis->SetTitle(fUpperGXaxis->GetTitle());

   // (re)set all the axes properties to what we want them
   fUpperGXaxis->SetTitle("");

   fUpperGXaxis->SetX1(upLM);
   fUpperGXaxis->SetX2(1-upRM);
   fUpperGXaxis->SetY1(upBM*(1-sf)+sf);
   fUpperGXaxis->SetY2(upBM*(1-sf)+sf);
   fUpperGXaxis->SetWmin(first);
   fUpperGXaxis->SetWmax(last);

   fUpperGYaxis->SetX1(upLM);
   fUpperGYaxis->SetX2(upLM);
   fUpperGYaxis->SetY1(upBM*(1-sf)+sf);
   fUpperGYaxis->SetY2( (1-upTM)*(1-sf)+sf );
   fUpperGYaxis->SetWmin(upYFirst);
   fUpperGYaxis->SetWmax(upYLast);

   fLowerGXaxis->SetX1(lowLM);
   fLowerGXaxis->SetX2(1-lowRM);
   fLowerGXaxis->SetY1(lowBM*sf);
   fLowerGXaxis->SetY2(lowBM*sf);
   fLowerGXaxis->SetWmin(first);
   fLowerGXaxis->SetWmax(last);

   fLowerGYaxis->SetX1(lowLM);
   fLowerGYaxis->SetX2(lowLM);
   fLowerGYaxis->SetY1(lowBM*sf);
   fLowerGYaxis->SetY2((1-lowTM)*sf);
   fLowerGYaxis->SetWmin(lowYFirst);
   fLowerGYaxis->SetWmax(lowYLast);

   fUpperGXaxis->SetNdivisions(fSharedXAxis->GetNdivisions());
   fUpperGYaxis->SetNdivisions(fUpYaxis->GetNdivisions());
   fLowerGXaxis->SetNdivisions(fSharedXAxis->GetNdivisions());
   fLowerGYaxis->SetNdivisions(fLowYaxis->GetNdivisions());

   fUpperGXaxis->SetOption("+U"+xopt);
   fUpperGYaxis->SetOption("S"+upyopt);
   fLowerGXaxis->SetOption("+S"+xopt);
   fLowerGYaxis->SetOption("-S"+lowyopt);

   // normalize the tick sizes. y axis ticks should be consistent
   // even if their length is different
   Double_t ratio = ( (upBM-(1-upTM))*(1-sf) ) / ( (lowBM-(1-lowTM))*sf ) ;
   fUpperGXaxis->SetLabelSize(0.);
   Double_t ticksize = fUpperGYaxis->GetTickSize()*ratio;
   fLowerGYaxis->SetTickSize(ticksize);

   if (fHideLabelMode == TRatioPlot::HideLabelMode::kForceHideUp) {

      fUpperGYaxis->ChangeLabel(1, -1, 0);

   } else if (fHideLabelMode == TRatioPlot::HideLabelMode::kForceHideLow) {

      fLowerGYaxis->ChangeLabel(-1, -1, 0);

   } else {
      if (GetSeparationMargin() < 0.025) {

         if (fHideLabelMode != TRatioPlot::HideLabelMode::kNoHide) {
            if (fHideLabelMode == TRatioPlot::HideLabelMode::kHideUp) {
               fUpperGYaxis->ChangeLabel(1, -1, 0);
            } else if (fHideLabelMode == TRatioPlot::HideLabelMode::kHideLow) {
               fLowerGYaxis->ChangeLabel(-1, -1, 0);
            }
         }

      } else {
         // reset
         if (fHideLabelMode == TRatioPlot::HideLabelMode::kHideUp) {
            fUpperGYaxis->ChangeLabel(0);
         } else if (fHideLabelMode == TRatioPlot::HideLabelMode::kHideLow) {
            fLowerGYaxis->ChangeLabel(0);
         }

      }
   }

   // Create the axes on the other sides of the graphs
   // This is steered by an option on the containing pad or self
   if (axistop || axisright) {

      // only actually create them once, reuse otherwise b/c memory
      if (fUpperGXaxisMirror == 0) {
         fUpperGXaxisMirror = (TGaxis*)fUpperGXaxis->Clone();
         if (axistop) fUpperGXaxisMirror->Draw();
      }

      if (fLowerGXaxisMirror == 0) {
         fLowerGXaxisMirror = (TGaxis*)fLowerGXaxis->Clone();
         if (axistop) fLowerGXaxisMirror->Draw();
      }

      if (fUpperGYaxisMirror == 0) {
         fUpperGYaxisMirror = (TGaxis*)fUpperGYaxis->Clone();
         if (axisright) fUpperGYaxisMirror->Draw();
      }

      if (fLowerGYaxisMirror == 0) {
         fLowerGYaxisMirror = (TGaxis*)fLowerGYaxis->Clone();
         if (axisright) fLowerGYaxisMirror->Draw();
      }

      // import attributes from shared axes
      ImportAxisAttributes(fUpperGXaxisMirror, GetUpperRefXaxis());
      ImportAxisAttributes(fUpperGYaxisMirror, GetUpperRefYaxis());
      ImportAxisAttributes(fLowerGXaxisMirror, GetLowerRefXaxis());
      ImportAxisAttributes(fLowerGYaxisMirror, GetLowerRefYaxis());

      // remove titles
      fUpperGXaxisMirror->SetTitle("");
      fUpperGYaxisMirror->SetTitle("");
      fLowerGXaxisMirror->SetTitle("");
      fLowerGYaxisMirror->SetTitle("");

      // move them about and set required positions
      fUpperGXaxisMirror->SetX1(upLM);
      fUpperGXaxisMirror->SetX2(1-upRM);
      fUpperGXaxisMirror->SetY1((1-upTM)*(1-sf)+sf);
      fUpperGXaxisMirror->SetY2((1-upTM)*(1-sf)+sf);
      fUpperGXaxisMirror->SetWmin(first);
      fUpperGXaxisMirror->SetWmax(last);

      fUpperGYaxisMirror->SetX1(1-upRM);
      fUpperGYaxisMirror->SetX2(1-upRM);
      fUpperGYaxisMirror->SetY1(upBM*(1-sf)+sf);
      fUpperGYaxisMirror->SetY2( (1-upTM)*(1-sf)+sf );
      fUpperGYaxisMirror->SetWmin(upYFirst);
      fUpperGYaxisMirror->SetWmax(upYLast);

      fLowerGXaxisMirror->SetX1(lowLM);
      fLowerGXaxisMirror->SetX2(1-lowRM);
      fLowerGXaxisMirror->SetY1((1-lowTM)*sf);
      fLowerGXaxisMirror->SetY2((1-lowTM)*sf);
      fLowerGXaxisMirror->SetWmin(first);
      fLowerGXaxisMirror->SetWmax(last);

      fLowerGYaxisMirror->SetX1(1-lowRM);
      fLowerGYaxisMirror->SetX2(1-lowRM);
      fLowerGYaxisMirror->SetY1(lowBM*sf);
      fLowerGYaxisMirror->SetY2((1-lowTM)*sf);
      fLowerGYaxisMirror->SetWmin(lowYFirst);
      fLowerGYaxisMirror->SetWmax(lowYLast);

      // also needs normalized tick size
      fLowerGYaxisMirror->SetTickSize(ticksize);

      fUpperGXaxisMirror->SetOption("-S"+xopt);
      fUpperGYaxisMirror->SetOption("+S"+upyopt);
      fLowerGXaxisMirror->SetOption("-S"+xopt);
      fLowerGYaxisMirror->SetOption("+S"+lowyopt);

      fUpperGXaxisMirror->SetNdivisions(fSharedXAxis->GetNdivisions());
      fUpperGYaxisMirror->SetNdivisions(fUpYaxis->GetNdivisions());
      fLowerGXaxisMirror->SetNdivisions(fSharedXAxis->GetNdivisions());
      fLowerGYaxisMirror->SetNdivisions(fLowYaxis->GetNdivisions());

      fUpperGXaxisMirror->SetLabelSize(0.);
      fLowerGXaxisMirror->SetLabelSize(0.);
      fUpperGYaxisMirror->SetLabelSize(0.);
      fLowerGYaxisMirror->SetLabelSize(0.);
   }

   padsav->cd();

}

////////////////////////////////////////////////////////////////////////////////
/// Sets the margins of all the pads to the value specified in class members.
/// This one is called whenever those are changed, e.g. in setters

void TRatioPlot::SetPadMargins()
{
   fUpperPad->SetTopMargin(fUpTopMargin);
   fUpperPad->SetBottomMargin(fUpBottomMargin);
   fUpperPad->SetLeftMargin(fLeftMargin);
   fUpperPad->SetRightMargin(fRightMargin);
   fLowerPad->SetTopMargin(fLowTopMargin);
   fLowerPad->SetBottomMargin(fLowBottomMargin);
   fLowerPad->SetLeftMargin(fLeftMargin);
   fLowerPad->SetRightMargin(fRightMargin);
}

////////////////////////////////////////////////////////////////////////////////
/// Figures out which pad margin has deviated from the stored ones,
/// to figure out what the new nominal is and set the other pad to it
/// subsequently.

Bool_t TRatioPlot::SyncPadMargins()
{

   Bool_t changed = kFALSE;

   if (fUpperPad->GetLeftMargin() != fLeftMargin) {
      fLeftMargin = fUpperPad->GetLeftMargin();
      changed = kTRUE;
   }
   else if (fLowerPad->GetLeftMargin() != fLeftMargin) {
      fLeftMargin = fLowerPad->GetLeftMargin();
      changed = kTRUE;
   }

   if (fUpperPad->GetRightMargin() != fRightMargin) {
      fRightMargin = fUpperPad->GetRightMargin();
      changed = kTRUE;
   }
   else if (fLowerPad->GetRightMargin() != fRightMargin) {
      fRightMargin = fLowerPad->GetRightMargin();
      changed = kTRUE;
   }

   // only reset margins, if any of the margins changed
   if (changed) {
      SetPadMargins();
   }

   Bool_t verticalChanged = kFALSE;

   if (fUpperPad->GetBottomMargin() != fUpBottomMargin) {

      verticalChanged = kTRUE;
      fUpBottomMargin = fUpperPad->GetBottomMargin();

   }

   if (fLowerPad->GetTopMargin() != fLowTopMargin) {

      verticalChanged = kTRUE;
      fLowTopMargin = fLowerPad->GetTopMargin();

   }

   if (fLowerPad->GetBottomMargin() != fLowBottomMargin) {

      fLowBottomMargin = fLowerPad->GetBottomMargin();

   }

   if (fUpperPad->GetTopMargin() != fUpTopMargin) {

      fUpTopMargin = fUpperPad->GetTopMargin();

   }

   // only reset margins, if any of the margins changed
   if (verticalChanged) {
      SetPadMargins();
   }

   return changed || verticalChanged;

}

////////////////////////////////////////////////////////////////////////////////
/// Slot that receives the RangeAxisChanged signal from any of the pads and
/// reacts correspondingly.

void TRatioPlot::RangeAxisChanged()
{
   // check if the ratio plot is already drawn.
   if (!IsDrawn()) {
      // not drawn yet
       return;
   }

   // Only run this concurrently once, in case it's called async
   if (fIsUpdating) {
      return;
   }

   fIsUpdating = kTRUE;

   // find out if logx has changed
   if (fParentPad->GetLogx()) {
      if (!fUpperPad->GetLogx() || !fLowerPad->GetLogx()) {
         fParentPad->SetLogx(kFALSE);
      }
   } else {
      if (fUpperPad->GetLogx() || fLowerPad->GetLogx()) {
         fParentPad->SetLogx(kTRUE);
      }
   }

   // set log to pad
   fUpperPad->SetLogx(fParentPad->GetLogx());
   fLowerPad->SetLogx(fParentPad->GetLogx());

   // get axis ranges for upper and lower
   TAxis *uprefx = GetUpperRefXaxis();
   Double_t upFirst = uprefx->GetBinLowEdge(uprefx->GetFirst());
   Double_t upLast  = uprefx->GetBinUpEdge(uprefx->GetLast());

   TAxis *lowrefx = GetLowerRefXaxis();
   Double_t lowFirst = lowrefx->GetBinLowEdge(lowrefx->GetFirst());
   Double_t lowLast = lowrefx->GetBinUpEdge(lowrefx->GetLast());

   Double_t globFirst = fSharedXAxis->GetBinLowEdge(fSharedXAxis->GetFirst());
   Double_t globLast = fSharedXAxis->GetBinUpEdge(fSharedXAxis->GetLast());

   Bool_t upChanged = kFALSE;
   Bool_t lowChanged = kFALSE;

   // determine which one has changed
   if (upFirst != globFirst || upLast != globLast) {
      fSharedXAxis->SetRangeUser(upFirst, upLast);
      upChanged = kTRUE;
   }
   else if (lowFirst != globFirst || lowLast != globLast) {
      fSharedXAxis->SetRangeUser(lowFirst, lowLast);
      lowChanged = kTRUE;
   }

   if (upChanged || lowChanged) {
      SyncAxesRanges();
      CreateVisualAxes();
      CreateGridline();

      // @TODO: Updating is not working when zooming on the lower plot. Axes update, but upper hist only on resize
      fUpperPad->Modified();
      fLowerPad->Modified();
      fTopPad->Modified();
      fParentPad->Modified();
   }

   // sync the margins in case the user has dragged one of them
   Bool_t marginsChanged = SyncPadMargins();

   if (marginsChanged) {
      fUpperPad->Modified();
      fLowerPad->Modified();
      fTopPad->Modified();
      fParentPad->Modified();
   }

   CreateVisualAxes();
   CreateGridline();
   fIsUpdating = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for the UnZoomed signal that was introduced to TAxis.
/// Unzoom both pads

void TRatioPlot::UnZoomed()
{
   // this is what resets the range
   fSharedXAxis->SetRange(0, 0);
   SyncAxesRanges();

   // Flushing
   fUpperPad->Modified();
   fLowerPad->Modified();
   fTopPad->Modified();
   fParentPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot that handles common resizing of upper and lower pad.

void TRatioPlot::SubPadResized()
{

   if (fIsPadUpdating) {
      return;
   }

   fIsPadUpdating = kTRUE;

   Float_t upylow = fUpperPad->GetYlowNDC();
   Float_t lowylow = fLowerPad->GetYlowNDC();
   Float_t lowh = fLowerPad->GetHNDC();
   Float_t lowyup = lowylow + lowh;

   Bool_t changed = kFALSE;

   if (upylow != fSplitFraction) {
      // up changed
      SetSplitFraction(upylow);
      changed = kTRUE;
   }
   else if (lowyup != fSplitFraction) {
      // low changed
      SetSplitFraction(lowyup);
      changed = kTRUE;
   }

   if (changed) {
      CreateVisualAxes();
   }

   fIsPadUpdating = kFALSE;

}

////////////////////////////////////////////////////////////////////////////////
/// Check if ... is drawn.

Bool_t TRatioPlot::IsDrawn()
{
   TList *siblings = fParentPad->GetListOfPrimitives();
   return siblings->FindObject(this) != 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the fraction of the parent pad, at which the to sub pads should meet

void TRatioPlot::SetSplitFraction(Float_t sf)
{
   if (fParentPad == 0) {
      Warning("SetSplitFraction", "Can only be used after TRatioPlot has been drawn.");
      return;
   }

   fSplitFraction = sf;
   double pm = fInsetWidth;
   double width = fParentPad->GetWNDC();
   double height = fParentPad->GetHNDC();
   double f = height/width;

   fUpperPad->SetPad(pm*f, fSplitFraction, 1.-pm*f, 1.-pm);
   fLowerPad->SetPad(pm*f, pm, 1.-pm*f, fSplitFraction);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the inset on the outer sides of all the pads. It's used to make the outer
/// pad draggable.

void TRatioPlot::SetInsetWidth(Double_t width)
{
   if (fParentPad == 0) {
      Warning("SetInsetWidth", "Can only be used after TRatioPlot has been drawn.");
      return;
   }

   fInsetWidth = width;
   SetSplitFraction(fSplitFraction);

   double pm = fInsetWidth;
   double w = fParentPad->GetWNDC();
   double h = fParentPad->GetHNDC();
   double f = h/w;
   fTopPad->SetPad(pm*f, pm, 1-pm*f, 1-pm);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the confidence levels used to calculate the bands in the fit residual
/// case. Defaults to 1 and 2 sigma.

void TRatioPlot::SetConfidenceLevels(Double_t c1, Double_t c2)
{
   fCl1 = c1;
   fCl2 = c2;
   if (!BuildLowerPlot()) return;
}

////////////////////////////////////////////////////////////////////////////////
/// Set where horizontal, dashed lines are drawn on the lower pad.
/// Can be used to override existing default lines (or disable them).
///
/// \param gridlines Vector of y positions for the dashes lines
///
/// Begin_Macro(source)
/// ../../../tutorials/hist/ratioplot4.C
/// End_Macro

void TRatioPlot::SetGridlines(std::vector<double> gridlines)
{
   fGridlinePositions = gridlines;
}

////////////////////////////////////////////////////////////////////////////////
/// Set where horizontal, dashed lines are drawn on the lower pad.
/// Can be used to override existing default lines (or disable them).
///
/// \param gridlines Double_t array of y positions for the dashed lines
/// \param numGridlines Length of gridlines

void TRatioPlot::SetGridlines(Double_t *gridlines, Int_t numGridlines)
{
   fGridlinePositions.clear();

   for (Int_t i=0;i<numGridlines;++i) {
      fGridlinePositions.push_back(gridlines[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the confidence interval colors.
///
/// \param ci1 Color of the 1 sigma band
/// \param ci2 Color of the 2 sigma band
/// Sets the color of the 1 and 2 sigma bands in the fit residual case.
/// Begin_Macro(source)
/// ../../../tutorials/hist/ratioplot5.C
/// End_Macro

void TRatioPlot::SetConfidenceIntervalColors(Color_t ci1, Color_t ci2)
{
   fCi1Color = ci1;
   fCi2Color = ci2;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal method to import TAxis attributes to a TGaxis. Copied from
/// `TGaxis::ImportAxisAttributes`

void TRatioPlot::ImportAxisAttributes(TGaxis *gaxis, TAxis *axis)
{
   gaxis->SetLineColor(axis->GetAxisColor());
   gaxis->SetTextColor(axis->GetTitleColor());
   gaxis->SetTextFont(axis->GetTitleFont());
   gaxis->SetLabelColor(axis->GetLabelColor());
   gaxis->SetLabelFont(axis->GetLabelFont());
   gaxis->SetLabelSize(axis->GetLabelSize());
   gaxis->SetLabelOffset(axis->GetLabelOffset());
   gaxis->SetTickSize(axis->GetTickLength());
   gaxis->SetTitle(axis->GetTitle());
   gaxis->SetTitleOffset(axis->GetTitleOffset());
   gaxis->SetTitleSize(axis->GetTitleSize());
   gaxis->SetBit(TAxis::kCenterTitle,   axis->TestBit(TAxis::kCenterTitle));
   gaxis->SetBit(TAxis::kCenterLabels,  axis->TestBit(TAxis::kCenterLabels));
   gaxis->SetBit(TAxis::kRotateTitle,   axis->TestBit(TAxis::kRotateTitle));
   gaxis->SetBit(TAxis::kNoExponent,    axis->TestBit(TAxis::kNoExponent));
   gaxis->SetBit(TAxis::kTickPlus,      axis->TestBit(TAxis::kTickPlus));
   gaxis->SetBit(TAxis::kTickMinus,     axis->TestBit(TAxis::kTickMinus));
   gaxis->SetBit(TAxis::kMoreLogLabels, axis->TestBit(TAxis::kMoreLogLabels));
   if (axis->GetDecimals())      gaxis->SetBit(TAxis::kDecimals); //the bit is in TAxis::fAxis2
   gaxis->SetTimeFormat(axis->GetTimeFormat());
}
