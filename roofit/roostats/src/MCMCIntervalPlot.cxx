// @(#)root/roostats:$Id$
// Authors: Kevin Belasco        17/06/2009
// Authors: Kyle Cranmer         17/06/2009
/*************************************************************************
 * Project: RooStats                                                     *
 * Package: RooFit/RooStats                                              *
 *************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class RooStats::MCMCIntervalPlot
    \ingroup Roostats

This class provides simple and straightforward utilities to plot a MCMCInterval
object.  Basic use only requires a few lines once you have an MCMCInterval*:

~~~ {.cpp}
      MCMCIntervalPlot plot(*interval);
      plot.Draw();
~~~

The standard Draw() function will currently draw the confidence interval
range with bars if 1-D and a contour if 2-D.  The MCMC posterior will also be
plotted for the 1-D case.

*/

#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MarkovChain.h"

#include "TLine.h"
#include "TList.h"
#include "TGraph.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "TH2.h"
#include "TH1F.h"
#include "RooArgList.h"
#include "TAxis.h"
#include "RooGlobalFunc.h"

#include <iostream>

// Extra draw commands
//static const char* POSTERIOR_HIST = "posterior_hist";
//static const char* POSTERIOR_KEYS_PDF = "posterior_keys_pdf";
//static const char* POSTERIOR_KEYS_PRODUCT = "posterior_keys_product";
//static const char* HIST_INTERVAL = "hist_interval";
//static const char* KEYS_PDF_INTERVAL = "keys_pdf_interval";
//static const char* TAIL_FRACTION_INTERVAL = "tail_fraction_interval";
//static const char* OPTION_SEP = ":";

ClassImp(RooStats::MCMCIntervalPlot);

using namespace std;
using namespace RooStats;

////////////////////////////////////////////////////////////////////////////////

MCMCIntervalPlot::MCMCIntervalPlot()
{
   fInterval = NULL;
   fParameters = NULL;
   fPosteriorHist = NULL;
   fPosteriorKeysPdf = NULL;
   fPosteriorKeysProduct = NULL;
   fDimension = 0;
   fLineColor = kBlack;
   fShadeColor = kGray;
   fLineWidth = 1;
   //fContourColor = kBlack;
   fShowBurnIn = kTRUE;
   fWalk = NULL;
   fBurnIn = NULL;
   fFirst = NULL;
   fParamGraph = NULL;
   fNLLGraph = NULL;
   fNLLHist = NULL;
   fWeightHist = NULL;
   fPosteriorHistHistCopy = NULL;
   fPosteriorHistTFCopy = NULL;
}

////////////////////////////////////////////////////////////////////////////////

MCMCIntervalPlot::MCMCIntervalPlot(MCMCInterval& interval)
{
   SetMCMCInterval(interval);
   fPosteriorHist = NULL;
   fPosteriorKeysPdf = NULL;
   fPosteriorKeysProduct = NULL;
   fLineColor = kBlack;
   fShadeColor = kGray;
   fLineWidth = 1;
   //fContourColor = kBlack;
   fShowBurnIn = kTRUE;
   fWalk = NULL;
   fBurnIn = NULL;
   fFirst = NULL;
   fParamGraph = NULL;
   fNLLGraph = NULL;
   fNLLHist = NULL;
   fWeightHist = NULL;
   fPosteriorHistHistCopy = NULL;
   fPosteriorHistTFCopy = NULL;
}

////////////////////////////////////////////////////////////////////////////////

MCMCIntervalPlot::~MCMCIntervalPlot()
{
   delete fParameters;
   // kbelasco: why does deleting fPosteriorHist remove the graphics
   // but deleting TGraphs doesn't?
   //delete fPosteriorHist;
   // can we delete fNLLHist and fWeightHist?
   //delete fNLLHist;
   //delete fWeightHist;

   // kbelasco: should we delete fPosteriorKeysPdf and fPosteriorKeysProduct?
   delete fPosteriorKeysPdf;
   delete fPosteriorKeysProduct;

   delete fWalk;
   delete fBurnIn;
   delete fFirst;
   delete fParamGraph;
   delete fNLLGraph;
}

////////////////////////////////////////////////////////////////////////////////

void MCMCIntervalPlot::SetMCMCInterval(MCMCInterval& interval)
{
   fInterval = &interval;
   fDimension = fInterval->GetDimension();
   fParameters = fInterval->GetParameters();
}

////////////////////////////////////////////////////////////////////////////////

void MCMCIntervalPlot::Draw(const Option_t* options)
{
   DrawInterval(options);
}

////////////////////////////////////////////////////////////////////////////////

void MCMCIntervalPlot::DrawPosterior(const Option_t* options)
{
   if (fInterval->GetUseKeys())
      DrawPosteriorKeysPdf(options);
   else
      DrawPosteriorHist(options);
}

////////////////////////////////////////////////////////////////////////////////

void* MCMCIntervalPlot::DrawPosteriorHist(const Option_t* /*options*/,
      const char* title, Bool_t scale)
{
   if (fPosteriorHist == NULL)
      fPosteriorHist = fInterval->GetPosteriorHist();

   if (fPosteriorHist == NULL) {
      coutE(InputArguments) << "MCMCIntervalPlot::DrawPosteriorHist: "
         << "Couldn't get posterior histogram." << endl;
      return NULL;
   }

   // kbelasco: annoying hack because histogram drawing fails when it sees
   // an un-recognized option like POSTERIOR_HIST, etc.
   //const Option_t* myOpt = NULL;

   //TString tmpOpt(options);
   //if (tmpOpt.Contains("same"))
   //   myOpt = "same";

   // scale so highest bin has height 1
   if (scale)
      fPosteriorHist->Scale(1/fPosteriorHist->GetBinContent(
               fPosteriorHist->GetMaximumBin()));

   TString ourTitle(GetTitle());
   if (ourTitle.CompareTo("") == 0) {
      if (title)
         fPosteriorHist->SetTitle(title);
   } else
      fPosteriorHist->SetTitle(GetTitle());

   //fPosteriorHist->Draw(myOpt);

   return (void*)fPosteriorHist;
}

////////////////////////////////////////////////////////////////////////////////

void* MCMCIntervalPlot::DrawPosteriorKeysPdf(const Option_t* options)
{
   if (fPosteriorKeysPdf == NULL)
      fPosteriorKeysPdf = fInterval->GetPosteriorKeysPdf();

   if (fPosteriorKeysPdf == NULL) {
      coutE(InputArguments) << "MCMCIntervalPlot::DrawPosteriorKeysPdf: "
         << "Couldn't get posterior Keys PDF." << endl;
      return NULL;
   }

   TString title(GetTitle());
   Bool_t isEmpty = (title.CompareTo("") == 0);

   if (fDimension == 1) {
      RooRealVar* v = (RooRealVar*)fParameters->first();
      RooPlot* frame = v->frame();
      if (frame == NULL) {
         coutE(InputArguments) << "MCMCIntervalPlot::DrawPosteriorKeysPdf: "
                               << "Invalid parameter" << endl;
         return NULL;
      }
      if (isEmpty)
         frame->SetTitle(Form("Posterior Keys PDF for %s", v->GetName()));
      else
         frame->SetTitle(GetTitle());
      //fPosteriorKeysPdf->plotOn(frame);
      //fPosteriorKeysPdf->plotOn(frame,
      //      RooFit::Normalization(1, RooAbsReal::Raw));
      //frame->Draw(options);
      return (void*)frame;
   } else if (fDimension == 2) {
      RooArgList* axes = fInterval->GetAxes();
      RooRealVar* xVar = (RooRealVar*)axes->at(0);
      RooRealVar* yVar = (RooRealVar*)axes->at(1);
      TH2F* keysHist = (TH2F*)fPosteriorKeysPdf->createHistogram(
            "keysPlot2D", *xVar, RooFit::YVar(*yVar), RooFit::Scaling(kFALSE));
      if (isEmpty)
         keysHist->SetTitle(
               Form("MCMC histogram of posterior Keys PDF for %s, %s",
                  axes->at(0)->GetName(), axes->at(1)->GetName()));
      else
         keysHist->SetTitle(GetTitle());

      keysHist->Draw(options);
      delete axes;
      return NULL;
   }
   return NULL;
}

////////////////////////////////////////////////////////////////////////////////

void MCMCIntervalPlot::DrawInterval(const Option_t* options)
{
   switch (fInterval->GetIntervalType()) {
      case MCMCInterval::kShortest:
         DrawShortestInterval(options);
         break;
      case MCMCInterval::kTailFraction:
         DrawTailFractionInterval(options);
         break;
      default:
         coutE(InputArguments) << "MCMCIntervalPlot::DrawInterval(): " <<
            "Interval type not supported" << endl;
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////

void MCMCIntervalPlot::DrawShortestInterval(const Option_t* options)
{
   if (fInterval->GetUseKeys())
      DrawKeysPdfInterval(options);
   else
      DrawHistInterval(options);
}

////////////////////////////////////////////////////////////////////////////////

void MCMCIntervalPlot::DrawKeysPdfInterval(const Option_t* options)
{
   TString title(GetTitle());
   Bool_t isEmpty = (title.CompareTo("") == 0);

   if (fDimension == 1) {
      // Draw the posterior keys PDF as well so the user can see where the
      // limit bars line up
      // fDimension == 1, so we know we will receive a RooPlot
      RooPlot* frame = (RooPlot*)DrawPosteriorKeysPdf(options);

      //Double_t height = 1;
      //Double_t height = 2.0 * fInterval->GetKeysPdfCutoff();
      Double_t height = fInterval->GetKeysMax();

      RooRealVar* p = (RooRealVar*)fParameters->first();
      Double_t ul = fInterval->UpperLimitByKeys(*p);
      Double_t ll = fInterval->LowerLimitByKeys(*p);

      if (frame != NULL && fPosteriorKeysPdf != NULL) {
         // draw shading in interval
         if (isEmpty)
            frame->SetTitle(NULL);
         else
            frame->SetTitle(GetTitle());
         frame->GetYaxis()->SetTitle(Form("Posterior for parameter %s",
                  p->GetName()));
         fPosteriorKeysPdf->plotOn(frame,
               RooFit::Normalization(1, RooAbsReal::Raw),
               RooFit::Range(ll, ul, kFALSE),
               RooFit::VLines(),
               RooFit::DrawOption("F"),
               RooFit::MoveToBack(),
               RooFit::FillColor(fShadeColor));

         // hack - this is drawn twice now:
         // once by DrawPosteriorKeysPdf (which also configures things and sets
         // the title), and once again here so the shading shows up behind.
         fPosteriorKeysPdf->plotOn(frame,
               RooFit::Normalization(1, RooAbsReal::Raw));
      }
      if (frame) {
   frame->Draw(options);
      }

      TLine* llLine = new TLine(ll, 0, ll, height);
      TLine* ulLine = new TLine(ul, 0, ul, height);
      llLine->SetLineColor(fLineColor);
      ulLine->SetLineColor(fLineColor);
      llLine->SetLineWidth(fLineWidth);
      ulLine->SetLineWidth(fLineWidth);
      llLine->Draw(options);
      ulLine->Draw(options);
   } else if (fDimension == 2) {
      if (fPosteriorKeysPdf == NULL)
         fPosteriorKeysPdf = fInterval->GetPosteriorKeysPdf();

      if (fPosteriorKeysPdf == NULL) {
         coutE(InputArguments) << "MCMCIntervalPlot::DrawKeysPdfInterval: "
            << "Couldn't get posterior Keys PDF." << endl;
         return;
      }

      RooArgList* axes = fInterval->GetAxes();
      RooRealVar* xVar = (RooRealVar*)axes->at(0);
      RooRealVar* yVar = (RooRealVar*)axes->at(1);
      TH2F* contHist = (TH2F*)fPosteriorKeysPdf->createHistogram(
          "keysContour2D", *xVar, RooFit::YVar(*yVar), RooFit::Scaling(kFALSE));
      //if (isEmpty)
      //   contHist->SetTitle(Form("MCMC Keys conf. interval for %s, %s",
      //            axes->at(0)->GetName(), axes->at(1)->GetName()));
      //else
      //   contHist->SetTitle(GetTitle());
      if (!isEmpty)
         contHist->SetTitle(GetTitle());
      else
         contHist->SetTitle(NULL);

      contHist->SetStats(kFALSE);

      TString tmpOpt(options);
      if (!tmpOpt.Contains("CONT2")) tmpOpt.Append("CONT2");

      Double_t cutoff = fInterval->GetKeysPdfCutoff();
      contHist->SetContour(1, &cutoff);
      contHist->SetLineColor(fLineColor);
      contHist->SetLineWidth(fLineWidth);
      contHist->Draw(tmpOpt.Data());
      delete axes;
   } else {
      coutE(InputArguments) << "MCMCIntervalPlot::DrawKeysPdfInterval: "
         << " Sorry: " << fDimension << "-D plots not currently supported" << endl;
   }
}

////////////////////////////////////////////////////////////////////////////////

void MCMCIntervalPlot::DrawHistInterval(const Option_t* options)
{
   TString title(GetTitle());
   Bool_t isEmpty = (title.CompareTo("") == 0);

   if (fDimension == 1) {
      // draw lower and upper limits
      RooRealVar* p = (RooRealVar*)fParameters->first();
      Double_t ul = fInterval->UpperLimitByHist(*p);
      Double_t ll = fInterval->LowerLimitByHist(*p);

      // Draw the posterior histogram as well so the user can see where the
      // limit bars line up
      // fDimension == 1, so we know will get a TH1F*
      TH1F* hist = (TH1F*)DrawPosteriorHist(options, NULL, false);
      if (hist == NULL) return;
      if (isEmpty)
         hist->SetTitle(NULL);
      else
         hist->SetTitle(GetTitle());
      hist->GetYaxis()->SetTitle(Form("Posterior for parameter %s",
               p->GetName()));
      hist->SetStats(kFALSE);
      TH1F* copy = (TH1F*)hist->Clone(Form("%s_copy", hist->GetTitle()));
      Double_t histCutoff = fInterval->GetHistCutoff();

      Int_t i;
      Int_t nBins = copy->GetNbinsX();
      Double_t height;
      for (i = 1; i <= nBins; i++) {
         // remove bins with height < cutoff
         height = copy->GetBinContent(i);
         if (height < histCutoff) {
            copy->SetBinContent(i, 0);
            copy->SetBinError(i, 0);
         }
      }

      hist->Scale(1/hist->GetBinContent(hist->GetMaximumBin()));
      copy->Scale(1/copy->GetBinContent(hist->GetMaximumBin()));

      copy->SetFillStyle(1001);
      copy->SetFillColor(fShadeColor);
      hist->Draw(options);
      // to show the interval area
      copy->Draw("HIST SAME");

      fPosteriorHistHistCopy = copy;

      TLine* llLine = new TLine(ll, 0, ll, 1);
      TLine* ulLine = new TLine(ul, 0, ul, 1);
      llLine->SetLineColor(fLineColor);
      ulLine->SetLineColor(fLineColor);
      llLine->SetLineWidth(fLineWidth);
      ulLine->SetLineWidth(fLineWidth);
      llLine->Draw(options);
      ulLine->Draw(options);

   } else if (fDimension == 2) {
      if (fPosteriorHist == NULL)
         fPosteriorHist = fInterval->GetPosteriorHist();

      if (fPosteriorHist == NULL) {
         coutE(InputArguments) << "MCMCIntervalPlot::DrawHistInterval: "
            << "Couldn't get posterior histogram." << endl;
         return;
      }

      RooArgList* axes = fInterval->GetAxes();
      //if (isEmpty)
      //   fPosteriorHist->SetTitle(
      //         Form("MCMC histogram conf. interval for %s, %s",
      //            axes->at(0)->GetName(), axes->at(1)->GetName()));
      //else
      //   fPosteriorHist->SetTitle(GetTitle());
      if (!isEmpty)
         fPosteriorHist->SetTitle(GetTitle());
      else
         fPosteriorHist->SetTitle(NULL);
      delete axes;

      fPosteriorHist->SetStats(kFALSE);

      TString tmpOpt(options);
      if (!tmpOpt.Contains("CONT2")) tmpOpt.Append("CONT2");

      Double_t cutoff = fInterval->GetHistCutoff();
      fPosteriorHist->SetContour(1, &cutoff);
      fPosteriorHist->SetLineColor(fLineColor);
      fPosteriorHist->SetLineWidth(fLineWidth);
      fPosteriorHist->Draw(tmpOpt.Data());
   } else {
      coutE(InputArguments) << "MCMCIntervalPlot::DrawHistInterval: "
         << " Sorry: " << fDimension << "-D plots not currently supported" << endl;
   }
}

////////////////////////////////////////////////////////////////////////////////

void MCMCIntervalPlot::DrawTailFractionInterval(const Option_t* options)
{
   TString title(GetTitle());
   Bool_t isEmpty = (title.CompareTo("") == 0);

   if (fDimension == 1) {
      // Draw the posterior histogram as well so the user can see where the
      // limit bars line up
      RooRealVar* p = (RooRealVar*)fParameters->first();
      Double_t ul = fInterval->UpperLimitTailFraction(*p);
      Double_t ll = fInterval->LowerLimitTailFraction(*p);

      TH1F* hist = (TH1F*)DrawPosteriorHist(options, NULL, false);
      if (hist == NULL) return;
      if (isEmpty)
         hist->SetTitle(NULL);
      else
         hist->SetTitle(GetTitle());
      hist->GetYaxis()->SetTitle(Form("Posterior for parameter %s",
               p->GetName()));
      hist->SetStats(kFALSE);
      TH1F* copy = (TH1F*)hist->Clone(Form("%s_copy", hist->GetTitle()));

      Int_t i;
      Int_t nBins = copy->GetNbinsX();
      Double_t center;
      for (i = 1; i <= nBins; i++) {
         // remove bins outside interval
         center = copy->GetBinCenter(i);
         if (center < ll || center > ul) {
            copy->SetBinContent(i, 0);
            copy->SetBinError(i, 0);
         }
      }

      hist->Scale(1/hist->GetBinContent(hist->GetMaximumBin()));
      copy->Scale(1/copy->GetBinContent(hist->GetMaximumBin()));

      copy->SetFillStyle(1001);
      copy->SetFillColor(fShadeColor);
      hist->Draw(options);
      copy->Draw("hist same");

      // draw lower and upper limits
      TLine* llLine = new TLine(ll, 0, ll, 1);
      TLine* ulLine = new TLine(ul, 0, ul, 1);
      llLine->SetLineColor(fLineColor);
      ulLine->SetLineColor(fLineColor);
      llLine->SetLineWidth(fLineWidth);
      ulLine->SetLineWidth(fLineWidth);
      llLine->Draw(options);
      ulLine->Draw(options);
   } else {
      coutE(InputArguments) << "MCMCIntervalPlot::DrawTailFractionInterval: "
         << " Sorry: " << fDimension << "-D plots not currently supported"
         << endl;
   }
}

////////////////////////////////////////////////////////////////////////////////

void* MCMCIntervalPlot::DrawPosteriorKeysProduct(const Option_t* options)
{
   if (fPosteriorKeysProduct == NULL)
      fPosteriorKeysProduct = fInterval->GetPosteriorKeysProduct();

   if (fPosteriorKeysProduct == NULL) {
      coutE(InputArguments) << "MCMCIntervalPlot::DrawPosteriorKeysProduct: "
         << "Couldn't get posterior Keys product." << endl;
      return NULL;
   }

   RooArgList* axes = fInterval->GetAxes();

   TString title(GetTitle());
   Bool_t isEmpty = (title.CompareTo("") == 0);

   if (fDimension == 1) {
      RooPlot* frame = ((RooRealVar*)fParameters->first())->frame();
      if (!frame) return NULL;
      if (isEmpty)
         frame->SetTitle(Form("Posterior Keys PDF * Heaviside product for %s",
                  axes->at(0)->GetName()));
      else
         frame->SetTitle(GetTitle());
      //fPosteriorKeysProduct->plotOn(frame);
      fPosteriorKeysProduct->plotOn(frame,
            RooFit::Normalization(1, RooAbsReal::Raw));
      frame->Draw(options);
      return (void*)frame;
   } else if (fDimension == 2) {
      RooRealVar* xVar = (RooRealVar*)axes->at(0);
      RooRealVar* yVar = (RooRealVar*)axes->at(1);
      TH2F* productHist = (TH2F*)fPosteriorKeysProduct->createHistogram(
            "prodPlot2D", *xVar, RooFit::YVar(*yVar), RooFit::Scaling(kFALSE));
      if (isEmpty)
         productHist->SetTitle(
               Form("MCMC Posterior Keys Product Hist. for %s, %s",
                  axes->at(0)->GetName(), axes->at(1)->GetName()));
      else
         productHist->SetTitle(GetTitle());
      productHist->Draw(options);
      return NULL;
   }
   delete axes;
   return NULL;
}

////////////////////////////////////////////////////////////////////////////////

void MCMCIntervalPlot::DrawChainScatter(RooRealVar& xVar, RooRealVar& yVar)
{
   const MarkovChain* markovChain = fInterval->GetChain();

   Int_t size = markovChain->Size();
   Int_t burnInSteps;
   if (fShowBurnIn)
      burnInSteps = fInterval->GetNumBurnInSteps();
   else
      burnInSteps = 0;

   Double_t* x = new Double_t[size - burnInSteps];
   Double_t* y = new Double_t[size - burnInSteps];
   Double_t* burnInX = NULL;
   Double_t* burnInY = NULL;
   if (burnInSteps > 0) {
      burnInX = new Double_t[burnInSteps];
      burnInY = new Double_t[burnInSteps];
   }
   Double_t firstX;
   Double_t firstY;

   for (Int_t i = burnInSteps; i < size; i++) {
      x[i - burnInSteps] = markovChain->Get(i)->getRealValue(xVar.GetName());
      y[i - burnInSteps] = markovChain->Get(i)->getRealValue(yVar.GetName());
   }

   for (Int_t i = 0; i < burnInSteps; i++) {
      burnInX[i] = markovChain->Get(i)->getRealValue(xVar.GetName());
      burnInY[i] = markovChain->Get(i)->getRealValue(yVar.GetName());
   }

   firstX = markovChain->Get(0)->getRealValue(xVar.GetName());
   firstY = markovChain->Get(0)->getRealValue(yVar.GetName());

   TString title(GetTitle());
   Bool_t isEmpty = (title.CompareTo("") == 0);

   TGraph* walk = new TGraph(size - burnInSteps, x, y);
   if (isEmpty)
      walk->SetTitle(Form("2-D Scatter Plot of Markov chain for %s, %s",
               xVar.GetName(), yVar.GetName()));
   else
      walk->SetTitle(GetTitle());
   // kbelasco: figure out how to set TGraph variable ranges
   walk->GetXaxis()->Set(xVar.numBins(), xVar.getMin(), xVar.getMax());
   walk->GetXaxis()->SetTitle(xVar.GetName());
   walk->GetYaxis()->Set(yVar.numBins(), yVar.getMin(), yVar.getMax());
   walk->GetYaxis()->SetTitle(yVar.GetName());
   walk->SetLineColor(kGray);
   walk->SetMarkerStyle(6);
   walk->SetMarkerColor(kViolet);
   walk->Draw("A,L,P,same");

   TGraph* burnIn = NULL;
   if (burnInX != NULL && burnInY != NULL) {
      burnIn = new TGraph(burnInSteps - 1, burnInX, burnInY);
      burnIn->SetLineColor(kPink);
      burnIn->SetMarkerStyle(6);
      burnIn->SetMarkerColor(kPink);
      burnIn->Draw("L,P,same");
   }

   TGraph* first = new TGraph(1, &firstX, &firstY);
   first->SetLineColor(kGreen);
   first->SetMarkerStyle(3);
   first->SetMarkerSize(2);
   first->SetMarkerColor(kGreen);
   first->Draw("L,P,same");

   //walkCanvas->Update();
   delete [] x;
   delete [] y;
   if (burnInX != NULL) delete [] burnInX;
   if (burnInY != NULL) delete [] burnInY;
   //delete walk;
   //delete burnIn;
   //delete first;
}

////////////////////////////////////////////////////////////////////////////////

void MCMCIntervalPlot::DrawParameterVsTime(RooRealVar& param)
{
   const MarkovChain* markovChain = fInterval->GetChain();
   Int_t size = markovChain->Size();
   Int_t numEntries = 2 * size;
   Double_t* value = new Double_t[numEntries];
   Double_t* time = new Double_t[numEntries];
   Double_t val;
   Int_t weight;
   Int_t t = 0;
   for (Int_t i = 0; i < size; i++) {
      val = markovChain->Get(i)->getRealValue(param.GetName());
      weight = (Int_t)markovChain->Weight();
      value[2*i] = val;
      value[2*i + 1] = val;
      time[2*i] = t;
      t += weight;
      time[2*i + 1] = t;
   }

   TString title(GetTitle());
   Bool_t isEmpty = (title.CompareTo("") == 0);

   TGraph* paramGraph = new TGraph(numEntries, time, value);
   if (isEmpty)
      paramGraph->SetTitle(Form("%s vs. time in Markov chain",param.GetName()));
   else
      paramGraph->SetTitle(GetTitle());
   paramGraph->GetXaxis()->SetTitle("Time (discrete steps)");
   paramGraph->GetYaxis()->SetTitle(param.GetName());
   //paramGraph->SetLineColor(fLineColor);
   paramGraph->Draw("A,L,same");
   delete [] value;
   delete [] time;
   //gPad->Update();
}

////////////////////////////////////////////////////////////////////////////////

void MCMCIntervalPlot::DrawNLLVsTime()
{
   const MarkovChain* markovChain = fInterval->GetChain();
   Int_t size = markovChain->Size();
   Int_t numEntries = 2 * size;
   Double_t* nllValue = new Double_t[numEntries];
   Double_t* time = new Double_t[numEntries];
   Double_t nll;
   Int_t weight;
   Int_t t = 0;
   for (Int_t i = 0; i < size; i++) {
      nll = markovChain->NLL(i);
      weight = (Int_t)markovChain->Weight();
      nllValue[2*i] = nll;
      nllValue[2*i + 1] = nll;
      time[2*i] = t;
      t += weight;
      time[2*i + 1] = t;
   }

   TString title(GetTitle());
   Bool_t isEmpty = (title.CompareTo("") == 0);

   TGraph* nllGraph = new TGraph(numEntries, time, nllValue);
   if (isEmpty)
      nllGraph->SetTitle("NLL value vs. time in Markov chain");
   else
      nllGraph->SetTitle(GetTitle());
   nllGraph->GetXaxis()->SetTitle("Time (discrete steps)");
   nllGraph->GetYaxis()->SetTitle("NLL (-log(likelihood))");
   //nllGraph->SetLineColor(fLineColor);
   nllGraph->Draw("A,L,same");
   //gPad->Update();
   delete [] nllValue;
   delete [] time;
}

////////////////////////////////////////////////////////////////////////////////

void MCMCIntervalPlot::DrawNLLHist(const Option_t* options)
{
   if (fNLLHist == NULL) {
      const MarkovChain* markovChain = fInterval->GetChain();
      // find the max NLL value
      Double_t maxNLL = 0;
      Int_t size = markovChain->Size();
      for (Int_t i = 0; i < size; i++)
         if (markovChain->NLL(i) > maxNLL)
            maxNLL = markovChain->NLL(i);
      RooRealVar* nllVar = fInterval->GetNLLVar();
      fNLLHist = new TH1F("mcmc_nll_hist", "MCMC NLL Histogram",
            nllVar->getBins(), 0, maxNLL);
      TString title(GetTitle());
      Bool_t isEmpty = (title.CompareTo("") == 0);
      if (!isEmpty)
         fNLLHist->SetTitle(GetTitle());
      fNLLHist->GetXaxis()->SetTitle("-log(likelihood)");
      for (Int_t i = 0; i < size; i++)
         fNLLHist->Fill(markovChain->NLL(i), markovChain->Weight());
   }
   fNLLHist->Draw(options);
}

////////////////////////////////////////////////////////////////////////////////

void MCMCIntervalPlot::DrawWeightHist(const Option_t* options)
{
   if (fWeightHist == NULL) {
      const MarkovChain* markovChain = fInterval->GetChain();
      // find the max weight value
      Double_t maxWeight = 0;
      Int_t size = markovChain->Size();
      for (Int_t i = 0; i < size; i++)
         if (markovChain->Weight(i) > maxWeight)
            maxWeight = markovChain->Weight(i);
      fWeightHist = new TH1F("mcmc_weight_hist", "MCMC Weight Histogram",
            (Int_t)(maxWeight + 1), 0, maxWeight * 1.02);
      for (Int_t i = 0; i < size; i++)
         fWeightHist->Fill(markovChain->Weight(i));
   }
   fWeightHist->Draw(options);
}

/*
/////////////////////////////////////////////////////////////////////
  // 3-d plot of the parameter points
  dataCanvas->cd(2);
  // also plot the points in the markov chain
  RooDataSet* markovChainData = ((MCMCInterval*)mcmcint)->GetChainAsDataSet();

  TTree& chain =  ((RooTreeDataStore*) markovChainData->store())->tree();
  chain.SetMarkerStyle(6);
  chain.SetMarkerColor(kRed);
  chain.Draw("s:ratioSigEff:ratioBkgEff","","box"); // 3-d box proportional to posterior

  // the points used in the profile construction
  TTree& parameterScan =  ((RooTreeDataStore*) fc.GetPointsToScan()->store())->tree();
  parameterScan.SetMarkerStyle(24);
  parameterScan.Draw("s:ratioSigEff:ratioBkgEff","","same");

  chain.SetMarkerStyle(6);
  chain.SetMarkerColor(kRed);
  //chain.Draw("s:ratioSigEff:ratioBkgEff", "_MarkovChain_local_nll","box");
  //chain.Draw("_MarkovChain_local_nll");
////////////////////////////////////////////////////////////////////
*/
