// @(#)root/roostats:$Id$
// Authors: Sven Kreiss    June 2010
// Authors: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class RooStats::SamplingDistPlot
    \ingroup Roostats

This class provides simple and straightforward utilities to plot SamplingDistribution
objects.
*/

#include "RooStats/SamplingDistPlot.h"

#include "RooStats/SamplingDistribution.h"

#include "RooRealVar.h"
#include "TStyle.h"
#include "TLine.h"
#include "TFile.h"
#include "TVirtualPad.h"  // for gPad

#include <algorithm>
#include <iostream>


#include "RooMsgService.h"

#include <limits>
#define NaN std::numeric_limits<float>::quiet_NaN()
#include "TMath.h"
#define IsNaN(a) TMath::IsNaN(a)

ClassImp(RooStats::SamplingDistPlot);

using namespace RooStats;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
/// SamplingDistPlot default constructor with bin size

SamplingDistPlot::SamplingDistPlot(Int_t nbins) :
   fHist(0),
   fLegend(NULL),
   fItems(),
   fOtherItems(),
   fRooPlot(NULL),
   fLogXaxis(kFALSE),
   fLogYaxis(kFALSE),
   fXMin(NaN), fXMax(NaN), fYMin(NaN), fYMax(NaN),
   fApplyStyle(kTRUE),
   fFillStyle(3004)
{
  fIterator = fItems.MakeIterator();
  fIsWeighted = kFALSE;
  fBins = nbins;
  fMarkerType = 20;
  fColor = 1;
}

////////////////////////////////////////////////////////////////////////////////
/// SamplingDistPlot constructor with bin size, min and max

SamplingDistPlot::SamplingDistPlot(Int_t nbins, Double_t min, Double_t max) :
   fHist(0),
   fLegend(NULL),
   fItems(),
   fOtherItems(),
   fRooPlot(NULL),
   fLogXaxis(kFALSE),
   fLogYaxis(kFALSE),
   fXMin(NaN), fXMax(NaN), fYMin(NaN), fYMax(NaN),
   fApplyStyle(kTRUE),
   fFillStyle(3004)
{
  fIterator = fItems.MakeIterator();
  fIsWeighted = kFALSE;
  fBins = nbins;
  fMarkerType = 20;
  fColor = 1;

  SetXRange( min, max );
}

////////////////////////////////////////////////////////////////////////////////
/// destructors - delete objects contained in the list

SamplingDistPlot::~SamplingDistPlot() {
   fItems.Delete();
   fOtherItems.Delete();
   if (fRooPlot) delete fRooPlot;
}


////////////////////////////////////////////////////////////////////////////////
/// adds sampling distribution (and normalizes if "NORMALIZE" is given as an option)

Double_t SamplingDistPlot::AddSamplingDistribution(const SamplingDistribution *samplingDist, Option_t *drawOptions) {
   fSamplingDistr = samplingDist->GetSamplingDistribution();
   if( fSamplingDistr.empty() ) {
      coutW(Plotting) << "Empty sampling distribution given to plot. Skipping." << endl;
      return 0.0;
   }
   SetSampleWeights(samplingDist);

   TString options(drawOptions);
   options.ToUpper();

   Double_t xmin(TMath::Infinity()), xmax(-TMath::Infinity());
   // remove cases where xmin and xmax are +/- inf
   for( unsigned int i=0; i < fSamplingDistr.size(); i++ ) {
      if( fSamplingDistr[i] < xmin  &&  fSamplingDistr[i] != -TMath::Infinity() ) {
         xmin = fSamplingDistr[i];
      }
      if( fSamplingDistr[i] > xmax  &&  fSamplingDistr[i] != TMath::Infinity() ) {
         xmax = fSamplingDistr[i];
      }
   }
   if( xmin >= xmax ) {
      coutW(Plotting) << "Could not determine xmin and xmax of sampling distribution that was given to plot." << endl;
      xmin = -1.0;
      xmax = 1.0;
   }


   // add 1.5 bins left and right
   assert(fBins > 1);
   double binWidth = (xmax-xmin)/(fBins);
   Double_t xlow = xmin - 1.5*binWidth;
   Double_t xup  = xmax + 1.5*binWidth;
   if( !IsNaN(fXMin) ) xlow = fXMin;
   if( !IsNaN(fXMax) ) xup = fXMax;

   fHist = new TH1F(samplingDist->GetName(), samplingDist->GetTitle(), fBins, xlow, xup);
   fHist->SetDirectory(0);  // make the object managed by this class

   if( fVarName.Length() == 0 ) fVarName = samplingDist->GetVarName();
   fHist->GetXaxis()->SetTitle(fVarName.Data());


   std::vector<Double_t>::iterator valuesIt = fSamplingDistr.begin();
   for (int w_idx = 0; valuesIt != fSamplingDistr.end(); ++valuesIt, ++w_idx) {
      if (fIsWeighted) fHist->Fill(*valuesIt, fSampleWeights[w_idx]);
      else fHist->Fill(*valuesIt);
   }

   // NORMALIZATION
   fHist->Sumw2();
   double weightSum = 1.0;
   if(options.Contains("NORMALIZE")) {
      weightSum = fHist->Integral("width");
      fHist->Scale(1./weightSum);

      options.ReplaceAll("NORMALIZE", "");
      options.Strip();
   }


   //some basic aesthetics
   fHist->SetMarkerStyle(fMarkerType);
   fHist->SetMarkerColor(fColor);
   fHist->SetLineColor(fColor);

   fMarkerType++;
   fColor++;

   fHist->SetStats(kFALSE);

   addObject(fHist, options.Data());

   TString title = samplingDist->GetTitle();
   if(fLegend  &&  title.Length() > 0) fLegend->AddEntry(fHist, title, "L");

   return 1./weightSum;
}

////////////////////////////////////////////////////////////////////////////////

Double_t SamplingDistPlot::AddSamplingDistributionShaded(const SamplingDistribution *samplingDist, Double_t minShaded, Double_t maxShaded, Option_t *drawOptions) {
   if( samplingDist->GetSamplingDistribution().empty() ) {
      coutW(Plotting) << "Empty sampling distribution given to plot. Skipping." << endl;
      return 0.0;
   }
   Double_t scaleFactor = AddSamplingDistribution(samplingDist, drawOptions);

   TH1F *shaded = (TH1F*)fHist->Clone((string(samplingDist->GetName())+string("_shaded")).c_str());
   shaded->SetDirectory(0);
   shaded->SetFillStyle(fFillStyle++);
   shaded->SetLineWidth(1);

   for (int i=0; i<shaded->GetNbinsX(); ++i) {
      if (shaded->GetBinCenter(i) < minShaded || shaded->GetBinCenter(i) > maxShaded){
         shaded->SetBinContent(i,0);
      }
   }

   TString options(drawOptions);
   options.ToUpper();
   if(options.Contains("NORMALIZE")) {
      options.ReplaceAll("NORMALIZE", "");
      options.Strip();
   }
   addObject(shaded, options.Data());

   return scaleFactor;
}

////////////////////////////////////////////////////////////////////////////////

void SamplingDistPlot::AddLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2, const char* title) {
   TLine *line = new TLine(x1, y1, x2, y2);
   line->SetLineWidth(3);
   line->SetLineColor(kBlack);

   if(fLegend  &&  title) fLegend->AddEntry(line, title, "L");

   addOtherObject(line, ""); // no options
}

////////////////////////////////////////////////////////////////////////////////
/// add an histogram (it will be cloned);

void SamplingDistPlot::AddTH1(TH1* h, Option_t *drawOptions) {
   if(fLegend  &&  h->GetTitle()) fLegend->AddEntry(h, h->GetTitle(), "L");
   TH1 * hcopy = (TH1*) h->Clone();
   hcopy->SetDirectory(0);
   addObject(hcopy, drawOptions);
}
void SamplingDistPlot::AddTF1(TF1* f, const char* title, Option_t *drawOptions) {
   if(fLegend  &&  title) fLegend->AddEntry(f, title, "L");
   addOtherObject(f->Clone(), drawOptions);
}

////////////////////////////////////////////////////////////////////////////////
///Determine if the sampling distribution has weights and store them

void SamplingDistPlot::SetSampleWeights(const SamplingDistribution* samplingDist)
{
  fIsWeighted = kFALSE;

  if(samplingDist->GetSampleWeights().size() != 0){
    fIsWeighted = kTRUE;
    fSampleWeights = samplingDist->GetSampleWeights();
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a generic object to this plot. The specified options will be
/// used to Draw() this object later. The caller transfers ownership
/// of the object with this call, and the object will be deleted
/// when its containing plot object is destroyed.

void SamplingDistPlot::addObject(TObject *obj, Option_t *drawOptions)
{

  if(0 == obj) {
    std::cerr << fName << "::addObject: called with a null pointer" << std::endl;
    return;
  }

  fItems.Add(obj,drawOptions);

  return;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a generic object to this plot. The specified options will be
/// used to Draw() this object later. The caller transfers ownership
/// of the object with this call, and the object will be deleted
/// when its containing plot object is destroyed.

void SamplingDistPlot::addOtherObject(TObject *obj, Option_t *drawOptions)
{
  if(0 == obj) {
     oocoutE(this,InputArguments) << fName << "::addOtherObject: called with a null pointer" << std::endl;
     return;
  }

  fOtherItems.Add(obj,drawOptions);

  return;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this plot and all of the elements it contains. The specified options
/// only apply to the drawing of our frame. The options specified in our add...()
/// methods will be used to draw each object we contain.

void SamplingDistPlot::Draw(Option_t * /*options */) {
   ApplyDefaultStyle();

   Double_t theMin(0.), theMax(0.), theYMin(NaN), theYMax(0.);
   GetAbsoluteInterval(theMin, theMax, theYMax);
   if( !IsNaN(fXMin) ) theMin = fXMin;
   if( !IsNaN(fXMax) ) theMax = fXMax;
   if( !IsNaN(fYMin) ) theYMin = fYMin;
   if( !IsNaN(fYMax) ) theYMax = fYMax;

   RooRealVar xaxis("xaxis", fVarName.Data(), theMin, theMax);

   //L.M. by drawing many times we create a memory leak ???
   if (fRooPlot) delete fRooPlot;

   bool dirStatus = RooPlot::addDirectoryStatus();
   // make the RooPlot managed by this class
   RooPlot::setAddDirectoryStatus(false);
   fRooPlot = xaxis.frame();
   RooPlot::setAddDirectoryStatus(dirStatus);
   if (!fRooPlot) {
     oocoutE(this,InputArguments) << "invalid variable to plot" << std::endl;
     return;
   }
   fRooPlot->SetTitle("");
   if( !IsNaN(theYMax) ) {
      //coutI(InputArguments) << "Setting maximum to " << theYMax << endl;
      fRooPlot->SetMaximum(theYMax);
   }
   if( !IsNaN(theYMin) ) {
      //coutI(InputArguments) << "Setting minimum to " << theYMin << endl;
      fRooPlot->SetMinimum(theYMin);
   }

   fIterator->Reset();
   TH1F *obj = 0;
   while ((obj = (TH1F*) fIterator->Next())) {
      //obj->Draw(fIterator->GetOption());
      // add cloned objects to avoid mem leaks
      TH1 * cloneObj = (TH1*)obj->Clone();
      if( !IsNaN(theYMax) ) {
         //coutI(InputArguments) << "Setting maximum of TH1 to " << theYMax << endl;
         cloneObj->SetMaximum(theYMax);
      }
      if( !IsNaN(theYMin) ) {
         //coutI(InputArguments) << "Setting minimum of TH1 to " << theYMin << endl;
         cloneObj->SetMinimum(theYMin);
      }
      cloneObj->SetDirectory(0);  // transfer ownership of the object
      fRooPlot->addTH1(cloneObj, fIterator->GetOption());
   }

   TIterator *otherIt = fOtherItems.MakeIterator();
   TObject *otherObj = NULL;
   while ((otherObj = otherIt->Next())) {
      TObject * cloneObj = otherObj->Clone();
      fRooPlot->addObject(cloneObj, otherIt->GetOption());
   }
   delete otherIt;


   if(fLegend) fRooPlot->addObject(fLegend);

   if(bool(gStyle->GetOptLogx()) != fLogXaxis) {
      if(!fApplyStyle) coutW(Plotting) << "gStyle will be changed to adjust SetOptLogx(...)" << endl;
      gStyle->SetOptLogx(fLogXaxis);
   }
   if(bool(gStyle->GetOptLogy()) != fLogYaxis) {
      if(!fApplyStyle) coutW(Plotting) << "gStyle will be changed to adjust SetOptLogy(...)" << endl;
      gStyle->SetOptLogy(fLogYaxis);
   }
   fRooPlot->Draw();

   // apply this since gStyle does not work for RooPlot
   if (gPad) {
      gPad->SetLogx(fLogXaxis);
      gPad->SetLogy(fLogYaxis);
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////

void SamplingDistPlot::ApplyDefaultStyle(void) {
   if(fApplyStyle) {
      // use plain black on white colors
      Int_t icol = 0;
      gStyle->SetFrameBorderMode( icol );
      gStyle->SetCanvasBorderMode( icol );
      gStyle->SetPadBorderMode( icol );
      gStyle->SetPadColor( icol );
      gStyle->SetCanvasColor( icol );
      gStyle->SetStatColor( icol );
      gStyle->SetFrameFillStyle( 0 );

      // set the paper & margin sizes
      gStyle->SetPaperSize( 20, 26 );

      if(fLegend) {
         fLegend->SetFillColor(0);
         fLegend->SetBorderSize(1);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

void SamplingDistPlot::GetAbsoluteInterval(Double_t &theMin, Double_t &theMax, Double_t &theYMax) const
{
   Double_t tmpmin = TMath::Infinity();
   Double_t tmpmax = -TMath::Infinity();
   Double_t tmpYmax = -TMath::Infinity();

  fIterator->Reset();
  TH1F *obj = 0;
  while((obj = (TH1F*)fIterator->Next())) {
    if(obj->GetXaxis()->GetXmin() < tmpmin) tmpmin = obj->GetXaxis()->GetXmin();
    if(obj->GetXaxis()->GetXmax() > tmpmax) tmpmax = obj->GetXaxis()->GetXmax();
    if(obj->GetMaximum() > tmpYmax) tmpYmax = obj->GetMaximum() + 0.1*obj->GetMaximum();
  }

  theMin = tmpmin;
  theMax = tmpmax;
  theYMax = tmpYmax;

  return;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets line color for given sampling distribution and
/// fill color for the associated shaded TH1F.

void SamplingDistPlot::SetLineColor(Color_t color, const SamplingDistribution *samplDist) {
   if (samplDist == 0) {
      fHist->SetLineColor(color);

      fIterator->Reset();
      TH1F *obj = 0;

      TString shadedName(fHist->GetName());
      shadedName += "_shaded";

      while ((obj = (TH1F*) fIterator->Next())) {
         if (!strcmp(obj->GetName(), shadedName.Data())) {
            obj->SetLineColor(color);
            obj->SetFillColor(color);
            //break;
         }
      }
   } else {
      fIterator->Reset();
      TH1F *obj = 0;

      TString shadedName(samplDist->GetName());
      shadedName += "_shaded";

      while ((obj = (TH1F*) fIterator->Next())) {
         if (!strcmp(obj->GetName(), samplDist->GetName())) {
            obj->SetLineColor(color);
            //break;
         }
         if (!strcmp(obj->GetName(), shadedName.Data())) {
            obj->SetLineColor(color);
            obj->SetFillColor(color);
            //break;
         }
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////

void SamplingDistPlot::SetLineWidth(Width_t lwidth, const SamplingDistribution *samplDist)
{
  if(samplDist == 0){
    fHist->SetLineWidth(lwidth);
  }
  else{
    fIterator->Reset();
    TH1F *obj = 0;
    while((obj = (TH1F*)fIterator->Next())) {
      if(!strcmp(obj->GetName(),samplDist->GetName())){
   obj->SetLineWidth(lwidth);
   break;
      }
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////

void SamplingDistPlot::SetLineStyle(Style_t style, const SamplingDistribution *samplDist)
{
  if(samplDist == 0){
    fHist->SetLineStyle(style);
  }
  else{
    fIterator->Reset();
    TH1F *obj = 0;
    while((obj = (TH1F*)fIterator->Next())) {
      if(!strcmp(obj->GetName(),samplDist->GetName())){
   obj->SetLineStyle(style);
   break;
      }
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////

void SamplingDistPlot::SetMarkerStyle(Style_t style, const SamplingDistribution *samplDist)
{
  if(samplDist == 0){
    fHist->SetMarkerStyle(style);
  }
  else{
    fIterator->Reset();
    TH1F *obj = 0;
    while((obj = (TH1F*)fIterator->Next())) {
      if(!strcmp(obj->GetName(),samplDist->GetName())){
   obj->SetMarkerStyle(style);
   break;
      }
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////

void SamplingDistPlot::SetMarkerColor(Color_t color, const SamplingDistribution *samplDist)
{
  if(samplDist == 0){
    fHist->SetMarkerColor(color);
  }
  else{
    fIterator->Reset();
    TH1F *obj = 0;
    while((obj = (TH1F*)fIterator->Next())) {
      if(!strcmp(obj->GetName(),samplDist->GetName())){
   obj->SetMarkerColor(color);
   break;
      }
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////

void SamplingDistPlot::SetMarkerSize(Size_t size, const SamplingDistribution *samplDist)
{
  if(samplDist == 0){
    fHist->SetMarkerSize(size);
  }
  else{
    fIterator->Reset();
    TH1F *obj = 0;
    while((obj = (TH1F*)fIterator->Next())) {
      if(!strcmp(obj->GetName(),samplDist->GetName())){
   obj->SetMarkerSize(size);
   break;
      }
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////

TH1F* SamplingDistPlot::GetTH1F(const SamplingDistribution *samplDist)
{
  if(samplDist == NULL){
    return fHist;
  }else{
    fIterator->Reset();
    TH1F *obj = 0;
    while((obj = (TH1F*)fIterator->Next())) {
      if(!strcmp(obj->GetName(),samplDist->GetName())){
        return obj;
      }
    }
  }

  return NULL;
}

////////////////////////////////////////////////////////////////////////////////

void SamplingDistPlot::RebinDistribution(Int_t rebinFactor, const SamplingDistribution *samplDist)
{
  if(samplDist == 0){
    fHist->Rebin(rebinFactor);
  }
  else{
    fIterator->Reset();
    TH1F *obj = 0;
    while((obj = (TH1F*)fIterator->Next())) {
      if(!strcmp(obj->GetName(),samplDist->GetName())){
   obj->Rebin(rebinFactor);
   break;
      }
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
/// TODO test

void SamplingDistPlot::DumpToFile(const char* RootFileName, Option_t *option, const char *ftitle, Int_t compress) {
   // All the objects are written to rootfile

   if(!fRooPlot) {
      cout << "Plot was not drawn yet. Dump can only be saved after it was drawn with Draw()." << endl;
      return;
   }

   TFile ofile(RootFileName, option, ftitle, compress);
   ofile.cd();
   fRooPlot->Write();
   ofile.Close();
}
