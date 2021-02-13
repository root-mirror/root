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

#ifndef ROOSTATS_SamplingDistPlot
#define ROOSTATS_SamplingDistPlot

#include "Compression.h"
#include "RooList.h"
#include "RooPrintable.h"
#include "TNamed.h"
#include "TIterator.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include <vector>

#include "RooStats/SamplingDistribution.h"

#include "RooPlot.h"


namespace RooStats {

 class SamplingDistPlot : public TNamed, public RooPrintable {

   public:
    /// Constructors for SamplingDistribution
    SamplingDistPlot(Int_t nbins = 100);
    SamplingDistPlot(Int_t nbins, Double_t min, Double_t max);
//    SamplingDistPlot(const char* name, const char* title, Int_t nbins, Double_t xmin, Double_t xmax);

    /// Destructor of SamplingDistribution
    virtual ~SamplingDistPlot();

    /// adds the sampling distribution and returns the scale factor
    Double_t AddSamplingDistribution(const SamplingDistribution *samplingDist, Option_t *drawOptions="NORMALIZE HIST");
    /// Like AddSamplingDistribution, but also sets a shaded area in the
    /// minShaded and maxShaded boundaries.
    Double_t AddSamplingDistributionShaded(const SamplingDistribution *samplingDist, Double_t minShaded, Double_t maxShaded, Option_t *drawOptions="NORMALIZE HIST");

    /// add a line
    void AddLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2, const char* title = NULL);
    /// add a TH1
    void AddTH1(TH1* h, Option_t *drawOptions="");
    /// add a TF1
    void AddTF1(TF1* f, const char* title = NULL, Option_t *drawOptions="SAME");
    /// set legend
    void SetLegend(TLegend* l){ fLegend = l; }

    void Draw(Option_t *options=0);

    /// Applies a predefined style if fApplyStyle is kTRUE (default).
    void ApplyDefaultStyle(void);

    void SetLineColor(Color_t color, const SamplingDistribution *samplDist = 0);
    void SetLineWidth(Width_t lwidth, const SamplingDistribution *samplDist = 0);
    void SetLineStyle(Style_t style, const SamplingDistribution *samplDist = 0);

    void SetMarkerColor(Color_t color, const SamplingDistribution *samplDist = 0);
    void SetMarkerStyle(Style_t style, const SamplingDistribution *samplDist = 0);
    void SetMarkerSize(Size_t size, const SamplingDistribution *samplDist = 0);

    void RebinDistribution(Int_t rebinFactor, const SamplingDistribution *samplDist = 0);

    void SetAxisTitle(char *varName) { fVarName = TString(varName); }

    /// If you do not want SamplingDistPlot to interfere with your style settings, call this
    /// function with "false" before Draw().
    void SetApplyStyle(Bool_t s) { fApplyStyle = s; }

    /// Returns the TH1F associated with the give SamplingDistribution.
    /// Intended use: Access to member functions of TH1F like GetMean(),
    /// GetRMS() etc.
    /// The return objects is managed by  SamplingDistPlot
    TH1F* GetTH1F(const SamplingDistribution *samplDist = NULL);
    TH1 * GetHistogram(const SamplingDistribution *samplDist = NULL) { return GetTH1F(samplDist); }

    /// return plotter class used to draw the sampling distribution histograms
    /// object is managed by SamplingDistPlot
    RooPlot * GetPlot() { return fRooPlot; }

    /// changes plot to log scale on x axis
    void SetLogXaxis(Bool_t lx) { fLogXaxis = lx; }
    /// changes plot to log scale on y axis
    void SetLogYaxis(Bool_t ly) { fLogYaxis = ly; }

    /// change x range
    void SetXRange( double mi, double ma ) { fXMin = mi; fXMax = ma; }
    /// change y range
    void SetYRange( double mi, double ma ) { fYMin = mi; fYMax = ma; }

    /// write to Root file
    void DumpToFile(const char* RootFileName, Option_t *option="", const char *ftitle="", Int_t compress = ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault);

  private:
    std::vector<Double_t> fSamplingDistr;
    std::vector<Double_t> fSampleWeights;

    Bool_t fIsWeighted;

    Int_t fBins;
    Int_t fMarkerType;
    Int_t fColor;

    TString fVarName;

  protected:

    TH1F* fHist;
    TLegend *fLegend;

    RooList fItems; /// holds TH1Fs only
    RooList fOtherItems; /// other objects to be drawn like TLine etc.
    TIterator* fIterator; /// TODO remove class variable and instantiate locally as necessary
    RooPlot* fRooPlot;
    Bool_t fLogXaxis;
    Bool_t fLogYaxis;

    double fXMin, fXMax, fYMin, fYMax;

    Bool_t fApplyStyle;
    Style_t fFillStyle;

    void SetSampleWeights(const SamplingDistribution *samplingDist);

    void addObject(TObject *obj, Option_t *drawOptions=0); // for TH1Fs only
    void addOtherObject(TObject *obj, Option_t *drawOptions=0);
    void GetAbsoluteInterval(Double_t &theMin, Double_t &theMax, Double_t &theYMax) const;

    ClassDef(SamplingDistPlot,1)  /// Class containing the results of the HybridCalculator
  };
}

#endif
