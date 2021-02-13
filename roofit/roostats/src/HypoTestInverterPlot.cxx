// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class RooStats::HypoTestInverterPlot
    \ingroup Roostats

Class to plot a HypoTestInverterResult, the output of the HypoTestInverter calculator.

It can be used to plot the obtained p-values ( CLb, CLs+b or CLs) for each scanned point, as well as
the test statistic distributions (when a calculator based on pseudo-experiments is used) for the two
hypotheses.

*/

// include header file of this class
#include "RooStats/HypoTestInverterPlot.h"

// include other header files
#include "RooStats/HybridResult.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/SamplingDistribution.h"

#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TLine.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TVirtualPad.h"
#include "TError.h"
#include "Math/DistFuncMathCore.h"

#include <cmath>

using namespace std;

ClassImp(RooStats::HypoTestInverterPlot);

using namespace RooStats;

////////////////////////////////////////////////////////////////////////////////
/// constructor from a HypoTestInverterResult class
/// name and title are taken from the result class

HypoTestInverterPlot::HypoTestInverterPlot(HypoTestInverterResult* results ) :
   TNamed( results->GetName(), results->GetTitle() ),
   fResults(results)
{
}

////////////////////////////////////////////////////////////////////////////////
/// constructor with name and title from a HypoTestInverterResult class

HypoTestInverterPlot::HypoTestInverterPlot( const char* name,
                   const char* title,
                   HypoTestInverterResult* results ) :
   TNamed( TString(name), TString(title) ),
   fResults(results)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Make the plot of the result of the scan using the observed data.
/// By default plot CLs or CLsb depending if the flag UseCLs is set for the results
/// that are passed to this instance.
///
/// \param opt Options according to following list:
///  - Empty: Return CLs or CLs+b depending on the value of UseCLs.ƒ
///  - "CLB":  return  CLb plot
///  - "CLS+B" / "CLSPLUSB": return  CLs+b plot  independently of the flag
///  - "CLS":  return  CLs plot  independently of the flag

TGraphErrors* HypoTestInverterPlot::MakePlot(Option_t * opt)
{
   TString option(opt);
   option.ToUpper();
   enum Mode_t { Default, CLb, CLsPlusb, CLs };
   Mode_t type = Default;
   if (option.Contains("CLB")) type = CLb;
   else if (option.Contains("CLS+B") || option.Contains("CLSPLUSB")) type = CLsPlusb;
   else if (option.Contains("CLS" )) type = CLs;

   const int nEntries = fResults->ArraySize();

   // sort the arrays based on the x values
   std::vector<unsigned int> index(nEntries);
   TMath::SortItr(fResults->fXValues.begin(), fResults->fXValues.end(), index.begin(), false);

   // copy result in sorted arrays
   std::vector<Double_t> xArray;
   std::vector<Double_t> yArray;
   std::vector<Double_t> yErrArray;

   for (int i=0; i<nEntries; i++) {
      double CLVal = 0., CLErr = 0.;
      if (type == Default) {
         CLVal = fResults->GetYValue(index[i]);
         CLErr = fResults->GetYError(index[i]);
      } else if (type == CLb) {
         CLVal = fResults->CLb(index[i]);
         CLErr = fResults->CLbError(index[i]);
      } else if (type == CLsPlusb) {
         CLVal = fResults->CLsplusb(index[i]);
         CLErr = fResults->CLsplusbError(index[i]);
      } else if (type == CLs) {
         CLVal = fResults->CLs(index[i]);
         CLErr = fResults->CLsError(index[i]);
      }

      if (CLVal < 0. || !std::isfinite(CLVal)) {
         Warning("HypoTestInverterPlot::MakePlot", "Got a confidence level of %f at x=%f (failed fit?). Skipping this point.", CLVal, fResults->GetXValue(index[i]));
         continue;
      }

      yArray.push_back(CLVal);
      yErrArray.push_back(CLErr);
      xArray.push_back(fResults->GetXValue(index[i]));
   }

   TGraphErrors* graph = new TGraphErrors(static_cast<Int_t>(xArray.size()),&xArray.front(),&yArray.front(),0,&yErrArray.front());
   TString pValueName = "CLs";
   if (type == CLb ) pValueName = "CLb";
   if (type == CLsPlusb || (type == Default && !fResults->fUseCLs) ) pValueName = "CLs+b";
   TString name = pValueName + TString("_observed");
   TString title = TString("Observed ") + pValueName;
   graph->SetName(name);
   graph->SetTitle(title);
   graph->SetMarkerStyle(20);
   graph->SetLineWidth(2);
   return graph;
}

////////////////////////////////////////////////////////////////////////////////
/// Make the expected plot and the bands
/// nsig1 and nsig2 indicates the n-sigma value for the bands
/// if nsig1 = 0 no band is drawn (only expected value)
/// if nsig2 > nsig1 (default is nsig1=1 and nsig2=2) the second band is also drawn
/// The first band is drawn in green while the second in yellow
/// THe return result is a TMultiGraph object

TMultiGraph* HypoTestInverterPlot::MakeExpectedPlot(double nsig1, double nsig2 )
{
   const int nEntries = fResults->ArraySize();
   bool doFirstBand = (nsig1 > 0);
   bool doSecondBand = (nsig2 > nsig1);

   nsig1 = std::abs(nsig1);
   nsig2 = std::abs(nsig2);

   // sort the arrays based on the x values
   std::vector<unsigned int> index(nEntries);
   TMath::SortItr(fResults->fXValues.begin(), fResults->fXValues.end(), index.begin(), false);

   // create the graphs
   TGraph * g0 = new TGraph;
   TString pValueName = "CLs";
   if (!fResults->fUseCLs) pValueName = "CLs+b";
   g0->SetTitle(TString::Format("Expected %s - Median",pValueName.Data()) );
   TGraphAsymmErrors * g1 = 0;
   TGraphAsymmErrors * g2 = 0;
   if (doFirstBand) {
      g1 = new TGraphAsymmErrors;
      if (nsig1 - int(nsig1) < 0.01)
         g1->SetTitle(TString::Format("Expected %s #pm %d #sigma",pValueName.Data(),int(nsig1)) );
      else
         g1->SetTitle(TString::Format("Expected %s #pm %3.1f #sigma",pValueName.Data(),nsig1) );
   }
   if (doSecondBand) {
      g2 = new TGraphAsymmErrors;
      if (nsig2 - int(nsig2) < 0.01)
         g2->SetTitle(TString::Format("Expected %s #pm %d #sigma",pValueName.Data(),int(nsig2)) );
      else
         g2->SetTitle(TString::Format("Expected %s #pm %3.1f #sigma",pValueName.Data(),nsig2) );
   }
   double p[7];
   double q[7];
   p[0] = ROOT::Math::normal_cdf(-nsig2);
   p[1] = ROOT::Math::normal_cdf(-nsig1);
   p[2] = 0.5;
   p[3] = ROOT::Math::normal_cdf(nsig1);
   p[4] = ROOT::Math::normal_cdf(nsig2);

   bool resultIsAsymptotic = ( !fResults->GetNullTestStatDist(0) && !fResults->GetAltTestStatDist(0) );
   int np = 0;
   for (int j=0; j<nEntries; ++j) {
      int i = index[j]; // i is the order index
      SamplingDistribution * s = fResults->GetExpectedPValueDist(i);
      if ( !s)  continue;
      const std::vector<double> & values = s->GetSamplingDistribution();
      // special case for asymptotic results (cannot use TMath::quantile in that case)
      if (resultIsAsymptotic) {
         double maxSigma = fResults->fgAsymptoticMaxSigma;
         double dsig = 2* maxSigma/ (values.size() -1) ;
         int  i0 = (int) TMath::Floor ( ( -nsig2 +  maxSigma )/dsig + 0.5);
         int  i1 = (int) TMath::Floor ( (-nsig1 +  maxSigma )/dsig + 0.5);
         int  i2 = (int) TMath::Floor ( ( maxSigma)/dsig + 0.5);
         int  i3 = (int) TMath::Floor ( ( nsig1 + maxSigma)/dsig + 0.5);
         int  i4 = (int) TMath::Floor ( ( nsig2 + maxSigma)/dsig + 0.5);
         q[0] = values[i0];
         q[1] = values[i1];
         q[2] = values[i2];
         q[3] = values[i3];
         q[4] = values[i4];
      }
      else {
         double * x = const_cast<double *>(&values[0]); // need to change TMath::Quantiles
         TMath::Quantiles(values.size(), 5, x,q,p,false);
      }

      g0->SetPoint(np, fResults->GetXValue(i),  q[2]);
      if (g1) {
         g1->SetPoint(np, fResults->GetXValue(i),  q[2]);
         g1->SetPointEYlow(np, q[2] - q[1]); // -1 sigma errorr
         g1->SetPointEYhigh(np, q[3] - q[2]);//+1 sigma error
      }
      if (g2) {
         g2->SetPoint(np, fResults->GetXValue(i), q[2]);

         g2->SetPointEYlow(np, q[2]-q[0]);   // -2 -- -1 sigma error
         g2->SetPointEYhigh(np, q[4]-q[2]);
      }
      if (s) delete s;
      np++;
   }

   TString name = GetName() + TString("_expected");
   TString title = TString("Expected ") + GetTitle();
   TMultiGraph* graph = new TMultiGraph(name,title);

   // set the graphics options and add in multi graph
   // orderof adding is drawing order
   if (g2) {
      g2->SetFillColor(kYellow);
      graph->Add(g2,"3");
   }
   if (g1) {
      g1->SetFillColor(kGreen);
      graph->Add(g1,"3");
   }
   g0->SetLineStyle(2);
   g0->SetLineWidth(2);
   graph->Add(g0,"L");

   return graph;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

HypoTestInverterPlot::~HypoTestInverterPlot()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the result in the current canvas
/// Possible options:
///  - SAME : draw in the current axis
///  - OBS  :  draw only the observed plot
///  - EXP  :  draw only the expected plot
///
///  - CLB  : draw also the CLB
///  - 2CL  : draw both clsplusb and cls
///
/// default draw observed + expected with 1 and 2 sigma bands

void HypoTestInverterPlot::Draw(Option_t * opt) {

   TString option(opt);
   option.ToUpper();
   bool drawAxis = !option.Contains("SAME");
   bool drawObs = option.Contains("OBS") || !option.Contains("EXP");
   bool drawExp = option.Contains("EXP") || !option.Contains("OBS");
   bool drawCLb = option.Contains("CLB");
   bool draw2CL = option.Contains("2CL");

   TGraphErrors * gobs = 0;
   TGraph * gplot = 0;
   if (drawObs) {
      gobs = MakePlot();
      // add object to top-level directory to avoid mem leak
      if (gROOT) gROOT->Add(gobs);
      if (drawAxis) {
         gobs->Draw("APL");
         gplot = gobs;
         gplot->GetHistogram()->SetTitle( GetTitle() );
      }
      else gobs->Draw("PL");

   }
   TMultiGraph * gexp = 0;
   if (drawExp) {
      gexp = MakeExpectedPlot();
      // add object to current directory to avoid mem leak
      if (gROOT) gROOT->Add(gexp);
      if (drawAxis && !drawObs) {
         gexp->Draw("A");
         if (gexp->GetHistogram()) gexp->GetHistogram()->SetTitle( GetTitle() );
         gplot = (TGraph*) gexp->GetListOfGraphs()->First();
      }
      else
         gexp->Draw();

   }

   // draw also an horizontal  line at the desired conf level
   if (gplot) {
      double alpha = 1.-fResults->ConfidenceLevel();
      double x1 = gplot->GetXaxis()->GetXmin();
      double x2 = gplot->GetXaxis()->GetXmax();
      TLine * line = new TLine(x1, alpha, x2,alpha);
      line->SetLineColor(kRed);
      line->Draw();
      // put axis labels
      RooAbsArg * arg = fResults->fParameters.first();
      if (arg) gplot->GetXaxis()->SetTitle(arg->GetName());
      gplot->GetYaxis()->SetTitle("p value");
   }


   TGraph *gclb = 0;
   if (drawCLb) {
      gclb = MakePlot("CLb");
      if (gROOT) gROOT->Add(gclb);
      gclb->SetMarkerColor(kBlue+4);
      gclb->Draw("PL");
      // draw in red observed cls or clsb
      if (gobs) gobs->SetMarkerColor(kRed);
   }
   TGraph * gclsb = 0;
   TGraph * gcls = 0;
   if (draw2CL) {
      if (fResults->fUseCLs) {
         gclsb = MakePlot("CLs+b");
         if (gROOT) gROOT->Add(gclsb);
         gclsb->SetMarkerColor(kBlue);
         gclsb->Draw("PL");
         gclsb->SetLineStyle(3);
      }
      else {
         gcls = MakePlot("CLs");
         if (gROOT) gROOT->Add(gcls);
         gcls->SetMarkerColor(kBlue);
         gcls->Draw("PL");
         gcls->SetLineStyle(3);
      }
   }
   // draw again observed values otherwise will be covered by the bands
   if (gobs) {
      gobs->Draw("PL");
   }


   double y0 = 0.6;
   double verticalSize = (gexp || draw2CL || drawCLb ) ? 0.3 : 0.15;
   double y1 = y0 + verticalSize;
   TLegend * l = new TLegend(0.6,y0,0.9,y1);
   if (gobs) l->AddEntry(gobs,"","PEL");
   if (gclsb) l->AddEntry(gclsb,"","PEL");
   if (gcls) l->AddEntry(gcls,"","PEL");
   if (gclb) l->AddEntry(gclb,"","PEL");
   if (gexp) {
      // loop in reverse order (opposite to drawing one)
      int ngraphs =  gexp->GetListOfGraphs()->GetSize();
      for (int i = ngraphs-1; i>=0; --i) {
         TObject * obj =  gexp->GetListOfGraphs()->At(i);
         TString lopt = "F";
         if (i == ngraphs-1) lopt = "L";
         if (obj)  l->AddEntry(obj,"",lopt);
      }
   }
   l->Draw();
   // redraw the axis
   if (gPad) gPad->RedrawAxis();

}

////////////////////////////////////////////////////////////////////////////////
/// \param index Index of the result stored in HypoTestInverterResult
/// \param type Type of the test (see below)
/// \param nbins Number of bins
///  - type =0  null and alt
///  - type = 1 only null (S+B)
///  - type = 2 only alt  (B)

SamplingDistPlot * HypoTestInverterPlot::MakeTestStatPlot(int index, int type, int nbins) {
   SamplingDistPlot * pl = 0;
   if (type == 0) {
      HypoTestResult * result = (HypoTestResult*) fResults->fYObjects.At(index);
      if (result)
         pl = new HypoTestPlot(*result, nbins );
      return pl;
   }
   if (type == 1) {
      SamplingDistribution * sbDist = fResults->GetSignalAndBackgroundTestStatDist(index);
      if (sbDist) {
         pl = new SamplingDistPlot( nbins);
         pl->AddSamplingDistribution(sbDist);
         return pl;
      }
   }
   if (type == 2) {
      SamplingDistribution * bDist = fResults->GetBackgroundTestStatDist(index);
      if (bDist) {
         pl = new SamplingDistPlot( nbins);
         pl->AddSamplingDistribution(bDist);
         return pl;
      }
   }
   return 0;
}
