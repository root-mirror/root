// @(#)root/roostats:$Id$

/*************************************************************************
 * Project: RooStats                                                     *
 * Package: RooFit/RooStats                                              *
 * Authors:                                                              *
 *   Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke       *
 *************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


/** \class RooStats::LikelihoodIntervalPlot
    \ingroup Roostats

 This class provides simple and straightforward utilities to plot a LikelihoodInterval
 object.

*/

#include "RooStats/LikelihoodIntervalPlot.h"

#include <algorithm>
#include <iostream>
#include <cmath>

#include "TROOT.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TList.h"
#include "TGraph.h"
#include "TPad.h"
#include "TCanvas.h"
// need chisquare_quantile function - can use mathcore implementation
// for plotting not crucial that is less precise
#include "Math/QuantFuncMathCore.h"


#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include "RooProfileLL.h"
#include "TF1.h"

/// ClassImp for building the THtml documentation of the class
using namespace std;

ClassImp(RooStats::LikelihoodIntervalPlot);

using namespace RooStats;

////////////////////////////////////////////////////////////////////////////////
/// LikelihoodIntervalPlot default constructor
/// with default parameters

LikelihoodIntervalPlot::LikelihoodIntervalPlot()
{
  fInterval = 0;
  fNdimPlot = 0;
  fParamsPlot = 0;
  fColor = 0;
  fFillStyle = 4050; // half transparent
  fLineColor = 0;
  fMaximum = -1;
  fNPoints = 0;  // default depends if 1D or 2D
  // default is variable range
  fXmin = 0;
  fXmax = -1;
  fYmin = 0;
  fYmax = -1;
  fPrecision = -1; // use default
  fPlotObject = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// LikelihoodIntervalPlot copy constructor

LikelihoodIntervalPlot::LikelihoodIntervalPlot(LikelihoodInterval* theInterval)
{
  fInterval = theInterval;
  fParamsPlot = fInterval->GetParameters();
  fNdimPlot = fParamsPlot->getSize();
  fColor = 0;
  fLineColor = 0;
  fFillStyle = 4050; // half transparent
  fMaximum = -1;
  fNPoints = 0;  // default depends if 1D or 2D
  // default is variable range
  fXmin = 0;
  fXmax = -1;
  fYmin = 0;
  fYmax = -1;
  fPrecision = -1; // use default
  fPlotObject = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// LikelihoodIntervalPlot destructor

LikelihoodIntervalPlot::~LikelihoodIntervalPlot()
{
}

////////////////////////////////////////////////////////////////////////////////

void LikelihoodIntervalPlot::SetLikelihoodInterval(LikelihoodInterval* theInterval)
{
  fInterval = theInterval;
  fParamsPlot = fInterval->GetParameters();
  fNdimPlot = fParamsPlot->getSize();

  return;
}

////////////////////////////////////////////////////////////////////////////////

void LikelihoodIntervalPlot::SetPlotParameters(const RooArgSet *params)
{
  fNdimPlot = params->getSize();
  fParamsPlot = (RooArgSet*) params->clone((std::string(params->GetName())+"_clone").c_str());

  return;
}


////////////////////////////////////////////////////////////////////////////////
/// draw the log of the profiled likelihood function in 1D with the interval or
/// as a 2D plot with the contours.
/// Higher dimensional intervals cannot be drawn. One needs to call
/// SetPlotParameters to project interval in 1 or 2dim
///
/// ### Options for drawing 1D intervals
///
/// For 1D problem the log of the profiled likelihood function is drawn by default in a RooPlot as a
/// RooCurve
/// The plotting range (default is the full parameter range) and the precision of the RooCurve
/// can be specified by using SetRange(x1,x2) and SetPrecision(eps).
/// SetNPoints(npoints) can also be used  (default is npoints=100)
/// Optionally the function can be drawn as a TF1 (option="tf1") obtained by sampling the given npoints
/// in the given range
///
/// ### Options for drawing 2D intervals
///
/// For 2D case, a contour and optionally the profiled likelihood function is drawn by sampling npoints in
/// the given range. A 2d histogram of nbinsX=nbinsY = sqrt(npoints) is used for sampling the profiled likelihood.
/// The contour can be obtained by using Minuit or by the sampled histogram,
/// If using Minuit, the number of points specifies the number of contour points. If using an histogram the number of
/// points is approximately the total number of bins of the histogram.
/// Possible options:
/// -  minuit/nominuit:     use minuit for computing the contour
/// -   hist/nohist   :     sample in an histogram the profiled likelihood
///
/// Note that one can have both a drawing of the sampled likelihood and of the contour using minuit.
/// The default options is "minuit nohist"
/// The sampled histogram is drawn first by default using the option "colz" and then 8 probability contours at
/// these CL are drawn:  { 0.1,0.3,0.5,0.683,0.95,0.9973,0.9999366575,0.9999994267} re-drawing the histogram with the
/// option "cont3"
///
/// The drawn object (RooPlot or sampled histogram) is saved in the class and can be retrieved using GetPlottedObject()
/// In this way the user can eventually customize further the plot.
/// Note that the class does not delete the plotted object. It needs, if needed, to be deleted by the user

void LikelihoodIntervalPlot::Draw(const Option_t *options)
{
   TIter it = fParamsPlot->createIterator();
   // we need to check if parameters to plot is different than parameters of interval
   RooArgSet* intervalParams = fInterval->GetParameters();
   RooAbsArg * arg = 0;
   RooArgSet extraParams;
   while((arg=(RooAbsArg*)it.Next())) {
      if (!intervalParams->contains(*arg) ) {
         ccoutE(InputArguments) << "Parameter " << arg->GetName() << "is not in the list of LikelihoodInterval parameters "
                                << " - do not use for plotting " << std::endl;
         fNdimPlot--;
         extraParams.add(*arg);
      }
   }
   if (extraParams.getSize() > 0)
      fParamsPlot->remove(extraParams,true,true);

   if(fNdimPlot > 2){
      ccoutE(InputArguments) << "LikelihoodIntervalPlot::Draw(" << GetName()
                << ") ERROR: contours for more than 2 dimensions not implemented!" << std::endl;
      return;
   }

   // if the number of parameters to plot is less to the number of parameters of the LikelihoodInterval
   // we need to re-do the profile likelihood function, otherwise those parameters will not be profiled
   // when plotting
   RooAbsReal* newProfile = 0;
   RooAbsReal* oldProfile = fInterval->GetLikelihoodRatio();
   if (fNdimPlot != intervalParams->getSize() ) {
      RooProfileLL * profilell = dynamic_cast<RooProfileLL*>(oldProfile);
      if (!profilell) return;
      RooAbsReal & nll =  profilell->nll();
      newProfile = nll.createProfile(*fParamsPlot);
   }
   else {
      newProfile = oldProfile;
   }

   it.Reset();
   RooRealVar *myparam = (RooRealVar*) it.Next();

   // do a dummy evaluation around minimum to be sure profile has right minimum
   if (fInterval->GetBestFitParameters() ) {
      *fParamsPlot = *fInterval->GetBestFitParameters();
      newProfile->getVal();
   }

   // analyze options
   TString opt = options;
   opt.ToLower();

   TString title = GetTitle();
   int nPoints = fNPoints;

   if(fNdimPlot == 1) {

      // 1D drawing options
      // use RooPLot for drawing the 1D PL
      // if option is TF1 use TF1 for drawing
      bool useRooPlot = opt.Contains("rooplot") ||  ! (opt.Contains("tf1"));
      opt.ReplaceAll("rooplot","");
      opt.ReplaceAll("tf1","");


     //      if (title.Length() == 0)
     //         title = "- log profile likelihood ratio";

      if (nPoints <=0) nPoints = 100; // default in 1D

      const Double_t xcont_min = fInterval->LowerLimit(*myparam);
      const Double_t xcont_max = fInterval->UpperLimit(*myparam);

      RooRealVar* myarg = (RooRealVar *) newProfile->getVariables()->find(myparam->GetName());
      double x1 = myarg->getMin();
      double x2 = myarg->getMax();

      // default color values
      if (fColor == 0) fColor = kBlue;
      if (fLineColor == 0) fLineColor = kGreen;

      RooPlot * frame = 0;

      // use TF1 for drawing the function
      if (!useRooPlot) {

         // set a first estimate of range including 2 times upper and lower limit
         double xmin = std::max( x1, 2*xcont_min - xcont_max);
         double xmax = std::min( x2, 2*xcont_max - xcont_min);
         if (fXmin < fXmax) { xmin = fXmin; xmax = fXmax; }

         TF1 * tmp = newProfile->asTF(*myarg);
         assert(tmp != 0);
         tmp->SetRange(xmin, xmax);
         tmp->SetNpx(nPoints);

         // clone the function to avoid later to sample it
         TF1 * f1 = (TF1*) tmp->Clone();
         delete tmp;

         f1->SetTitle(title);
         TString name = TString(GetName()) + TString("_PLL_") + TString(myarg->GetName());
         f1->SetName(name);

         // set range for displaying x values where function <=  fMaximum
         // if no range is set amd
         // if no reasonable value found maintain first estimate
         x1 = xmin; x2 = xmax;
         if (fMaximum > 0 && fXmin >= fXmax ) {
            double x0 = f1->GetX(0, xmin, xmax);
            // check that minimum is between xmin and xmax
            if ( x0 > x1 && x0 < x2) {
               x1 = f1->GetX(fMaximum, xmin, x0);
               x2 = f1->GetX(fMaximum, x0, xmax);
               f1->SetMaximum(fMaximum);
            //std::cout << "setting range to " << x1 << " , " << x2 << " x0 = " << x0 << std::endl;
            }
         }

         f1->SetRange(x1,x2);


         f1->SetLineColor(kBlue);
         f1->GetXaxis()->SetTitle(myarg->GetName());
         f1->GetYaxis()->SetTitle(Form("- log #lambda(%s)",myparam->GetName()));
         f1->Draw(opt);
         fPlotObject = f1->GetHistogram();

      }
      else {
         // use a RooPlot for drawing the PL function
         double xmin = myparam->getMin(); double xmax =  myparam->getMax();
         if (fXmin < fXmax) { xmin = fXmin; xmax = fXmax; }

         // set nbins (must be used in combination with precision )
         // the curve will evaluate 2 * nbins if precision is > 1
         int prevBins = myarg->getBins();
         myarg->setBins(fNPoints);

         // want to set range on frame not function
         frame = myarg->frame(xmin,xmax,nPoints);
         // for ycutoff line
         x1= xmin;
         x2=xmax;
         frame->SetTitle(title);
         frame->GetYaxis()->SetTitle(Form("- log #lambda(%s)",myparam->GetName()));
         //    frame->GetYaxis()->SetTitle("- log profile likelihood ratio");


         // plot
         RooCmdArg cmd;
         if (fPrecision > 0) cmd = RooFit::Precision(fPrecision);
         newProfile->plotOn(frame,cmd,RooFit::LineColor(fColor));

         frame->SetMaximum(fMaximum);
         frame->SetMinimum(0.);

         myarg->setBins(prevBins);
         fPlotObject = frame;
      }


      //myarg->setVal(xcont_max);
      //const Double_t Yat_Xmax = newProfile->getVal();
      Double_t Yat_Xmax = 0.5*ROOT::Math::chisquared_quantile(fInterval->ConfidenceLevel(),1);

      TLine *Yline_cutoff = new TLine(x1,Yat_Xmax,x2,Yat_Xmax);
      TLine *Yline_min = new TLine(xcont_min,0.,xcont_min,Yat_Xmax);
      TLine *Yline_max = new TLine(xcont_max,0.,xcont_max,Yat_Xmax);

      Yline_cutoff->SetLineColor(fLineColor);
      Yline_min->SetLineColor(fLineColor);
      Yline_max->SetLineColor(fLineColor);

      if (!useRooPlot) {
         // need to draw the line
         Yline_cutoff->Draw();
         Yline_min->Draw();
         Yline_max->Draw();
      }
      else {
         // add line in the RooPlot
         frame->addObject(Yline_min);
         frame->addObject(Yline_max);
         frame->addObject(Yline_cutoff);
         frame->Draw(opt);
      }


      return;
   }

   // case of 2 dimensions

   else if(fNdimPlot == 2){

      //2D drawing options

      // use Minuit for drawing the contours of the PL (default case)
      bool useMinuit = !opt.Contains("nominuit");
      // plot histogram in 2D
      bool plotHist = !opt.Contains("nohist");
      opt.ReplaceAll("nominuit","");
      opt.ReplaceAll("nohist","");
      if (opt.Contains("minuit") ) useMinuit= true;
      if (useMinuit) plotHist = false; // switch off hist by default in case of Minuit
      if (opt.Contains("hist") ) plotHist= true;
      opt.ReplaceAll("minuit","");
      opt.ReplaceAll("hist","");

      RooRealVar *myparamY = (RooRealVar*)it.Next();

      Double_t cont_level = ROOT::Math::chisquared_quantile(fInterval->ConfidenceLevel(),fNdimPlot); // level for -2log LR
      cont_level = cont_level/2; // since we are plotting -log LR

      RooArgList params(*newProfile->getVariables());
      // set values and error for the POI to the best fit values
      for (int i = 0; i < params.getSize(); ++i) {
         RooRealVar & par =  (RooRealVar &) params[i];
         RooRealVar * fitPar =  (RooRealVar *) (fInterval->GetBestFitParameters()->find(par.GetName() ) );
         if (fitPar) {
            par.setVal( fitPar->getVal() );
         }
      }
      // do a profile evaluation to start from the best fit values of parameters
      newProfile->getVal();

      if (title.Length() == 0)
         title = TString("Contour of ") + TString(myparamY->GetName() ) + TString(" vs ") + TString(myparam->GetName() );
      // add also labels
      title = TString::Format("%s;%s;%s",title.Data(),myparam->GetName(),myparamY->GetName());

      if (nPoints <=0) nPoints = 40; // default in 2D

      double xmin = myparam->getMin(); double xmax =  myparam->getMax();
      double ymin = myparamY->getMin(); double ymax =  myparamY->getMax();
      if (fXmin < fXmax) { xmin = fXmin; xmax = fXmax; }
      if (fYmin < fYmax) { ymin = fYmin; ymax = fYmax; }


      if (!useMinuit || plotHist) {

         // find contour from a scanned histogram of points

         // draw directly the TH2 from the profile LL
         TString histName = TString::Format("_hist2D__%s_%s",myparam->GetName(),myparamY->GetName() );
         int nBins = int( std::sqrt(double(nPoints)) + 0.5 );
         TH2* hist2D = new TH2D(histName, title, nBins, xmin, xmax, nBins, ymin, ymax );
         newProfile->fillHistogram(hist2D, RooArgList(*myparam,*myparamY), 1, 0, false, 0, false);

         hist2D->SetTitle(title);
         hist2D->SetStats(kFALSE);

         //need many color levels for drawing with option colz
         if (plotHist) {

            const int nLevels = 51;
            double contLevels[nLevels];
            contLevels[0] = 0.01;
            double maxVal = (fMaximum > 0) ? fMaximum : hist2D->GetMaximum();
            for (int k = 1; k < nLevels; ++k) {
               contLevels[k] = k*maxVal/double(nLevels-1);
            }
            hist2D->SetContour(nLevels,contLevels);

            if (fMaximum>0) hist2D->SetMaximum(fMaximum);

            hist2D->DrawClone("COLZ");
         }


         //need now less contours for drawing with option cont

         const int nLevels = 8;
         double contLevels[nLevels];
         // last 3 are the 3,4,5 sigma levels
         double confLevels[nLevels] = { 0.1,0.3,0.5,0.683,0.95,0.9973,0.9999366575,0.9999994267};
         for (int k = 0; k < nLevels; ++k) {
            //contLevels[k] = 0.5*ROOT::Math::chisquared_quantile(1.-2.*ROOT::Math::normal_cdf_c(nSigmaLevels[k],1),2);
            contLevels[k] = 0.5*ROOT::Math::chisquared_quantile(confLevels[k],2);
         }
         hist2D->SetContour(nLevels,contLevels);
         if (fLineColor) hist2D->SetLineColor(fLineColor);

         // default options for drawing a second histogram
         TString tmpOpt = opt;
         tmpOpt.ReplaceAll("same","");
         if (tmpOpt.Length() < 3) opt += "cont3";
         // if histo is plotted draw on top
         if (plotHist) opt += TString(" same");
         hist2D->Draw(opt.Data());
         gPad->Update();

         // case of plotting contours without minuit
         if (!useMinuit) {

            // set levels of contours if make contours without minuit
            TH2 * h = (TH2*) hist2D->Clone();
            h->SetContour(1,&cont_level);

            TVirtualPad * currentPad = gPad;
            // o a temporary draw to get the contour graph
            TCanvas * tmpCanvas = new TCanvas("tmpCanvas","tmpCanvas");
            h->Draw("CONT LIST");
            gPad->Update();

            // get graphs from the contours
            TObjArray *contoursOrig = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
            // CLONE THE LIST IN CASE IT GETS DELETED
            TObjArray *contours = 0;
            if (contoursOrig) contours = (TObjArray*) contoursOrig->Clone();

            delete tmpCanvas;
            delete h;
            gPad = currentPad;


            // in case of option CONT4 I need to re-make the Pad
            if (tmpOpt.Contains("cont4")) {
               Double_t bm = gPad->GetBottomMargin();
               Double_t lm = gPad->GetLeftMargin();
               Double_t rm = gPad->GetRightMargin();
               Double_t tm = gPad->GetTopMargin();
               Double_t x1 = hist2D->GetXaxis()->GetXmin();
               Double_t y1 = hist2D->GetYaxis()->GetXmin();
               Double_t x2 = hist2D->GetXaxis()->GetXmax();
               Double_t y2 = hist2D->GetYaxis()->GetXmax();

               TPad *null=new TPad("null","null",0,0,1,1);
               null->SetFillStyle(0);
               null->SetFrameFillStyle(0);
               null->Draw();
               null->cd();
               null->Range(x1-(x2-x1)*(lm/(1-rm-lm)),
                           y1-(y2-y1)*(bm/(1-tm-lm)),
                           x2+(x2-x1)*(rm/(1-rm-lm)),
                           y2+(y2-y1)*(tm/(1-tm-lm)));

               gPad->Update();
            }


            if (contours) {
               int ncontours = contours->GetSize();
               for (int icont = 0; icont < ncontours; ++icont) {
                  TList *  contourList = (TList*)contours->At(icont);
                  if (contourList && contourList->GetSize() > 0) {
                     TIterator * itgr = contourList->MakeIterator();
                     TGraph * gr = 0;
                     while( (gr = dynamic_cast<TGraph*>(itgr->Next()) ) ){
                        if (fLineColor) gr->SetLineColor(fLineColor);
                        gr->SetLineStyle(kDashed);
                        gr->SetLineWidth(3);
                        if (fColor) {
                           gr->SetFillColor(fColor);
                           gr->Draw("FL");
                        }
                        else
                           gr->Draw("L");
                     }
                     delete itgr;
                  }
               }
            }
            else {
               ccoutE(InputArguments) << "LikelihoodIntervalPlot::Draw(" << GetName()
                                      << ") ERROR: no contours found in ListOfSpecial" << std::endl;
            }

            fPlotObject = hist2D;

         }
      }
      if (useMinuit) {

         // find contours  using Minuit
         TGraph * gr = new TGraph(nPoints+1);

         int ncp = fInterval->GetContourPoints(*myparam, *myparamY, gr->GetX(), gr->GetY(),nPoints);

         if (int(ncp) < nPoints) {
            std::cout << "Warning - Less points calculated in contours np = " << ncp << " / " << nPoints << std::endl;
            for (int i = ncp; i < nPoints; ++i) gr->RemovePoint(i);
         }
         // add last point to same as first one to close the contour
         gr->SetPoint(ncp, gr->GetX()[0], gr->GetY()[0] );
         if (!opt.Contains("c")) opt.Append("L");  // use by default option L if C is not specified
         // draw first a dummy 2d histogram gfor the axis
         if (!opt.Contains("same") && !plotHist) {

            TH2F* hist2D = new TH2F("_hist2D",title, nPoints, xmin, xmax, nPoints, ymin, ymax );
            hist2D->GetXaxis()->SetTitle(myparam->GetName());
            hist2D->GetYaxis()->SetTitle(myparamY->GetName());
            hist2D->SetBit(TH1::kNoStats); // do not draw statistics
            hist2D->SetFillStyle(fFillStyle);
            hist2D->SetMaximum(1);  // to avoid problem with subsequents draws
            hist2D->Draw("AXIS");
         }
         if (fLineColor) gr->SetLineColor(fLineColor);
         if (fColor) {
            // draw contour as filled area (add option "F")
            gr->SetFillColor(fColor);
            opt.Append("F");
         }
         gr->SetLineWidth(3);
         if (opt.Contains("same"))  gr->SetFillStyle(fFillStyle); // put transparent
         gr->Draw(opt);
         TString name = TString("Graph_of_") + TString(fInterval->GetName());
         gr->SetName(name);

         if (!fPlotObject) fPlotObject = gr;
         else if (fPlotObject->IsA() != TH2D::Class() ) fPlotObject = gr;

      }

      // draw also the minimum
      const RooArgSet * bestFitParams = fInterval->GetBestFitParameters();
      if (bestFitParams) {
         TGraph * gr0 = new TGraph(1);
         double x0 = bestFitParams->getRealValue(myparam->GetName());
         double y0 = bestFitParams->getRealValue(myparamY->GetName());
         gr0->SetPoint(0,x0,y0);
         gr0->SetMarkerStyle(33);
         if (fColor)  {
            if (fColor != kBlack) gr0->SetMarkerColor(fColor+4);
            else  gr0->SetMarkerColor(kGray);
         }
         gr0->Draw("P");
         delete bestFitParams;
      }



   }

   // need to delete if a new profileLL was made
   if (newProfile != oldProfile) delete newProfile;

   return;
}
