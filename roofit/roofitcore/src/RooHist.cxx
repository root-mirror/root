/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
\file RooHist.cxx
\class RooHist
\ingroup Roofitcore

A RooHist is a graphical representation of binned data based on the
TGraphAsymmErrors class. Error bars are calculated using either Poisson
or Binomial statistics. A RooHist is used to represent histograms in
a RooPlot.
**/

#include "RooFit.h"

#include "RooHist.h"
#include "RooHistError.h"
#include "RooCurve.h"
#include "RooScaledFunc.h"
#include "RooMsgService.h"

#include "TH1.h"
#include "TClass.h"
#include "Riostream.h"
#include <iomanip>

using namespace std;

ClassImp(RooHist);
  ;


////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooHist::RooHist() :
  _nominalBinWidth(1),
  _nSigma(1),
  _entries(0),
  _rawEntries(0)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Create an empty histogram that can be filled with the addBin()
/// and addAsymmetryBin() methods. Use the optional parameter to
/// specify the confidence level in units of sigma to use for
/// calculating error bars. The nominal bin width specifies the
/// default used by addBin(), and is used to set the relative
/// normalization of bins with different widths.

  RooHist::RooHist(Double_t nominalBinWidth, Double_t nSigma, Double_t /*xErrorFrac*/, Double_t /*scaleFactor*/) :
    TGraphAsymmErrors(), _nominalBinWidth(nominalBinWidth), _nSigma(nSigma), _rawEntries(-1)
{
  initialize();
}


////////////////////////////////////////////////////////////////////////////////
/// Create a histogram from the contents of the specified TH1 object
/// which may have fixed or variable bin widths. Error bars are
/// calculated using Poisson statistics. Prints a warning and rounds
/// any bins with non-integer contents. Use the optional parameter to
/// specify the confidence level in units of sigma to use for
/// calculating error bars. The nominal bin width specifies the
/// default used by addBin(), and is used to set the relative
/// normalization of bins with different widths. If not set, the
/// nominal bin width is calculated as range/nbins.

RooHist::RooHist(const TH1 &data, Double_t nominalBinWidth, Double_t nSigma, RooAbsData::ErrorType etype, Double_t xErrorFrac,
       Bool_t correctForBinWidth, Double_t scaleFactor) :
  TGraphAsymmErrors(), _nominalBinWidth(nominalBinWidth), _nSigma(nSigma), _rawEntries(-1)
{
  initialize();
  // copy the input histogram's name and title
  SetName(data.GetName());
  SetTitle(data.GetTitle());
  // calculate our nominal bin width if necessary
  if(_nominalBinWidth == 0) {
    const TAxis *axis= ((TH1&)data).GetXaxis();
    if(axis->GetNbins() > 0) _nominalBinWidth= (axis->GetXmax() - axis->GetXmin())/axis->GetNbins();
  }
  setYAxisLabel(data.GetYaxis()->GetTitle());

  // initialize our contents from the input histogram's contents
  Int_t nbin= data.GetNbinsX();
  for(Int_t bin= 1; bin <= nbin; bin++) {
    Axis_t x= data.GetBinCenter(bin);
    Stat_t y= data.GetBinContent(bin);
    Stat_t dy = data.GetBinError(bin) ;
    if (etype==RooAbsData::Poisson) {
      addBin(x,y,data.GetBinWidth(bin),xErrorFrac,scaleFactor);
    } else if (etype==RooAbsData::SumW2) {
      addBinWithError(x,y,dy,dy,data.GetBinWidth(bin),xErrorFrac,correctForBinWidth,scaleFactor);
    } else {
      addBinWithError(x,y,0,0,data.GetBinWidth(bin),xErrorFrac,correctForBinWidth,scaleFactor);
    }
  }
  // add over/underflow bins to our event count
  _entries+= data.GetBinContent(0) + data.GetBinContent(nbin+1);
}



////////////////////////////////////////////////////////////////////////////////
/// Create a histogram from the asymmetry between the specified TH1 objects
/// which may have fixed or variable bin widths, but which must both have
/// the same binning. The asymmetry is calculated as (1-2)/(1+2). Error bars are
/// calculated using Binomial statistics. Prints a warning and rounds
/// any bins with non-integer contents. Use the optional parameter to
/// specify the confidence level in units of sigma to use for
/// calculating error bars. The nominal bin width specifies the
/// default used by addAsymmetryBin(), and is used to set the relative
/// normalization of bins with different widths. If not set, the
/// nominal bin width is calculated as range/nbins.

RooHist::RooHist(const TH1 &data1, const TH1 &data2, Double_t nominalBinWidth, Double_t nSigma,
       RooAbsData::ErrorType etype, Double_t xErrorFrac, Bool_t efficiency, Double_t scaleFactor) :
  TGraphAsymmErrors(), _nominalBinWidth(nominalBinWidth), _nSigma(nSigma), _rawEntries(-1)
{
  initialize();
  // copy the first input histogram's name and title
  SetName(data1.GetName());
  SetTitle(data1.GetTitle());
  // calculate our nominal bin width if necessary
  if(_nominalBinWidth == 0) {
    const TAxis *axis= ((TH1&)data1).GetXaxis();
    if(axis->GetNbins() > 0) _nominalBinWidth= (axis->GetXmax() - axis->GetXmin())/axis->GetNbins();
  }

  if (!efficiency) {
    setYAxisLabel(Form("Asymmetry (%s - %s)/(%s + %s)",
           data1.GetName(),data2.GetName(),data1.GetName(),data2.GetName()));
  } else {
    setYAxisLabel(Form("Efficiency (%s)/(%s + %s)",
           data1.GetName(),data1.GetName(),data2.GetName()));
  }
  // initialize our contents from the input histogram contents
  Int_t nbin= data1.GetNbinsX();
  if(data2.GetNbinsX() != nbin) {
    coutE(InputArguments) << "RooHist::RooHist: histograms have different number of bins" << endl;
    return;
  }
  for(Int_t bin= 1; bin <= nbin; bin++) {
    Axis_t x= data1.GetBinCenter(bin);
    if(fabs(data2.GetBinCenter(bin)-x)>1e-10) {
      coutW(InputArguments) << "RooHist::RooHist: histograms have different centers for bin " << bin << endl;
    }
    Stat_t y1= data1.GetBinContent(bin);
    Stat_t y2= data2.GetBinContent(bin);
    if (!efficiency) {

      if (etype==RooAbsData::Poisson) {
   addAsymmetryBin(x,roundBin(y1),roundBin(y2),data1.GetBinWidth(bin),xErrorFrac,scaleFactor);
      } else if (etype==RooAbsData::SumW2) {
   Stat_t dy1= data1.GetBinError(bin);
   Stat_t dy2= data2.GetBinError(bin);
   addAsymmetryBinWithError(x,y1,y2,dy1,dy2,data1.GetBinWidth(bin),xErrorFrac,scaleFactor);
      } else {
   addAsymmetryBinWithError(x,y1,y2,0,0,data1.GetBinWidth(bin),xErrorFrac,scaleFactor);
      }

    } else {

      if (etype==RooAbsData::Poisson) {
   addEfficiencyBin(x,roundBin(y1),roundBin(y2),data1.GetBinWidth(bin),xErrorFrac,scaleFactor);
      } else if (etype==RooAbsData::SumW2) {
   Stat_t dy1= data1.GetBinError(bin);
   Stat_t dy2= data2.GetBinError(bin);
   addEfficiencyBinWithError(x,y1,y2,dy1,dy2,data1.GetBinWidth(bin),xErrorFrac,scaleFactor);
      } else {
   addEfficiencyBinWithError(x,y1,y2,0,0,data1.GetBinWidth(bin),xErrorFrac,scaleFactor);
      }

    }

  }
  // we do not have a meaningful number of entries
  _entries= -1;
}



////////////////////////////////////////////////////////////////////////////////
/// Create histogram as sum of two existing histograms. If Poisson errors are selected the histograms are
/// added and Poisson confidence intervals are calculated for the summed content. If wgt1 and wgt2 are not
/// 1 in this mode, a warning message is printed. If SumW2 errors are selected the histograms are added
/// and the histograms errors are added in quadrature, taking the weights into account.

RooHist::RooHist(const RooHist& hist1, const RooHist& hist2, Double_t wgt1, Double_t wgt2,
       RooAbsData::ErrorType etype, Double_t xErrorFrac) : _rawEntries(-1)
{
  // Initialize the histogram
  initialize() ;

  // Copy all non-content properties from hist1
  SetName(hist1.GetName()) ;
  SetTitle(hist1.GetTitle()) ;
  _nominalBinWidth=hist1._nominalBinWidth ;
  _nSigma=hist1._nSigma ;
  setYAxisLabel(hist1.getYAxisLabel()) ;

  if (!hist1.hasIdenticalBinning(hist2)) {
    coutE(InputArguments) << "RooHist::RooHist input histograms have incompatible binning, combined histogram will remain empty" << endl ;
    return ;
  }

  if (etype==RooAbsData::Poisson) {
    // Add histograms with Poisson errors

    // Issue warning if weights are not 1
    if (wgt1!=1.0 || wgt2 != 1.0) {
      coutW(InputArguments) << "RooHist::RooHist: WARNING: Poisson errors of weighted sum of two histograms is not well defined! " << endl
             << "                  Summed histogram bins will rounded to nearest integer for Poisson confidence interval calculation" << endl ;
    }

    // Add histograms, calculate Poisson confidence interval on sum value
    Int_t i,n=hist1.GetN() ;
    for(i=0 ; i<n ; i++) {
      Double_t x1,y1,x2,y2,dx1 ;
      hist1.GetPoint(i,x1,y1) ;
      dx1 = hist1.GetErrorX(i) ;
      hist2.GetPoint(i,x2,y2) ;
      addBin(x1,roundBin(wgt1*y1+wgt2*y2),2*dx1/xErrorFrac,xErrorFrac) ;
    }

  } else {
    // Add histograms with SumW2 errors

    // Add histograms, calculate combined sum-of-weights error
    Int_t i,n=hist1.GetN() ;
    for(i=0 ; i<n ; i++) {
      Double_t x1,y1,x2,y2,dx1,dy1,dy2 ;
      hist1.GetPoint(i,x1,y1) ;
      dx1 = hist1.GetErrorX(i) ;
      dy1 = hist1.GetErrorY(i) ;
      dy2 = hist2.GetErrorY(i) ;
      hist2.GetPoint(i,x2,y2) ;
      Double_t dy = sqrt(wgt1*wgt1*dy1*dy1+wgt2*wgt2*dy2*dy2) ;
      addBinWithError(x1,wgt1*y1+wgt2*y2,dy,dy,2*dx1/xErrorFrac,xErrorFrac) ;
    }
  }

}


////////////////////////////////////////////////////////////////////////////////
/// Create histogram from a pdf or function. Errors are computed based on the fit result provided.
///
/// This signature is intended for unfolding/deconvolution scenarios,
/// where a pdf is constructed as "data minus background" and is thus
/// intended to be displayed as "data" (or at least data-like).
/// Usage of this signature is triggered by the draw style "P" in RooAbsReal::plotOn.
/// 
/// More details.
/// \param[in] f The function to be plotted.
/// \param[in] x The variable on the x-axis
/// \param[in] xErrorFrac Size of the errror in x as a fraction of the bin width
/// \param[in] scaleFactor arbitrary scaling of the y-values
/// \param[in] normVars variables over which to normalize
RooHist::RooHist(const RooAbsReal &f, RooAbsRealLValue &x, Double_t xErrorFrac, Double_t scaleFactor, const RooArgSet *normVars, const RooFitResult* fr) :
  TGraphAsymmErrors(), _nSigma(1), _rawEntries(-1)
{
  // grab the function's name and title
  TString name(f.GetName());
  SetName(name.Data());
  TString title(f.GetTitle());
  SetTitle(title.Data());
  // append " ( [<funit> ][/ <xunit> ])" to our y-axis label if necessary
  if(0 != strlen(f.getUnit()) || 0 != strlen(x.getUnit())) {
    title.Append(" ( ");
    if(0 != strlen(f.getUnit())) {
      title.Append(f.getUnit());
      title.Append(" ");
    }
    if(0 != strlen(x.getUnit())) {
      title.Append("/ ");
      title.Append(x.getUnit());
      title.Append(" ");
    }
    title.Append(")");
  }
  setYAxisLabel(title.Data());

  RooAbsFunc *funcPtr = nullptr;
  RooAbsFunc *rawPtr  = nullptr;
  funcPtr= f.bindVars(x,normVars,kTRUE);

  // apply a scale factor if necessary
  if(scaleFactor != 1) {
    rawPtr= funcPtr;
    funcPtr= new RooScaledFunc(*rawPtr,scaleFactor);
  }
  
  // apply a scale factor if necessary
  assert(funcPtr);

  // calculate the points to add to our curve
  int xbins = x.numBins();
  RooArgSet nset;
  if(normVars) nset.add(*normVars);
  for(int i=0; i<xbins; ++i){
    double xval = x.getBinning().binCenter(i);
    double xwidth = x.getBinning().binWidth(i);
    Axis_t xval_ax = xval;
    double yval = (*funcPtr)(&xval);
    double yerr = sqrt(yval);
    if(fr) yerr = f.getPropagatedError(*fr,nset);
    addBinWithError(xval_ax,yval,yerr,yerr,xwidth,xErrorFrac,false,scaleFactor) ;
    _entries += yval;
  }
  _nominalBinWidth = 1.;
  
  // cleanup
  delete funcPtr;
  if(rawPtr) delete rawPtr;
}


////////////////////////////////////////////////////////////////////////////////
/// Perform common initialization for all constructors.

void RooHist::initialize()
{
  SetMarkerStyle(8);
  _entries= 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the number of events of the dataset associated with this RooHist.
/// This is the number of events in the RooHist itself, unless a different
/// value was specified through setRawEntries()

Double_t RooHist::getFitRangeNEvt() const
{
  return (_rawEntries==-1 ? _entries : _rawEntries) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate integral of histogram in given range

Double_t RooHist::getFitRangeNEvt(Double_t xlo, Double_t xhi) const
{
  Double_t sum(0) ;
  for (int i=0 ; i<GetN() ; i++) {
    Double_t x,y ;

    GetPoint(i,x,y) ;

    if (x>=xlo && x<=xhi) {
      sum += y ;
    }
  }

  if (_rawEntries!=-1) {
    coutW(Plotting) << "RooHist::getFitRangeNEvt() WARNING: The number of normalisation events associated to histogram " << GetName() << " is not equal to number of events in this histogram."
          << "\n\t\t This is due a cut being applied while plotting the data. Automatic normalisation over a sub-range of a plot variable assumes"
          << "\n\t\t that the effect of that cut is uniform across the plot, which may be an incorrect assumption. To obtain a correct normalisation, it needs to be passed explicitly:"
          << "\n\t\t\t data->plotOn(frame01,CutRange(\"SB1\"));"
          << "\n\t\t\t const double nData = data->sumEntries(\"\", \"SB1\"); //or the cut string such as sumEntries(\"x > 0.\");"
          << "\n\t\t\t model.plotOn(frame01, RooFit::Normalization(nData, RooAbsReal::NumEvent), ProjectionRange(\"SB1\"));" << endl ;
    sum *= _rawEntries / _entries ;
  }

  return sum ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return (average) bin width of this RooHist

Double_t RooHist::getFitRangeBinW() const
{
  return _nominalBinWidth ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the nearest positive integer to the input value
/// and print a warning if an adjustment is required.

Int_t RooHist::roundBin(Double_t y)
{
  if(y < 0) {
    coutW(Plotting) << fName << "::roundBin: rounding negative bin contents to zero: " << y << endl;
    return 0;
  }
  Int_t n= (Int_t)(y+0.5);
  if(fabs(y-n)>1e-6) {
    coutW(Plotting) << fName << "::roundBin: rounding non-integer bin contents: " << y << endl;
  }
  return n;
}



////////////////////////////////////////////////////////////////////////////////
/// Add a bin to this histogram with the specified integer bin contents
/// and using an error bar calculated with Poisson statistics. The bin width
/// is used to set the relative scale of bins with different widths.

void RooHist::addBin(Axis_t binCenter, Double_t n, Double_t binWidth, Double_t xErrorFrac, Double_t scaleFactor)
{
  if (n<0) {
    coutW(Plotting) << "RooHist::addBin(" << GetName() << ") WARNING: negative entry set to zero when Poisson error bars are requested" << endl ;
  }

  Double_t scale= 1;
  if(binWidth > 0) {
    scale= _nominalBinWidth/binWidth;
  }
  _entries+= n;
  Int_t index= GetN();

  // calculate Poisson errors for this bin
  Double_t ym,yp,dx(0.5*binWidth);

  if (fabs((double)((n-Int_t(n))>1e-5))) {
    // need interpolation
    Double_t ym1(0),yp1(0),ym2(0),yp2(0) ;
    Int_t n1 = Int_t(n) ;
    Int_t n2 = n1+1 ;
    if(!RooHistError::instance().getPoissonInterval(n1,ym1,yp1,_nSigma) ||
       !RooHistError::instance().getPoissonInterval(n2,ym2,yp2,_nSigma)) {
      coutE(Plotting) << "RooHist::addBin: unable to add bin with " << n << " events" << endl;
    }
    ym = ym1 + (n-n1)*(ym2-ym1) ;
    yp = yp1 + (n-n1)*(yp2-yp1) ;
    coutW(Plotting) << "RooHist::addBin(" << GetName()
          << ") WARNING: non-integer bin entry " << n << " with Poisson errors, interpolating between Poisson errors of adjacent integer" << endl ;
  } else {
  // integer case
  if(!RooHistError::instance().getPoissonInterval(Int_t(n),ym,yp,_nSigma)) {
      coutE(Plotting) << "RooHist::addBin: unable to add bin with " << n << " events" << endl;
      return;
    }
  }

  SetPoint(index,binCenter,n*scale*scaleFactor);
  SetPointError(index,dx*xErrorFrac,dx*xErrorFrac,scale*(n-ym)*scaleFactor,scale*(yp-n)*scaleFactor);
  updateYAxisLimits(scale*yp);
  updateYAxisLimits(scale*ym);
}



////////////////////////////////////////////////////////////////////////////////
/// Add a bin to this histogram with the specified bin contents
/// and error. The bin width is used to set the relative scale of
/// bins with different widths.

void RooHist::addBinWithError(Axis_t binCenter, Double_t n, Double_t elow, Double_t ehigh, Double_t binWidth,
               Double_t xErrorFrac, Bool_t correctForBinWidth, Double_t scaleFactor)
{
  Double_t scale= 1;
  if(binWidth > 0 && correctForBinWidth) {
    scale= _nominalBinWidth/binWidth;
  }
  _entries+= n;
  Int_t index= GetN();

  Double_t dx(0.5*binWidth) ;
  SetPoint(index,binCenter,n*scale*scaleFactor);
  SetPointError(index,dx*xErrorFrac,dx*xErrorFrac,elow*scale*scaleFactor,ehigh*scale*scaleFactor);
  updateYAxisLimits(scale*(n-elow));
  updateYAxisLimits(scale*(n+ehigh));
}




////////////////////////////////////////////////////////////////////////////////
/// Add a bin to this histogram with the specified bin contents
/// and error. The bin width is used to set the relative scale of
/// bins with different widths.

void RooHist::addBinWithXYError(Axis_t binCenter, Double_t n, Double_t exlow, Double_t exhigh, Double_t eylow, Double_t eyhigh,
            Double_t scaleFactor)
{
  _entries+= n;
  Int_t index= GetN();

  SetPoint(index,binCenter,n*scaleFactor);
  SetPointError(index,exlow,exhigh,eylow*scaleFactor,eyhigh*scaleFactor);
  updateYAxisLimits(scaleFactor*(n-eylow));
  updateYAxisLimits(scaleFactor*(n+eyhigh));
}





////////////////////////////////////////////////////////////////////////////////
/// Add a bin to this histogram with the value (n1-n2)/(n1+n2)
/// using an error bar calculated with Binomial statistics.

void RooHist::addAsymmetryBin(Axis_t binCenter, Int_t n1, Int_t n2, Double_t binWidth, Double_t xErrorFrac, Double_t scaleFactor)
{
  Double_t scale= 1;
  if(binWidth > 0) scale= _nominalBinWidth/binWidth;
  Int_t index= GetN();

  // calculate Binomial errors for this bin
  Double_t ym,yp,dx(0.5*binWidth);
  if(!RooHistError::instance().getBinomialIntervalAsym(n1,n2,ym,yp,_nSigma)) {
    coutE(Plotting) << "RooHist::addAsymmetryBin: unable to calculate binomial error for bin with " << n1 << "," << n2 << " events" << endl;
    return;
  }

  Double_t a= (Double_t)(n1-n2)/(n1+n2);
  SetPoint(index,binCenter,a*scaleFactor);
  SetPointError(index,dx*xErrorFrac,dx*xErrorFrac,(a-ym)*scaleFactor,(yp-a)*scaleFactor);
  updateYAxisLimits(scale*yp);
  updateYAxisLimits(scale*ym);
}



////////////////////////////////////////////////////////////////////////////////
/// Add a bin to this histogram with the value (n1-n2)/(n1+n2)
/// using an error bar calculated with Binomial statistics.

void RooHist::addAsymmetryBinWithError(Axis_t binCenter, Double_t n1, Double_t n2, Double_t en1, Double_t en2, Double_t binWidth, Double_t xErrorFrac, Double_t scaleFactor)
{
  Double_t scale= 1;
  if(binWidth > 0) scale= _nominalBinWidth/binWidth;
  Int_t index= GetN();

  // calculate Binomial errors for this bin
  Double_t ym,yp,dx(0.5*binWidth);
  Double_t a= (Double_t)(n1-n2)/(n1+n2);

  Double_t error = 2*sqrt( pow(en1,2)*pow(n2,2) + pow(en2,2)*pow(n1,2) ) / pow(n1+n2,2) ;
  ym=a-error ;
  yp=a+error ;

  SetPoint(index,binCenter,a*scaleFactor);
  SetPointError(index,dx*xErrorFrac,dx*xErrorFrac,(a-ym)*scaleFactor,(yp-a)*scaleFactor);
  updateYAxisLimits(scale*yp);
  updateYAxisLimits(scale*ym);
}



////////////////////////////////////////////////////////////////////////////////
/// Add a bin to this histogram with the value n1/(n1+n2)
/// using an error bar calculated with Binomial statistics.

void RooHist::addEfficiencyBin(Axis_t binCenter, Int_t n1, Int_t n2, Double_t binWidth, Double_t xErrorFrac, Double_t scaleFactor)
{
  Double_t scale= 1;
  if(binWidth > 0) scale= _nominalBinWidth/binWidth;
  Int_t index= GetN();

  Double_t a= (Double_t)(n1)/(n1+n2);

  // calculate Binomial errors for this bin
  Double_t ym,yp,dx(0.5*binWidth);
  if(!RooHistError::instance().getBinomialIntervalEff(n1,n2,ym,yp,_nSigma)) {
    coutE(Plotting) << "RooHist::addEfficiencyBin: unable to calculate binomial error for bin with " << n1 << "," << n2 << " events" << endl;
    return;
  }

  SetPoint(index,binCenter,a*scaleFactor);
  SetPointError(index,dx*xErrorFrac,dx*xErrorFrac,(a-ym)*scaleFactor,(yp-a)*scaleFactor);
  updateYAxisLimits(scale*yp);
  updateYAxisLimits(scale*ym);
}



////////////////////////////////////////////////////////////////////////////////
/// Add a bin to this histogram with the value n1/(n1+n2)
/// using an error bar calculated with Binomial statistics.

void RooHist::addEfficiencyBinWithError(Axis_t binCenter, Double_t n1, Double_t n2, Double_t en1, Double_t en2, Double_t binWidth, Double_t xErrorFrac, Double_t scaleFactor)
{
  Double_t scale= 1;
  if(binWidth > 0) scale= _nominalBinWidth/binWidth;
  Int_t index= GetN();

  Double_t a= (Double_t)(n1)/(n1+n2);

  Double_t error = sqrt( pow(en1,2)*pow(n2,2) + pow(en2,2)*pow(n1,2) ) / pow(n1+n2,2) ;

  // calculate Binomial errors for this bin
  Double_t ym,yp,dx(0.5*binWidth);
  ym=a-error ;
  yp=a+error ;


  SetPoint(index,binCenter,a*scaleFactor);
  SetPointError(index,dx*xErrorFrac,dx*xErrorFrac,(a-ym)*scaleFactor,(yp-a)*scaleFactor);
  updateYAxisLimits(scale*yp);
  updateYAxisLimits(scale*ym);
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooHist::~RooHist()
{
}



////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if binning of this RooHist is identical to that of 'other'

Bool_t RooHist::hasIdenticalBinning(const RooHist& other) const
{
  // First check if number of bins is the same
  if (GetN() != other.GetN()) {
    return kFALSE ;
  }

  // Next require that all bin centers are the same
  Int_t i ;
  for (i=0 ; i<GetN() ; i++) {
    Double_t x1,x2,y1,y2 ;

    GetPoint(i,x1,y1) ;
    other.GetPoint(i,x2,y2) ;

    if (fabs(x1-x2)>1e-10) {
      return kFALSE ;
    }

  }

  return kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if contents of this RooHist is identical within given
/// relative tolerance to that of 'other'

Bool_t RooHist::isIdentical(const RooHist& other, Double_t tol) const
{
  // Make temporary TH1s output of RooHists to perform Kolmogorov test
  TH1::AddDirectory(kFALSE) ;
  TH1F h_self("h_self","h_self",GetN(),0,1) ;
  TH1F h_other("h_other","h_other",GetN(),0,1) ;
  TH1::AddDirectory(kTRUE) ;

  for (Int_t i=0 ; i<GetN() ; i++) {
    h_self.SetBinContent(i+1,GetY()[i]) ;
    h_other.SetBinContent(i+1,other.GetY()[i]) ;
  }

  Double_t M = h_self.KolmogorovTest(&h_other,"M") ;
  if (M>tol) {
    Double_t kprob = h_self.KolmogorovTest(&h_other) ;
    cout << "RooHist::isIdentical() tolerance exceeded M=" << M << " (tol=" << tol << "), corresponding prob = " << kprob << endl ;
    return kFALSE ;
  }

  return kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print info about this histogram to the specified output stream.
///
///   Standard: number of entries
///      Shape: error CL and maximum value
///    Verbose: print our bin contents and errors

void RooHist::printMultiline(ostream& os, Int_t contents, Bool_t verbose, TString indent) const
{
  RooPlotable::printMultiline(os,contents,verbose,indent);
  os << indent << "--- RooHist ---" << endl;
  Int_t n= GetN();
  os << indent << "  Contains " << n << " bins" << endl;
  if(verbose) {
    os << indent << "  Errors calculated at" << _nSigma << "-sigma CL" << endl;
    os << indent << "  Bin Contents:" << endl;
    for(Int_t i= 0; i < n; i++) {
      os << indent << setw(3) << i << ") x= " <<  fX[i];
      if(fEXhigh[i] > 0 || fEXlow[i] > 0) {
   os << " +" << fEXhigh[i] << " -" << fEXlow[i];
      }
      os << " , y = " << fY[i] << " +" << fEYhigh[i] << " -" << fEYlow[i] << endl;
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Print name of RooHist

void RooHist::printName(ostream& os) const
{
  os << GetName() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print title of RooHist

void RooHist::printTitle(ostream& os) const
{
  os << GetTitle() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print class name of RooHist

void RooHist::printClassName(ostream& os) const
{
  os << IsA()->GetName() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Create and return RooHist containing  residuals w.r.t to given curve.
/// If normalize is true, the residuals are normalized by the histogram
/// errors creating a RooHist with pull values

RooHist* RooHist::makeResidHist(const RooCurve& curve, bool normalize, bool useAverage) const
{

  // Copy all non-content properties from hist1
  RooHist* hist = new RooHist(_nominalBinWidth) ;
  if (normalize) {
    hist->SetName(Form("pull_%s_%s",GetName(),curve.GetName())) ;
    hist->SetTitle(Form("Pull of %s and %s",GetTitle(),curve.GetTitle())) ;
  } else {
    hist->SetName(Form("resid_%s_%s",GetName(),curve.GetName())) ;
    hist->SetTitle(Form("Residual of %s and %s",GetTitle(),curve.GetTitle())) ;
  }

  // Determine range of curve
  Double_t xstart,xstop,y ;
  curve.GetPoint(0,xstart,y) ;
  curve.GetPoint(curve.GetN()-1,xstop,y) ;

  // Add histograms, calculate Poisson confidence interval on sum value
  for(Int_t i=0 ; i<GetN() ; i++) {
    Double_t x,point;
    GetPoint(i,x,point) ;

    // Only calculate pull for bins inside curve range
    if (x<xstart || x>xstop) continue ;

    Double_t yy ;
    if (useAverage) {
      Double_t exl = GetErrorXlow(i);
      Double_t exh = GetErrorXhigh(i) ;
      if (exl<=0 ) exl = GetErrorX(i);
      if (exh<=0 ) exh = GetErrorX(i);
      if (exl<=0 ) exl = 0.5*getNominalBinWidth();
      if (exh<=0 ) exh = 0.5*getNominalBinWidth();
      yy = point - curve.average(x-exl,x+exh) ;
    } else {
      yy = point - curve.interpolate(x) ;
    }

    Double_t dyl = GetErrorYlow(i) ;
    Double_t dyh = GetErrorYhigh(i) ;
    if (normalize) {
        Double_t norm = (yy>0?dyl:dyh);
   if (norm==0.) {
     coutW(Plotting) << "RooHist::makeResisHist(" << GetName() << ") WARNING: point " << i << " has zero error, setting residual to zero" << endl ;
     yy=0 ;
     dyh=0 ;
     dyl=0 ;
   } else {
     yy   /= norm;
     dyh /= norm;
     dyl /= norm;
   }
    }
    hist->addBinWithError(x,yy,dyl,dyh);
  }
  return hist ;
}
