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
\file RooCurve.cxx
\class RooCurve
\ingroup Roofitcore

A RooCurve is a one-dimensional graphical representation of a real-valued function.
A curve is approximated by straight line segments with end points chosen to give
a "good" approximation to the true curve. The goodness of the approximation is
controlled by a precision and a resolution parameter.

A RooCurve derives from TGraph, so it can either be drawn as a line (default) or
as points:
```
RooPlot *p = y.plotOn(x.frame());
p->getAttMarker("curve_y")->SetMarkerStyle(20);
p->setDrawOptions("curve_y","PL");
p->Draw();
```

To retrieve a RooCurve from a RooPlot, use RooPlot::getCurve().
**/

#include "RooFit.h"

#include "RooCurve.h"
#include "RooHist.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooRealIntegral.h"
#include "RooRealBinding.h"
#include "RooScaledFunc.h"
#include "RooMsgService.h"

#include "Riostream.h"
#include "TClass.h"
#include "TMath.h"
#include "TAxis.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include <iomanip>
#include <deque>
#include <algorithm>

using namespace std ;

ClassImp(RooCurve);


////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

RooCurve::RooCurve() : _showProgress(kFALSE)
{
  initialize();
}


////////////////////////////////////////////////////////////////////////////////
/// Create a 1-dim curve of the value of the specified real-valued expression
/// as a function of x. Use the optional precision parameter to control
/// how precisely the smooth curve is rasterized. Use the optional argument set
/// to specify how the expression should be normalized. Use the optional scale
/// factor to rescale the expression after normalization.
/// If shiftToZero is set, the entire curve is shifted down to make the lowest
/// point of the curve go through zero.
RooCurve::RooCurve(const RooAbsReal &f, RooAbsRealLValue &x, Double_t xlo, Double_t xhi, Int_t xbins,
    Double_t scaleFactor, const RooArgSet *normVars, Double_t prec, Double_t resolution,
    Bool_t shiftToZero, WingMode wmode, Int_t nEvalError, Int_t doEEVal, Double_t eeVal,
    Bool_t showProg) :
        TGraph(),
        RooPlotable(),
        _showProgress(showProg)
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

  RooAbsFunc *funcPtr = 0;
  RooAbsFunc *rawPtr  = 0;
  funcPtr= f.bindVars(x,normVars,kTRUE);

  // apply a scale factor if necessary
  if(scaleFactor != 1) {
    rawPtr= funcPtr;
    funcPtr= new RooScaledFunc(*rawPtr,scaleFactor);
  }
  assert(0 != funcPtr);

  // calculate the points to add to our curve
  Double_t prevYMax = getYAxisMax() ;
  if(xbins > 0){
    // regular mode - use the sampling hint to decide where to evaluate the pdf
    list<Double_t>* hint = f.plotSamplingHint(x,xlo,xhi) ;
    addPoints(*funcPtr,xlo,xhi,xbins+1,prec,resolution,wmode,nEvalError,doEEVal,eeVal,hint);
    if (_showProgress) {
      ccoutP(Plotting) << endl ;
    }
    if (hint) {
      delete hint ;
    }
  } else {
    // if number of bins is set to <= 0, skip any interpolation and just evaluate the pdf at the bin centers
    // this is useful when plotting a pdf like a histogram
    int nBinsX = x.numBins();
    for(int i=0; i<nBinsX; ++i){
      double xval = x.getBinning().binCenter(i);
      addPoint(xval,(*funcPtr)(&xval)) ;      
    }
  }
  initialize();

  // cleanup
  delete funcPtr;
  if(rawPtr) delete rawPtr;
  if (shiftToZero) shiftCurveToZero(prevYMax) ;

  // Adjust limits
  for (int i=0 ; i<GetN() ; i++) {
    updateYAxisLimits(fY[i]);
  }
  this->Sort();
}



////////////////////////////////////////////////////////////////////////////////
/// Create a 1-dim curve of the value of the specified real-valued
/// expression as a function of x. Use the optional precision
/// parameter to control how precisely the smooth curve is
/// rasterized. If shiftToZero is set, the entire curve is shifted
/// down to make the lowest point in of the curve go through zero.

RooCurve::RooCurve(const char *name, const char *title, const RooAbsFunc &func,
		   Double_t xlo, Double_t xhi, UInt_t minPoints, Double_t prec, Double_t resolution,
		   Bool_t shiftToZero, WingMode wmode, Int_t nEvalError, Int_t doEEVal, Double_t eeVal) :
  _showProgress(kFALSE)
{
  SetName(name);
  SetTitle(title);
  Double_t prevYMax = getYAxisMax() ;
  addPoints(func,xlo,xhi,minPoints+1,prec,resolution,wmode,nEvalError,doEEVal,eeVal);  
  initialize();
  if (shiftToZero) shiftCurveToZero(prevYMax) ;

  // Adjust limits
  for (int i=0 ; i<GetN() ; i++) {
    updateYAxisLimits(fY[i]);
  }
  this->Sort();
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor of a curve as sum of two other curves.
/// \f[
/// C_\mathrm{sum} = \mathrm{scale1}*c1 + \mathrm{scale2}*c2
/// \f]
///
/// \param[in] name Name of the curve (to retrieve it from a plot)
/// \param[in] title Title (for plotting).
/// \param[in] c1 First curve.
/// \param[in] c2 Second curve.
/// \param[in] scale1 Scale y values for c1 by this factor.
/// \param[in] scale1 Scale y values for c2 by this factor.

RooCurve::RooCurve(const char* name, const char* title, const RooCurve& c1, const RooCurve& c2, Double_t scale1, Double_t scale2) :
  _showProgress(kFALSE)
{
  initialize() ;
  SetName(name) ;
  SetTitle(title) ;

  // Make deque of points in X
  deque<Double_t> pointList ;
  Double_t x,y ;

  // Add X points of C1
  Int_t i1,n1 = c1.GetN() ;
  for (i1=0 ; i1<n1 ; i1++) {
    c1.GetPoint(i1,x,y) ;
    pointList.push_back(x) ;
  }

  // Add X points of C2
  Int_t i2,n2 = c2.GetN() ;
  for (i2=0 ; i2<n2 ; i2++) {
    c2.GetPoint(i2,x,y) ;
    pointList.push_back(x) ;
  }
  
  // Sort X points
  sort(pointList.begin(),pointList.end()) ;

  // Loop over X points
  Double_t last(-RooNumber::infinity()) ;
  for (auto point : pointList) {

    if ((point-last)>1e-10) {
      // Add OR of points to new curve, skipping duplicate points within tolerance
      addPoint(point,scale1*c1.interpolate(point)+scale2*c2.interpolate(point)) ;
    }
    last = point ;
  }

  this->Sort();
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooCurve::~RooCurve() 
{
}



////////////////////////////////////////////////////////////////////////////////
/// Perform initialization that is common to all curves

void RooCurve::initialize() 
{
  // set default line width in pixels
  SetLineWidth(3);
  // set default line color
  SetLineColor(kBlue);
}



////////////////////////////////////////////////////////////////////////////////
/// Find lowest point in curve and move all points in curve so that
/// lowest point will go exactly through zero

void RooCurve::shiftCurveToZero(Double_t prevYMax) 
{
  Int_t i ;
  Double_t minVal(1e30) ;
  Double_t maxVal(-1e30) ;

  // First iteration, find current lowest point
  for (i=1 ; i<GetN()-1 ; i++) {
    Double_t x,y ;
    GetPoint(i,x,y) ;
    if (y<minVal) minVal=y ;
    if (y>maxVal) maxVal=y ;
  }

  // Second iteration, lower all points by minVal
  for (i=1 ; i<GetN()-1 ; i++) {
    Double_t x,y ;
    GetPoint(i,x,y) ;
    SetPoint(i,x,y-minVal) ;
  }

  // Check if y-axis range needs readjustment
  if (getYAxisMax()>prevYMax) {
    Double_t newMax = maxVal - minVal ;
    setYAxisLimits(getYAxisMin(), newMax<prevYMax ? prevYMax : newMax) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Add points calculated with the specified function, over the range (xlo,xhi).
/// Add at least minPoints equally spaced points, and add sufficient points so that
/// the maximum deviation from the final straight-line segments is prec*(ymax-ymin),
/// down to a minimum horizontal spacing of resolution*(xhi-xlo).

void RooCurve::addPoints(const RooAbsFunc &func, Double_t xlo, Double_t xhi,
			 Int_t minPoints, Double_t prec, Double_t resolution, WingMode wmode,
			 Int_t numee, Bool_t doEEVal, Double_t eeVal, list<Double_t>* samplingHint) 
{
  // check the inputs
  if(!func.isValid()) {
    coutE(InputArguments) << fName << "::addPoints: input function is not valid" << endl;
    return;
  }
  if(minPoints <= 0 || xhi <= xlo) {
    coutE(InputArguments) << fName << "::addPoints: bad input (nothing added)" << endl;
    return;
  }

  // Perform a coarse scan of the function to estimate its y range.
  // Save the results so we do not have to re-evaluate at the scan points.

  // Adjust minimum number of points to external sampling hint if used
  if (samplingHint) {
    minPoints = samplingHint->size() ;
  }

  Double_t dx= (xhi-xlo)/(minPoints-1.);
  std::vector<double> yval(minPoints);
  
  // Get list of initial x values. If function provides sampling hint use that,
  // otherwise use default binning of frame
  std::vector<double> xval;
  if (!samplingHint) {
    for(int step= 0; step < minPoints; step++) {
      xval.push_back(xlo + step*dx) ;
    }    
  } else {
    std::copy(samplingHint->begin(), samplingHint->end(), std::back_inserter(xval));
  }

  for (unsigned int step=0; step < xval.size(); ++step) {
    Double_t xx = xval[step];
    if (step == static_cast<unsigned int>(minPoints-1))
      xx -= 1e-15;

    yval[step]= func(&xx);
    if (_showProgress) {
      ccoutP(Plotting) << "." ;
      cout.flush() ;
    }

    if (RooAbsReal::numEvalErrors()>0) {
      if (numee>=0) {
        coutW(Plotting) << "At observable [x]=" << xx <<  " " ;
        RooAbsReal::printEvalErrors(ccoutW(Plotting),numee) ;
      }
      if (doEEVal) {
        yval[step]=eeVal ;
      }
    }
    RooAbsReal::clearEvalErrorLog() ;
  }

  const double ymax = *std::max_element(yval.begin(), yval.end());
  const double ymin = *std::min_element(yval.begin(), yval.end());
  Double_t yrangeEst=(ymax-ymin) ;

  // store points of the coarse scan and calculate any refinements necessary
  Double_t minDx= resolution*(xhi-xlo);
  Double_t x1,x2= xlo;

  if (wmode==Extended) {
    // Add two points to make curve jump from 0 to yval at the left end of the plotting range.
    // This ensures that filled polygons are drawn properly. The first point needs to be to the
    // left of the second. Since points are sorted later, its x coordinate is shifted by 1/1000 dx.
    addPoint(xlo-dx*1.001, 0);
    addPoint(xlo-dx,yval[0]) ;
  } else if (wmode==Straight) {
    addPoint(xlo,0) ;
  }

  addPoint(xlo,yval[0]);

  auto iter2 = xval.begin() ;
  x1 = *iter2 ;
  int step=1 ;
  while(true) {
    x1= x2;
    ++iter2 ;
    if (iter2==xval.end()) {
      break ;
    }
    x2= *iter2 ;
    if (prec<0) {
      // If precision is <0, no attempt at recursive interpolation is made
      addPoint(x2,yval[step]) ;
    } else {
      addRange(func,x1,x2,yval[step-1],yval[step],prec*yrangeEst,minDx,numee,doEEVal,eeVal);
    }
    step++ ;
  }
  addPoint(xhi,yval[minPoints-1]) ;

  if (wmode==Extended) {
    // Add two points to close polygon. The order matters. Since they are sorted in x later, the second
    // point is shifted by 1/1000 * dx.
    addPoint(xhi+dx,yval[minPoints-1]) ;
    addPoint(xhi+dx*1.001, 0);
  } else if (wmode==Straight) {
    addPoint(xhi,0) ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Fill the range (x1,x2) with points calculated using func(&x). No point will
/// be added at x1, and a point will always be added at x2. The density of points
/// will be calculated so that the maximum deviation from a straight line
/// approximation is prec*(ymax-ymin) down to the specified minimum horizontal spacing.

void RooCurve::addRange(const RooAbsFunc& func, Double_t x1, Double_t x2,
			Double_t y1, Double_t y2, Double_t minDy, Double_t minDx,
			Int_t numee, Bool_t doEEVal, Double_t eeVal) 
{
  // Explicitly skip empty ranges to eliminate point duplication
  if (fabs(x2-x1)<1e-20) {
    return ;
  }

  // calculate our value at the midpoint of this range
  Double_t xmid= 0.5*(x1+x2);
  Double_t ymid= func(&xmid);
  if (_showProgress) {
    ccoutP(Plotting) << "." ;
    cout.flush() ;
  }

  if (RooAbsReal::numEvalErrors()>0) {
    if (numee>=0) {
      coutW(Plotting) << "At observable [x]=" << xmid <<  " " ;
      RooAbsReal::printEvalErrors(ccoutW(Plotting),numee) ;
    }
    if (doEEVal) {
      ymid=eeVal ;
    }
  }
  RooAbsReal::clearEvalErrorLog() ;

  // test if the midpoint is sufficiently close to a straight line across this interval
  Double_t dy= ymid - 0.5*(y1+y2);
  if((xmid - x1 >= minDx) && fabs(dy)>0 && fabs(dy) >= minDy) {
    // fill in each subrange
    addRange(func,x1,xmid,y1,ymid,minDy,minDx,numee,doEEVal,eeVal);
    addRange(func,xmid,x2,ymid,y2,minDy,minDx,numee,doEEVal,eeVal);
  }
  else {
    // add the endpoint
    addPoint(x2,y2);
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Add a point with the specified coordinates. Update our y-axis limits.

void RooCurve::addPoint(Double_t x, Double_t y) 
{
//   cout << "RooCurve("<< GetName() << ") adding point at (" << x << "," << y << ")" << endl ;
  Int_t next= GetN();
  SetPoint(next, x, y);
  updateYAxisLimits(y) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the number of events associated with the plotable object,
/// it is always 1 for curves

Double_t RooCurve::getFitRangeNEvt() const {
  return 1;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the number of events associated with the plotable object,
/// in the given range. It is always 1 for curves

Double_t RooCurve::getFitRangeNEvt(Double_t, Double_t) const 
{
  return 1 ;
}


////////////////////////////////////////////////////////////////////////////////
/// Get the bin width associated with this plotable object.
/// It is alwats zero for curves

Double_t RooCurve::getFitRangeBinW() const {
  return 0 ;
}



////////////////////////////////////////////////////////////////////////////////

void RooCurve::printName(ostream& os) const 
// 
{
  // Print the name of this curve
  os << GetName() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Print the title of this curve

void RooCurve::printTitle(ostream& os) const 
{
  os << GetTitle() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Print the class name of this curve

void RooCurve::printClassName(ostream& os) const 
{
  os << IsA()->GetName() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print the details of this curve

void RooCurve::printMultiline(ostream& os, Int_t /*contents*/, Bool_t /*verbose*/, TString indent) const
{
  os << indent << "--- RooCurve ---" << endl ;
  Int_t n= GetN();
  os << indent << "  Contains " << n << " points" << endl;
  os << indent << "  Graph points:" << endl;
  for(Int_t i= 0; i < n; i++) {
    os << indent << setw(3) << i << ") x = " << fX[i] << " , y = " << fY[i] << endl;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate the chi^2/NDOF of this curve with respect to the histogram
/// 'hist' accounting nFitParam floating parameters in case the curve
/// was the result of a fit

Double_t RooCurve::chiSquare(const RooHist& hist, Int_t nFitParam) const 
{
  Int_t i,np = hist.GetN() ;
  Double_t x,y,eyl,eyh,exl,exh ;

  // Find starting and ending bin of histogram based on range of RooCurve
  Double_t xstart,xstop ;

  GetPoint(0,xstart,y) ;
  GetPoint(GetN()-1,xstop,y) ;

  Int_t nbin(0) ;

  Double_t chisq(0) ;
  for (i=0 ; i<np ; i++) {   

    // Retrieve histogram contents
    hist.GetPoint(i,x,y) ;

    // Check if point is in range of curve
    if (x<xstart || x>xstop) continue ;

    eyl = hist.GetEYlow()[i] ;
    eyh = hist.GetEYhigh()[i] ;
    exl = hist.GetEXlow()[i] ;
    exh = hist.GetEXhigh()[i] ;

    // Integrate function over this bin
    Double_t avg = average(x-exl,x+exh) ;

    // Add pull^2 to chisq
    if (y!=0) {      
      Double_t pull = (y>avg) ? ((y-avg)/eyl) : ((y-avg)/eyh) ;
      chisq += pull*pull ;
      nbin++ ;
    }
  }

  // Return chisq/nDOF 
  return chisq / (nbin-nFitParam) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return average curve value in [xFirst,xLast] by integrating curve between points
/// and dividing by xLast-xFirst

Double_t RooCurve::average(Double_t xFirst, Double_t xLast) const
{
  if (xFirst>=xLast) {
    coutE(InputArguments) << "RooCurve::average(" << GetName() 
			  << ") invalid range (" << xFirst << "," << xLast << ")" << endl ;
    return 0 ;
  }

  // Find Y values and begin and end points
  Double_t yFirst = interpolate(xFirst,1e-10) ;
  Double_t yLast = interpolate(xLast,1e-10) ;

  // Find first and last mid points
  Int_t ifirst = findPoint(xFirst,1e10) ;
  Int_t ilast  = findPoint(xLast,1e10) ;
  Double_t xFirstPt,yFirstPt,xLastPt,yLastPt ;
  GetPoint(ifirst,xFirstPt,yFirstPt) ;
  GetPoint(ilast,xLastPt,yLastPt) ;

  Double_t tolerance=1e-3*(xLast-xFirst) ;

  // Handle trivial scenario -- no midway points, point only at or outside given range
  if (ilast-ifirst==1 &&(xFirstPt-xFirst)<-1*tolerance && (xLastPt-xLast)>tolerance) {
    return 0.5*(yFirst+yLast) ;
  }
 
  // If first point closest to xFirst is at xFirst or before xFirst take the next point
  // as the first midway point   
  if ((xFirstPt-xFirst)<-1*tolerance) {
    ifirst++ ;
    const_cast<RooCurve&>(*this).GetPoint(ifirst,xFirstPt,yFirstPt) ;
  }
  
  // If last point closest to yLast is at yLast or beyond yLast the the previous point
  // as the last midway point
  if ((xLastPt-xLast)>tolerance) {
    ilast-- ;
    const_cast<RooCurve&>(*this).GetPoint(ilast,xLastPt,yLastPt) ;
  }

  Double_t sum(0),x1,y1,x2,y2 ;

  // Trapezoid integration from lower edge to first midpoint
  sum += (xFirstPt-xFirst)*(yFirst+yFirstPt)/2 ;

  // Trapezoid integration between midpoints
  Int_t i ;
  for (i=ifirst ; i<ilast ; i++) {
    const_cast<RooCurve&>(*this).GetPoint(i,x1,y1) ;
    const_cast<RooCurve&>(*this).GetPoint(i+1,x2,y2) ;
    sum += (x2-x1)*(y1+y2)/2 ;
  }

  // Trapezoid integration from last midpoint to upper edge 
  sum += (xLast-xLastPt)*(yLastPt+yLast)/2 ;
  return sum/(xLast-xFirst) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Find the nearest point to xvalue. Return -1 if distance
/// exceeds tolerance

Int_t RooCurve::findPoint(Double_t xvalue, Double_t tolerance) const
{
  Double_t delta(std::numeric_limits<double>::max()),x,y ;
  Int_t i,n = GetN() ;
  Int_t ibest(-1) ;
  for (i=0 ; i<n ; i++) {
    GetPoint(i,x,y);
    if (fabs(xvalue-x)<delta) {
      delta = fabs(xvalue-x) ;
      ibest = i ;
    }
  }

  return (delta<tolerance)?ibest:-1 ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return linearly interpolated value of curve at xvalue. If distance
/// to nearest point is less than tolerance, return nearest point value
/// instead

Double_t RooCurve::interpolate(Double_t xvalue, Double_t tolerance) const
{
  // Find best point
  int n = GetN() ;
  int ibest = findPoint(xvalue,1e10) ;
  
  // Get position of best point
  Double_t xbest, ybest ;
  const_cast<RooCurve*>(this)->GetPoint(ibest,xbest,ybest) ;

  // Handle trivial case of being dead on
  if (fabs(xbest-xvalue)<tolerance) {
    return ybest ;
  }

  // Get nearest point on other side w.r.t. xvalue
  Double_t xother,yother, retVal(0) ;
  if (xbest<xvalue) {
    if (ibest==n-1) {
      // Value beyond end requested -- return value of last point
      return ybest ;
    }
    const_cast<RooCurve*>(this)->GetPoint(ibest+1,xother,yother) ;        
    if (xother==xbest) return ybest ;
    retVal = ybest + (yother-ybest)*(xvalue-xbest)/(xother-xbest) ; 

  } else {
    if (ibest==0) {
      // Value before 1st point requested -- return value of 1st point
      return ybest ;
    }
    const_cast<RooCurve*>(this)->GetPoint(ibest-1,xother,yother) ;    
    if (xother==xbest) return ybest ;
    retVal = yother + (ybest-yother)*(xvalue-xother)/(xbest-xother) ;
  }
 
  return retVal ;
}




////////////////////////////////////////////////////////////////////////////////
/// Construct filled RooCurve represented error band that captures alpha% of the variations
/// of the curves passed through argument variations, where the percentage alpha corresponds to
/// the central interval fraction of a significance Z

RooCurve* RooCurve::makeErrorBand(const vector<RooCurve*>& variations, Double_t Z) const
{
  RooCurve* band = new RooCurve ;
  band->SetName(Form("%s_errorband",GetName())) ;
  band->SetLineWidth(1) ;
  band->SetFillColor(kCyan) ;
  band->SetLineColor(kCyan) ;

  vector<double> bandLo(GetN()) ;
  vector<double> bandHi(GetN()) ;
  for (int i=0 ; i<GetN() ; i++) {
    calcBandInterval(variations,i,Z,bandLo[i],bandHi[i],kFALSE) ;
  }
  
  for (int i=0 ; i<GetN() ; i++) {
    band->addPoint(GetX()[i],bandLo[i]) ;
  }
  for (int i=GetN()-1 ; i>=0 ; i--) {
    band->addPoint(GetX()[i],bandHi[i]) ;
  }	 
  // if the axis of the old graph is alphanumeric, copy the labels to the new one as well
  if(this->GetXaxis() && this->GetXaxis()->IsAlphanumeric()){
    band->GetXaxis()->Set(this->GetXaxis()->GetNbins(),this->GetXaxis()->GetXmin(),this->GetXaxis()->GetXmax());
    for(int i=0; i<this->GetXaxis()->GetNbins(); ++i){
      band->GetXaxis()->SetBinLabel(i+1,this->GetXaxis()->GetBinLabel(i+1));
    }
  }
  
  return band ;
}




////////////////////////////////////////////////////////////////////////////////
/// Construct filled RooCurve represented error band represent the error added in quadrature defined by the curves arguments
/// plusVar and minusVar corresponding to one-sigma variations of each parameter. The resulting error band, combined used the correlation matrix C
/// is multiplied with the significance parameter Z to construct the equivalent of a Z sigma error band (in Gaussian approximation)

RooCurve* RooCurve::makeErrorBand(const vector<RooCurve*>& plusVar, const vector<RooCurve*>& minusVar, const TMatrixD& C, Double_t Z) const
{
  
  RooCurve* band = new RooCurve ;
  band->SetName(Form("%s_errorband",GetName())) ;
  band->SetLineWidth(1) ;
  band->SetFillColor(kCyan) ;
  band->SetLineColor(kCyan) ;

  vector<double> bandLo(GetN()) ;
  vector<double> bandHi(GetN()) ;
  for (int i=0 ; i<GetN() ; i++) {
    calcBandInterval(plusVar,minusVar,i,C,Z,bandLo[i],bandHi[i]) ;
  }
  
  for (int i=0 ; i<GetN() ; i++) {
    band->addPoint(GetX()[i],bandLo[i]) ;
  }
  for (int i=GetN()-1 ; i>=0 ; i--) {
    band->addPoint(GetX()[i],bandHi[i]) ;
  }	   
  
  // if the axis of the old graph is alphanumeric, copy the labels to the new one as well
  if(this->GetXaxis() && this->GetXaxis()->IsAlphanumeric()){
    band->GetXaxis()->Set(this->GetXaxis()->GetNbins(),this->GetXaxis()->GetXmin(),this->GetXaxis()->GetXmax());
    for(int i=0; i<this->GetXaxis()->GetNbins(); ++i){
      band->GetXaxis()->SetBinLabel(i+1,this->GetXaxis()->GetBinLabel(i+1));
    }
  }

  return band ;
}





////////////////////////////////////////////////////////////////////////////////
/// Retrieve variation points from curves

void RooCurve::calcBandInterval(const vector<RooCurve*>& plusVar, const vector<RooCurve*>& minusVar,Int_t i, const TMatrixD& C, Double_t /*Z*/, Double_t& lo, Double_t& hi) const
{
  vector<double> y_plus(plusVar.size()), y_minus(minusVar.size()) ;
  Int_t j(0) ;
  for (vector<RooCurve*>::const_iterator iter=plusVar.begin() ; iter!=plusVar.end() ; ++iter) {
    y_plus[j++] = (*iter)->interpolate(GetX()[i]) ;    
  }
  j=0 ;
  for (vector<RooCurve*>::const_iterator iter=minusVar.begin() ; iter!=minusVar.end() ; ++iter) {
    y_minus[j++] = (*iter)->interpolate(GetX()[i]) ;
  }
  Double_t y_cen = GetY()[i] ;
  Int_t n = j ;

  // Make vector of variations
  TVectorD F(plusVar.size()) ;
  for (j=0 ; j<n ; j++) {
    F[j] = (y_plus[j]-y_minus[j])/2 ;
  }

  // Calculate error in linear approximation from variations and correlation coefficient
  Double_t sum = F*(C*F) ;

  lo= y_cen + sqrt(sum) ;
  hi= y_cen - sqrt(sum) ;
}



////////////////////////////////////////////////////////////////////////////////

void RooCurve::calcBandInterval(const vector<RooCurve*>& variations,Int_t i,Double_t Z, Double_t& lo, Double_t& hi, Bool_t approxGauss) const
{
  vector<double> y(variations.size()) ;
  Int_t j(0) ;
  for (vector<RooCurve*>::const_iterator iter=variations.begin() ; iter!=variations.end() ; ++iter) {
    y[j++] = (*iter)->interpolate(GetX()[i]) ;
}

  if (!approxGauss) {
    // Construct central 68% interval from variations collected at each point
    Double_t pvalue = TMath::Erfc(Z/sqrt(2.)) ;
    Int_t delta = Int_t( y.size()*(pvalue)/2 + 0.5) ;
    sort(y.begin(),y.end()) ;    
    lo = y[delta] ;
    hi = y[y.size()-delta] ;  
  } else {
    // Estimate R.M.S of variations at each point and use that as Gaussian sigma
    Double_t sum_y(0), sum_ysq(0) ;
    for (unsigned int k=0 ; k<y.size() ; k++) {
      sum_y   += y[k] ;
      sum_ysq += y[k]*y[k] ;
    }
    sum_y /= y.size() ;
    sum_ysq /= y.size() ;

    Double_t rms = sqrt(sum_ysq - (sum_y*sum_y)) ;
    lo = GetY()[i] - Z*rms ;
    hi = GetY()[i] + Z*rms ;    
  }
}




////////////////////////////////////////////////////////////////////////////////
/// Return true if curve is identical to other curve allowing for given
/// absolute tolerance on each point compared point.

Bool_t RooCurve::isIdentical(const RooCurve& other, Double_t tol) const 
{
  // Determine X range and Y range
  Int_t n= min(GetN(),other.GetN());
  Double_t xmin(1e30), xmax(-1e30), ymin(1e30), ymax(-1e30) ;
  for(Int_t i= 0; i < n; i++) {
    if (fX[i]<xmin) xmin=fX[i] ;
    if (fX[i]>xmax) xmax=fX[i] ;
    if (fY[i]<ymin) ymin=fY[i] ;
    if (fY[i]>ymax) ymax=fY[i] ;
  }
  const double Yrange=ymax-ymin ;

  Bool_t ret(kTRUE) ;
  for(Int_t i= 2; i < n-2; i++) {
    Double_t yTest = interpolate(other.fX[i],1e-10) ;
    Double_t rdy = fabs(yTest-other.fY[i])/Yrange ;
    if (rdy>tol) {
      ret = false;
      cout << "RooCurve::isIdentical[" << std::setw(3) << i << "] Y tolerance exceeded (" << std::setprecision(5) << std::setw(10) << rdy << ">" << tol << "),";
      cout << "  x,y=(" << std::right << std::setw(10) << fX[i] << "," << std::setw(10) << fY[i] << ")\tref: y="
          << std::setw(10) << other.interpolate(fX[i], 1.E-15) << ". [Nearest point from ref: ";
      auto j = other.findPoint(fX[i], 1.E10);
      std::cout << "j=" << j << "\tx,y=(" << std::setw(10) << other.fX[j] << "," << std::setw(10) << other.fY[j] << ") ]" << "\trange=" << Yrange << std::endl;
    }
  }
      
  return ret ;
}


