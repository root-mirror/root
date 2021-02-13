/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *   AB, Adrian Bevan, Liverpool University, bevan@slac.stanford.edu         *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California,         *
 *                          Liverpool University,                            *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/** \class Roo2DKeysPdf
    \ingroup Roofit

Two-dimensional kernel estimation PDF.

<b>This function has been superseded by the more general RooNDKeysPdf.</b>
*/

#include "RooFit.h"

#include "Roo2DKeysPdf.h"
#include "RooRealVar.h"
#include "TTree.h"
#include "TH2.h"
#include "TFile.h"
#include "TBranch.h"
#include "TMath.h"

//#include <math.h>

using namespace std;

ClassImp(Roo2DKeysPdf);


////////////////////////////////////////////////////////////////////////////////
/// Constructor.
/// \param[in] name
/// \param[in] title
/// \param[in] xx
/// \param[in] yy
/// \param[in] data
/// \param[in] options
/// \param[in] widthScaleFactor

Roo2DKeysPdf::Roo2DKeysPdf(const char *name, const char *title,
                       RooAbsReal& xx, RooAbsReal & yy, RooDataSet& data,  TString options, Double_t widthScaleFactor):
  RooAbsPdf(name,title),
  x("x", "x dimension",this, xx),
  y("y", "y dimension",this, yy)
{
  setWidthScaleFactor(widthScaleFactor);
  loadDataSet(data, options);
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.
/// \param[in] other
/// \param[in] name

Roo2DKeysPdf::Roo2DKeysPdf(const Roo2DKeysPdf & other, const char* name) :
  RooAbsPdf(other,name),
  x("x", this, other.x),
  y("y", this, other.y)
{
  if(_verbosedebug) { cout << "Roo2DKeysPdf::Roo2DKeysPdf copy ctor" << endl; }

  _xMean   = other._xMean;
  _xSigma  = other._xSigma;
  _yMean   = other._yMean;
  _ySigma  = other._ySigma;
  _n       = other._n;

  _BandWidthType    = other._BandWidthType;
  _MirrorAtBoundary = other._MirrorAtBoundary;
  _widthScaleFactor = other._widthScaleFactor;

  _2pi     = other._2pi;
  _sqrt2pi = other._sqrt2pi;
  _nEvents = other._nEvents;
  _n16     = other._n16;
  _debug   = other._debug;
  _verbosedebug   = other._verbosedebug;
  _vverbosedebug  = other._vverbosedebug;

  _lox       = other._lox;
  _hix       = other._hix;
  _loy       = other._loy;
  _hiy       = other._hiy;
  _xoffset   = other._xoffset;
  _yoffset   = other._yoffset;

  _x  = new Double_t[_nEvents];
  _y  = new Double_t[_nEvents];
  _hx = new Double_t[_nEvents];
  _hy = new Double_t[_nEvents];

  //copy the data and bandwidths
  for(Int_t iEvt = 0; iEvt< _nEvents; iEvt++)
  {
    _x[iEvt]  = other._x[iEvt];
    _y[iEvt]  = other._y[iEvt];
    _hx[iEvt] = other._hx[iEvt];
    _hy[iEvt] = other._hy[iEvt];
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

Roo2DKeysPdf::~Roo2DKeysPdf() {
  if(_verbosedebug) { cout << "Roo2DKeysPdf::Roo2KeysPdf dtor" << endl; }
    delete[] _x;
    delete[] _hx;
    delete[] _y;
    delete[] _hy;
}


////////////////////////////////////////////////////////////////////////////////
/// Loads a new data set into the class instance.
/// Returns 1 in case of error, 0 otherwise.
/// \param[in] data
/// \param[in] options

Int_t Roo2DKeysPdf::loadDataSet(RooDataSet& data, TString options)
{
  if(_verbosedebug) { cout << "Roo2DKeysPdf::loadDataSet" << endl; }

  setOptions(options);

  if(_verbosedebug) { cout << "Roo2DKeysPdf::loadDataSet(RooDataSet& data, TString options)" << endl; }

  _2pi       = 2.0*TMath::Pi() ;   //use pi from math.h
  _sqrt2pi   = sqrt(_2pi);
  _nEvents   = (Int_t)data.numEntries();
  if(_nEvents == 0)
  {
    cout << "ERROR:  Roo2DKeysPdf::loadDataSet The input data set is empty.  Unable to begin generating the PDF" << endl;
    return 1;
  }
  _n16       =  TMath::Power(_nEvents, -0.166666666); // = (4/[n(dim(R) + 2)])^1/(dim(R)+4); dim(R) = 2

  _lox       = x.min();
  _hix       = x.max();
  _loy       = y.min();
  _hiy       = y.max();

  _x  = new Double_t[_nEvents];
  _y  = new Double_t[_nEvents];
  _hx = new Double_t[_nEvents];
  _hy = new Double_t[_nEvents];

  Double_t x0 = 0.0;
  Double_t x1 = 0.0;
  Double_t x_2 = 0.0;
  Double_t y0 = 0.0;
  Double_t y1 = 0.0;
  Double_t y_2 = 0.0;

  //check that the data contain the variable we are interested in
  Int_t bad = 0;
  const RooAbsReal & xx = x.arg();
  const RooAbsReal & yy = y.arg();
  if(! (RooRealVar*)( (RooArgSet*)data.get(0) )->find( xx.GetName() ) )
  {
    cout << "Roo2DKeysPdf::Roo2DKeysPdf invalid RooAbsReal name: "<<xx.GetName()<<" not in the data set" <<endl;
    bad = 1;
  }
  if(! (RooRealVar*)( (RooArgSet*)data.get(0) )->find( yy.GetName() ) )
  {
    cout << "Roo2DKeysPdf::Roo2DKeysPdf invalid RooAbsReal name: "<<yy.GetName()<<" not in the data set" << endl;
    bad = 1;
  }
  if(bad)
  {
    cout << "Roo2DKeysPdf::Roo2DKeysPdf Unable to initilize object; incompatible RooDataSet doesn't contain"<<endl;
    cout << "                           all of the RooAbsReal arguments"<<endl;
    return 1;
  }

  //copy the data into local arrays
  const RooArgSet * values = data.get();
  const RooRealVar* X = ((RooRealVar*)(values->find(xx.GetName())) ) ;
  const RooRealVar* Y = ((RooRealVar*)(values->find(yy.GetName())) ) ;

  for (Int_t j=0;j<_nEvents;++j)
  {
    data.get(j) ;

    _x[j] = X->getVal() ;
    _y[j] = Y->getVal() ;

    x0+=1; x1+=_x[j]; x_2+=_x[j]*_x[j];
    y0+=1; y1+=_y[j]; y_2+=_y[j]*_y[j];
  }

  //==========================================//
  //calculate the mean and sigma for the data //
  //==========================================//
  if(_nEvents == 0)
  {
    cout << "Roo2DKeysPdf::Roo2DKeysPdf Empty data set was used; can't generate a PDF"<<endl;
  }

  _xMean  = x1/x0;
  _xSigma = sqrt(x_2/_nEvents-_xMean*_xMean);

  _yMean  = y1/y0;
  _ySigma = sqrt(y_2/_nEvents-_yMean*_yMean);

  _n = Double_t(1)/(_2pi*_nEvents*_xSigma*_ySigma);

  //calculate the PDF
  return calculateBandWidth(_BandWidthType);
}


////////////////////////////////////////////////////////////////////////////////

void Roo2DKeysPdf::setOptions(TString options)
{
  if(_verbosedebug) { cout << "Roo2DKeysPdf::setOptions" << endl; }

  options.ToLower();
  if( options.Contains("a") )   _BandWidthType    = 0;
  else                          _BandWidthType    = 1;
  if( options.Contains("n") )   _BandWidthType    = 1;
  else                          _BandWidthType    = 0;
  if( options.Contains("m") )   _MirrorAtBoundary = 1;
  else                          _MirrorAtBoundary = 0;
  if( options.Contains("d") )   _debug            = 1;
  else                          _debug            = 0;
  if( options.Contains("v") )   { _debug         = 1; _verbosedebug = 1; }
  else                            _verbosedebug  = 0;
  if( options.Contains("vv") )  { _vverbosedebug = 1; }
  else                            _vverbosedebug = 0;

  if( _debug )
  {
    cout << "Roo2DKeysPdf::setOptions(TString options)    options = "<< options << endl;
    cout << "\t_BandWidthType    = " << _BandWidthType    << endl;
    cout << "\t_MirrorAtBoundary = " << _MirrorAtBoundary << endl;
    cout << "\t_debug            = " << _debug            << endl;
    cout << "\t_verbosedebug     = " << _verbosedebug     << endl;
    cout << "\t_vverbosedebug    = " << _vverbosedebug    << endl;
  }
}


////////////////////////////////////////////////////////////////////////////////

void Roo2DKeysPdf::getOptions(void) const
{
  cout << "Roo2DKeysPdf::getOptions(void)" << endl;
  cout << "\t_BandWidthType                           = " << _BandWidthType    << endl;
  cout << "\t_MirrorAtBoundary                        = " << _MirrorAtBoundary << endl;
  cout << "\t_debug                                   = " << _debug            << endl;
  cout << "\t_verbosedebug                            = " << _verbosedebug     << endl;
  cout << "\t_vverbosedebug                           = " << _vverbosedebug    << endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Calculates the kernel bandwidth for x & y and the probability look up table _p[i][j]
/// \param[in] kernel

Int_t Roo2DKeysPdf::calculateBandWidth(Int_t kernel)
{
  if(_verbosedebug) { cout << "Roo2DKeysPdf::calculateBandWidth(Int_t kernel)" << endl; }
  if(kernel != -999)
  {
    _BandWidthType = kernel;
  }

  Double_t h = 0.0;

  Double_t sigSum       = _xSigma*_xSigma + _ySigma*_ySigma;
  Double_t sqrtSum      = sqrt( sigSum );
  Double_t sigProd      = _ySigma*_xSigma;
  if(sigProd != 0.0)  h = _n16*sqrt( sigSum/sigProd );
  if(sqrtSum == 0)
  {
    cout << "Roo2DKeysPdf::calculateBandWidth The sqr(variance sum) == 0.0. " << " Your dataset represents a delta function."<<endl;
    return 1;
  }

  Double_t hXSigma = h * _xSigma;
  Double_t hYSigma = h * _ySigma;
  Double_t xhmin   = hXSigma * sqrt(2.)/10;  //smallest anticipated bandwidth
  Double_t yhmin   = hYSigma * sqrt(2.)/10;

  //////////////////////////////////////
  //calculate bandwidths from the data//
  //////////////////////////////////////
  if(_BandWidthType == 1)  //calculate a trivial bandwidth
  {
    cout << "Roo2DKeysPdf::calculateBandWidth Using a normal bandwidth (same for a given dimension) based on"<<endl;
    cout << "                                 h_j = n^{-1/6}*sigma_j for the j^th dimension and n events * "<<_widthScaleFactor<<endl;
    Double_t hxGaussian = _n16 * _xSigma * _widthScaleFactor;
    Double_t hyGaussian = _n16 * _ySigma * _widthScaleFactor;
    for(Int_t j=0;j<_nEvents;++j)
    {
      _hx[j] = hxGaussian;
      _hy[j] = hyGaussian;
      if(_hx[j]<xhmin) _hx[j] = xhmin;
      if(_hy[j]<yhmin) _hy[j] = yhmin;
     }
  }
  else //use an adaptive bandwidth to reduce the dependence on global data distribution
  {
    cout << "Roo2DKeysPdf::calculateBandWidth Using an adaptive bandwidth (in general different for all events) [default]"<<endl;
    cout << "                                 scaled by a factor of "<<_widthScaleFactor<<endl;
    Double_t xnorm   = h * TMath::Power(_xSigma/sqrtSum, 1.5) * _widthScaleFactor;
    Double_t ynorm   = h * TMath::Power(_ySigma/sqrtSum, 1.5) * _widthScaleFactor;
    for(Int_t j=0;j<_nEvents;++j)
    {
      Double_t f_ti =  TMath::Power( g(_x[j], _x, hXSigma, _y[j], _y, hYSigma), -0.25 ) ;
      _hx[j] = xnorm * f_ti;
      _hy[j] = ynorm * f_ti;
      if(_hx[j]<xhmin) _hx[j] = xhmin;
      if(_hy[j]<yhmin) _hy[j] = yhmin;
    }
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Evaluates the kernel estimation for x,y, interpolating between the points if necessary
///
/// Uses the caching intrinsic in RFC to bypass the grid and remove
/// the grid and extrapolation approximation in the kernel estimation method
/// implementation.

Double_t Roo2DKeysPdf::evaluate() const
{
  if(_vverbosedebug) { cout << "Roo2DKeysPdf::evaluate()" << endl; }
  return evaluateFull(x,y);
}


////////////////////////////////////////////////////////////////////////////////
/// Evaluates the sum of the product of the 2D kernels
/// for use in calculating the fixed kernel estimate, f,
/// given the bandwidths _hx[j] and _hy[j].
///
/// _n is calculated once in the constructor.
/// \param[in] thisX
/// \param[in] thisY

Double_t Roo2DKeysPdf::evaluateFull(Double_t thisX, Double_t thisY) const
{
  if( _vverbosedebug ) { cout << "Roo2DKeysPdf::evaluateFull()" << endl; }

  Double_t f=0.0;

  Double_t rx2, ry2, zx, zy;
  if( _MirrorAtBoundary )
  {
    for (Int_t j = 0; j < _nEvents; ++j)
    {
      rx2 = 0.0; ry2 = 0.0; zx = 0.0; zy = 0.0;
      if(_hx[j] != 0.0) rx2 = (thisX - _x[j])/_hx[j];
      if(_hy[j] != 0.0) ry2 = (thisY - _y[j])/_hy[j];

      if(_hx[j] != 0.0) zx = exp(-0.5*rx2*rx2)/_hx[j];
      if(_hy[j] != 0.0) zy = exp(-0.5*ry2*ry2)/_hy[j];

      zx += highBoundaryCorrection(thisX, _hx[j], x.max(), _x[j])
 	 +   lowBoundaryCorrection(thisX, _hx[j], x.min(), _x[j]);
      zy += highBoundaryCorrection(thisY, _hy[j], y.max(), _y[j])
 	 +   lowBoundaryCorrection(thisY, _hy[j], y.min(), _y[j]);
      f += zy * zx;
      //      f += _n * zy * zx; // ooops this is a normalisation factor :(
    }
  }
  else
  {
    for (Int_t j = 0; j < _nEvents; ++j)
    {
      rx2 = 0.0; ry2 = 0.0; zx = 0.0; zy = 0.0;
      if(_hx[j] != 0.0) rx2 = (thisX - _x[j])/_hx[j];
      if(_hy[j] != 0.0) ry2 = (thisY - _y[j])/_hy[j];

      if(_hx[j] != 0.0) zx = exp(-0.5*rx2*rx2)/_hx[j];
      if(_hy[j] != 0.0) zy = exp(-0.5*ry2*ry2)/_hy[j];
      f += zy * zx;
      //      f += _n * zy * zx; // ooops this is a normalisation factor  :(
    }
  }
  return f;
}


////////////////////////////////////////////////////////////////////////////////
/// Apply the mirror at boundary correction to a dimension given the space position to evaluate
/// at (thisVar), the bandwidth at this position (thisH), the boundary (high/low) and the
/// value of the data kernel that this correction is being applied to tVar (i.e. the _x[ix] etc.).
/// \param[in] thisVar
/// \param[in] thisH
/// \param[in] high
/// \param[in] tVar

Double_t Roo2DKeysPdf::highBoundaryCorrection(Double_t thisVar, Double_t thisH, Double_t high, Double_t tVar) const
{
  if(_vverbosedebug) { cout << "Roo2DKeysPdf::highBoundaryCorrection" << endl; }

  if(thisH == 0.0) return 0.0;
  Double_t correction = (thisVar + tVar - 2.0* high )/thisH;
  return exp(-0.5*correction*correction)/thisH;
}


////////////////////////////////////////////////////////////////////////////////

Double_t Roo2DKeysPdf::lowBoundaryCorrection(Double_t thisVar, Double_t thisH, Double_t low, Double_t tVar) const
{
  if(_vverbosedebug) { cout << "Roo2DKeysPdf::lowBoundaryCorrection" << endl; }

  if(thisH == 0.0) return 0.0;
  Double_t correction = (thisVar + tVar - 2.0* low )/thisH;
  return exp(-0.5*correction*correction)/thisH;
}


////////////////////////////////////////////////////////////////////////////////
/// Calculates f(t_i) for the bandwidths.
/// \f$ g = 1/(N_{evt} * \sigma_j * \sqrt{2\pi})*\sum_{all evts}{\prod d K[ \exp({-(xd - ti)/\sigma_{j}d^2}) ]}\f$
/// \param[in] varMean1
/// \param[in] _var1
/// \param[in] sigma1
/// \param[in] varMean2
/// \param[in] _var2
/// \param[in] sigma2

Double_t Roo2DKeysPdf::g(Double_t varMean1, Double_t * _var1, Double_t sigma1, Double_t varMean2, Double_t * _var2, Double_t sigma2) const
{
  if((_nEvents == 0.0) || (sigma1 == 0.0) || (sigma2 == 0)) return 0.0;

  Double_t c1 = -1.0/(2.0*sigma1*sigma1);
  Double_t c2 = -1.0/(2.0*sigma2*sigma2);
  Double_t d  = 4.0*c1*c2  /(_sqrt2pi*_nEvents);
  Double_t z  = 0.0;

  for (Int_t i = 0; i < _nEvents; ++i)
  {
    Double_t r1 =  _var1[i] - varMean1;
    Double_t r2 =  _var2[i] - varMean2;
    z          += exp( c1 * r1*r1 ) * exp( c2 * r2*r2 );
  }
  z = z*d;
  return z;
}


////////////////////////////////////////////////////////////////////////////////

Int_t Roo2DKeysPdf::getBandWidthType() const
{
  if(_BandWidthType == 1)  cout << "The Bandwidth Type selected is Trivial" << endl;
  else                     cout << "The Bandwidth Type selected is Adaptive" << endl;

  return _BandWidthType;
}


////////////////////////////////////////////////////////////////////////////////

Double_t Roo2DKeysPdf::getMean(const char * axis) const
{
  if(!strcmp(axis,x.GetName()) || !strcmp(axis,"x") || !strcmp(axis,"X"))      return _xMean;
  else if(!strcmp(axis,y.GetName()) || !strcmp(axis,"y") || !strcmp(axis,"Y")) return _yMean;
  else
  {
    cout << "Roo2DKeysPdf::getMean unknown axis "<<axis<<endl;
  }
  return 0.0;
}


////////////////////////////////////////////////////////////////////////////////

Double_t Roo2DKeysPdf::getSigma(const char * axis) const
{
  if(!strcmp(axis,x.GetName()) || !strcmp(axis,"x") || !strcmp(axis,"X"))      return _xSigma;
  else if(!strcmp(axis,y.GetName()) || !strcmp(axis,"y") || !strcmp(axis,"Y")) return _ySigma;
  else
  {
    cout << "Roo2DKeysPdf::getSigma unknown axis "<<axis<<endl;
  }
  return 0.0;
}



////////////////////////////////////////////////////////////////////////////////

void Roo2DKeysPdf::writeToFile(char * outputFile, const char * name) const
{
  TString histName = name;
  histName        += "_hist";
  TString nName    = name;
  nName           += "_Ntuple";
  writeHistToFile( outputFile,    histName);
  writeNTupleToFile( outputFile,  nName);
}


////////////////////////////////////////////////////////////////////////////////
/// Plots the PDF as a histogram and saves it to a file, so that it can be loaded in
/// as a Roo2DHist PDF in the future to save on calculation time.
/// \param[in] outputFile Name of the file where to store the PDF
/// \param[in] histName PDF histogram name

void Roo2DKeysPdf::writeHistToFile(char * outputFile, const char * histName) const
{
  TFile * file = 0;
  cout << "Roo2DKeysPdf::writeHistToFile This member function is temporarily disabled" <<endl;
  //make sure that any existing file is not over written
  file = new TFile(outputFile, "UPDATE");
  if (!file)
  {
    cout << "Roo2DKeysPdf::writeHistToFile unable to open file "<< outputFile <<endl;
    return;
  }


  const RooAbsReal & xx = x.arg();
  const RooAbsReal & yy = y.arg();
  RooArgSet values( RooArgList( xx, yy ));
  RooRealVar * xArg = ((RooRealVar*)(values.find(xx.GetName())) ) ;
  RooRealVar * yArg = ((RooRealVar*)(values.find(yy.GetName())) ) ;

  TH2F * hist = (TH2F*)xArg->createHistogram("hist", *yArg);
  hist = (TH2F*)this->fillHistogram(hist, RooArgList(*xArg, *yArg) );
  hist->SetName(histName);

  file->Write();
  file->Close();
}


////////////////////////////////////////////////////////////////////////////////
/// Saves the data and calculated bandwidths to a file,
/// as a record of what produced the PDF and to give a reduced
/// data set in order to facilitate re-calculation in the future.
/// \param[in] outputFile Name of the file where to store the data
/// \param[in] name Name of the tree which will contain the data

void Roo2DKeysPdf::writeNTupleToFile(char * outputFile, const char * name) const
{
  TFile * file = 0;

  //make sure that any existing file is not over written
  file = new TFile(outputFile, "UPDATE");
  if (!file)
  {
    cout << "Roo2DKeysPdf::writeNTupleToFile unable to open file "<< outputFile <<endl;
    return;
  }
  RooAbsReal & xArg = (RooAbsReal&)x.arg();
  RooAbsReal & yArg = (RooAbsReal&)y.arg();

  Double_t theX, theY, hx/*, hy*/;
  TString label = name;
  label += " the source data for 2D Keys PDF";
  TTree * _theTree =  new TTree(name, label);
  if(!_theTree) { cout << "Unable to get a TTree for output" << endl; return; }
  _theTree->SetAutoSave(1000000000);  // autosave when 1 Gbyte written

  //name the TBranches the same as the RooAbsReal's
  const char * xname = xArg.GetName();
  const char * yname = yArg.GetName();
  if (!strcmp(xname,"")) xname = "x";
  if (!strcmp(yname,"")) yname = "y";

  _theTree->Branch(xname, &theX, " x/D");
  _theTree->Branch(yname, &theY, " y/D");
  _theTree->Branch("hx",  &hx,   " hx/D");
  _theTree->Branch("hy",  &hx,  " hy/D");

  for(Int_t iEvt = 0; iEvt < _nEvents; iEvt++)
  {
    theX = _x[iEvt];
    theY = _y[iEvt];
    hx   = _hx[iEvt];
    hx   = _hy[iEvt];
    _theTree->Fill();
  }
  file->Write();
  file->Close();
}


////////////////////////////////////////////////////////////////////////////////
/// Prints out _p[_nPoints][_nPoints] indicating the domain limits.
/// \param[out] out Output stream where to print

void Roo2DKeysPdf::PrintInfo(ostream & out) const
{
  out << "Roo2DKeysPDF instance domain information:"<<endl;
  out << "\tX_min          = " << _lox <<endl;
  out << "\tX_max          = " << _hix <<endl;
  out << "\tY_min          = " << _loy <<endl;
  out << "\tY_max          = " << _hiy <<endl;

  out << "Data information:" << endl;
  out << "\t<x>             = " << _xMean <<endl;
  out << "\tsigma(x)       = " << _xSigma <<endl;
  out << "\t<y>             = " << _yMean <<endl;
  out << "\tsigma(y)       = " << _ySigma <<endl;

  out << "END of info for Roo2DKeys pdf instance"<< endl;
}
