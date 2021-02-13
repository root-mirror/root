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
\file RooParamBinning.cxx
\class RooParamBinning
\ingroup Roofitcore

Class RooParamBinning is an implementation of RooAbsBinning that constructs
a binning with a range definition that depends on external RooAbsReal objects.
The external RooAbsReal definitions are explicitly allowed to depend on other
observables and parameters, and make it possible to define non-rectangular
range definitions in RooFit. Objects of class RooParamBinning are made
by the RooRealVar::setRange() that takes RooAbsReal references as arguments
**/

#include "RooFit.h"

#include "RooParamBinning.h"
#include "RooMsgService.h"

#include "Riostream.h"


using namespace std;

ClassImp(RooParamBinning);
;


////////////////////////////////////////////////////////////////////////////////
/// Default constructor
///   cout << "RooParamBinning(" << this << ") default ctor" << endl ;

RooParamBinning::RooParamBinning(const char* name) : 
  RooAbsBinning(name), _xlo(0), _xhi(0), _nbins(100), _binw(0), _lp(0), _owner(0)
{  
  _array = 0 ;
}


////////////////////////////////////////////////////////////////////////////////
/// Construct binning with 'nBins' bins and with a range
/// parameterized by external RooAbsReals xloIn and xhiIn.

RooParamBinning::RooParamBinning(RooAbsReal& xloIn, RooAbsReal& xhiIn, Int_t nBins, const char* name) :
  RooAbsBinning(name),
  _array(0), 
  _xlo(&xloIn),
  _xhi(&xhiIn),
  _nbins(nBins),
  _binw(0),
  _lp(0),
  _owner(0)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooParamBinning::~RooParamBinning() 
{
  if (_array) delete[] _array ;
  if (_lp) delete _lp ;
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor
///   cout << "RooParamBinning::cctor(" << this << ") orig = " << &other << endl ;

RooParamBinning::RooParamBinning(const RooParamBinning& other, const char* name) :
  RooAbsBinning(name), _binw(0), _owner(0)
{
  _array = 0 ;

  if (other._lp) {
//     cout << "RooParamBinning::cctor(this = " << this << ") taking addresses from orig  ListProxy" << endl ;
    _xlo = (RooAbsReal*) other._lp->at(0) ;
    _xhi = (RooAbsReal*) other._lp->at(1) ;

  } else {

//     cout << "RooParamBinning::cctor(this = " << this << ") taking addresses from orig pointers " << other._xlo << " " << other._xhi << endl ;

    _xlo   = other._xlo ;
    _xhi   = other._xhi ;
  }

  _nbins = other._nbins ;
  _lp = 0 ;

  //cout << "RooParamBinning::cctor(this = " << this << " xlo = " << &_xlo << " xhi = " << &_xhi << " _lp = " << _lp << " owner = " << _owner << ")" << endl ;
}



////////////////////////////////////////////////////////////////////////////////
/// Hook function called by RooAbsRealLValue when this binning
/// is inserted as binning for into given owner. Create
/// list proxy registered with owner that will track and implement
/// server directs to external RooAbsReals of this binning

void RooParamBinning::insertHook(RooAbsRealLValue& owner) const  
{
  _owner = &owner ;

  // If list proxy already exists update pointers from proxy
//   cout << "RooParamBinning::insertHook(" << this << "," << GetName() << ") _lp at beginning = " << _lp << endl ;
  if (_lp) {
//     cout << "updating raw pointers from list proxy contents" << endl ;
    _xlo = xlo() ;
    _xhi = xhi() ;
    delete _lp ;
  }
//   cout << "_xlo = " << _xlo << " _xhi = " << _xhi << endl ;

  // If list proxy does not exist, create it now
  _lp = new RooListProxy(Form("range::%s",GetName()),"lp",&owner,kFALSE,kTRUE) ;
  _lp->add(*_xlo) ;
  _lp->add(*_xhi) ;
  _xlo = 0 ;
  _xhi = 0 ;


}


////////////////////////////////////////////////////////////////////////////////
/// Hook function called by RooAbsRealLValue when this binning
/// is removed as binning for into given owner. Delete list
/// proxy that was inserted in owner

void RooParamBinning::removeHook(RooAbsRealLValue& /*owner*/) const  
{
  _owner = 0 ;
  
  // Remove list proxy from owner
  if (_lp) {
    _xlo = xlo() ;
    _xhi = xhi() ;
    delete _lp ;
    _lp = 0 ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Adjust range by adjusting values of external RooAbsReal values
/// Only functional when external representations are lvalues

void RooParamBinning::setRange(Double_t newxlo, Double_t newxhi) 
{
  if (newxlo>newxhi) {
    coutE(InputArguments) << "RooParamBinning::setRange: ERROR low bound > high bound" << endl ;
    return ;
  }

  RooAbsRealLValue* xlolv = dynamic_cast<RooAbsRealLValue*>(xlo()) ;
  if (xlolv) {
    xlolv->setVal(newxlo) ;
  } else {
    coutW(InputArguments) << "RooParamBinning::setRange: WARNING lower bound not represented by lvalue, cannot set lower bound value through setRange()" << endl ;
  }

  RooAbsRealLValue* xhilv = dynamic_cast<RooAbsRealLValue*>(xhi()) ;
  if (xhilv) {
    xhilv->setVal(newxhi) ;
  } else {
    coutW(InputArguments) << "RooParamBinning::setRange: WARNING upper bound not represented by lvalue, cannot set upper bound value through setRange()" << endl ;
  }

}



////////////////////////////////////////////////////////////////////////////////
/// Return the fit bin index for the current value

Int_t RooParamBinning::binNumber(Double_t x) const  
{
  if (x >= xhi()->getVal()) return _nbins-1 ;
  if (x < xlo()->getVal()) return 0 ;

  return Int_t((x - xlo()->getVal())/averageBinWidth()) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the central value of the 'i'-th fit bin

Double_t RooParamBinning::binCenter(Int_t i) const 
{
  if (i<0 || i>=_nbins) {
    coutE(InputArguments) << "RooParamBinning::binCenter ERROR: bin index " << i 
			  << " is out of range (0," << _nbins-1 << ")" << endl ;
    return 0 ;
  }

  return xlo()->getVal() + (i + 0.5)*averageBinWidth() ;  
}




////////////////////////////////////////////////////////////////////////////////
/// Return average bin width

Double_t RooParamBinning::binWidth(Int_t /*bin*/) const 
{
  return (xhi()->getVal()-xlo()->getVal())/_nbins ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the low edge of the 'i'-th fit bin

Double_t RooParamBinning::binLow(Int_t i) const 
{
  if (i<0 || i>=_nbins) {
    coutE(InputArguments) << "RooParamBinning::binLow ERROR: bin index " << i 
			  << " is out of range (0," << _nbins-1 << ")" << endl ;
    return 0 ;
  }

  return xlo()->getVal() + i*binWidth(i) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the high edge of the 'i'-th fit bin

Double_t RooParamBinning::binHigh(Int_t i) const 
{
  if (i<0 || i>=_nbins) {
    coutE(InputArguments) << "RooParamBinning::fitBinHigh ERROR: bin index " << i 
			  << " is out of range (0," << _nbins-1 << ")" << endl ;
    return 0 ;
  }

  return xlo()->getVal() + (i + 1)*binWidth(i) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return array of bin boundaries

Double_t* RooParamBinning::array() const 
{
  if (_array) delete[] _array ;
  _array = new Double_t[_nbins+1] ;

  Int_t i ;
  for (i=0 ; i<=_nbins ; i++) {
    _array[i] = xlo()->getVal() + i*binWidth(i) ;
  }
  return _array ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print details of binning

void RooParamBinning::printMultiline(ostream &os, Int_t /*content*/, Bool_t /*verbose*/, TString indent) const
{
  os << indent << "_xlo = " << _xlo << endl ;
  os << indent << "_xhi = " << _xhi << endl ;
  if (_lp) {
    os << indent << "xlo() = " << xlo() << endl ;
    os << indent << "xhi() = " << xhi() << endl ;
  }  
  if (xlo()) {
    xlo()->Print("t") ;
  }
  if (xhi()) {
    xhi()->Print("t") ;
  }
}


