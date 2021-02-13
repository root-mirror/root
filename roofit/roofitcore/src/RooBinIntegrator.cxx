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
\file RooBinIntegrator.cxx
\class RooBinIntegrator
\ingroup Roofitcore

RooBinIntegrator computes the integral over a binned distribution by summing the bin
contents of all bins.
**/


#include "RooFit.h"
#include "Riostream.h"

#include "TClass.h"
#include "RooBinIntegrator.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooNumber.h"
#include "RooIntegratorBinding.h"
#include "RooNumIntConfig.h"
#include "RooNumIntFactory.h"
#include "RooMsgService.h"

#include <assert.h>



using namespace std;

ClassImp(RooBinIntegrator);
;

// Register this class with RooNumIntConfig

////////////////////////////////////////////////////////////////////////////////
/// Register RooBinIntegrator, is parameters and capabilities with RooNumIntFactory

void RooBinIntegrator::registerIntegrator(RooNumIntFactory& fact)
{
  RooRealVar numBins("numBins","Number of bins in range",100) ;
  RooBinIntegrator* proto = new RooBinIntegrator() ;
  fact.storeProtoIntegrator(proto,RooArgSet(numBins)) ;
  RooNumIntConfig::defaultConfig().method1D().setLabel(proto->IsA()->GetName()) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooBinIntegrator::RooBinIntegrator() : _numBins(0), _useIntegrandLimits(kFALSE), _x(0)
{
}


////////////////////////////////////////////////////////////////////////////////
/// Construct integrator on given function binding binding

RooBinIntegrator::RooBinIntegrator(const RooAbsFunc& function) : 
  RooAbsIntegrator(function)
{
  _useIntegrandLimits= kTRUE;
  assert(0 != integrand() && integrand()->isValid());

  // Allocate coordinate buffer size after number of function dimensions
  _x = new Double_t[_function->getDimension()] ;
  _numBins = 100 ;

  _xmin.resize(_function->getDimension()) ;
  _xmax.resize(_function->getDimension()) ;

  for (UInt_t i=0 ; i<_function->getDimension() ; i++) {
    _xmin[i]= integrand()->getMinLimit(i);
    _xmax[i]= integrand()->getMaxLimit(i);

    // Retrieve bin configuration from integrand
    list<Double_t>* tmp = integrand()->binBoundaries(i) ;
    if (!tmp) {
      oocoutW((TObject*)0,Integration) << "RooBinIntegrator::RooBinIntegrator WARNING: integrand provide no binning definition observable #" 
				     << i << " substituting default binning of " << _numBins << " bins" << endl ;
      tmp = new list<Double_t> ;
      for (Int_t j=0 ; j<=_numBins ; j++) {
	tmp->push_back(_xmin[i]+j*(_xmax[i]-_xmin[i])/_numBins) ;
      }
    }
    _binb.push_back(tmp) ;
  }
  checkLimits();

} 


////////////////////////////////////////////////////////////////////////////////
/// Construct integrator on given function binding binding

RooBinIntegrator::RooBinIntegrator(const RooAbsFunc& function, const RooNumIntConfig& config) : 
  RooAbsIntegrator(function), _binb(0)
{
  const RooArgSet& configSet = config.getConfigSection(IsA()->GetName()) ;  
  _useIntegrandLimits= kTRUE;
  _numBins = (Int_t) configSet.getRealValue("numBins") ;
  assert(0 != integrand() && integrand()->isValid());
  
  // Allocate coordinate buffer size after number of function dimensions
  _x = new Double_t[_function->getDimension()] ;

  for (UInt_t i=0 ; i<_function->getDimension() ; i++) {
    _xmin.push_back(integrand()->getMinLimit(i));
    _xmax.push_back(integrand()->getMaxLimit(i));
    
    // Retrieve bin configuration from integrand
    list<Double_t>* tmp = integrand()->binBoundaries(i) ;
    if (!tmp) {
      oocoutW((TObject*)0,Integration) << "RooBinIntegrator::RooBinIntegrator WARNING: integrand provide no binning definition observable #" 
				     << i << " substituting default binning of " << _numBins << " bins" << endl ;
      tmp = new list<Double_t> ;
      for (Int_t j=0 ; j<=_numBins ; j++) {
	tmp->push_back(_xmin[i]+j*(_xmax[i]-_xmin[i])/_numBins) ;
      }
    }
    _binb.push_back(tmp) ;
  }

  checkLimits();
} 


////////////////////////////////////////////////////////////////////////////////
/// Clone integrator with new function binding and configuration. Needed by RooNumIntFactory

RooAbsIntegrator* RooBinIntegrator::clone(const RooAbsFunc& function, const RooNumIntConfig& config) const
{
  return new RooBinIntegrator(function,config) ;
}





////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooBinIntegrator::~RooBinIntegrator()
{
  if(_x) delete[] _x;
  for (vector<list<Double_t>*>::iterator iter = _binb.begin() ; iter!=_binb.end() ; ++iter) {
    delete (*iter) ;
  }

}


////////////////////////////////////////////////////////////////////////////////
/// Change our integration limits. Return kTRUE if the new limits are
/// ok, or otherwise kFALSE. Always returns kFALSE and does nothing
/// if this object was constructed to always use our integrand's limits.

Bool_t RooBinIntegrator::setLimits(Double_t *xmin, Double_t *xmax) 
{
  if(_useIntegrandLimits) {
    oocoutE((TObject*)0,Integration) << "RooBinIntegrator::setLimits: cannot override integrand's limits" << endl;
    return kFALSE;
  }
  _xmin[0]= *xmin;
  _xmax[0]= *xmax;
  return checkLimits();
}


////////////////////////////////////////////////////////////////////////////////
/// Check that our integration range is finite and otherwise return kFALSE.
/// Update the limits from the integrand if requested.

Bool_t RooBinIntegrator::checkLimits() const 
{
  if(_useIntegrandLimits) {
    assert(0 != integrand() && integrand()->isValid());
    _xmin.resize(_function->getDimension()) ;
    _xmax.resize(_function->getDimension()) ;
    for (UInt_t i=0 ; i<_function->getDimension() ; i++) {
      _xmin[i]= integrand()->getMinLimit(i);
      _xmax[i]= integrand()->getMaxLimit(i);
    }
  }
  for (UInt_t i=0 ; i<_function->getDimension() ; i++) {
    if (_xmax[i]<=_xmin[i]) {
      oocoutE((TObject*)0,Integration) << "RooBinIntegrator::checkLimits: bad range with min >= max (_xmin = " << _xmin[i] << " _xmax = " << _xmax[i] << ")" << endl;
      return kFALSE;
    }
    if (RooNumber::isInfinite(_xmin[i]) || RooNumber::isInfinite(_xmax[i])) {
      return kFALSE ;
    }
  }
  
  return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate numeric integral at given set of function binding parameters

Double_t RooBinIntegrator::integral(const Double_t *) 
{
  assert(isValid());

  double sum = 0. ;

  if (_function->getDimension()==1) {
    list<Double_t>::iterator iter = _binb[0]->begin() ;
    Double_t xlo = *iter ; ++iter ;
    for (; iter!=_binb[0]->end() ; ++iter) {
      Double_t xhi = *iter ;
      Double_t xcenter = (xhi+xlo)/2 ;
      Double_t binInt = integrand(xvec(xcenter))*(xhi-xlo) ;
      sum += binInt ;
      //cout << "RBI::integral over " << _function->getName() << " 1D binInt[" << xcenter << "] = " << binInt << " running sum = " << sum << endl ;
      xlo=xhi ;
    }
  }

  if (_function->getDimension()==2) {

    list<Double_t>::iterator iter1 = _binb[0]->begin() ;

    Double_t x1lo = *iter1 ; ++iter1 ;
    for (; iter1!=_binb[0]->end() ; ++iter1) {

      Double_t x1hi = *iter1 ;
      Double_t x1center = (x1hi+x1lo)/2 ;
      
      list<Double_t>::iterator iter2 = _binb[1]->begin() ;
      Double_t x2lo = *iter2 ; ++iter2 ;
      for (; iter2!=_binb[1]->end() ; ++iter2) {

	Double_t x2hi = *iter2 ;
	Double_t x2center = (x2hi+x2lo)/2 ;
      	
	Double_t binInt = integrand(xvec(x1center,x2center))*(x1hi-x1lo)*(x2hi-x2lo) ;
	//cout << "RBI::integral 2D binInt[" << x1center << "," << x2center << "] = " << binInt << " binv = " << (x1hi-x1lo) << "*" << (x2hi-x2lo) << endl ;
	sum += binInt ;
	x2lo=x2hi ;
      }
      x1lo=x1hi ;
    }    
  }

  if (_function->getDimension()==3) {

    list<Double_t>::iterator iter1 = _binb[0]->begin() ;

    Double_t x1lo = *iter1 ; ++iter1 ;
    for (; iter1!=_binb[0]->end() ; ++iter1) {

      Double_t x1hi = *iter1 ;
      Double_t x1center = (x1hi+x1lo)/2 ;
      
      list<Double_t>::iterator iter2 = _binb[1]->begin() ;
      Double_t x2lo = *iter2 ; ++iter2 ;
      for (; iter2!=_binb[1]->end() ; ++iter2) {

	Double_t x2hi = *iter2 ;
	Double_t x2center = (x2hi+x2lo)/2 ;

	list<Double_t>::iterator iter3 = _binb[2]->begin() ;
	Double_t x3lo = *iter3 ; ++iter3 ;
	for (; iter3!=_binb[2]->end() ; ++iter3) {

	  Double_t x3hi = *iter3 ;
	  Double_t x3center = (x3hi+x3lo)/2 ;
	  
	  Double_t binInt = integrand(xvec(x1center,x2center,x3center))*(x1hi-x1lo)*(x2hi-x2lo)*(x3hi-x3lo) ;
	  //cout << "RBI::integral 3D binInt[" << x1center << "," << x2center << "," << x3center << "] = " << binInt << endl ;
	  sum += binInt ;
	  
	  x3lo=x3hi ;
	}
	x2lo=x2hi ;
      }
      x1lo=x1hi ;
    }    
  }

  return sum;
}


