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
\file RooSecondMoment.cxx
\class RooSecondMoment
\ingroup Roofitcore

RooSecondMoment represents the first, second, or third order derivative
of any RooAbsReal as calculated (numerically) by the MathCore Richardson
derivator class.
**/


#include "RooFit.h"

#include "Riostream.h"
#include <math.h>

#include "RooSecondMoment.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooErrorHandler.h"
#include "RooArgSet.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooFunctor.h"
#include "RooGlobalFunc.h"
#include "RooConstVar.h"
#include "RooRealIntegral.h"
#include "RooNumIntConfig.h"
#include "RooFormulaVar.h"
#include "RooLinearVar.h"
#include "RooProduct.h"
#include <string>
using namespace std;


ClassImp(RooSecondMoment);



////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooSecondMoment::RooSecondMoment() 
{
}



////////////////////////////////////////////////////////////////////////////////

RooSecondMoment::RooSecondMoment(const char* name, const char* title, RooAbsReal& func, RooRealVar& x, Bool_t centr, Bool_t takeRoot) :
  RooAbsMoment(name, title,func,x,2,takeRoot),
  _xf("!xf","xf",this,kFALSE,kFALSE),
  _ixf("!ixf","ixf",this),
  _if("!if","if",this),
  _xfOffset(0)
{
  setExpensiveObjectCache(func.expensiveObjectCache()) ;
  
  RooAbsReal* XF(0) ;
  if (centr) {

    string m1name=Form("%s_moment1",GetName()) ;
    RooAbsReal* mom1 = func.mean(x) ;
    _mean.setArg(*mom1) ;
    
    string pname=Form("%s_product",name) ;
    _xfOffset = mom1->getVal() ;
    XF = new RooFormulaVar(pname.c_str(),Form("pow((@0-%f),2)*@1",_xfOffset),RooArgList(x,func)) ;  
        
  } else {

    string pname=Form("%s_product",name) ;
    XF = new RooProduct(pname.c_str(),pname.c_str(),RooArgList(x,x,func)) ;
  }
    
  XF->setExpensiveObjectCache(func.expensiveObjectCache()) ;

  if (func.isBinnedDistribution(x)) {
    XF->specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator");
  }

  RooRealIntegral* intXF = (RooRealIntegral*) XF->createIntegral(x) ;
  RooRealIntegral* intF =  (RooRealIntegral*) func.createIntegral(x) ;
  intXF->setCacheNumeric(kTRUE) ;
  intF->setCacheNumeric(kTRUE) ;

  _xf.setArg(*XF) ;
  _ixf.setArg(*intXF) ;
  _if.setArg(*intF) ;
  addOwnedComponents(RooArgSet(*XF,*intXF,*intF)) ;
}

////////////////////////////////////////////////////////////////////////////////

RooSecondMoment::RooSecondMoment(const char* name, const char* title, RooAbsReal& func, RooRealVar& x, const RooArgSet& nset, 
		     Bool_t centr, Bool_t takeRoot, Bool_t intNSet) :
  RooAbsMoment(name, title,func,x,2,takeRoot),
  _xf("!xf","xf",this,kFALSE,kFALSE),
  _ixf("!ixf","ixf",this),
  _if("!if","if",this),
  _xfOffset(0)
{
  setExpensiveObjectCache(func.expensiveObjectCache()) ;

  _nset.add(nset) ;

  RooAbsReal* XF(0) ;
  if (centr) {

    string m1name=Form("%s_moment1",GetName()) ;
    RooAbsReal* mom1 = func.mean(x,nset) ;
    _mean.setArg(*mom1) ;
    
    string pname=Form("%s_product",name) ;
    _xfOffset = mom1->getVal() ;
    XF = new RooFormulaVar(pname.c_str(),Form("pow((@0-%f),2)*@1",_xfOffset),RooArgList(x,func)) ;  


  } else {

    string pname=Form("%s_product",name) ;
    XF = new RooProduct(pname.c_str(),pname.c_str(),RooArgList(x,x,func)) ;

  }
  
  XF->setExpensiveObjectCache(func.expensiveObjectCache()) ;

  if (func.isBinnedDistribution(x)) {
    XF->specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator");
  }
  if (intNSet && _nset.getSize()>0 && func.isBinnedDistribution(_nset)) {
      XF->specialIntegratorConfig(kTRUE)->method2D().setLabel("RooBinIntegrator");
      XF->specialIntegratorConfig(kTRUE)->methodND().setLabel("RooBinIntegrator");
  }

  RooArgSet intSet(x) ;
  if (intNSet) intSet.add(_nset,kTRUE) ;
  RooRealIntegral* intXF = (RooRealIntegral*) XF->createIntegral(intSet,&_nset) ;
  RooRealIntegral* intF =  (RooRealIntegral*) func.createIntegral(intSet,&_nset) ;
  intXF->setCacheNumeric(kTRUE) ;
  intF->setCacheNumeric(kTRUE) ;

  _xf.setArg(*XF) ;
  _ixf.setArg(*intXF) ;
  _if.setArg(*intF) ;
  addOwnedComponents(RooArgSet(*XF,*intXF,*intF)) ;
}



////////////////////////////////////////////////////////////////////////////////

RooSecondMoment::RooSecondMoment(const RooSecondMoment& other, const char* name) :
  RooAbsMoment(other, name), 
  _xf("xf",this,other._xf),
  _ixf("ixf",this,other._ixf),
  _if("if",this,other._if),
  _xfOffset(other._xfOffset)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooSecondMoment::~RooSecondMoment() 
{
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate value  

Double_t RooSecondMoment::evaluate() const 
{
  Double_t ratio = _ixf / _if ;

  if (_mean.absArg()) {
    ratio -= (_mean - _xfOffset)*(_mean-_xfOffset) ;
  }  

  Double_t ret =  _takeRoot ? sqrt(ratio) : ratio ;
  return ret ;
}


