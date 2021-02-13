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
\file RooConstraintSum.cxx
\class RooConstraintSum
\ingroup Roofitcore

RooConstraintSum calculates the sum of the -(log) likelihoods of
a set of RooAbsPfs that represent constraint functions. This class
is used to calculate the composite -log(L) of constraints to be
added to the regular -log(L) in RooAbsPdf::fitTo() with Constrain(..)
arguments.
**/


#include "RooFit.h"

#include "Riostream.h"
#include <math.h>

#include "RooConstraintSum.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooErrorHandler.h"
#include "RooArgSet.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMsgService.h"

using namespace std;

ClassImp(RooConstraintSum);



////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooConstraintSum::RooConstraintSum()
{

}




////////////////////////////////////////////////////////////////////////////////
/// Constructor with set of constraint p.d.f.s. All elements in constraintSet must inherit from RooAbsPdf

RooConstraintSum::RooConstraintSum(const char* name, const char* title, const RooArgSet& constraintSet, const RooArgSet& normSet) :
  RooAbsReal(name, title),
  _set1("set1","First set of components",this),
  _paramSet("paramSet","Set of parameters",this)
{
  for (const auto comp : constraintSet) {
    if (!dynamic_cast<RooAbsPdf*>(comp)) {
      coutE(InputArguments) << "RooConstraintSum::ctor(" << GetName() << ") ERROR: component " << comp->GetName() 
			    << " is not of type RooAbsPdf" << endl ;
      RooErrorHandler::softAbort() ;
    }
    _set1.add(*comp) ;
  }

  _paramSet.add(normSet) ;
}





////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooConstraintSum::RooConstraintSum(const RooConstraintSum& other, const char* name) :
  RooAbsReal(other, name), 
  _set1("set1",this,other._set1),
  _paramSet("paramSet",this,other._paramSet)
{

}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooConstraintSum::~RooConstraintSum() 
{

}



////////////////////////////////////////////////////////////////////////////////
/// Return sum of -log of constraint p.d.f.s

Double_t RooConstraintSum::evaluate() const 
{
  Double_t sum(0);

  for (const auto comp : _set1) {
    sum -= static_cast<RooAbsPdf*>(comp)->getLogVal(&_paramSet);
  }
  
  return sum;
}

