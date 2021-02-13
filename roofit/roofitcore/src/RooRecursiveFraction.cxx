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
\file RooRecursiveFraction.cxx
\class RooRecursiveFraction
\ingroup Roofitcore

Class RooRecursiveFraction is a RooAbsReal implementation that
calculates the plain fraction of sum of RooAddPdf components
from a set of recursive fractions: for a given set of input fractions
\f$ {a_i} \f$, it returns \f$ a_n * \prod_{i=0}^{n-1} (1 - a_i) \f$.
**/


#include "RooFit.h"

#include "Riostream.h"
#include <math.h>

#include "RooRecursiveFraction.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooErrorHandler.h"
#include "RooArgSet.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMsgService.h"

using namespace std;

ClassImp(RooRecursiveFraction);



////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooRecursiveFraction::RooRecursiveFraction()
{

}



////////////////////////////////////////////////////////////////////////////////
/// Constructor of plain RooAddPdf fraction from list of recursive fractions

RooRecursiveFraction::RooRecursiveFraction(const char* name, const char* title, const RooArgList& fracList) :
  RooAbsReal(name, title),
  _list("list","First set of components",this)
{
  for (Int_t ifrac=fracList.getSize()-1 ; ifrac>=0 ; ifrac--) {
    RooAbsArg* comp = fracList.at(ifrac) ;
    if (!dynamic_cast<RooAbsReal*>(comp)) {
      std::stringstream errorMsg;
      errorMsg << "RooRecursiveFraction::ctor(" << GetName() << ") ERROR: component " << comp->GetName()
			    << " is not of type RooAbsReal" << endl ;
      coutE(InputArguments) << errorMsg.str();
      throw std::invalid_argument(errorMsg.str());
    }

    _list.add(*comp) ;    
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooRecursiveFraction::RooRecursiveFraction(const RooRecursiveFraction& other, const char* name) :
  RooAbsReal(other, name), 
  _list("list",this,other._list)
{

}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooRecursiveFraction::~RooRecursiveFraction() 
{

}



////////////////////////////////////////////////////////////////////////////////
/// Calculate and return value of \f$ a_n * \prod_{i=0}^{n-1} (1 - a_i) \f$.
Double_t RooRecursiveFraction::evaluate() const 
{
  const RooArgSet* nset = _list.nset() ;

  // Note that input coefficients are saved in reverse in this list.
  Double_t prod = static_cast<RooAbsReal&>(_list[0]).getVal(nset);

  for (unsigned int i=1; i < _list.size(); ++i) {
    prod *= (1 - static_cast<RooAbsReal&>(_list[i]).getVal(nset));
  }
    
  return prod ;
}

