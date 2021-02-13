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
\file RooRealConstant.cxx
\class RooRealConstant
\ingroup Roofitcore

RooRealConstant provides static functions to create and keep track
of RooRealVar constants. Instead of creating such constants by
hand (e.g. RooRealVar one("one","one",1)), simply use
~~~{.cpp}
 RooRealConstant::value(1.0)
~~~
whenever a reference to RooRealVar with constant value 1.0 is needed.
RooRealConstant keeps an internal database of previously created
RooRealVar objects and will recycle them as appropriate.
**/

#include "RooFit.h"

#include <math.h>
#include <sstream>
#include "RooRealConstant.h"
#include "RooConstVar.h"
#include "RooArgList.h"

using namespace std;

ClassImp(RooRealConstant);



////////////////////////////////////////////////////////////////////////////////
/// Return a constant value object with given value.
/// Return previously created object if available,
/// otherwise create a new one on the fly.

RooConstVar& RooRealConstant::value(Double_t value)
{
  // Lookup existing constant
  for (auto varArg : constDB()) {
    auto var = static_cast<RooConstVar*>(varArg);
    if ((var->getVal()==value) && (!var->getAttribute("REMOVAL_DUMMY"))) return *var ;
  }

  // Create new constant
  std::ostringstream s ;
  s << value ;

  auto var = new RooConstVar(s.str().c_str(),s.str().c_str(),value) ;
  var->setAttribute("RooRealConstant_Factory_Object",kTRUE) ;
  constDB().addOwned(*var) ;

  return *var ;
}


////////////////////////////////////////////////////////////////////////////////
/// Create a dummy node used in node-removal operations

RooConstVar& RooRealConstant::removalDummy()
{
  RooConstVar* var = new RooConstVar("REMOVAL_DUMMY","REMOVAL_DUMMY",1) ;
  var->setAttribute("RooRealConstant_Factory_Object",kTRUE) ;
  var->setAttribute("REMOVAL_DUMMY") ;
  constDB().addOwned(*var) ;

  return *var ;
}



////////////////////////////////////////////////////////////////////////////////
/// One-time initialization of constants database

RooArgList& RooRealConstant::constDB()
{
  static RooArgList constDB("RooRealVar Constants Database");
  return constDB;
}
