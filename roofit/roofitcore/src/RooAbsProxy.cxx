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

#include "RooFit.h"

#include "Riostream.h"
#include "RooAbsProxy.h"
#include "RooArgSet.h"
#include "RooAbsArg.h"

/**
\file RooAbsProxy.cxx
\class RooAbsProxy
\ingroup Roofitcore

RooAbsProxy is the abstact interface for proxy classes.
Proxy classes hold pointers to other Roofit objects 
and process serverRedirect changes so that the proxied
pointers are updated accordingly on a clone or copy of
of the owning class
**/


using namespace std;

ClassImp(RooAbsProxy);
;


////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooAbsProxy::RooAbsProxy() : _nset(0)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooAbsProxy::RooAbsProxy(const char* /*name*/, const RooAbsProxy& other) : 
  _nset(other._nset)
{
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor

void RooAbsProxy::changeNormSet(const RooArgSet* newNormSet) 
{
  _nset = (RooArgSet*) newNormSet ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print proxy name

void RooAbsProxy::print(ostream& os, Bool_t /*addContents*/) const 
{ 
  os << name() << endl ; 
}
