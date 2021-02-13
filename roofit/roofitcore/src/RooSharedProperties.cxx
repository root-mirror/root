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
\file RooSharedProperties.cxx
\class RooSharedProperties
\ingroup Roofitcore

Class RooSharedProperties is the base class for shared properties
that can be stored in RooSharedPropertiesList.
**/

#include "RooFit.h"
#include "RooSharedProperties.h"
#include "RooMsgService.h"
#include "RooTrace.h"

#include "Riostream.h"
using std::cout ;
using std::endl ;

using namespace std;

ClassImp(RooSharedProperties);
;



////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooSharedProperties::RooSharedProperties() : _refCount(0), _inSharedList(kFALSE)
{
  RooTrace::create(this) ;
} 


////////////////////////////////////////////////////////////////////////////////
/// Constructor with unique-id string

RooSharedProperties::RooSharedProperties(const char* uuidstr) : _uuid(uuidstr), _refCount(0), _inSharedList(kFALSE)
{
  RooTrace::create(this) ;
} 



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooSharedProperties::~RooSharedProperties() 
{
  RooTrace::destroy(this) ;
} 



////////////////////////////////////////////////////////////////////////////////
/// Return true of unique id of this property is equal to that of other

Bool_t RooSharedProperties::operator==(const RooSharedProperties& other) const
{
  return (_uuid==other._uuid) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Printing interface

void RooSharedProperties::Print(Option_t* /*opts*/) const 
{
  cout << "RooSharedProperties(" << this << ") UUID = " << _uuid.AsString() << endl ;
}
