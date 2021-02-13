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
\file RooInt.cxx
\class RooInt
\ingroup Roofitcore

RooInt is a minimal implementation of a TObject holding a Int_t
value.
**/

#include "RooFit.h"

#include "RooInt.h"

using namespace std;

ClassImp(RooInt);
;



////////////////////////////////////////////////////////////////////////////////
/// Facilitate sorting of RooInts in ROOT container classes
/// Return -1 or +1 if 'other' is a RooInt with value
/// greater or lesser than self. Return zero if other
/// object is not a RooInt

Int_t RooInt::Compare(const TObject* other) const 
{
  const RooInt* otherD = dynamic_cast<const RooInt*>(other) ;
  if (!otherD) return 0 ;
  return (_value>otherD->_value) ? 1 : -1 ;
}
