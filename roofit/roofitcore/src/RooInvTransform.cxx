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
\file RooInvTransform.cxx
\class RooInvTransform
\ingroup Roofitcore

Lightweight function binding that returns the inverse of an input function binding.
Apply the change of variables transformation x -> 1/x to the input
function and its range. The function must be one dimensional and its
range cannot include zero.
**/


#include "RooFit.h"

#include "RooInvTransform.h"

#include "Riostream.h"
#include <math.h>

using namespace std;

ClassImp(RooInvTransform);
;


////////////////////////////////////////////////////////////////////////////////
/// Constructor. 

RooInvTransform::RooInvTransform(const RooAbsFunc &func) :
  RooAbsFunc(func.getDimension()), _func(&func)
{
}
