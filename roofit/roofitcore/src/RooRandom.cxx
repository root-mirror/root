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
\file RooRandom.cxx
\class RooRandom
\ingroup Roofitcore

This class provides a static interface for generating random numbers.
By default a private copy of TRandom3 is used to generate all random numbers.
**/
#include <cassert>

#include "RooFit.h"

#include "RooRandom.h"
#include "RooQuasiRandomGenerator.h"

#include "TRandom3.h"

using namespace std;

ClassImp(RooRandom);
  ;


TRandom* RooRandom::_theGenerator = 0;
RooQuasiRandomGenerator* RooRandom::_theQuasiGenerator = 0;
RooRandom::Guard RooRandom::guard;

////////////////////////////////////////////////////////////////////////////////

RooRandom::Guard::~Guard()
{ delete RooRandom::_theGenerator; delete RooRandom::_theQuasiGenerator; }

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to a singleton random-number generator
/// implementation. Creates the object the first time it is called.

TRandom *RooRandom::randomGenerator() 
{
  if (!_theGenerator) _theGenerator= new TRandom3();
  return _theGenerator;
}


////////////////////////////////////////////////////////////////////////////////
/// set the random number generator; takes ownership of the object passed as parameter

void RooRandom::setRandomGenerator(TRandom* gen)
{
  if (_theGenerator) delete _theGenerator;
  _theGenerator = gen;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to a singleton quasi-random generator
/// implementation. Creates the object the first time it is called.

RooQuasiRandomGenerator *RooRandom::quasiGenerator() 
{
  if(!_theQuasiGenerator) _theQuasiGenerator= new RooQuasiRandomGenerator();
  return _theQuasiGenerator;
}


////////////////////////////////////////////////////////////////////////////////
/// Return a number uniformly distributed from (0,1)

Double_t RooRandom::uniform(TRandom *generator) 
{
  return generator->Rndm();
}


////////////////////////////////////////////////////////////////////////////////
/// Fill the vector provided with random numbers uniformly distributed from (0,1)

void RooRandom::uniform(UInt_t dimension, Double_t vector[], TRandom *generator) 
{
  generator->RndmArray(dimension, vector);
}


////////////////////////////////////////////////////////////////////////////////
/// Return an integer uniformly distributed from [0,n-1]

UInt_t RooRandom::integer(UInt_t n, TRandom *generator) 
{
  return generator->Integer(n);
}


////////////////////////////////////////////////////////////////////////////////
/// Return a Gaussian random variable with mean 0 and variance 1.

Double_t RooRandom::gaussian(TRandom *generator) 
{
  return generator->Gaus();
}


////////////////////////////////////////////////////////////////////////////////
/// Return a quasi-random number in the range (0,1) using the
/// Niederreiter base 2 generator described in Bratley, Fox, Niederreiter,
/// ACM Trans. Model. Comp. Sim. 2, 195 (1992).

Bool_t RooRandom::quasi(UInt_t dimension, Double_t vector[], RooQuasiRandomGenerator *generator) 
{
  return generator->generate(dimension,vector);
}
