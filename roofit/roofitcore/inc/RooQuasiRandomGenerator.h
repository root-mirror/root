/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooQuasiRandomGenerator.h,v 1.7 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_QUASI_RANDOM_GENERATOR
#define ROO_QUASI_RANDOM_GENERATOR

#include "Rtypes.h"

class RooQuasiRandomGenerator {
public:
  RooQuasiRandomGenerator();
  virtual ~RooQuasiRandomGenerator();
  void reset();
  Bool_t generate(UInt_t dimension, Double_t vector[]);
  enum { MaxDimension = 12 , NBits = 31 , MaxDegree = 50 , MaxPrimitiveDegree = 5 };
protected:
  void calculateCoefs(UInt_t dimension);
  void calculateV(const int px[], int px_degree,
        int pb[], int * pb_degree, int v[], int maxv);
  void polyMultiply(const int pa[], int pa_degree, const int pb[],
          int pb_degree, int pc[], int  * pc_degree);
  // Z_2 field operations
  inline Int_t add(Int_t x, Int_t y) const { return (x+y)%2; }
  inline Int_t mul(Int_t x, Int_t y) const { return (x*y)%2; }
  inline Int_t sub(Int_t x, Int_t y) const { return add(x,y); }
private:
  Int_t *_nextq;
  Int_t _sequenceCount;

  static Bool_t _coefsCalculated;
  static Int_t _cj[NBits][MaxDimension];
  static const Int_t _primitivePoly[MaxDimension+1][MaxPrimitiveDegree+1];
  static const Int_t _polyDegree[MaxDimension+1];

  ClassDef(RooQuasiRandomGenerator,0) // Quasi-random number generator
};

#endif


