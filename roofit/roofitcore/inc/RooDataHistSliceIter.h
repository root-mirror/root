/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
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
#ifndef ROO_DATA_HIST_SLICE_ITER
#define ROO_DATA_HIST_SLICE_ITER

#include "TIterator.h"
#include "RooArgSet.h"

class RooDataHist;

typedef TIterator* pTIterator ;

class RooDataHistSliceIter : public TIterator {
public:
  // Constructors, assignment etc.
  RooDataHistSliceIter(const RooDataHistSliceIter& other) ;
  virtual ~RooDataHistSliceIter() ;

  // Iterator implementation
  virtual const TCollection* GetCollection() const ;
  virtual TObject* Next() ;
  virtual void Reset() ;
  virtual bool operator!=(const TIterator &aIter) const ;
  virtual TObject *operator*() const ;
protected:
  
  friend class RooDataHist ;
  RooDataHistSliceIter(RooDataHist& hist, RooAbsArg& sliceArg) ;

  RooDataHist* _hist ;
  RooAbsArg* _sliceArg ;  
  Int_t      _baseIndex ;
  Int_t      _stepSize ;
  Int_t      _nStep ;
  Int_t      _curStep ;
  
  TIterator& operator=(const TIterator&) { 
    // Prohibit iterator assignment 
    return *this ; 
  }

  ClassDef(RooDataHistSliceIter,0) // Iterator over a one-dimensional slice of a RooDataHist
};

#endif
