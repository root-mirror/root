/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooSetPair.h,v 1.4 2007/08/09 19:55:47 wouter Exp $
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
#ifndef ROO_SET_PAIR
#define ROO_SET_PAIR

#include <iostream>
#include "TObject.h"
#include "RooArgSet.h"

class RooLinkedListElem ;
class TBuffer ;

class RooSetPair : public TObject {
public:

  // Initial element ctor
  RooSetPair(const RooArgSet* set1=0, const RooArgSet* set2=0) :
    _set1(const_cast<RooArgSet*>(set1)), _set2(const_cast<RooArgSet*>(set2)) {
  }

  // Destructor
  virtual ~RooSetPair();

  RooArgSet* _set1 ;
  RooArgSet* _set2 ;

  virtual ULong_t Hash() const {
    return TString::Hash((void*)&_set1,2*sizeof(void*)) ;
  }

protected:


  // Forbidden
  RooSetPair(const RooSetPair&) ;

  ClassDef(RooSetPair,0) // Utility class holding a pair of RooArgSet pointers
} ;



#endif
