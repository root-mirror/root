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
#ifndef ROO_CHANGE_TRACKER
#define ROO_CHANGE_TRACKER

#include "RooAbsReal.h"
#include "RooListProxy.h"
#include <vector>

class RooRealVar;
class RooArgList ;

class RooChangeTracker : public RooAbsReal {
public:

  RooChangeTracker() ;
  RooChangeTracker(const char *name, const char *title, const RooArgSet& trackSet, Bool_t checkValues=kFALSE) ;
  virtual ~RooChangeTracker() ;

  RooChangeTracker(const RooChangeTracker& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooChangeTracker(*this, newname); }

  Bool_t hasChanged(Bool_t clearState) ;

  RooArgSet parameters() const ;


protected:

  RooListProxy     _realSet ;        // List of reals to track
  RooListProxy     _catSet ;         // List of categories to check
  std::vector<Double_t> _realRef ;   // Reference values for reals
  std::vector<Int_t>    _catRef ;    // Reference values for categories
  Bool_t       _checkVal ;           // Check contents as well if true

  Bool_t        _init ; //!

  Double_t evaluate() const { return 1 ; }

  ClassDef(RooChangeTracker,1) // Meta object that tracks changes in set of other arguments
};

#endif
