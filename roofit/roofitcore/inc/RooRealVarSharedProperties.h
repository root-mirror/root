/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooRealVarSharedProperties.h,v 1.3 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_REAL_VAR_SHARED_PROPERTY
#define ROO_REAL_VAR_SHARED_PROPERTY

#include "RooSharedProperties.h"

#include <memory>
#include <unordered_map>
#include <string>

class RooAbsBinning;

class RooRealVarSharedProperties : public RooSharedProperties {
public:

  RooRealVarSharedProperties() ;
  RooRealVarSharedProperties(const RooRealVarSharedProperties&);
  RooRealVarSharedProperties(const char* uuidstr) ;
  virtual ~RooRealVarSharedProperties() ;
  void disownBinnings() {
    _ownBinnings = false;
  }

  RooSharedProperties* clone() {
    auto tmp = new RooRealVarSharedProperties(*this);
    tmp->disownBinnings();
    return tmp;
  }

protected:

  friend class RooRealVar ;

  std::unordered_map<std::string,RooAbsBinning*> _altBinning ;  // Optional alternative ranges and binnings
  bool _ownBinnings{true}; //!
  ClassDef(RooRealVarSharedProperties,2) // Shared properties of a RooRealVar clone set
};


#endif
