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
#ifndef ROO_NUM_GEN_FACTORY
#define ROO_NUM_GEN_FACTORY

#include <map>
#include <string>
#include "TObject.h"
#include "RooLinkedList.h"
#include "RooAbsNumGenerator.h"
class RooNumGenConfig ;
class RooAbsReal ;

class RooNumGenFactory ;
typedef void (*RooNumGenInitializerFunc)(RooNumGenFactory&) ;

class RooNumGenFactory : public TObject {
public:

  static RooNumGenFactory& instance() ;
  virtual ~RooNumGenFactory();

  Bool_t storeProtoSampler(RooAbsNumGenerator* proto, const RooArgSet& defConfig) ;
  const RooAbsNumGenerator* getProtoSampler(const char* name) ;

  RooAbsNumGenerator* createSampler(RooAbsReal& func, const RooArgSet& genVars, const RooArgSet& condVars,
                const RooNumGenConfig& config, Bool_t verbose=kFALSE, RooAbsReal* maxFuncVal=0) ;


protected:

  friend class RooNumGenConfig ;

  std::map<std::string,RooAbsNumGenerator*> _map ;

  RooNumGenFactory();
  RooNumGenFactory(const RooNumGenFactory& other) ;


  ClassDef(RooNumGenFactory,1) // Numeric Generator factory
};

#endif


