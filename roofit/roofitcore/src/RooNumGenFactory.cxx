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
\file RooNumGenFactory.cxx
\class RooNumGenFactory
\ingroup Roofitcore

RooNumGenFactory is a factory to instantiate numeric integrators
from a given function binding and a given configuration. The factory
searches for a numeric integrator registered with the factory that
has the ability to perform the numeric integration. The choice of
method may depend on the number of dimensions integrated,
the nature of the integration limits (closed or open ended) and
the preference of the caller as encoded in the configuration object.
**/

#include "TClass.h"
#include "Riostream.h"

#include "RooFit.h"

#include "RooNumGenFactory.h"
#include "RooArgSet.h"
#include "RooAbsFunc.h"
#include "RooNumGenConfig.h"
#include "RooNumber.h"

#include "RooAcceptReject.h"
#include "RooFoamGenerator.h"


#include "RooMsgService.h"

using namespace std ;

ClassImp(RooNumGenFactory);



////////////////////////////////////////////////////////////////////////////////
/// Constructor. Register all known integrators by calling
/// their static registration functions

RooNumGenFactory::RooNumGenFactory()
{
  RooAcceptReject::registerSampler(*this) ;
  RooFoamGenerator::registerSampler(*this) ;

  // Prepare default
  RooNumGenConfig::defaultConfig().method1D(kFALSE,kFALSE).setLabel("RooFoamGenerator") ;
  RooNumGenConfig::defaultConfig().method1D(kTRUE ,kFALSE).setLabel("RooAcceptReject") ;
  RooNumGenConfig::defaultConfig().method1D(kFALSE,kTRUE ).setLabel("RooAcceptReject") ;
  RooNumGenConfig::defaultConfig().method1D(kTRUE, kTRUE ).setLabel("RooAcceptReject") ;
  
  RooNumGenConfig::defaultConfig().method2D(kFALSE,kFALSE).setLabel("RooFoamGenerator") ;
  RooNumGenConfig::defaultConfig().method2D(kTRUE ,kFALSE).setLabel("RooAcceptReject") ;
  RooNumGenConfig::defaultConfig().method2D(kFALSE,kTRUE ).setLabel("RooAcceptReject") ;
  RooNumGenConfig::defaultConfig().method2D(kTRUE, kTRUE ).setLabel("RooAcceptReject") ;
  
  RooNumGenConfig::defaultConfig().methodND(kFALSE,kFALSE).setLabel("RooFoamGenerator") ;
  RooNumGenConfig::defaultConfig().methodND(kTRUE ,kFALSE).setLabel("RooAcceptReject") ;
  RooNumGenConfig::defaultConfig().methodND(kFALSE,kTRUE ).setLabel("RooAcceptReject") ;
  RooNumGenConfig::defaultConfig().methodND(kTRUE, kTRUE ).setLabel("RooAcceptReject") ;
  
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooNumGenFactory::~RooNumGenFactory()
{
  std::map<std::string,RooAbsNumGenerator*>::iterator iter = _map.begin() ;
  while (iter != _map.end()) {
    delete iter->second ;
    ++iter ;
  }  
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooNumGenFactory::RooNumGenFactory(const RooNumGenFactory& other) : TObject(other)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Static method returning reference to singleton instance of factory

RooNumGenFactory& RooNumGenFactory::instance()
{
  static RooNumGenFactory instance;
  return instance;
}


////////////////////////////////////////////////////////////////////////////////
/// Method accepting registration of a prototype numeric integrator along with a RooArgSet of its
/// default configuration options and an optional list of names of other numeric integrators
/// on which this integrator depends. Returns true if integrator was previously registered

Bool_t RooNumGenFactory::storeProtoSampler(RooAbsNumGenerator* proto, const RooArgSet& defConfig) 
{
  TString name = proto->IsA()->GetName() ;

  if (getProtoSampler(name)) {
    //cout << "RooNumGenFactory::storeSampler() ERROR: integrator '" << name << "' already registered" << endl ;
    return kTRUE ;
  }

  // Add to factory 
  _map[name.Data()] = proto ;

  // Add default config to master config
  RooNumGenConfig::defaultConfig().addConfigSection(proto,defConfig) ;
  
  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return prototype integrator with given (class) name

const RooAbsNumGenerator* RooNumGenFactory::getProtoSampler(const char* name) 
{
  if (_map.count(name)==0) {
    return 0 ;
  } 
  
  return _map[name] ;
}



////////////////////////////////////////////////////////////////////////////////
/// Construct a numeric integrator instance that operates on function 'func' and is configured
/// with 'config'. If ndimPreset is greater than zero that number is taken as the dimensionality
/// of the integration, otherwise it is queried from 'func'. This function iterators over list
/// of available prototype integrators and returns an clone attached to the given function of
/// the first class that matches the specifications of the requested integration considering
/// the number of dimensions, the nature of the limits (open ended vs closed) and the user
/// preference stated in 'config'

RooAbsNumGenerator* RooNumGenFactory::createSampler(RooAbsReal& func, const RooArgSet& genVars, const RooArgSet& condVars, const RooNumGenConfig& config, Bool_t verbose, RooAbsReal* maxFuncVal) 
{
  // Find method defined configuration
  Int_t ndim = genVars.getSize() ;
  Bool_t cond = (condVars.getSize() > 0) ? kTRUE : kFALSE ;

  Bool_t hasCat(kFALSE) ;
  for (const auto arg : genVars) {
    if (arg->IsA()==RooCategory::Class()) {
      hasCat=kTRUE ;
      break ;
    }
  }


  TString method ;
  switch(ndim) {
  case 1:
    method = config.method1D(cond,hasCat).getCurrentLabel() ;
    break ;

  case 2:
    method = config.method2D(cond,hasCat).getCurrentLabel() ;
    break ;

  default:
    method = config.methodND(cond,hasCat).getCurrentLabel() ;
    break ;
  }

  // Check that a method was defined for this case
  if (!method.CompareTo("N/A")) {
    oocoutE((TObject*)0,Integration) << "RooNumGenFactory::createSampler: No sampler method has been defined for " 
				     << (cond?"a conditional ":"a ") << ndim << "-dimensional p.d.f" << endl ;
    return 0 ;    
  }

  // Retrieve proto integrator and return clone configured for the requested integration task
  const RooAbsNumGenerator* proto = getProtoSampler(method) ;  
  RooAbsNumGenerator* engine =  proto->clone(func,genVars,condVars,config,verbose,maxFuncVal) ;
  return engine ;
}
