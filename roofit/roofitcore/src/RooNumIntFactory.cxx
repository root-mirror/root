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
\file RooNumIntFactory.cxx
\class RooNumIntFactory
\ingroup Roofitcore

RooNumIntFactory is a factory to instantiate numeric integrators
from a given function binding and a given configuration. The factory
searches for a numeric integrator registered with the factory that
has the ability to perform the numeric integration. The choice of
method may depend on the number of dimensions integrated,
the nature of the integration limits (closed or open ended) and
the preference of the caller as encoded in the configuration object.
**/

#include "TClass.h"
#include "TSystem.h"
#include "Riostream.h"

#include "RooFit.h"

#include "RooNumIntFactory.h"
#include "RooArgSet.h"
#include "RooAbsFunc.h"
#include "RooNumIntConfig.h"
#include "RooNumber.h"

#include "RooIntegrator1D.h"
#include "RooBinIntegrator.h"
#include "RooIntegrator2D.h"
#include "RooSegmentedIntegrator1D.h"
#include "RooSegmentedIntegrator2D.h"
#include "RooImproperIntegrator1D.h"
#include "RooMCIntegrator.h"
//#include "RooGaussKronrodIntegrator1D.h"
//#include "RooAdaptiveGaussKronrodIntegrator1D.h"
#include "RooAdaptiveIntegratorND.h"

#include "RooMsgService.h"

using namespace std ;

ClassImp(RooNumIntFactory)



////////////////////////////////////////////////////////////////////////////////
/// Register all known integrators by calling
/// their static registration functions
void RooNumIntFactory::init() {
  RooBinIntegrator::registerIntegrator(*this) ;
  RooIntegrator1D::registerIntegrator(*this) ;
  RooIntegrator2D::registerIntegrator(*this) ;
  RooSegmentedIntegrator1D::registerIntegrator(*this) ;
  RooSegmentedIntegrator2D::registerIntegrator(*this) ;
  RooImproperIntegrator1D::registerIntegrator(*this) ;
  RooMCIntegrator::registerIntegrator(*this) ;
  // GSL integrator is now in RooFitMore and it register itself
  //RooAdaptiveGaussKronrodIntegrator1D::registerIntegrator(*this) ;
  //RooGaussKronrodIntegrator1D::registerIntegrator(*this) ;  
  RooAdaptiveIntegratorND::registerIntegrator(*this) ;

  RooNumIntConfig::defaultConfig().method1D().setLabel("RooIntegrator1D") ;
  RooNumIntConfig::defaultConfig().method1DOpen().setLabel("RooImproperIntegrator1D") ;
  RooNumIntConfig::defaultConfig().method2D().setLabel("RooAdaptiveIntegratorND") ;
  RooNumIntConfig::defaultConfig().methodND().setLabel("RooAdaptiveIntegratorND") ;

  //if GSL is available load (and register GSL integrator)
#ifdef R__HAS_MATHMORE
  int iret = gSystem->Load("libRooFitMore");
  if (iret < 0) {
     oocoutE((TObject*)nullptr, Integration) << " RooNumIntFactory::Init : libRooFitMore cannot be loaded. GSL integrators will not beavailable ! " << std::endl;
  }
#endif
}


////////////////////////////////////////////////////////////////////////////////
/// Static method returning reference to singleton instance of factory

RooNumIntFactory& RooNumIntFactory::instance()
{
  static unique_ptr<RooNumIntFactory> instance;

  if (!instance) {
    // This is needed to break a deadlock. During init(),
    // other functions may call back to this one here. So we need to construct first,
    // and ensure that we can return an instance while we are waiting for init()
    // to finish.
    instance.reset(new RooNumIntFactory);
    instance->init();
  }

  return *instance;
}



////////////////////////////////////////////////////////////////////////////////
/// Method accepting registration of a prototype numeric integrator along with a RooArgSet of its
/// default configuration options and an optional list of names of other numeric integrators
/// on which this integrator depends. Returns true if integrator was previously registered

Bool_t RooNumIntFactory::storeProtoIntegrator(RooAbsIntegrator* proto, const RooArgSet& defConfig, const char* depName) 
{
  TString name = proto->IsA()->GetName() ;

  if (getProtoIntegrator(name)) {
    //cout << "RooNumIntFactory::storeIntegrator() ERROR: integrator '" << name << "' already registered" << endl ;
    return kTRUE ;
  }

  // Add to factory
  _map[name.Data()] = std::make_pair(unique_ptr<RooAbsIntegrator>(proto), std::string(depName));

  // Add default config to master config
  RooNumIntConfig::defaultConfig().addConfigSection(proto,defConfig) ;
  
  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return prototype integrator with given (class) name

const RooAbsIntegrator* RooNumIntFactory::getProtoIntegrator(const char* name) const
{
  auto item = _map.find(name);
  
  return item == _map.end() ? nullptr : item->second.first.get();
}



////////////////////////////////////////////////////////////////////////////////
/// Get list of class names of integrators needed by integrator named 'name'

const char* RooNumIntFactory::getDepIntegratorName(const char* name) const
{
  auto item = _map.find(name);

  return item == _map.end() ? nullptr : item->second.second.c_str();
}



////////////////////////////////////////////////////////////////////////////////
/// Construct a numeric integrator instance that operates on function 'func' and is configured
/// with 'config'. If ndimPreset is greater than zero that number is taken as the dimensionality
/// of the integration, otherwise it is queried from 'func'. This function iterators over list
/// of available prototype integrators and returns an clone attached to the given function of
/// the first class that matches the specifications of the requested integration considering
/// the number of dimensions, the nature of the limits (open ended vs closed) and the user
/// preference stated in 'config'

RooAbsIntegrator* RooNumIntFactory::createIntegrator(RooAbsFunc& func, const RooNumIntConfig& config, Int_t ndimPreset, Bool_t isBinned) const
{
  // First determine dimensionality and domain of integrand  
  Int_t ndim = ndimPreset>0 ? ndimPreset : ((Int_t)func.getDimension()) ;

  Bool_t openEnded = kFALSE ;
  Int_t i ;
  for (i=0 ; i<ndim ; i++) {
    if(RooNumber::isInfinite(func.getMinLimit(i)) ||
       RooNumber::isInfinite(func.getMaxLimit(i))) {
      openEnded = kTRUE ;
    }
  }

  // Find method defined configuration
  TString method ;
  switch(ndim) {
  case 1:
    method = openEnded ? config.method1DOpen().getCurrentLabel() : config.method1D().getCurrentLabel() ;
    break ;

  case 2:
    method = openEnded ? config.method2DOpen().getCurrentLabel() : config.method2D().getCurrentLabel() ;
    break ;

  default:
    method = openEnded ? config.methodNDOpen().getCurrentLabel() : config.methodND().getCurrentLabel() ;
    break ;
  }

  // If distribution is binned and not open-ended override with bin integrator
  if (isBinned & !openEnded) {
    method = "RooBinIntegrator" ;
  }

  // Check that a method was defined for this case
  if (!method.CompareTo("N/A")) {
    oocoutE((TObject*)0,Integration) << "RooNumIntFactory::createIntegrator: No integration method has been defined for " 
				     << (openEnded?"an open ended ":"a ") << ndim << "-dimensional integral" << endl ;
    return 0 ;    
  }

  // Retrieve proto integrator and return clone configured for the requested integration task
  const RooAbsIntegrator* proto = getProtoIntegrator(method) ;  
  RooAbsIntegrator* engine =  proto->clone(func,config) ;
  if (config.printEvalCounter()) {
    engine->setPrintEvalCounter(kTRUE) ;
  }
  return engine ;
}
