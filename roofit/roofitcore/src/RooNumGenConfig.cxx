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
\file RooNumGenConfig.cxx
\class RooNumGenConfig
\ingroup Roofitcore

RooNumGenConfig holds the configuration parameters of the various
numeric integrators used by RooRealIntegral. RooRealIntegral and RooAbsPdf
use this class in the (normalization) integral configuration interface
**/

#include "RooFit.h"
#include "Riostream.h"

#include "RooNumGenConfig.h"
#include "RooArgSet.h"
#include "RooAbsNumGenerator.h"
#include "RooNumGenFactory.h"
#include "RooMsgService.h"

#include "TClass.h"



using namespace std;

ClassImp(RooNumGenConfig);


////////////////////////////////////////////////////////////////////////////////
/// Return reference to instance of default numeric integrator configuration object

RooNumGenConfig& RooNumGenConfig::defaultConfig() 
{
  static RooNumGenConfig defaultConfig;
  return defaultConfig;
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor 

RooNumGenConfig::RooNumGenConfig() : 
  _method1D("method1D","1D sampling method"),
  _method1DCat("method1DCat","1D sampling method for pdfs with categories"),
  _method1DCond("method1DCond","1D sampling method for conditional pfs"),
  _method1DCondCat("method1DCond","1D sampling method for conditional pfs with categories"),
  _method2D("method2D","2D sampling method"),
  _method2DCat("method2DCat","2D sampling method for pdfs with categories"),
  _method2DCond("method2DCond","2D sampling method for conditional pfs"),
  _method2DCondCat("method2DCond","2D sampling method for conditional pfs with categories"),
  _methodND("methodND","ND sampling method"),
  _methodNDCat("methodNDCat","ND sampling method for pdfs with categories"),
  _methodNDCond("methodNDCond","ND sampling method for conditional pfs"),
  _methodNDCondCat("methodNDCond","ND sampling method for conditional pfs with categories")
{
  // Set all methods to undefined
  // Defined methods will be registered by static initialization routines
  // of the various numeric integrator engines
  _method1D.defineType("N/A",0) ;
  _method1DCat.defineType("N/A",0) ;
  _method1DCond.defineType("N/A",0) ;
  _method1DCondCat.defineType("N/A",0) ;

  _method2D.defineType("N/A",0) ;
  _method2DCat.defineType("N/A",0) ;
  _method2DCond.defineType("N/A",0) ;
  _method2DCondCat.defineType("N/A",0) ;

  _methodND.defineType("N/A",0) ;
  _methodNDCat.defineType("N/A",0) ;
  _methodNDCond.defineType("N/A",0) ;
  _methodNDCondCat.defineType("N/A",0) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooNumGenConfig::~RooNumGenConfig()
{
  // Delete all configuration data
  _configSets.Delete() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooNumGenConfig::RooNumGenConfig(const RooNumGenConfig& other) :
  TObject(other), RooPrintable(other),
  _method1D(other._method1D),
  _method1DCat(other._method1DCat),
  _method1DCond(other._method1DCond),
  _method1DCondCat(other._method1DCondCat),
  _method2D(other._method2D),
  _method2DCat(other._method2DCat),
  _method2DCond(other._method2DCond),
  _method2DCondCat(other._method2DCondCat),
  _methodND(other._methodND),
  _methodNDCat(other._methodNDCat),
  _methodNDCond(other._methodNDCond),
  _methodNDCondCat(other._methodNDCondCat)
{
  // Clone all configuration dat
  TIterator* iter = other._configSets.MakeIterator() ;
  RooArgSet* set ;
  while((set=(RooArgSet*)iter->Next())) {
    RooArgSet* setCopy = (RooArgSet*) set->snapshot() ;
    setCopy->setName(set->GetName()) ;
   _configSets.Add(setCopy);
  }
  delete iter ;
}


////////////////////////////////////////////////////////////////////////////////
/// Assignment operator from other RooNumGenConfig

RooNumGenConfig& RooNumGenConfig::operator=(const RooNumGenConfig& other) 
{
  // Prevent self-assignment 
  if (&other==this) {
    return *this ;
  }

  // Copy common properties
  _method1D.setIndex(other._method1D.getCurrentIndex()) ;
  _method1DCat.setIndex(other._method1DCat.getCurrentIndex()) ;
  _method1DCond.setIndex(other._method1DCond.getCurrentIndex()) ;
  _method1DCondCat.setIndex(other._method1DCondCat.getCurrentIndex()) ;

  _method2D.setIndex(other._method2D.getCurrentIndex()) ;
  _method2DCat.setIndex(other._method2DCat.getCurrentIndex()) ;
  _method2DCond.setIndex(other._method2DCond.getCurrentIndex()) ;
  _method2DCondCat.setIndex(other._method2DCondCat.getCurrentIndex()) ;

  _methodND.setIndex(other._methodND.getCurrentIndex()) ;
  _methodNDCat.setIndex(other._methodNDCat.getCurrentIndex()) ;
  _methodNDCond.setIndex(other._methodNDCond.getCurrentIndex()) ;
  _methodNDCondCat.setIndex(other._methodNDCondCat.getCurrentIndex()) ;

  // Delete old integrator-specific configuration data
  _configSets.Delete() ;

  // Copy new integrator-specific data
  TIterator* iter = other._configSets.MakeIterator() ;
  RooArgSet* set ;
  while((set=(RooArgSet*)iter->Next())) {
    RooArgSet* setCopy = (RooArgSet*) set->snapshot() ;
    setCopy->setName(set->GetName()) ;
   _configSets.Add(setCopy);
  }
  delete iter ;

  return *this ;
}




////////////////////////////////////////////////////////////////////////////////

RooCategory& RooNumGenConfig::method1D(Bool_t cond, Bool_t cat) 
{
  if (cond && cat) return _method1DCondCat ;
  if (cond) return _method1DCond ;
  if (cat) return _method1DCat ;
  return _method1D ;
}



////////////////////////////////////////////////////////////////////////////////

RooCategory& RooNumGenConfig::method2D(Bool_t cond, Bool_t cat) 
{
  if (cond && cat) return _method2DCondCat ;
  if (cond) return _method2DCond ;
  if (cat) return _method2DCat ;
  return _method2D ;
}



////////////////////////////////////////////////////////////////////////////////

RooCategory& RooNumGenConfig::methodND(Bool_t cond, Bool_t cat) 
{
  if (cond && cat) return _methodNDCondCat ;
  if (cond) return _methodNDCond ;
  if (cat) return _methodNDCat ;
  return _methodND ;
}



////////////////////////////////////////////////////////////////////////////////

const RooCategory& RooNumGenConfig::method1D(Bool_t cond, Bool_t cat) const 
{
  return const_cast<RooNumGenConfig*>(this)->method1D(cond,cat) ;
}



////////////////////////////////////////////////////////////////////////////////

const RooCategory& RooNumGenConfig::method2D(Bool_t cond, Bool_t cat) const 
{
  return const_cast<RooNumGenConfig*>(this)->method2D(cond,cat) ;
}



////////////////////////////////////////////////////////////////////////////////

const RooCategory& RooNumGenConfig::methodND(Bool_t cond, Bool_t cat) const 
{
  return const_cast<RooNumGenConfig*>(this)->methodND(cond,cat) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Add a configuration section for a particular integrator. Integrator name and capabilities are
/// automatically determined from instance passed as 'proto'. The defaultConfig object is associated
/// as the default configuration for the integrator. 

Bool_t RooNumGenConfig::addConfigSection(const RooAbsNumGenerator* proto, const RooArgSet& inDefaultConfig)
{
  std::string name = proto->IsA()->GetName();

  // Register integrator for appropriate dimensionalities
  
  _method1D.defineType(name) ;
  _method2D.defineType(name) ; 
  _methodND.defineType(name) ;

  if (proto->canSampleConditional()) {
    _method1DCond.defineType(name) ;
    _method2DCond.defineType(name) ;
    _methodNDCond.defineType(name) ;
  }
  if (proto->canSampleCategories()) {
    _method1DCat.defineType(name) ;
    _method2DCat.defineType(name) ;
    _methodNDCat.defineType(name) ;
  }

  if (proto->canSampleConditional() && proto->canSampleCategories()) {
    _method1DCondCat.defineType(name) ;
    _method2DCondCat.defineType(name) ;
    _methodNDCondCat.defineType(name) ;
  }
  
  // Store default configuration parameters
  RooArgSet* config = (RooArgSet*) inDefaultConfig.snapshot() ;
  config->setName(name.c_str());
  _configSets.Add(config) ;

  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return section with configuration parameters for integrator with given (class) name

RooArgSet& RooNumGenConfig::getConfigSection(const char* name)  
{
  return const_cast<RooArgSet&>((const_cast<const RooNumGenConfig*>(this)->getConfigSection(name))) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Retrieve configuration information specific to integrator with given name

const RooArgSet& RooNumGenConfig::getConfigSection(const char* name) const
{
  static RooArgSet dummy ;
  RooArgSet* config = (RooArgSet*) _configSets.FindObject(name) ;
  if (!config) {
    oocoutE((TObject*)0,InputArguments) << "RooNumGenConfig::getIntegrator: ERROR: no configuration stored for integrator '" << name << "'" << endl ;
    return dummy ;
  }
  return *config ;
}


////////////////////////////////////////////////////////////////////////////////

RooPrintable::StyleOption RooNumGenConfig::defaultPrintStyle(Option_t* opt) const
{ 
  if (!opt) {
    return kStandard ;
  }

  TString o(opt) ;
  o.ToLower() ;

  if (o.Contains("v")) {
    return kVerbose ;
  } 
  return kStandard ;
}



////////////////////////////////////////////////////////////////////////////////
/// Detailed printing interface

void RooNumGenConfig::printMultiline(ostream &os, Int_t /*content*/, Bool_t verbose, TString indent) const
{
  os << endl ;
  os << indent << "1-D sampling method: " << _method1D.getCurrentLabel() << endl ;
  if (_method1DCat.getCurrentIndex()!=_method1D.getCurrentIndex()) {
    os << " (" << _method1DCat.getCurrentLabel() << " if with categories)" << endl ;
  }
  if (_method1DCond.getCurrentIndex()!=_method1D.getCurrentIndex()) {
    os << " (" << _method1DCond.getCurrentLabel() << " if conditional)" << endl ;
  }
  if (_method1DCondCat.getCurrentIndex()!=_method1D.getCurrentIndex()) {
    os << " (" << _method1DCondCat.getCurrentLabel() << " if conditional with categories)" << endl ;    
  }
  os << endl ;

  os << indent << "2-D sampling method: " << _method2D.getCurrentLabel() << endl ;
  if (_method2DCat.getCurrentIndex()!=_method2D.getCurrentIndex()) {
    os << " (" << _method2DCat.getCurrentLabel() << " if with categories)" << endl ;
  }
  if (_method2DCond.getCurrentIndex()!=_method2D.getCurrentIndex()) {
    os << " (" << _method2DCond.getCurrentLabel() << " if conditional)" << endl ;
  }
  if (_method2DCondCat.getCurrentIndex()!=_method2D.getCurrentIndex()) {
    os << " (" << _method2DCondCat.getCurrentLabel() << " if conditional with categories)" << endl ;
  }
  os << endl ;

  os << indent << "N-D sampling method: " << _methodND.getCurrentLabel() << endl ;
  if (_methodNDCat.getCurrentIndex()!=_methodND.getCurrentIndex()) {
    os << " (" << _methodNDCat.getCurrentLabel() << " if with categories)" << endl ;
  }
  if (_methodNDCond.getCurrentIndex()!=_methodND.getCurrentIndex()) {
    os << " (" << _methodNDCond.getCurrentLabel() << " if conditional)" << endl ;
  }
  if (_methodNDCondCat.getCurrentIndex()!=_methodND.getCurrentIndex()) {
    os << " (" << _methodNDCondCat.getCurrentLabel() << " if conditional with categories)" << endl ;
  }
  os << endl ;
   
  if (verbose) {

    os << endl << "Available sampling methods:" << endl << endl ;
    TIterator* cIter = _configSets.MakeIterator() ;
    RooArgSet* configSet ;
    while ((configSet=(RooArgSet*)cIter->Next())) {

      os << indent << "*** " << configSet->GetName() << " ***" << endl ;
      os << indent << "Capabilities: " ;
      const RooAbsNumGenerator* proto = RooNumGenFactory::instance().getProtoSampler(configSet->GetName()) ;
     if (proto->canSampleConditional()) os << "[Conditional] " ;
     if (proto->canSampleCategories()) os << "[Categories] " ;
      os << endl ;

      os << "Configuration: " << endl ;
      configSet->printMultiline(os,kName|kValue|kTitle) ;
      os << endl ;

    }

    delete cIter ;
  }
}
