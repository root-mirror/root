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
\file RooFoamGenerator.cxx
\class RooFoamGenerator
\ingroup Roofitcore

Class RooFoamGenerator is a generic toy monte carlo generator that implement
the TFOAM sampling technique on any positively valued function.
The RooFoamGenerator generator is used by the various generator context
classes to take care of generation of observables for which p.d.fs
do not define internal methods.

The foam generator reacts to the following config options:
- nCell[123N]D
- nSample
- chatLevel
Access those using:
    myPdf->specialGeneratorConfig()->getConfigSection("RooFoamGenerator").setRealValue("nSample",1e4);

\see rf902_numgenconfig.C
**/


#include "RooFit.h"
#include "Riostream.h"

#include "RooFoamGenerator.h"
#include "RooAbsReal.h"
#include "RooCategory.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooRandom.h"
#include "RooErrorHandler.h"

#include "TString.h"
#include "TIterator.h"
#include "RooMsgService.h"
#include "TClass.h"
#include "TFoam.h"
#include "RooTFoamBinding.h"
#include "RooNumGenFactory.h"
#include "RooNumGenConfig.h"

#include <assert.h>

using namespace std;

ClassImp(RooFoamGenerator);
  ;


////////////////////////////////////////////////////////////////////////////////
/// Register RooIntegrator1D, is parameters and capabilities with RooNumIntFactory

void RooFoamGenerator::registerSampler(RooNumGenFactory& fact)
{
  // Register RooIntegrator1D, is parameters and capabilities with RooNumIntFactory
  RooRealVar nSample("nSample","Number of samples per cell",200,0,1e6) ;
  RooRealVar nCell1D("nCell1D","Number of cells for 1-dim generation",30,0,1e6) ;
  RooRealVar nCell2D("nCell2D","Number of cells for 2-dim generation",500,0,1e6) ;
  RooRealVar nCell3D("nCell3D","Number of cells for 3-dim generation",5000,0,1e6) ;
  RooRealVar nCellND("nCellND","Number of cells for N-dim generation",10000,0,1e6) ;
  RooRealVar chatLevel("chatLevel","TFOAM 'chat level' (verbosity)",0,0,2) ;

  RooFoamGenerator* proto = new RooFoamGenerator ;
  fact.storeProtoSampler(proto,RooArgSet(nSample,nCell1D,nCell2D,nCell3D,nCellND,chatLevel)) ;
}




////////////////////////////////////////////////////////////////////////////////

RooFoamGenerator::RooFoamGenerator(const RooAbsReal &func, const RooArgSet &genVars, const RooNumGenConfig& config, Bool_t verbose, const RooAbsReal* maxFuncVal) :
  RooAbsNumGenerator(func,genVars,verbose,maxFuncVal)
{
  _binding = new RooTFoamBinding(*_funcClone,_realVars) ;
 
  _tfoam = new TFoam("TFOAM") ;
  _tfoam->SetkDim(_realVars.getSize()) ;
  _tfoam->SetRho(_binding) ;
  _tfoam->SetPseRan(RooRandom::randomGenerator()) ;
  switch(_realVars.getSize()) {
  case 1:_tfoam->SetnCells((Int_t)config.getConfigSection("RooFoamGenerator").getRealValue("nCell1D")) ; break ;
  case 2:_tfoam->SetnCells((Int_t)config.getConfigSection("RooFoamGenerator").getRealValue("nCell2D")) ; break ;
  case 3:_tfoam->SetnCells((Int_t)config.getConfigSection("RooFoamGenerator").getRealValue("nCell3D")) ; break ;
  default:_tfoam->SetnCells((Int_t)config.getConfigSection("RooFoamGenerator").getRealValue("nCellND")) ; break ;
  }
  _tfoam->SetnSampl((Int_t)config.getConfigSection("RooFoamGenerator").getRealValue("nSample")) ;
  _tfoam->SetPseRan(RooRandom::randomGenerator()) ;
  _tfoam->SetChat((Int_t)config.getConfigSection("RooFoamGenerator").getRealValue("chatLevel")) ;
  _tfoam->Initialize() ;

  _vec = new Double_t[_realVars.getSize()] ;
  _xmin  = new Double_t[_realVars.getSize()] ;
  _range = new Double_t[_realVars.getSize()] ;
  
  Int_t i(0) ;
  for (const auto arg : _realVars) {
    auto var = static_cast<const RooRealVar*>(arg);
    _xmin[i] = var->getMin() ;
    _range[i] = var->getMax() - var->getMin() ;
    i++ ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooFoamGenerator::~RooFoamGenerator() 
{
  delete[] _vec ;
  delete[] _xmin ;
  delete[] _range ;
  delete _tfoam ;
  delete _binding ;
}



////////////////////////////////////////////////////////////////////////////////
/// are we actually generating anything? (the cache always contains at least our function value)

const RooArgSet *RooFoamGenerator::generateEvent(UInt_t /*remaining*/, Double_t& /*resampleRatio*/) 
{
  const RooArgSet *event= _cache->get();
  if(event->getSize() == 1) return event;

  _tfoam->MakeEvent() ;
  _tfoam->GetMCvect(_vec) ;
  
  // Transfer contents to dataset
  Int_t i(0) ;
  for (auto arg : _realVars) {
    auto var = static_cast<RooRealVar*>(arg);
    var->setVal(_xmin[i] + _range[i]*_vec[i]) ;
    i++ ;
  }
  return &_realVars ;
}
