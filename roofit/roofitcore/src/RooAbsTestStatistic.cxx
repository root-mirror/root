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
\file RooAbsTestStatistic.cxx
\class RooAbsTestStatistic
\ingroup Roofitcore

RooAbsTestStatistic is the abstract base class for all test
statistics. Test statistics that evaluate the PDF at each data
point should inherit from the RooAbsOptTestStatistic class which
implements several generic optimizations that can be done for such
quantities.

This test statistic base class organizes calculation of test
statistic values for RooSimultaneous PDF as a combination of test
statistic values for the PDF components of the simultaneous PDF and
organizes multi-processor parallel calculation of test statistic
values. For the latter, the test statistic value is calculated in
partitions in parallel executing processes and a posteriori
combined in the main thread.
**/

#include "RooAbsTestStatistic.h"

#include "RooFit.h"
#include "RooAbsPdf.h"
#include "RooSimultaneous.h"
#include "RooAbsData.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooNLLVar.h"
#include "RooRealMPFE.h"
#include "RooErrorHandler.h"
#include "RooMsgService.h"
#include "RooProdPdf.h"
#include "RooRealSumPdf.h"
#include "RooConstVar.h"
#include "RooRealIntegral.h"
#include "RooAbsCategoryLValue.h"

#include "TTimeStamp.h"
#include "TClass.h"
#include <string>
#include <fstream>
#include <sstream>

// timing
#include "RooTimer.h"
// getpid and getppid:
#include "unistd.h"
#include <stdexcept>

using namespace std;

ClassImp(RooAbsTestStatistic);


////////////////////////////////////////////////////////////////////////////////
/// Create a test statistic from the given function and the data.
/// \param[in] name Name of the test statistic
/// \param[in] title Title (for plotting)
/// \param[in] real Function to be used for tests
/// \param[in] data Data to fit function to
/// \param[in] projDeps A set of projected observables
/// \param[in] rangeName Fit data only in range with given name
/// \param[in] addCoefRangeName If not null, all RooAddPdf components of `real` will be instructed to fix their fraction definitions to the given named range.
/// \param[in] nCPU If larger than one, the test statistic calculation will be parallelized over multiple processes.
/// By default the data is split with 'bulk' partitioning (each process calculates a contigious block of fraction 1/nCPU
/// of the data). For binned data this approach may be suboptimal as the number of bins with >0 entries
/// in each processing block many vary greatly thereby distributing the workload rather unevenly.
/// \param[in] interleave is set to true, the interleave partitioning strategy is used where each partition
/// i takes all bins for which (ibin % ncpu == i) which is more likely to result in an even workload.
/// \param[in] verbose Be more verbose.
/// \param[in] splitCutRange If true, a different rangeName constructed as rangeName_{catName} will be used
/// as range definition for each index state of a RooSimultaneous. This means that a different range can be defined
/// for each category such as
/// ```
/// myVariable.setRange("range_pi0", 135, 210);
/// myVariable.setRange("range_gamma", 50, 210);
/// ```
/// if the categories are called "pi0" and "gamma".

RooAbsTestStatistic::RooAbsTestStatistic(const char *name, const char *title, RooAbsReal& real, RooAbsData& data,
                                         const RooArgSet& projDeps, RooAbsTestStatistic::Configuration const& cfg) :
  RooAbsReal(name,title),
  _paramSet("paramSet","Set of parameters",this),
  _func(&real),
  _data(&data),
  _projDeps((RooArgSet*)projDeps.Clone()),
  _rangeName(cfg.rangeName),
  _addCoefRangeName(cfg.addCoefRangeName),
  _splitRange(cfg.splitCutRange),
  _verbose(cfg.verbose),
  // Determine if RooAbsReal is a RooSimultaneous
  _gofOpMode{(cfg.nCPU>1 || cfg.nCPU==-1) ? MPMaster : (dynamic_cast<RooSimultaneous*>(_func) ? SimMaster : Slave)},
  _nEvents{data.numEntries()},
  _nCPU(cfg.nCPU != -1 ? cfg.nCPU : 1),
  _mpinterl(cfg.interleave),
  _CPUAffinity(cfg.CPUAffinity)
{
  // Register all parameters as servers
  _paramSet.add(*std::unique_ptr<RooArgSet>{real.getParameters(&data)});

  if (cfg.rangeName.find(',') != std::string::npos) {
    auto errorMsg = std::string("Ranges ") + cfg.rangeName
            + " were passed to the RooAbsTestStatistic with name \"" + name + "\", "
            + "but it doesn't support multiple comma-separated fit ranges!\n" +
            + "Instead, one should combine multiple RooAbsTestStatistic objects "
            + "(see RooAbsPdf::createNLL for an example with RooNLLVar).";
    coutE(InputArguments) <<  errorMsg << std::endl;
    throw std::invalid_argument(errorMsg);
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooAbsTestStatistic::RooAbsTestStatistic(const RooAbsTestStatistic& other, const char* name) :
  RooAbsReal(other,name),
  _paramSet("paramSet","Set of parameters",this),
  _func(other._func),
  _data(other._data),
  _projDeps((RooArgSet*)other._projDeps->Clone()),
  _rangeName(other._rangeName),
  _addCoefRangeName(other._addCoefRangeName),
  _splitRange(other._splitRange),
  _verbose(other._verbose),
  // Determine if RooAbsReal is a RooSimultaneous
  _gofOpMode{(other._nCPU>1 || other._nCPU==-1) ? MPMaster : (dynamic_cast<RooSimultaneous*>(_func) ? SimMaster : Slave)},
  _nEvents{_data->numEntries()},
  _gofSplitMode(other._gofSplitMode),
  _nCPU(other._nCPU != -1 ? other._nCPU : 1),
  _mpinterl(other._mpinterl),
  _CPUAffinity(other._CPUAffinity),
  _doOffset(other._doOffset),
  _offset(other._offset),
  _evalCarry(other._evalCarry)
{
  // Our parameters are those of original
  _paramSet.add(other._paramSet) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooAbsTestStatistic::~RooAbsTestStatistic()
{
  if (MPMaster == _gofOpMode && _init) {
    for (Int_t i = 0; i < _nCPU; ++i) delete _mpfeArray[i];
    delete[] _mpfeArray ;
  }

  if (SimMaster == _gofOpMode && _init) {
    for (Int_t i = 0; i < _nGof; ++i) delete _gofArray[i];
    delete[] _gofArray ;
  }

  delete _projDeps ;

}



////////////////////////////////////////////////////////////////////////////////
/// Calculate and return value of test statistic. If the test statistic
/// is calculated from a RooSimultaneous, the test statistic calculation
/// is performed separately on each simultaneous p.d.f component and associated
/// data, and then combined. If the test statistic calculation is parallelized,
/// partitions are calculated in nCPU processes and combined a posteriori.

Double_t RooAbsTestStatistic::evaluate() const
{
  // One-time Initialization
  if (!_init) {
    const_cast<RooAbsTestStatistic*>(this)->initialize() ;
  }

  if (SimMaster == _gofOpMode) {
    // Evaluate array of owned GOF objects
    Double_t ret = 0.;

    if (_mpinterl == RooFit::BulkPartition || _mpinterl == RooFit::Interleave ) {
      ret = combinedValue((RooAbsReal**)_gofArray,_nGof);
    } else {
      Double_t sum = 0., carry = 0.;
      for (Int_t i = 0 ; i < _nGof; ++i) {
	if (i % _numSets == _setNum || (_mpinterl==RooFit::Hybrid && _gofSplitMode[i] != RooFit::SimComponents )) {
	  Double_t y = _gofArray[i]->getValV();
	  carry += _gofArray[i]->getCarry();
	  y -= carry;
	  const Double_t t = sum + y;
	  carry = (t - sum) - y;
	  sum = t;
	}
      }
      ret = sum ;
      _evalCarry = carry;
    }

    // Only apply global normalization if SimMaster doesn't have MP master
    if (numSets()==1) {
      const Double_t norm = globalNormalization();
      ret /= norm;
      _evalCarry /= norm;
    }

    return ret ;

  } else if (MPMaster == _gofOpMode) {

    // Start calculations in parallel
    for (Int_t i = 0; i < _nCPU; ++i) _mpfeArray[i]->calculate();


    Double_t sum(0), carry = 0.;
    for (Int_t i = 0; i < _nCPU; ++i) {
      Double_t y = _mpfeArray[i]->getValV();
      carry += _mpfeArray[i]->getCarry();
      y -= carry;
      const Double_t t = sum + y;
      carry = (t - sum) - y;
      sum = t;
    }

    Double_t ret = sum ;
    _evalCarry = carry;

    const Double_t norm = globalNormalization();
    ret /= norm;
    _evalCarry /= norm;

    return ret ;

  } else {

    // Evaluate as straight FUNC
    Int_t nFirst(0), nLast(_nEvents), nStep(1) ;
    
    switch (_mpinterl) {
    case RooFit::BulkPartition:
      nFirst = _nEvents * _setNum / _numSets ;
      nLast  = _nEvents * (_setNum+1) / _numSets ;
      nStep  = 1 ;
      break;
      
    case RooFit::Interleave:
      nFirst = _setNum ;
      nLast  = _nEvents ;
      nStep  = _numSets ;
      break ;
      
    case RooFit::SimComponents:
      nFirst = 0 ;
      nLast  = _nEvents ;
      nStep  = 1 ;
      break ;
      
    case RooFit::Hybrid:
      throw std::logic_error("this should never happen");
      break ;
    }

    Double_t ret = evaluatePartition(nFirst,nLast,nStep);

    if (numSets()==1) {
      const Double_t norm = globalNormalization();
      ret /= norm;
      _evalCarry /= norm;
    }

    return ret ;

  }
}



////////////////////////////////////////////////////////////////////////////////
/// One-time initialization of the test statistic. Setup
/// infrastructure for simultaneous p.d.f processing and/or
/// parallelized processing if requested

Bool_t RooAbsTestStatistic::initialize()
{
  if (_init) return kFALSE;

  if (MPMaster == _gofOpMode) {
    initMPMode(_func,_data,_projDeps,_rangeName,_addCoefRangeName) ;
  } else if (SimMaster == _gofOpMode) {
    initSimMode((RooSimultaneous*)_func,_data,_projDeps,_rangeName,_addCoefRangeName) ;
  }
  _init = kTRUE;
  return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Forward server redirect calls to component test statistics

Bool_t RooAbsTestStatistic::redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t)
{
  if (SimMaster == _gofOpMode && _gofArray) {
    // Forward to slaves
    for (Int_t i = 0; i < _nGof; ++i) {
      if (_gofArray[i]) {
	_gofArray[i]->recursiveRedirectServers(newServerList,mustReplaceAll,nameChange);
      }
    }
  } else if (MPMaster == _gofOpMode&& _mpfeArray) {
    // Forward to slaves
    for (Int_t i = 0; i < _nCPU; ++i) {
      if (_mpfeArray[i]) {
	_mpfeArray[i]->recursiveRedirectServers(newServerList,mustReplaceAll,nameChange);
// 	cout << "redirecting servers on " << _mpfeArray[i]->GetName() << endl;
      }
    }
  }
  return kFALSE;
}



////////////////////////////////////////////////////////////////////////////////
/// Add extra information on component test statistics when printing
/// itself as part of a tree structure

void RooAbsTestStatistic::printCompactTreeHook(ostream& os, const char* indent)
{
  if (SimMaster == _gofOpMode) {
    // Forward to slaves
    os << indent << "RooAbsTestStatistic begin GOF contents" << endl ;
    for (Int_t i = 0; i < _nGof; ++i) {
      if (_gofArray[i]) {
	TString indent2(indent);
	indent2 += Form("[%d] ",i);
	_gofArray[i]->printCompactTreeHook(os,indent2);
      }
    }
    os << indent << "RooAbsTestStatistic end GOF contents" << endl;
  } else if (MPMaster == _gofOpMode) {
    // WVE implement this
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Forward constant term optimization management calls to component
/// test statistics

void RooAbsTestStatistic::constOptimizeTestStatistic(ConstOpCode opcode, Bool_t doAlsoTrackingOpt)
{
  initialize();
  if (SimMaster == _gofOpMode) {
    // Forward to slaves
    for (Int_t i = 0; i < _nGof; ++i) {
      // In SimComponents Splitting strategy only constOptimize the terms that are actually used
      RooFit::MPSplit effSplit = (_mpinterl!=RooFit::Hybrid) ? _mpinterl : _gofSplitMode[i];
      if ( (i % _numSets == _setNum) || (effSplit != RooFit::SimComponents) ) {
	if (_gofArray[i]) _gofArray[i]->constOptimizeTestStatistic(opcode,doAlsoTrackingOpt);
      }
    }
  } else if (MPMaster == _gofOpMode) {
    for (Int_t i = 0; i < _nCPU; ++i) {
      _mpfeArray[i]->constOptimizeTestStatistic(opcode,doAlsoTrackingOpt);
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Set MultiProcessor set number identification of this instance

void RooAbsTestStatistic::setMPSet(Int_t inSetNum, Int_t inNumSets)
{
  _setNum = inSetNum; _numSets = inNumSets;
  _extSet = _mpinterl==RooFit::SimComponents ? _setNum : (_numSets - 1);
  
  if (SimMaster == _gofOpMode) {
    // Forward to slaves
    initialize();
    for (Int_t i = 0; i < _nGof; ++i) {
      if (_gofArray[i]) _gofArray[i]->setMPSet(inSetNum,inNumSets);
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Initialize multi-processor calculation mode. Create component test statistics in separate
/// processed that are connected to this process through a RooAbsRealMPFE front-end class.

void RooAbsTestStatistic::initMPMode(RooAbsReal* real, RooAbsData* data, const RooArgSet* projDeps, std::string const& rangeName, std::string const& addCoefRangeName)
{
  _mpfeArray = new pRooRealMPFE[_nCPU];

  // Create proto-goodness-of-fit
  Configuration cfg;
  cfg.rangeName = rangeName;
  cfg.addCoefRangeName = addCoefRangeName;
  cfg.nCPU = 1;
  cfg.interleave = _mpinterl;
  cfg.CPUAffinity=_CPUAffinity;
  cfg.verbose = _verbose;
  cfg.splitCutRange = _splitRange;
  // This configuration parameter is stored in the RooAbsOptTestStatistic.
  // It would have been cleaner to move the member variable into RooAbsTestStatistic,
  // but to avoid incrementing the class version we do the dynamic_cast trick.
  if(auto thisAsRooAbsOptTestStatistic = dynamic_cast<RooAbsOptTestStatistic const*>(this)) {
    cfg.integrateOverBinsPrecision = thisAsRooAbsOptTestStatistic->_integrateBinsPrecision;
  }
  RooAbsTestStatistic* gof = create(GetName(),GetTitle(),*real,*data,*projDeps,cfg);
  gof->recursiveRedirectServers(_paramSet);

  for (Int_t i = 0; i < _nCPU; ++i) {
    gof->setMPSet(i,_nCPU);
    gof->SetName(Form("%s_GOF%d",GetName(),i));
    gof->SetTitle(Form("%s_GOF%d",GetTitle(),i));

    ccoutD(Eval) << "RooAbsTestStatistic::initMPMode: starting remote server process #" << i << endl;
    _mpfeArray[i] = new RooRealMPFE(Form("%s_%lx_MPFE%d",GetName(),(ULong_t)this,i),Form("%s_%lx_MPFE%d",GetTitle(),(ULong_t)this,i),*gof, i, _nCPU,false);
    //_mpfeArray[i]->setVerbose(kTRUE,kTRUE);
    _mpfeArray[i]->initialize();
    if (i > 0) {
      _mpfeArray[i]->followAsSlave(*_mpfeArray[0]);
    }

    if (_CPUAffinity == kTRUE) {
      _mpfeArray[i]->setCpuAffinity(i);
    }
  }
  _mpfeArray[_nCPU - 1]->addOwnedComponents(*gof);
  coutI(Eval) << "RooAbsTestStatistic::initMPMode: started " << _nCPU << " remote server process." << endl;
  return ;
}


////////////////////////////////////////////////////////////////////////////////
/// Activate timing of numerical integral normalization terms in the pdf.
/// This function should be called from the process that evaluates the pdf.
/// Using RooRealMPFE this means it should be called from the serverLoop().

void RooAbsTestStatistic::_setNumIntTimingInPdfs(Bool_t flag) {
  // find all pdf nodes with a normalization integral and set the activate timing flag on them
  // Get list of branch nodes in expression
  RooArgSet blist;
  RooAbsArg* node;

  _func->branchNodeServerList(&blist);

  // Iterator over branch nodes
  RooFIter iter = blist.fwdIterator();
  while((node = iter.next())) {
    RooAbsPdf *pdfNode = dynamic_cast<RooAbsPdf *>(node);
    if (!pdfNode) continue;
    // Skip self-normalized nodes
    if (pdfNode->selfNormalized()) continue;

    // TODO: move everything below to RooAbsPdf::setNumIntTiming(RooArgSet& obsSet, Bool_t flag) and only call that here? Or do we not want to pass around the observables?
    // set attribute on or off
    pdfNode->set_num_int_timing_flag(flag);

    // Retrieve normalization integral object for branch nodes that are pdfs
    RooRealIntegral *normint = const_cast<RooRealIntegral *>(dynamic_cast<const RooRealIntegral *>(pdfNode->getNormIntegral(*(_data->get()))));
    if (!normint) continue;

    normint->activateTimingNumInts();
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize simultaneous p.d.f processing mode. Strip simultaneous
/// p.d.f into individual components, split dataset in subset
/// matching each component and create component test statistics for
/// each of them.

void RooAbsTestStatistic::initSimMode(RooSimultaneous* simpdf, RooAbsData* data,
                                      const RooArgSet* projDeps,
                                      std::string const& rangeName, std::string const& addCoefRangeName)
{

  RooAbsCategoryLValue& simCat = const_cast<RooAbsCategoryLValue&>(simpdf->indexCat());

  TString simCatName(simCat.GetName());
  std::unique_ptr<TList> dsetList{const_cast<RooAbsData*>(data)->split(simCat, processEmptyDataSets())};
  if (!dsetList) {
    coutE(Fitting) << "RooAbsTestStatistic::initSimMode(" << GetName() << ") ERROR: index category of simultaneous pdf is missing in dataset, aborting" << endl;
    throw std::runtime_error("RooAbsTestStatistic::initSimMode() ERROR, index category of simultaneous pdf is missing in dataset, aborting");
  }

  // Count number of used states
  Int_t n = 0;
  _nGof = 0;

  for (const auto& catState : simCat) {
    // Retrieve the PDF for this simCat state
    RooAbsPdf* pdf = simpdf->getPdf(catState.first.c_str());
    RooAbsData* dset = (RooAbsData*) dsetList->FindObject(catState.first.c_str());

    if (pdf && dset && (0. != dset->sumEntries() || processEmptyDataSets())) {
      ++_nGof;
    }
  }

  // Allocate arrays
  _gofArray = new pRooAbsTestStatistic[_nGof];
  _gofSplitMode.resize(_nGof);

  // Create array of regular fit contexts, containing subset of data and single fitCat PDF
  for (const auto& catState : simCat) {
    const std::string& catName = catState.first;
    // Retrieve the PDF for this simCat state
    RooAbsPdf* pdf = simpdf->getPdf(catName.c_str());
    RooAbsData* dset = (RooAbsData*) dsetList->FindObject(catName.c_str());

    if (pdf && dset && (0. != dset->sumEntries() || processEmptyDataSets())) {
      ccoutI(Fitting) << "RooAbsTestStatistic::initSimMode: creating slave calculator #" << n << " for state " << catName
          << " (" << dset->numEntries() << " dataset entries)" << endl;


      // *** START HERE
      // WVE HACK determine if we have a RooRealSumPdf and then treat it like a binned likelihood
      RooAbsPdf* binnedPdf = 0 ;
      Bool_t binnedL = kFALSE ;
      if (pdf->getAttribute("BinnedLikelihood") && pdf->IsA()->InheritsFrom(RooRealSumPdf::Class())) {
        // Simplest case: top-level of component is a RRSP
        binnedPdf = pdf ;
        binnedL = kTRUE ;
      } else if (pdf->IsA()->InheritsFrom(RooProdPdf::Class())) {
        // Default case: top-level pdf is a product of RRSP and other pdfs
        RooFIter iter = ((RooProdPdf*)pdf)->pdfList().fwdIterator() ;
        RooAbsArg* component ;
        while ((component = iter.next())) {
          if (component->getAttribute("BinnedLikelihood") && component->IsA()->InheritsFrom(RooRealSumPdf::Class())) {
            binnedPdf = (RooAbsPdf*) component ;
            binnedL = kTRUE ;
          }
          if (component->getAttribute("MAIN_MEASUREMENT")) {
            // not really a binned pdf, but this prevents a (potentially) long list of subsidiary measurements to be passed to the slave calculator
            binnedPdf = (RooAbsPdf*) component ;
          }
        }
      }
      // WVE END HACK
      // Below here directly pass binnedPdf instead of PROD(binnedPdf,constraints) as constraints are evaluated elsewhere anyway
      // and omitting them reduces model complexity and associated handling/cloning times
      Configuration cfg;
      cfg.addCoefRangeName = addCoefRangeName;
      cfg.interleave = _mpinterl;
      cfg.CPUAffinity=_CPUAffinity;
      cfg.verbose = _verbose;
      cfg.splitCutRange = _splitRange;
      cfg.binnedL = binnedL;
      // This configuration parameter is stored in the RooAbsOptTestStatistic.
      // It would have been cleaner to move the member variable into RooAbsTestStatistic,
      // but to avoid incrementing the class version we do the dynamic_cast trick.
      if(auto thisAsRooAbsOptTestStatistic = dynamic_cast<RooAbsOptTestStatistic const*>(this)) {
        cfg.integrateOverBinsPrecision = thisAsRooAbsOptTestStatistic->_integrateBinsPrecision;
      }
      if (_splitRange && !rangeName.empty()) {
        cfg.rangeName = rangeName + "_" + catName;
        cfg.nCPU = _nCPU*(_mpinterl?-1:1);
      } else {
        cfg.rangeName = rangeName;
        cfg.nCPU = _nCPU;
      }
      _gofArray[n] = create(catName.c_str(),catName.c_str(),(binnedPdf?*binnedPdf:*pdf),*dset,*projDeps,cfg);
      _gofArray[n]->setSimCount(_nGof);
      // *** END HERE

      // Fill per-component split mode with Bulk Partition for now so that Auto will map to bulk-splitting of all components
      if (_mpinterl==RooFit::Hybrid) {
        if (dset->numEntries()<10) {
          //cout << "RAT::initSim("<< GetName() << ") MP mode is auto, setting split mode for component "<< n << " to SimComponents"<< endl ;
          _gofSplitMode[n] = RooFit::SimComponents;
          _gofArray[n]->_mpinterl = RooFit::SimComponents;
        } else {
          //cout << "RAT::initSim("<< GetName() << ") MP mode is auto, setting split mode for component "<< n << " to BulkPartition"<< endl ;
          _gofSplitMode[n] = RooFit::BulkPartition;
          _gofArray[n]->_mpinterl = RooFit::BulkPartition;
        }
      }

      // Servers may have been redirected between instantiation and (deferred) initialization

      std::unique_ptr<RooArgSet> actualParams{binnedPdf ? binnedPdf->getParameters(dset) : pdf->getParameters(dset)};
      std::unique_ptr<RooArgSet> selTargetParams{(RooArgSet*) _paramSet.selectCommon(*actualParams)};

      _gofArray[n]->recursiveRedirectServers(*selTargetParams);

       ++n;

    } else {
      if ((!dset || (0. != dset->sumEntries() && !processEmptyDataSets())) && pdf) {
        if (_verbose) {
          ccoutD(Fitting) << "RooAbsTestStatistic::initSimMode: state " << catName
              << " has no data entries, no slave calculator created" << endl;
        }
      }
    }
  }
  coutI(Fitting) << "RooAbsTestStatistic::initSimMode: created " << n << " slave calculators." << endl;
  
//  // Delete datasets by hand as TList::Delete() doesn't see our datasets as 'on the heap'...
//  TIterator* iter = dsetList->MakeIterator();
//  TObject* ds;
//  while((ds = iter->Next())) {
//    delete ds;
//  }
//  delete iter;
   // NOTE: commented out, because it actually causes two error messages, "Error in TList::Clear: A list is accessing an object (0x7f8320ae8f64) already deleted (list name = TList)" and another like that
}


////////////////////////////////////////////////////////////////////////////////
/// Change dataset that is used to given one. If cloneData is kTRUE, a clone of
/// in the input dataset is made.  If the test statistic was constructed with
/// a range specification on the data, the cloneData argument is ignored and
/// the data is always cloned.
Bool_t RooAbsTestStatistic::setData(RooAbsData& indata, Bool_t cloneData) 
{ 
  // Trigger refresh of likelihood offsets 
  if (isOffsetting()) {
    enableOffsetting(kFALSE);
    enableOffsetting(kTRUE);
  }

  switch(operMode()) {
  case Slave:
    // Delegate to implementation
    return setDataSlave(indata, cloneData);
  case SimMaster:
    // Forward to slaves
    //     cout << "RATS::setData(" << GetName() << ") SimMaster, calling setDataSlave() on slave nodes" << endl;
    if (indata.canSplitFast()) {
      for (Int_t i = 0; i < _nGof; ++i) {
	RooAbsData* compData = indata.getSimData(_gofArray[i]->GetName());
	_gofArray[i]->setDataSlave(*compData, cloneData);
      }
    } else if (0 == indata.numEntries()) {
      // For an unsplit empty dataset, simply assign empty dataset to each component
      for (Int_t i = 0; i < _nGof; ++i) {
	_gofArray[i]->setDataSlave(indata, cloneData);
      }
    } else {
      const RooAbsCategoryLValue& indexCat = static_cast<RooSimultaneous*>(_func)->indexCat();
      TList* dlist = indata.split(indexCat, kTRUE);
      if (!dlist) {
        coutF(DataHandling) << "Tried to split '" << indata.GetName() << "' into categories of '" << indexCat.GetName()
            << "', but splitting failed. Input data:" << std::endl;
        indata.Print("V");
        throw std::runtime_error("Error when setting up test statistic: dataset couldn't be split into categories.");
      }

      for (Int_t i = 0; i < _nGof; ++i) {
        RooAbsData* compData = (RooAbsData*) dlist->FindObject(_gofArray[i]->GetName());
        // 	cout << "component data for index " << _gofArray[i]->GetName() << " is " << compData << endl;
        if (compData) {
          _gofArray[i]->setDataSlave(*compData,kFALSE,kTRUE);
        } else {
          coutE(DataHandling) << "RooAbsTestStatistic::setData(" << GetName() << ") ERROR: Cannot find component data for state " << _gofArray[i]->GetName() << endl;
        }
      }
      delete dlist; // delete only list, data will be used
    }
    break;
  case MPMaster:
    // Not supported
    coutF(DataHandling) << "RooAbsTestStatistic::setData(" << GetName() << ") FATAL: setData() is not supported in multi-processor mode" << endl;
    throw std::runtime_error("RooAbsTestStatistic::setData is not supported in MPMaster mode");
    break;
  }

  return kTRUE;
}



void RooAbsTestStatistic::enableOffsetting(Bool_t flag) 
{
  // Apply internal value offsetting to control numeric precision
  if (!_init) {
    const_cast<RooAbsTestStatistic*>(this)->initialize() ;
  }
  
  switch(operMode()) {
  case Slave:
    _doOffset = flag ;
    // Clear offset if feature is disabled to that it is recalculated next time it is enabled
    if (!_doOffset) {
      _offset = 0 ;
    }
    setValueDirty() ;
    break ;
  case SimMaster:
    _doOffset = flag;
    for (Int_t i = 0; i < _nGof; ++i) {
      _gofArray[i]->enableOffsetting(flag);
    }
    break ;
  case MPMaster:    
    _doOffset = flag;
    for (Int_t i = 0; i < _nCPU; ++i) {
      _mpfeArray[i]->enableOffsetting(flag);
    }
    break;
  }
}


Double_t RooAbsTestStatistic::getCarry() const
{ return _evalCarry; }
