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
\file RooMCStudy.cxx
\class RooMCStudy
\ingroup Roofitcore

RooMCStudy is a helper class to facilitate Monte Carlo studies
such as 'goodness-of-fit' studies, that involve fitting a PDF 
to multiple toy Monte Carlo sets. These may be generated from either same PDF
or from a different PDF with similar parameters.

Given a fit and a generator PDF (they might be identical), RooMCStudy can produce
toyMC samples and/or fit these.
It accumulates the post-fit parameters of each iteration in a dataset. These can be
retrieved using fitParams() or fitParDataSet(). This dataset additionally contains the
variables
- NLL: The value of the negative log-likelihood for each run.
- ngen: The number of events generated for each run.

Additional plotting routines simplify the task of plotting
the distribution of the minimized likelihood, the fitted parameter values,
fitted error and pull distribution.

RooMCStudy provides the option to insert add-in modules
that modify the generate-and-fit cycle and allow to perform
extra steps in the cycle. Output of these modules can be stored
alongside the fit results in the aggregate results dataset.
These study modules should derive from the class RooAbsMCStudyModule.

Check the RooFit tutorials
- rf801_mcstudy.C
- rf802_mcstudy_addons.C
- rf803_mcstudy_addons2.C
- rf804_mcstudy_constr.C
for usage examples.
**/



#include "RooFit.h"
#include "snprintf.h"
#include <iostream>

#include "RooMCStudy.h"
#include "RooAbsMCStudyModule.h"

#include "RooGenContext.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooErrorVar.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooRandom.h"
#include "RooCmdConfig.h"
#include "RooGlobalFunc.h"
#include "RooPullVar.h"
#include "RooMsgService.h"
#include "RooProdPdf.h"

using namespace std ;

ClassImp(RooMCStudy);
  ;


/**
Construct Monte Carlo Study Manager. This class automates generating data from a given PDF,
fitting the PDF to data and accumulating the fit statistics.

\param[in] model The PDF to be studied
\param[in] observables The variables of the PDF to be considered observables
\param[in] argX Arguments from the table below

<table>
<tr><th> Optional arguments <th>
<tr><td> Silence()                         <td> Suppress all RooFit messages during running below PROGRESS level
<tr><td> FitModel(const RooAbsPdf&)        <td> The PDF for fitting if it is different from the PDF for generating.
<tr><td> ConditionalObservables(const RooArgSet& set)  <td> The set of observables that the PDF should _not_ be normalized over
<tr><td> Binned(Bool_t flag)               <td> Bin the dataset before fitting it. Speeds up fitting of large data samples
<tr><td> FitOptions(const char*)           <td> Classic fit options, provided for backward compatibility
<tr><td> FitOptions(....)                  <td> Options to be used for fitting. All named arguments inside FitOptions() are passed to RooAbsPdf::fitTo().
                                                `Save()` is especially interesting to be able to retrieve fit results of each run using fitResult().
<tr><td> Verbose(Bool_t flag)              <td> Activate informational messages in event generation phase
<tr><td> Extended(Bool_t flag)             <td> Determine number of events for each sample anew from a Poisson distribution
<tr><td> Constrain(const RooArgSet& pars)  <td> Apply internal constraints on given parameters in fit and sample constrained parameter values from constraint p.d.f for each toy.
<tr><td> ExternalConstraints(const RooArgSet& ) <td> Apply internal constraints on given parameters in fit and sample constrained parameter values from constraint p.d.f for each toy.
<tr><td> ProtoData(const RooDataSet&, Bool_t randOrder)
         <td> Prototype data for the event generation. If the randOrder flag is set, the order of the dataset will be re-randomized for each generation
              cycle to protect against systematic biases if the number of generated events does not exactly match the number of events in the prototype dataset
              at the cost of reduced precision with mu equal to the specified number of events
</table>
*/
RooMCStudy::RooMCStudy(const RooAbsPdf& model, const RooArgSet& observables,
   		       const RooCmdArg& arg1, const RooCmdArg& arg2,
   		       const RooCmdArg& arg3,const RooCmdArg& arg4,const RooCmdArg& arg5,
   		       const RooCmdArg& arg6,const RooCmdArg& arg7,const RooCmdArg& arg8) : TNamed("mcstudy","mcstudy")

{
  // Stuff all arguments in a list
  RooLinkedList cmdList;
  cmdList.Add(const_cast<RooCmdArg*>(&arg1)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg2)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg3)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg4)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg5)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg6)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg7)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg8)) ;

  // Select the pdf-specific commands 
  RooCmdConfig pc(Form("RooMCStudy::RooMCStudy(%s)",model.GetName())) ;
  
  pc.defineObject("fitModel","FitModel",0,0) ;
  pc.defineObject("condObs","ProjectedDependents",0,0) ;
  pc.defineObject("protoData","PrototypeData",0,0) ;
  pc.defineSet("cPars","Constrain",0,0) ;
  pc.defineSet("extCons","ExternalConstraints",0,0) ;
  pc.defineInt("silence","Silence",0,0) ;
  pc.defineInt("randProtoData","PrototypeData",0,0) ;
  pc.defineInt("verboseGen","Verbose",0,0) ;
  pc.defineInt("extendedGen","Extended",0,0) ;
  pc.defineInt("binGenData","Binned",0,0) ;
  pc.defineString("fitOpts","FitOptions",0,"") ;
  pc.defineInt("dummy","FitOptArgs",0,0) ;
  pc.defineMutex("FitOptions","FitOptArgs") ; // can have either classic or new-style fit options
  pc.defineMutex("Constrain","FitOptions") ; // constraints only work with new-style fit options
  pc.defineMutex("ExternalConstraints","FitOptions") ; // constraints only work with new-style fit options
  
  // Process and check varargs 
  pc.process(cmdList) ;
  if (!pc.ok(kTRUE)) {
    // WVE do something here
    throw std::string("RooMCStudy::RooMCStudy() Error in parsing arguments passed to contructor") ;
    return ;
  }

  // Save fit command options
  if (pc.hasProcessed("FitOptArgs")) {
    RooCmdArg* fitOptArg = static_cast<RooCmdArg*>(cmdList.FindObject("FitOptArgs")) ;
    for (Int_t i=0 ; i<fitOptArg->subArgs().GetSize() ;i++) {
      _fitOptList.Add(new RooCmdArg(static_cast<RooCmdArg&>(*fitOptArg->subArgs().At(i)))) ;
    }
  }

  // Decode command line arguments
  _silence = pc.getInt("silence") ;
  _verboseGen = pc.getInt("verboseGen") ;
  _extendedGen = pc.getInt("extendedGen") ;
  _binGenData = pc.getInt("binGenData") ;
  _randProto = pc.getInt("randProtoData") ;

  // Process constraints specifications
  const RooArgSet* cParsTmp = pc.getSet("cPars") ;
  const RooArgSet* extCons = pc.getSet("extCons") ;

  RooArgSet* cPars = new RooArgSet ;
  if (cParsTmp) {
    cPars->add(*cParsTmp) ;
  }
  
  // If constraints are specified, add to fit options
  if (cPars) {
    _fitOptList.Add(RooFit::Constrain(*cPars).Clone()) ;
  }
  if (extCons) {
    _fitOptList.Add(RooFit::ExternalConstraints(*extCons).Clone()) ;
  }

  // Make list of all constraints
  RooArgSet allConstraints ;
  RooArgSet consPars ;
  if (cPars) {
    RooArgSet* constraints = model.getAllConstraints(observables,*cPars,kTRUE) ;
    if (constraints) {
      allConstraints.add(*constraints) ;
      delete constraints ;
    }
  }
  
  // Construct constraint p.d.f
  if (allConstraints.getSize()>0) {
    _constrPdf = new RooProdPdf("mcs_constr_prod","RooMCStudy constraints product",allConstraints) ;

    if (cPars) {
      consPars.add(*cPars) ;
    } else {
      RooArgSet* params = model.getParameters(observables) ;
      RooArgSet* cparams = _constrPdf->getObservables(*params) ;
      consPars.add(*cparams) ;
      delete params ;
      delete cparams ;
    }
    _constrGenContext = _constrPdf->genContext(consPars,0,0,_verboseGen) ;

    _perExptGenParams = kTRUE ;

    coutI(Generation) << "RooMCStudy::RooMCStudy: INFO have pdf with constraints, will generate parameters from constraint pdf for each experiment" << endl ;


  } else {
    _constrPdf = 0 ;
    _constrGenContext=0 ;

    _perExptGenParams = kFALSE ;
  }

  
  // Extract generator and fit models
  _genModel = const_cast<RooAbsPdf*>(&model) ;
  _genSample = 0 ;
  RooAbsPdf* fitModel = static_cast<RooAbsPdf*>(pc.getObject("fitModel",0)) ;
  _fitModel = fitModel ? fitModel : _genModel ;
  
  // Extract conditional observables and prototype data
  _genProtoData = static_cast<RooDataSet*>(pc.getObject("protoData",0)) ;
  if (pc.getObject("condObs",0)) {
    _projDeps.add(static_cast<RooArgSet&>(*pc.getObject("condObs",0))) ;
  }
  
  _dependents.add(observables) ;
     
  _allDependents.add(_dependents) ;  
  _fitOptions = pc.getString("fitOpts") ;
  _canAddFitResults = kTRUE ;
  
  if (_extendedGen && _genProtoData && !_randProto) {
    oocoutW(_fitModel,Generation) << "RooMCStudy::RooMCStudy: WARNING Using generator option 'e' (Poisson distribution of #events) together " << endl
				  << "                        with a prototype dataset implies incomplete sampling or oversampling of proto data." << endl
				  << "                        Use option \"r\" to randomize prototype dataset order and thus to randomize" << endl
				  << "                        the set of over/undersampled prototype events for each generation cycle." << endl ;
  }
  
  _genParams = _genModel->getParameters(&_dependents) ;
  if (!_binGenData) {
    _genContext = _genModel->genContext(_dependents,_genProtoData,0,_verboseGen) ;
    _genContext->attach(*_genParams) ;
  } else {
    _genContext = 0 ;
  }

  _genInitParams = (RooArgSet*) _genParams->snapshot(kFALSE) ;

  // Store list of parameters and save initial values separately
  _fitParams = _fitModel->getParameters(&_dependents) ;
  _fitInitParams = (RooArgSet*) _fitParams->snapshot(kTRUE) ;
  
  _nExpGen = _extendedGen ? _genModel->expectedEvents(&_dependents) : 0 ;
  
  // Place holder for NLL
  _nllVar = new RooRealVar("NLL","-log(Likelihood)",0) ;

  // Place holder for number of generated events
  _ngenVar = new RooRealVar("ngen","number of generated events",0) ;
  
  // Create data set containing parameter values, errors and pulls
  RooArgSet tmp2(*_fitParams) ;
  tmp2.add(*_nllVar) ;
  tmp2.add(*_ngenVar) ;

  // Mark all variable to store their errors in the dataset
  tmp2.setAttribAll("StoreError",kTRUE) ;
  tmp2.setAttribAll("StoreAsymError",kTRUE) ;
  TString fpdName ;
  if (_fitModel==_genModel) {
    fpdName = Form("fitParData_%s",_fitModel->GetName()) ;
  } else {
    fpdName= Form("fitParData_%s_%s",_fitModel->GetName(),_genModel->GetName()) ;
  }

  _fitParData = new RooDataSet(fpdName.Data(),"Fit Parameters DataSet",tmp2) ;
  tmp2.setAttribAll("StoreError",kFALSE) ;
  tmp2.setAttribAll("StoreAsymError",kFALSE) ;

  if (_perExptGenParams) {
    _genParData = new RooDataSet("genParData","Generated Parameters dataset",*_genParams) ;
  } else {
    _genParData = 0 ;
  }
  
  // Append proto variables to allDependents
  if (_genProtoData) {
    _allDependents.add(*_genProtoData->get(),kTRUE) ;
  }

  // Call module initializers
  list<RooAbsMCStudyModule*>::iterator iter ;
  for (iter=_modList.begin() ; iter!= _modList.end() ; ++iter) {
    Bool_t ok = (*iter)->doInitializeInstance(*this) ;
    if (!ok) {
      oocoutE(_fitModel,Generation) << "RooMCStudy::ctor: removing study module " << (*iter)->GetName() << " from analysis chain because initialization failed" << endl ;
      iter = _modList.erase(iter) ;
    }
  }
  
}


////////////////////////////////////////////////////////////////////////////////
/// \deprecated PLEASE USE CONSTRUCTOR WITH NAMED ARGUMENTS. RETAINED FOR BACKWARD COMPATIBILY.
///
/// Constructor with a generator and fit model. Both models may point
/// to the same object. The 'dependents' set of variables is generated 
/// in the generator phase. The optional prototype dataset is passed to
/// the generator
///
/// Available generator options
///  v  - Verbose
///  e  - Extended: use Poisson distribution for Nevts generated
///
/// Available fit options
///  See RooAbsPdf::fitTo()
///

RooMCStudy::RooMCStudy(const RooAbsPdf& genModel, const RooAbsPdf& fitModel, 
   		       const RooArgSet& dependents, const char* genOptions, 
   		       const char* fitOptions, const RooDataSet* genProtoData, 
   		       const RooArgSet& projDeps) :
  TNamed("mcstudy","mcstudy"),
  _genModel((RooAbsPdf*)&genModel), 
  _genProtoData(genProtoData),
  _projDeps(projDeps),
  _constrPdf(0),
  _constrGenContext(0),
  _dependents(dependents), 
  _allDependents(dependents), 
  _fitModel((RooAbsPdf*)&fitModel), 
  _nllVar(0),
  _ngenVar(0),
  _genParData(0),
  _fitOptions(fitOptions),
  _canAddFitResults(kTRUE),
  _perExptGenParams(0),
  _silence(kFALSE)
{
  // Decode generator options
  TString genOpt(genOptions) ;
  genOpt.ToLower() ;
  _verboseGen = genOpt.Contains("v") ;
  _extendedGen = genOpt.Contains("e") ;
  _binGenData = genOpt.Contains("b") ;
  _randProto = genOpt.Contains("r") ;
  
  if (_extendedGen && genProtoData && !_randProto) {
    oocoutE(_fitModel,Generation) << "RooMCStudy::RooMCStudy: WARNING Using generator option 'e' (Poisson distribution of #events) together " << endl
				  << "                        with a prototype dataset implies incomplete sampling or oversampling of proto data." << endl
				  << "                        Use option \"r\" to randomize prototype dataset order and thus to randomize" << endl
				  << "                        the set of over/undersampled prototype events for each generation cycle." << endl ;
  }
  
  if (!_binGenData) {
    _genContext = genModel.genContext(dependents,genProtoData,0,_verboseGen) ;
  } else {
    _genContext = 0 ;
  }
  _genParams = _genModel->getParameters(&_dependents) ;
  _genSample = 0 ;
  RooArgSet* tmp = genModel.getParameters(&dependents) ;
  _genInitParams = (RooArgSet*) tmp->snapshot(kFALSE) ;
  delete tmp ;
  
  // Store list of parameters and save initial values separately
  _fitParams = fitModel.getParameters(&dependents) ;
  _fitInitParams = (RooArgSet*) _fitParams->snapshot(kTRUE) ;
  
  _nExpGen = _extendedGen ? genModel.expectedEvents(&dependents) : 0 ;
  
  // Place holder for NLL
  _nllVar = new RooRealVar("NLL","-log(Likelihood)",0) ;
  
  // Place holder for number of generated events
  _ngenVar = new RooRealVar("ngen","number of generated events",0) ;
  
  // Create data set containing parameter values, errors and pulls
  RooArgSet tmp2(*_fitParams) ;
  tmp2.add(*_nllVar) ;
  tmp2.add(*_ngenVar) ;
  
  // Mark all variable to store their errors in the dataset
  tmp2.setAttribAll("StoreError",kTRUE) ;
  tmp2.setAttribAll("StoreAsymError",kTRUE) ;
  _fitParData = new RooDataSet("fitParData","Fit Parameters DataSet",tmp2) ;
  tmp2.setAttribAll("StoreError",kFALSE) ;
  tmp2.setAttribAll("StoreAsymError",kFALSE) ;
  
  // Append proto variables to allDependents
  if (genProtoData) {
    _allDependents.add(*genProtoData->get(),kTRUE) ;
  }

  // Call module initializers
  list<RooAbsMCStudyModule*>::iterator iter ;
  for (iter=_modList.begin() ; iter!= _modList.end() ; ++iter) {
    Bool_t ok = (*iter)->doInitializeInstance(*this) ;
    if (!ok) {
      oocoutE(_fitModel,Generation) << "RooMCStudy::ctor: removing study module " << (*iter)->GetName() << " from analysis chain because initialization failed" << endl ;
      iter = _modList.erase(iter) ;
    }
  }
  
}



////////////////////////////////////////////////////////////////////////////////

RooMCStudy::~RooMCStudy() 
{  
  _genDataList.Delete() ;
  _fitOptList.Delete() ;
  delete _ngenVar ;
  delete _fitParData ;
  delete _genParData ;
  delete _fitInitParams ;
  delete _fitParams ;
  delete _genInitParams ;
  delete _genParams ;
  delete _genContext ;
  delete _nllVar ;
  delete _constrPdf ;
  delete _constrGenContext ;
}



////////////////////////////////////////////////////////////////////////////////
/// Insert given RooMCStudy add-on module to the processing chain
/// of this MCStudy object

void RooMCStudy::addModule(RooAbsMCStudyModule& module) 
{
  module.doInitializeInstance(*this) ;
  _modList.push_back(&module) ;        
}



////////////////////////////////////////////////////////////////////////////////
/// Run engine method. Generate and/or fit, according to flags, 'nSamples' samples of 'nEvtPerSample' events.
/// If keepGenData is set, all generated data sets will be kept in memory and can be accessed
/// later via genData().
///
/// When generating, data sets will be written out in ascii form if the pattern string is supplied
/// The pattern, which is a template for snprintf, should look something like "data/toymc_%04d.dat"
/// and should contain one integer field that encodes the sample serial number.
///
/// When fitting only, data sets may optionally be read from ascii files, using the same file
/// pattern.
///

Bool_t RooMCStudy::run(Bool_t doGenerate, Bool_t DoFit, Int_t nSamples, Int_t nEvtPerSample, Bool_t keepGenData, const char* asciiFilePat) 
{
  RooFit::MsgLevel oldLevel(RooFit::FATAL) ;
  if (_silence) {
    oldLevel = RooMsgService::instance().globalKillBelow() ;
    RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS) ;
  }

  list<RooAbsMCStudyModule*>::iterator iter ;
  for (iter=_modList.begin() ; iter!= _modList.end() ; ++iter) {
    (*iter)->initializeRun(nSamples) ;
  }  
  
  Int_t prescale = nSamples>100 ? Int_t(nSamples/100) : 1 ;

  while(nSamples--) {
    
    if (nSamples%prescale==0) {
      oocoutP(_fitModel,Generation) << "RooMCStudy::run: " ;
      if (doGenerate) ooccoutI(_fitModel,Generation) << "Generating " ;
      if (doGenerate && DoFit) ooccoutI(_fitModel,Generation) << "and " ;
      if (DoFit) ooccoutI(_fitModel,Generation) << "fitting " ;
      ooccoutP(_fitModel,Generation) << "sample " << nSamples << endl ;
    }

    _genSample = 0;
    Bool_t existingData = kFALSE ;
    if (doGenerate) {
      // Generate sample
      Int_t nEvt(nEvtPerSample) ;

      // Reset generator parameters to initial values
      *_genParams = *_genInitParams ;

      // If constraints are present, sample generator values from constraints
      if (_constrPdf) {
	RooDataSet* tmp = _constrGenContext->generate(1) ;
	*_genParams = *tmp->get() ;
	delete tmp ;
      }

      // Save generated parameters if required
      if (_genParData) {
	_genParData->add(*_genParams) ;
      }

      // Call module before-generation hook
      list<RooAbsMCStudyModule*>::iterator iter2 ;
      for (iter2=_modList.begin() ; iter2!= _modList.end() ; ++iter2) {
	(*iter2)->processBeforeGen(nSamples) ;
      }  

      if (_binGenData) {

	// Calculate the number of (extended) events for this run
	if (_extendedGen) {
	  _nExpGen = _genModel->expectedEvents(&_dependents) ;
	  nEvt = RooRandom::randomGenerator()->Poisson(nEvtPerSample==0?_nExpGen:nEvtPerSample) ;
	}	

	// Binned generation
	_genSample = _genModel->generateBinned(_dependents,nEvt) ;

      } else {

	// Calculate the number of (extended) events for this run
	if (_extendedGen) {
	  _nExpGen = _genModel->expectedEvents(&_dependents) ;
	  nEvt = RooRandom::randomGenerator()->Poisson(nEvtPerSample==0?_nExpGen:nEvtPerSample) ;
	}
	
	// Optional randomization of protodata for this run
	if (_randProto && _genProtoData && _genProtoData->numEntries()!=nEvt) {
	  oocoutI(_fitModel,Generation) << "RooMCStudy: (Re)randomizing event order in prototype dataset (Nevt=" << nEvt << ")" << endl ;
	  Int_t* newOrder = _genModel->randomizeProtoOrder(_genProtoData->numEntries(),nEvt) ;
	  _genContext->setProtoDataOrder(newOrder) ;
	  delete[] newOrder ;
	}

	coutP(Generation) << "RooMCStudy: now generating " << nEvt << " events" << endl ;
	
	// Actual generation of events
	if (nEvt>0) {
	  _genSample = _genContext->generate(nEvt) ;
	} else {
	  // Make empty dataset
	  _genSample = new RooDataSet("emptySample","emptySample",_dependents) ;
	}	
      } 

	
    //} else if (asciiFilePat && &asciiFilePat) { //warning: the address of 'asciiFilePat' will always evaluate as 'true'
    } else if (asciiFilePat) {

      // Load sample from ASCII file
      char asciiFile[1024] ;
      snprintf(asciiFile,1024,asciiFilePat,nSamples) ;
      RooArgList depList(_allDependents) ;
      _genSample = RooDataSet::read(asciiFile,depList,"q") ;      
      
    } else {
      
      // Load sample from internal list
      _genSample = (RooDataSet*) _genDataList.At(nSamples) ;
      existingData = kTRUE ;
      if (!_genSample) {
   	oocoutW(_fitModel,Generation) << "RooMCStudy::run: WARNING: Sample #" << nSamples << " not loaded, skipping" << endl ;
   	continue ;
      }
    }

    // Save number of generated events
    _ngenVar->setVal(_genSample->sumEntries()) ;

    // Call module between generation and fitting hook
    list<RooAbsMCStudyModule*>::iterator iter3 ;
    for (iter3=_modList.begin() ; iter3!= _modList.end() ; ++iter3) {
      (*iter3)->processBetweenGenAndFit(nSamples) ;
    }  
    
    if (DoFit) fitSample(_genSample) ;

    // Call module between generation and fitting hook
    for (iter3=_modList.begin() ; iter3!= _modList.end() ; ++iter3) {
      (*iter3)->processAfterFit(nSamples) ;
    }  
    
    // Optionally write to ascii file
    if (doGenerate && asciiFilePat && *asciiFilePat) {
      char asciiFile[1024] ;
      snprintf(asciiFile,1024,asciiFilePat,nSamples) ;
      RooDataSet* unbinnedData = dynamic_cast<RooDataSet*>(_genSample) ;
      if (unbinnedData) {
	unbinnedData->write(asciiFile) ;
      } else {
	coutE(InputArguments) << "RooMCStudy::run(" << GetName() << ") ERROR: ASCII writing of binned datasets is not supported" << endl ;
      }
    }
    
    // Add to list or delete
    if (!existingData) {
      if (keepGenData) {
	_genDataList.Add(_genSample) ;
      } else {
	delete _genSample ;
      }
    }
  }

  for (iter=_modList.begin() ; iter!= _modList.end() ; ++iter) {
    RooDataSet* auxData = (*iter)->finalizeRun() ;
    if (auxData) {
      _fitParData->merge(auxData) ;
    }
  }  

  _canAddFitResults = kFALSE ;

  if (_genParData) {
    const RooArgSet* genPars = _genParData->get() ;
    TIterator* iter2 = genPars->createIterator() ;
    RooAbsArg* arg ;
    while((arg=(RooAbsArg*)iter2->Next())) {
      _genParData->changeObservableName(arg->GetName(),Form("%s_gen",arg->GetName())) ;
    }
    delete iter2 ;
    
    _fitParData->merge(_genParData) ;
  }

  if (DoFit) calcPulls() ;

  if (_silence) {
    RooMsgService::instance().setGlobalKillBelow(oldLevel) ;
  }

  return kFALSE ;
}






////////////////////////////////////////////////////////////////////////////////
/// Generate and fit 'nSamples' samples of 'nEvtPerSample' events.
/// If keepGenData is set, all generated data sets will be kept in memory and can be accessed
/// later via genData().
///
/// Data sets will be written out in ascii form if the pattern string is supplied.
/// The pattern, which is a template for snprintf, should look something like "data/toymc_%04d.dat"
/// and should contain one integer field that encodes the sample serial number.
///

Bool_t RooMCStudy::generateAndFit(Int_t nSamples, Int_t nEvtPerSample, Bool_t keepGenData, const char* asciiFilePat) 
{
  // Clear any previous data in memory
  _fitResList.Delete() ; // even though the fit results are owned by gROOT, we still want to scratch them here.
  _genDataList.Delete() ;
  _fitParData->reset() ;
  
  return run(kTRUE,kTRUE,nSamples,nEvtPerSample,keepGenData,asciiFilePat) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Generate 'nSamples' samples of 'nEvtPerSample' events.
/// If keepGenData is set, all generated data sets will be kept in memory 
/// and can be accessed later via genData().
///
/// Data sets will be written out in ascii form if the pattern string is supplied.
/// The pattern, which is a template for snprintf, should look something like "data/toymc_%04d.dat"
/// and should contain one integer field that encodes the sample serial number.
///

Bool_t RooMCStudy::generate(Int_t nSamples, Int_t nEvtPerSample, Bool_t keepGenData, const char* asciiFilePat) 
{
  // Clear any previous data in memory
  _genDataList.Delete() ;
  
  return run(kTRUE,kFALSE,nSamples,nEvtPerSample,keepGenData,asciiFilePat) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Fit 'nSamples' datasets, which are read from ASCII files.
///
/// The ascii file pattern, which is a template for snprintf, should look something like "data/toymc_%04d.dat"
/// and should contain one integer field that encodes the sample serial number.
///

Bool_t RooMCStudy::fit(Int_t nSamples, const char* asciiFilePat) 
{
  // Clear any previous data in memory
  _fitResList.Delete() ; // even though the fit results are owned by gROOT, we still want to scratch them here.
  _fitParData->reset() ;
  
  return run(kFALSE,kTRUE,nSamples,0,kFALSE,asciiFilePat) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Fit 'nSamples' datasets, as supplied in 'dataSetList'
/// 

Bool_t RooMCStudy::fit(Int_t nSamples, TList& dataSetList) 
{
  // Clear any previous data in memory
  _fitResList.Delete() ; // even though the fit results are owned by gROOT, we still want to scratch them here.
  _genDataList.Delete() ;
  _fitParData->reset() ;
  
  // Load list of data sets
  TIterator* iter = dataSetList.MakeIterator() ;
  RooAbsData* gset ;
  while((gset=(RooAbsData*)iter->Next())) {
    _genDataList.Add(gset) ;
  }
  delete iter ;
  
  return run(kFALSE,kTRUE,nSamples,0,kTRUE,0) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Reset all fit parameters to the initial model
/// parameters at the time of the RooMCStudy constructor

void RooMCStudy::resetFitParams()
{
  *_fitParams = *_fitInitParams ;
}



////////////////////////////////////////////////////////////////////////////////
/// Internal function. Performs actual fit according to specifications

RooFitResult* RooMCStudy::doFit(RooAbsData* genSample)
{
  // Fit model to data set
  TString fitOpt2(_fitOptions) ; fitOpt2.Append("r") ;
  if (_silence) {
    fitOpt2.Append("b") ;
  }
  
  // Optionally bin dataset before fitting
  RooAbsData* data ;
  if (_binGenData) {    
    RooArgSet* depList = _fitModel->getObservables(genSample) ;
    data = new RooDataHist(genSample->GetName(),genSample->GetTitle(),*depList,*genSample) ;
    delete depList ;
  } else {
    data = genSample ;
  }
  
  RooFitResult* fr ;
  if (_fitOptList.GetSize()==0) {
    if (_projDeps.getSize()>0) {
      fr = (RooFitResult*) _fitModel->fitTo(*data,RooFit::ConditionalObservables(_projDeps),RooFit::FitOptions(fitOpt2)) ;
    } else {
      fr = (RooFitResult*) _fitModel->fitTo(*data,RooFit::FitOptions(fitOpt2)) ;
    }
  } else {
    RooCmdArg save  = RooFit::Save() ;
    RooCmdArg condo = RooFit::ConditionalObservables(_projDeps) ;
    RooCmdArg plevel = RooFit::PrintLevel(-1) ;
    RooLinkedList fitOptList(_fitOptList) ;
    fitOptList.Add(&save) ;
    if (_projDeps.getSize()>0) {
      fitOptList.Add(&condo) ;
    }
    if (_silence) {
      fitOptList.Add(&plevel) ;
    }
    fr = (RooFitResult*) _fitModel->fitTo(*data,fitOptList) ;
  }

  if (_binGenData) delete data ;

  return fr ;
}



////////////////////////////////////////////////////////////////////////////////
/// Redo fit on 'current' toy sample, or if genSample is not NULL
/// do fit on given sample instead

RooFitResult* RooMCStudy::refit(RooAbsData* genSample) 
{
  if (!genSample) {
    genSample = _genSample ;
  }

  RooFitResult* fr(0) ;
  if (genSample->sumEntries()>0) {
    fr = doFit(genSample) ;
  }

  return fr ;
}



////////////////////////////////////////////////////////////////////////////////
/// Internal method. Fit given dataset with fit model. If fit
/// converges (TMinuit status code zero) The fit results are appended
/// to the fit results dataset
///
/// If the fit option "r" is supplied, the RooFitResult
/// objects will always be saved, regardless of the
/// fit status. RooFitResults objects can be retrieved
/// later via fitResult().
///  

Bool_t RooMCStudy::fitSample(RooAbsData* genSample) 
{  
  // Reset all fit parameters to their initial values  
  resetFitParams() ;

  // Perform actual fit
  Bool_t ok ;
  RooFitResult* fr(0) ;
  if (genSample->sumEntries()>0) {
    fr = doFit(genSample) ;
    ok = (fr->status()==0) ;
  } else {
    ok = kFALSE ;
  }

  // If fit converged, store parameters and NLL
  if (ok) {
    _nllVar->setVal(fr->minNll()) ;
    RooArgSet tmp(*_fitParams) ;
    tmp.add(*_nllVar) ;
    tmp.add(*_ngenVar) ;

    _fitParData->add(tmp) ;
  }
  
  // Store fit result if requested by user
  Bool_t userSaveRequest = kFALSE ;
  if (_fitOptList.GetSize()>0) {
    if (_fitOptList.FindObject("Save")) userSaveRequest = kTRUE ;
  } else {
    if (_fitOptions.Contains("r")) userSaveRequest = kTRUE ;
  }

  if (userSaveRequest) {
    _fitResList.Add(fr) ;
  } else {
    delete fr ;
  }
    
  return !ok ;
}



////////////////////////////////////////////////////////////////////////////////
/// Utility function to add fit result from external fit to this RooMCStudy
/// and process its results through the standard RooMCStudy statistics gathering tools.
/// This function allows users to run the toy MC generation and/or fitting
/// in a distributed way and to collect and analyze the results in a RooMCStudy
/// as if they were run locally.
///
/// This method is only functional if this RooMCStudy object is cleanm, i.e. it was not used
/// to generate and/or fit any samples.

Bool_t RooMCStudy::addFitResult(const RooFitResult& fr) 
{  
  if (!_canAddFitResults) {
    oocoutE(_fitModel,InputArguments) << "RooMCStudy::addFitResult: ERROR cannot add fit results in current state" << endl ;
    return kTRUE ;
  }
  
  // Transfer contents of fit result to fitParams ;
  *_fitParams = RooArgSet(fr.floatParsFinal()) ;
  
  // If fit converged, store parameters and NLL
  Bool_t ok = (fr.status()==0) ;
  if (ok) {
    _nllVar->setVal(fr.minNll()) ;
    RooArgSet tmp(*_fitParams) ;
    tmp.add(*_nllVar) ;
    tmp.add(*_ngenVar) ;
    _fitParData->add(tmp) ;
  }
  
  // Store fit result if requested by user
  if (_fitOptions.Contains("r")) {
    _fitResList.Add((TObject*)&fr) ;
  }  
  
  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate the pulls for all fit parameters in
/// the fit results data set, and add them to that dataset.

void RooMCStudy::calcPulls() 
{
  for (auto it = _fitParams->begin(); it != _fitParams->end(); ++it) {
    const auto par = static_cast<RooRealVar*>(*it);
    RooErrorVar* err = par->errorVar();
    _fitParData->addColumn(*err);
    delete err;

    TString name(par->GetName()), title(par->GetTitle()) ;
    name.Append("pull") ;
    title.Append(" Pull") ;

    if (!par->hasError(false)) {
      coutW(Generation) << "Fit parameter '" << par->GetName() << "' does not have an error."
          " A pull distribution cannot be generated. This might be caused by the parameter being constant or"
          " because the fits were not run." << std::endl;
      continue;
    }

    // First look in fitParDataset to see if per-experiment generated value has been stored
    auto genParOrig = static_cast<RooAbsReal*>(_fitParData->get()->find(Form("%s_gen",par->GetName())));
    if (genParOrig && _perExptGenParams) {

      RooPullVar pull(name,title,*par,*genParOrig) ;
      _fitParData->addColumn(pull,kFALSE) ;

    } else {
      // If not use fixed generator value
      genParOrig = static_cast<RooAbsReal*>(_genInitParams->find(par->GetName()));

      if (!genParOrig) {
        std::size_t index = it - _fitParams->begin();
        genParOrig = index < _genInitParams->size() ?
            static_cast<RooAbsReal*>((*_genInitParams)[index]) :
            nullptr;

        if (genParOrig) {
          coutW(Generation) << "The fit parameter '" << par->GetName() << "' is not in the model that was used to generate toy data. "
              "The parameter '" << genParOrig->GetName() << "'=" << genParOrig->getVal() << " was found at the same position in the generator model."
              " It will be used to compute pulls."
              "\nIf this is not desired, the parameters of the generator model need to be renamed or reordered." << std::endl;
        }
      }

      if (genParOrig) {
        std::unique_ptr<RooAbsReal> genPar(static_cast<RooAbsReal*>(genParOrig->Clone("truth")));
        RooPullVar pull(name,title,*par,*genPar);

        _fitParData->addColumn(pull,kFALSE) ;
      } else {
        coutE(Generation) << "Cannot generate pull distribution for the fit parameter '" << par->GetName() << "'."
            "\nNo similar parameter was found in the set of parameters that were used to generate toy data." << std::endl;
      }
    }
  }
}




////////////////////////////////////////////////////////////////////////////////
/// Return a RooDataSet containing the post-fit parameters of each toy cycle.
/// This dataset also contains any additional output that was generated
/// by study modules that were added to this RooMCStudy.
/// By default, the two following variables are added (apart from fit parameters):
/// - NLL: The value of the negative log-likelihood for each run.
/// - ngen: Number of events generated for each run.
const RooDataSet& RooMCStudy::fitParDataSet()
{
  if (_canAddFitResults) {
    calcPulls() ;  
    _canAddFitResults = kFALSE ; 
  }
  
  return *_fitParData ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return an argset with the fit parameters for the given sample number
///
/// NB: The fit parameters are only stored for successfull fits,
///     thus the maximum sampleNum can be less that the number
///     of generated samples and if so, the indeces will
///     be out of synch with genData() and fitResult()

const RooArgSet* RooMCStudy::fitParams(Int_t sampleNum) const 
{
  // Check if sampleNum is in range
  if (sampleNum<0 || sampleNum>=_fitParData->numEntries()) {
    oocoutE(_fitModel,InputArguments) << "RooMCStudy::fitParams: ERROR, invalid sample number: " << sampleNum << endl ;    
    return 0 ;
  }
  
  return _fitParData->get(sampleNum) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the RooFitResult of the fit with the given run number.
///
/// \note Fit results are not saved by default. This requires passing `FitOptions(Save(), ...)`
/// to the constructor.
const RooFitResult* RooMCStudy::fitResult(Int_t sampleNum) const
{
  // Check if sampleNum is in range
  if (sampleNum<0 || sampleNum>=_fitResList.GetSize()) {
    oocoutE(_fitModel,InputArguments) << "RooMCStudy::fitResult: ERROR, invalid sample number: " << sampleNum << endl ;    
    return 0 ;
  }
  
  // Retrieve fit result object
  const RooFitResult* fr = (RooFitResult*) _fitResList.At(sampleNum) ;
  if (fr) {
    return fr ;
  } else {
    oocoutE(_fitModel,InputArguments) << "RooMCStudy::fitResult: ERROR, no fit result saved for sample " 
			  << sampleNum << ", did you use the 'r; fit option?" << endl ;
  }
  return 0 ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the given generated dataset. This method will only return datasets
/// if during the run cycle it was indicated that generator data should be saved.

RooAbsData* RooMCStudy::genData(Int_t sampleNum) const 
{
  // Check that generated data was saved
  if (_genDataList.GetSize()==0) {
    oocoutE(_fitModel,InputArguments) << "RooMCStudy::genData() ERROR, generated data was not saved" << endl ;
    return 0 ;
  }
  
  // Check if sampleNum is in range
  if (sampleNum<0 || sampleNum>=_genDataList.GetSize()) {
    oocoutE(_fitModel,InputArguments) << "RooMCStudy::genData() ERROR, invalid sample number: " << sampleNum << endl ;    
    return 0 ;
  }
  
  return  (RooAbsData*) _genDataList.At(sampleNum) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Plot the distribution of fitted values of a parameter. The parameter shown is the one from which the RooPlot
/// was created, e.g.
///
/// RooPlot* frame = param.frame(100,-10,10) ;
/// mcstudy.paramOn(frame,LineStyle(kDashed)) ;
/// 
/// Any named arguments passed to plotParamOn() are forwarded to the underlying plotOn() call

RooPlot* RooMCStudy::plotParamOn(RooPlot* frame, const RooCmdArg& arg1, const RooCmdArg& arg2, const RooCmdArg& arg3, const RooCmdArg& arg4, 
   				 const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8) 
{
  _fitParData->plotOn(frame,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8) ;
  return frame ;
}



////////////////////////////////////////////////////////////////////////////////
/// Plot the distribution of the fitted value of the given parameter on a newly created frame.
///
/// <table>
/// <tr><th> Optional arguments <th>
/// <tr><td> FrameRange(double lo, double hi) <td> Set range of frame to given specification
/// <tr><td> FrameBins(int bins)              <td> Set default number of bins of frame to given number
/// <tr><td> Frame()                       <td> Pass supplied named arguments to RooAbsRealLValue::frame() function. See there
///     for list of allowed arguments
/// </table>
/// If no frame specifications are given, the AutoRange() feature will be used to set the range
/// Any other named argument is passed to the RooAbsData::plotOn() call. See that function for allowed options

RooPlot* RooMCStudy::plotParam(const char* paramName, const RooCmdArg& arg1, const RooCmdArg& arg2, const RooCmdArg& arg3, const RooCmdArg& arg4, 
   			       const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8) 
{

  // Find parameter in fitParDataSet
  RooRealVar* param = static_cast<RooRealVar*>(_fitParData->get()->find(paramName)) ;
  if (!param) {
    oocoutE(_fitModel,InputArguments) << "RooMCStudy::plotParam: ERROR: no parameter defined with name " << paramName << endl ;  
    return 0 ;
  }

  // Forward to implementation below
  return plotParam(*param,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Plot the distribution of the fitted value of the given parameter on a newly created frame.
/// \copydetails RooMCStudy::plotParam(const char* paramName, const RooCmdArg& arg1, const RooCmdArg& arg2, const RooCmdArg& arg3, const RooCmdArg& arg4,
/// const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8)

RooPlot* RooMCStudy::plotParam(const RooRealVar& param, const RooCmdArg& arg1, const RooCmdArg& arg2, const RooCmdArg& arg3, const RooCmdArg& arg4, 
   			       const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8) 
{
  // Stuff all arguments in a list
  RooLinkedList cmdList;
  cmdList.Add(const_cast<RooCmdArg*>(&arg1)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg2)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg3)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg4)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg5)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg6)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg7)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg8)) ;
  
  RooPlot* frame = makeFrameAndPlotCmd(param, cmdList) ;
  if (frame) {
    _fitParData->plotOn(frame, cmdList) ;
  }

  return frame ;
}



////////////////////////////////////////////////////////////////////////////////
/// Plot the distribution of the -log(L) values on a newly created frame.
///
/// <table>
/// <tr><th> Optional arguments <th>
/// <tr><td> FrameRange(double lo, double hi) <td> Set range of frame to given specification
/// <tr><td> FrameBins(int bins)              <td> Set default number of bins of frame to given number
/// <tr><td> Frame()                       <td> Pass supplied named arguments to RooAbsRealLValue::frame() function. See there
///     for list of allowed arguments
/// </table>
///
/// If no frame specifications are given, the AutoRange() feature will be used to set the range.
/// Any other named argument is passed to the RooAbsData::plotOn() call. See that function for allowed options

RooPlot* RooMCStudy::plotNLL(const RooCmdArg& arg1, const RooCmdArg& arg2,
                     const RooCmdArg& arg3, const RooCmdArg& arg4,
                     const RooCmdArg& arg5, const RooCmdArg& arg6,
                     const RooCmdArg& arg7, const RooCmdArg& arg8) 
{
  return plotParam(*_nllVar,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Plot the distribution of the fit errors for the specified parameter on a newly created frame.
///
/// <table>
/// <tr><th> Optional arguments <th>
/// <tr><td> FrameRange(double lo, double hi) <td> Set range of frame to given specification
/// <tr><td> FrameBins(int bins)              <td> Set default number of bins of frame to given number
/// <tr><td> Frame()                       <td> Pass supplied named arguments to RooAbsRealLValue::frame() function. See there
///     for list of allowed arguments
/// </table>
///
/// If no frame specifications are given, the AutoRange() feature will be used to set a default range.
/// Any other named argument is passed to the RooAbsData::plotOn() call. See that function for allowed options.

RooPlot* RooMCStudy::plotError(const RooRealVar& param, const RooCmdArg& arg1, const RooCmdArg& arg2,
                     const RooCmdArg& arg3, const RooCmdArg& arg4,
                     const RooCmdArg& arg5, const RooCmdArg& arg6,
                     const RooCmdArg& arg7, const RooCmdArg& arg8) 
{
  if (_canAddFitResults) {
    calcPulls() ;
    _canAddFitResults=kFALSE ;
  }

  RooErrorVar* evar = param.errorVar() ;
  RooRealVar* evar_rrv = static_cast<RooRealVar*>(evar->createFundamental()) ;
  RooPlot* frame = plotParam(*evar_rrv,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8) ;
  delete evar_rrv ;
  delete evar ;
  return frame ;
}



////////////////////////////////////////////////////////////////////////////////
/// Plot the distribution of pull values for the specified parameter on a newly created frame. If asymmetric
/// errors are calculated in the fit (by MINOS) those will be used in the pull calculation.
///
/// If the parameters of the models for generation and fit differ, simple heuristics are used to find the
/// corresponding parameters:
/// - Parameters have the same name: They will be used to compute pulls.
/// - Parameters have different names: The position of the fit parameter in the set of fit parameters will be
///   computed. The parameter at the same position in the set of generator parameters will be used.
///
/// Further options:
/// <table>
/// <tr><th> Arguments <th> Effect
/// <tr><td> FrameRange(double lo, double hi) <td> Set range of frame to given specification
/// <tr><td> FrameBins(int bins)              <td> Set default number of bins of frame to given number
/// <tr><td> Frame()                       <td> Pass supplied named arguments to RooAbsRealLValue::frame() function. See there
///     for list of allowed arguments
/// <tr><td> FitGauss(Bool_t flag)            <td> Add a gaussian fit to the frame
/// </table>
///
/// If no frame specifications are given, the AutoSymRange() feature will be used to set a default range.
/// Any other named argument is passed to the RooAbsData::plotOn(). See that function for allowed options.

RooPlot* RooMCStudy::plotPull(const RooRealVar& param, const RooCmdArg& arg1, const RooCmdArg& arg2,
                     const RooCmdArg& arg3, const RooCmdArg& arg4,
                     const RooCmdArg& arg5, const RooCmdArg& arg6,
                     const RooCmdArg& arg7, const RooCmdArg& arg8) 
{
  // Stuff all arguments in a list
  RooLinkedList cmdList;
  cmdList.Add(const_cast<RooCmdArg*>(&arg1)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg2)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg3)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg4)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg5)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg6)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg7)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg8)) ;

  TString name(param.GetName()), title(param.GetTitle()) ;
  name.Append("pull") ; title.Append(" Pull") ;
  RooRealVar pvar(name,title,-100,100) ;
  pvar.setBins(100) ;


  RooPlot* frame = makeFrameAndPlotCmd(pvar, cmdList, kTRUE) ;
  if (frame) {

    // Pick up optonal FitGauss command from list
    RooCmdConfig pc(Form("RooMCStudy::plotPull(%s)",_genModel->GetName())) ;
    pc.defineInt("fitGauss","FitGauss",0,0) ;
    pc.allowUndefined() ;
    pc.process(cmdList) ;
    Bool_t fitGauss=pc.getInt("fitGauss") ;

    // Pass stripped command list to plotOn()
    pc.stripCmdList(cmdList,"FitGauss") ;
    const bool success = _fitParData->plotOn(frame,cmdList) ;

    if (!success) {
      coutF(Plotting) << "No pull distribution for the parameter '" << param.GetName() << "'. Check logs for errors." << std::endl;
      return frame;
    }

    // Add Gaussian fit if requested
    if (fitGauss) {
      RooRealVar pullMean("pullMean","Mean of pull",0,-10,10) ;
      RooRealVar pullSigma("pullSigma","Width of pull",1,0.1,5) ;
      RooGenericPdf pullGauss("pullGauss","Gaussian of pull",
			      "exp(-0.5*(@0-@1)*(@0-@1)/(@2*@2))",
			      RooArgSet(pvar,pullMean,pullSigma)) ;
      pullGauss.fitTo(*_fitParData,RooFit::Minos(0),RooFit::PrintLevel(-1)) ;
      pullGauss.plotOn(frame) ;
      pullGauss.paramOn(frame,_fitParData) ;
    }
  }
  return frame;
}



////////////////////////////////////////////////////////////////////////////////
/// Internal function. Construct RooPlot from given parameter and modify the list of named
/// arguments 'cmdList' to only contain the plot arguments that should be forwarded to 
/// RooAbsData::plotOn()

RooPlot* RooMCStudy::makeFrameAndPlotCmd(const RooRealVar& param, RooLinkedList& cmdList, Bool_t symRange) const 
{
  // Select the frame-specific commands 
  RooCmdConfig pc(Form("RooMCStudy::plotParam(%s)",_genModel->GetName())) ;
  pc.defineInt("nbins","Bins",0,0) ;
  pc.defineDouble("xlo","Range",0,0) ;
  pc.defineDouble("xhi","Range",1,0) ;
  pc.defineInt("dummy","FrameArgs",0,0) ;
  pc.defineMutex("Bins","FrameArgs") ;
  pc.defineMutex("Range","FrameArgs") ;

  // Process and check varargs 
  pc.allowUndefined() ;
  pc.process(cmdList) ;
  if (!pc.ok(kTRUE)) {
    return 0 ;
  }
  
  // Make frame according to specs
  Int_t nbins = pc.getInt("nbins") ;
  Double_t xlo = pc.getDouble("xlo") ;
  Double_t xhi = pc.getDouble("xhi") ;
  RooPlot* frame ; 

  if (pc.hasProcessed("FrameArgs")) {
    // Explicit frame arguments are given, pass them on
    RooCmdArg* frameArg = static_cast<RooCmdArg*>(cmdList.FindObject("FrameArgs")) ;
    frame = param.frame(frameArg->subArgs()) ;
  } else {
    // FrameBins, FrameRange or none are given, build custom frame command list
    RooCmdArg bins = RooFit::Bins(nbins) ;
    RooCmdArg range = RooFit::Range(xlo,xhi) ;
    RooCmdArg autor = symRange ? RooFit::AutoSymRange(*_fitParData,0.2) : RooFit::AutoRange(*_fitParData,0.2) ;
    RooLinkedList frameCmdList ;

    if (pc.hasProcessed("Bins")) frameCmdList.Add(&bins) ;
    if (pc.hasProcessed("Range")) {
      frameCmdList.Add(&range) ;
    } else {
      frameCmdList.Add(&autor) ;
    }
    frame = param.frame(frameCmdList) ;
  }
  
  // Filter frame command from list and pass on to plotOn() 
  pc.stripCmdList(cmdList,"FrameArgs,Bins,Range") ;

  return frame ;
}



////////////////////////////////////////////////////////////////////////////////
/// Create a RooPlot of the -log(L) distribution in the range lo-hi
/// with 'nBins' bins

RooPlot* RooMCStudy::plotNLL(Double_t lo, Double_t hi, Int_t nBins) 
{
  RooPlot* frame = _nllVar->frame(lo,hi,nBins) ;
  
  _fitParData->plotOn(frame) ;
  return frame ;
}



////////////////////////////////////////////////////////////////////////////////
/// Create a RooPlot of the distribution of the fitted errors of the given parameter. 
/// The frame is created with a range [lo,hi] and plotted data will be binned in 'nbins' bins

RooPlot* RooMCStudy::plotError(const RooRealVar& param, Double_t lo, Double_t hi, Int_t nbins) 
{
  if (_canAddFitResults) {
    calcPulls() ;
    _canAddFitResults=kFALSE ;
  }

  RooErrorVar* evar = param.errorVar() ;
  RooPlot* frame = evar->frame(lo,hi,nbins) ;
  _fitParData->plotOn(frame) ;

  delete evar ;
  return frame ;
}



////////////////////////////////////////////////////////////////////////////////
/// Create a RooPlot of the pull distribution for the given
/// parameter.  The range lo-hi is plotted in nbins.  If fitGauss is
/// set, an unbinned ML fit of the distribution to a Gaussian p.d.f
/// is performed. The fit result is overlaid on the returned RooPlot
/// and a box with the fitted mean and sigma is added.
///
/// If the parameters of the models for generation and fit differ, simple heuristics are used to find the
/// corresponding parameters:
/// - Parameters have the same name: They will be used to compute pulls.
/// - Parameters have different names: The position of the fit parameter in the set of fit parameters will be
///   computed. The parameter at the same position in the set of generator parameters will be used.

RooPlot* RooMCStudy::plotPull(const RooRealVar& param, Double_t lo, Double_t hi, Int_t nbins, Bool_t fitGauss) 
{
  if (_canAddFitResults) {
    calcPulls() ;
    _canAddFitResults=kFALSE ;
  }


  TString name(param.GetName()), title(param.GetTitle()) ;
  name.Append("pull") ; title.Append(" Pull") ;
  RooRealVar pvar(name,title,lo,hi) ;
  pvar.setBins(nbins) ;

  RooPlot* frame = pvar.frame() ;
  const bool success = _fitParData->plotOn(frame);

  if (!success) {
    coutF(Plotting) << "No pull distribution for the parameter '" << param.GetName() << "'. Check logs for errors." << std::endl;
    return frame;
  }

  if (fitGauss) {
    RooRealVar pullMean("pullMean","Mean of pull",0,lo,hi) ;
    RooRealVar pullSigma("pullSigma","Width of pull",1,0,5) ;
    RooGenericPdf pullGauss("pullGauss","Gaussian of pull",
			    "exp(-0.5*(@0-@1)*(@0-@1)/(@2*@2))",
			    RooArgSet(pvar,pullMean,pullSigma)) ;
    pullGauss.fitTo(*_fitParData,RooFit::Minos(0),RooFit::PrintLevel(-1)) ;
    pullGauss.plotOn(frame) ;
    pullGauss.paramOn(frame,_fitParData) ;
  }

  return frame ;
}


////////////////////////////////////////////////////////////////////////////////
/// If one of the TObject we have a referenced to is deleted, remove the
/// reference.

void RooMCStudy::RecursiveRemove(TObject *obj)
{
   _fitResList.RecursiveRemove(obj);
   _genDataList.RecursiveRemove(obj);
   _fitOptList.RecursiveRemove(obj);
   if (_ngenVar == obj) _ngenVar = nullptr;

   if (_fitParData) _fitParData->RecursiveRemove(obj);
   if (_fitParData == obj) _fitParData = nullptr;

   if (_genParData) _genParData->RecursiveRemove(obj);
   if (_genParData == obj) _genParData = nullptr;
}

