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
\file RooRealIntegral.cxx
\class RooRealIntegral
\ingroup Roofitcore

RooRealIntegral performs hybrid numerical/analytical integrals of RooAbsReal objects.
The class performs none of the actual integration, but only manages the logic
of what variables can be integrated analytically, accounts for eventual jacobian
terms and defines what numerical integrations needs to be done to complement the
analytical integral.
The actual analytical integrations (if any) are done in the PDF themselves, the numerical
integration is performed in the various implementations of the RooAbsIntegrator base class.
**/

#include "RooRealIntegral.h"

#include "RooFit.h"

#include "RooMsgService.h"
#include "RooArgSet.h"
#include "RooAbsRealLValue.h"
#include "RooAbsCategoryLValue.h"
#include "RooRealBinding.h"
#include "RooRealAnalytic.h"
#include "RooInvTransform.h"
#include "RooSuperCategory.h"
#include "RooNumIntFactory.h"
#include "RooNumIntConfig.h"
#include "RooNameReg.h"
#include "RooExpensiveObjectCache.h"
#include "RooConstVar.h"
#include "RooDouble.h"
#include "RooTrace.h"
#include "RooHelpers.h"

#include "ROOT/RMakeUnique.hxx"
#include "RooTimer.h"
// getpid and getppid:
#include "unistd.h"

#include "TClass.h"

#include <iostream>

using namespace std;

ClassImp(RooRealIntegral);


Int_t RooRealIntegral::_cacheAllNDim(2) ;


////////////////////////////////////////////////////////////////////////////////

RooRealIntegral::RooRealIntegral() : 
  _valid(kFALSE),
  _respectCompSelect(true),
  _funcNormSet(0),
  _iconfig(0),
  _mode(0),
  _intOperMode(Hybrid),
  _restartNumIntEngine(kFALSE),
  _numIntEngine(0),
  _numIntegrand(0),
  _rangeName(0),
  _params(0),
  _cacheNum(kFALSE),
  _timeNumInt(kFALSE)
{
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Construct integral of 'function' over observables in 'depList'
/// in range 'rangeName'  with normalization observables 'funcNormSet' 
/// (for p.d.f.s). In the integral is performed to the maximum extent
/// possible the internal (analytical) integrals advertised by function.
/// The other integrations are performed numerically. The optional
/// config object prescribes how these numeric integrations are configured.
///

RooRealIntegral::RooRealIntegral(const char *name, const char *title, 
				 const RooAbsReal& function, const RooArgSet& depList,
				 const RooArgSet* funcNormSet, const RooNumIntConfig* config,
				 const char* rangeName) :
  RooAbsReal(name,title), 
  _valid(kTRUE), 
  _respectCompSelect(true),
  _sumList("!sumList","Categories to be summed numerically",this,kFALSE,kFALSE), 
  _intList("!intList","Variables to be integrated numerically",this,kFALSE,kFALSE), 
  _anaList("!anaList","Variables to be integrated analytically",this,kFALSE,kFALSE), 
  _jacList("!jacList","Jacobian product term",this,kFALSE,kFALSE), 
  _facList("!facList","Variables independent of function",this,kFALSE,kTRUE),
  _function("!func","Function to be integrated",this,
	    const_cast<RooAbsReal&>(function),kFALSE,kFALSE), 
  _iconfig((RooNumIntConfig*)config),
  _sumCat("!sumCat","SuperCategory for summation",this,kFALSE,kFALSE),
  _mode(0),
  _intOperMode(Hybrid), 
  _restartNumIntEngine(kFALSE),
  _numIntEngine(0), 
  _numIntegrand(0),
  _rangeName((TNamed*)RooNameReg::ptr(rangeName)),
  _params(0),
  _cacheNum(kFALSE)
{
  //   A) Check that all dependents are lvalues 
  //
  //   B) Check if list of dependents can be re-expressed in        
  //      lvalues that are higher in the expression tree            
  //
  //   C) Check for dependents that the PDF insists on integrating  
  //      analytically iself                                        
  //
  //   D) Make list of servers that can be integrated analytically  
  //      Add all parameters/dependents as value/shape servers      
  //
  //   E) Interact with function to make list of objects actually integrated analytically   
  //
  //   F) Make list of numerical integration variables consisting of:               
  //     - Category dependents of RealLValues in analytical integration             
  //     - Leaf nodes server lists of function server that are not analytically integrated   
  //     - Make Jacobian list for analytically integrated RealLValues            
  //
  //   G) Split numeric list in integration list and summation list   
  //

  oocxcoutI(&function,Integration) << "RooRealIntegral::ctor(" << GetName() << ") Constructing integral of function " 
				     << function.GetName() << " over observables" << depList << " with normalization " 
				     << (funcNormSet?*funcNormSet:RooArgSet()) << " with range identifier " 
				     << (rangeName?rangeName:"<none>") << endl ;

  
  // Choose same expensive object cache as integrand
  setExpensiveObjectCache(function.expensiveObjectCache()) ;
//   cout << "RRI::ctor(" << GetName() << ") setting expensive object cache to " << &expensiveObjectCache() << " as taken from " << function.GetName() << endl ;

  // Use objects integrator configuration if none is specified
  if (!_iconfig) _iconfig = (RooNumIntConfig*) function.getIntegratorConfig() ;

  // Save private copy of funcNormSet, if supplied, excluding factorizing terms
  if (funcNormSet) {
    _funcNormSet = new RooArgSet ;
    for (const auto nArg : *funcNormSet) {
      if (function.dependsOn(*nArg)) {
        _funcNormSet->addClone(*nArg) ;
      }
    }
  } else {
    _funcNormSet = 0 ;
  }

  //_funcNormSet = funcNormSet ? (RooArgSet*)funcNormSet->snapshot(kFALSE) : 0 ;
  
  // Make internal copy of dependent list
  RooArgSet intDepList(depList) ;

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  // * A) Check that all dependents are lvalues and filter out any
  //      dependents that the PDF doesn't explicitly depend on
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  
  for (auto arg : intDepList) {
    if(!arg->isLValue()) {
      coutE(InputArguments) << ClassName() << "::" << GetName() << ": cannot integrate non-lvalue ";
      arg->Print("1");
      _valid= kFALSE;
    }
    if (!function.dependsOn(*arg)) {
      RooAbsArg* argClone = (RooAbsArg*) arg->Clone() ;
      _facListOwned.addOwned(*argClone) ;
      _facList.add(*argClone) ;
      addServer(*argClone,kFALSE,kTRUE) ;
    }
  }

  if (_facList.getSize()>0) {
    oocxcoutI(&function,Integration) << function.GetName() << ": Factorizing obserables are " << _facList << endl ;
  }


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  // * B) Check if list of dependents can be re-expressed in       *
  // *    lvalues that are higher in the expression tree           *
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


  // Initial fill of list of LValue branches
  RooArgSet exclLVBranches("exclLVBranches") ;
  RooArgSet branchList,branchListVD ;
  function.branchNodeServerList(&branchList) ;

  for (auto branch: branchList) {
    RooAbsRealLValue    *realArgLV = dynamic_cast<RooAbsRealLValue*>(branch) ;
    RooAbsCategoryLValue *catArgLV = dynamic_cast<RooAbsCategoryLValue*>(branch) ;
    if ((realArgLV && (realArgLV->isJacobianOK(intDepList)!=0)) || catArgLV) {
      exclLVBranches.add(*branch) ;
//       cout << "exclv branch = " << endl ;
//       branch->printCompactTree() ;
    }
    if (dependsOnValue(*branch)) {
      branchListVD.add(*branch) ;
    } else {
//       cout << "value of self does not depend on branch " << branch->GetName() << endl ;
    }
  }
  exclLVBranches.remove(depList,kTRUE,kTRUE) ;

  // Initial fill of list of LValue leaf servers (put in intDepList)
  RooArgSet exclLVServers("exclLVServers") ;
  exclLVServers.add(intDepList) ;
  
  // Obtain mutual exclusive dependence by iterative reduction
  Bool_t converged(kFALSE) ;
  while(!converged) {
    converged=kTRUE ;

    // Reduce exclLVServers to only those serving exclusively exclLVBranches
    std::vector<RooAbsArg*> toBeRemoved;
    for (auto server : exclLVServers) {
      if (!servesExclusively(server,exclLVBranches,branchListVD)) {
        toBeRemoved.push_back(server);
        converged=kFALSE ;
      }
    }
    exclLVServers.remove(toBeRemoved.begin(), toBeRemoved.end());
    
    // Reduce exclLVBranches to only those depending exclusively on exclLVservers
    // Attention: counting loop, since erasing from container
    for (std::size_t i=0; i < exclLVBranches.size(); ++i) {
      const RooAbsArg* branch = exclLVBranches[i];
      RooArgSet* brDepList = branch->getObservables(&intDepList) ;
      RooArgSet bsList(*brDepList,"bsList") ;
      delete brDepList ;
      bsList.remove(exclLVServers,kTRUE,kTRUE) ;
      if (bsList.getSize()>0) {
        exclLVBranches.remove(*branch,kTRUE,kTRUE) ;
        --i;
        converged=kFALSE ;
      }
    }
  }

  // Eliminate exclLVBranches that do not depend on any LVServer
  // Attention: Counting loop, since modifying container
  for (std::size_t i=0; i < exclLVBranches.size(); ++i) {
    const RooAbsArg* branch = exclLVBranches[i];
    if (!branch->dependsOnValue(exclLVServers)) {
      exclLVBranches.remove(*branch,kTRUE,kTRUE) ;
      --i;
    }
  }
     
  // Replace exclusive lvalue branch servers with lvalue branches
  // WVE Don't do this for binned distributions - deal with this using numeric integration with transformed bin boundaroes
  if (exclLVServers.getSize()>0 && !function.isBinnedDistribution(exclLVBranches)) {
//     cout << "activating LVservers " << exclLVServers << " for use in integration " << endl ;
    intDepList.remove(exclLVServers) ;
    intDepList.add(exclLVBranches) ;

    //cout << "intDepList removing exclLVServers " << exclLVServers << endl ;
    //cout << "intDepList adding exclLVBranches " << exclLVBranches << endl ;

  }

     
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  // * C) Check for dependents that the PDF insists on integrating *
  //      analytically iself                                       *
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  RooArgSet anIntOKDepList ;
  for (auto arg : intDepList) {
    if (function.forceAnalyticalInt(*arg)) {
      anIntOKDepList.add(*arg) ;
    }
  }
  
  if (anIntOKDepList.getSize()>0) {
    oocxcoutI(&function,Integration) << function.GetName() << ": Observables that function forcibly requires to be integrated internally " << anIntOKDepList << endl ;
  }


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  // * D) Make list of servers that can be integrated analytically *
  //      Add all parameters/dependents as value/shape servers     *
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  for (const auto arg : function.servers()) {

    //cout << "considering server" << arg->GetName() << endl ;

    // Dependent or parameter?
    if (!arg->dependsOnValue(intDepList)) {

      //cout << " server does not depend on observables, adding server as value server to integral" << endl ;

      if (function.dependsOnValue(*arg)) {
        addServer(*arg,kTRUE,kFALSE) ;
      }

      continue ;

    } else {

      // Add final dependents of arg as shape servers
      RooArgSet argLeafServers ;
      arg->leafNodeServerList(&argLeafServers,0,kFALSE) ;

      //arg->printCompactTree() ;
      //cout << "leaf nodes of server are " << argLeafServers << " depList = " << depList << endl ;

      // Skip arg if it is neither value or shape server
      if (!arg->isValueServer(function) && !arg->isShapeServer(function)) {
        //cout << " server is neither value not shape server of function, ignoring" << endl ;
        continue ;
      }

      for (const auto leaf : argLeafServers) {

        //cout << " considering leafnode " << leaf->GetName() << " of server " << arg->GetName() << endl ;

        if (depList.find(leaf->GetName()) && function.dependsOnValue(*leaf)) {

          RooAbsRealLValue* leaflv = dynamic_cast<RooAbsRealLValue*>(leaf) ;
          if (leaflv && leaflv->getBinning(rangeName).isParameterized()) {
            oocxcoutD(&function,Integration) << function.GetName() << " : Observable " << leaf->GetName() << " has parameterized binning, add value dependence of boundary objects rather than shape of leaf" << endl ;
            if (leaflv->getBinning(rangeName).lowBoundFunc()) {
              addServer(*leaflv->getBinning(rangeName).lowBoundFunc(),kTRUE,kFALSE) ;
            }
            if(leaflv->getBinning(rangeName).highBoundFunc()) {
              addServer(*leaflv->getBinning(rangeName).highBoundFunc(),kTRUE,kFALSE) ;
            }
          } else {
            oocxcoutD(&function,Integration) << function.GetName() << ": Adding observable " << leaf->GetName() << " of server "
                << arg->GetName() << " as shape dependent" << endl ;
            addServer(*leaf,kFALSE,kTRUE) ;
          }
        } else if (!depList.find(leaf->GetName())) {

          if (function.dependsOnValue(*leaf)) {
            oocxcoutD(&function,Integration) << function.GetName() << ": Adding parameter " << leaf->GetName() << " of server " << arg->GetName() << " as value dependent" << endl ;
            addServer(*leaf,kTRUE,kFALSE) ;
          } else {
            oocxcoutD(&function,Integration) << function.GetName() << ": Adding parameter " << leaf->GetName() << " of server " << arg->GetName() << " as shape dependent" << endl ;
            addServer(*leaf,kFALSE,kTRUE) ;
          }
        }
      }
    }

    // If this dependent arg is self-normalized, stop here
    //if (function.selfNormalized()) continue ;

    Bool_t depOK(kFALSE) ;
    // Check for integratable AbsRealLValue

    //cout << "checking server " << arg->IsA()->GetName() << "::" << arg->GetName() << endl ;

    if (arg->isDerived()) {
      RooAbsRealLValue    *realArgLV = dynamic_cast<RooAbsRealLValue*>(arg) ;
      RooAbsCategoryLValue *catArgLV = dynamic_cast<RooAbsCategoryLValue*>(arg) ;
      //cout << "realArgLV = " << realArgLV << " intDepList = " << intDepList << endl ;
      if ((realArgLV && intDepList.find(realArgLV->GetName()) && (realArgLV->isJacobianOK(intDepList)!=0)) || catArgLV) {	

 	//cout  << " arg " << arg->GetName() << " is derived LValue with valid jacobian" << endl ;

	// Derived LValue with valid jacobian
	depOK = kTRUE ;
	
	// Now, check for overlaps
	Bool_t overlapOK = kTRUE ;
	for (const auto otherArg : function.servers()) {
	  // skip comparison with self
	  if (arg==otherArg) continue ;
	  if (otherArg->IsA()==RooConstVar::Class()) continue ;
	  if (arg->overlaps(*otherArg,kTRUE)) {
	    //cout << "arg " << arg->GetName() << " overlaps with " << otherArg->GetName() << endl ;
	    //overlapOK=kFALSE ;
	  }
	}      	
	// coverity[DEADCODE]
	if (!overlapOK) depOK=kFALSE ;      

	//cout << "overlap check returns OK=" << (depOK?"T":"F") << endl ;

      }
    } else {
      // Fundamental types are always OK
      depOK = kTRUE ;
    }
    
    // Add server to list of dependents that are OK for analytical integration
    if (depOK) {
      anIntOKDepList.add(*arg,kTRUE) ;      
      oocxcoutI(&function,Integration) << function.GetName() << ": Observable " << arg->GetName() << " is suitable for analytical integration (if supported by p.d.f)" << endl ;
    }
  }
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  // * E) interact with function to make list of objects actually integrated analytically  *
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  RooArgSet anIntDepList ;

  RooArgSet *anaSet = new RooArgSet( _anaList, Form("UniqueCloneOf_%s",_anaList.GetName()));
  _mode = ((RooAbsReal&)_function.arg()).getAnalyticalIntegralWN(anIntOKDepList,*anaSet,_funcNormSet,RooNameReg::str(_rangeName)) ;    
  _anaList.removeAll() ;
  _anaList.add(*anaSet);    
  delete anaSet;

  // Avoid confusion -- if mode is zero no analytical integral is defined regardless of contents of _anaListx
  if (_mode==0) {
    _anaList.removeAll() ;
  }

  if (_mode!=0) {
    oocxcoutI(&function,Integration) << function.GetName() << ": Function integrated observables " << _anaList << " internally with code " << _mode << endl ;
  }

  // when _funcNormSet is a nullptr a warning message appears for RooAddPdf functions
  // This is not a problem since we do noty use the returned value from getVal()
  // we then disable the produced warning message in the RooFit::Eval topic
  std::unique_ptr<RooHelpers::LocalChangeMsgLevel> msgChanger;
  if (_funcNormSet == nullptr) {
     // remove only the RooFit::Eval message topic from current active streams
     // passed level can be whatever if we provide a false as last argument   
     msgChanger = std::make_unique<RooHelpers::LocalChangeMsgLevel>(RooFit::WARNING, 0u, RooFit::Eval, false);
  }

  // WVE kludge: synchronize dset for use in analyticalIntegral
  // LM : is this really needed ??
  function.getVal(_funcNormSet) ;
  // delete LocalChangeMsgLevel which will restore previous message level
  msgChanger.reset(nullptr); 

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  // * F) Make list of numerical integration variables consisting of:            *  
  // *   - Category dependents of RealLValues in analytical integration          *  
  // *   - Expanded server lists of server that are not analytically integrated  *
  // *    Make Jacobian list with analytically integrated RealLValues            *
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  RooArgSet numIntDepList ;


  // Loop over actually analytically integrated dependents
  for (const auto arg : _anaList) {

    // Process only derived RealLValues
    if (arg->IsA()->InheritsFrom(RooAbsRealLValue::Class()) && arg->isDerived() && !arg->isFundamental()) {

      // Add to list of Jacobians to calculate
      _jacList.add(*arg) ;

      // Add category dependent of LValueReal used in integration
      auto argDepList = std::unique_ptr<RooArgSet>(arg->getObservables(&intDepList));
      for (const auto argDep : *argDepList) {
        if (argDep->IsA()->InheritsFrom(RooAbsCategoryLValue::Class()) && intDepList.contains(*argDep)) {
          numIntDepList.add(*argDep,kTRUE) ;
        }
      }
    }
  }

  
  // If nothing was integrated analytically, swap back LVbranches for LVservers for subsequent numeric integration
  if (_anaList.getSize()==0) {    
    if (exclLVServers.getSize()>0) {
      //cout << "NUMINT phase analList is empty. exclLVServers = " << exclLVServers << endl ;
      intDepList.remove(exclLVBranches) ;      
      intDepList.add(exclLVServers) ;
     }            
  }
  //cout << "NUMINT intDepList = " << intDepList << endl ;

  // Loop again over function servers to add remaining numeric integrations
  for (const auto arg : function.servers()) {

    //cout << "processing server for numeric integration " << arg->IsA()->GetName() << "::" << arg->GetName() << endl ;

    // Process only servers that are not treated analytically
    if (!_anaList.find(arg->GetName()) && arg->dependsOn(intDepList)) {

      // Process only derived RealLValues
      if (dynamic_cast<RooAbsLValue*>(arg) && arg->isDerived() && intDepList.contains(*arg)) {
        numIntDepList.add(*arg,kTRUE) ;
      } else {

        // WVE this will only get the observables, but not l-value transformations
        // Expand server in final dependents
        auto argDeps = std::unique_ptr<RooArgSet>(arg->getObservables(&intDepList));

        if (argDeps->getSize()>0) {

          // Add final dependents, that are not forcibly integrated analytically,
          // to numerical integration list
          for (const auto dep : *argDeps) {
            if (!_anaList.find(dep->GetName())) {
              numIntDepList.add(*dep,kTRUE) ;
            }
          }
        }
      }
    }
  }

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  // * G) Split numeric list in integration list and summation list  *
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  // Split numeric integration list in summation and integration lists
  for (const auto arg : numIntDepList) {
    if (arg->IsA()->InheritsFrom(RooAbsRealLValue::Class())) {
      _intList.add(*arg) ;
    } else if (arg->IsA()->InheritsFrom(RooAbsCategoryLValue::Class())) {
      _sumList.add(*arg) ;
    }
  }

  if (_anaList.getSize()>0) {
    oocxcoutI(&function,Integration) << function.GetName() << ": Observables " << _anaList << " are analytically integrated with code " << _mode << endl ;
  }
  if (_intList.getSize()>0) {
    oocxcoutI(&function,Integration) << function.GetName() << ": Observables " << _intList << " are numerically integrated" << endl ;
  }
  if (_sumList.getSize()>0) {
    oocxcoutI(&function,Integration) << function.GetName() << ": Observables " << _sumList << " are numerically summed" << endl ;
  }
  

  // Determine operating mode
  if (numIntDepList.getSize()>0) {
    // Numerical and optional Analytical integration
    _intOperMode = Hybrid ;
  } else if (_anaList.getSize()>0) {
    // Purely analytical integration
    _intOperMode = Analytic ;    
  } else {
    // No integration performed
    _intOperMode = PassThrough ;
  }

  // Determine auto-dirty status  
  autoSelectDirtyMode() ;

  // Create value caches for _intList and _sumList
  _intList.snapshot(_saveInt) ;
  _sumList.snapshot(_saveSum) ;

  
  if (_sumList.getSize()>0) {
    RooSuperCategory *sumCat = new RooSuperCategory(Form("%s_sumCat",GetName()),"sumCat",_sumList) ;
    _sumCat.addOwned(*sumCat) ;
  }

  // activate timing on numerical integrals
  if (RooTrace::time_numInts()) {
    activateTimingNumInts();
  }

  TRACE_CREATE
}


void RooRealIntegral::activateTimingNumInts() {
  // activate timing on numerical integrals
//  const RooAbsArg& pdfNode = _function.arg();
  const RooAbsPdf& pdfNode = dynamic_cast<const RooAbsPdf&>(_function.arg());
  // TODO: check whether indeed only RooAbsPdf _functions have integral components
  Bool_t timing_flag = pdfNode.num_int_timing_flag();

  RooFIter ni_iter = _intList.fwdIterator();
  while (RooAbsArg *normint = ni_iter.next()) {
    // Integral expressions can be composite objects (in case of disjoint normalization ranges)
    // Therefore: retrieve list of branch nodes of integral expression
    RooArgList bi;
    normint->branchNodeServerList(&bi);
    RooFIter ib_iter = bi.fwdIterator();
    RooAbsArg *inode;
    while ((inode = ib_iter.next())) {
      // If a RooRealIntegal component is found...
      if (inode->IsA() == RooRealIntegral::Class()) {
        // Retrieve the number of real dimensions that is integrated numerically,
        RooRealIntegral *rri = (RooRealIntegral *) inode;
        rri->setNumIntTiming(timing_flag);
      }
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Set appropriate cache operation mode for integral depending on cache operation
/// mode of server objects

void RooRealIntegral::autoSelectDirtyMode() 
{
  // If any of our servers are is forcedDirty or a projectedDependent, then we need to be ADirty
  for (const auto server : _serverList) {
    if (server->isValueServer(*this)) {
      RooArgSet leafSet ;
      server->leafNodeServerList(&leafSet) ;
      for (const auto leaf : leafSet) {
        if (leaf->operMode()==ADirty && leaf->isValueServer(*this)) {
          setOperMode(ADirty) ;
          break ;
        }
        if (leaf->getAttribute("projectedDependent")) {
          setOperMode(ADirty) ;
          break ;
        }
      }
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Utility function that returns true if 'object server' is a server
/// to exactly one of the RooAbsArgs in 'exclLVBranches'

Bool_t RooRealIntegral::servesExclusively(const RooAbsArg* server,const RooArgSet& exclLVBranches, const RooArgSet& allBranches) const
{  
  // Determine if given server serves exclusively exactly one of the given nodes in exclLVBranches

  // Special case, no LV servers available
  if (exclLVBranches.getSize()==0) return kFALSE ;

  // If server has no clients and is not an LValue itself, return false
   if (server->_clientList.empty() && exclLVBranches.find(server->GetName())) {
     return kFALSE ;
   }

   // WVE must check for value relations only here!!!!

   // Loop over all clients
   Int_t numLVServ(0) ;
   for (const auto client : server->valueClients()) {
     // If client is not an LValue, recurse
     if (!(exclLVBranches.find(client->GetName())==client)) {
       if (allBranches.find(client->GetName())==client) {
         if (!servesExclusively(client,exclLVBranches,allBranches)) {
           // Client is a non-LValue that doesn't have an exclusive LValue server
           return kFALSE ;
         }
       }
     } else {
       // Client is an LValue       
       numLVServ++ ;
     }
   }

   return (numLVServ==1) ;
}




////////////////////////////////////////////////////////////////////////////////
/// (Re)Initialize numerical integration engine if necessary. Return kTRUE if
/// successful, or otherwise kFALSE.

Bool_t RooRealIntegral::initNumIntegrator() const
{
  // if we already have an engine, check if it still works for the present limits.
  if(0 != _numIntEngine) {
    if(_numIntEngine->isValid() && _numIntEngine->checkLimits() && !_restartNumIntEngine ) return kTRUE;
    // otherwise, cleanup the old engine
    delete _numIntEngine ;
    _numIntEngine= 0;
    if(0 != _numIntegrand) {
      delete _numIntegrand;
      _numIntegrand= 0;
    }
  }

  // All done if there are no arguments to integrate numerically
  if(0 == _intList.getSize()) return kTRUE;
  
  // Bind the appropriate analytic integral (specified by _mode) of our RooRealVar object to
  // those of its arguments that will be integrated out numerically.
  if(_mode != 0) {
    _numIntegrand= new RooRealAnalytic(_function.arg(),_intList,_mode,_funcNormSet,_rangeName);
  }
  else {
    _numIntegrand= new RooRealBinding(_function.arg(),_intList,_funcNormSet,kFALSE,_rangeName);
  }
  if(0 == _numIntegrand || !_numIntegrand->isValid()) {
    coutE(Integration) << ClassName() << "::" << GetName() << ": failed to create valid integrand." << endl;
    return kFALSE;
  }

  // Create appropriate numeric integrator using factory
  Bool_t isBinned = _function.arg().isBinnedDistribution(_intList) ;
  _numIntEngine = RooNumIntFactory::instance().createIntegrator(*_numIntegrand,*_iconfig,0,isBinned) ;

  if(0 == _numIntEngine || !_numIntEngine->isValid()) {
    coutE(Integration) << ClassName() << "::" << GetName() << ": failed to create valid integrator." << endl;
    return kFALSE;
  }

  cxcoutI(NumIntegration) << "RooRealIntegral::init(" << GetName() << ") using numeric integrator " 
			  << _numIntEngine->IsA()->GetName() << " to calculate Int" << _intList << endl ;

  if (_intList.getSize()>3) {
    cxcoutI(NumIntegration) << "RooRealIntegral::init(" << GetName() << ") evaluation requires " << _intList.getSize() << "-D numeric integration step. Evaluation may be slow, sufficient numeric precision for fitting & minimization is not guaranteed" << endl ;
  }

  _restartNumIntEngine = kFALSE ;
  return kTRUE;
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooRealIntegral::RooRealIntegral(const RooRealIntegral& other, const char* name) : 
  RooAbsReal(other,name), 
  _valid(other._valid),
  _respectCompSelect(other._respectCompSelect),
  _sumList("!sumList",this,other._sumList),
  _intList("!intList",this,other._intList), 
  _anaList("!anaList",this,other._anaList),
  _jacList("!jacList",this,other._jacList),
  _facList("!facList","Variables independent of function",this,kFALSE,kTRUE),
  _function("!func",this,other._function), 
  _iconfig(other._iconfig),
  _sumCat("!sumCat",this,other._sumCat),
  _mode(other._mode),
  _intOperMode(other._intOperMode), 
  _restartNumIntEngine(kFALSE),
  _numIntEngine(0), 
  _numIntegrand(0),
  _rangeName(other._rangeName),
  _params(0),
  _cacheNum(kFALSE),
  _timeNumInt(other._timeNumInt)
{
 _funcNormSet = other._funcNormSet ? (RooArgSet*)other._funcNormSet->snapshot(kFALSE) : 0 ;

 for (const auto arg : other._facList) {
   RooAbsArg* argClone = (RooAbsArg*) arg->Clone() ;
   _facListOwned.addOwned(*argClone) ;
   _facList.add(*argClone) ;
   addServer(*argClone,kFALSE,kTRUE) ;
 }

 other._intList.snapshot(_saveInt) ;
 other._sumList.snapshot(_saveSum) ;

  // activate timing on numerical integrals
  activateTimingNumInts();

  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////

RooRealIntegral::~RooRealIntegral()
  // Destructor
{
  if (_numIntEngine) delete _numIntEngine ;
  if (_numIntegrand) delete _numIntegrand ;
  if (_funcNormSet) delete _funcNormSet ;
  if (_params) delete _params ;

  TRACE_DESTROY
}





////////////////////////////////////////////////////////////////////////////////

RooAbsReal* RooRealIntegral::createIntegral(const RooArgSet& iset, const RooArgSet* nset, const RooNumIntConfig* cfg, const char* rangeName) const 
{
  // Handle special case of no integration with default algorithm
  if (iset.getSize()==0) {
    return RooAbsReal::createIntegral(iset,nset,cfg,rangeName) ;
  }

  // Special handling of integral of integral, return RooRealIntegral that represents integral over all dimensions in one pass
  RooArgSet isetAll(iset) ;
  isetAll.add(_sumList) ;
  isetAll.add(_intList) ;
  isetAll.add(_anaList) ;
  isetAll.add(_facList) ;

  const RooArgSet* newNormSet(0) ;
  RooArgSet* tmp(0) ;
  if (nset && !_funcNormSet) {
    newNormSet = nset ;
  } else if (!nset && _funcNormSet) {
    newNormSet = _funcNormSet ;
  } else if (nset && _funcNormSet) {
    tmp = new RooArgSet ;
    tmp->add(*nset) ;
    tmp->add(*_funcNormSet,kTRUE) ;
    newNormSet = tmp ;
  } 
  RooAbsReal* ret =  _function.arg().createIntegral(isetAll,newNormSet,cfg,rangeName) ;

  if (tmp) {
    delete tmp ;
  }

  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return value of object. If the cache is clean, return the
/// cached value, otherwise recalculate on the fly and refill
/// the cache

Double_t RooRealIntegral::getValV(const RooArgSet* nset) const
{
//   // fast-track clean-cache processing
//   if (_operMode==AClean) {
//     return _value ;
//   }

  if (nset && nset!=_lastNSet) {
    ((RooAbsReal*) this)->setProxyNormSet(nset) ;    
    _lastNSet = (RooArgSet*) nset ;
  }

  if (isValueOrShapeDirtyAndClear()) {
    _value = traceEval(nset) ;
  } 

  return _value ;
}



////////////////////////////////////////////////////////////////////////////////
/// Perform the integration and return the result

Double_t RooRealIntegral::evaluate() const 
{
  GlobalSelectComponentRAII selCompRAII(_globalSelectComp || !_respectCompSelect);
  
  Double_t retVal(0) ;
  switch (_intOperMode) {    
    
  case Hybrid: 
    {
      // Cache numeric integrals in >1d expensive object cache
      RooDouble* cacheVal(0) ;
      if ((_cacheNum && _intList.getSize()>0) || _intList.getSize()>=_cacheAllNDim) {
        cacheVal = (RooDouble*) expensiveObjectCache().retrieveObject(GetName(),RooDouble::Class(),parameters())  ;
      }

      if (cacheVal) {
        retVal = *cacheVal ;
	//	cout << "using cached value of integral" << GetName() << endl ;
      } else {


        // Find any function dependents that are AClean
        // and switch them temporarily to ADirty
        Bool_t origState = inhibitDirty() ;
        setDirtyInhibit(kTRUE) ;

        // try to initialize our numerical integration engine
        if(!(_valid= initNumIntegrator())) {
          coutE(Integration) << ClassName() << "::" << GetName()
                             << ":evaluate: cannot initialize numerical integrator" << endl;
          return 0;
        }

        // Save current integral dependent values
        _saveInt = _intList ;
        _saveSum = _sumList ;

        // Evaluate sum/integral
        retVal = sum() ;


        // This must happen BEFORE restoring dependents, otherwise no dirty state propagation in restore step
        setDirtyInhibit(origState) ;

        // Restore integral dependent values
        _intList=_saveInt ;
        _sumList=_saveSum ;

        // Cache numeric integrals in >1d expensive object cache
        if ((_cacheNum && _intList.getSize()>0) || _intList.getSize()>=_cacheAllNDim) {
          RooDouble* val = new RooDouble(retVal) ;
          expensiveObjectCache().registerObject(_function.arg().GetName(),GetName(),*val,parameters())  ;
          //  	  cout << "### caching value of integral" << GetName() << " in " << &expensiveObjectCache() << endl ;
        }

      }
      break ;
    }
  case Analytic:
    {
      retVal =  ((RooAbsReal&)_function.arg()).analyticalIntegralWN(_mode,_funcNormSet,RooNameReg::str(_rangeName)) / jacobianProduct() ;
      cxcoutD(Tracing) << "RooRealIntegral::evaluate_analytic(" << GetName() 
		       << ")func = " << _function.arg().IsA()->GetName() << "::" << _function.arg().GetName()
		       << " raw = " << retVal << " _funcNormSet = " << (_funcNormSet?*_funcNormSet:RooArgSet()) << endl ;

      
      break ;
    }

  case PassThrough:
    {
      //setDirtyInhibit(kTRUE) ;
      retVal= _function.arg().getVal(_funcNormSet) ;      
      //setDirtyInhibit(kFALSE) ;
      break ;
    }
  }
  

  // Multiply answer with integration ranges of factorized variables
  if (_facList.getSize()>0) {
    for (const auto arg : _facList) {
      // Multiply by fit range for 'real' dependents
      if (arg->IsA()->InheritsFrom(RooAbsRealLValue::Class())) {
        RooAbsRealLValue* argLV = (RooAbsRealLValue*)arg ;
        retVal *= (argLV->getMax() - argLV->getMin()) ;
      }
      // Multiply by number of states for category dependents
      if (arg->IsA()->InheritsFrom(RooAbsCategoryLValue::Class())) {
        RooAbsCategoryLValue* argLV = (RooAbsCategoryLValue*)arg ;
        retVal *= argLV->numTypes() ;
      }    
    } 
  }


  if (dologD(Tracing)) {
    cxcoutD(Tracing) << "RooRealIntegral::evaluate(" << GetName() << ") anaInt = " << _anaList << " numInt = " << _intList << _sumList << " mode = " ;
    switch(_mode) {
    case Hybrid: ccoutD(Tracing) << "Hybrid" ; break ;
    case Analytic: ccoutD(Tracing) << "Analytic" ; break ;
    case PassThrough: ccoutD(Tracing) << "PassThrough" ; break ;
    }

    ccxcoutD(Tracing) << "raw*fact = " << retVal << endl ;
  }

  return retVal ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return product of jacobian terms originating from analytical integration

Double_t RooRealIntegral::jacobianProduct() const 
{
  if (_jacList.getSize()==0) {
    return 1 ;
  }

  Double_t jacProd(1) ;
  for (const auto elm : _jacList) {
    auto arg = static_cast<const RooAbsRealLValue*>(elm);
    jacProd *= arg->jacobian() ;
  }

  // Take fabs() here: if jacobian is negative, min and max are swapped and analytical integral
  // will be positive, so must multiply with positive jacobian.
  return fabs(jacProd) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Perform summation of list of category dependents to be integrated

Double_t RooRealIntegral::sum() const
{  
  if (_sumList.getSize()!=0) {
    // Add integrals for all permutations of categories summed over
    Double_t total(0) ;

    RooSuperCategory* sumCat = (RooSuperCategory*) _sumCat.first() ;
    for (const auto& nameIdx : *sumCat) {
      sumCat->setIndex(nameIdx);
      if (!_rangeName || sumCat->inRange(RooNameReg::str(_rangeName))) {
        total += integrate() / jacobianProduct() ;
      }
    }

    return total ;

  } else {
    // Simply return integral 
    Double_t ret = integrate() / jacobianProduct() ;
    return ret ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Perform hybrid numerical/analytical integration over all real-valued dependents

Double_t RooRealIntegral::integrate() const
{
  if (!_numIntEngine) {
    // Trivial case, fully analytical integration
    return ((RooAbsReal&)_function.arg()).analyticalIntegralWN(_mode,_funcNormSet,RooNameReg::str(_rangeName)) ;
  } else {
    return _numIntEngine->calculate()  ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Intercept server redirects and reconfigure internal object accordingly

Bool_t RooRealIntegral::redirectServersHook(const RooAbsCollection& /*newServerList*/, 
					    Bool_t /*mustReplaceAll*/, Bool_t /*nameChange*/, Bool_t /*isRecursive*/) 
{
  _restartNumIntEngine = kTRUE ;

  autoSelectDirtyMode() ;

  // Update contents value caches for _intList and _sumList
  _saveInt.removeAll() ;
  _saveSum.removeAll() ;
  _intList.snapshot(_saveInt) ;
  _sumList.snapshot(_saveSum) ;

  // Delete parameters cache if we have one
  if (_params) {
    delete _params ;
    _params = 0 ;
  }

  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////

const RooArgSet& RooRealIntegral::parameters() const
{
  if (!_params) {
    _params = new RooArgSet("params") ;
    
    RooArgSet params ;
    for (const auto server : _serverList) {
      if (server->isValueServer(*this)) _params->add(*server) ;
    }
  }

  return *_params ;
}


////////////////////////////////////////////////////////////////////////////////
/// Check if current value is valid

Bool_t RooRealIntegral::isValidReal(Double_t /*value*/, Bool_t /*printError*/) const 
{
  return kTRUE ;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if component selection is allowed

Bool_t RooRealIntegral::getAllowComponentSelection() const {
  return _respectCompSelect;
}

////////////////////////////////////////////////////////////////////////////////
/// Set component selection to be allowed/forbidden

void RooRealIntegral::setAllowComponentSelection(Bool_t allow){
  _respectCompSelect = allow;
}

////////////////////////////////////////////////////////////////////////////////
/// Customized printing of arguments of a RooRealIntegral to more intuitively reflect the contents of the
/// integration operation

void RooRealIntegral::printMetaArgs(ostream& os) const
{

  if (intVars().getSize()!=0) {
    os << "Int " ;
  }
  os << _function.arg().GetName() ;
  if (_funcNormSet) {
    os << "_Norm" ;
    os << *_funcNormSet ;
    os << " " ;
  }
 
  // List internally integrated observables and factorizing observables as analytically integrated
  RooArgSet tmp(_anaList) ;
  tmp.add(_facList) ;
  if (tmp.getSize()>0) {
    os << "d[Ana]" ;
    os << tmp ;
    os << " " ;
  }
  
  // List numerically integrated and summed observables as numerically integrated
  RooArgSet tmp2(_intList) ;
  tmp2.add(_sumList) ;
  if (tmp2.getSize()>0) {
    os << " d[Num]" ;
    os << tmp2 ;
    os << " " ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Print the state of this object to the specified output stream.

void RooRealIntegral::printMultiline(ostream& os, Int_t contents, Bool_t verbose, TString indent) const
{
  RooAbsReal::printMultiline(os,contents,verbose,indent) ;
  os << indent << "--- RooRealIntegral ---" << endl; 
  os << indent << "  Integrates ";
  _function.arg().printStream(os,kName|kArgs,kSingleLine,indent);
  TString deeper(indent);
  deeper.Append("  ");
  os << indent << "  operating mode is " 
     << (_intOperMode==Hybrid?"Hybrid":(_intOperMode==Analytic?"Analytic":"PassThrough")) << endl ;
  os << indent << "  Summed discrete args are " << _sumList << endl ;
  os << indent << "  Numerically integrated args are " << _intList << endl;
  os << indent << "  Analytically integrated args using mode " << _mode << " are " << _anaList << endl ;
  os << indent << "  Arguments included in Jacobian are " << _jacList << endl ;
  os << indent << "  Factorized arguments are " << _facList << endl ;
  os << indent << "  Function normalization set " ;
  if (_funcNormSet) 
    _funcNormSet->Print("1") ; 
  else
    os << "<none>" ;
  
  os << endl ;
} 



////////////////////////////////////////////////////////////////////////////////
/// Global switch to cache all integral values that integrate at least ndim dimensions numerically 

void RooRealIntegral::setCacheAllNumeric(Int_t ndim) {
  _cacheAllNDim = ndim ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return minimum dimensions of numeric integration for which values are cached. 

Int_t RooRealIntegral::getCacheAllNumeric() 
{
  return _cacheAllNDim ;
}

void RooRealIntegral::setNumIntTiming(Bool_t flag) {
  Int_t numIntDim = this->numIntRealVars().getSize();
  // .. and activate timing if numeric integration occurs
  if (numIntDim > 0) {
    _timeNumInt = flag;
  }
}