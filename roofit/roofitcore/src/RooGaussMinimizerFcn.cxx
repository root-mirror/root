/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   AL, Alfio Lazzaro,   INFN Milan,        alfio.lazzaro@mi.infn.it        *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *   VC, Vince Croft,     DIANA / NYU,        vincent.croft@cern.ch          *
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef __ROOFIT_NOROOMINIMIZER

//////////////////////////////////////////////////////////////////////////////
//
// RooGaussMinimizerFcn is am interface class to the ROOT::Math function
// for minization. See RooGaussMinimizer.cxx for more information.
//                                                                                   

#include <iostream>

#include "RooFit.h"
#include "RooMinimizerFcn.h"

#include "Riostream.h"

#include "TIterator.h"
#include "TClass.h"

#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooAbsRealLValue.h"
#include "RooMsgService.h"

#include "RooMinimizer.h"
#include "RooGaussMinimizer.h"

using namespace std;

RooGaussMinimizerFcn::RooGaussMinimizerFcn(RooAbsReal *funct, RooGaussMinimizer* context,
			   bool verbose) :
  _funct(funct), _context(context),
  // Reset the *largest* negative log-likelihood value we have seen so far
  _maxFCN(-1e30), _numBadNLL(0),  
  _printEvalErrors(10), _doEvalErrorWall(kTRUE),
  _nDim(0), _logfile(0),
  _verbose(verbose)
{ 

  _evalCounter = 0 ;
  
  // Examine parameter list
  RooArgSet* paramSet = _funct->getParameters(RooArgSet());
  RooArgList paramList(*paramSet);
  delete paramSet;

  _floatParamList = (RooArgList*) paramList.selectByAttrib("Constant",kFALSE);
  if (_floatParamList->getSize()>1) {
    _floatParamList->sort();
  }
  _floatParamList->setName("floatParamList");

  _constParamList = (RooArgList*) paramList.selectByAttrib("Constant",kTRUE);
  if (_constParamList->getSize()>1) {
    _constParamList->sort();
  }
  _constParamList->setName("constParamList");

  // Remove all non-RooRealVar parameters from list (MINUIT cannot handle them)
  TIterator* pIter = _floatParamList->createIterator();
  RooAbsArg* arg;
  while ((arg=(RooAbsArg*)pIter->Next())) {
    if (!arg->IsA()->InheritsFrom(RooAbsRealLValue::Class())) {
      oocoutW(_context,Minimization) << "RooGaussMinimizerFcn::RooGaussMinimizerFcn: removing parameter " 
				     << arg->GetName()
                                     << " from list because it is not of type RooRealVar" << endl;
      _floatParamList->remove(*arg);
    }
  }
  delete pIter;

  _nDim = _floatParamList->getSize();

  updateFloatVec() ;
  
  // Save snapshot of initial lists
  _initFloatParamList = (RooArgList*) _floatParamList->snapshot(kFALSE) ;
  _initConstParamList = (RooArgList*) _constParamList->snapshot(kFALSE) ;

}



RooGaussMinimizerFcn::RooGaussMinimizerFcn(const RooGaussMinimizerFcn& other) : ROOT::Math::IMultiGradFunction(other), 
  _evalCounter(other._evalCounter),
  _funct(other._funct),
  _context(other._context),
  _maxFCN(other._maxFCN),
  _numBadNLL(other._numBadNLL),
  _printEvalErrors(other._printEvalErrors),
  _doEvalErrorWall(other._doEvalErrorWall),
  _nDim(other._nDim),
  _logfile(other._logfile),
  _verbose(other._verbose),
  _floatParamVec(other._floatParamVec)
{  
  _floatParamList = new RooArgList(*other._floatParamList) ;
  _constParamList = new RooArgList(*other._constParamList) ;
  _initFloatParamList = (RooArgList*) other._initFloatParamList->snapshot(kFALSE) ;
  _initConstParamList = (RooArgList*) other._initConstParamList->snapshot(kFALSE) ;  
}


RooGaussMinimizerFcn::~RooGaussMinimizerFcn()
{
  delete _floatParamList;
  delete _initFloatParamList;
  delete _constParamList;
  delete _initConstParamList;
}


ROOT::Math::IMultiGradFunction* RooGaussMinimizerFcn::Clone() const 
{  
  return new RooGaussMinimizerFcn(*this) ;
}


Bool_t RooGaussMinimizerFcn::Synchronize(std::vector<ROOT::Fit::ParameterSettings>& parameters, 
				 Bool_t optConst, Bool_t verbose)
{

  // Internal function to synchronize TMinimizer with current
  // information in RooAbsReal function parameters
  
  Bool_t constValChange(kFALSE) ;
  Bool_t constStatChange(kFALSE) ;
  
  Int_t index(0) ;
  
  // Handle eventual migrations from constParamList -> floatParamList
  for(index= 0; index < _constParamList->getSize() ; index++) {

    RooRealVar *par= dynamic_cast<RooRealVar*>(_constParamList->at(index)) ;
    if (!par) continue ;

    RooRealVar *oldpar= dynamic_cast<RooRealVar*>(_initConstParamList->at(index)) ;
    if (!oldpar) continue ;

    // Test if constness changed
    if (!par->isConstant()) {      
    
      // Remove from constList, add to floatList
      _constParamList->remove(*par) ;
      _floatParamList->add(*par) ;
      _initFloatParamList->addClone(*oldpar) ;      
      _initConstParamList->remove(*oldpar) ;
      constStatChange=kTRUE ;
      _nDim++ ;

      if (verbose) {
	oocoutI(_context,Minimization) << "RooGaussMinimizerFcn::synchronize: parameter " 
				     << par->GetName() << " is now floating." << endl ;
      }
    } 

    // Test if value changed
    if (par->getVal()!= oldpar->getVal()) {
      constValChange=kTRUE ;      
      if (verbose) {
	oocoutI(_context,Minimization) << "RooGaussMinimizerFcn::synchronize: value of constant parameter " 
				       << par->GetName() 
				       << " changed from " << oldpar->getVal() << " to " 
				       << par->getVal() << endl ;
      }
    }

  }

  // Update reference list
  *_initConstParamList = *_constParamList ;
  
  // Synchronize MINUIT with function state
  // Handle floatParamList
  for(index= 0; index < _floatParamList->getSize(); index++) {
    RooRealVar *par= dynamic_cast<RooRealVar*>(_floatParamList->at(index)) ;
    
    if (!par) continue ;

    Double_t pstep(0) ;
    Double_t pmin(0) ;
    Double_t pmax(0) ;

    if(!par->isConstant()) {

      // Verify that floating parameter is indeed of type RooRealVar 
      if (!par->IsA()->InheritsFrom(RooRealVar::Class())) {
	oocoutW(_context,Minimization) << "RooGaussMinimizerFcn::fit: Error, non-constant parameter " 
				       << par->GetName() 
				       << " is not of type RooRealVar, skipping" << endl ;
	_floatParamList->remove(*par);
	index--;
	_nDim--;
	continue ;
      }

      // Set the limits, if not infinite
      if (par->hasMin() )
	pmin = par->getMin();
      if (par->hasMax() )
	pmax = par->getMax();
      
      // Calculate step size
      pstep = par->getError();
      if(pstep <= 0) {
	// Floating parameter without error estitimate
	if (par->hasMin() && par->hasMax()) {
	  pstep= 0.1*(pmax-pmin);

	  // Trim default choice of error if within 2 sigma of limit
	  if (pmax - par->getVal() < 2*pstep) {
	    pstep = (pmax - par->getVal())/2 ;
	  } else if (par->getVal() - pmin < 2*pstep) {
	    pstep = (par->getVal() - pmin )/2 ;	    
	  }	  

	  // If trimming results in zero error, restore default
	  if (pstep==0) {
	    pstep= 0.1*(pmax-pmin);
	  }

	} else {
	  pstep=1 ;
	}						  
	if(verbose) {
	  oocoutW(_context,Minimization) << "RooGaussMinimizerFcn::synchronize: WARNING: no initial error estimate available for "
					 << par->GetName() << ": using " << pstep << endl;
	}
      }       
    } else {
      pmin = par->getVal() ;
      pmax = par->getVal() ;      
    }

    // new parameter
    if (index>=Int_t(parameters.size())) {

      if (par->hasMin() && par->hasMax()) {
	parameters.push_back(ROOT::Fit::ParameterSettings(par->GetName(),
							  par->getVal(),
							  pstep,
							  pmin,pmax));
      }
      else {
	parameters.push_back(ROOT::Fit::ParameterSettings(par->GetName(),
							  par->getVal(),
							  pstep));
        if (par->hasMin() ) 
           parameters.back().SetLowerLimit(pmin);
        else if (par->hasMax() ) 
           parameters.back().SetUpperLimit(pmax);
      }
      
      continue;

    }

    Bool_t oldFixed = parameters[index].IsFixed();
    Double_t oldVar = parameters[index].Value();
    Double_t oldVerr = parameters[index].StepSize();
    Double_t oldVlo = parameters[index].LowerLimit();
    Double_t oldVhi = parameters[index].UpperLimit();

    if (par->isConstant() && !oldFixed) {

      // Parameter changes floating -> constant : update only value if necessary
      if (oldVar!=par->getVal()) {
	parameters[index].SetValue(par->getVal());
	if (verbose) {
	  oocoutI(_context,Minimization) << "RooGaussMinimizerFcn::synchronize: value of parameter " 
					 << par->GetName() << " changed from " << oldVar 
					 << " to " << par->getVal() << endl ;
	}
      }
      parameters[index].Fix();
      constStatChange=kTRUE ;
      if (verbose) {
	oocoutI(_context,Minimization) << "RooGaussMinimizerFcn::synchronize: parameter " 
				       << par->GetName() << " is now fixed." << endl ;
      }

    } else if (par->isConstant() && oldFixed) {
      
      // Parameter changes constant -> constant : update only value if necessary
      if (oldVar!=par->getVal()) {
	parameters[index].SetValue(par->getVal());
	constValChange=kTRUE ;

	if (verbose) {
	  oocoutI(_context,Minimization) << "RooGaussMinimizerFcn::synchronize: value of fixed parameter " 
					 << par->GetName() << " changed from " << oldVar 
					 << " to " << par->getVal() << endl ;
	}
      }

    } else {
      // Parameter changes constant -> floating
      if (!par->isConstant() && oldFixed) {
	parameters[index].Release();
	constStatChange=kTRUE ;
	
	if (verbose) {
	  oocoutI(_context,Minimization) << "RooGaussMinimizerFcn::synchronize: parameter " 
					 << par->GetName() << " is now floating." << endl ;
	}
      } 

      // Parameter changes constant -> floating : update all if necessary      
      if (oldVar!=par->getVal() || oldVlo!=pmin || oldVhi != pmax || oldVerr!=pstep) {
	parameters[index].SetValue(par->getVal()); 
	parameters[index].SetStepSize(pstep);
        if (par->hasMin() && par->hasMax() ) 
           parameters[index].SetLimits(pmin,pmax);  
        else if (par->hasMin() )
           parameters[index].SetLowerLimit(pmin);  
        else if (par->hasMax() )
           parameters[index].SetUpperLimit(pmax);  
      }

      // Inform user about changes in verbose mode
      if (verbose) {
	// if ierr<0, par was moved from the const list and a message was already printed

	if (oldVar!=par->getVal()) {
	  oocoutI(_context,Minimization) << "RooGaussMinimizerFcn::synchronize: value of parameter " 
					 << par->GetName() << " changed from " << oldVar << " to " 
					 << par->getVal() << endl ;
	}
	if (oldVlo!=pmin || oldVhi!=pmax) {
	  oocoutI(_context,Minimization) << "RooGaussMinimizerFcn::synchronize: limits of parameter " 
					 << par->GetName() << " changed from [" << oldVlo << "," << oldVhi 
					 << "] to [" << pmin << "," << pmax << "]" << endl ;
	}

	// If oldVerr=0, then parameter was previously fixed
	if (oldVerr!=pstep && oldVerr!=0) {
	  oocoutI(_context,Minimization) << "RooGaussMinimizerFcn::synchronize: error/step size of parameter " 
					 << par->GetName() << " changed from " << oldVerr << " to " << pstep << endl ;
	}
      }      
    }
  }

  if (optConst) {
    if (constStatChange) {

      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors) ;

      oocoutI(_context,Minimization) << "RooGaussMinimizerFcn::synchronize: set of constant parameters changed, rerunning const optimizer" << endl ;
      _funct->constOptimizeTestStatistic(RooAbsArg::ConfigChange) ;
    } else if (constValChange) {
      oocoutI(_context,Minimization) << "RooGaussMinimizerFcn::synchronize: constant parameter values changed, rerunning const optimizer" << endl ;
      _funct->constOptimizeTestStatistic(RooAbsArg::ValueChange) ;
    }
    
    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors) ;  

  }

  updateFloatVec() ;

  return 0 ;  

}

Double_t RooGaussMinimizerFcn::GetPdfParamVal(Int_t index)
{
  // Access PDF parameter value by ordinal index (needed by MINUIT)

  return ((RooRealVar*)_floatParamList->at(index))->getVal() ;
}

Double_t RooGaussMinimizerFcn::GetPdfParamErr(Int_t index)
{
  // Access PDF parameter error by ordinal index (needed by MINUIT)
  return ((RooRealVar*)_floatParamList->at(index))->getError() ;
}


void RooGaussMinimizerFcn::SetPdfParamErr(Int_t index, Double_t value)
{
  // Modify PDF parameter error by ordinal index (needed by MINUIT)

  ((RooRealVar*)_floatParamList->at(index))->setError(value) ;
}



void RooGaussMinimizerFcn::ClearPdfParamAsymErr(Int_t index)
{
  // Modify PDF parameter error by ordinal index (needed by MINUIT)

  ((RooRealVar*)_floatParamList->at(index))->removeAsymError() ;
}


void RooGaussMinimizerFcn::SetPdfParamErr(Int_t index, Double_t loVal, Double_t hiVal)
{
  // Modify PDF parameter error by ordinal index (needed by MINUIT)

  ((RooRealVar*)_floatParamList->at(index))->setAsymError(loVal,hiVal) ;
}


void RooGaussMinimizerFcn::BackProp(const ROOT::Fit::FitResult &results)
{
  // Transfer MINUIT fit results back into RooFit objects

  for (Int_t index= 0; index < _nDim; index++) {
    Double_t value = results.Value(index);
    SetPdfParamVal(index, value);

    // Set the parabolic error    
    Double_t err = results.Error(index);
    SetPdfParamErr(index, err);
    
    Double_t eminus = results.LowerError(index);
    Double_t eplus = results.UpperError(index);

    if(eplus > 0 || eminus < 0) {
      // Store the asymmetric error, if it is available
      SetPdfParamErr(index, eminus,eplus);
    } else {
      // Clear the asymmetric error
      ClearPdfParamAsymErr(index) ;
    }
  }

}

Bool_t RooGaussMinimizerFcn::SetLogFile(const char* inLogfile) 
{
  // Change the file name for logging of a RooGaussMinimizer of all MINUIT steppings
  // through the parameter space. If inLogfile is null, the current log file
  // is closed and logging is stopped.

  if (_logfile) {
    oocoutI(_context,Minimization) << "RooGaussMinimizerFcn::setLogFile: closing previous log file" << endl ;
    _logfile->close() ;
    delete _logfile ;
    _logfile = 0 ;
  }
  _logfile = new ofstream(inLogfile) ;
  if (!_logfile->good()) {
    oocoutI(_context,Minimization) << "RooGaussMinimizerFcn::setLogFile: cannot open file " << inLogfile << endl ;
    _logfile->close() ;
    delete _logfile ;
    _logfile= 0;
  }  
  
  return kFALSE ;

}


void RooGaussMinimizerFcn::ApplyCovarianceMatrix(TMatrixDSym& V) 
{
  // Apply results of given external covariance matrix. i.e. propagate its errors
  // to all RRV parameter representations and give this matrix instead of the
  // HESSE matrix at the next save() call

  for (Int_t i=0 ; i<_nDim ; i++) {
    // Skip fixed parameters
    if (_floatParamList->at(i)->isConstant()) {
      continue ;
    }
    SetPdfParamErr(i, sqrt(V(i,i))) ;		  
  }

}


Bool_t RooGaussMinimizerFcn::SetPdfParamVal(const Int_t &index, const Double_t &value) const
{
  //RooRealVar* par = (RooRealVar*)_floatParamList->at(index);
  RooRealVar* par = (RooRealVar*)_floatParamVec[index] ;

  if (par->getVal()!=value) {
    if (_verbose) cout << par->GetName() << "=" << value << ", " ;
    
    par->setVal(value);
    return kTRUE;
  }

  return kFALSE;
}



////////////////////////////////////////////////////////////////////////////////

void RooGaussMinimizerFcn::updateFloatVec() 
{
  _floatParamVec.clear() ;
  RooFIter iter = _floatParamList->fwdIterator() ;
  RooAbsArg* arg ;
  _floatParamVec = std::vector<RooAbsArg*>(_floatParamList->getSize()) ;
  Int_t i(0) ;
  while((arg=iter.next())) {
    _floatParamVec[i++] = arg ;
  }
}



double RooGaussMinimizerFcn::DoEval(const double *x) const 
{

  // Set the parameter values for this iteration
  for (int index = 0; index < _nDim; index++) {
    if (_logfile && _logfile->is_open()) (*_logfile) << x[index] << " " ;
    SetPdfParamVal(index,x[index]);
  }

  // Calculate the function for these parameters  
  RooAbsReal::setHideOffset(kFALSE) ;
  double fvalue = _funct->getVal();
  RooAbsReal::setHideOffset(kTRUE) ;

  if (!std::isfinite(fvalue) || RooAbsReal::numEvalErrors() > 0 || fvalue > 1e30) {

    if (_printEvalErrors>=0) {

      if (_doEvalErrorWall) {
        oocoutW(_context,Minimization) << "RooGaussMinimizerFcn: Minimized function has error status." << endl
				       << "Returning maximum FCN so far (" << _maxFCN
				       << ") to force MIGRAD to back out of this region. Error log follows" << endl ;
      } else {
        oocoutW(_context,Minimization) << "RooGaussMinimizerFcn: Minimized function has error status but is ignored" << endl ;
      }

      TIterator* iter = _floatParamList->createIterator() ;
      RooRealVar* var ;
      Bool_t first(kTRUE) ;
      ooccoutW(_context,Minimization) << "Parameter values: " ;
      while((var=(RooRealVar*)iter->Next())) {
        if (first) { first = kFALSE ; } else ooccoutW(_context,Minimization) << ", " ;
        ooccoutW(_context,Minimization) << var->GetName() << "=" << var->getVal() ;
      }
      delete iter ;
      ooccoutW(_context,Minimization) << endl ;

      RooAbsReal::printEvalErrors(ooccoutW(_context,Minimization),_printEvalErrors) ;
      ooccoutW(_context,Minimization) << endl ;
    }

    if (_doEvalErrorWall) {
      fvalue = _maxFCN+1 ;
    }

    RooAbsReal::clearEvalErrorLog() ;
    _numBadNLL++ ;
  } else if (fvalue>_maxFCN) {
    _maxFCN = fvalue ;
  }
      
  // Optional logging
  if (_logfile && _logfile->is_open())
    (*_logfile) << setprecision(15) << fvalue << setprecision(4) << endl;
  if (_verbose) {
    cout << "\nprevFCN" << (_funct->isOffsetting()?"-offset":"") << " = " << setprecision(10) 
         << fvalue << setprecision(4) << "  " ;
    cout.flush() ;
  }

  _evalCounter++ ;
  cout<< "func eval "<<fvalue<<endl;
  return fvalue;
}



double RooGaussMinimizerFcn::DoDerivative(const double *x, unsigned int icoord) const {
  // Set the parameter values for this iteration
  // EGP TODO: this is already done in DoEval as well; find efficient way to do only once
  for (int index = 0; index < _nDim; index++) {
    if (_logfile && _logfile->is_open()) (*_logfile) << x[index] << " " ;
    SetPdfParamVal(index,x[index]);
  }

  // Calculate the function for these parameters
  RooAbsReal::setHideOffset(kFALSE) ; // EGP TODO: check whether this is necessary

  ///// EGP TODO: REPLACE BELOW DERIVATIVE CALCULATION WITH THE FANCY MINUIT TYPE STUFF
  double dx = max(1e-5 * x[icoord], 1e-8);
  double fvalue_0 = _funct->getVal();

  if (_logfile) (*_logfile) << x[icoord] << " " ;
  SetPdfParamVal(icoord,x[icoord] + dx);

  double fvalue_dx = _funct->getVal();

  double derivative_i_value = (fvalue_dx - fvalue_0) / dx; //######## OI THIS IS WHERE WE COMPUTE THE GRADIENT####
  ///// EGP TODO: REPLACE ABOVE DERIVATIVE CALCULATION WITH THE FANCY MINUIT TYPE STUFF

  RooAbsReal::setHideOffset(kTRUE) ; // EGP TODO: check whether this is necessary

  // EGP TODO: decide whether to do error handling and logging, like in DoEval

  // EGP TOOO: update this when changing the derivative algorithm
  // Count the function calls necessary for this derivative and use that.
  // Except when the derivative itself calls DoEval where the counter is already updated!
  _evalCounter += 2;

  cout << "grad value " << derivative_i_value << endl;
  return derivative_i_value;
}


#endif

