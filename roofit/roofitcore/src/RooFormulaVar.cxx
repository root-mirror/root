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

//////////////////////////////////////////////////////////////////////////////
/// \class RooFormulaVar
///
/// A RooFormulaVar is a generic implementation of a real-valued object,
/// which takes a RooArgList of servers and a C++ expression string defining how
/// its value should be calculated from the given list of servers.
/// RooFormulaVar uses a RooFormula object to perform the expression evaluation.
///
/// If RooAbsPdf objects are supplied to RooFormulaVar as servers, their
/// raw (unnormalized) values will be evaluated. Use RooGenericPdf, which
/// constructs generic PDF functions, to access their properly normalized
/// values.
///
/// The string expression can be any valid TFormula expression referring to the
/// listed servers either by name or by their ordinal list position. These three are
/// equivalent:
/// ```
///   RooFormulaVar("gen", "x*y", RooArgList(x,y))       // reference by name
///   RooFormulaVar("gen", "@0*@1", RooArgList(x,y))     // reference by ordinal with @
///   RooFormulaVar("gen", "x[0]*x[1]", RooArgList(x,y)) // TFormula-builtin reference by ordinal
/// ```
/// Note that `x[i]` is an expression reserved for TFormula. All variable references
/// are automatically converted to the TFormula-native format. If a variable with
/// the name `x` is given, the RooFormula interprets `x[i]` as a list position,
/// but `x` without brackets as the name of a RooFit object.
///
/// The last two versions, while slightly less readable, are more versatile because
/// the names of the arguments are not hard coded.
///


#include "RooFit.h"
#include "Riostream.h"

#include "RooFormulaVar.h"
#include "RooStreamParser.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMsgService.h"
#include "RooTrace.h"


using namespace std;

ClassImp(RooFormulaVar);



////////////////////////////////////////////////////////////////////////////////
/// Constructor with formula expression and list of input variables.
/// \param[in] name Name of the formula.
/// \param[in] title Title of the formula.
/// \param[in] inFormula Expression to be evaluated.
/// \param[in] dependents Variables that should be passed to the formula.
/// \param[in] checkVariables Check that all variables from `dependents` are used in the expression.
RooFormulaVar::RooFormulaVar(const char *name, const char *title, const char* inFormula, const RooArgList& dependents,
    bool checkVariables) :
  RooAbsReal(name,title), 
  _actualVars("actualVars","Variables used by formula expression",this),
  _formExpr(inFormula)
{  
  _actualVars.add(dependents) ; 

  if (_actualVars.getSize()==0) {
    _value = traceEval(0);
  } else {
    _formula.reset(new RooFormula(GetName(), _formExpr, _actualVars, checkVariables));
    _formExpr = _formula->formulaString().c_str();
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor with formula expression, title and list of input variables.
/// \param[in] name Name of the formula.
/// \param[in] title Formula expression. Will also be used as the title.
/// \param[in] dependents Variables that should be passed to the formula.
/// \param[in] checkVariables Check that all variables from `dependents` are used in the expression.
RooFormulaVar::RooFormulaVar(const char *name, const char *title, const RooArgList& dependents,
    bool checkVariables) :
  RooAbsReal(name,title),
  _actualVars("actualVars","Variables used by formula expression",this),
  _formExpr(title)
{  
  _actualVars.add(dependents) ; 

  if (_actualVars.getSize()==0) {
    _value = traceEval(0);
  } else {
    _formula.reset(new RooFormula(GetName(), _formExpr, _actualVars, checkVariables));
    _formExpr = _formula->formulaString().c_str();
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooFormulaVar::RooFormulaVar(const RooFormulaVar& other, const char* name) : 
  RooAbsReal(other, name), 
  _actualVars("actualVars",this,other._actualVars),
  _formExpr(other._formExpr)
{
  if (other._formula && other._formula->ok()) {
    _formula.reset(new RooFormula(GetName(), _formExpr, _actualVars, /*checkVariables=*/false));
    _formExpr = _formula->formulaString().c_str();
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Return reference to internal RooFormula object.
/// If it doesn't exist, create it on the fly.
RooFormula& RooFormulaVar::getFormula() const
{
  if (!_formula) {
    // After being read from file, the formula object might not exist, yet:
    auto theFormula = new RooFormula(GetName(), _formExpr, _actualVars);
    const_cast<std::unique_ptr<RooFormula>&>(this->_formula).reset(theFormula);
    const_cast<TString&>(_formExpr) = _formula->formulaString().c_str();
  }

  return *_formula;
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate current value of object from internal formula

Double_t RooFormulaVar::evaluate() const
{
  return getFormula().eval(_lastNSet);
}


////////////////////////////////////////////////////////////////////////////////
/// Evaluate the formula for all entries of our servers found in `inputData`.
RooSpan<double> RooFormulaVar::evaluateSpan(RooBatchCompute::RunContext& inputData, const RooArgSet* normSet) const {
  if (normSet != _lastNSet) {
    // TODO: Remove dependence on _lastNSet
    // See also comment in RooAbsReal::getValBatch().
    std::cerr << "Formula " << GetName() << " " << GetTitle() << "\n\tBeing evaluated with normSet " << normSet << "\n";
    normSet->Print("V");
    std::cerr << "\tHowever, _lastNSet = " << _lastNSet << "\n";
    if (_lastNSet) _lastNSet->Print("V");

    throw std::logic_error("Got conflicting norm sets. This shouldn't happen.");
  }

  return formula().evaluateSpan(this, inputData, normSet);
}


////////////////////////////////////////////////////////////////////////////////
/// Propagate server change information to embedded RooFormula object

Bool_t RooFormulaVar::redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t /*isRecursive*/)
{
  bool success = getFormula().changeDependents(newServerList,mustReplaceAll,nameChange);

  _formExpr = getFormula().GetTitle();
  return success;
}



////////////////////////////////////////////////////////////////////////////////
/// Print info about this object to the specified stream.   

void RooFormulaVar::printMultiline(ostream& os, Int_t contents, Bool_t verbose, TString indent) const
{
  RooAbsReal::printMultiline(os,contents,verbose,indent);
  if(verbose) {
    indent.Append("  ");
    os << indent;
    getFormula().printMultiline(os,contents,verbose,indent);
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Add formula expression as meta argument in printing interface

void RooFormulaVar::printMetaArgs(ostream& os) const 
{
  os << "formula=\"" << _formExpr << "\" " ;
}




////////////////////////////////////////////////////////////////////////////////
/// Read object contents from given stream

Bool_t RooFormulaVar::readFromStream(istream& /*is*/, Bool_t /*compact*/, Bool_t /*verbose*/)
{
  coutE(InputArguments) << "RooFormulaVar::readFromStream(" << GetName() << "): can't read" << endl ;
  return kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Write object contents to given stream

void RooFormulaVar::writeToStream(ostream& os, Bool_t compact) const
{
  if (compact) {
    cout << getVal() << endl ;
  } else {
    os << GetTitle() ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Forward the plot sampling hint from the p.d.f. that defines the observable obs  

std::list<Double_t>* RooFormulaVar::binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const
{
  for (const auto par : _actualVars) {
    auto func = static_cast<const RooAbsReal*>(par);
    list<Double_t>* binb = nullptr;

    if (func && (binb = func->binBoundaries(obs,xlo,xhi)) ) {
      return binb;
    }
  }
  
  return nullptr;
}



////////////////////////////////////////////////////////////////////////////////
/// Forward the plot sampling hint from the p.d.f. that defines the observable obs  

std::list<Double_t>* RooFormulaVar::plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const
{
  for (const auto par : _actualVars) {
    auto func = dynamic_cast<const RooAbsReal*>(par);
    list<Double_t>* hint = nullptr;

    if (func && (hint = func->plotSamplingHint(obs,xlo,xhi)) ) {
      return hint;
    }
  }

  return nullptr;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the default error level for MINUIT error analysis
/// If the formula contains one or more RooNLLVars and 
/// no RooChi2Vars, return the defaultErrorLevel() of
/// RooNLLVar. If the addition contains one ore more RooChi2Vars
/// and no RooNLLVars, return the defaultErrorLevel() of
/// RooChi2Var. If the addition contains neither or both
/// issue a warning message and return a value of 1

Double_t RooFormulaVar::defaultErrorLevel() const 
{
  RooAbsReal* nllArg(0) ;
  RooAbsReal* chi2Arg(0) ;

  for (const auto arg : _actualVars) {
    if (dynamic_cast<RooNLLVar*>(arg)) {
      nllArg = (RooAbsReal*)arg ;
    }
    if (dynamic_cast<RooChi2Var*>(arg)) {
      chi2Arg = (RooAbsReal*)arg ;
    }
  }

  if (nllArg && !chi2Arg) {
    coutI(Minimization) << "RooFormulaVar::defaultErrorLevel(" << GetName() 
			<< ") Formula contains a RooNLLVar, using its error level" << endl ;
    return nllArg->defaultErrorLevel() ;
  } else if (chi2Arg && !nllArg) {
    coutI(Minimization) << "RooFormulaVar::defaultErrorLevel(" << GetName() 
	 << ") Formula contains a RooChi2Var, using its error level" << endl ;
    return chi2Arg->defaultErrorLevel() ;
  } else if (!nllArg && !chi2Arg) {
    coutI(Minimization) << "RooFormulaVar::defaultErrorLevel(" << GetName() << ") WARNING: "
		      << "Formula contains neither RooNLLVar nor RooChi2Var server, using default level of 1.0" << endl ;
  } else {
    coutI(Minimization) << "RooFormulaVar::defaultErrorLevel(" << GetName() << ") WARNING: "
			<< "Formula contains BOTH RooNLLVar and RooChi2Var server, using default level of 1.0" << endl ;
  }

  return 1.0 ;
}




