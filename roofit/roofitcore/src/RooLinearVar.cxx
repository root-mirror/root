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
/// \class RooLinearVar
/// RooLinearVar is the most general form of a derived real-valued object that can
/// be used by RooRealIntegral to integrate over. The requirements for this are
/// * Can be modified directly (i.e. invertible formula)
/// * Jacobian term in integral is constant (but not necessarily 1)
///
/// This class implements the most general form that satisfies these requirements
/// \f[
///    RLV = \mathrm{slope} \cdot x + \mathrm{offset}
/// \f]
/// \f$ x \f$ is required to be a RooRealVar to meet the invertibility criterium,
/// `slope` and `offset` are RooAbsReals, but cannot overlap with \f$ x \f$,
/// *i.e.*, \f$ x \f$ may not be a server of `slope` and `offset`.
///
/// In the context of a dataset, `slope` may not contain any real-valued dependents
/// (to satisfyt the constant Jacobian requirement). This check cannot be enforced at
/// construction time, but can be performed at run time through the isJacobianOK(depList)
/// member function.
///
///

#include "RooFit.h"

#include <cmath>

#include "TClass.h"
#include "RooLinearVar.h"
#include "RooStreamParser.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooNumber.h"
#include "RooBinning.h"
#include "RooMsgService.h"



using namespace std;

ClassImp(RooLinearVar);


////////////////////////////////////////////////////////////////////////////////
/// Constructor with RooAbsRealLValue variable and RooAbsReal slope and offset

RooLinearVar::RooLinearVar(const char *name, const char *title, RooAbsRealLValue& variable, 
			   const RooAbsReal& slope, const RooAbsReal& offs, const char *unit) :
  RooAbsRealLValue(name, title, unit), 
  _binning(variable.getBinning(),slope.getVal(),offs.getVal()),
  _var("var","variable",this,variable,kTRUE,kTRUE),
  _slope("slope","slope",this,(RooAbsReal&)slope),
  _offset("offset","offset",this,(RooAbsReal&)offs)
{
  // Slope and offset may not depend on variable
  if (slope.dependsOnValue(variable) || offs.dependsOnValue(variable)) {
    coutE(InputArguments) << "RooLinearVar::RooLinearVar(" << GetName() 
			  << "): ERROR, slope(" << slope.GetName() << ") and offset(" 
			  << offs.GetName() << ") may not depend on variable(" 
			  << variable.GetName() << ")" << endl ;
    assert(0) ;
  }

  // Initial plot range and number of bins from dependent variable
//   setPlotRange(variable.getPlotMin()*_slope + _offset,
//                variable.getPlotMax()*_slope + _offset) ;
//   setPlotBins(variable.getPlotBins()) ;
	       
}  



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooLinearVar::RooLinearVar(const RooLinearVar& other, const char* name) :
  RooAbsRealLValue(other,name), 
  _binning(other._binning),
  _var("var",this,other._var),
  _slope("slope",this,other._slope),
  _offset("offset",this,other._offset)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooLinearVar::~RooLinearVar() 
{
  _altBinning.Delete() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate current value of this object  

Double_t RooLinearVar::evaluate() const
{
  return _offset + _var * _slope ;
}



////////////////////////////////////////////////////////////////////////////////
/// Assign given value to linear transformation: sets input variable to (value-offset)/slope
/// If slope is zerom an error message is printed and no assignment is made

void RooLinearVar::setVal(Double_t value) 
{
  //cout << "RooLinearVar::setVal(" << GetName() << "): new value = " << value << endl ;

  // Prevent DIV0 problems
  if (_slope == 0.) {
    coutE(Eval) << "RooLinearVar::setVal(" << GetName() << "): ERROR: slope is zero, cannot invert relation" << endl ;
    return ;
  }

  // Invert formula 'value = offset + slope*var'
  ((RooRealVar&)_var.arg()).setVal((value - _offset) / _slope) ;

}



////////////////////////////////////////////////////////////////////////////////
/// Returns true if Jacobian term associated with current
/// expression tree is indeed constant.

Bool_t RooLinearVar::isJacobianOK(const RooArgSet& depList) const
{
  if (!((RooAbsRealLValue&)_var.arg()).isJacobianOK(depList)) {
    return kFALSE ;
  }

  // Check if jacobian has no real-valued dependents
  RooAbsArg* arg ;
  TIterator* dIter = depList.createIterator() ;
  while ((arg=(RooAbsArg*)dIter->Next())) {
    if (arg->IsA()->InheritsFrom(RooAbsReal::Class())) {
      if (_slope.arg().dependsOnValue(*arg)) {
// 	cout << "RooLinearVar::isJacobianOK(" << GetName() << ") return kFALSE because slope depends on value of " << arg->GetName() << endl ;
	return kFALSE ;
      }
    }
  }
  delete dIter ;
//   cout << "RooLinearVar::isJacobianOK(" << GetName() << ") return kTRUE" << endl ;
  return kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return value of Jacobian associated with the transformation

Double_t RooLinearVar::jacobian() const 
{
  return _slope*((RooAbsRealLValue&)_var.arg()).jacobian() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Read object contents from stream

Bool_t RooLinearVar::readFromStream(istream& /*is*/, Bool_t /*compact*/, Bool_t /*verbose*/) 
{
  return kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Write object contents to stream

void RooLinearVar::writeToStream(ostream& os, Bool_t compact) const
{
  if (compact) {
    os << getVal() ;
  } else {
    os << _slope.arg().GetName() << " * " << _var.arg().GetName() << " + " << _offset.arg().GetName() ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Retrieve binning of this linear transformation. A RooLinearVar does not have its own
/// binnings but uses linearly transformed binnings of the input variable. If a given
/// binning exists on the input variable, it will also exist on this linear transformation,
/// and a binning adaptor object is created on the fly.

 RooAbsBinning& RooLinearVar::getBinning(const char* name, Bool_t verbose, Bool_t createOnTheFly) 
{
  // Normalization binning
  if (name==0) {
    _binning.updateInput(((RooAbsRealLValue&)_var.arg()).getBinning(),_slope,_offset) ;
    return _binning ;
  } 

  // Alternative named range binnings, look for existing translator binning first
  RooLinTransBinning* altBinning = (RooLinTransBinning*) _altBinning.FindObject(name) ;
  if (altBinning) {
    altBinning->updateInput(((RooAbsRealLValue&)_var.arg()).getBinning(name,verbose),_slope,_offset) ;
    return *altBinning ;
  }

  // If binning is not found return default binning, if creation is not requested
  if (!_var.arg().hasRange(name) && !createOnTheFly) {
    return _binning ;
  }

  // Create translator binning on the fly
  RooAbsBinning& sourceBinning = ((RooAbsRealLValue&)_var.arg()).getBinning(name,verbose) ;
  RooLinTransBinning* transBinning = new RooLinTransBinning(sourceBinning,_slope,_offset) ;
  _altBinning.Add(transBinning) ;

  return *transBinning ;
}


////////////////////////////////////////////////////////////////////////////////
/// Const version of getBinning()

const RooAbsBinning& RooLinearVar::getBinning(const char* name, Bool_t verbose, Bool_t createOnTheFly) const
{
  return const_cast<RooLinearVar*>(this)->getBinning(name,verbose,createOnTheFly) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a list of all binning names. An empty name implies the default binning.
/// A 0 pointer should be passed to getBinning in this case.

std::list<std::string> RooLinearVar::getBinningNames() const
{
  std::list<std::string> binningNames(1, "");

  RooFIter iter = _altBinning.fwdIterator();
  const RooAbsArg* binning = 0;
  while((binning = iter.next())) {
    const char* name = binning->GetName();
    binningNames.push_back(name);
  }

  return binningNames;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if binning with given name exists.If a given binning
/// exists on the input variable, it will also exists on this linear
/// transformation.

Bool_t RooLinearVar::hasBinning(const char* name) const 
{
  return ((RooAbsRealLValue&)_var.arg()).hasBinning(name) ;
}
