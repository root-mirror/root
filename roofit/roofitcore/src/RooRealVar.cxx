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
\file RooRealVar.cxx
\class RooRealVar
\ingroup Roofitcore

RooRealVar represents a variable that can be changed from the outside.
For example by the user or a fitter.

It can be written into datasets, can hold a (possibly asymmetric) error, and
can have several ranges. These can be accessed with names, to e.g. limit fits
or integrals to sub ranges. The range without any name is used as default range.
**/

#include "RooRealVar.h"

#include "RooStreamParser.h"
#include "RooErrorVar.h"
#include "RooRangeBinning.h"
#include "RooCmdConfig.h"
#include "RooMsgService.h"
#include "RooParamBinning.h"
#include "RooVectorDataStore.h"
#include "RooTrace.h"
#include "RooRealVarSharedProperties.h"
#include "RooUniformBinning.h"
#include "RunContext.h"

#include "TTree.h"
#include "TBuffer.h"
#include "TBranch.h"
#include "snprintf.h"

using namespace std;

ClassImp(RooRealVar);


Bool_t RooRealVar::_printScientific(kFALSE) ;
Int_t  RooRealVar::_printSigDigits(5) ;

/// Return a reference to a map of weak pointers to RooRealVarSharedProperties.
std::map<std::string,std::weak_ptr<RooRealVarSharedProperties>>& RooRealVar::_sharedPropList() 
{
  static std::map<std::string,std::weak_ptr<RooRealVarSharedProperties>> sharedPropList;
  return sharedPropList; 
}

/// Return a dummy object to use when properties are not initialised.
RooRealVarSharedProperties& RooRealVar::_nullProp()
{
  static const std::unique_ptr<RooRealVarSharedProperties> nullProp(new RooRealVarSharedProperties("00000000-0000-0000-0000-000000000000"));
  return *nullProp; 
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

RooRealVar::RooRealVar()  :  _error(0), _asymErrLo(0), _asymErrHi(0), _binning(new RooUniformBinning())
{
  _fast = kTRUE ;
  TRACE_CREATE
}


////////////////////////////////////////////////////////////////////////////////
/// Create a constant variable with a value and optional unit.
RooRealVar::RooRealVar(const char *name, const char *title,
		       Double_t value, const char *unit) :
  RooAbsRealLValue(name, title, unit), _error(-1), _asymErrLo(1), _asymErrHi(-1),
  _binning(new RooUniformBinning(-1,1,100))
{
  _value = value ;
  _fast = kTRUE ;
  removeRange();
  setConstant(kTRUE) ;
  TRACE_CREATE
}


////////////////////////////////////////////////////////////////////////////////
/// Create a variable allowed to float in the given range.
/// The initial value will be set to the center of the range.
RooRealVar::RooRealVar(const char *name, const char *title,
		       Double_t minValue, Double_t maxValue,
		       const char *unit) :
  RooAbsRealLValue(name, title, unit), _error(-1), _asymErrLo(1), _asymErrHi(-1),
  _binning(new RooUniformBinning(minValue,maxValue,100))
{
  _fast = kTRUE ;

  if (RooNumber::isInfinite(minValue)) {
    if (RooNumber::isInfinite(maxValue)) {
      // [-inf,inf]
      _value = 0 ;
    } else {
      // [-inf,X]
      _value= maxValue ;
    }
  } else {
    if (RooNumber::isInfinite(maxValue)) {
      // [X,inf]
      _value = minValue ;
    } else {
      // [X,X]
      _value= 0.5*(minValue + maxValue);
    }
  }

  //   setPlotRange(minValue,maxValue) ;
  setRange(minValue,maxValue) ;
  TRACE_CREATE
}


////////////////////////////////////////////////////////////////////////////////
/// Create a variable with the given starting value. It is allowed to float
/// within the defined range. Optionally, a unit can be specified for axis labels.
RooRealVar::RooRealVar(const char *name, const char *title,
		       Double_t value, Double_t minValue, Double_t maxValue,
		       const char *unit) :
  RooAbsRealLValue(name, title, unit), _error(-1), _asymErrLo(1), _asymErrHi(-1),
  _binning(new RooUniformBinning(minValue,maxValue,100))
{
    _fast = kTRUE ;
    setRange(minValue,maxValue) ;

    Double_t clipValue ;
    inRange(value,0,&clipValue) ;
    _value = clipValue ;

    TRACE_CREATE
}


////////////////////////////////////////////////////////////////////////////////
/// Copy Constructor

RooRealVar::RooRealVar(const RooRealVar& other, const char* name) :
  RooAbsRealLValue(other,name),
  _error(other._error),
  _asymErrLo(other._asymErrLo),
  _asymErrHi(other._asymErrHi)
{
  _sharedProp = other.sharedProp();
  if (other._binning) {
     _binning.reset(other._binning->clone());
     _binning->insertHook(*this) ;
  }
  _fast = kTRUE ;

  for (const auto& item : other._altNonSharedBinning) {
    std::unique_ptr<RooAbsBinning> abc( item.second->clone() );
    abc->insertHook(*this) ;
    _altNonSharedBinning[item.first] = std::move(abc);
  }

  TRACE_CREATE

}

/// Assign the values of another RooRealVar to this instance.
RooRealVar& RooRealVar::operator=(const RooRealVar& other) {
  RooAbsRealLValue::operator=(other);

  _error = other._error;
  _asymErrLo = other._asymErrLo;
  _asymErrHi = other._asymErrHi;

  _binning.reset();
  if (other._binning) {
    _binning.reset(other._binning->clone());
    _binning->insertHook(*this) ;
  }

  _altNonSharedBinning.clear();
  for (const auto& item : other._altNonSharedBinning) {
    RooAbsBinning* abc = item.second->clone();
    _altNonSharedBinning[item.first].reset(abc);
    abc->insertHook(*this);
  }

  _sharedProp = other.sharedProp();

  return *this;
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooRealVar::~RooRealVar()
{
  TRACE_DESTROY
}


////////////////////////////////////////////////////////////////////////////////
/// Return value of variable

Double_t RooRealVar::getValV(const RooArgSet*) const
{
  return _value ;
}


////////////////////////////////////////////////////////////////////////////////
/// Retrieve data column of this variable.
/// \param inputData Struct with data arrays.
/// \param normSet Ignored.
/// 1. Check if `inputData` has a column of data registered for this variable (checks the pointer).
/// 2. If not, check if there's an object with the same name, and use this object's values.
/// 3. If there is no such object, return a batch of size one with the current value of the variable.
/// For cases 2. and 3., the data column in `inputData` is associated to this object, so the next call can return it immediately.
RooSpan<const double> RooRealVar::getValues(RooBatchCompute::RunContext& inputData, const RooArgSet*) const {
  auto item = inputData.spans.find(this);
  if (item != inputData.spans.end()) {
    return item->second;
  }

  for (const auto& var_span : inputData.spans) {
    auto var = var_span.first;
    if (var->namePtr() == namePtr()) {
      // A variable with the same name exists in the input data. Use their values as ours.
      inputData.spans[this] = var_span.second;
      return var_span.second;
    }
  }

  auto output = inputData.makeBatch(this, 1);
  output[0] = _value;

  return output;
}


////////////////////////////////////////////////////////////////////////////////
/// Set value of variable to 'value'. If 'value' is outside
/// range of object, clip value into range

#include <unistd.h> // getpid, pid_t


void RooRealVar::setVal(Double_t value)
{
  Double_t clipValue ;
  inRange(value,0,&clipValue) ;

//  printf("setVal: %s %p on PID %d\n", GetName(), this, getpid());

  if (clipValue != _value) {
    setValueDirty() ;
    _value = clipValue;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Set value of variable to `value`. If `value` is outside of the
/// range named `rangeName`, clip value into that range.
void RooRealVar::setVal(Double_t value, const char* rangeName)
{
  Double_t clipValue ;
  inRange(value,rangeName,&clipValue) ;

  if (clipValue != _value) {
    setValueDirty() ;
    _value = clipValue;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Return a RooAbsRealLValue representing the error associated
/// with this variable. The callers takes ownership of the
/// return object

RooErrorVar* RooRealVar::errorVar() const
{
  TString name(GetName()), title(GetTitle()) ;
  name.Append("err") ;
  title.Append(" Error") ;

  return new RooErrorVar(name,title,*this) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Returns true if variable has a binning named 'name'.

Bool_t RooRealVar::hasBinning(const char* name) const
{
  return sharedProp()->_altBinning.find(name) != sharedProp()->_altBinning.end();
}



////////////////////////////////////////////////////////////////////////////////
/// Return binning definition with name. If binning with 'name' is not found it is created
/// on the fly as a clone of the default binning if createOnTheFly is true, otherwise
/// a reference to the default binning is returned. If verbose is true a message
/// is printed if a binning is created on the fly.

const RooAbsBinning& RooRealVar::getBinning(const char* name, Bool_t verbose, Bool_t createOnTheFly) const
{
  return const_cast<RooRealVar*>(this)->getBinning(name, verbose, createOnTheFly) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return binning definition with name. If binning with 'name' is not found it is created
/// on the fly as a clone of the default binning if createOnTheFly is true, otherwise
/// a reference to the default binning is returned. If verbose is true a message
/// is printed if a binning is created on the fly.

RooAbsBinning& RooRealVar::getBinning(const char* name, Bool_t verbose, Bool_t createOnTheFly)
{
  // Return default (normalization) binning and range if no name is specified
  if (name==0) {
    return *_binning ;
  }

  if (strchr(name, ',')) {
    coutW(InputArguments) << "Asking variable " << GetName() << "for binning '" << name
        << "', but comma in binning names is not supported." << std::endl;
  }

  // Check if non-shared binning with this name has been created already
  auto item = _altNonSharedBinning.find(name);
  if (item != _altNonSharedBinning.end()) {
    return *item->second;
  }

  // Check if binning with this name has been created already
  auto item2 = sharedProp()->_altBinning.find(name);
  if (item2 != sharedProp()->_altBinning.end()) {
    return *item2->second;
  }


  // Return default binning if requested binning doesn't exist
  if (!createOnTheFly) {
    return *_binning ;
  }

  // Create a new RooRangeBinning with this name with default range
  auto binning = new RooRangeBinning(getMin(),getMax(),name) ;
  if (verbose) {
    coutI(Eval) << "RooRealVar::getBinning(" << GetName() << ") new range named '"
		<< name << "' created with default bounds" << endl ;
  }
  sharedProp()->_altBinning[name] = binning;

  return *binning ;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a list of all binning names. An empty name implies the default binning and
/// a NULL pointer should be passed to getBinning in this case.

std::list<std::string> RooRealVar::getBinningNames() const
{
  std::list<std::string> binningNames;
  if (_binning) {
    binningNames.push_back("");
  }

  for (const auto& item : _altNonSharedBinning) {
    binningNames.push_back(item.first);
  }
  for (const auto& item : sharedProp()->_altBinning) {
    binningNames.push_back(item.first);
  }

  return binningNames;
}

void RooRealVar::removeMin(const char* name) {
  getBinning(name).setMin(-RooNumber::infinity());
}
void RooRealVar::removeMax(const char* name) {
  getBinning(name).setMax(RooNumber::infinity());
}
void RooRealVar::removeRange(const char* name) {
  getBinning(name).setRange(-RooNumber::infinity(),RooNumber::infinity());
}


////////////////////////////////////////////////////////////////////////////////
/// Create a uniform binning under name 'name' for this variable.
/// \param[in] nBins Number of bins. The limits are taken from the currently set limits.
/// \param[in] name Optional name. If name is null, install as default binning.
void RooRealVar::setBins(Int_t nBins, const char* name) {
  setBinning(RooUniformBinning(getMin(name),getMax(name),nBins),name);
}

////////////////////////////////////////////////////////////////////////////////
/// Add given binning under name 'name' with this variable. If name is null,
/// the binning is installed as the default binning.
void RooRealVar::setBinning(const RooAbsBinning& binning, const char* name)
{
  std::unique_ptr<RooAbsBinning> newBinning( binning.clone() );

  // Process insert hooks required for parameterized binnings
  if (!name || name[0] == 0) {
    if (_binning) {
      _binning->removeHook(*this) ;
    }
    newBinning->insertHook(*this) ;
    _binning = std::move(newBinning);
  } else {
    // Remove any old binning with this name
    auto sharedProps = sharedProp();
    auto item = sharedProps->_altBinning.find(name);
    if (item != sharedProps->_altBinning.end()) {
      item->second->removeHook(*this);
      if (sharedProps->_ownBinnings)
        delete item->second;

      sharedProps->_altBinning.erase(item);
    }
    auto item2 = _altNonSharedBinning.find(name);
    if (item2 != _altNonSharedBinning.end()) {
      item2->second->removeHook(*this);
      _altNonSharedBinning.erase(item2);
    }

    // Install new
    newBinning->SetName(name) ;
    newBinning->SetTitle(name) ;
    newBinning->insertHook(*this) ;
    if (newBinning->isShareable()) {
      sharedProp()->_altBinning[name] = newBinning.release();
    } else {
      _altNonSharedBinning[name] = std::move(newBinning);
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Set minimum of name range to given value. If name is null
/// minimum of default range is set

void RooRealVar::setMin(const char* name, Double_t value)
{
  // Set new minimum of fit range
  RooAbsBinning& binning = getBinning(name,kTRUE,kTRUE) ;

  // Check if new limit is consistent
  if (value >= getMax()) {
    coutW(InputArguments) << "RooRealVar::setMin(" << GetName()
			  << "): Proposed new fit min. larger than max., setting min. to max." << endl ;
    binning.setMin(getMax()) ;
  } else {
    binning.setMin(value) ;
  }

  // Clip current value in window if it fell out
  if (!name) {
    Double_t clipValue ;
    if (!inRange(_value,0,&clipValue)) {
      setVal(clipValue) ;
    }
  }

  setShapeDirty() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Set maximum of name range to given value. If name is null
/// maximum of default range is set

void RooRealVar::setMax(const char* name, Double_t value)
{
  // Set new maximum of fit range
  RooAbsBinning& binning = getBinning(name,kTRUE,kTRUE) ;

  // Check if new limit is consistent
  if (value < getMin()) {
    coutW(InputArguments) << "RooRealVar::setMax(" << GetName()
			  << "): Proposed new fit max. smaller than min., setting max. to min." << endl ;
    binning.setMax(getMin()) ;
  } else {
    binning.setMax(value) ;
  }

  // Clip current value in window if it fell out
  if (!name) {
    Double_t clipValue ;
    if (!inRange(_value,0,&clipValue)) {
      setVal(clipValue) ;
    }
  }

  setShapeDirty() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Set a fit or plotting range.
/// Ranges can be selected for e.g. fitting, plotting or integration. Note that multiple
/// variables can have ranges with the same name, so multi-dimensional PDFs can be sliced.
/// See also the tutorial rf203_ranges.C
/// \param[in] name Name this range (so it can be selected later for fitting or
/// plotting). If the name is `nullptr`, the function sets the limits of the default range.
/// \param[in] min Miniminum of the range.
/// \param[in] max Maximum of the range.
void RooRealVar::setRange(const char* name, Double_t min, Double_t max)
{
  Bool_t exists = name == nullptr || sharedProp()->_altBinning.count(name) > 0;

  // Set new fit range
  RooAbsBinning& binning = getBinning(name,kFALSE,kTRUE) ;

  // Check if new limit is consistent
  if (min>max) {
    coutW(InputArguments) << "RooRealVar::setRange(" << GetName()
			  << "): Proposed new fit max. smaller than min., setting max. to min." << endl ;
    binning.setRange(min,min) ;
  } else {
    binning.setRange(min,max) ;
  }

  if (!exists) {
    coutI(Eval) << "RooRealVar::setRange(" << GetName()
		<< ") new range named '" << name << "' created with bounds ["
		<< min << "," << max << "]" << endl ;
  }

  setShapeDirty() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Set or modify a parameterised range, i.e., a range the varies in dependence
/// of parameters.
/// See setRange() for more details.
void RooRealVar::setRange(const char* name, RooAbsReal& min, RooAbsReal& max)
{
  RooParamBinning pb(min,max,100) ;
  setBinning(pb,name) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Read object contents from given stream

Bool_t RooRealVar::readFromStream(istream& is, Bool_t compact, Bool_t verbose)
{
  TString token,errorPrefix("RooRealVar::readFromStream(") ;
  errorPrefix.Append(GetName()) ;
  errorPrefix.Append(")") ;
  RooStreamParser parser(is,errorPrefix) ;
  Double_t value(0) ;

  if (compact) {
    // Compact mode: Read single token
    if (parser.readDouble(value,verbose)) return kTRUE ;
    if (isValidReal(value,verbose)) {
      setVal(value) ;
      return kFALSE ;
    } else {
      return kTRUE ;
    }

  } else {
    // Extended mode: Read multiple tokens on a single line
    Bool_t haveValue(kFALSE) ;
    Bool_t haveConstant(kFALSE) ;
    removeError() ;
    removeAsymError() ;

    Bool_t reprocessToken = kFALSE ;
    while(1) {
      if (parser.atEOL() || parser.atEOF()) break ;

      if (!reprocessToken) {
	token=parser.readToken() ;
      }
      reprocessToken = kFALSE ;

      if (!token.CompareTo("+")) {

	// Expect +/- as 3-token sequence
	if (parser.expectToken("/",kTRUE) ||
	    parser.expectToken("-",kTRUE)) {
	  break ;
	}

	// Next token is error or asymmetric error, check if first char of token is a '('
	TString tmp = parser.readToken() ;
	if (tmp.CompareTo("(")) {
	  // Symmetric error, convert token do double

	  Double_t error ;
	  parser.convertToDouble(tmp,error) ;
	  setError(error) ;

	} else {
	  // Have error
	  Double_t asymErrLo=0., asymErrHi=0.;
	  if (parser.readDouble(asymErrLo,kTRUE) ||
	      parser.expectToken(",",kTRUE) ||
	      parser.readDouble(asymErrHi,kTRUE) ||
	      parser.expectToken(")",kTRUE)) break ;
	  setAsymError(asymErrLo,asymErrHi) ;
	}

      } else if (!token.CompareTo("C")) {

	// Set constant
	setConstant(kTRUE) ;
	haveConstant = kTRUE ;

      } else if (!token.CompareTo("P")) {

	// Next tokens are plot limits
	Double_t plotMin(0), plotMax(0) ;
        Int_t plotBins(0) ;
	if (parser.expectToken("(",kTRUE) ||
	    parser.readDouble(plotMin,kTRUE) ||
	    parser.expectToken("-",kTRUE) ||
	    parser.readDouble(plotMax,kTRUE) ||
            parser.expectToken(":",kTRUE) ||
            parser.readInteger(plotBins,kTRUE) ||
	    parser.expectToken(")",kTRUE)) break ;
//   	setPlotRange(plotMin,plotMax) ;
	coutW(Eval) << "RooRealVar::readFromStrem(" << GetName()
	     << ") WARNING: plot range deprecated, removed P(...) token" << endl ;

      } else if (!token.CompareTo("F")) {

	// Next tokens are fit limits
	Double_t fitMin, fitMax ;
	Int_t fitBins ;
	if (parser.expectToken("(",kTRUE) ||
	    parser.readDouble(fitMin,kTRUE) ||
	    parser.expectToken("-",kTRUE) ||
	    parser.readDouble(fitMax,kTRUE) ||
	    parser.expectToken(":",kTRUE) ||
	    parser.readInteger(fitBins,kTRUE) ||
	    parser.expectToken(")",kTRUE)) break ;
	//setBins(fitBins) ;
	//setRange(fitMin,fitMax) ;
	coutW(Eval) << "RooRealVar::readFromStream(" << GetName()
	     << ") WARNING: F(lo-hi:bins) token deprecated, use L(lo-hi) B(bins)" << endl ;
	if (!haveConstant) setConstant(kFALSE) ;

      } else if (!token.CompareTo("L")) {

	// Next tokens are fit limits
	Double_t fitMin = 0.0, fitMax = 0.0;
//	Int_t fitBins ;
	if (parser.expectToken("(",kTRUE) ||
	    parser.readDouble(fitMin,kTRUE) ||
	    parser.expectToken("-",kTRUE) ||
	    parser.readDouble(fitMax,kTRUE) ||
	    parser.expectToken(")",kTRUE)) break ;
	setRange(fitMin,fitMax) ;
	if (!haveConstant) setConstant(kFALSE) ;

      } else if (!token.CompareTo("B")) {

	// Next tokens are fit limits
	Int_t fitBins = 0;
	if (parser.expectToken("(",kTRUE) ||
	    parser.readInteger(fitBins,kTRUE) ||
	    parser.expectToken(")",kTRUE)) break ;
	setBins(fitBins) ;

      } else {
	// Token is value
	if (parser.convertToDouble(token,value)) { parser.zapToEnd() ; break ; }
	haveValue = kTRUE ;
	// Defer value assignment to end
      }
    }
    if (haveValue) setVal(value) ;
    return kFALSE ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Write object contents to given stream

void RooRealVar::writeToStream(ostream& os, Bool_t compact) const
{
  if (compact) {
    // Write value only
    os << getVal() ;
  } else {

    // Write value with error (if not zero)
    if (_printScientific) {
      char fmtVal[16], fmtErr[16] ;
      snprintf(fmtVal,16,"%%.%de",_printSigDigits) ;
      snprintf(fmtErr,16,"%%.%de",(_printSigDigits+1)/2) ;
      if (_value>=0) os << " " ;
      os << Form(fmtVal,_value) ;

      if (hasAsymError()) {
	os << " +/- (" << Form(fmtErr,getAsymErrorLo())
	   << ", " << Form(fmtErr,getAsymErrorHi()) << ")" ;
      } else  if (hasError()) {
	os << " +/- " << Form(fmtErr,getError()) ;
      }

      os << " " ;
    } else {
      TString* tmp = format(_printSigDigits,"EFA") ;
      os << tmp->Data() << " " ;
      delete tmp ;
    }

    // Append limits if not constants
    if (isConstant()) {
      os << "C " ;
    }

    // Append fit limits
    os << "L(" ;
    if(hasMin()) {
      os << getMin();
    }
    else {
      os << "-INF";
    }
    if(hasMax()) {
      os << " - " << getMax() ;
    }
    else {
      os << " - +INF";
    }
    os << ") " ;

    if (getBins()!=100) {
      os << "B(" << getBins() << ") " ;
    }

    // Add comment with unit, if unit exists
    if (!_unit.IsNull())
      os << "// [" << getUnit() << "]" ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Print value of variable

void RooRealVar::printValue(ostream& os) const
{
  os << getVal() ;

  if(hasError() && !hasAsymError()) {
    os << " +/- " << getError() ;
  } else if (hasAsymError()) {
    os << " +/- (" << getAsymErrorLo() << "," << getAsymErrorHi() << ")" ;
  }

}


////////////////////////////////////////////////////////////////////////////////
/// Print extras of variable: (asymmetric) error, constant flag, limits and binning

void RooRealVar::printExtras(ostream& os) const
{
  // Append limits if not constants
  if (isConstant()) {
    os << "C " ;
  }

  // Append fit limits
  os << " L(" ;
  if(hasMin()) {
    os << getMin();
  }
  else {
    os << "-INF";
  }
  if(hasMax()) {
    os << " - " << getMax() ;
  }
  else {
    os << " - +INF";
  }
  os << ") " ;

  if (getBins()!=100) {
    os << "B(" << getBins() << ") " ;
  }

  // Add comment with unit, if unit exists
  if (!_unit.IsNull())
    os << "// [" << getUnit() << "]" ;

//   cout << " _value = " << &_value << " _error = " << &_error ;


}


////////////////////////////////////////////////////////////////////////////////
/// Mapping of Print() option string to RooPrintable contents specifications

Int_t RooRealVar::defaultPrintContents(Option_t* opt) const
{
  if (opt && TString(opt)=="I") {
    return kName|kClassName|kValue ;
  }
  return kName|kClassName|kValue|kExtras ;
}


////////////////////////////////////////////////////////////////////////////////
/// Detailed printing interface

void RooRealVar::printMultiline(ostream& os, Int_t contents, Bool_t verbose, TString indent) const
{
  RooAbsRealLValue::printMultiline(os,contents,verbose,indent);
  os << indent << "--- RooRealVar ---" << endl;
  TString unit(_unit);
  if(!unit.IsNull()) unit.Prepend(' ');
  os << indent << "  Error = " << getError() << unit << endl;
}



////////////////////////////////////////////////////////////////////////////////
/// Format contents of RooRealVar for pretty printing on RooPlot
/// parameter boxes. This function processes the named arguments
/// taken by paramOn() and translates them to an option string
/// parsed by RooRealVar::format(Int_t sigDigits, const char *options)

TString* RooRealVar::format(const RooCmdArg& formatArg) const
{
  RooCmdArg tmp(formatArg) ;
  tmp.setProcessRecArgs(kTRUE) ;

  RooCmdConfig pc(Form("RooRealVar::format(%s)",GetName())) ;
  pc.defineString("what","FormatArgs",0,"") ;
  pc.defineInt("autop","FormatArgs::AutoPrecision",0,2) ;
  pc.defineInt("fixedp","FormatArgs::FixedPrecision",0,2) ;
  pc.defineInt("tlatex","FormatArgs::TLatexStyle",0,0) ;
  pc.defineInt("latex","FormatArgs::LatexStyle",0,0) ;
  pc.defineInt("latext","FormatArgs::LatexTableStyle",0,0) ;
  pc.defineInt("verbn","FormatArgs::VerbatimName",0,0) ;
  pc.defineMutex("FormatArgs::TLatexStyle","FormatArgs::LatexStyle","FormatArgs::LatexTableStyle") ;
  pc.defineMutex("FormatArgs::AutoPrecision","FormatArgs::FixedPrecision") ;

  // Process & check varargs
  pc.process(tmp) ;
  if (!pc.ok(kTRUE)) {
    return 0 ;
  }

  // Extract values from named arguments
  TString options ;
  options = pc.getString("what") ;

  if (pc.getInt("tlatex")) {
    options += "L" ;
  } else if (pc.getInt("latex")) {
    options += "X" ;
  } else if (pc.getInt("latext")) {
    options += "Y" ;
  }

  if (pc.getInt("verbn")) options += "V" ;
  Int_t sigDigits = 2 ;
  if (pc.hasProcessed("FormatArgs::AutoPrecision")) {
    options += "P" ;
    sigDigits = pc.getInt("autop") ;
  } else if (pc.hasProcessed("FormatArgs::FixedPrecision")) {
    options += "F" ;
    sigDigits = pc.getInt("fixedp") ;
  }

  return format(sigDigits,options) ;
}




////////////////////////////////////////////////////////////////////////////////
/// Format numeric value of RooRealVar and its error in a variety of ways
///
/// To control what is shown use the following options
/// N = show name
/// H = hide value
/// E = show error
/// A = show asymmetric error instead of parabolic error (if available)
/// U = show unit
///
/// To control how it is shown use these options
/// L = TLatex mode
/// X = Latex mode
/// Y = Latex table mode ( '=' replaced by '&' )
/// V = Make name \\verbatim in Latex mode
/// P = use error to control shown precision
/// F = force fixed precision
///

TString *RooRealVar::format(Int_t sigDigits, const char *options) const
{
  //cout << "format = " << options << endl ;

  // parse the options string
  TString opts(options);
  opts.ToLower();
  Bool_t showName= opts.Contains("n");
  Bool_t hideValue= opts.Contains("h");
  Bool_t showError= opts.Contains("e");
  Bool_t showUnit= opts.Contains("u");
  Bool_t tlatexMode= opts.Contains("l");
  Bool_t latexMode= opts.Contains("x");
  Bool_t latexTableMode = opts.Contains("y") ;
  Bool_t latexVerbatimName = opts.Contains("v") ;

  if (latexTableMode) latexMode = kTRUE ;
  Bool_t asymError= opts.Contains("a") ;
  Bool_t useErrorForPrecision= (((showError && hasError(kFALSE) && !isConstant()) || opts.Contains("p")) && !opts.Contains("f")) ;
  // calculate the precision to use
  if(sigDigits < 1) sigDigits= 1;
  Int_t leadingDigitVal = 0;
  if (useErrorForPrecision) {
    leadingDigitVal = (Int_t)floor(log10(fabs(_error+1e-10)));
    if (_value==0&&_error==0) leadingDigitVal=0 ;
  } else {
    leadingDigitVal = (Int_t)floor(log10(fabs(_value+1e-10)));
    if (_value==0) leadingDigitVal=0 ;
  }
  Int_t leadingDigitErr= (Int_t)floor(log10(fabs(_error+1e-10)));
  Int_t whereVal= leadingDigitVal - sigDigits + 1;
  Int_t whereErr= leadingDigitErr - sigDigits + 1;
  char fmtVal[16], fmtErr[16];

  if (_value<0) whereVal -= 1 ;
  snprintf(fmtVal,16,"%%.%df", whereVal < 0 ? -whereVal : 0);
  snprintf(fmtErr,16,"%%.%df", whereErr < 0 ? -whereErr : 0);
  TString *text= new TString();
  if(latexMode) text->Append("$");
  // begin the string with "<name> = " if requested
  if(showName) {
    if (latexTableMode && latexVerbatimName) {
      text->Append("\\verb+") ;
    }
    text->Append(getPlotLabel());
    if (latexVerbatimName) text->Append("+") ;

    if (!latexTableMode) {
      text->Append(" = ");
    } else {
      text->Append(" $ & $ ");
    }
  }

  // Add leading space if value is positive
  if (_value>=0) text->Append(" ") ;

  // append our value if requested
  char buffer[256];
  if(!hideValue) {
    chopAt(_value, whereVal);
    snprintf(buffer, 256,fmtVal, _value);
    text->Append(buffer);
  }

  // append our error if requested and this variable is not constant
  if(hasError(kFALSE) && showError && !(asymError && hasAsymError(kFALSE))) {
    if(tlatexMode) {
      text->Append(" #pm ");
    }
    else if(latexMode) {
      text->Append("\\pm ");
    }
    else {
      text->Append(" +/- ");
    }
    snprintf(buffer, 256,fmtErr, getError());
    text->Append(buffer);
  }

  if (asymError && hasAsymError() && showError) {
    if(tlatexMode) {
      text->Append(" #pm ");
      text->Append("_{") ;
      snprintf(buffer, 256,fmtErr, getAsymErrorLo());
      text->Append(buffer);
      text->Append("}^{+") ;
      snprintf(buffer, 256,fmtErr, getAsymErrorHi());
      text->Append(buffer);
      text->Append("}") ;
    }
    else if(latexMode) {
      text->Append("\\pm ");
      text->Append("_{") ;
      snprintf(buffer, 256,fmtErr, getAsymErrorLo());
      text->Append(buffer);
      text->Append("}^{+") ;
      snprintf(buffer, 256,fmtErr, getAsymErrorHi());
      text->Append(buffer);
      text->Append("}") ;
    }
    else {
      text->Append(" +/- ");
      text->Append(" (") ;
      snprintf(buffer, 256, fmtErr, getAsymErrorLo());
      text->Append(buffer);
      text->Append(", ") ;
      snprintf(buffer, 256, fmtErr, getAsymErrorHi());
      text->Append(buffer);
      text->Append(")") ;
    }

  }

  // append our units if requested
  if(!_unit.IsNull() && showUnit) {
    text->Append(' ');
    text->Append(_unit);
  }
  if(latexMode) text->Append("$");
  return text;
}



////////////////////////////////////////////////////////////////////////////////
/// Utility to calculate number of decimals to show
/// based on magnitude of error

Double_t RooRealVar::chopAt(Double_t what, Int_t where) const
{
  Double_t scale= pow(10.0,where);
  Int_t trunc= (Int_t)floor(what/scale + 0.5);
  return (Double_t)trunc*scale;
}



////////////////////////////////////////////////////////////////////////////////
/// Overload RooAbsReal::attachToTree to also attach
/// branches for errors  and/or asymmetric errors
/// attribute StoreError and/or StoreAsymError are set

void RooRealVar::attachToVStore(RooVectorDataStore& vstore)
{
  // Follow usual procedure for value

  if (getAttribute("StoreError") || getAttribute("StoreAsymError") || vstore.isFullReal(this) ) {

    RooVectorDataStore::RealFullVector* rfv = vstore.addRealFull(this) ;
    rfv->setBuffer(this,&_value);

    // Attach/create additional branch for error
    if (getAttribute("StoreError") || vstore.hasError(this) ) {
      rfv->setErrorBuffer(&_error) ;
    }

    // Attach/create additional branches for asymmetric error
    if (getAttribute("StoreAsymError") || vstore.hasAsymError(this)) {
      rfv->setAsymErrorBuffer(&_asymErrLo,&_asymErrHi) ;
    }

  } else {

    RooAbsReal::attachToVStore(vstore) ;

  }
}



////////////////////////////////////////////////////////////////////////////////
/// Overload RooAbsReal::attachToTree to also attach
/// branches for errors  and/or asymmetric errors
/// attribute StoreError and/or StoreAsymError are set

void RooRealVar::attachToTree(TTree& t, Int_t bufSize)
{
  // Follow usual procedure for value
  RooAbsReal::attachToTree(t,bufSize) ;
//   cout << "RooRealVar::attachToTree(" << this << ") name = " << GetName()
//        << " StoreError = " << (getAttribute("StoreError")?"T":"F") << endl ;

  // Attach/create additional branch for error
  if (getAttribute("StoreError")) {
    TString errName(GetName()) ;
    errName.Append("_err") ;
    TBranch* branch = t.GetBranch(errName) ;
    if (branch) {
      t.SetBranchAddress(errName,&_error) ;
    } else {
      TString format2(errName);
      format2.Append("/D");
      t.Branch(errName, &_error, (const Text_t*)format2, bufSize);
    }
  }

  // Attach/create additional branches for asymmetric error
  if (getAttribute("StoreAsymError")) {
    TString loName(GetName()) ;
    loName.Append("_aerr_lo") ;
    TBranch* lobranch = t.GetBranch(loName) ;
    if (lobranch) {
      t.SetBranchAddress(loName,&_asymErrLo) ;
    } else {
      TString format2(loName);
      format2.Append("/D");
      t.Branch(loName, &_asymErrLo, (const Text_t*)format2, bufSize);
    }

    TString hiName(GetName()) ;
    hiName.Append("_aerr_hi") ;
    TBranch* hibranch = t.GetBranch(hiName) ;
    if (hibranch) {
      t.SetBranchAddress(hiName,&_asymErrHi) ;
    } else {
      TString format2(hiName);
      format2.Append("/D");
      t.Branch(hiName, &_asymErrHi, (const Text_t*)format2, bufSize);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Overload RooAbsReal::fillTreeBranch to also
/// fill tree branches with (asymmetric) errors
/// if requested.

void RooRealVar::fillTreeBranch(TTree& t)
{
  // First determine if branch is taken
  TString cleanName(cleanBranchName()) ;
  TBranch* valBranch = t.GetBranch(cleanName) ;
  if (!valBranch) {
    coutE(Eval) << "RooAbsReal::fillTreeBranch(" << GetName() << ") ERROR: not attached to tree" << endl ;
    assert(0) ;
  }
  valBranch->Fill() ;

  if (getAttribute("StoreError")) {
    TString errName(GetName()) ;
    errName.Append("_err") ;
    TBranch* errBranch = t.GetBranch(errName) ;
    if (errBranch) errBranch->Fill() ;
  }

  if (getAttribute("StoreAsymError")) {
    TString loName(GetName()) ;
    loName.Append("_aerr_lo") ;
    TBranch* loBranch = t.GetBranch(loName) ;
    if (loBranch) loBranch->Fill() ;

    TString hiName(GetName()) ;
    hiName.Append("_aerr_hi") ;
    TBranch* hiBranch = t.GetBranch(hiName) ;
    if (hiBranch) hiBranch->Fill() ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Copy the cached value of another RooAbsArg to our cache
/// Warning: This function copies the cached values of source,
///          it is the callers responsibility to make sure the cache is clean

void RooRealVar::copyCache(const RooAbsArg* source, Bool_t valueOnly, Bool_t setValDirty)
{
  // Follow usual procedure for valueklog
  RooAbsReal::copyCache(source,valueOnly,setValDirty) ;

  if (valueOnly) return ;

  // Copy error too, if source has one
  RooRealVar* other = dynamic_cast<RooRealVar*>(const_cast<RooAbsArg*>(source)) ;
  if (other) {
    // Copy additional error value
    _error = other->_error ;
    _asymErrLo = other->_asymErrLo ;
    _asymErrHi = other->_asymErrHi ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class RooRealVar.

void RooRealVar::Streamer(TBuffer &R__b)
{
  UInt_t R__s, R__c;
  if (R__b.IsReading()) {

    Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
    RooAbsRealLValue::Streamer(R__b);
    if (R__v==1) {
      coutI(Eval) << "RooRealVar::Streamer(" << GetName() << ") converting version 1 data format" << endl ;
      Double_t fitMin, fitMax ;
      Int_t fitBins ;
      R__b >> fitMin;
      R__b >> fitMax;
      R__b >> fitBins;
      _binning.reset(new RooUniformBinning(fitMin,fitMax,fitBins));
    }
    R__b >> _error;
    R__b >> _asymErrLo;
    R__b >> _asymErrHi;
    if (R__v>=2) {
      RooAbsBinning* binning;
      R__b >> binning;
      _binning.reset(binning);
    }
    if (R__v==3) {
      // In v3, properties were written as pointers, so read now and install:
      RooRealVarSharedProperties* tmpProp;
      R__b >> tmpProp;
      installSharedProp(std::shared_ptr<RooRealVarSharedProperties>(tmpProp));
    }
    if (R__v>=4) {
      // In >= v4, properties were written directly, but they might be the "_nullProp"
      auto tmpProp = std::make_shared<RooRealVarSharedProperties>();
      tmpProp->Streamer(R__b);
      installSharedProp(std::move(tmpProp));
    }

    R__b.CheckByteCount(R__s, R__c, RooRealVar::IsA());

  } else {

    R__c = R__b.WriteVersion(RooRealVar::IsA(), kTRUE);
    RooAbsRealLValue::Streamer(R__b);
    R__b << _error;
    R__b << _asymErrLo;
    R__b << _asymErrHi;
    R__b << _binning.get();
    if (_sharedProp) {
      _sharedProp->Streamer(R__b) ;
    } else {
      _nullProp().Streamer(R__b) ;
    }
    R__b.SetByteCount(R__c, kTRUE);

  }
}

/// Hand out our shared property, create on the fly and register
/// in shared map if necessary.
std::shared_ptr<RooRealVarSharedProperties> RooRealVar::sharedProp() const {
  if (!_sharedProp) {
    const_cast<RooRealVar*>(this)->installSharedProp(std::make_shared<RooRealVarSharedProperties>());
  }

  return _sharedProp;
}


////////////////////////////////////////////////////////////////////////////////
/// Install the shared property into the member _sharedProp.
/// If a property with same name already exists, discard the incoming one,
/// and share the existing.
/// `nullptr` and properties equal to the RooRealVar::_nullProp will not be installed.
void RooRealVar::installSharedProp(std::shared_ptr<RooRealVarSharedProperties>&& prop) {
  if (prop == nullptr || (*prop == _nullProp())) {
    _sharedProp = nullptr;
    return;
  }


  auto& weakPtr = _sharedPropList()[prop->asString().Data()];
  std::shared_ptr<RooRealVarSharedProperties> existingProp;
  if ( (existingProp = weakPtr.lock()) ) {
    // Property exists, discard incoming
    _sharedProp = std::move(existingProp);
    // Incoming is not allowed to delete the binnings now - they are owned by the other instance
    prop->disownBinnings();
  } else {
    // Doesn't exist. Install, register weak pointer for future sharing
    _sharedProp = std::move(prop);
    weakPtr = _sharedProp;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Stop sharing properties.
void RooRealVar::deleteSharedProperties()
{
  _sharedProp.reset();

  for (auto it = _sharedPropList().begin(); it != _sharedPropList().end();) {
    if (it->second.expired()) {
      it = _sharedPropList().erase(it);
    } else {
      ++it;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
/// If true, contents of RooRealVars will be printed in scientific notation

void RooRealVar::printScientific(Bool_t flag)
{
  _printScientific = flag ;
}


////////////////////////////////////////////////////////////////////////////////
/// Set number of digits to show when printing RooRealVars

void RooRealVar::printSigDigits(Int_t ndig)
{
  _printSigDigits = ndig>1?ndig:1 ;
}
