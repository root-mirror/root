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
\file RooAbsCategoryLValue.cxx
\class RooAbsCategoryLValue
\ingroup Roofitcore

RooAbsCategoryLValue is the common abstract base class for objects that represent a
discrete value that can be set from the outside, *i.e.* that may appear on the left
hand side of an assignment ("*lvalue*").

Each implementation must provide the functions setIndex()/setLabel() to allow direct modification
of the value. RooAbsCategoryLValue may be derived, but its functional relation
to other RooAbsArgs must be invertible.
**/

#include "RooAbsCategoryLValue.h"

#include "RooFit.h"
#include "RooArgSet.h"
#include "RooRandom.h"
#include "RooMsgService.h"

#include "TString.h"

using namespace std;

ClassImp(RooAbsCategoryLValue);



////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooAbsCategoryLValue::RooAbsCategoryLValue(const char *name, const char *title) : 
  RooAbsCategory(name,title)
{
  setValueDirty() ;  
  setShapeDirty() ;  
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooAbsCategoryLValue::RooAbsCategoryLValue(const RooAbsCategoryLValue& other, const char* name) :
  RooAbsCategory(other, name), RooAbsLValue(other)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooAbsCategoryLValue::~RooAbsCategoryLValue()
{
}



////////////////////////////////////////////////////////////////////////////////
/// Assignment operator from integer index number

RooAbsArg& RooAbsCategoryLValue::operator=(Int_t index) 
{
  setIndex(index,kTRUE) ;
  return *this ;
}



////////////////////////////////////////////////////////////////////////////////
/// Assignment operator from string pointer

RooAbsArg& RooAbsCategoryLValue::operator=(const char *label) 
{
  setLabel(label) ;
  return *this ;
}



////////////////////////////////////////////////////////////////////////////////
/// Assignment from another RooAbsCategory. This will use the *state name*
/// of the other object to set the corresponding state. This is less efficient
/// then directly assigning the state index.
RooAbsArg& RooAbsCategoryLValue::operator=(const RooAbsCategory& other) 
{
  if (&other==this) return *this ;

  const auto index = lookupIndex(other.getCurrentLabel());
  if (index == std::numeric_limits<value_type>::min()) {
    coutE(ObjectHandling) << "Trying to assign the label '" << other.getCurrentLabel() << "' to category'"
        << GetName() << "', but such a label is not defined." << std::endl;
    return *this;
  }

  _currentIndex = index;
  setValueDirty();

  return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// Set our state to our `n`th defined type.
/// \return true in case of an error.
Bool_t RooAbsCategoryLValue::setOrdinal(UInt_t n)
{
  return setIndex(getOrdinal(n).second, true);
}



////////////////////////////////////////////////////////////////////////////////
/// Copy the cached value from given source and raise dirty flag.
/// It is the callers responsability to ensure that the sources
/// cache is clean(valid) before this function is called, e.g. by
/// calling syncCache() on the source.

void RooAbsCategoryLValue::copyCache(const RooAbsArg* source, Bool_t valueOnly, Bool_t setValDirty) 
{
  RooAbsCategory::copyCache(source,valueOnly,setValDirty) ;

  if (isValid()) {
    setIndex(_currentIndex); // force back-propagation
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Randomize current value.
/// If the result is not in the range, the randomisation is repeated.
void RooAbsCategoryLValue::randomize(const char* rangeName) 
{
  const auto& theStateNames = stateNames();

  if (_insertionOrder.size() == theStateNames.size()) {
    // If users didn't manipulate the state map directly, the order of insertion has to be respected to
    // ensure backward compatibility.
    // This heavily uses strings, though.
    do {
      const UInt_t ordinal = RooRandom::integer(theStateNames.size());
      const auto item = theStateNames.find(_insertionOrder[ordinal]);
      setIndex(item->second);
    } while (!inRange(rangeName));
  } else {
    // When not having to respect the insertion order, can just advance the iterator
    do {
      const UInt_t ordinal = RooRandom::integer(theStateNames.size());
      const auto it = std::next(theStateNames.begin(), ordinal);
      setIndex(it->second);
    } while (!inRange(rangeName));
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Set category to i-th fit bin, which is the i-th registered state.

void RooAbsCategoryLValue::setBin(Int_t ibin, const char* rangeName) 
{
  // Check validity of ibin
  if (ibin<0 || ibin>=numBins(rangeName)) {
    coutE(InputArguments) << "RooAbsCategoryLValue::setBin(" << GetName() << ") ERROR: bin index " << ibin
			  << " is out of range (0," << numBins(rangeName)-1 << ")" << endl ;
    return ;
  }

  if (rangeName) {
    coutF(InputArguments) << "RooAbsCategoryLValue::setBin(" << GetName() << ") ERROR: ranges not implemented"
        " for setting bins in categories." << std::endl;
    throw std::logic_error("Ranges not implemented for setting bins in categories.");
  }

  // Retrieve state corresponding to bin
  const auto& type = getOrdinal(ibin);
  assert(type.second != std::numeric_limits<value_type>::min());

  // Set value to requested state
  setIndex(type.second);
}


////////////////////////////////////////////////////////////////////////////////
/// Return the number of fit bins ( = number of types )
Int_t RooAbsCategoryLValue::numBins(const char* rangeName) const
{
  return numTypes(rangeName) ;
}
