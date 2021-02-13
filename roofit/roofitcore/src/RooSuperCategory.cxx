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
\file RooSuperCategory.cxx
\class RooSuperCategory
\ingroup Roofitcore

The RooSuperCategory can join several RooAbsCategoryLValue objects into
a single category. For this, it uses a RooMultiCategory, which takes care
of enumerating all the permutations of possible states.
In addition, the super category derives from RooAbsCategoryLValue, *i.e.*, it allows for
setting its state (opposed to the RooMultiCategory, which just reacts
to the states of its subcategories). This requires that all input categories
are lvalues as well. This is because a modification of the state of the
supercategory will propagate to its input categories.
**/

#include "RooSuperCategory.h"

#include "RooFit.h"
#include "Riostream.h"
#include "RooStreamParser.h"
#include "RooArgSet.h"
#include "RooFitLegacy/RooMultiCatIter.h"
#include "RooAbsCategoryLValue.h"
#include "RooMsgService.h"

#include "TString.h"
#include "TClass.h"

using namespace std;

ClassImp(RooSuperCategory);

RooSuperCategory::RooSuperCategory() :
  RooAbsCategoryLValue(),
  _multiCat("MultiCatProxy", "Stores a RooMultiCategory", this, true, true, true) { }

////////////////////////////////////////////////////////////////////////////////
/// Construct a super category from other categories.
/// \param[in] name Name of this object
/// \param[in] title Title (for e.g. printing)
/// \param[in] inputCategories RooArgSet with category objects. These all need to derive from RooAbsCategoryLValue, *i.e.*
/// one needs to be able to assign to them.
RooSuperCategory::RooSuperCategory(const char *name, const char *title, const RooArgSet& inputCategories) :
  RooAbsCategoryLValue(name, title),
  _multiCat("MultiCatProxy", "Stores a RooMultiCategory", this,
      *new RooMultiCategory((std::string(name) + "_internalMultiCat").c_str(), title, inputCategories), true, true, true)
{  
  // Check category list
  for (const auto arg : inputCategories) {
    if (!arg->IsA()->InheritsFrom(RooAbsCategoryLValue::Class())) {
      coutE(InputArguments) << "RooSuperCategory::RooSuperCategory(" << GetName() << "): input category " << arg->GetName() 
			    << " is not an lvalue. Use RooMultiCategory instead." << endl ;
      throw std::invalid_argument("Arguments of RooSuperCategory must be lvalues.");
    }
  }
  setShapeDirty();
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooSuperCategory::RooSuperCategory(const RooSuperCategory& other, const char *name) :
    RooAbsCategoryLValue(other, name),
    _multiCat("MultiCatProxy", this, other._multiCat)
{
  RooSuperCategory::setIndex(other.getCurrentIndex(), true);
  setShapeDirty();
}


////////////////////////////////////////////////////////////////////////////////
/// Make an iterator over all state permutations of 
/// the input categories of this supercategory.
/// The iterator just generates state names, it does not set them.
TIterator* RooSuperCategory::MakeIterator() const 
{
  return new RooMultiCatIter(_multiCat->inputCatList());
}


////////////////////////////////////////////////////////////////////////////////
/// Set the value of the super category to the specified index.
/// This will propagate to the sub-categories, and set their state accordingly.
bool RooSuperCategory::setIndex(Int_t index, Bool_t printError)
{
  if (index < 0) {
    if (printError)
      coutE(InputArguments) << "RooSuperCategory can only have positive index states. Got " << index << std::endl;
    return true;
  }

  bool error = false;
  for (auto arg : _multiCat->_catSet) {
    auto cat = static_cast<RooAbsCategoryLValue*>(arg);
    if (cat->size() == 0) {
      if (printError)
        coutE(InputArguments) << __func__ << ": Found a category with zero states. Cannot set state for '"
            << cat->GetName() << "'." << std::endl;
      continue;
    }
    const value_type thisIndex = index % cat->size();
    error |= cat->setOrdinal(thisIndex);
    index = (index - thisIndex) / cat->size();
  }

  return error;
}



////////////////////////////////////////////////////////////////////////////////
/// Set the value of the super category by specifying the state name.
/// This looks up the corresponding index number, and calls setIndex().
Bool_t RooSuperCategory::setLabel(const char* label, Bool_t printError)
{
  const value_type index = _multiCat->lookupIndex(label);
  return setIndex(index, printError);
}


////////////////////////////////////////////////////////////////////////////////
/// Print the state of this object to the specified output stream.

void RooSuperCategory::printMultiline(ostream& os, Int_t content, Bool_t verbose, TString indent) const
{
  RooAbsCategory::printMultiline(os,content,verbose,indent) ;
  
  if (verbose) {     
    os << indent << "--- RooSuperCategory ---" << '\n';
    os << indent << "  Internal RooMultiCategory:" << '\n';
    _multiCat->printMultiline(os, content, verbose, indent+"  ");

    os << std::endl;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Check that all input category states are in the given range.
Bool_t RooSuperCategory::inRange(const char* rangeName) const 
{
  for (const auto c : _multiCat->inputCatList()) {
    auto cat = static_cast<RooAbsCategoryLValue*>(c);
    if (!cat->inRange(rangeName)) {
      return false;
    }
  }

  return true;
}


////////////////////////////////////////////////////////////////////////////////
/// Check that any of the input categories has a range with the given name.
Bool_t RooSuperCategory::hasRange(const char* rangeName) const 
{
  for (const auto c : _multiCat->inputCatList()) {
    auto cat = static_cast<RooAbsCategoryLValue*>(c);
    if (cat->hasRange(rangeName)) return true;
  }

  return false;
}
