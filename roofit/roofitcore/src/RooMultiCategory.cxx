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
\file RooMultiCategory.cxx
\class RooMultiCategory
\ingroup Roofitcore

RooMultiCategory connects several RooAbsCategory objects into
a single category. The states of the multi-category consist of all the permutations
of the input categories. 
RooMultiCategory states are automatically defined and updated whenever an input
category modifies its list of states.

A RooMultiCategory is not an lvalue, *i.e.* one cannot set its states. Its state simply follows
as a computation from the states of the input categories. This is because the input categories
don't need to be lvalues, so their states cannot be set by the MultiCategory. If all input categories
are lvalues, the RooSuperCategory can be used. It works like RooMultiCategory, but allows for
setting the states.
**/

#include "RooMultiCategory.h"

#include "RooFit.h"
#include "RooStreamParser.h"
#include "RooArgSet.h"
#include "RooAbsCategory.h"
#include "RooMsgService.h"

#include "TString.h"

using namespace std;

ClassImp(RooMultiCategory);



////////////////////////////////////////////////////////////////////////////////
/// Construct a product of the given set of input RooAbsCategories in `inInputCatList`.
/// The state names of this product category are {S1;S2,S3,...Sn} where Si are the state names
/// of the input categories.

RooMultiCategory::RooMultiCategory(const char *name, const char *title, const RooArgSet& inputCategories) :
  RooAbsCategory(name, title), _catSet("input","Input category set",this,kTRUE,kTRUE)
{  
  // Copy category list
  for (const auto arg : inputCategories) {
    if (!dynamic_cast<RooAbsCategory*>(arg)) {
      coutE(InputArguments) << "RooMultiCategory::RooMultiCategory(" << GetName() << "): input argument " << arg->GetName() 
			    << " is not a RooAbsCategory" << endl ;
    }
    _catSet.add(*arg) ;
  }
  setShapeDirty();
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooMultiCategory::RooMultiCategory(const RooMultiCategory& other, const char *name) :
  RooAbsCategory(other,name), _catSet("input",this,other._catSet)
{
  setShapeDirty();
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooMultiCategory::~RooMultiCategory() 
{
}


////////////////////////////////////////////////////////////////////////////////
/// Compile a string with all the labels of the serving categories,
/// such as `{1Jet;1Lepton;2Tag}`.
std::string RooMultiCategory::createLabel() const
{
  // Construct composite label name
  std::string label;
  Bool_t first = true;
  for (const auto arg : _catSet) {
    auto cat = static_cast<const RooAbsCategory*>(arg);

    label += first ? '{' : ';';
    label += cat->getCurrentLabel();
    first = false;
  }
  label += '}';

  return label ;
}


#ifndef NDEBUG

#include "RooFitLegacy/RooMultiCatIter.h"
namespace {
/// Check that root-6.22 redesign of category interfaces yields same labels
std::string computeLabelOldStyle(const RooArgSet& catSet, unsigned int index) {
  RooMultiCatIter iter(catSet) ;
  TObjString* obj ;
  for (unsigned int i=0; (obj=(TObjString*)iter.Next()); ++i) {
    if (i == index) {
      return obj->String().Data();
    }
  }

  return {};
}
}
#endif


////////////////////////////////////////////////////////////////////////////////
/// Calculate the current value.
/// This enumerates the states of each serving category, and calculates a unique
/// state number. The first category occupies the state numbers \f$ 0, \ldots \mathrm{size}_\mathrm{first}-1 \f$,
/// the second category \f$ (0, \ldots \mathrm{size}_\mathrm{second}-1) * \mathrm{size}_\mathrm{first} \f$ etc.
RooAbsCategory::value_type RooMultiCategory::evaluate() const
{
  value_type computedStateIndex = 0;
  value_type multiplier = 1;
  for (const auto arg : _catSet) {
    auto cat = static_cast<const RooAbsCategory*>(arg);
    if (cat->size() == 0) {
      coutW(InputArguments) << __func__ << " Trying to build a multi-category state based on "
          "a category with zero states. Fix '" << cat->GetName() << "'." << std::endl;
      continue;
    }
    computedStateIndex += cat->getCurrentOrdinalNumber() * multiplier;
    multiplier *= cat->size();
  }

#ifndef NDEBUG
  assert(hasIndex(computedStateIndex));
  _currentIndex = computedStateIndex;
  assert(createLabel() == computeLabelOldStyle(_catSet, computedStateIndex));
#endif

  return computedStateIndex;
}



////////////////////////////////////////////////////////////////////////////////
/// Print the state of this object to the specified output stream.

void RooMultiCategory::printMultiline(ostream& os, Int_t content, Bool_t verbose, TString indent) const
{
  RooAbsCategory::printMultiline(os,content,verbose,indent) ;
  
  if (verbose) {     
    os << indent << "--- RooMultiCategory ---" << endl;
    os << indent << "  Input category list:" << endl ;
    TString moreIndent(indent) ;
    moreIndent.Append("   ") ;
    _catSet.printStream(os,kName|kValue,kStandard,moreIndent.Data()) ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Write object contents to given stream

void RooMultiCategory::writeToStream(ostream& os, Bool_t compact) const
{
  RooAbsCategory::writeToStream(os,compact) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Get current label. If labels haven't been computed, yet, or if the shape is
/// dirty, a recomputation is triggered.
const char* RooMultiCategory::getCurrentLabel() const {
  for (const auto& item : stateNames()) {
    if (item.second == getCurrentIndex())
      return item.first.c_str();
  }

  return "";
}


////////////////////////////////////////////////////////////////////////////////
/// Inspect all the subcategories, and enumerate and name their states.
void RooMultiCategory::recomputeShape() {
  // Propagate up:
  setShapeDirty();

  clearTypes();

  unsigned int totalSize = 1;
  for (const auto arg : _catSet) {
    auto cat = static_cast<const RooAbsCategory*>(arg);
    totalSize *= cat->size();
  }

  for (unsigned int i=0; i < totalSize; ++i) {
    unsigned int workingIndex = i;
    std::string catName = "{";
    for (const auto arg : _catSet) {
      auto cat = static_cast<const RooAbsCategory*>(arg);
      unsigned int thisStateOrdinal = workingIndex % cat->size();
      const auto& thisState = cat->getOrdinal(thisStateOrdinal);
      catName += thisState.first + ';';
      workingIndex = (workingIndex - thisStateOrdinal) / cat->size();
    }
    catName[catName.size()-1] = '}';

    // It's important that we define the states unchecked, because for checking that name
    // or index are available, recomputeShape() would be called.
    defineStateUnchecked(catName, i);
  }
  assert(_stateNames.size() == totalSize);
  assert(_insertionOrder.size() == totalSize);

  // Possible new state numbers will invalidate all cached numbers
  setValueDirty();
}
