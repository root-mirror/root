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
\class RooCategory
\ingroup Roofitcore

RooCategory is an object to represent discrete states.
States have names and index numbers, and the index numbers can be written into datasets and
used in calculations.
A category is "fundamental", i.e., its value doesn't depend on the value of other objects.
(Objects in datasets cannot depend on other objects' values, they need to be self-consistent.)

A category object can be used to *e.g.* conduct a simultaneous fit of
the same observable in multiple categories.

### Setting up a category
1. A category can be set up like this:
~~~{.cpp}
RooCategory myCat("myCat", "Lepton multiplicity category", {
                  {"0Lep", 0},
                  {"1Lep", 1},
                  {"2Lep", 2},
                  {"3Lep", 3}
});
~~~
2. Like this:
~~~{.cpp}
RooCategory myCat("myCat", "Asymmetry");
myCat["left"]  = -1;
myCat["right"] =  1;
~~~
3. Or like this:
~~~{.cpp}
RooCategory myCat("myCat", "Asymmetry");
myCat.defineType("left", -1);
myCat.defineType("right", 1);
~~~
Inspect the pairs of state names and state numbers like this:
~~~{.cpp}
for (const auto& nameIdx : myCat) {
  std::cout << nameIdx.first << " --> " << nameIdx.second << std::endl;
}
~~~

### Changing category states
Category states can be modified either by using the index state (faster) or state names.
For example:
~~~{.cpp}
myCat.setIndex(5);
myCat.setLabel("left");
for (const auto& otherNameIdx : otherCat) {
  myCat.setIndex(otherNameIdx);
}
~~~

Also refer to \ref tutorial_roofit, especially rf404_categories.C for an introduction, and to rf405_realtocatfuncs.C and rf406_cattocatfuncs.C
for advanced uses of categories.
**/

#include "RooCategory.h"

#include "RooFit.h"
#include "RooArgSet.h"
#include "RooStreamParser.h"
#include "RooMsgService.h"
#include "RooTrace.h"
#include "RooHelpers.h"
#include "RooFitLegacy/RooCategorySharedProperties.h"
#include "RooFitLegacy/RooCatTypeLegacy.h"

#include "TBuffer.h"
#include "TString.h"
#include "ROOT/RMakeUnique.hxx"
#include "TList.h"

#include <iostream>
#include <cstdlib>

using namespace std;

ClassImp(RooCategory);

std::map<std::string, std::weak_ptr<RooCategory::RangeMap_t>> RooCategory::_uuidToSharedRangeIOHelper; // Helper for restoring shared properties
std::map<std::string, std::weak_ptr<RooCategory::RangeMap_t>> RooCategory::_sharedRangeIOHelper;


////////////////////////////////////////////////////////////////////////////////

RooCategory::RooCategory()
{
  TRACE_CREATE 
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor. Types must be defined using defineType() before variable can be used
RooCategory::RooCategory(const char *name, const char *title) : 
  RooAbsCategoryLValue(name,title),
  _ranges(new RangeMap_t())
{
  setValueDirty() ;  
  setShapeDirty() ;  
  TRACE_CREATE 
}


////////////////////////////////////////////////////////////////////////////////
/// Create a new category and define allowed states.
/// \param[in] name Name used to refer to this object.
/// \param[in] title Title for e.g. plotting.
/// \param[in] allowedStates Map of allowed states. Pass e.g. `{ {"0Lep", 0}, {"1Lep:, 1} }`
RooCategory::RooCategory(const char* name, const char* title, const std::map<std::string, int>& allowedStates) :
  RooAbsCategoryLValue(name,title),
  _ranges(new RangeMap_t())
{
  defineTypes(allowedStates);
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooCategory::RooCategory(const RooCategory& other, const char* name) :
  RooAbsCategoryLValue(other, name),
  _ranges(other._ranges)
{
  TRACE_CREATE   
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooCategory::~RooCategory()
{
  TRACE_DESTROY
}




////////////////////////////////////////////////////////////////////////////////
/// Set value by specifying the index code of the desired state.
/// If printError is set, a message will be printed if
/// the specified index does not represent a valid state.
/// \return bool signalling if an error occurred.
Bool_t RooCategory::setIndex(Int_t index, Bool_t printError) 
{
  if (!hasIndex(index)) {
    if (printError) {
      coutE(InputArguments) << "RooCategory: Trying to set invalid state " << index << " for category " << GetName() << std::endl;
    }
    return true;
  }

  _currentIndex = index;
  setValueDirty();

  return false;
}



////////////////////////////////////////////////////////////////////////////////
/// Set value by specifying the name of the desired state.
/// If printError is set, a message will be printed if
/// the specified label does not represent a valid state.
/// \return false on success.
Bool_t RooCategory::setLabel(const char* label, Bool_t printError) 
{
  const auto item = stateNames().find(label);
  if (item != stateNames().end()) {
    _currentIndex = item->second;
    setValueDirty();
    return false;
  }

  if (printError) {
    coutE(InputArguments) << "Trying to set invalid state label '" << label << "' for category " << GetName() << std::endl;
  }

  return true;
}



////////////////////////////////////////////////////////////////////////////////
/// Define a state with given name.
/// The lowest available positive integer is assigned as index. Category
/// state labels may not contain semicolons.
/// \return True in case of an error.
bool RooCategory::defineType(const std::string& label)
{ 
  if (label.find(';') != std::string::npos) {
    coutE(InputArguments) << "RooCategory::defineType(" << GetName()
        << "): semicolons not allowed in label name" << endl ;
    return true;
  }

  return RooAbsCategory::defineState(label) == invalidCategory();
}


////////////////////////////////////////////////////////////////////////////////
/// Define a state with given name and index. Category
/// state labels may not contain semicolons.
/// \return True in case of error.
bool RooCategory::defineType(const std::string& label, Int_t index)
{
  if (label.find(';') != std::string::npos) {
    coutE(InputArguments) << "RooCategory::defineType(" << GetName()
			<< "): semicolons not allowed in label name" << endl ;
    return true;
  }

  return RooAbsCategory::defineState(label, index) == invalidCategory();
}


////////////////////////////////////////////////////////////////////////////////
/// Define multiple states in a single call. Use like:
/// ```
/// myCat.defineTypes({ {"0Lep", 0}, {"1Lep", 1}, {"2Lep", 2}, {"3Lep", 3} });
/// ```
/// Note: When labels or indices are defined multiple times, an error message is printed,
/// and the corresponding state is ignored.
void RooCategory::defineTypes(const std::map<std::string, int>& allowedStates) {
  for (const auto& nameAndIdx : allowedStates) {
    defineType(nameAndIdx.first, nameAndIdx.second);
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Access a named state. If a state with this name doesn't exist yet, the state is
/// assigned the next available positive integer.
/// \param[in] stateName Name of the state to be accessed.
/// \return Reference to the category index. If no state exists, it will be created on the fly.
RooAbsCategory::value_type& RooCategory::operator[](const std::string& stateName) {
  setShapeDirty();
  if (stateNames().count(stateName) == 0) {
    _insertionOrder.push_back(stateName);
    return stateNames()[stateName] = nextAvailableStateIndex();

  }

  return stateNames()[stateName];
}


////////////////////////////////////////////////////////////////////////////////
/// Return a reference to the map of state names to index states.
/// This can be used to manipulate the category.
/// \note Calling this function will **always** trigger recomputations of
/// of **everything** that depends on this category, since in case the map gets
/// manipulated, names or indices might change. Also, the order that states have
/// been inserted in gets lost. This changes what is returned by getOrdinal().
std::map<std::string, RooAbsCategory::value_type>& RooCategory::states() {
  auto& theStates = stateNames();
  setValueDirty();
  setShapeDirty();
  _insertionOrder.clear();
  return theStates;
}


////////////////////////////////////////////////////////////////////////////////
/// Read object contents from given stream. If token is a decimal digit, try to
/// find a corresponding state for it. If that succeeds, the state denoted by this
/// index is used. Otherwise, interpret it as a label.
Bool_t RooCategory::readFromStream(istream& is, Bool_t /*compact*/, Bool_t verbose) 
{
  // Read single token
  RooStreamParser parser(is) ;
  TString token = parser.readToken() ;

  if (token.IsDec() && hasIndex(std::stoi(token.Data()))) {
    return setIndex(std::stoi(token.Data()), verbose);
  } else {
    return setLabel(token,verbose) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// compact only at the moment

void RooCategory::writeToStream(ostream& os, Bool_t compact) const
{
  if (compact) {
    os << getCurrentIndex() ;
  } else {
    os << getCurrentLabel() ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Clear the named range.
/// \note This affects **all** copies of this category, because they are sharing
/// range definitions. This ensures that categories inside a dataset and their
/// counterparts on the outside will both see a modification of the range.
void RooCategory::clearRange(const char* name, Bool_t silent)
{
  std::map<std::string, std::vector<value_type>>::iterator item = _ranges->find(name);
  if (item == _ranges->end()) {
    if (!silent)
      coutE(InputArguments) << "RooCategory::clearRange(" << GetName() << ") ERROR: must specify valid range name" << endl ;
    return;
  }
  
  _ranges->erase(item);
}


////////////////////////////////////////////////////////////////////////////////

void RooCategory::setRange(const char* name, const char* stateNameList) 
{
  clearRange(name,kTRUE) ;
  addToRange(name,stateNameList) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Add the given state to the given range.
/// \note This creates or accesses a **shared** map with allowed ranges. All copies of this
/// category will share this range such that a category inside a dataset and its
/// counterpart on the outside will both see a modification of the range.
void RooCategory::addToRange(const char* name, RooAbsCategory::value_type stateIndex) {
  auto item = _ranges->find(name);
  if (item == _ranges->end()) {
    if (!name) {
      coutE(Contents) << "RooCategory::addToRange(" << GetName()
          << "): Need valid range name." << std::endl;
      return;
    }

    item = _ranges->emplace(name, std::vector<value_type>()).first;
    coutI(Contents) << "RooCategory::setRange(" << GetName() 
        << ") new range named '" << name << "' created for state " << stateIndex << endl ;
  }

  item->second.push_back(stateIndex);
}


////////////////////////////////////////////////////////////////////////////////
/// Add the list of state names to the given range. State names can be separated
/// with ','.
/// \note This creates or accesses a **shared** map with allowed ranges. All copies of this
/// category will share this range such that a category inside a dataset and its
/// counterpart on the outside will both see a modification of the range.
void RooCategory::addToRange(const char* name, const char* stateNameList)
{
  if (!stateNameList) {
    coutE(InputArguments) << "RooCategory::setRange(" << GetName() << ") ERROR: must specify valid name and state name list" << endl ;
    return;
  }

  // Parse list of state names, verify that each is valid and add them to the list
  for (const auto& token : RooHelpers::tokenise(stateNameList, ",")) {
    const value_type idx = lookupIndex(token);
    if (idx != invalidCategory().second) {
      addToRange(name, idx);
    } else {
      coutW(InputArguments) << "RooCategory::setRange(" << GetName() << ") WARNING: Ignoring invalid state name '" 
			    << token << "' in state name list" << endl ;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Check if the state is in the given range.
/// If no range is specified either as argument or if no range has been defined for this category
/// (*i.e.*, the default range is meant), all category states count as being in range.
bool RooCategory::isStateInRange(const char* rangeName, RooAbsCategory::value_type stateIndex) const {
  if (rangeName == nullptr || _ranges->empty())
    return true;

  const auto item = _ranges->find(rangeName);
  if (item == _ranges->end())
    return false;

  const std::vector<value_type>& vec = item->second;
  return std::find(vec.begin(), vec.end(), stateIndex) != vec.end();
}


////////////////////////////////////////////////////////////////////////////////
/// Check if the state is in the given range.
/// If no range is specified (*i.e.*, the default range), all category states count as being in range.
/// This overload requires a name lookup. Recommend to use the category index with
/// RooCategory::isStateInRange(const char*, RooAbsCategory::value_type) const.
bool RooCategory::isStateInRange(const char* rangeName, const char* stateName) const
{
  // Check that both input arguments are not null pointers
  if (!rangeName) {
    return true;
  }

  if (!stateName) {
    coutE(InputArguments) << "RooCategory::isStateInRange(" << GetName() << ") ERROR: must specify valid state name" << endl ;
    return false;
  }

  return isStateInRange(rangeName, lookupIndex(stateName));
}


////////////////////////////////////////////////////////////////////////////////

void RooCategory::Streamer(TBuffer &R__b)
{
  UInt_t R__s, R__c;
  if (R__b.IsReading()) {
    
    Version_t R__v = R__b.ReadVersion(&R__s, &R__c);

    if (R__v==1) {
      RooAbsCategoryLValue::Streamer(R__b);

      // In v1, properties were a direct pointer:
      RooCategorySharedProperties* props = nullptr;
      R__b >> props;
      installLegacySharedProp(props);
      // props was allocated by I/O system, we cannot delete here in case it gets reused

    } else if (R__v == 2) {
      RooAbsCategoryLValue::Streamer(R__b);

      // In v2, properties were written directly into the class buffer
      auto props = std::make_unique<RooCategorySharedProperties>();
      props->Streamer(R__b);
      installLegacySharedProp(props.get());

    } else {
      // Starting at v3, ranges are shared using a shared pointer, which cannot be read by ROOT's I/O.
      // Instead, ranges are written as a normal pointer, and here we restore the sharing.
      R__b.ReadClassBuffer(RooCategory::Class(), this, R__v, R__s, R__c);
      installSharedRange(std::unique_ptr<RangeMap_t>(_rangesPointerForIO));
      _rangesPointerForIO = nullptr;
    }

    R__b.CheckByteCount(R__s, R__c, RooCategory::IsA());
    
  } else {
    // Since we cannot write shared pointers yet, assign the shared ranges to a normal pointer,
    // write, and restore.
    if (_ranges)
      _rangesPointerForIO = _ranges.get();

    R__b.WriteClassBuffer(RooCategory::Class(), this);
    _rangesPointerForIO = nullptr;
  }
}


/// When reading old versions of the class, we get instances of shared properties.
/// Since these only contain ranges with numbers, just convert to vectors of numbers.
void RooCategory::installLegacySharedProp(const RooCategorySharedProperties* props) {
  if (props == nullptr || (*props == RooCategorySharedProperties("00000000-0000-0000-0000-000000000000")))
    return;

  auto& weakPtr = _uuidToSharedRangeIOHelper[props->asString().Data()];
  if (auto existingObject = weakPtr.lock()) {
    // We know this range, start sharing
    _ranges = std::move(existingObject);
  } else {
    // This range is unknown, make a new object
    _ranges = std::make_shared<std::map<std::string, std::vector<value_type>>>();
    auto& rangesMap = *_ranges;

    // Copy the data:
    std::unique_ptr<TIterator> iter(props->_altRanges.MakeIterator());
    while (TList* olist = (TList*)iter->Next()) {
      std::vector<value_type>& vec = rangesMap[olist->GetName()];


      std::unique_ptr<TIterator> citer(olist->MakeIterator());
      while (RooCatType* ctype = (RooCatType*)citer->Next()) {
        vec.push_back(ctype->getVal());
      }
    }

    // Register the shared_ptr for future sharing
    weakPtr = _ranges;
  }
}


/// In current versions of the class, a map with ranges can be shared between instances.
/// If an instance with the same name alreday uses the same map, the instances will start sharing.
/// Otherwise, this instance will be registered, and future copies being read will share with this
/// one.
void RooCategory::installSharedRange(std::unique_ptr<RangeMap_t>&& rangeMap) {
  if (rangeMap == nullptr)
    return;

  auto checkRangeMapsEqual = [](const RooCategory::RangeMap_t& a, const RooCategory::RangeMap_t& b) {
    if (&a == &b)
      return true;

    if (a.size() != b.size())
      return false;

    auto vecsEqual = [](const std::vector<RooAbsCategory::value_type>& aa, const std::vector<RooAbsCategory::value_type>& bb) {
      return aa.size() == bb.size() && std::equal(aa.begin(), aa.end(), bb.begin());
    };

    for (const auto& itemA : a) {
      const auto itemB = b.find(itemA.first);
      if (itemB == b.end())
        return false;

      if (!vecsEqual(itemA.second, itemB->second))
        return false;
    }

    return true;
  };


  auto& weakPtr = _sharedRangeIOHelper[GetName()];
  auto existingMap = weakPtr.lock();
  if (existingMap && checkRangeMapsEqual(*rangeMap, *existingMap)) {
    // We know this map, use the shared one.
    _ranges = std::move(existingMap);
    if (rangeMap.get() == existingMap.get()) {
      // This happens when ROOT's IO has written the same pointer twice. We cannot delete now.
      (void) rangeMap.release(); // clang-tidy is normally right that this leaks. Here, we need to leave the result unused, though.
    }
  } else {
    // We don't know this map. Register for sharing.
    _ranges = std::move(rangeMap);
    weakPtr = _ranges;
  }
}
