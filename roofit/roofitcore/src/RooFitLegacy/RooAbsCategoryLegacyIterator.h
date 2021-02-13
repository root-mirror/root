// Author: Stephan Hageboeck, CERN  3 Feb 2020

/*****************************************************************************
 * RooFit
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2019, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOABSCATEGORYLEGACYITERATOR_H_
#define ROOABSCATEGORYLEGACYITERATOR_H_

#include "TIterator.h"
#include "RooAbsCategory.h"
#include <map>
#include <vector>
#include <string>

#include "RooFitLegacy/RooCatTypeLegacy.h"

/// \deprecated Legacy class to iterate through legacy RooAbsCategory states.
/// Use RooAbsCategory::begin(), RooAbsCategory::end() or range-based for loops instead.
/// \ingroup Roofitlegacy
class RooAbsCategoryLegacyIterator : public TIterator {
  public:
    RooAbsCategoryLegacyIterator(const std::map<std::string, RooAbsCategory::value_type>& stateNames) :
          _origStateNames(&stateNames), index(-1) {
      populate();
    }

    virtual const TCollection *GetCollection() const override { return nullptr; }

    TObject* Next() override {
      ++index;
      return this->operator*();
    }

    void Reset() override { populate(); index = -1; }

    TObject* operator*() const override {
      if (! (0 <= index && index < (int)_origStateNames->size()) )
        return nullptr;

      return const_cast<RooCatType*>(&_legacyStates[index]);
    }

    RooAbsCategoryLegacyIterator& operator=(const RooAbsCategoryLegacyIterator&) = default;

    TIterator& operator=(const TIterator&) override {
      throw std::logic_error("Assigning from another iterator is not supported for the RooAbsCategoryLegacyIterator.");
    }

  private:
    void populate() {
      _legacyStates.clear();

      for (const auto& item : *_origStateNames) {
        _legacyStates.emplace_back(item.first.c_str(), item.second);
      }
      std::sort(_legacyStates.begin(), _legacyStates.end(), [](const RooCatType& left, const RooCatType& right){
        return left.getVal() < right.getVal();
      });
    }

    const std::map<std::string, RooAbsCategory::value_type>* _origStateNames;
    std::vector<RooCatType> _legacyStates;
    int index;
};


#endif /* ROOABSCATEGORYLEGACYITERATOR_H_ */
