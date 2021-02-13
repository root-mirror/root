/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooCategory.h,v 1.27 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_CATEGORY
#define ROO_CATEGORY

#include "RooAbsCategoryLValue.h"

#include <vector>
#include <map>
#include <string>

class RooCategorySharedProperties;

class RooCategory final : public RooAbsCategoryLValue {
public:
  // Constructor, assignment etc.
  RooCategory() ;
  RooCategory(const char *name, const char *title);
  RooCategory(const char* name, const char* title, const std::map<std::string, int>& allowedStates);
  RooCategory(const RooCategory& other, const char* name=0) ;
  RooCategory& operator=(const RooCategory&) = delete;
  virtual ~RooCategory();
  virtual TObject* clone(const char* newname) const override { return new RooCategory(*this,newname); }

  /// Return current index.
  virtual value_type getCurrentIndex() const override final {
    return RooCategory::evaluate();
  }

  virtual Bool_t setIndex(Int_t index, bool printError = true) override;
  using RooAbsCategoryLValue::setIndex;
  virtual Bool_t setLabel(const char* label, bool printError = true) override;
  using RooAbsCategoryLValue::setLabel;
  
  // I/O streaming interface (machine readable)
  virtual Bool_t readFromStream(std::istream& is, Bool_t compact, Bool_t verbose=kFALSE) override;
  virtual void writeToStream(std::ostream& os, Bool_t compact) const override ;

  bool defineType(const std::string& label);
  bool defineType(const std::string& label, Int_t index);
  void defineTypes(const std::map<std::string, int>& allowedStates);
  value_type& operator[](const std::string& stateName);
  std::map<std::string, RooAbsCategory::value_type>& states();

  /// \cond LEGACY
  bool defineType(const char* label) {
    return defineType(std::string(label));
  }
  bool defineType(const char* label, Int_t index) {
    return defineType(std::string(label), index);
  }
  /// \endcond

  /// Clear all defined category states.
  void clear() {
    clearTypes();
  }

  void clearRange(const char* name, Bool_t silent) ;
  void setRange(const char* rangeName, const char* stateNameList) ;
  void addToRange(const char* rangeName, RooAbsCategory::value_type stateIndex);
  void addToRange(const char* rangeName, const char* stateNameList) ;


  /// \name RooFit interface
  /// @{

  /// Tell whether we can be stored in a dataset. Always true for RooCategory.
  inline virtual Bool_t isFundamental() const override {
    return true;
  }

  /// Does our value or shape depend on any other arg? Always false for RooCategory.
  virtual Bool_t isDerived() const override {
    return false;
  }

  Bool_t isStateInRange(const char* rangeName, RooAbsCategory::value_type stateIndex) const ;
  Bool_t isStateInRange(const char* rangeName, const char* stateName) const ;
  /// Check if the currently defined category state is in the range with the given name.
  /// If no ranges are defined, the state counts as being in range.
  virtual Bool_t inRange(const char* rangeName) const override {
    return isStateInRange(rangeName, RooCategory::evaluate());
  }
  /// Returns true if category has a range with given name defined.
  virtual bool hasRange(const char* rangeName) const override {
    return _ranges->find(rangeName) != _ranges->end();
  }

  ///@}

protected:
  /// \copydoc RooAbsCategory::evaluate() const
  /// Returns the currently set state index. If this is invalid,
  /// returns the first-set index.
  virtual value_type evaluate() const override {
    if (hasIndex(_currentIndex))
      return _currentIndex;

    if (_insertionOrder.empty()) {
      return invalidCategory().second;
    } else {
      auto item = stateNames().find(_insertionOrder.front());
      assert(item != stateNames().end());
      return item->second;
    }
  }

  /// This category's shape does not depend on others, and does not need recomputing.
  void recomputeShape() override { };

private:

  using RangeMap_t = std::map<std::string, std::vector<value_type>>;
  /// Map range names to allowed category states. Note that this must be shared between copies,
  /// so categories in datasets have the same ranges as their counterparts outside of the dataset.
  std::shared_ptr<RangeMap_t> _ranges{new RangeMap_t()}; //!
  RangeMap_t* _rangesPointerForIO{nullptr}; // Pointer to the same object as _ranges, but not shared for I/O.

  void installLegacySharedProp(const RooCategorySharedProperties* sp);
  void installSharedRange(std::unique_ptr<RangeMap_t>&& rangeMap);
  // Helper for restoring shared ranges from old versions of this class read from files. Maps TUUID names to shared ranges.
  static std::map<std::string, std::weak_ptr<RangeMap_t>> _uuidToSharedRangeIOHelper;
  // Helper for restoring shared ranges from current versions of this class read from files. Maps category names to shared ranges.
  static std::map<std::string, std::weak_ptr<RangeMap_t>> _sharedRangeIOHelper;

  ClassDefOverride(RooCategory, 3) // Discrete valued variable type
};

#endif
