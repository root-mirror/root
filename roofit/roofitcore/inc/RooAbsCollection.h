/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooAbsCollection.h,v 1.26 2007/08/09 19:55:47 wouter Exp $
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
#ifndef ROO_ABS_COLLECTION
#define ROO_ABS_COLLECTION

#include "TObject.h"
#include "TString.h"
#include "RooAbsArg.h"
#include "RooPrintable.h"
#include "RooCmdArg.h"
#include "RooLinkedListIter.h"
#include <string>
#include <unordered_map>
#include <vector>

class RooCmdArg;

class RooAbsCollection : public TObject, public RooPrintable {
public:
  using Storage_t = std::vector<RooAbsArg*>;
  using const_iterator = Storage_t::const_iterator;


  // Constructors, assignment etc.
  RooAbsCollection();
  RooAbsCollection(const char *name);
  virtual TObject* clone(const char* newname) const = 0 ;
  virtual TObject* create(const char* newname) const = 0 ;
  virtual TObject* Clone(const char* newname=0) const { 
    return clone(newname?newname:GetName()) ; 
  }
  virtual ~RooAbsCollection();

  // Create a copy of an existing list. New variables cannot be added
  // to a copied list. The variables in the copied list are independent
  // of the original variables.
  RooAbsCollection(const RooAbsCollection& other, const char *name="");
  RooAbsCollection& operator=(const RooAbsCollection& other);
  RooAbsCollection& assignValueOnly(const RooAbsCollection& other, Bool_t oneSafe=kFALSE) ;
  void assignFast(const RooAbsCollection& other, Bool_t setValDirty=kTRUE) ;

  // Copy list and contents (and optionally 'deep' servers)
  RooAbsCollection *snapshot(Bool_t deepCopy=kTRUE) const ;
  Bool_t snapshot(RooAbsCollection& output, Bool_t deepCopy=kTRUE) const ;

  /// Set the size at which the collection will automatically start using an extra
  /// lookup table instead of performing a linear search.
  void setHashTableSize(Int_t number) {
    _sizeThresholdForMapSearch = number;
  }
  /// Query the size at which the collection will automatically start using an extra
  /// lookup table instead of performing a linear search.
  Int_t getHashTableSize() const { 
    return _sizeThresholdForMapSearch;
  }

  // List content management
  virtual Bool_t add(const RooAbsArg& var, Bool_t silent=kFALSE) ;
  virtual Bool_t addOwned(RooAbsArg& var, Bool_t silent=kFALSE);
  virtual RooAbsArg *addClone(const RooAbsArg& var, Bool_t silent=kFALSE) ;
  virtual Bool_t replace(const RooAbsArg& var1, const RooAbsArg& var2) ;
  virtual Bool_t remove(const RooAbsArg& var, Bool_t silent=kFALSE, Bool_t matchByNameOnly=kFALSE) ;
  virtual void removeAll() ;

  virtual Bool_t add(const RooAbsCollection& list, Bool_t silent=kFALSE) ;
  virtual Bool_t addOwned(const RooAbsCollection& list, Bool_t silent=kFALSE);
  virtual void   addClone(const RooAbsCollection& list, Bool_t silent=kFALSE);
  Bool_t replace(const RooAbsCollection &other);
  Bool_t remove(const RooAbsCollection& list, Bool_t silent=kFALSE, Bool_t matchByNameOnly=kFALSE) ;
  template<class forwardIt>
  void remove(forwardIt rangeBegin, forwardIt rangeEnd, Bool_t silent = kFALSE, Bool_t matchByNameOnly = kFALSE) {
      for (forwardIt it = rangeBegin; it != rangeEnd; ++it) {
        static_assert(std::is_same<
            typename std::iterator_traits<forwardIt>::value_type,
            RooAbsArg*>::value, "Can only remove lists of RooAbsArg*.");
        auto castedElm = static_cast<RooAbsArg*>(*it);
        remove(*castedElm, silent, matchByNameOnly);
      }
  }

  // Group operations on AbsArgs
  void setAttribAll(const Text_t* name, Bool_t value=kTRUE) ;

  // List search methods
  RooAbsArg *find(const char *name) const ;
  RooAbsArg *find(const RooAbsArg&) const ;

  /// Check if collection contains an argument with the same name as var.
  /// To check for a specific instance, use containsInstance().
  Bool_t contains(const RooAbsArg& var) const { 
    return find(var) != nullptr;
  }
  /// Check if this exact instance is in this collection.
  Bool_t containsInstance(const RooAbsArg& var) const { 
    return std::find(_list.begin(), _list.end(), &var) != _list.end();
  }
  RooAbsCollection* selectByAttrib(const char* name, Bool_t value) const ;
  RooAbsCollection* selectCommon(const RooAbsCollection& refColl) const ;
  RooAbsCollection* selectByName(const char* nameList, Bool_t verbose=kFALSE) const ;
  Bool_t equals(const RooAbsCollection& otherColl) const ; 
  Bool_t overlaps(const RooAbsCollection& otherColl) const ;

  /// TIterator-style iteration over contained elements.
  /// \note These iterators are slow. Use begin() and end() or
  /// range-based for loop instead.
  inline TIterator* createIterator(Bool_t dir = kIterForward) const
  R__SUGGEST_ALTERNATIVE("begin(), end() and range-based for loops.") {
    // Create and return an iterator over the elements in this collection
    return new RooLinkedListIter(makeLegacyIterator(dir));
  }

  /// TIterator-style iteration over contained elements.
  /// \note This iterator is slow. Use begin() and end() or range-based for loop instead.
  RooLinkedListIter iterator(Bool_t dir = kIterForward) const
  R__SUGGEST_ALTERNATIVE("begin(), end() and range-based for loops.") {
    return RooLinkedListIter(makeLegacyIterator(dir));
  }

  /// One-time forward iterator.
  /// \note Use begin() and end() or range-based for loop instead.
  RooFIter fwdIterator() const
  R__SUGGEST_ALTERNATIVE("begin(), end() and range-based for loops.") {
    return RooFIter(makeLegacyIterator());
  }

  const_iterator begin() const {
    return _list.begin();
  }
  
  const_iterator end() const {
    return _list.end();
  }

  Storage_t::const_reverse_iterator rbegin() const {
    return _list.rbegin();
  }

  Storage_t::const_reverse_iterator rend() const {
      return _list.rend();
    }

  Storage_t::size_type size() const {
    return _list.size();
  }

  bool empty() const {
    return _list.empty();
  }

  void reserve(Storage_t::size_type count) {
    _list.reserve(count);
  }

  /// Clear contents. If the collection is owning, it will also delete the contents.
  void clear() {
    removeAll();
  }

  inline Int_t getSize() const { 
    // Return the number of elements in the collection
    return _list.size();
  }
  
  inline RooAbsArg *first() const { 
    // Return the first element in this collection
    return _list.front();
  }

  RooAbsArg * operator[](Storage_t::size_type i) const {
    return _list[i];
  }


  /// Returns index of given arg, or -1 if arg is not in the collection.
  inline Int_t index(const RooAbsArg* arg) const {
    auto item = std::find(_list.begin(), _list.end(), arg);
    return item != _list.end() ? item - _list.begin() : -1;
  }

  /// Returns index of given arg, or -1 if arg is not in the collection.
  inline Int_t index(const RooAbsArg& arg) const {
    return index(&arg);
  }

  Int_t index(const char* name) const;

  inline virtual void Print(Option_t *options= 0) const {
    // Printing interface (human readable)
    printStream(defaultPrintStream(),defaultPrintContents(options),defaultPrintStyle(options));
  }
  std::string contentsString() const ;


  virtual void printName(std::ostream& os) const ;
  virtual void printTitle(std::ostream& os) const ;
  virtual void printClassName(std::ostream& os) const ;
  virtual void printValue(std::ostream& os) const ;
  virtual void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose=kFALSE, TString indent="") const ;

  virtual Int_t defaultPrintContents(Option_t* opt) const ;

  // Latex printing methods
  void printLatex(const RooCmdArg& arg1=RooCmdArg(), const RooCmdArg& arg2=RooCmdArg(),	
		  const RooCmdArg& arg3=RooCmdArg(), const RooCmdArg& arg4=RooCmdArg(),	
		  const RooCmdArg& arg5=RooCmdArg(), const RooCmdArg& arg6=RooCmdArg(),	
		  const RooCmdArg& arg7=RooCmdArg(), const RooCmdArg& arg8=RooCmdArg()) const ;
  void printLatex(std::ostream& ofs, Int_t ncol, const char* option="NEYU", Int_t sigDigit=1, 
                  const RooLinkedList& siblingLists=RooLinkedList(), const RooCmdArg* formatCmd=0) const ;

  void setName(const char *name) {
    // Set name of collection
    _name= name; 
  }
  const char* GetName() const { 
    // Return namer of collection
    return _name.Data() ; 
  }
  Bool_t isOwning() const { 
    // Does collection own contents?
    return _ownCont ; 
  }

  Bool_t allInRange(const char* rangeSpec) const ;

  void dump() const ;

  void releaseOwnership() { _ownCont = kFALSE ; }
  void takeOwnership() { _ownCont = kTRUE ; }

  void sort(Bool_t reverse = false);

  virtual void RecursiveRemove(TObject *obj);

  void useHashMapForFind(bool flag) const;

protected:
  Storage_t _list; // Actual object storage
  using LegacyIterator_t = TIteratorToSTLInterface<Storage_t>;

  Bool_t _ownCont;  // Flag to identify a list that owns its contents.
  TString _name;    // Our name.
  Bool_t _allRRV ; // All contents are RRV

  void safeDeleteList() ;

  // Support for snapshot method 
  Bool_t addServerClonesToList(const RooAbsArg& var) ;

  inline TNamed* structureTag() { if (_structureTag==0) makeStructureTag() ; return _structureTag ; }
  inline TNamed* typedStructureTag() { if (_typedStructureTag==0) makeTypedStructureTag() ; return _typedStructureTag ; }

  mutable TNamed* _structureTag{nullptr}; //! Structure tag
  mutable TNamed* _typedStructureTag{nullptr}; //! Typed structure tag
  
  inline void clearStructureTags() { _structureTag=0 ; _typedStructureTag = 0 ; }

  void makeStructureTag() ;
  void makeTypedStructureTag() ;
  
private:
  std::unique_ptr<LegacyIterator_t> makeLegacyIterator (bool forward = true) const;
  mutable std::unique_ptr<std::unordered_map<const TNamed*, Storage_t::value_type>> _nameToItemMap; //!
  std::size_t _sizeThresholdForMapSearch; //!
  void insert(RooAbsArg*);
  RooAbsArg* tryFastFind(const TNamed* namePtr) const;

  ClassDef(RooAbsCollection,3) // Collection of RooAbsArg objects
};

#endif
