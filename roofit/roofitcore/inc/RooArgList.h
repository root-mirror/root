/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooArgList.h,v 1.14 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_ARG_LIST
#define ROO_ARG_LIST

#include "RooAbsCollection.h"

class RooArgList : public RooAbsCollection {
public:

  // Constructors, assignment etc.
  RooArgList();
  RooArgList(const RooArgSet& set) ;
  explicit RooArgList(const TCollection& tcoll, const char* name="") ;
  explicit RooArgList(const char *name);
  /// Construct a (non-owning) RooArgList from one or more
  /// RooFit objects.
  /// \param arg A RooFit object to be put in the set.
  /// \param varsOrName Arbitrary number of
  ///   - RooFit objects deriving from RooAbsArg.
  ///   - A c-string to name the set.
  template<typename... Arg_t>
  RooArgList(const RooAbsArg& arg, const Arg_t&... argsOrName)
  /*NB: Making this a delegating constructor led to linker errors with MSVC*/ {
    processArg(arg);
    // Expand parameter pack in C++ 11 way:
    int dummy[] = { 0, (processArg(argsOrName), 0) ... };
    (void)dummy;
  };

  /// Construct from iterators.
  /// \tparam Iterator_t An iterator pointing to RooFit objects or references thereof.
  /// \param beginIt Iterator to first element to add.
  /// \param end Iterator to end of range to be added.
  /// \param name Optional name of the collection.
  template<typename Iterator_t,
      typename value_type = typename std::iterator_traits<Iterator_t>::value_type,
      typename = std::enable_if<std::is_convertible<const value_type*, const RooAbsArg*>::value> >
  RooArgList(Iterator_t beginIt, Iterator_t endIt, const char* name="") :
  RooArgList(name) {
    for (auto it = beginIt; it != endIt; ++it) {
      add(*it);
    }
  }

  virtual ~RooArgList();
  // Create a copy of an existing list. New variables cannot be added
  // to a copied list. The variables in the copied list are independent
  // of the original variables.
  RooArgList(const RooArgList& other, const char *name="");
  virtual TObject* clone(const char* newname) const { return new RooArgList(*this,newname); }
  virtual TObject* create(const char* newname) const { return new RooArgList(newname); }
  RooArgList& operator=(const RooArgList& other) { RooAbsCollection::operator=(other) ; return *this ; }


  /// Return object at given index, or nullptr if index is out of range
  inline RooAbsArg* at(Int_t idx) const { 

    if (idx >= static_cast<Int_t>(_list.size()))
      return nullptr;

    return _list[idx];
  }

  // I/O streaming interface (machine readable)
  virtual Bool_t readFromStream(std::istream& is, Bool_t compact, Bool_t verbose=kFALSE) ;
  virtual void writeToStream(std::ostream& os, Bool_t compact) ;  

  RooAbsArg& operator[](Int_t idx) const ; 

private:
  void processArg(const RooAbsArg& arg) { add(arg); }
  void processArg(const char* name) { _name = name; }

  ClassDef(RooArgList,1) // Ordered list of RooAbsArg objects
};

#endif
