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
  ///            Note that you can also pass a `double` as first argument
  ///            when constructing a RooArgList, and another templated
  ///            constructor will be used where a RooConstVar is implicitly
  ///            created from the `double` value.
  /// \param varsOrName Arbitrary number of
  ///   - RooFit objects deriving from RooAbsArg.
  ///   - `double`s from which a RooConstVar is implicitly created via `RooFit::RooConst`.
  ///   - A c-string to name the set.
  template<typename... Ts>
  RooArgList(const RooAbsArg& arg, const Ts&... moreArgsOrName)
  /*NB: Making this a delegating constructor led to linker errors with MSVC*/ {
    processArgs(arg, moreArgsOrName...);
  }

  template<typename... Ts>
  explicit RooArgList(double arg, const Ts&... moreArgsOrName) {
    processArgs(arg, moreArgsOrName...);
  }

  /// Construct from iterators.
  /// \tparam Iterator_t An iterator pointing to RooFit objects or pointers/references thereof.
  /// \param beginIt Iterator to first element to add.
  /// \param endIt Iterator to end of range to be added.
  /// \param name Optional name of the collection.
  template<typename Iterator_t,
      typename value_type = typename std::remove_pointer<typename std::iterator_traits<Iterator_t>::value_type>,
      typename = std::enable_if<std::is_convertible<const value_type*, const RooAbsArg*>::value> >
  RooArgList(Iterator_t beginIt, Iterator_t endIt, const char* name="") :
  RooArgList(name) {
    for (auto it = beginIt; it != endIt; ++it) {
      processArg(*it);
    }
  }

  /// Construct a non-owning RooArgList from a vector of RooAbsArg pointers.
  /// This constructor is mainly intended for pyROOT. With cppyy, a Python list
  /// or tuple can be implicitly converted to an std::vector, and by enabling
  /// implicit construction of a RooArgList from a std::vector, we indirectly
  /// enable implicit conversion from a Python list/tuple to RooArgLists.
  /// \param vec A vector with pointers to the arguments.
  RooArgList(std::vector<RooAbsArg*> const& vec) {
    for(auto const& arg : vec) {
      add(*arg);
    }
  }

  virtual ~RooArgList();
  // Create a copy of an existing list. New variables cannot be added
  // to a copied list. The variables in the copied list are independent
  // of the original variables.
  RooArgList(const RooArgList& other, const char *name="");
  /// Move constructor.
  RooArgList(RooArgList && other) : RooAbsCollection(std::move(other)) {}
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
  virtual bool readFromStream(std::istream& is, bool compact, bool verbose=false);
  virtual void writeToStream(std::ostream& os, bool compact); 

  /// Access element by index.
  RooAbsArg& operator[](Int_t idx) const {
    assert(0 <= idx && idx < static_cast<Int_t>(_list.size()));
    return *_list[idx];
  }

protected:
  virtual bool canBeAdded(RooAbsArg const&, bool) const  { return true; }

private:
  template<typename... Args_t>
  void processArgs(Args_t &&... args) {
    // Expand parameter pack in C++ 11 way:
    int dummy[] = { 0, (processArg(std::forward<Args_t>(args)), 0) ... };
    (void)dummy;
  }
  void processArg(const RooAbsArg& arg) { add(arg); }
  void processArg(const RooAbsArg* arg) { add(*arg); }
  void processArg(const char* name) { _name = name; }
  void processArg(double value);

  ClassDef(RooArgList,1) // Ordered list of RooAbsArg objects
};


namespace RooFitShortHand {

template<class... Args_t>
RooArgList L(Args_t&&... args) {
  return {std::forward<Args_t>(args)...};
}

} // namespace RooFitShortHand


#endif
