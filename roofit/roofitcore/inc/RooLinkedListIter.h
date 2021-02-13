/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooLinkedListIter.h,v 1.11 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_LINKED_LIST_ITER
#define ROO_LINKED_LIST_ITER

#include "TIterator.h"
#include "RooLinkedList.h"

#include <memory>
#include <stdexcept>
#include <assert.h>

/// Interface for RooFIter-compatible iterators
class GenericRooFIter
{
  public:
  /// Return next element or nullptr if at end.
  virtual RooAbsArg * next() = 0;
  virtual ~GenericRooFIter() {}
};

////////////////////////////////////////////////////////////////////////////////////////////
/// A one-time forward iterator working on RooLinkedList or RooAbsCollection.
/// This wrapper separates the interface visible to the outside from the actual
/// implementation of the iterator.
class RooFIter final
{
  public:
  RooFIter(std::unique_ptr<GenericRooFIter> && itImpl) : fIterImpl{std::move(itImpl)} {}
  RooFIter(const RooFIter &) = delete;
  RooFIter(RooFIter &&) = default;
  RooFIter & operator=(const RooFIter &) = delete;
  RooFIter & operator=(RooFIter &&) = default;

  /// Return next element or nullptr if at end.
  RooAbsArg *next() {
    return fIterImpl->next();
  }

  private:
  std::unique_ptr<GenericRooFIter> fIterImpl;
};

////////////////////////////////////////////////////////////////////////////////////////////
/// Implementation of the GenericRooFIter interface for the RooLinkedList
class RooFIterForLinkedList final : public GenericRooFIter
{
  public:
  RooFIterForLinkedList() {}
  RooFIterForLinkedList(const RooLinkedList* list) : fPtr (list->_first) {}

  /// Return next element in collection
  RooAbsArg *next() override {
    if (!fPtr) return nullptr ;
    TObject* arg = fPtr->_arg ;
    fPtr = fPtr->_next;
    return (RooAbsArg*) arg ;
  }

 private:
    const RooLinkedListElem * fPtr{nullptr};  //! Next link element
};




////////////////////////////////////////////////////////////////////////////////////////////
/// TIterator and GenericRooFIter front end with STL back end.
///
/// By default, this iterators counts, at which position the current element should be.
/// On request, it does an index access to the underlying collection, and returns the element.
/// This happens because the RooLinkedList, which used to be the default collection in RooFit,
/// will not invalidate iterators when inserting elements. Since the default is now an STL collection,
/// reallocations might invalidate the iterator.
///
/// With an iterator that counts, only inserting before or at the iterator position will create problems.
/// deal with reallocations while iterating. Therefore, this iterator will also check that the last element
/// it was pointing to is the the current element when it is invoked again. This ensures that
/// inserting or removing before this iterator does not happen, which was possible with
/// the linked list iterators of RooFit.
/// When NDEBUG is defined, these checks will disappear.
/// \note This is a legacy iterator that only exists to not break old code. Use begin(), end() and
/// range-based for loops with RooArgList and RooArgSet.
template<class STLContainer>
class TIteratorToSTLInterface final : public TIterator , public GenericRooFIter {
public:

  TIteratorToSTLInterface(const STLContainer & container) :
    TIterator(),
    GenericRooFIter(),
    fSTLContainer(container),
    fIndex(0)
#ifdef NDEBUG
    ,fCurrentElem{nullptr}
#else
    ,fCurrentElem{fSTLContainer.empty() ? nullptr : fSTLContainer.front()}
#endif
  {

  }

  TIterator & operator=(const TIterator &) override {
    throw;
  }

  const TCollection *GetCollection() const override {
    return nullptr;
  }

  RooAbsArg * next() override {
    if (atEnd())
      return nullptr;
#ifdef NDEBUG
    return fSTLContainer[fIndex++];
#else
    return nextChecked();
#endif
  }


  TObject * Next() override {
    return static_cast<TObject*>(next());
  }

  void Reset() override {
    fIndex = 0;
#ifndef NDEBUG
    fCurrentElem = fSTLContainer.empty() ? nullptr : fSTLContainer.front();
#endif

  }

  Bool_t operator!=(const TIterator & other) const override {
    const auto * castedOther =
        dynamic_cast<const TIteratorToSTLInterface<STLContainer>*>(&other);
    return !castedOther || &fSTLContainer != &(castedOther->fSTLContainer)
        || fIndex == castedOther->fIndex;
  }

  TObject * operator*() const override {
      if (atEnd())
        return nullptr;

    #ifndef NDEBUG
      assert(fCurrentElem == fSTLContainer[fIndex]);
    #endif

      return static_cast<TObject*>(fSTLContainer[fIndex]);
    }


private:
  bool atEnd() const {
    return fSTLContainer.empty()
        || fIndex >= fSTLContainer.size();
  }


  RooAbsArg * nextChecked() {
    RooAbsArg * ret = fSTLContainer.at(fIndex);
    if (fCurrentElem != nullptr && ret != fCurrentElem) {
      throw std::logic_error("A RooCollection should not be modified while iterating. "
          "Only inserting at end is acceptable.");
    }
    fCurrentElem = ++fIndex < fSTLContainer.size() ? fSTLContainer[fIndex] : nullptr;

    return ret;
  }


  const STLContainer & fSTLContainer; //!
  std::size_t fIndex; //!
  const RooAbsArg * fCurrentElem; //!
};




////////////////////////////////////////////////////////////////////////////////////////////
/// A wrapper around TIterator derivatives.
///
/// It is called RooLinkedListIter because all classes assume that the RooAbsCollections use
/// a RooLinkedList, which is not true, any more.
/// The purpose of this wrapper is to act on the outside like a RooLinkedListIter, even though
/// the underlying implementation may work an a different container, like e.g.
/// an STL container. This is needed to not break user code that is using a RooLinkedList or
/// a RooAbsCollection.
///
/// \note All code using this iterator as an iterator over a RooAbsCollection should move
/// to begin() and end() or range-based for loops. These are faster.
class RooLinkedListIter final : public TIterator {

  public:
  RooLinkedListIter(std::shared_ptr<TIterator> iterImpl) :
    fIterImpl{std::move(iterImpl)} {

  }

  RooLinkedListIter(const RooLinkedListIter &) = delete;
  RooLinkedListIter(RooLinkedListIter &&) = default;
  RooLinkedListIter & operator=(const RooLinkedListIter &) = delete;
  RooLinkedListIter & operator=(RooLinkedListIter &&) = default;

  TIterator &operator=(const TIterator & other) override {fIterImpl->operator=(other); return *this;}
  const TCollection *GetCollection() const override {return nullptr;}

  TObject * Next() override {return fIterImpl->Next();}
  void Reset() override {fIterImpl->Reset();}
  Bool_t operator!=(const TIterator & other) const override {return fIterImpl->operator!=(other);}
  TObject * operator*() const override {return fIterImpl->operator*();}

  private:
  std::shared_ptr<TIterator> fIterImpl; //!
};


////////////////////////////////////////////////////////////////////////////////////////////
/// Implementation of the actual iterator on RooLinkedLists.
///
class RooLinkedListIterImpl final : public TIterator {
public:

  RooLinkedListIterImpl() {
    // coverity[UNINIT_CTOR]
  } ;


  RooLinkedListIterImpl(const RooLinkedList* list, Bool_t forward) :
    TIterator(), _list(list), _ptr(forward ? _list->_first : _list->_last),
      _forward(forward)
  { }

  RooLinkedListIterImpl(const RooLinkedListIterImpl& other) :
    TIterator(other), _list(other._list), _ptr(other._ptr),
    _forward(other._forward)
  {
    // Copy constructor
  }

  virtual ~RooLinkedListIterImpl() { }

  TIterator& operator=(const TIterator& other) {

    // Iterator assignment operator

    if (&other==this) return *this ;
    const RooLinkedListIterImpl* iter = dynamic_cast<const RooLinkedListIterImpl*>(&other) ;
    if (iter) {
      _list = iter->_list ;
      _ptr = iter->_ptr ;
      _forward = iter->_forward ;
    }
    return *this ;
  }

  virtual const TCollection *GetCollection() const {
    // Dummy
    return 0 ;
  }

  virtual TObject *Next() {
    // Return next element in collection
    if (!_ptr) return 0 ;
    TObject* arg = _ptr->_arg ;
    _ptr = _forward ? _ptr->_next : _ptr->_prev ;
    return arg ;
  }

  TObject *NextNV() {
    // Return next element in collection
    if (!_ptr) return 0 ;
    TObject* arg = _ptr->_arg ;
    _ptr = _forward ? _ptr->_next : _ptr->_prev ;
    return arg ;
  }


  virtual void Reset() {
    // Return iterator to first element in collection
    _ptr = _forward ? _list->_first : _list->_last ;
  }

  bool operator!=(const TIterator &aIter) const {
    const RooLinkedListIterImpl *iter(dynamic_cast<const RooLinkedListIterImpl*>(&aIter));
    if (iter) return (_ptr != iter->_ptr);
    return false; // for base class we don't implement a comparison
  }

  bool operator!=(const RooLinkedListIterImpl &aIter) const {
    return (_ptr != aIter._ptr);
  }

  virtual TObject *operator*() const {
    // Return element iterator points to
    return (_ptr ? _ptr->_arg : nullptr);
  }

protected:
  const RooLinkedList* _list ;     //! Collection iterated over
  const RooLinkedListElem* _ptr ;  //! Next link element
  Bool_t _forward ;                //!  Iterator direction
};




#endif
