/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooLinkedListElem.h,v 1.11 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_LINKED_LIST_ELEM
#define ROO_LINKED_LIST_ELEM

#include "Rtypes.h"

class TObject ;
class RooLinkedListElem ;
class TBuffer ;

namespace RooLinkedListImplDetails {
    class Chunk;
    class Pool;
}

class RooLinkedListElem {
public:
  // Initial element ctor
  RooLinkedListElem() :
    _prev(0), _next(0), _arg(0), _refCount(0) {
  }

  void init(TObject* arg, RooLinkedListElem* after=0) {
   _arg = arg ;
   _refCount = 1 ;

   if (after) {
     _prev = after ;
     _next = after->_next ;
     after->_next = this ;
     if (_next) {
       _next->_prev = this ;
     }
   }
 }

 void release() {
   if (_prev) _prev->_next = _next ;
   if (_next) _next->_prev = _prev ;
   _prev = 0 ;
   _next = 0 ;
 }

  RooLinkedListElem(TObject* arg) :
    // Constructor with payload
    _prev(0), _next(0), _arg(arg), _refCount(1) {
  }

  RooLinkedListElem(TObject* arg, RooLinkedListElem* after) :
    // Constructor with payload and next chain element
    _prev(after), _next(after->_next), _arg(arg), _refCount(1) {

    // Insert self in link
    after->_next = this ;
    if (_next) _next->_prev = this ;
  }

  // Destructor
  virtual ~RooLinkedListElem() {
    // Remove self from link
    if (_prev) _prev->_next = _next ;
    if (_next) _next->_prev = _prev ;
  }

  Int_t refCount() const { return _refCount ; }
  Int_t incRefCount() { return ++_refCount ; }
  Int_t decRefCount() { return --_refCount ; }

protected:
  friend class RooHashTable ;
  friend class RooLinkedList ;
  friend class RooLinkedListImplDetails::Pool;
  friend class RooLinkedListImplDetails::Chunk;
  friend class RooLinkedListIterImpl ;
  friend class RooFIterForLinkedList ;
  RooLinkedListElem* _prev ; // Link to previous element in list
  RooLinkedListElem* _next ; // Link to next element in list
  TObject*   _arg ;          // Link to contents
  Int_t      _refCount ;     //! Reference count

protected:

  // Forbidden
  RooLinkedListElem(const RooLinkedListElem&) ;

  ClassDef(RooLinkedListElem,1) // Element of RooLinkedList container class
} ;



#endif
