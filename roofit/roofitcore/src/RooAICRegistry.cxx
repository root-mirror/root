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
\file RooAICRegistry.cxx
\class RooAICRegistry
\ingroup Roofitcore

RooAICRegistry is a utility class for operator p.d.f
classes that keeps track of analytical integration codes and
associated normalization and integration sets.  
**/

#include "RooFit.h"

#include "RooAICRegistry.h"
#include "RooMsgService.h"
#include "RooArgSet.h"

#include "Riostream.h"


using namespace std;

ClassImp(RooAICRegistry);
;

////////////////////////////////////////////////////////////////////////////////

RooAICRegistry::RooAICRegistry(UInt_t size) 
  : _clArr(0), _asArr1(0), _asArr2(0), _asArr3(0), _asArr4(0)
{
  _clArr.reserve(size);
  _asArr1.reserve(size);
  _asArr2.reserve(size);
  _asArr3.reserve(size);
  _asArr4.reserve(size);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooAICRegistry::RooAICRegistry(const RooAICRegistry& other)
  : _clArr(other._clArr), _asArr1(other._clArr.size(), 0), _asArr2(other._clArr.size(), 0),
    _asArr3(other._clArr.size(), 0), _asArr4(other._clArr.size(), 0)
{
  // Copy code-list array if other PDF has one
  UInt_t size = other._clArr.size();
  if (size) {
    _asArr1.resize(size, 0);
    _asArr2.resize(size, 0);
    _asArr3.resize(size, 0);
    _asArr4.resize(size, 0);
    for(UInt_t i = 0; i < size; ++i) {
      _asArr1[i] = other._asArr1[i] ? ((RooArgSet*)other._asArr1[i]->snapshot(kFALSE)) : 0; 
      _asArr2[i] = other._asArr2[i] ? ((RooArgSet*)other._asArr2[i]->snapshot(kFALSE)) : 0;
      _asArr3[i] = other._asArr3[i] ? ((RooArgSet*)other._asArr3[i]->snapshot(kFALSE)) : 0;
      _asArr4[i] = other._asArr4[i] ? ((RooArgSet*)other._asArr4[i]->snapshot(kFALSE)) : 0;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooAICRegistry::~RooAICRegistry() 
{
  // Delete code list array, if allocated
  for (unsigned int i = 0; i < _clArr.size(); ++i) {
    if (_asArr1[i]) delete   _asArr1[i];
    if (_asArr2[i]) delete   _asArr2[i];
    if (_asArr3[i]) delete   _asArr3[i];
    if (_asArr4[i]) delete   _asArr4[i];
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Store given arrays of integer codes, and up to four RooArgSets in
/// the registry (each setX pointer may be null). The registry
/// clones all RooArgSets internally so the RooArgSets passed as
/// arguments do not need to live beyond the store() call. The return
/// value is a unique master code for the given configuration of
/// integers and RooArgSets. If an identical combination is
/// previously stored in the registry no objects are stored and the
/// unique code of the existing entry is returned.

Int_t RooAICRegistry::store(const std::vector<Int_t>& codeList, RooArgSet* set1,
                            RooArgSet* set2, RooArgSet* set3, RooArgSet* set4)
{
  // Loop over code-list array  
  for (UInt_t i = 0; i < _clArr.size(); ++i) {
    // Existing slot, compare with current list, if matched return index
    Bool_t match(kTRUE) ;
    
    // Check that array contents is identical
    match &= _clArr[i] == codeList;

    // Check that supplied configuration of lists is identical
    if (_asArr1[i] && !set1) match=kFALSE ;
    if (!_asArr1[i] && set1) match=kFALSE ;
    if (_asArr2[i] && !set2) match=kFALSE ;
    if (!_asArr2[i] && set2) match=kFALSE ;
    if (_asArr3[i] && !set3) match=kFALSE ;
    if (!_asArr3[i] && set3) match=kFALSE ;
    if (_asArr4[i] && !set4) match=kFALSE ;
    if (!_asArr4[i] && set4) match=kFALSE ;
    
    // Check that contents of arrays is identical
    if (_asArr1[i] && set1 && !set1->equals(*_asArr1[i])) match=kFALSE ;
    if (_asArr2[i] && set2 && !set2->equals(*_asArr2[i])) match=kFALSE ;	
    if (_asArr3[i] && set3 && !set3->equals(*_asArr3[i])) match=kFALSE ;	
    if (_asArr4[i] && set4 && !set4->equals(*_asArr4[i])) match=kFALSE ;	
    
    if (match) {
      if (set1) delete set1 ;
      if (set2) delete set2 ;
      if (set3) delete set3 ;
      if (set4) delete set4 ;
      return i ;
    }
  }

  // Store code list and return index
  _clArr.push_back(codeList);
  _asArr1.push_back(set1 ? (RooArgSet*)set1->snapshot(kFALSE) : 0);
  _asArr2.push_back(set2 ? (RooArgSet*)set2->snapshot(kFALSE) : 0);
  _asArr3.push_back(set3 ? (RooArgSet*)set3->snapshot(kFALSE) : 0);
  _asArr4.push_back(set4 ? (RooArgSet*)set4->snapshot(kFALSE) : 0);

  if (set1) delete set1 ;
  if (set2) delete set2 ;
  if (set3) delete set3 ;
  if (set4) delete set4 ;
  return _clArr.size() - 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve the array of integer codes associated with the given master code

const std::vector<Int_t>& RooAICRegistry::retrieve(Int_t masterCode) const 
{
  return _clArr[masterCode] ;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve the array of integer codes associated with the given master code
/// and set the passed set1 pointer to the first RooArgSet associated with this master code

const std::vector<Int_t>& RooAICRegistry::retrieve(Int_t masterCode, pRooArgSet& set1) const 
{
  set1 = _asArr1[masterCode] ;
  return _clArr[masterCode] ;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve the array of integer codes associated with the given master code
/// and set the passed set1,set2 pointers to the first and second  RooArgSets associated with this 
/// master code respectively

const std::vector<Int_t>& RooAICRegistry::retrieve
(Int_t masterCode, pRooArgSet& set1, pRooArgSet& set2) const 
{
  set1 = _asArr1[masterCode] ;
  set2 = _asArr2[masterCode] ;
  return _clArr[masterCode] ;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve the array of integer codes associated with the given master code
/// and set the passed set1-4 pointers to the four  RooArgSets associated with this 
/// master code respectively

const std::vector<Int_t>& RooAICRegistry::retrieve
(Int_t masterCode, pRooArgSet& set1, pRooArgSet& set2, pRooArgSet& set3, pRooArgSet& set4) const 
{
  set1 = _asArr1[masterCode] ;
  set2 = _asArr2[masterCode] ;
  set3 = _asArr3[masterCode] ;
  set4 = _asArr4[masterCode] ;
  return _clArr[masterCode] ;
}
