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
\file RooMultiCatIter.cxx
\class RooMultiCatIter
\ingroup Roofitcore

RooMultiCatIter iterators over all state permutations of a list of categories.
It serves as the state iterator for a RooSuperCategory or a RooMultiCategory.
Since this iterator only constructs state labels and does not change the value
of its input categories, it is not required that its inputs are LValues.  
For cases where all inputs are LValues (such as for RooSuperCategory) the
values of the input can be changes by assigning the super category the
string label generated by this iterator
**/

#include "RooFit.h"

#include "RooAbsCategoryLValue.h"
#include "RooAbsCategoryLValue.h"
#include "RooMultiCatIter.h"

using namespace std;

ClassImp(RooMultiCatIter);
;



////////////////////////////////////////////////////////////////////////////////
/// Construct iterator over all permutations of states of categories in catList.
/// If rangeName is not null, iteration is restricted to states that are selected
/// in the given range name

RooMultiCatIter::RooMultiCatIter(const RooArgSet& catList, const char* rangeName) : _catList("catList") 
{
  if (rangeName) {
    _rangeName = rangeName ;
  }
  initialize(catList) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooMultiCatIter::RooMultiCatIter(const RooMultiCatIter& other) : TIterator(other), _catList("catList")
{
  initialize(other._catList) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Build iterator array for given catList

void RooMultiCatIter::initialize(const RooArgSet& catList) 
{
  // Copy RooCategory list into internal argset
  TIterator* catIter = catList.createIterator() ;
  TObject* obj ;
  while ((obj = catIter->Next())) {
    _catList.add((RooAbsArg&)*obj) ;
  }
  delete catIter ;
  
  // Allocate storage for component iterators
  _nIter = catList.getSize() ;
  _iterList   = new pTIterator[_nIter] ;
  _catPtrList = new pRooCategory[_nIter] ;
  _curTypeList = new RooCatType[_nIter] ;

  // Construct component iterators
  _curIter = 0 ;
  _curItem = nullptr;
  TIterator* cIter = _catList.createIterator() ;
  RooAbsCategoryLValue* cat ;
  while((cat=(RooAbsCategoryLValue*)cIter->Next())) {
    _catPtrList[_curIter] = cat ;
    _iterList[_curIter] = cat->typeIterator() ;
    _iterList[_curIter]->Next() ;
//     _curTypeList[_curIter] = *first ;
//     _curTypeList[_curIter].SetName(first->GetName()) ;
//     cout << "init: _curTypeList[" << _curIter << "] set to " << first->GetName() << endl ;
//     _iterList[_curIter]->Reset() ;    
    _curIter++ ;
  }
  delete cIter ;

  Reset() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooMultiCatIter::~RooMultiCatIter() 
{
  for (_curIter=0 ; _curIter<_nIter ; _curIter++) {
    delete _iterList[_curIter] ;
  }
  delete[] _iterList ;
  delete[] _catPtrList ;
  delete[] _curTypeList ;
}



////////////////////////////////////////////////////////////////////////////////
/// Dummy implementation, always returns zero

const TCollection* RooMultiCatIter::GetCollection() const 
{
  //return &_catList.getCollection() ;
  return nullptr;
}



////////////////////////////////////////////////////////////////////////////////
/// Construct string with composite object
/// label corresponding to the state name
/// of a RooMultiCategory or RooSuperCategory
/// constructed from this set of input categories

TObjString* RooMultiCatIter::compositeLabel() 
{
  TString& str = _compositeLabel.String() ;

  str = "{" ;
  Int_t i ;
  for (i=0 ; i<_nIter ; i++) {
    if (i>0) str.Append(";") ;
    str.Append(_curTypeList[i].GetName()) ;
  }
  str.Append("}") ;

  return &_compositeLabel ;
}



////////////////////////////////////////////////////////////////////////////////
/// Iterator increment operator

TObject* RooMultiCatIter::Next() 
{
  // Check for end
  if (_curIter==_nIter) {
     _curItem = nullptr;
     return nullptr;
  }

  RooCatType* next = (RooCatType*) _iterList[_curIter]->Next() ;
  if (next) { 

    // Increment current iterator
    _curTypeList[_curIter] = *next ;
    _curTypeList[_curIter].SetName(next->GetName()) ;
    
    // If higher order increment was successful, reset master iterator
    if (_curIter>0) _curIter=0 ;

    _curItem = compositeLabel() ;
    return _curItem ;    
  } else {

    // Reset current iterator
    _iterList[_curIter]->Reset() ;
    next = (RooCatType*) _iterList[_curIter]->Next() ;
    if (next) {
      _curTypeList[_curIter] = *next ; 
      _curTypeList[_curIter].SetName(next->GetName()) ;
    }
    //if (next) _catPtrList[_curIter]->setIndex(next->getVal()) ;

    // Increment next iterator 
    _curIter++ ;
    _curItem = Next() ;
    return _curItem ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Rewind master iterator

void RooMultiCatIter::Reset() 
{
  for (_curIter=0 ; _curIter<_nIter ; _curIter++) {
    TIterator* cIter = _iterList[_curIter] ;
    cIter->Reset() ;
    RooCatType* first = (RooCatType*) cIter->Next() ;
    if (first) {
      if (_curIter==0) cIter->Reset() ;
      _curTypeList[_curIter] = *first ;
      _curTypeList[_curIter].SetName(first->GetName()) ;
    }
  }
  _curIter=0 ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return current item (dummy)

TObject *RooMultiCatIter::operator*() const
{
  return _curItem ;
}


////////////////////////////////////////////////////////////////////////////////
/// Comparison operator to other iterator
/// Returns true if both iterator iterate over the
/// same set of input categories and are not at the
/// same sequential position

bool RooMultiCatIter::operator!=(const TIterator &aIter) const
{
   if ((aIter.IsA() == RooMultiCatIter::Class())) {
      const RooMultiCatIter &iter(dynamic_cast<const RooMultiCatIter &>(aIter));
      return (_curItem != iter._curItem);
   }
   
   return false;
}

