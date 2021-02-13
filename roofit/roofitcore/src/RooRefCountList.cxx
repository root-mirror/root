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
\file RooRefCountList.cxx
\class RooRefCountList
\ingroup Roofitcore

A RooRefCountList is a RooLinkedList that keeps a reference counter
with each added node. Multiple Add()s of the same object will increase
the counter instead of adding multiple copies. Remove() decrements the 
reference count until zero, when the object is actually removed.
**/

#include "RooFit.h"

#include "RooRefCountList.h"

#include "Riostream.h"
#include <stdlib.h>

using namespace std;

ClassImp(RooRefCountList);
  ;



////////////////////////////////////////////////////////////////////////////////
/// Default constructor construct lists with initial hash table size of 17

RooRefCountList::RooRefCountList()
  : RooLinkedList(0) 
{ 
}



////////////////////////////////////////////////////////////////////////////////
/// Add object to list with given reference count increment
/// List takes ownership of object.

void RooRefCountList::Add(TObject* obj, Int_t count) 
{
  // Check if we already have it
  TObject* listObj = FindObject(obj) ;
  if (!listObj) {
    // Add to list with reference count 
    RooLinkedList::Add(obj, count) ;
    //cout << "RooRefCountList::AddLast(" << obj << ") adding object" << endl ;
  } else {
    RooLinkedListElem* link = findLink(obj) ;
    if(link) {
      while(count--) link->incRefCount() ;    
    }
    //cout << "RooRefCountList::AddLast(" << obj << ") incremented reference count to " << link->refCount() << endl ;
  }

}



////////////////////////////////////////////////////////////////////////////////
/// Remove object from list and if reference count 
/// reaches zero delete object itself as well.

Bool_t RooRefCountList::Remove(TObject* obj) 
{
  RooLinkedListElem* link = findLink(obj) ;
  if (!link) {
    return 0 ;
  } else {
    if (link->decRefCount()==0) {
      //cout << "RooRefCountList::AddLast(" << obj << ") removed object" << endl ;
      return RooLinkedList::Remove(obj) ;
    }
    //cout << "RooRefCountList::AddLast(" << obj << ") decremented reference count to " << link->refCount() << endl ;
  }
  return 0 ;
}



////////////////////////////////////////////////////////////////////////////////
/// Remove object from list and delete object itself
/// regardless of reference count

Bool_t RooRefCountList::RemoveAll(TObject* obj)
{
  return RooLinkedList::Remove(obj) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return reference count associated with 'obj'

Int_t RooRefCountList::refCount(TObject* obj) const
{
  RooLinkedListElem* link = findLink(obj) ;
  if (!link) {
    return 0 ;
  } else {
    return link->refCount() ;
  }
}
