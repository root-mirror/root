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
\class RooCategorySharedProperties
\ingroup Roofitlegacy

RooCategorySharedProperties is the container for all properties
that are shared between instance of RooCategory objects that
are clones of each other. At present the only property that is
shared in this way is the list of alternate named range definitions
**/

#include "RooFitLegacy/RooCategorySharedProperties.h"

#include "RooFit.h"
#include "RooFitLegacy/RooCatTypeLegacy.h"
#include "TList.h"
#include <iostream>


using namespace std ;

ClassImp(RooCategorySharedProperties);
;


////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooCategorySharedProperties::RooCategorySharedProperties()
{
} 


////////////////////////////////////////////////////////////////////////////////
/// Constructor with unique-id string

RooCategorySharedProperties::RooCategorySharedProperties(const char* uuidstr) : RooSharedProperties(uuidstr)
{
} 





////////////////////////////////////////////////////////////////////////////////

RooCategorySharedProperties::RooCategorySharedProperties(const RooCategorySharedProperties& other) :
  RooSharedProperties(other)
{
  cout << "RooCategorySharedProperties::cctor()" << endl ;
  // Copy constructor
  TIterator* iter = other._altRanges.MakeIterator() ;
  TList* olist ;
  while((olist=(TList*)iter->Next())) {
    TList* mylist = new TList ; 
    mylist->SetName(olist->GetName()) ;
    RooCatType* ctype ;
    TIterator* citer = olist->MakeIterator() ;    
    while ((ctype=(RooCatType*)citer->Next())) {
      mylist->Add(new RooCatType(*ctype)) ;
    }
    delete citer ;
    mylist->SetOwner(kTRUE) ;
    _altRanges.Add(mylist) ;
  }
  delete iter ;
}





////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooCategorySharedProperties::~RooCategorySharedProperties() 
{
  _altRanges.Delete() ;
} 


