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

//////////////////////////////////////////////////////////////////////////////
/// \class RooArgList
/// RooArgList is a container object that can hold multiple RooAbsArg objects.
/// The container has list semantics which means that:
///
///  - Contained objects are ordered, The iterator
///    follows the object insertion order.
///
///  - Objects can be retrieved by name and index
///
///  - Multiple objects with the same name are allowed
///
/// Ownership of contents.
///
/// Unowned objects are inserted with the add() method. Owned objects
/// are added with addOwned() or addClone(). A RooArgSet either owns all
/// of it contents, or none, which is determined by the first <add>
/// call. Once an ownership status is selected, inappropriate <add> calls
/// will return error status. Clearing the list via removeAll() resets the
/// ownership status. Arguments supplied in the constructor are always added
/// as unowned elements.
///
///

#include "RooArgList.h"

#include "RooStreamParser.h"
#include "RooAbsRealLValue.h"
#include "RooAbsCategoryLValue.h"
#include "RooTrace.h"
#include "RooMsgService.h"

#include <stdexcept>

using namespace std;

ClassImp(RooArgList);


////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooArgList::RooArgList() :
  RooAbsCollection()
{
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor from a RooArgSet. 

RooArgList::RooArgList(const RooArgSet& set) :
  RooAbsCollection(set.GetName())
{
  add(set) ;
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Empty list constructor

RooArgList::RooArgList(const char *name) :
  RooAbsCollection(name)
{
  TRACE_CREATE
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor from a root TCollection. Elements in the collection that
/// do not inherit from RooAbsArg will be skipped. A warning message
/// will be printed for every skipped item.

RooArgList::RooArgList(const TCollection& tcoll, const char* name) :
  RooAbsCollection(name)
{
  TIterator* iter = tcoll.MakeIterator() ;
  TObject* obj ;
  while((obj=iter->Next())) {
    if (!dynamic_cast<RooAbsArg*>(obj)) {
      coutW(InputArguments) << "RooArgList::RooArgList(TCollection) element " << obj->GetName() 
			    << " is not a RooAbsArg, ignored" << endl ;
      continue ;
    }
    add(*(RooAbsArg*)obj) ;
  }
  delete iter ;
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor. Note that a copy of a list is always non-owning,
/// even the source list is owning. To create an owning copy of
/// a list (owning or not), use the snaphot() method.

RooArgList::RooArgList(const RooArgList& other, const char *name) 
  : RooAbsCollection(other,name)
{
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooArgList::~RooArgList() 
{
  TRACE_DESTROY
}



////////////////////////////////////////////////////////////////////////////////
/// Array operator. Element in slot `idx` must exist.
/// \throws std::invalid_argument if `idx` is out of range.
///
/// When used as
/// ```
///  myArgList[4] = x;
/// ```
/// note that the element contained in the list will not be
/// replaced! Instead, `operator=` of the existing element is called.
RooAbsArg& RooArgList::operator[](Int_t idx) const 
{     
  RooAbsArg* arg = at(idx) ;
  if (!arg) {
    coutE(InputArguments) << "RooArgList::operator[](" << GetName() << ") ERROR: index " 
			  << idx << " out of range (0," << getSize() << ")" << endl ;
    throw std::invalid_argument(std::string("Index ") + to_string(idx) + " is out of range.");
  }
  return *arg ; 
}



////////////////////////////////////////////////////////////////////////////////
/// Write the contents of the argset in ASCII form to given stream.
/// 
/// All elements will be printed on a single line separated by a single 
/// white space. The contents of each element is written by the arguments' 
/// writeToStream() function

void RooArgList::writeToStream(ostream& os, Bool_t compact) 
{
  if (!compact) {
    coutE(InputArguments) << "RooArgList::writeToStream(" << GetName() << ") non-compact mode not supported" << endl ;
    return ;
  }

  for (const auto obj : _list) {
    obj->writeToStream(os,kTRUE);
    os << " " ;
  }
  os << endl ;
}



////////////////////////////////////////////////////////////////////////////////
/// Read the contents of the argset in ASCII form from given stream.
/// 
/// A single line is read, and all elements are assumed to be separated 
/// by white space. The value of each argument is read by the arguments 
/// readFromStream function.

Bool_t RooArgList::readFromStream(istream& is, Bool_t compact, Bool_t verbose) 
{
  if (!compact) {
    coutE(InputArguments) << "RooArgList::readFromStream(" << GetName() << ") non-compact mode not supported" << endl ;
    return kTRUE ;
  }    

  RooStreamParser parser(is) ;
  for (auto next : _list) {
    if (!next->getAttribute("Dynamic")) {
      if (next->readFromStream(is,kTRUE,verbose)) {
        parser.zapToEnd() ;

        return kTRUE ;
      }	
    } else {
    }
  }
  
  if (!parser.atEOL()) {
    TString rest = parser.readLine() ;
    if (verbose) {
      coutW(InputArguments) << "RooArgSet::readFromStream(" << GetName() 
			    << "): ignoring extra characters at end of line: '" << rest << "'" << endl ;
    }
  }

  return kFALSE ;  
}

