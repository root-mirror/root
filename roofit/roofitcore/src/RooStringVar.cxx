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
\file RooStringVar.cxx
\class RooStringVar
\ingroup Roofitcore

RooStringVar is a RooAbsArg implementing string values.
**/

#include "RooStringVar.h"

#include "RooFit.h"
#include "Riostream.h"
#include "TTree.h"
#include "RooStreamParser.h"
#include "RooMsgService.h"
#include "TBranch.h"


////////////////////////////////////////////////////////////////////////////////
/// Constructor with initial value. The size argument is ignored.
RooStringVar::RooStringVar(const char *name, const char *title, const char* value, Int_t) :
  RooAbsArg(name, title),
  _string(value)
{
  setValueDirty();
}  



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooStringVar::RooStringVar(const RooStringVar& other, const char* name) :
  RooAbsArg(other, name),
  _string(other._string)
{
  setValueDirty();
}



////////////////////////////////////////////////////////////////////////////////
/// Read object contents from given stream
bool RooStringVar::readFromStream(std::istream& is, Bool_t compact, Bool_t)
{
  TString token,errorPrefix("RooStringVar::readFromStream(") ;
  errorPrefix.Append(GetName()) ;
  errorPrefix.Append(")") ;
  RooStreamParser parser(is,errorPrefix) ;

  TString newValue ;

  if (compact) {
    parser.readString(newValue,kTRUE) ;
  } else {
    newValue = parser.readLine() ;
  }

  _string = newValue;
  setValueDirty();

  return false;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy cache of another RooAbsArg to our cache
///
/// Warning: This function copies the cached values of source,
///          it is the callers responsibility to make sure the cache is clean

void RooStringVar::copyCache(const RooAbsArg* source, Bool_t /*valueOnly*/, Bool_t setValDirty)
{
  auto other = dynamic_cast<const RooStringVar*>(source) ;
  assert(other);

  _string = other->_string;
  if (setValDirty) {
    setValueDirty() ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Attach object to a branch of given TTree

void RooStringVar::attachToTree(TTree& t, Int_t)
{
  // First determine if branch is taken
  TBranch* branch ;
  if ((branch = t.GetBranch(GetName()))) {
    t.SetBranchAddress(GetName(), &_string);
  } else {
    t.Branch(GetName(), &_string);
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Fill tree branch associated with this object

void RooStringVar::fillTreeBranch(TTree& t)
{
  // First determine if branch is taken
  TBranch* branch = t.GetBranch(GetName()) ;
  if (!branch) {
    coutE(DataHandling) << "RooAbsString::fillTreeBranch(" << GetName() << ") ERROR: not attached to tree" << std::endl;
    assert(false);
    return;
  }
  branch->Fill() ;
}



////////////////////////////////////////////////////////////////////////////////
/// (De)Activate associated tree branch

void RooStringVar::setTreeBranchStatus(TTree& t, Bool_t active)
{
  TBranch* branch = t.GetBranch(GetName()) ;
  if (branch) {
    t.SetBranchStatus(GetName(),active?1:0) ;
  }
}


