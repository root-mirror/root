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
\file RooList.cxx
\class RooList
\ingroup Roofitcore

A RooList is a TList with extra support for working with options
that are associated with each node. This is a utility class for RooPlot
**/

#include "RooFit.h"

#include "RooList.h"
#include "RooMsgService.h"

#include "Riostream.h"

using namespace std;

ClassImp(RooList);



////////////////////////////////////////////////////////////////////////////////
/// Find the link corresponding to the named object in this list.
/// Return 0 if the object is not found or does not have an Option_t
/// string associated with its link. Also print a warning message
/// if caller is non-zero.

TObjOptLink *RooList::findLink(const char *name, const char *caller) const 
{
  if(0 == strlen(name)) return 0;
  TObjLink *link = FirstLink();
  while (link) {
    TObject *obj= link->GetObject();
    if (obj->GetName() && !strcmp(name, obj->GetName())) break;
    link = link->Next();
  }
  if(0 == link) {
    if(strlen(caller)) {
      coutE(InputArguments) << caller << ": cannot find object named \"" << name << "\"" << endl;
    }
    return 0;
  }
  return dynamic_cast<TObjOptLink*>(link);
}


////////////////////////////////////////////////////////////////////////////////
/// Move the target object immediately before the specified object,
/// preserving any Option_t associated with the target link.

Bool_t RooList::moveBefore(const char *before, const char *target, const char *caller) 
{
  // Find the target object's link
  TObjOptLink *targetLink= findLink(target,caller);
  if(0 == targetLink) return kFALSE;

  // Find the insert-before object's link
  TObjOptLink *beforeLink= findLink(before,caller);
  if(0 == beforeLink) return kFALSE;

  // Remember the target link's object and options
  TObject *obj= targetLink->GetObject();
  TString opt= targetLink->GetOption();

  // Remove the target object in its present position
  Remove(targetLink);

  // Add it back in its new position
  if(beforeLink == fFirst) {
    RooList::AddFirst(obj, opt.Data());
  }
  else {
    // coverity[RESOURCE_LEAK]
    NewOptLink(obj, opt.Data(), beforeLink->Prev());
    fSize++;
    Changed();
  }
  return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Move the target object immediately after the specified object,
/// preserving any Option_t associated with the target link.

Bool_t RooList::moveAfter(const char *after, const char *target, const char *caller) 
{
  // Find the target object's link
  TObjOptLink *targetLink= findLink(target,caller);
  if(0 == targetLink) return kFALSE;

  // Find the insert-after object's link
  TObjOptLink *afterLink= findLink(after,caller);
  if(0 == afterLink) return kFALSE;

  // Remember the target link's object and options
  TObject *obj= targetLink->GetObject();
  TString opt= targetLink->GetOption();

  // Remove the target object in its present position
  Remove(targetLink);

  // Add it back in its new position
  if(afterLink == fLast) {
    RooList::AddLast(obj, opt.Data());
  }
  else {
    NewOptLink(obj, opt.Data(), afterLink);
    fSize++;
    Changed();
  }
  return kTRUE;
}
