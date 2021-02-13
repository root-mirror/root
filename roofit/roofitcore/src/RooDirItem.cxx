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
\file RooDirItem.cxx
\class RooDirItem
\ingroup Roofitcore

RooDirItem is a utility base class for RooFit objects that are to be attached
to ROOT directories. Concrete classes inherit the appendToDir and removeToDir
methods that can be used to safely attach and detach one self from a TDirectory
**/

#include "RooFit.h"

#include <iostream>
#include "TList.h"
#include "TDirectoryFile.h"
#include "TString.h"
#include "RooDirItem.h"

using namespace std;

ClassImp(RooDirItem);


////////////////////////////////////////////////////////////////////////////////
/// Remove object from directory it was added to

void RooDirItem::removeFromDir(TObject* obj) 
{
  if (_dir) {
    _dir->GetList()->Remove(obj) ;
  }
} 



////////////////////////////////////////////////////////////////////////////////
/// Append object to directory. If forceMemoryResident is
/// true, nothing happens.
void RooDirItem::appendToDir(TObject* obj, Bool_t forceMemoryResident) 
{
  if (forceMemoryResident) {
    // If we are not going into a file, appending to a directory
    // doesn't make sense. It only creates global state and congestion.
    return;
  } else {
    // Append self to present gDirectory
    _dir = gDirectory ;
    gDirectory->Append(obj) ;
  }
}

