// Author:  Sergey Linev, GSI  6/04/2017

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TWebSnapshot.h"

#include "TString.h"

#include <ROOT/RMakeUnique.hxx>

/** \class TWebSnapshot
\ingroup webgui6

Paint state of object to transfer to JavaScript side

*/

///////////////////////////////////////////////////////////////////////////////////////////
/// destructor

TWebSnapshot::~TWebSnapshot()
{
   SetSnapshot(kNone, nullptr);
}

///////////////////////////////////////////////////////////////////////////////////////////
/// SetUse pointer to assign object id - TString::Hash

void TWebSnapshot::SetSnapshot(Int_t kind, TObject *snapshot, Bool_t owner)
{
   if (fSnapshot && fOwner) delete fSnapshot;
   fKind = kind;
   fSnapshot = snapshot;
   fOwner = owner;
}

///////////////////////////////////////////////////////////////////////////////////////////
/// Use pointer to assign object id - TString::Hash

void TWebSnapshot::SetObjectIDAsPtr(void *ptr)
{
   UInt_t hash = TString::Hash(&ptr, sizeof(ptr));
   SetObjectID(std::to_string(hash));
}

///////////////////////////////////////////////////////////////////////////////////////////
/// Create new entry in list of primitives

TWebSnapshot &TPadWebSnapshot::NewPrimitive(TObject *obj, const std::string &opt)
{
   fPrimitives.emplace_back(std::make_unique<TWebSnapshot>());
   if (obj) {
      fPrimitives.back()->SetObjectIDAsPtr(obj);
      fPrimitives.back()->SetOption(opt);
   }
   return *(fPrimitives.back());
}

///////////////////////////////////////////////////////////////////////////////////////////
/// Create new entry for subpad

TPadWebSnapshot &TPadWebSnapshot::NewSubPad()
{
   auto res = new TPadWebSnapshot(IsReadOnly());
   fPrimitives.emplace_back(res);
   return *res;
}

///////////////////////////////////////////////////////////////////////////////////////////
/// Create new entry in list of primitives in the front

TWebSnapshot &TPadWebSnapshot::NewSpecials()
{
   fPrimitives.emplace(fPrimitives.begin(), std::make_unique<TWebSnapshot>());
   return *(fPrimitives.front());
}
