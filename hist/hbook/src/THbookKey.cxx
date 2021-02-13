// @(#)root/hbook:$Id$
// Author: Rene Brun   20/02/2002

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "THbookKey.h"
#include "THbookTree.h"
#include "TBrowser.h"
#include "snprintf.h"

ClassImp(THbookKey);

/** \class THbookKey
    \ingroup Hist
    \brief HBOOK Key
*/


////////////////////////////////////////////////////////////////////////////////
///constructor

THbookKey::THbookKey(Int_t id, THbookFile *file)
{
   fDirectory = file;
   fID = id;
   char name[10];
   snprintf(name,10,"h%d",id);
   SetName(name);
}


////////////////////////////////////////////////////////////////////////////////

THbookKey::~THbookKey()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Read object from disk and call its Browse() method.
/// If object with same name already exist in memory delete it (like
/// TDirectory::Get() is doing), except when the key references a
/// folder in which case we don't want to re-read the folder object
/// since it might contain new objects not yet saved.

void THbookKey::Browse(TBrowser *b)
{
   fDirectory->cd();

   TObject *obj = fDirectory->GetList()->FindObject(GetName());
   if (obj && !obj->IsFolder()) {
      if (obj->InheritsFrom(TCollection::Class()))
         obj->Delete();   // delete also collection elements
      delete obj;
      obj = 0;
   }

   if (!obj)
      obj = fDirectory->Get(fID);

   if (b && obj) {
      obj->Browse(b);
      b->SetRefreshFlag(kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
///an hbook key is not a folder

Bool_t THbookKey::IsFolder() const
{
   Bool_t ret = kFALSE;


   return( ret );
}
