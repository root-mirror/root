// @(#)root/tree:$Id$
// Author: Rene Brun   27/01/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TLeafObject
\ingroup tree

A TLeaf for a general object derived from TObject.
*/

#include "TLeafObject.h"
#include "TBranch.h"
#include "TBuffer.h"
#include "TClass.h"
#include "TMethodCall.h"
#include "TDataType.h"

ClassImp(TLeafObject);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for LeafObject.

TLeafObject::TLeafObject(): TLeaf()
{
   fClass      = 0;
   fObjAddress = 0;
   fVirtual    = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a LeafObject.

TLeafObject::TLeafObject(TBranch *parent, const char *name, const char *type)
   :TLeaf(parent, name,type)
{
   SetTitle(type);
   fClass      = TClass::GetClass(type);
   fObjAddress = 0;
   fVirtual    = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor for a LeafObject.

TLeafObject::~TLeafObject()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Pack leaf elements in Basket output buffer.

void TLeafObject::FillBasket(TBuffer &b)
{
   if (!fObjAddress) return;
   TObject *object  = GetObject();
   if (object) {
      if (fVirtual) {
         UChar_t n = (UChar_t) strlen(object->ClassName());
         b << n;
         b.WriteFastArray(object->ClassName(),n+1);
      }
      object->Streamer(b);
   } else {
      if (fClass) {
         if (fClass->Property() & kIsAbstract) object = new TObject;
         else                                  object = (TObject *)fClass->New();
         object->SetBit(kInvalidObject);
         object->SetUniqueID(123456789);
         object->Streamer(b);
         if (fClass->Property() & kIsAbstract) delete object;
         else                                  fClass->Destructor(object);
      } else {
         Error("FillBasket","Attempt to write a NULL object in leaf:%s",GetName());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns pointer to method corresponding to name.
///
/// name is a string with the general form  "method(list of params)"
/// If list of params is omitted, () is assumed;

TMethodCall *TLeafObject::GetMethodCall(const char *name)
{
   char *namecpy = new char[strlen(name)+1];
   strcpy(namecpy,name);
   char *params = strchr(namecpy,'(');
   if (params) { *params = 0; params++; }
   else params = (char *) ")";

   if (!fClass) fClass      = TClass::GetClass(GetTitle());
   TMethodCall *m = new TMethodCall(fClass, namecpy, params);
   delete [] namecpy;
   if (m->GetMethod()) return m;
   Error("GetMethodCall","Unknown method:%s",name);
   delete m;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns name of leaf type.

const char *TLeafObject::GetTypeName() const
{
   return fTitle.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// This method must be overridden to handle object notification.

Bool_t TLeafObject::Notify()
{
   fClass      = TClass::GetClass(GetTitle());
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Prints leaf value.

void TLeafObject::PrintValue(Int_t) const
{
   printf("%zx\n",(size_t)GetValuePointer());
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer.

void TLeafObject::ReadBasket(TBuffer &b)
{
   char classname[128];
   UChar_t n;
   if (fVirtual) {
      b >> n;
      b.ReadFastArray(classname,n+1);
      fClass      = TClass::GetClass(classname);
   }
   if (fClass) {
      TObject *object;
      if (!fObjAddress) {
         Long_t *voidobj = new Long_t[1];
         fObjAddress  = (void **)voidobj;
         *fObjAddress = (TObject *)fClass->New();
      }
      object = (TObject*)(*fObjAddress);
      if (fBranch->IsAutoDelete()) {
         fClass->Destructor(object);
         object = (TObject *)fClass->New();
      }
      if (!object) return;

      if (fClass->GetState() > TClass::kEmulated) {
         object->Streamer(b);
      } else {
         //emulated class has no Streamer
         if (!TestBit(kWarn)) {
            Warning("ReadBasket","%s::Streamer not available, using TClass::ReadBuffer instead",fClass->GetName());
            SetBit(kWarn);
         }
         fClass->ReadBuffer(b,object);
      }
      // in case we had written a null pointer a Zombie object was created
      // we must delete it
      if (object->TestBit(kInvalidObject)) {
         if (object->GetUniqueID() == 123456789) {
            fClass->Destructor(object);
            object = 0;
         }
      }
      *fObjAddress = object;
   } else GetBranch()->SetAddress(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set leaf buffer data address.

void TLeafObject::SetAddress(void *add)
{
   fObjAddress = (void **)add;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TLeafObject.

void TLeafObject::Streamer(TBuffer &b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      if (R__v > 3 || R__v == 2) {
         b.ReadClassBuffer(TLeafObject::Class(), this, R__v, R__s, R__c);
         if (R__v == 2) fVirtual = kTRUE;
         fObjAddress = 0;
         fClass  = TClass::GetClass(fTitle.Data());
         if (!fClass) Warning("Streamer","Cannot find class:%s",fTitle.Data());

         // We should rewarn in this process.
         ResetBit(kWarn);
         ResetBit(kOldWarn);

         return;
      }
      //====process old versions before automatic schema evolution
      TLeaf::Streamer(b);
      fObjAddress = 0;
      fClass  = TClass::GetClass(fTitle.Data());
      if (!fClass) Warning("Streamer","Cannot find class:%s",fTitle.Data());
      if (R__v  < 1) fVirtual = kFALSE;
      if (R__v == 1) fVirtual = kTRUE;
      if (R__v == 3) b >> fVirtual;
      // We should rewarn in this process.
      ResetBit(kOldWarn);
      //====end of old versions

   } else {
      b.WriteClassBuffer(TLeafObject::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if this leaf is does not have any sub-branch/leaf.

Bool_t TLeafObject::IsOnTerminalBranch() const
{
   if (fBranch->GetListOfBranches()->GetEntriesFast()) return kFALSE;
   return kTRUE;
}
