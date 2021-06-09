// @(#)root/cont
// Author: Bianca-Cristina Cristescu February 2014

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TListOfEnums
A collection of TEnum objects designed for fast access given a
DeclId_t and for keep track of TEnum that were described
unloaded enum.
*/

#include <forward_list>

#include "TListOfEnums.h"
#include "TClass.h"
#include "TExMap.h"
#include "TEnum.h"
#include "TGlobal.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"

const unsigned int listSize=3;

ClassImp(TListOfEnums);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TListOfEnums::TListOfEnums(TClass *cl /*=0*/) :
   THashList(listSize), fClass(cl), fIds(0), fUnloaded(0), fIsLoaded(kFALSE), fLastLoadMarker(0)
{
   fIds = new TExMap(listSize);
   fUnloaded = new THashList(listSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TListOfEnums::~TListOfEnums()
{
   THashList::Delete();
   delete fIds;
   fUnloaded->Delete();
   delete fUnloaded;
}

////////////////////////////////////////////////////////////////////////////////
/// Add pair<id, object> to the map of functions and their ids.

void TListOfEnums::MapObject(TObject *obj)
{
   TEnum *e = dynamic_cast<TEnum *>(obj);
   if (e && e->GetDeclId()) {
      fIds->Add((Long64_t)e->GetDeclId(), (Long64_t)e);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add object at the beginning of the list.

void TListOfEnums::AddFirst(TObject *obj)
{
   THashList::AddFirst(obj);
   MapObject(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Add object at the beginning of the list and also store option.
/// Storing an option is useful when one wants to change the behaviour
/// of an object a little without having to create a complete new
/// copy of the object. This feature is used, for example, by the Draw()
/// method. It allows the same object to be drawn in different ways.

void TListOfEnums::AddFirst(TObject *obj, Option_t *opt)
{
   THashList::AddFirst(obj, opt);
   MapObject(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Add object at the end of the list.

void TListOfEnums::AddLast(TObject *obj)
{
   THashList::AddLast(obj);
   MapObject(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Add object at the end of the list and also store option.
/// Storing an option is useful when one wants to change the behaviour
/// of an object a little without having to create a complete new
/// copy of the object. This feature is used, for example, by the Draw()
/// method. It allows the same object to be drawn in different ways.

void TListOfEnums::AddLast(TObject *obj, Option_t *opt)
{
   THashList::AddLast(obj, opt);
   MapObject(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Insert object at location idx in the list.

void TListOfEnums::AddAt(TObject *obj, Int_t idx)
{
   THashList::AddAt(obj, idx);
   MapObject(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Insert object after object after in the list.

void TListOfEnums::AddAfter(const TObject *after, TObject *obj)
{
   THashList::AddAfter(after, obj);
   MapObject(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Insert object after object after in the list.

void TListOfEnums::AddAfter(TObjLink *after, TObject *obj)
{
   THashList::AddAfter(after, obj);
   MapObject(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Insert object before object before in the list.

void TListOfEnums::AddBefore(const TObject *before, TObject *obj)
{
   THashList::AddBefore(before, obj);
   MapObject(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Insert object before object before in the list.

void TListOfEnums::AddBefore(TObjLink *before, TObject *obj)
{
   THashList::AddBefore(before, obj);
   MapObject(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all objects from the list. Does not delete the objects unless
/// the THashList is the owner (set via SetOwner()).

void TListOfEnums::Clear(Option_t *option)
{
   fUnloaded->Clear(option);
   fIds->Clear();
   THashList::Clear(option);
   fIsLoaded = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete all TDataMember object files.

void TListOfEnums::Delete(Option_t *option /* ="" */)
{
   fUnloaded->Delete(option);
   THashList::Delete(option);
   fIsLoaded = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the TEnum corresponding to the Decl 'id' or NULL if it does not
/// exist.

TEnum *TListOfEnums::Find(DeclId_t id) const
{
   if (!id) return 0;

   return (TEnum *)fIds->GetValue((Long64_t)id);
}

////////////////////////////////////////////////////////////////////////////////
/// Return (after creating it if necessary) the TEnum
/// describing the enum corresponding to the Decl 'id'.

TEnum *TListOfEnums::Get(DeclId_t id, const char *name)
{
   if (!id) return 0;

   TEnum *e = Find(id);
   if (e) return e;

   // If this declID is not found as key, we look for the enum by name.
   // Indeed it could have been generated by protoclasses.
#if defined(R__MUST_REVISIT)
# if R__MUST_REVISIT(6,4)
   "This special case can be removed once PCMs are available."
# endif
#endif
   e = static_cast<TEnum*>(THashList::FindObject(name));
   if (e) {
      // In this case, we update the declId, update its constants and add the enum to the ids map and return.
      // At this point it is like it came from the interpreter.
      if (0 == e->GetDeclId()){
         e->Update(id);
         fIds->Add((Long64_t)id, (Long64_t)e);
         gInterpreter->UpdateEnumConstants(e, fClass);
      }
      return e;
   }

   if (fClass) {
      if (!fClass->HasInterpreterInfoInMemory()) {
         // The interpreter does not know about this class yet (or a problem
         // occurred that prevented the proper updating of fClassInfo).
         // So this decl can not possibly be part of this class.
         // [In addition calling GetClassInfo would trigger a late parsing
         //  of the header which we want to avoid].
         return 0;
      }
      if (!gInterpreter->ClassInfo_Contains(fClass->GetClassInfo(), id)) return 0;
   } else {
      if (!gInterpreter->ClassInfo_Contains(0, id)) return 0;
   }

   R__LOCKGUARD(gInterpreterMutex);

   // Let's see if this is a reload ...
   // can we check for reloads for enums?
   e = (TEnum *)fUnloaded->FindObject(name);
   if (e) {
      e->Update(id);
      gInterpreter->UpdateEnumConstants(e, fClass);
   } else {
      e = gInterpreter->CreateEnum((void *)id, fClass);
   }
   // Calling 'just' THahList::Add would turn around and call
   // TListOfEnums::AddLast which should *also* do the fIds->Add.
   THashList::AddLast(e);
   fIds->Add((Long64_t)id, (Long64_t)e);

   return e;
}

////////////////////////////////////////////////////////////////////////////////
/// Return an object from the list of enums *if and only if* is has already
/// been loaded in the list.  This is an internal routine.

TEnum *TListOfEnums::GetObject(const char *name) const
{
   return (TEnum*)THashList::FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a pair<id, object> from the map of functions and their ids.

void TListOfEnums::UnmapObject(TObject *obj)
{
   TEnum *e = dynamic_cast<TEnum *>(obj);
   if (e && e->GetDeclId()) {
      fIds->Remove((Long64_t)e->GetDeclId());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object from this collection and recursively remove the object
/// from all other objects (and collections).
/// This function overrides TCollection::RecursiveRemove that calls
/// the Remove function. THashList::Remove cannot be called because
/// it uses the hash value of the hash table. This hash value
/// is not available anymore when RecursiveRemove is called from
/// the TObject destructor.

void TListOfEnums::RecursiveRemove(TObject *obj)
{
   if (!obj) return;

   THashList::RecursiveRemove(obj);
   fUnloaded->RecursiveRemove(obj);
   UnmapObject(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object from the list.

TObject *TListOfEnums::Remove(TObject *obj)
{
   Bool_t found;

   found = THashList::Remove(obj);
   if (!found) {
      found = fUnloaded->Remove(obj);
   }
   UnmapObject(obj);
   if (found) return obj;
   else return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object via its objlink from the list.

TObject *TListOfEnums::Remove(TObjLink *lnk)
{
   if (!lnk) return 0;

   TObject *obj = lnk->GetObject();

   THashList::Remove(lnk);
   fUnloaded->Remove(obj);
   UnmapObject(obj);
   return obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Load all the DataMembers known to the interpreter for the scope 'fClass'
/// into this collection.

void TListOfEnums::Load()
{
   if (fClass && fClass->Property() & (kIsClass | kIsStruct | kIsUnion)) {
      // Class and union are not extendable, if we already
      // loaded all the data member there is no need to recheck
      if (fIsLoaded) return;
   }

   // This will provoke the parsing of the headers if need be.
   if (fClass && fClass->GetClassInfo() == 0) return;

   R__LOCKGUARD(gInterpreterMutex);

   ULong64_t currentTransaction = gInterpreter->GetInterpreterStateMarker();
   if (currentTransaction == fLastLoadMarker) {
      return;
   }
   fLastLoadMarker = currentTransaction;

   // In the case of namespace, even if we have loaded before we need to
   // load again in case there was new data member added.

   // Mark the list as loaded to avoid an infinite recursion in the case
   // where we have a data member that is a variable size array.  In that
   // case TDataMember::Init needs to get/load the list to find the data
   // member used as the array size.
   fIsLoaded = kTRUE;

   // Respawn the unloaded enums if they come from protoclasses, i.e. they
   // have a 0 declId.
#if defined(R__MUST_REVISIT)
# if R__MUST_REVISIT(6,4)
   "This special case can be removed once PCMs are available."
# endif
#endif

   std::forward_list<TEnum*> respownedEnums;
   for (auto enumAsObj : *fUnloaded){
      TEnum* en = static_cast<TEnum*>(enumAsObj);
      if (0 == en->GetDeclId()){
         THashList::AddLast(en);
         respownedEnums.push_front(en);
      }
   }

   for (auto en : respownedEnums)
      fUnloaded->Remove(en);

   // We cannot clear the whole unloaded list. It is too much.
//   fUnloaded->Clear();

   gInterpreter->LoadEnums(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Mark 'all func' as being unloaded.
/// After the unload, the data member can no longer be found directly,
/// until the decl can be found again in the interpreter (in which
/// the func object will be reused.

void TListOfEnums::Unload()
{
   TObjLink *lnk = FirstLink();
   while (lnk) {
      TEnum *data = (TEnum *)lnk->GetObject();

      if (data->GetDeclId())
         fIds->Remove((Long64_t)data->GetDeclId());
      fUnloaded->Add(data);

      lnk = lnk->Next();
   }

   THashList::Clear();
   fIsLoaded = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Mark enum 'e' as being unloaded.
/// After the unload, the data member can no longer be found directly,
/// until the decl can be found again in the interpreter (in which
/// the func object will be reused.

void TListOfEnums::Unload(TEnum *e)
{
   if (THashList::Remove(e)) {
      // We contains the object, let remove it from the other internal
      // list and move it to the list of unloaded objects.
      if (e->GetDeclId())
         fIds->Remove((Long64_t)e->GetDeclId());
      fUnloaded->Add(e);
   }
}
