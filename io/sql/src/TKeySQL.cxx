// @(#)root/sql:$Id$
// Author: Sergey Linev  20/11/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**
\class TKeySQL
\ingroup IO

TKeySQL represents meta-information about object, which was written to
SQL database. It keeps object id, which used to locate object data
from database tables.
*/

#include "TKeySQL.h"

#include "TROOT.h"
#include "TClass.h"

#include "TSQLResult.h"
#include "TBufferSQL2.h"
#include "TSQLStructure.h"
#include "TSQLFile.h"
#include <cstdlib>

ClassImp(TKeySQL);

////////////////////////////////////////////////////////////////////////////////
/// Creates TKeySQL and convert obj data to TSQLStructure via TBufferSQL2

TKeySQL::TKeySQL(TDirectory *mother, const TObject *obj, const char *name, const char *title)
   : TKey(mother)
{
   if (name)
      SetName(name);
   else if (obj) {
      SetName(obj->GetName());
      fClassName = obj->ClassName();
   } else
      SetName("Noname");

   if (title)
      SetTitle(title);

   StoreKeyObject((void *)obj, obj ? obj->IsA() : 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Creates TKeySQL and convert obj data to TSQLStructure via TBufferSQL2

TKeySQL::TKeySQL(TDirectory *mother, const void *obj, const TClass *cl, const char *name, const char *title)
   : TKey(mother)
{
   if (name && *name)
      SetName(name);
   else
      SetName(cl ? cl->GetName() : "Noname");

   if (title)
      SetTitle(title);

   StoreKeyObject(obj, cl);
}

////////////////////////////////////////////////////////////////////////////////
/// Create TKeySQL object, which corresponds to single entry in keys table

TKeySQL::TKeySQL(TDirectory *mother, Long64_t keyid, Long64_t objid, const char *name, const char *title,
                 const char *keydatetime, Int_t cycle, const char *classname)
   : TKey(mother), fKeyId(keyid), fObjId(objid)
{
   SetName(name);
   if (title)
      SetTitle(title);
   TDatime dt(keydatetime);
   fDatime = dt;
   fCycle = cycle;
   fClassName = classname;
}

////////////////////////////////////////////////////////////////////////////////
/// Compares keydata with provided and return kTRUE if key was modified
/// Used in TFile::StreamKeysForDirectory() method to verify data for that keys
/// should be updated

Bool_t TKeySQL::IsKeyModified(const char *keyname, const char *keytitle, const char *keydatime, Int_t cycle,
                              const char *classname)
{
   Int_t len1 = (GetName() == 0) ? 0 : strlen(GetName());
   Int_t len2 = (keyname == 0) ? 0 : strlen(keyname);
   if (len1 != len2)
      return kTRUE;
   if ((len1 > 0) && (strcmp(GetName(), keyname) != 0))
      return kTRUE;

   len1 = (GetTitle() == 0) ? 0 : strlen(GetTitle());
   len2 = (keytitle == 0) ? 0 : strlen(keytitle);
   if (len1 != len2)
      return kTRUE;
   if ((len1 > 0) && (strcmp(GetTitle(), keytitle) != 0))
      return kTRUE;

   const char *tm = GetDatime().AsSQLString();
   len1 = (tm == 0) ? 0 : strlen(tm);
   len2 = (keydatime == 0) ? 0 : strlen(keydatime);
   if (len1 != len2)
      return kTRUE;
   if ((len1 > 0) && (strcmp(tm, keydatime) != 0))
      return kTRUE;

   if (cycle != GetCycle())
      return kTRUE;

   len1 = (GetClassName() == 0) ? 0 : strlen(GetClassName());
   len2 = (classname == 0) ? 0 : strlen(classname);
   if (len1 != len2)
      return kTRUE;
   if ((len1 > 0) && (strcmp(GetClassName(), classname) != 0))
      return kTRUE;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Removes key from current directory
/// Note: TKeySQL object is not deleted. You still have to call "delete key"

void TKeySQL::Delete(Option_t * /*option*/)
{
   TSQLFile *f = (TSQLFile *)GetFile();

   if (f != 0)
      f->DeleteKeyFromDB(GetDBKeyId());

   fMotherDir->GetListOfKeys()->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// return sql id of parent directory

Long64_t TKeySQL::GetDBDirId() const
{
   return GetMotherDir() ? GetMotherDir()->GetSeekDir() : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Stores object, associated with key, into data tables

void TKeySQL::StoreKeyObject(const void *obj, const TClass *cl)
{
   TSQLFile *f = (TSQLFile *)GetFile();

   fCycle = GetMotherDir()->AppendKey(this);

   fKeyId = f->DefineNextKeyId();

   fObjId = f->StoreObjectInTables(fKeyId, obj, cl);

   if (cl)
      fClassName = cl->GetName();

   if (GetDBObjId() >= 0) {
      fDatime.Set();
      if (!f->WriteKeyData(this)) {
         // cannot add entry to keys table
         Error("StoreKeyObject", "Cannot write data to key tables");
         // delete everything relevant for that key
         f->DeleteKeyFromDB(GetDBKeyId());
         fObjId = -1;
      }
   }

   if (GetDBObjId() < 0)
      GetMotherDir()->GetListOfKeys()->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// To read an object from the file.
/// The object associated to this key is read from the file into memory.
/// Before invoking this function, obj has been created via the
/// default constructor.

Int_t TKeySQL::Read(TObject *tobj)
{
   if (tobj == 0)
      return 0;

   void *res = ReadKeyObject(tobj, 0);

   return res == 0 ? 0 : 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Read object derived from TObject class
/// If it is not TObject or in case of error, return 0

TObject *TKeySQL::ReadObj()
{
   TObject *tobj = (TObject *)ReadKeyObject(0, TObject::Class());

   if (tobj) {
      if (gROOT->GetForceStyle())
         tobj->UseCurrentStyle();
      if (tobj->IsA() == TDirectoryFile::Class()) {
         TDirectoryFile *dir = (TDirectoryFile *)tobj;
         dir->SetName(GetName());
         dir->SetTitle(GetTitle());
         dir->SetSeekDir(GetDBKeyId());
         dir->SetMother(fMotherDir);
         dir->ReadKeys();
         fMotherDir->Append(dir);
      }
   }

   return tobj;
}

////////////////////////////////////////////////////////////////////////////////
/// Read object derived from TObject class
/// If it is not TObject or in case of error, return 0

TObject *TKeySQL::ReadObjWithBuffer(char * /*bufferRead*/)
{
   TObject *tobj = (TObject *)ReadKeyObject(0, TObject::Class());

   if (tobj) {
      if (gROOT->GetForceStyle())
         tobj->UseCurrentStyle();
      if (tobj->IsA() == TDirectoryFile::Class()) {
         TDirectoryFile *dir = (TDirectoryFile *)tobj;
         dir->SetName(GetName());
         dir->SetTitle(GetTitle());
         dir->SetSeekDir(GetDBKeyId());
         dir->SetMother(fMotherDir);
         dir->ReadKeys();
         fMotherDir->Append(dir);
      }
   }

   return tobj;
}

////////////////////////////////////////////////////////////////////////////////
/// Read object of any type from SQL database

void *TKeySQL::ReadObjectAny(const TClass *expectedClass)
{
   void *res = ReadKeyObject(0, expectedClass);

   if (res && (expectedClass == TDirectoryFile::Class())) {
      TDirectoryFile *dir = (TDirectoryFile *)res;
      dir->SetName(GetName());
      dir->SetTitle(GetTitle());
      dir->SetSeekDir(GetDBKeyId());
      dir->SetMother(fMotherDir);
      dir->ReadKeys();
      fMotherDir->Append(dir);
   }

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Read object, associated with key, from database

void *TKeySQL::ReadKeyObject(void *obj, const TClass *expectedClass)
{
   TSQLFile *f = (TSQLFile *)GetFile();

   if ((GetDBKeyId() <= 0) || !f)
      return obj;

   TBufferSQL2 buffer(TBuffer::kRead, f);

   buffer.InitMap();

   TClass *cl = nullptr;

   void *res = buffer.SqlReadAny(GetDBKeyId(), GetDBObjId(), &cl, obj);

   if (!cl || !res)
      return nullptr;

   Int_t delta = 0;

   if (expectedClass != 0) {
      delta = cl->GetBaseClassOffset(expectedClass);
      if (delta < 0) {
         if (obj == 0)
            cl->Destructor(res);
         return 0;
      }
      if (cl->GetState() > TClass::kEmulated && expectedClass->GetState() <= TClass::kEmulated) {
         // we cannot mix a compiled class with an emulated class in the inheritance
         Warning("XmlReadAny", "Trying to read an emulated class (%s) to store in a compiled pointer (%s)",
                 cl->GetName(), expectedClass->GetName());
      }
   }

   return ((char *)res) + delta;
}
