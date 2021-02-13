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
\class TBufferSQL2
\ingroup IO

Converts data to SQL statements or read data from SQL tables.

Class for serializing/deserializing object to/from SQL data base.
It redefines most of TBuffer class function to convert simple types,
array of simple types and objects to/from TSQLStructure objects.
TBufferSQL2 class uses streaming mechanism, provided by ROOT system,
therefore most of ROOT and user classes can be stored. There are
limitations for complex objects like TTree, TClonesArray, TDirectory and
few other, which can not be converted to SQL (yet).
*/

#include "TBufferSQL2.h"

#include "TROOT.h"
#include "TDataType.h"
#include "TClass.h"
#include "TClassTable.h"
#include "TMap.h"
#include "TStreamerInfo.h"
#include "TStreamerElement.h"
#include "TMemberStreamer.h"
#include "TStreamer.h"
#include "TStreamerInfoActions.h"
#include "snprintf.h"

#include <iostream>
#include <cstdlib>
#include <string>

#include "TSQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"
#include "TSQLStructure.h"
#include "TSQLObjectData.h"
#include "TSQLFile.h"
#include "TSQLClassInfo.h"

ClassImp(TBufferSQL2);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor, should not be used

TBufferSQL2::TBufferSQL2()
   : TBufferText(), fSQL(nullptr), fIOVersion(1), fStructure(nullptr), fStk(0), fReadBuffer(), fErrorFlag(0),
     fCompressLevel(ROOT::RCompressionSetting::EAlgorithm::kUseGlobal), fReadVersionBuffer(-1), fObjIdCounter(1), fIgnoreVerification(kFALSE),
     fCurrentData(nullptr), fObjectsInfos(nullptr), fFirstObjId(0), fLastObjId(0), fPoolsMap(nullptr)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Creates buffer object to serialize/deserialize data to/from sql.
/// This constructor should be used, if data from buffer supposed to be stored in file.
/// Mode should be either TBuffer::kRead or TBuffer::kWrite.

TBufferSQL2::TBufferSQL2(TBuffer::EMode mode, TSQLFile *file)
   : TBufferText(mode, file), fSQL(nullptr), fIOVersion(1), fStructure(nullptr), fStk(0), fReadBuffer(), fErrorFlag(0),
     fCompressLevel(ROOT::RCompressionSetting::EAlgorithm::kUseGlobal), fReadVersionBuffer(-1), fObjIdCounter(1), fIgnoreVerification(kFALSE),
     fCurrentData(nullptr), fObjectsInfos(nullptr), fFirstObjId(0), fLastObjId(0), fPoolsMap(nullptr)
{
   fSQL = file;
   if (file) {
      SetCompressionLevel(file->GetCompressionLevel());
      fIOVersion = file->GetIOVersion();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy sql buffer.

TBufferSQL2::~TBufferSQL2()
{
   if (fStructure)
      delete fStructure;

   if (fObjectsInfos) {
      fObjectsInfos->Delete();
      delete fObjectsInfos;
   }

   if (fPoolsMap) {
      fPoolsMap->DeleteValues();
      delete fPoolsMap;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Convert object of any class to sql structures
/// Return pointer on created TSQLStructure
/// TSQLStructure object will be owned by TBufferSQL2

TSQLStructure *TBufferSQL2::SqlWriteAny(const void *obj, const TClass *cl, Long64_t objid)
{
   fErrorFlag = 0;

   fStructure = nullptr;

   fFirstObjId = objid;
   fObjIdCounter = objid;

   SqlWriteObject(obj, cl, kTRUE);

   if (gDebug > 3)
      if (fStructure) {
         std::cout << "==== Printout of Sql structures ===== " << std::endl;
         fStructure->Print("*");
         std::cout << "=========== End printout ============ " << std::endl;
      }

   return fStructure;
}

////////////////////////////////////////////////////////////////////////////////
/// Recreate object from sql structure.
/// Return pointer to read object.
/// if (cl!=0) returns pointer to class of object

void *TBufferSQL2::SqlReadAny(Long64_t keyid, Long64_t objid, TClass **cl, void *obj)
{
   if (cl)
      *cl = nullptr;
   if (!fSQL)
      return nullptr;

   fCurrentData = 0;
   fErrorFlag = 0;

   fReadVersionBuffer = -1;

   fObjectsInfos = fSQL->SQLObjectsInfo(keyid);
   fFirstObjId = objid;
   fLastObjId = objid;
   if (fObjectsInfos) {
      TSQLObjectInfo *objinfo = (TSQLObjectInfo *)fObjectsInfos->Last();
      if (objinfo)
         fLastObjId = objinfo->GetObjId();
   }

   return SqlReadObjectDirect(obj, cl, objid);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns object info like classname and version
/// Should be taken from buffer, which is produced in the beginning

Bool_t TBufferSQL2::SqlObjectInfo(Long64_t objid, TString &clname, Version_t &version)
{
   if ((objid < 0) || !fObjectsInfos)
      return kFALSE;

   // suppose that objects info are sorted out

   Long64_t shift = objid - fFirstObjId;

   TSQLObjectInfo *info = nullptr;
   if ((shift >= 0) && (shift <= fObjectsInfos->GetLast())) {
      info = (TSQLObjectInfo *)fObjectsInfos->At(shift);
      if (info->GetObjId() != objid)
         info = nullptr;
   }

   if (!info) {
      // I hope, i will never get inside it
      Info("SqlObjectInfo", "Standard not works %lld", objid);
      for (Int_t n = 0; n <= fObjectsInfos->GetLast(); n++) {
         info = (TSQLObjectInfo *)fObjectsInfos->At(n);
         if (info->GetObjId() == objid)
            break;
         info = nullptr;
      }
   }

   if (!info)
      return kFALSE;

   clname = info->GetObjClassName();
   version = info->GetObjVersion();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates TSQLObjectData for specified object id and specified class
///
/// Object data for each class can be stored in two different tables.
/// First table contains data in column-wise form for simple types like integer,
/// strings and so on when second table contains any other data which cannot
/// be converted into column-wise representation.
/// TSQLObjectData will contain results of the requests to both such tables for
/// concrete object id.

TSQLObjectData *TBufferSQL2::SqlObjectData(Long64_t objid, TSQLClassInfo *sqlinfo)
{
   TSQLResult *classdata = nullptr;
   TSQLRow *classrow = nullptr;

   if (sqlinfo->IsClassTableExist()) {

      TSQLObjectDataPool *pool = nullptr;

      if (fPoolsMap)
         pool = (TSQLObjectDataPool *)fPoolsMap->GetValue(sqlinfo);

      if (pool && (fLastObjId >= fFirstObjId)) {
         if (gDebug > 4)
            Info("SqlObjectData", "Before request to %s", sqlinfo->GetClassTableName());
         TSQLResult *alldata = fSQL->GetNormalClassDataAll(fFirstObjId, fLastObjId, sqlinfo);
         if (gDebug > 4)
            Info("SqlObjectData", "After request res = 0x%lx", (Long_t)alldata);
         if (!alldata) {
            Error("SqlObjectData", "Cannot get data from table %s", sqlinfo->GetClassTableName());
            return nullptr;
         }

         if (!fPoolsMap)
            fPoolsMap = new TMap();
         pool = new TSQLObjectDataPool(sqlinfo, alldata);
         fPoolsMap->Add(sqlinfo, pool);
      }

      if (!pool)
         return nullptr;

      if (pool->GetSqlInfo() != sqlinfo) {
         Error("SqlObjectData", "Missmatch in pools map !!! CANNOT BE !!!");
         return nullptr;
      }

      classdata = pool->GetClassData();

      classrow = pool->GetObjectRow(objid);
      if (!classrow) {
         Error("SqlObjectData", "Can not find row for objid = %lld in table %s", objid, sqlinfo->GetClassTableName());
         return nullptr;
      }
   }

   TSQLResult *blobdata = nullptr;
   TSQLStatement *blobstmt = fSQL->GetBlobClassDataStmt(objid, sqlinfo);

   if (!blobstmt)
      blobdata = fSQL->GetBlobClassData(objid, sqlinfo);

   return new TSQLObjectData(sqlinfo, objid, classdata, classrow, blobdata, blobstmt);
}

////////////////////////////////////////////////////////////////////////////////
/// Write object to buffer.
/// If object was written before, only pointer will be stored
/// Return id of saved object

Int_t TBufferSQL2::SqlWriteObject(const void *obj, const TClass *cl, Bool_t cacheReuse, TMemberStreamer *streamer,
                                  Int_t streamer_index)
{
   if (gDebug > 1)
      Info("SqlWriteObject", "Object: %p Class: %s", obj, (cl ? cl->GetName() : "null"));

   PushStack();

   Long64_t objid = -1;

   if (!cl)
      obj = nullptr;

   if (!obj) {
      objid = 0;
   } else {
      Long64_t value = GetObjectTag(obj);
      if (value > 0)
         objid = fFirstObjId + value - 1;
   }

   if (gDebug > 1)
      Info("SqlWriteObject", "Find objectid %ld", (long)objid);

   if (objid >= 0) {
      Stack()->SetObjectPointer(objid);
      PopStack();
      return objid;
   }

   objid = fObjIdCounter++;

   Stack()->SetObjectRef(objid, cl);

   if (cacheReuse)
      MapObject(obj, cl, objid - fFirstObjId + 1);

   if (streamer)
      (*streamer)(*this, (void *)obj, streamer_index);
   else
      ((TClass *)cl)->Streamer((void *)obj, *this);

   if (gDebug > 1)
      Info("SqlWriteObject", "Done write of %s", cl->GetName());

   PopStack();

   return objid;
}

////////////////////////////////////////////////////////////////////////////////
/// Read object from the buffer

void *TBufferSQL2::SqlReadObject(void *obj, TClass **cl, TMemberStreamer *streamer, Int_t streamer_index,
                                 const TClass *onFileClass)
{
   if (cl)
      *cl = nullptr;

   if (fErrorFlag > 0)
      return obj;

   Bool_t findptr = kFALSE;

   const char *refid = fCurrentData->GetValue();
   if ((refid == 0) || (strlen(refid) == 0)) {
      Error("SqlReadObject", "Invalid object reference value");
      fErrorFlag = 1;
      return obj;
   }

   Long64_t objid = (Long64_t)std::stoll(refid);

   if (gDebug > 2)
      Info("SqlReadObject", "Starting objid: %ld column: %s", (long)objid, fCurrentData->GetLocatedField());

   if (!fCurrentData->IsBlobData() || fCurrentData->VerifyDataType(sqlio::ObjectPtr, kFALSE)) {
      if (objid == 0) {
         obj = nullptr;
         findptr = kTRUE;
      } else if (objid == -1) {
         findptr = kTRUE;
      } else if (objid >= fFirstObjId) {
         void *obj1 = nullptr;
         TClass *cl1 = nullptr;
         GetMappedObject(objid - fFirstObjId + 1, obj1, cl1);
         if (obj1 && cl1) {
            obj = obj1;
            if (cl)
               *cl = cl1;
         }
      }
   }

   if ((gDebug > 3) && findptr)
      Info("SqlReadObject", "Found pointer %p cl %s", obj, ((cl && *cl) ? (*cl)->GetName() : "null"));

   if (findptr) {
      fCurrentData->ShiftToNextValue();
      return obj;
   }

   if (fCurrentData->IsBlobData())
      if (!fCurrentData->VerifyDataType(sqlio::ObjectRef)) {
         Error("SqlReadObject", "Object reference or pointer is not found in blob data");
         fErrorFlag = 1;
         return obj;
      }

   fCurrentData->ShiftToNextValue();

   if ((gDebug > 2) || (objid < 0))
      Info("SqlReadObject", "Found object reference %ld", (long)objid);

   return SqlReadObjectDirect(obj, cl, objid, streamer, streamer_index, onFileClass);
}

////////////////////////////////////////////////////////////////////////////////
/// Read object data.
/// Class name and version are taken from special objects table.

void *TBufferSQL2::SqlReadObjectDirect(void *obj, TClass **cl, Long64_t objid, TMemberStreamer *streamer,
                                       Int_t streamer_index, const TClass *onFileClass)
{
   TString clname;
   Version_t version;

   if (!SqlObjectInfo(objid, clname, version))
      return obj;

   if (gDebug > 2)
      Info("SqlReadObjectDirect", "objid = %lld clname = %s ver = %d", objid, clname.Data(), version);

   TSQLClassInfo *sqlinfo = fSQL->FindSQLClassInfo(clname.Data(), version);

   TClass *objClass = TClass::GetClass(clname);
   if (objClass == TDirectory::Class())
      objClass = TDirectoryFile::Class();

   if (!objClass || !sqlinfo) {
      Error("SqlReadObjectDirect", "Class %s is not known", clname.Data());
      return obj;
   }

   if (!obj)
      obj = objClass->New();

   MapObject(obj, objClass, objid - fFirstObjId + 1);

   PushStack()->SetObjectRef(objid, objClass);

   TSQLObjectData *olddata = fCurrentData;

   if (sqlinfo->IsClassTableExist()) {
      // TObject and TString classes treated differently
      if ((objClass == TObject::Class()) || (objClass == TString::Class())) {

         TSQLObjectData *objdata = new TSQLObjectData;
         if (objClass == TObject::Class())
            TSQLStructure::UnpackTObject(fSQL, this, objdata, objid, version);
         else if (objClass == TString::Class())
            TSQLStructure::UnpackTString(fSQL, this, objdata, objid, version);

         Stack()->AddObjectData(objdata);
         fCurrentData = objdata;
      } else
         // before normal streamer first version will be read and
         // then streamer functions of TStreamerInfo class
         fReadVersionBuffer = version;
   } else {
      TSQLObjectData *objdata = SqlObjectData(objid, sqlinfo);
      if (!objdata || !objdata->PrepareForRawData()) {
         Error("SqlReadObjectDirect", "No found raw data for obj %lld in class %s version %d table", objid,
               clname.Data(), version);
         fErrorFlag = 1;
         return obj;
      }

      Stack()->AddObjectData(objdata);

      fCurrentData = objdata;
   }

   if (streamer) {
      streamer->SetOnFileClass(onFileClass);
      (*streamer)(*this, (void *)obj, streamer_index);
   } else {
      objClass->Streamer((void *)obj, *this, onFileClass);
   }

   PopStack();

   if (gDebug > 1)
      Info("SqlReadObjectDirect", "Read object of class %s done", objClass->GetName());

   if (cl)
      *cl = objClass;

   fCurrentData = olddata;

   return obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Function is called from TStreamerInfo WriteBuffer and Readbuffer functions
/// and indent new level in data structure.
/// This call indicates, that TStreamerInfo functions starts streaming
/// object data of correspondent class

void TBufferSQL2::IncrementLevel(TVirtualStreamerInfo *info)
{
   if (!info)
      return;

   PushStack()->SetStreamerInfo((TStreamerInfo *)info);

   if (gDebug > 2)
      Info("IncrementLevel", "Info: %s", info->GetName());

   WorkWithClass(info->GetName(), info->GetClassVersion());
}

////////////////////////////////////////////////////////////////////////////////
/// Function is called from TStreamerInfo WriteBuffer and Readbuffer functions
/// and decrease level in sql structure.

void TBufferSQL2::DecrementLevel(TVirtualStreamerInfo *info)
{
   if (Stack()->GetElement())
      PopStack(); // for element
   PopStack();    // for streamerinfo

   // restore value of object data
   fCurrentData = Stack()->GetObjectData(kTRUE);

   if (gDebug > 2)
      Info("DecrementLevel", "Info: %s", info->GetName());
}

////////////////////////////////////////////////////////////////////////////////
/// Function is called from TStreamerInfo WriteBuffer and Readbuffer functions
/// and add/verify next element in sql tables
/// This calls allows separate data, correspondent to one class member, from another

void TBufferSQL2::SetStreamerElementNumber(TStreamerElement *elem, Int_t comp_type)
{
   if (Stack()->GetElement())
      PopStack(); // was with if (number > 0), i.e. not first element.
   TSQLStructure *curr = Stack();

   TStreamerInfo *info = curr->GetStreamerInfo();
   if (!info) {
      Error("SetStreamerElementNumber", "Error in structures stack");
      return;
   }

   WorkWithElement(elem, comp_type);
}

////////////////////////////////////////////////////////////////////////////////
/// This method inform buffer data of which class now
/// will be streamed. When reading, classversion should be specified
/// as was read by TBuffer::ReadVersion().
///
/// ClassBegin(), ClassEnd() & ClassMemeber() should be used in
/// custom class streamers to specify which kind of data are
/// now streamed to/from buffer. That information is used to correctly
/// convert class data to/from "normal" sql tables with meaningfull names
/// and correct datatypes. Without that functions data from custom streamer
/// will be saved as "raw" data in special _streamer_ table one value after another
/// Such MUST be used when object is written with standard ROOT streaming
/// procedure, but should be read back in custom streamer.
/// For example, custom streamer of TNamed class may look like:

void TBufferSQL2::ClassBegin(const TClass *cl, Version_t classversion)
{
   // void TNamed::Streamer(TBuffer &b)
   //   UInt_t R__s, R__c;
   //   if (b.IsReading()) {
   //      Version_t R__v = b.ReadVersion(&R__s, &R__c);
   //      b.ClassBegin(TNamed::Class(), R__v);
   //      b.ClassMember("TObject");
   //      TObject::Streamer(b);
   //      b.ClassMember("fName","TString");
   //      fName.Streamer(b);
   //      b.ClassMember("fTitle","TString");
   //      fTitle.Streamer(b);
   //      b.ClassEnd(TNamed::Class());
   //      b.SetBufferOffset(R__s+R__c+sizeof(UInt_t));
   //   } else {
   //      TNamed::Class()->WriteBuffer(b,this);
   //   }

   if (classversion < 0)
      classversion = cl->GetClassVersion();

   PushStack()->SetCustomClass(cl, classversion);

   if (gDebug > 2)
      Info("ClassBegin", "Class: %s", cl->GetName());

   WorkWithClass(cl->GetName(), classversion);
}

////////////////////////////////////////////////////////////////////////////////
/// Method indicates end of streaming of classdata in custom streamer.
/// See ClassBegin() method for more details.

void TBufferSQL2::ClassEnd(const TClass *cl)
{
   if (Stack()->GetType() == TSQLStructure::kSqlCustomElement)
      PopStack(); // for element
   PopStack();    // for streamerinfo

   // restore value of object data
   fCurrentData = Stack()->GetObjectData(kTRUE);

   if (gDebug > 2)
      Info("ClassEnd", "Class: %s", cl->GetName());
}

////////////////////////////////////////////////////////////////////////////////
/// Method indicates name and typename of class memeber,
/// which should be now streamed in custom streamer
/// Following combinations are supported:
/// see TBufferXML::ClassMember for the details.

void TBufferSQL2::ClassMember(const char *name, const char *typeName, Int_t arrsize1, Int_t arrsize2)
{
   if (!typeName)
      typeName = name;

   if (!name || (strlen(name) == 0)) {
      Error("ClassMember", "Invalid member name");
      fErrorFlag = 1;
      return;
   }

   TString tname = typeName;

   Int_t typ_id = -1;

   if (strcmp(typeName, "raw:data") == 0)
      typ_id = TStreamerInfo::kMissing;

   if (typ_id < 0) {
      TDataType *dt = gROOT->GetType(typeName);
      if (dt)
         if ((dt->GetType() > 0) && (dt->GetType() < 20))
            typ_id = dt->GetType();
   }

   if (typ_id < 0)
      if (strcmp(name, typeName) == 0) {
         TClass *cl = TClass::GetClass(tname.Data());
         if (cl)
            typ_id = TStreamerInfo::kBase;
      }

   if (typ_id < 0) {
      Bool_t isptr = kFALSE;
      if (tname[tname.Length() - 1] == '*') {
         tname.Resize(tname.Length() - 1);
         isptr = kTRUE;
      }
      TClass *cl = TClass::GetClass(tname.Data());
      if (!cl) {
         Error("ClassMember", "Invalid class specifier %s", typeName);
         fErrorFlag = 1;
         return;
      }

      if (cl->IsTObject())
         typ_id = isptr ? TStreamerInfo::kObjectp : TStreamerInfo::kObject;
      else
         typ_id = isptr ? TStreamerInfo::kAnyp : TStreamerInfo::kAny;

      if ((cl == TString::Class()) && !isptr)
         typ_id = TStreamerInfo::kTString;
   }

   TStreamerElement *elem = nullptr;

   if (typ_id == TStreamerInfo::kMissing) {
      elem = new TStreamerElement(name, "title", 0, typ_id, "raw:data");
   } else if (typ_id == TStreamerInfo::kBase) {
      TClass *cl = TClass::GetClass(tname.Data());
      if (cl) {
         TStreamerBase *b = new TStreamerBase(tname.Data(), "title", 0);
         b->SetBaseVersion(cl->GetClassVersion());
         elem = b;
      }
   } else if ((typ_id > 0) && (typ_id < 20)) {
      elem = new TStreamerBasicType(name, "title", 0, typ_id, typeName);
   } else if ((typ_id == TStreamerInfo::kObject) || (typ_id == TStreamerInfo::kTObject) ||
              (typ_id == TStreamerInfo::kTNamed)) {
      elem = new TStreamerObject(name, "title", 0, tname.Data());
   } else if (typ_id == TStreamerInfo::kObjectp) {
      elem = new TStreamerObjectPointer(name, "title", 0, tname.Data());
   } else if (typ_id == TStreamerInfo::kAny) {
      elem = new TStreamerObjectAny(name, "title", 0, tname.Data());
   } else if (typ_id == TStreamerInfo::kAnyp) {
      elem = new TStreamerObjectAnyPointer(name, "title", 0, tname.Data());
   } else if (typ_id == TStreamerInfo::kTString) {
      elem = new TStreamerString(name, "title", 0);
   }

   if (!elem) {
      Error("ClassMember", "Invalid combination name = %s type = %s", name, typeName);
      fErrorFlag = 1;
      return;
   }

   if (arrsize1 > 0) {
      elem->SetArrayDim(arrsize2 > 0 ? 2 : 1);
      elem->SetMaxIndex(0, arrsize1);
      if (arrsize2 > 0)
         elem->SetMaxIndex(1, arrsize2);
   }

   // return stack to CustomClass node
   if (Stack()->GetType() == TSQLStructure::kSqlCustomElement)
      PopStack();

   // we indicate that there is no streamerinfo
   WorkWithElement(elem, -1);
}

////////////////////////////////////////////////////////////////////////////////
/// This function is a part of IncrementLevel method.
/// Also used in StartClass method

void TBufferSQL2::WorkWithClass(const char *classname, Version_t classversion)
{
   if (IsReading()) {
      Long64_t objid = 0;

      //      if ((fCurrentData!=0) && fCurrentData->VerifyDataType(sqlio::ObjectInst, kFALSE))
      //        if (!fCurrentData->IsBlobData()) Info("WorkWithClass","Big problem %s", fCurrentData->GetValue());

      if (fCurrentData && fCurrentData->IsBlobData() && fCurrentData->VerifyDataType(sqlio::ObjectInst, kFALSE)) {
         objid = atoi(fCurrentData->GetValue());
         fCurrentData->ShiftToNextValue();
         TString sobjid;
         sobjid.Form("%lld", objid);
         Stack()->ChangeValueOnly(sobjid.Data());
      } else
         objid = Stack()->DefineObjectId(kTRUE);
      if (objid < 0) {
         Error("WorkWithClass", "cannot define object id");
         fErrorFlag = 1;
         return;
      }

      TSQLClassInfo *sqlinfo = fSQL->FindSQLClassInfo(classname, classversion);
      if (!sqlinfo) {
         Error("WorkWithClass", "Can not find table for class %s version %d", classname, classversion);
         fErrorFlag = 1;
         return;
      }

      TSQLObjectData *objdata = SqlObjectData(objid, sqlinfo);
      if (!objdata) {
         Error("WorkWithClass", "Request error for data of object %lld for class %s version %d", objid, classname,
               classversion);
         fErrorFlag = 1;
         return;
      }

      Stack()->AddObjectData(objdata);

      fCurrentData = objdata;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// This function is a part of SetStreamerElementNumber method.
/// It is introduced for reading of data for specified data member of class.
/// Used also in ReadFastArray methods to resolve problem of compressed data,
/// when several data members of the same basic type streamed with single ...FastArray call

void TBufferSQL2::WorkWithElement(TStreamerElement *elem, Int_t /* comp_type */)
{
   if (gDebug > 2)
      Info("WorkWithElement", "elem = %s", elem->GetName());

   TSQLStructure *stack = Stack(1);
   TStreamerInfo *info = stack->GetStreamerInfo();
   Int_t number = info ? info->GetElements()->IndexOf(elem) : -1;

   if (number >= 0)
      PushStack()->SetStreamerElement(elem, number);
   else
      PushStack()->SetCustomElement(elem);

   if (IsReading()) {

      if (!fCurrentData) {
         Error("WorkWithElement", "Object data is lost");
         fErrorFlag = 1;
         return;
      }

      fCurrentData = Stack()->GetObjectData(kTRUE);

      Int_t located = Stack()->LocateElementColumn(fSQL, this, fCurrentData);

      if (located == TSQLStructure::kColUnknown) {
         Error("WorkWithElement", "Cannot locate correct column in the table");
         fErrorFlag = 1;
         return;
      } else if ((located == TSQLStructure::kColObject) || (located == TSQLStructure::kColObjectArray) ||
                 (located == TSQLStructure::kColParent)) {
         // search again for object data while for BLOB it should be already assign
         fCurrentData = Stack()->GetObjectData(kTRUE);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Suppressed function of TBuffer

TClass *TBufferSQL2::ReadClass(const TClass *, UInt_t *)
{
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Suppressed function of TBuffer

void TBufferSQL2::WriteClass(const TClass *)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Read version value from buffer
/// actually version is normally defined by table name
/// and kept in intermediate variable fReadVersionBuffer

Version_t TBufferSQL2::ReadVersion(UInt_t *start, UInt_t *bcnt, const TClass *)
{
   Version_t res = 0;

   if (start)
      *start = 0;
   if (bcnt)
      *bcnt = 0;

   if (fReadVersionBuffer >= 0) {
      res = fReadVersionBuffer;
      fReadVersionBuffer = -1;
      if (gDebug > 3)
         Info("ReadVersion", "from buffer = %d", (int)res);
   } else if (fCurrentData && fCurrentData->IsBlobData() && fCurrentData->VerifyDataType(sqlio::Version)) {
      TString value = fCurrentData->GetValue();
      res = value.Atoi();
      if (gDebug > 3)
         Info("ReadVersion", "from blob %s = %d", fCurrentData->GetBlobPrefixName(), (int)res);
      fCurrentData->ShiftToNextValue();
   } else {
      Error("ReadVersion", "No correspondent tags to read version");
      fErrorFlag = 1;
   }

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Copies class version to buffer, but not writes it to sql immidiately
/// Version will be used to produce complete table
/// name, which will include class version

UInt_t TBufferSQL2::WriteVersion(const TClass *cl, Bool_t /* useBcnt */)
{
   if (gDebug > 2)
      Info("WriteVersion", "cl:%s ver:%d", (cl ? cl->GetName() : "null"), (int)(cl ? cl->GetClassVersion() : 0));

   if (cl)
      Stack()->AddVersion(cl);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read object from buffer. Only used from TBuffer.

void *TBufferSQL2::ReadObjectAny(const TClass *)
{
   return SqlReadObject(0);
}

////////////////////////////////////////////////////////////////////////////////
/// ?????? Skip any kind of object from buffer
/// !!!!!! fix me, not yet implemented
/// Should be just skip of current column later

void TBufferSQL2::SkipObjectAny()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Write object to buffer. Only used from TBuffer

void TBufferSQL2::WriteObjectClass(const void *actualObjStart, const TClass *actualClass, Bool_t cacheReuse)
{
   if (gDebug > 2)
      Info("WriteObjectClass", "class %s", (actualClass ? actualClass->GetName() : " null"));
   SqlWriteObject(actualObjStart, actualClass, cacheReuse);
}

////////////////////////////////////////////////////////////////////////////////
/// Template method to read array content

template <typename T>
R__ALWAYS_INLINE void TBufferSQL2::SqlReadArrayContent(T *arr, Int_t arrsize, Bool_t withsize)
{
   if (gDebug > 3)
      Info("SqlReadArrayContent", "size %d", (int)(arrsize));
   PushStack()->SetArray(withsize ? arrsize : -1);
   Int_t indx(0), first, last;
   if (fCurrentData->IsBlobData()) {
      while (indx < arrsize) {
         const char *name = fCurrentData->GetBlobPrefixName();
         if (strstr(name, sqlio::IndexSepar) == 0) {
            sscanf(name, "[%d", &first);
            last = first;
         } else {
            sscanf(name, "[%d..%d", &first, &last);
         }
         if ((first != indx) || (last < first) || (last >= arrsize)) {
            Error("SqlReadArrayContent", "Error reading array content %s", name);
            fErrorFlag = 1;
            break;
         }
         SqlReadBasic(arr[indx++]);
         while (indx <= last)
            arr[indx++] = arr[first];
      }
   } else {
      while (indx < arrsize)
         SqlReadBasic(arr[indx++]);
   }
   PopStack();
   if (gDebug > 3)
      Info("SqlReadArrayContent", "done");
}

template <typename T>
R__ALWAYS_INLINE Int_t TBufferSQL2::SqlReadArray(T *&arr, Bool_t is_static)
{
   Int_t n = SqlReadArraySize();
   if (n <= 0)
      return 0;
   if (!arr) {
      if (is_static)
         return 0;
      arr = new T[n];
   }
   SqlReadArrayContent(arr, n, kTRUE);
   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Bool_t from buffer

Int_t TBufferSQL2::ReadArray(Bool_t *&b)
{
   return SqlReadArray(b);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Char_t from buffer

Int_t TBufferSQL2::ReadArray(Char_t *&c)
{
   return SqlReadArray(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UChar_t from buffer

Int_t TBufferSQL2::ReadArray(UChar_t *&c)
{
   return SqlReadArray(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Short_t from buffer

Int_t TBufferSQL2::ReadArray(Short_t *&h)
{
   return SqlReadArray(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UShort_t from buffer

Int_t TBufferSQL2::ReadArray(UShort_t *&h)
{
   return SqlReadArray(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Int_t from buffer

Int_t TBufferSQL2::ReadArray(Int_t *&i)
{
   return SqlReadArray(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UInt_t from buffer

Int_t TBufferSQL2::ReadArray(UInt_t *&i)
{
   return SqlReadArray(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long_t from buffer

Int_t TBufferSQL2::ReadArray(Long_t *&l)
{
   return SqlReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong_t from buffer

Int_t TBufferSQL2::ReadArray(ULong_t *&l)
{
   return SqlReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long64_t from buffer

Int_t TBufferSQL2::ReadArray(Long64_t *&l)
{
   return SqlReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong64_t from buffer

Int_t TBufferSQL2::ReadArray(ULong64_t *&l)
{
   return SqlReadArray(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float_t from buffer

Int_t TBufferSQL2::ReadArray(Float_t *&f)
{
   return SqlReadArray(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double_t from buffer

Int_t TBufferSQL2::ReadArray(Double_t *&d)
{
   return SqlReadArray(d);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Bool_t from buffer

Int_t TBufferSQL2::ReadStaticArray(Bool_t *b)
{
   return SqlReadArray(b, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Char_t from buffer

Int_t TBufferSQL2::ReadStaticArray(Char_t *c)
{
   return SqlReadArray(c, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UChar_t from buffer

Int_t TBufferSQL2::ReadStaticArray(UChar_t *c)
{
   return SqlReadArray(c, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Short_t from buffer

Int_t TBufferSQL2::ReadStaticArray(Short_t *h)
{
   return SqlReadArray(h, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UShort_t from buffer

Int_t TBufferSQL2::ReadStaticArray(UShort_t *h)
{
   return SqlReadArray(h, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Int_t from buffer

Int_t TBufferSQL2::ReadStaticArray(Int_t *i)
{
   return SqlReadArray(i, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UInt_t from buffer

Int_t TBufferSQL2::ReadStaticArray(UInt_t *i)
{
   return SqlReadArray(i, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long_t from buffer

Int_t TBufferSQL2::ReadStaticArray(Long_t *l)
{
   return SqlReadArray(l, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong_t from buffer

Int_t TBufferSQL2::ReadStaticArray(ULong_t *l)
{
   return SqlReadArray(l, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long64_t from buffer

Int_t TBufferSQL2::ReadStaticArray(Long64_t *l)
{
   return SqlReadArray(l, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong64_t from buffer

Int_t TBufferSQL2::ReadStaticArray(ULong64_t *l)
{
   return SqlReadArray(l, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float_t from buffer

Int_t TBufferSQL2::ReadStaticArray(Float_t *f)
{
   return SqlReadArray(f, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double_t from buffer

Int_t TBufferSQL2::ReadStaticArray(Double_t *d)
{
   return SqlReadArray(d, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Template method to read content of array, which not include size of array

template <typename T>
R__ALWAYS_INLINE void TBufferSQL2::SqlReadFastArray(T *arr, Int_t arrsize)
{
   if (arrsize > 0)
      SqlReadArrayContent(arr, arrsize, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Bool_t from buffer

void TBufferSQL2::ReadFastArray(Bool_t *b, Int_t n)
{
   SqlReadFastArray(b, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Char_t from buffer
/// if nodename==CharStar, read all array as string

void TBufferSQL2::ReadFastArray(Char_t *c, Int_t n)
{
   if ((n > 0) && fCurrentData->IsBlobData() && fCurrentData->VerifyDataType(sqlio::CharStar, kFALSE)) {
      const char *buf = SqlReadCharStarValue();
      if (!buf || (n <= 0))
         return;
      Int_t size = strlen(buf);
      if (size < n)
         size = n;
      memcpy(c, buf, size);
   } else {
      SqlReadFastArray(c, n);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UChar_t from buffer

void TBufferSQL2::ReadFastArray(UChar_t *c, Int_t n)
{
   SqlReadFastArray(c, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Short_t from buffer

void TBufferSQL2::ReadFastArray(Short_t *h, Int_t n)
{
   SqlReadFastArray(h, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UShort_t from buffer

void TBufferSQL2::ReadFastArray(UShort_t *h, Int_t n)
{
   SqlReadFastArray(h, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Int_t from buffer

void TBufferSQL2::ReadFastArray(Int_t *i, Int_t n)
{
   SqlReadFastArray(i, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of UInt_t from buffer

void TBufferSQL2::ReadFastArray(UInt_t *i, Int_t n)
{
   SqlReadFastArray(i, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long_t from buffer

void TBufferSQL2::ReadFastArray(Long_t *l, Int_t n)
{
   SqlReadFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong_t from buffer

void TBufferSQL2::ReadFastArray(ULong_t *l, Int_t n)
{
   SqlReadFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Long64_t from buffer

void TBufferSQL2::ReadFastArray(Long64_t *l, Int_t n)
{
   SqlReadFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ULong64_t from buffer

void TBufferSQL2::ReadFastArray(ULong64_t *l, Int_t n)
{
   SqlReadFastArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float_t from buffer

void TBufferSQL2::ReadFastArray(Float_t *f, Int_t n)
{
   SqlReadFastArray(f, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n characters from the I/O buffer.
/// Used only from TLeafC, dummy implementation here

void TBufferSQL2::ReadFastArrayString(Char_t *c, Int_t n)
{
   SqlReadFastArray(c, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double_t from buffer

void TBufferSQL2::ReadFastArray(Double_t *d, Int_t n)
{
   SqlReadFastArray(d, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Same functionality as TBuffer::ReadFastArray(...) but
/// instead of calling cl->Streamer(obj,buf) call here
/// buf.StreamObject(obj, cl). In that case it is easy to understand where
/// object data is started and finished

void TBufferSQL2::ReadFastArray(void *start, const TClass *cl, Int_t n, TMemberStreamer *streamer,
                                const TClass *onFileClass)
{
   if (gDebug > 2)
      Info("ReadFastArray", "(void *");

   if (streamer) {
      StreamObjectExtra(start, streamer, cl, 0, onFileClass);
      //      (*streamer)(*this,start,0);
      return;
   }

   int objectSize = cl->Size();
   char *obj = (char *)start;
   char *end = obj + n * objectSize;

   for (; obj < end; obj += objectSize) {
      StreamObject(obj, cl, onFileClass);
   }
   //   TBuffer::ReadFastArray(start, cl, n, s);
}

////////////////////////////////////////////////////////////////////////////////
/// Same functionality as TBuffer::ReadFastArray(...) but
/// instead of calling cl->Streamer(obj,buf) call here
/// buf.StreamObject(obj, cl). In that case it is easy to understand where
/// object data is started and finished

void TBufferSQL2::ReadFastArray(void **start, const TClass *cl, Int_t n, Bool_t isPreAlloc, TMemberStreamer *streamer,
                                const TClass *onFileClass)
{
   if (gDebug > 2)
      Info("ReadFastArray", "(void **  pre = %d  n = %d", isPreAlloc, n);

   Bool_t oldStyle = kFALSE; // flag used to reproduce old-style I/O actions for kSTLp

   if ((fIOVersion < 2) && !isPreAlloc) {
      TStreamerElement *elem = Stack(0)->GetElement();
      if (elem && ((elem->GetType() == TStreamerInfo::kSTLp) ||
                   (elem->GetType() == TStreamerInfo::kSTLp + TStreamerInfo::kOffsetL)))
         oldStyle = kTRUE;
   }

   if (streamer) {
      if (isPreAlloc) {
         for (Int_t j = 0; j < n; j++) {
            if (!start[j])
               start[j] = ((TClass *)cl)->New();
         }
      }
      if (oldStyle)
         (*streamer)(*this, (void *)start, n);
      else
         StreamObjectExtra((void *)start, streamer, cl, 0, onFileClass);
      return;
   }

   if (!isPreAlloc) {

      for (Int_t j = 0; j < n; j++) {
         if (oldStyle) {
            if (!start[j])
               start[j] = ((TClass *)cl)->New();
            ((TClass *)cl)->Streamer(start[j], *this);
            continue;
         }

         // delete the object or collection
         if (start[j] && TStreamerInfo::CanDelete())
            ((TClass *)cl)->Destructor(start[j], kFALSE); // call delete and desctructor
         start[j] = ReadObjectAny(cl);
      }

   } else { // case //-> in comment

      for (Int_t j = 0; j < n; j++) {
         if (!start[j])
            start[j] = ((TClass *)cl)->New();
         StreamObject(start[j], cl, onFileClass);
      }
   }

   if (gDebug > 2)
      Info("ReadFastArray", "(void ** Done");

   //   TBuffer::ReadFastArray(startp, cl, n, isPreAlloc, s);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads array size, written in raw data table.
/// Used in ReadArray methods, where TBuffer need to read array size first.

Int_t TBufferSQL2::SqlReadArraySize()
{
   const char *value = SqlReadValue(sqlio::Array);
   if (!value || (strlen(value) == 0))
      return 0;
   Int_t sz = atoi(value);
   return sz;
}

template <typename T>
R__ALWAYS_INLINE void TBufferSQL2::SqlWriteArray(T *arr, Int_t arrsize, Bool_t withsize)
{
   if (!withsize && (arrsize <= 0))
      return;
   PushStack()->SetArray(withsize ? arrsize : -1);
   Int_t indx = 0;
   if (fCompressLevel > 0) {
      while (indx < arrsize) {
         Int_t curr = indx++;
         while ((indx < arrsize) && (arr[indx] == arr[curr]))
            indx++;
         SqlWriteBasic(arr[curr]);
         Stack()->ChildArrayIndex(curr, indx - curr);
      }
   } else {
      for (; indx < arrsize; indx++) {
         SqlWriteBasic(arr[indx]);
         Stack()->ChildArrayIndex(indx, 1);
      }
   }
   PopStack();
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Bool_t to buffer

void TBufferSQL2::WriteArray(const Bool_t *b, Int_t n)
{
   SqlWriteArray(b, n, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Char_t to buffer

void TBufferSQL2::WriteArray(const Char_t *c, Int_t n)
{
   SqlWriteArray(c, n, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UChar_t to buffer

void TBufferSQL2::WriteArray(const UChar_t *c, Int_t n)
{
   SqlWriteArray(c, n, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Short_t to buffer

void TBufferSQL2::WriteArray(const Short_t *h, Int_t n)
{
   SqlWriteArray(h, n, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UShort_t to buffer

void TBufferSQL2::WriteArray(const UShort_t *h, Int_t n)
{
   SqlWriteArray(h, n, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Int_ to buffer

void TBufferSQL2::WriteArray(const Int_t *i, Int_t n)
{
   SqlWriteArray(i, n, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UInt_t to buffer

void TBufferSQL2::WriteArray(const UInt_t *i, Int_t n)
{
   SqlWriteArray(i, n, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Long_t to buffer

void TBufferSQL2::WriteArray(const Long_t *l, Int_t n)
{
   SqlWriteArray(l, n, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of ULong_t to buffer

void TBufferSQL2::WriteArray(const ULong_t *l, Int_t n)
{
   SqlWriteArray(l, n, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Long64_t to buffer

void TBufferSQL2::WriteArray(const Long64_t *l, Int_t n)
{
   SqlWriteArray(l, n, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of ULong64_t to buffer

void TBufferSQL2::WriteArray(const ULong64_t *l, Int_t n)
{
   SqlWriteArray(l, n, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Float_t to buffer

void TBufferSQL2::WriteArray(const Float_t *f, Int_t n)
{
   SqlWriteArray(f, n, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Double_t to buffer

void TBufferSQL2::WriteArray(const Double_t *d, Int_t n)
{
   SqlWriteArray(d, n, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Bool_t to buffer

void TBufferSQL2::WriteFastArray(const Bool_t *b, Int_t n)
{
   SqlWriteArray(b, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Char_t to buffer
/// it will be reproduced as CharStar node with string as attribute

void TBufferSQL2::WriteFastArray(const Char_t *c, Int_t n)
{
   Bool_t usedefault = (n == 0);

   const Char_t *ccc = c;
   // check if no zeros in the array
   if (!usedefault)
      for (int i = 0; i < n; i++)
         if (*ccc++ == 0) {
            usedefault = kTRUE;
            break;
         }

   if (usedefault) {
      SqlWriteArray(c, n);
   } else {
      Char_t *buf = new Char_t[n + 1];
      memcpy(buf, c, n);
      buf[n] = 0;
      SqlWriteValue(buf, sqlio::CharStar);
      delete[] buf;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UChar_t to buffer

void TBufferSQL2::WriteFastArray(const UChar_t *c, Int_t n)
{
   SqlWriteArray(c, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Short_t to buffer

void TBufferSQL2::WriteFastArray(const Short_t *h, Int_t n)
{
   SqlWriteArray(h, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UShort_t to buffer

void TBufferSQL2::WriteFastArray(const UShort_t *h, Int_t n)
{
   SqlWriteArray(h, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Int_t to buffer

void TBufferSQL2::WriteFastArray(const Int_t *i, Int_t n)
{
   SqlWriteArray(i, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of UInt_t to buffer

void TBufferSQL2::WriteFastArray(const UInt_t *i, Int_t n)
{
   SqlWriteArray(i, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Long_t to buffer

void TBufferSQL2::WriteFastArray(const Long_t *l, Int_t n)
{
   SqlWriteArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of ULong_t to buffer

void TBufferSQL2::WriteFastArray(const ULong_t *l, Int_t n)
{
   SqlWriteArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Long64_t to buffer

void TBufferSQL2::WriteFastArray(const Long64_t *l, Int_t n)
{
   SqlWriteArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of ULong64_t to buffer

void TBufferSQL2::WriteFastArray(const ULong64_t *l, Int_t n)
{
   SqlWriteArray(l, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Float_t to buffer

void TBufferSQL2::WriteFastArray(const Float_t *f, Int_t n)
{
   SqlWriteArray(f, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of Double_t to buffer

void TBufferSQL2::WriteFastArray(const Double_t *d, Int_t n)
{
   SqlWriteArray(d, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n characters into the I/O buffer.
/// Used only by TLeafC, just dummy implementation here

void TBufferSQL2::WriteFastArrayString(const Char_t *c, Int_t n)
{
   SqlWriteArray(c, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Same functionality as TBuffer::WriteFastArray(...) but
/// instead of calling cl->Streamer(obj,buf) call here
/// buf.StreamObject(obj, cl). In that case it is easy to understand where
/// object data is started and finished

void TBufferSQL2::WriteFastArray(void *start, const TClass *cl, Int_t n, TMemberStreamer *streamer)
{
   if (streamer) {
      StreamObjectExtra(start, streamer, cl, 0);
      //      (*streamer)(*this, start, 0);
      return;
   }

   char *obj = (char *)start;
   if (!n)
      n = 1;
   int size = cl->Size();

   for (Int_t j = 0; j < n; j++, obj += size)
      StreamObject(obj, cl);
}

////////////////////////////////////////////////////////////////////////////////
/// Same functionality as TBuffer::WriteFastArray(...) but
/// instead of calling cl->Streamer(obj,buf) call here
/// buf.StreamObject(obj, cl). In that case it is easy to understand where
/// object data is started and finished

Int_t TBufferSQL2::WriteFastArray(void **start, const TClass *cl, Int_t n, Bool_t isPreAlloc, TMemberStreamer *streamer)
{

   Bool_t oldStyle = kFALSE; // flag used to reproduce old-style I/O actions for kSTLp

   if ((fIOVersion < 2) && !isPreAlloc) {
      TStreamerElement *elem = Stack(0)->GetElement();
      if (elem && ((elem->GetType() == TStreamerInfo::kSTLp) ||
                   (elem->GetType() == TStreamerInfo::kSTLp + TStreamerInfo::kOffsetL)))
         oldStyle = kTRUE;
   }

   if (streamer) {
      if (oldStyle)
         (*streamer)(*this, (void *)start, n);
      else
         StreamObjectExtra((void *)start, streamer, cl, 0);
      return 0;
   }

   int strInfo = 0;

   Int_t res = 0;

   if (!isPreAlloc) {

      for (Int_t j = 0; j < n; j++) {
         // must write StreamerInfo if pointer is null
         if (!strInfo && !start[j] && !oldStyle)
            ForceWriteInfo(((TClass *)cl)->GetStreamerInfo(), kFALSE);
         strInfo = 2003;
         if (oldStyle)
            ((TClass *)cl)->Streamer(start[j], *this);
         else
            res |= WriteObjectAny(start[j], cl);
      }

   } else {
      // case //-> in comment

      for (Int_t j = 0; j < n; j++) {
         if (!start[j])
            start[j] = ((TClass *)cl)->New();
         StreamObject(start[j], cl);
      }
   }
   return res;

   //   return TBuffer::WriteFastArray(startp, cl, n, isPreAlloc, s);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream object to/from buffer

void TBufferSQL2::StreamObject(void *obj, const TClass *cl, const TClass *onFileClass)
{
   if (fIOVersion < 2) {
      TStreamerElement *elem = Stack(0)->GetElement();
      if (elem && (elem->GetType() == TStreamerInfo::kTObject)) {
         ((TObject *)obj)->TObject::Streamer(*this);
         return;
      } else if (elem && (elem->GetType() == TStreamerInfo::kTNamed)) {
         ((TNamed *)obj)->TNamed::Streamer(*this);
         return;
      }
   }

   if (gDebug > 1)
      Info("StreamObject", "class  %s", (cl ? cl->GetName() : "none"));
   if (IsReading())
      SqlReadObject(obj, 0, nullptr, 0, onFileClass);
   else
      SqlWriteObject(obj, cl, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream object to/from buffer

void TBufferSQL2::StreamObjectExtra(void *obj, TMemberStreamer *streamer, const TClass *cl, Int_t n,
                                    const TClass *onFileClass)
{
   if (!streamer)
      return;

   if (gDebug > 1)
      Info("StreamObjectExtra", "class = %s", cl->GetName());
   //   (*streamer)(*this, obj, n);

   if (IsReading())
      SqlReadObject(obj, 0, streamer, n, onFileClass);
   else
      SqlWriteObject(obj, cl, kTRUE, streamer, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Bool_t value from buffer

void TBufferSQL2::ReadBool(Bool_t &b)
{
   SqlReadBasic(b);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Char_t value from buffer

void TBufferSQL2::ReadChar(Char_t &c)
{
   SqlReadBasic(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads UChar_t value from buffer

void TBufferSQL2::ReadUChar(UChar_t &c)
{
   SqlReadBasic(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Short_t value from buffer

void TBufferSQL2::ReadShort(Short_t &h)
{
   SqlReadBasic(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads UShort_t value from buffer

void TBufferSQL2::ReadUShort(UShort_t &h)
{
   SqlReadBasic(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Int_t value from buffer

void TBufferSQL2::ReadInt(Int_t &i)
{
   SqlReadBasic(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads UInt_t value from buffer

void TBufferSQL2::ReadUInt(UInt_t &i)
{
   SqlReadBasic(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Long_t value from buffer

void TBufferSQL2::ReadLong(Long_t &l)
{
   SqlReadBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads ULong_t value from buffer

void TBufferSQL2::ReadULong(ULong_t &l)
{
   SqlReadBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Long64_t value from buffer

void TBufferSQL2::ReadLong64(Long64_t &l)
{
   SqlReadBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads ULong64_t value from buffer

void TBufferSQL2::ReadULong64(ULong64_t &l)
{
   SqlReadBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Float_t value from buffer

void TBufferSQL2::ReadFloat(Float_t &f)
{
   SqlReadBasic(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads Double_t value from buffer

void TBufferSQL2::ReadDouble(Double_t &d)
{
   SqlReadBasic(d);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads array of characters from buffer

void TBufferSQL2::ReadCharP(Char_t *c)
{
   const char *buf = SqlReadCharStarValue();
   if (buf)
      strcpy(c, buf);
}

////////////////////////////////////////////////////////////////////////////////
/// Read a TString

void TBufferSQL2::ReadTString(TString &s)
{
   if (fIOVersion < 2) {
      // original TBufferFile method can not be used, while used TString methods are private
      // try to reimplement close to the original
      Int_t nbig;
      UChar_t nwh;
      *this >> nwh;
      if (nwh == 0) {
         s.Resize(0);
      } else {
         if (nwh == 255)
            *this >> nbig;
         else
            nbig = nwh;

         char *data = new char[nbig];
         data[nbig] = 0;
         ReadFastArray(data, nbig);
         s = data;
         delete[] data;
      }
   } else {
      // TODO: new code - direct reading of string
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write a TString

void TBufferSQL2::WriteTString(const TString &s)
{
   if (fIOVersion < 2) {
      // original TBufferFile method, keep for compatibility
      Int_t nbig = s.Length();
      UChar_t nwh;
      if (nbig > 254) {
         nwh = 255;
         *this << nwh;
         *this << nbig;
      } else {
         nwh = UChar_t(nbig);
         *this << nwh;
      }
      const char *data = s.Data();
      WriteFastArray(data, nbig);
   } else {
      // TODO: make writing of string directly
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read a std::string

void TBufferSQL2::ReadStdString(std::string *obj)
{
   if (fIOVersion < 2) {
      if (!obj) {
         Error("ReadStdString", "The std::string address is nullptr but should not");
         return;
      }
      Int_t nbig;
      UChar_t nwh;
      *this >> nwh;
      if (nwh == 0) {
         obj->clear();
      } else {
         if (obj->size()) {
            // Insure that the underlying data storage is not shared
            (*obj)[0] = '\0';
         }
         if (nwh == 255) {
            *this >> nbig;
            obj->resize(nbig, '\0');
            ReadFastArray((char *)obj->data(), nbig);
         } else {
            obj->resize(nwh, '\0');
            ReadFastArray((char *)obj->data(), nwh);
         }
      }
   } else {
      // TODO: direct reading of std string
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write a std::string

void TBufferSQL2::WriteStdString(const std::string *obj)
{
   if (fIOVersion < 2) {
      if (!obj) {
         *this << (UChar_t)0;
         WriteFastArray("", 0);
         return;
      }

      UChar_t nwh;
      Int_t nbig = obj->length();
      if (nbig > 254) {
         nwh = 255;
         *this << nwh;
         *this << nbig;
      } else {
         nwh = UChar_t(nbig);
         *this << nwh;
      }
      WriteFastArray(obj->data(), nbig);
   } else {
      // TODO: make writing of string directly
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read a char* string

void TBufferSQL2::ReadCharStar(char *&s)
{
   delete[] s;
   s = nullptr;

   Int_t nch;
   *this >> nch;
   if (nch > 0) {
      s = new char[nch + 1];
      ReadFastArray(s, nch);
      s[nch] = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write a char* string

void TBufferSQL2::WriteCharStar(char *s)
{
   Int_t nch = 0;
   if (s) {
      nch = strlen(s);
      *this << nch;
      WriteFastArray(s, nch);
   } else {
      *this << nch;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Bool_t value to buffer

void TBufferSQL2::WriteBool(Bool_t b)
{
   SqlWriteBasic(b);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Char_t value to buffer

void TBufferSQL2::WriteChar(Char_t c)
{
   SqlWriteBasic(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes UChar_t value to buffer

void TBufferSQL2::WriteUChar(UChar_t c)
{
   SqlWriteBasic(c);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Short_t value to buffer

void TBufferSQL2::WriteShort(Short_t h)
{
   SqlWriteBasic(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes UShort_t value to buffer

void TBufferSQL2::WriteUShort(UShort_t h)
{
   SqlWriteBasic(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Int_t value to buffer

void TBufferSQL2::WriteInt(Int_t i)
{
   SqlWriteBasic(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes UInt_t value to buffer

void TBufferSQL2::WriteUInt(UInt_t i)
{
   SqlWriteBasic(i);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Long_t value to buffer

void TBufferSQL2::WriteLong(Long_t l)
{
   SqlWriteBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes ULong_t value to buffer

void TBufferSQL2::WriteULong(ULong_t l)
{
   SqlWriteBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Long64_t value to buffer

void TBufferSQL2::WriteLong64(Long64_t l)
{
   SqlWriteBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes ULong64_t value to buffer

void TBufferSQL2::WriteULong64(ULong64_t l)
{
   SqlWriteBasic(l);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Float_t value to buffer

void TBufferSQL2::WriteFloat(Float_t f)
{
   SqlWriteBasic(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes Double_t value to buffer

void TBufferSQL2::WriteDouble(Double_t d)
{
   SqlWriteBasic(d);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes array of characters to buffer

void TBufferSQL2::WriteCharP(const Char_t *c)
{
   SqlWriteValue(c, sqlio::CharStar);
}

////////////////////////////////////////////////////////////////////////////////
/// converts Char_t to string and creates correspondent sql structure

Bool_t TBufferSQL2::SqlWriteBasic(Char_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%d", value);
   return SqlWriteValue(buf, sqlio::Char);
}

////////////////////////////////////////////////////////////////////////////////
/// converts Short_t to string and creates correspondent sql structure

Bool_t TBufferSQL2::SqlWriteBasic(Short_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%hd", value);
   return SqlWriteValue(buf, sqlio::Short);
}

////////////////////////////////////////////////////////////////////////////////
/// converts Int_t to string and creates correspondent sql structure

Bool_t TBufferSQL2::SqlWriteBasic(Int_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%d", value);
   return SqlWriteValue(buf, sqlio::Int);
}

////////////////////////////////////////////////////////////////////////////////
/// converts Long_t to string and creates correspondent sql structure

Bool_t TBufferSQL2::SqlWriteBasic(Long_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%ld", value);
   return SqlWriteValue(buf, sqlio::Long);
}

////////////////////////////////////////////////////////////////////////////////
/// converts Long64_t to string and creates correspondent sql structure

Bool_t TBufferSQL2::SqlWriteBasic(Long64_t value)
{
   std::string buf = std::to_string(value);
   return SqlWriteValue(buf.c_str(), sqlio::Long64);
}

////////////////////////////////////////////////////////////////////////////////
/// converts Float_t to string and creates correspondent sql structure

Bool_t TBufferSQL2::SqlWriteBasic(Float_t value)
{
   char buf[200];
   ConvertFloat(value, buf, sizeof(buf), kTRUE);
   return SqlWriteValue(buf, sqlio::Float);
}

////////////////////////////////////////////////////////////////////////////////
/// converts Double_t to string and creates correspondent sql structure

Bool_t TBufferSQL2::SqlWriteBasic(Double_t value)
{
   char buf[200];
   ConvertDouble(value, buf, sizeof(buf), kTRUE);
   return SqlWriteValue(buf, sqlio::Double);
}

////////////////////////////////////////////////////////////////////////////////
/// converts Bool_t to string and creates correspondent sql structure

Bool_t TBufferSQL2::SqlWriteBasic(Bool_t value)
{
   return SqlWriteValue(value ? sqlio::True : sqlio::False, sqlio::Bool);
}

////////////////////////////////////////////////////////////////////////////////
/// converts UChar_t to string and creates correspondent sql structure

Bool_t TBufferSQL2::SqlWriteBasic(UChar_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%u", value);
   return SqlWriteValue(buf, sqlio::UChar);
}

////////////////////////////////////////////////////////////////////////////////
/// converts UShort_t to string and creates correspondent sql structure

Bool_t TBufferSQL2::SqlWriteBasic(UShort_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%hu", value);
   return SqlWriteValue(buf, sqlio::UShort);
}

////////////////////////////////////////////////////////////////////////////////
/// converts UInt_t to string and creates correspondent sql structure

Bool_t TBufferSQL2::SqlWriteBasic(UInt_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%u", value);
   return SqlWriteValue(buf, sqlio::UInt);
}

////////////////////////////////////////////////////////////////////////////////
/// converts ULong_t to string and creates correspondent sql structure

Bool_t TBufferSQL2::SqlWriteBasic(ULong_t value)
{
   char buf[50];
   snprintf(buf, sizeof(buf), "%lu", value);
   return SqlWriteValue(buf, sqlio::ULong);
}

////////////////////////////////////////////////////////////////////////////////
/// converts ULong64_t to string and creates correspondent sql structure

Bool_t TBufferSQL2::SqlWriteBasic(ULong64_t value)
{
   std::string buf = std::to_string(value);
   return SqlWriteValue(buf.c_str(), sqlio::ULong64);
}

//______________________________________________________________________________

Bool_t TBufferSQL2::SqlWriteValue(const char *value, const char *tname)
{
   // create structure in stack, which holds specified value

   Stack()->AddValue(value, tname);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Read current value from table and convert it to Char_t value

void TBufferSQL2::SqlReadBasic(Char_t &value)
{
   const char *res = SqlReadValue(sqlio::Char);
   if (res) {
      int n;
      sscanf(res, "%d", &n);
      value = n;
   } else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read current value from table and convert it to Short_t value

void TBufferSQL2::SqlReadBasic(Short_t &value)
{
   const char *res = SqlReadValue(sqlio::Short);
   if (res)
      sscanf(res, "%hd", &value);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read current value from table and convert it to Int_t value

void TBufferSQL2::SqlReadBasic(Int_t &value)
{
   const char *res = SqlReadValue(sqlio::Int);
   if (res)
      sscanf(res, "%d", &value);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read current value from table and convert it to Long_t value

void TBufferSQL2::SqlReadBasic(Long_t &value)
{
   const char *res = SqlReadValue(sqlio::Long);
   if (res)
      sscanf(res, "%ld", &value);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read current value from table and convert it to Long64_t value

void TBufferSQL2::SqlReadBasic(Long64_t &value)
{
   const char *res = SqlReadValue(sqlio::Long64);
   if (res)
      value = (Long64_t)std::stoll(res);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read current value from table and convert it to Float_t value

void TBufferSQL2::SqlReadBasic(Float_t &value)
{
   const char *res = SqlReadValue(sqlio::Float);
   if (res)
      sscanf(res, "%f", &value);
   else
      value = 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// Read current value from table and convert it to Double_t value

void TBufferSQL2::SqlReadBasic(Double_t &value)
{
   const char *res = SqlReadValue(sqlio::Double);
   if (res)
      sscanf(res, "%lf", &value);
   else
      value = 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// Read current value from table and convert it to Bool_t value

void TBufferSQL2::SqlReadBasic(Bool_t &value)
{
   const char *res = SqlReadValue(sqlio::Bool);
   if (res)
      value = (strcmp(res, sqlio::True) == 0);
   else
      value = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Read current value from table and convert it to UChar_t value

void TBufferSQL2::SqlReadBasic(UChar_t &value)
{
   const char *res = SqlReadValue(sqlio::UChar);
   if (res) {
      unsigned int n;
      sscanf(res, "%ud", &n);
      value = n;
   } else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read current value from table and convert it to UShort_t value

void TBufferSQL2::SqlReadBasic(UShort_t &value)
{
   const char *res = SqlReadValue(sqlio::UShort);
   if (res)
      sscanf(res, "%hud", &value);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read current value from table and convert it to UInt_t value

void TBufferSQL2::SqlReadBasic(UInt_t &value)
{
   const char *res = SqlReadValue(sqlio::UInt);
   if (res)
      sscanf(res, "%u", &value);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read current value from table and convert it to ULong_t value

void TBufferSQL2::SqlReadBasic(ULong_t &value)
{
   const char *res = SqlReadValue(sqlio::ULong);
   if (res)
      sscanf(res, "%lu", &value);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read current value from table and convert it to ULong64_t value

void TBufferSQL2::SqlReadBasic(ULong64_t &value)
{
   const char *res = SqlReadValue(sqlio::ULong64);
   if (res)
      value = (ULong64_t)std::stoull(res);
   else
      value = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read string value from current stack node

const char *TBufferSQL2::SqlReadValue(const char *tname)
{
   if (fErrorFlag > 0)
      return 0;

   if (!fCurrentData) {
      Error("SqlReadValue", "No object data to read from");
      fErrorFlag = 1;
      return 0;
   }

   if (!fIgnoreVerification)
      if (!fCurrentData->VerifyDataType(tname)) {
         fErrorFlag = 1;
         return 0;
      }

   fReadBuffer = fCurrentData->GetValue();

   fCurrentData->ShiftToNextValue();

   if (gDebug > 4)
      Info("SqlReadValue", "%s = %s", tname, fReadBuffer.Data());

   return fReadBuffer.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Read CharStar value, if it has special code, request it from large table

const char *TBufferSQL2::SqlReadCharStarValue()
{
   const char *res = SqlReadValue(sqlio::CharStar);
   if (!res || !fSQL)
      return nullptr;

   Long64_t objid = Stack()->DefineObjectId(kTRUE);

   Int_t strid = fSQL->IsLongStringCode(objid, res);
   if (strid <= 0)
      return res;

   fSQL->GetLongString(objid, strid, fReadBuffer);

   return fReadBuffer.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Push stack with structural information about streamed object

TSQLStructure *TBufferSQL2::PushStack()
{
   TSQLStructure *res = new TSQLStructure;
   if (!fStk) {
      fStructure = res;
   } else {
      fStk->Add(res);
   }

   fStk = res; // add in the stack
   return fStk;
}

////////////////////////////////////////////////////////////////////////////////
/// Pop stack

TSQLStructure *TBufferSQL2::PopStack()
{
   if (!fStk)
      return nullptr;
   fStk = fStk->GetParent();
   return fStk;
}

////////////////////////////////////////////////////////////////////////////////
/// returns head of stack

TSQLStructure *TBufferSQL2::Stack(Int_t depth)
{
   TSQLStructure *curr = fStk;
   while ((depth-- > 0) && curr)
      curr = curr->GetParent();
   return curr;
}

////////////////////////////////////////////////////////////////////////////////
/// Return current streamer info element

TVirtualStreamerInfo *TBufferSQL2::GetInfo()
{
   return Stack()->GetStreamerInfo();
}
