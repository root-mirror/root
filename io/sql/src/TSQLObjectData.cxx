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
\class TSQLObjectData
\ingroup IO

TSQLObjectData is used in TBufferSQL2 class in reading procedure.
It contains data, request from database table for one specific
object for one specific class. For instance, when data for
class TH1 required, requests will be done to
TH1_ver4 and TH1_raw4 tables and result of these requests
will be kept in single TSQLObjectData instance.
*/

#include "TSQLObjectData.h"

#include "TObjArray.h"
#include "TNamed.h"
#include "TList.h"
#include "TSQLRow.h"
#include "TSQLResult.h"
#include "TSQLClassInfo.h"
#include "TSQLStructure.h"
#include "TSQLStatement.h"

/**
\class TSQLObjectInfo
\ingroup IO
Info (classname, version) about object in database
*/

ClassImp(TSQLObjectInfo);

////////////////////////////////////////////////////////////////////////////////

TSQLObjectInfo::TSQLObjectInfo() : TObject(), fObjId(0), fClassName(), fVersion(0)
{
}

////////////////////////////////////////////////////////////////////////////////

TSQLObjectInfo::TSQLObjectInfo(Long64_t objid, const char *classname, Version_t version)
   : TObject(), fObjId(objid), fClassName(classname), fVersion(version)
{
}

////////////////////////////////////////////////////////////////////////////////

TSQLObjectInfo::~TSQLObjectInfo()
{
}

/**
\class TSQLObjectData
\ingroup IO
Keeps the data requested from the SQL server for an object.
*/

ClassImp(TSQLObjectData);

////////////////////////////////////////////////////////////////////////////////
/// default contrsuctor

TSQLObjectData::TSQLObjectData()
   : TObject(), fInfo(nullptr), fObjId(0), fOwner(kFALSE), fClassData(nullptr), fBlobData(nullptr), fBlobStmt(nullptr), fLocatedColumn(-1),
     fClassRow(nullptr), fBlobRow(nullptr), fLocatedField(nullptr), fLocatedValue(nullptr), fCurrentBlob(kFALSE), fBlobPrefixName(nullptr),
     fBlobTypeName(nullptr), fUnpack(nullptr)
{
}

////////////////////////////////////////////////////////////////////////////////
/// normal contrsuctor,

TSQLObjectData::TSQLObjectData(TSQLClassInfo *sqlinfo, Long64_t objid, TSQLResult *classdata, TSQLRow *classrow,
                               TSQLResult *blobdata, TSQLStatement *blobstmt)
   : TObject(), fInfo(sqlinfo), fObjId(objid), fOwner(kFALSE), fClassData(classdata), fBlobData(blobdata),
     fBlobStmt(blobstmt), fLocatedColumn(-1), fClassRow(classrow), fBlobRow(nullptr), fLocatedField(nullptr), fLocatedValue(nullptr),
     fCurrentBlob(kFALSE), fBlobPrefixName(nullptr), fBlobTypeName(nullptr), fUnpack(nullptr)
{
   // take ownership if no special row from data pool is provided
   if ((fClassData != nullptr) && (fClassRow == nullptr)) {
      fOwner = kTRUE;
      fClassRow = fClassData->Next();
   }

   ShiftBlobRow();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor of TSQLObjectData object

TSQLObjectData::~TSQLObjectData()
{
   if ((fClassData != nullptr) && fOwner)
      delete fClassData;
   if (fClassRow != nullptr)
      delete fClassRow;
   if (fBlobRow != nullptr)
      delete fBlobRow;
   if (fBlobData != nullptr)
      delete fBlobData;
   if (fUnpack != nullptr) {
      fUnpack->Delete();
      delete fUnpack;
   }
   if (fBlobStmt != nullptr)
      delete fBlobStmt;
}

////////////////////////////////////////////////////////////////////////////////
/// return number of columns in class table result

Int_t TSQLObjectData::GetNumClassFields()
{
   if (fClassData != nullptr)
      return fClassData->GetFieldCount();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// get name of class table column

const char *TSQLObjectData::GetClassFieldName(Int_t n)
{
   if (fClassData != nullptr)
      return fClassData->GetFieldName(n);
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// locate column of that name in results

Bool_t TSQLObjectData::LocateColumn(const char *colname, Bool_t isblob)
{
   if (fUnpack != nullptr) {
      fUnpack->Delete();
      delete fUnpack;
      fUnpack = nullptr;
   }

   fLocatedField = nullptr;
   fLocatedValue = nullptr;
   fCurrentBlob = kFALSE;

   if ((fClassData == nullptr) || (fClassRow == nullptr))
      return kFALSE;

   //   Int_t numfields = GetNumClassFields();

   Int_t ncol = fInfo->FindColumn(colname, kFALSE);
   if (ncol > 0) {
      fLocatedColumn = ncol;
      fLocatedField = GetClassFieldName(ncol);
      fLocatedValue = fClassRow->GetField(ncol);
   }

   /*   for (Int_t ncol=1;ncol<numfields;ncol++) {
         const char* fieldname = GetClassFieldName(ncol);
         if (strcmp(colname, fieldname)==0) {
            fLocatedColumn = ncol;
            fLocatedField = fieldname;
            fLocatedValue = fClassRow->GetField(ncol);
            break;
         }
      }
   */

   if (fLocatedField == nullptr)
      return kFALSE;

   if (!isblob)
      return kTRUE;

   if ((fBlobRow == nullptr) && (fBlobStmt == nullptr))
      return kFALSE;

   fCurrentBlob = kTRUE;

   ExtractBlobValues();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// shift cursor to next blob value

Bool_t TSQLObjectData::ShiftBlobRow()
{
   if (fBlobStmt != nullptr) {
      Bool_t res = fBlobStmt->NextResultRow();
      if (!res) {
         delete fBlobStmt;
         fBlobStmt = nullptr;
      }
      return res;
   }

   delete fBlobRow;
   fBlobRow = fBlobData ? fBlobData->Next() : nullptr;
   return fBlobRow != nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// extract from curent blob row value and names identifiers

Bool_t TSQLObjectData::ExtractBlobValues()
{
   const char *name = nullptr;

   Bool_t hasdata = kFALSE;

   if (fBlobStmt != nullptr) {
      name = fBlobStmt->GetString(0);
      fLocatedValue = fBlobStmt->GetString(1);
      hasdata = kTRUE;
   }

   if (!hasdata) {
      if (fBlobRow != nullptr) {
         fLocatedValue = fBlobRow->GetField(1);
         name = fBlobRow->GetField(0);
      }
   }

   if (name == nullptr) {
      fBlobPrefixName = nullptr;
      fBlobTypeName = nullptr;
      return kFALSE;
   }

   const char *separ = strstr(name, ":"); // SQLNameSeparator()

   if (separ == nullptr) {
      fBlobPrefixName = nullptr;
      fBlobTypeName = name;
   } else {
      fBlobPrefixName = name;
      separ += strlen(":"); // SQLNameSeparator()
      fBlobTypeName = separ;
   }

   //   if (gDebug>4)
   //      Info("ExtractBlobValues","Prefix:%s Type:%s",
   //            (fBlobPrefixName ? fBlobPrefixName : "null"),
   //            (fBlobTypeName ? fBlobTypeName : "null"));

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// add emulated data
/// this used to place version or TObject raw data, read from normal tables

void TSQLObjectData::AddUnpack(const char *tname, const char *value)
{
   TNamed *str = new TNamed(tname, value);
   if (fUnpack == nullptr) {
      fUnpack = new TObjArray();
      fBlobPrefixName = nullptr;
      fBlobTypeName = str->GetName();
      fLocatedValue = str->GetTitle();
   }

   fUnpack->Add(str);
}

////////////////////////////////////////////////////////////////////////////////
/// emulate integer value in raw data

void TSQLObjectData::AddUnpackInt(const char *tname, Int_t value)
{
   TString sbuf;
   sbuf.Form("%d", value);
   AddUnpack(tname, sbuf.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// shift to next column or next row in blob data

void TSQLObjectData::ShiftToNextValue()
{
   Bool_t doshift = kTRUE;

   if (fUnpack != nullptr) {
      TObject *prev = fUnpack->First();
      fUnpack->Remove(prev);
      delete prev;
      fUnpack->Compress();
      if (fUnpack->GetLast() >= 0) {
         TNamed *curr = (TNamed *)fUnpack->First();
         fBlobPrefixName = nullptr;
         fBlobTypeName = curr->GetName();
         fLocatedValue = curr->GetTitle();
         return;
      }
      delete fUnpack;
      fUnpack = nullptr;
      doshift = kFALSE;
   }

   if (fCurrentBlob) {
      if (doshift)
         ShiftBlobRow();
      ExtractBlobValues();
   } else if (fClassData != nullptr) {
      if (doshift)
         fLocatedColumn++;
      if (fLocatedColumn < GetNumClassFields()) {
         fLocatedField = GetClassFieldName(fLocatedColumn);
         fLocatedValue = fClassRow->GetField(fLocatedColumn);
      } else {
         fLocatedField = nullptr;
         fLocatedValue = nullptr;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// checks if data type corresponds to that stored in raw table

Bool_t TSQLObjectData::VerifyDataType(const char *tname, Bool_t errormsg)
{
   if (tname == nullptr) {
      if (errormsg)
         Error("VerifyDataType", "Data type not specified");
      return kFALSE;
   }

   // here maybe type of column can be checked
   if (!IsBlobData())
      return kTRUE;

   if (gDebug > 4)
      if ((fBlobTypeName == nullptr) && errormsg) {
         Error("VerifyDataType", "fBlobTypeName is null");
         return kFALSE;
      }

   TString v1(fBlobTypeName);
   TString v2(tname);

   //   if (strcmp(fBlobTypeName,tname)!=0) {
   if (v1 != v2) {
      if (errormsg)
         Error("VerifyDataType", "Data type missmatch %s - %s", fBlobTypeName, tname);
      return kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// prepare to read data from raw table

Bool_t TSQLObjectData::PrepareForRawData()
{
   if (!ExtractBlobValues())
      return kFALSE;

   fCurrentBlob = kTRUE;

   return kTRUE;
}

//===================================================================================

//________________________________________________________________________
//
// TSQLObjectDataPool contains list (pool) of data from single class table
// for differents objects, all belonging to the same key.
// This is typical situation when list of objects stored as single key.
// To optimize reading of such data, one query is submitted and results of that
// query kept in TSQLObjectDataPool object
//
//________________________________________________________________________

ClassImp(TSQLObjectDataPool);

////////////////////////////////////////////////////////////////////////////////

TSQLObjectDataPool::TSQLObjectDataPool() : TObject(), fInfo(nullptr), fClassData(nullptr), fIsMoreRows(kTRUE), fRowsPool(nullptr)
{
}

////////////////////////////////////////////////////////////////////////////////

TSQLObjectDataPool::TSQLObjectDataPool(TSQLClassInfo *info, TSQLResult *data)
   : TObject(), fInfo(info), fClassData(data), fIsMoreRows(kTRUE), fRowsPool(nullptr)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor of TSQLObjectDataPool class
/// Deletes not used rows and class data table

TSQLObjectDataPool::~TSQLObjectDataPool()
{
   if (fClassData != nullptr)
      delete fClassData;
   if (fRowsPool != nullptr) {
      fRowsPool->Delete();
      delete fRowsPool;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns single sql row with object data for that class

TSQLRow *TSQLObjectDataPool::GetObjectRow(Long64_t objid)
{
   if (fClassData == nullptr)
      return nullptr;

   Long64_t rowid;

   if (fRowsPool != nullptr) {
      TObjLink *link = fRowsPool->FirstLink();
      while (link != nullptr) {
         TSQLRow *row = (TSQLRow *)link->GetObject();
         rowid = sqlio::atol64(row->GetField(0));
         if (rowid == objid) {
            fRowsPool->Remove(link);
            return row;
         }

         link = link->Next();
      }
   }

   while (fIsMoreRows) {
      TSQLRow *row = fClassData->Next();
      if (row == nullptr)
         fIsMoreRows = kFALSE;
      else {
         rowid = sqlio::atol64(row->GetField(0));
         if (rowid == objid)
            return row;
         if (fRowsPool == nullptr)
            fRowsPool = new TList();
         fRowsPool->Add(row);
      }
   }

   return nullptr;
}
