// @(#)root/oracle:$Id$
// Author: Yan Liu and Shaowen Wang   23/11/04

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TOracleResult.h"
#include "TOracleRow.h"
#include "TList.h"

using namespace std;
using namespace oracle::occi;

ClassImp(TOracleResult)

////////////////////////////////////////////////////////////////////////////////
/// Oracle query result.

void TOracleResult::initResultSet(Statement *stmt)
{
   if (!stmt) {
      Error("initResultSet", "construction: empty statement");
   } else {
      try {
         fStmt = stmt;
         if (stmt->status() == Statement::RESULT_SET_AVAILABLE) {
            fResultType  = 1;
            fResult      = stmt->getResultSet();
            fFieldInfo   = (fResult==0) ? 0 : new vector<MetaData>(fResult->getColumnListMetaData());
            fFieldCount  = (fFieldInfo==0) ? 0 : fFieldInfo->size();
         } else if (stmt->status() == Statement::UPDATE_COUNT_AVAILABLE) {
            fResultType  = 3; // this is update_count_available
            fResult      = 0;
            fFieldInfo   = 0;
            fFieldCount  = 0;
            fUpdateCount = stmt->getUpdateCount();
         }
      } catch (SQLException &oraex) {
         Error("initResultSet", "%s", (oraex.getMessage()).c_str());
         MakeZombie();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

TOracleResult::TOracleResult(Connection *conn, Statement *stmt)
{
   fConn        = conn;
   fResult      = 0;
   fStmt        = 0;
   fPool        = 0;
   fRowCount    = 0;
   fFieldInfo   = 0;
   fResultType  = 0;
   fUpdateCount = 0;

   initResultSet(stmt);

   if (fResult) ProducePool();
}

////////////////////////////////////////////////////////////////////////////////
/// This construction func is only used to get table metainfo.

TOracleResult::TOracleResult(Connection *conn, const char *tableName)
{
   fResult      = 0;
   fStmt        = 0;
   fConn        = 0;
   fPool        = 0;
   fRowCount    = 0;
   fFieldInfo   = 0;
   fResultType  = 0;
   fUpdateCount = 0;
   fFieldCount  = 0;

   if (!tableName || !conn) {
      Error("TOracleResult", "construction: empty input parameter");
   } else {
      MetaData connMD = conn->getMetaData(tableName, MetaData::PTYPE_TABLE);
      fFieldInfo   = new vector<MetaData>(connMD.getVector(MetaData::ATTR_LIST_COLUMNS));
      fFieldCount  = fFieldInfo->size();
      fResultType  = 2; // indicates that this is just an table metainfo
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup Oracle query result.

TOracleResult::~TOracleResult()
{
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Close query result.

void TOracleResult::Close(Option_t *)
{
   if (fConn && fStmt) {
      if (fResult) fStmt->closeResultSet(fResult);
      fConn->terminateStatement(fStmt);
   }

   if (fPool) {
      fPool->Delete();
      delete fPool;
   }

   if (fFieldInfo)
      delete fFieldInfo;

   fResultType = 0;

   fStmt = 0;
   fResult = 0;
   fFieldInfo = 0;
   fPool = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if result set is open and field index within range.

Bool_t TOracleResult::IsValid(Int_t field)
{
   if (field < 0 || field >= fFieldCount) {
      Error("IsValid", "field index out of bounds");
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get number of fields in result.

Int_t TOracleResult::GetFieldCount()
{
   return fFieldCount;
}

////////////////////////////////////////////////////////////////////////////////
/// Get name of specified field.

const char *TOracleResult::GetFieldName(Int_t field)
{
   if (!IsValid(field))
      return 0;
   fNameBuffer = (*fFieldInfo)[field].getString(MetaData::ATTR_NAME);
   return fNameBuffer.c_str();
}

////////////////////////////////////////////////////////////////////////////////
/// Get next query result row. The returned object must be
/// deleted by the user.

TSQLRow *TOracleResult::Next()
{
   if (!fResult || (fResultType!=1)) return 0;

   if (fPool!=0) {
      TSQLRow* row = (TSQLRow*) fPool->First();
      if (row!=0) fPool->Remove(row);
      return row;
   }

   // if select query,
   try {
      if (fResult->next() != oracle::occi::ResultSet::END_OF_FETCH) {
         fRowCount++;
         return new TOracleRow(fResult, fFieldInfo);
      } else
         return 0;
   } catch (SQLException &oraex) {
      Error("Next", "%s", (oraex.getMessage()).c_str());
      MakeZombie();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

Int_t TOracleResult::GetRowCount() const
{
   if (!fResult) return 0;

   if (fPool==0) ((TOracleResult*) this)->ProducePool();

   return fRowCount;
}

////////////////////////////////////////////////////////////////////////////////

void TOracleResult::ProducePool()
{
   if (fPool!=0) return;

   TList* pool = new TList;
   TSQLRow* res = 0;
   while ((res = Next()) !=0) {
      pool->Add(res);
   }

   fPool = pool;
}
