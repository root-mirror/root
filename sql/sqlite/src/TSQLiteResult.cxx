// @(#)root/sqlite:$Id$
// Author: o.freyermuth <o.f@cern.ch>, 01/06/2013

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TSQLiteResult.h"
#include "TSQLiteRow.h"

#include <sqlite3.h>

ClassImp(TSQLiteResult);

////////////////////////////////////////////////////////////////////////////////
/// SQLite query result.

TSQLiteResult::TSQLiteResult(void *result)
{
   fResult     = (sqlite3_stmt *) result;

   // RowCount is -1, as sqlite cannot determine RowCount beforehand:
   fRowCount = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup SQLite query result.

TSQLiteResult::~TSQLiteResult()
{
   if (fResult)
      Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Close query result.

void TSQLiteResult::Close(Option_t *)
{
   if (!fResult)
      return;

   sqlite3_finalize(fResult);
   fResult     = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if result set is open and field index within range.

Bool_t TSQLiteResult::IsValid(Int_t field)
{
   if (!fResult) {
      Error("IsValid", "result set closed");
      return kFALSE;
   }
   if (field < 0 || field >= GetFieldCount()) {
      Error("IsValid", "field index out of bounds");
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get number of fields in result.

Int_t TSQLiteResult::GetFieldCount()
{
   if (!fResult) {
      Error("GetFieldCount", "result set closed");
      return 0;
   }
   return sqlite3_column_count(fResult);
}

////////////////////////////////////////////////////////////////////////////////
/// Get name of specified field.

const char *TSQLiteResult::GetFieldName(Int_t field)
{
   if (!fResult) {
      Error("GetFieldName", "result set closed");
      return nullptr;
   }
   return sqlite3_column_name(fResult, field);
}

////////////////////////////////////////////////////////////////////////////////
/// SQLite can not determine the row count for a Query, return -1 instead.
/// For similar functionality, call Next() until it retruns nullptr.

Int_t TSQLiteResult::GetRowCount() const
{
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get next query result row. The returned object must be
/// deleted by the user.

TSQLRow *TSQLiteResult::Next()
{
   if (!fResult) {
      Error("Next", "result set closed");
      return nullptr;
   }

   int ret = sqlite3_step(fResult);
   if ((ret != SQLITE_DONE) && (ret != SQLITE_ROW)) {
      Error("Statement", "SQL Error: %d %s", ret, sqlite3_errmsg(sqlite3_db_handle(fResult)));
      return nullptr;
   }
   if (ret == SQLITE_DONE) {
      // Finished executing, no other row!
      return nullptr;
   }
   return new TSQLiteRow((void *) fResult, -1);
}

