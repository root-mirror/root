// @(#)root/sqlite:$Id$
// Author: o.freyermuth <o.f@cern.ch>, 01/06/2013

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TSQLiteServer.h"
#include "TSQLiteResult.h"
#include "TSQLiteStatement.h"
#include "TSQLColumnInfo.h"
#include "TList.h"
#include "TSQLTableInfo.h"
#include "TSQLRow.h"

#include <sqlite3.h>

ClassImp(TSQLiteServer);

////////////////////////////////////////////////////////////////////////////////
/// Open a connection to an SQLite DB server. The db arguments should be
/// of the form "sqlite://<database>", e.g.:
/// "sqlite://test.sqlite" or "sqlite://:memory:" for a temporary database
/// in memory.
/// Note that for SQLite versions >= 3.7.7 the full string behind
/// "sqlite://" is handed to sqlite3_open_v2() with SQLITE_OPEN_URI activated,
/// so all URI accepted by it can be used.

TSQLiteServer::TSQLiteServer(const char *db, const char* /*uid*/, const char* /*pw*/)
{
   fSQLite = nullptr;
   fSrvInfo = "SQLite ";
   fSrvInfo += sqlite3_libversion();

   if (strncmp(db, "sqlite://", 9)) {
      TString givenProtocol(db, 9); // this TString-constructor allocs len+1 and does \0 termination already.
      Error("TSQLiteServer", "protocol in db argument should be sqlite it is %s",
            givenProtocol.Data());
      MakeZombie();
      return;
   }

   const char *dbase = db + 9;

#ifndef SQLITE_OPEN_URI
#define SQLITE_OPEN_URI 0x00000000
#endif
#if SQLITE_VERSION_NUMBER >= 3005000
   Int_t error = sqlite3_open_v2(dbase, &fSQLite, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_URI, NULL);
#else
   Int_t error = sqlite3_open(dbase, &fSQLite);
#endif

   if (error == 0) {
      // Set members of the abstract interface
      fType = "SQLite";
      fHost = "";
      fDB = dbase;
      // fPort != -1 means we are 'connected'
      fPort = 0;
   } else {
      Error("TSQLiteServer", "opening of %s failed with error: %d %s", dbase, sqlite3_errcode(fSQLite), sqlite3_errmsg(fSQLite));
      sqlite3_close(fSQLite);
      MakeZombie();
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Close SQLite DB.

TSQLiteServer::~TSQLiteServer()
{
   if (IsConnected()) {
      sqlite3_close(fSQLite);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Close connection to SQLite DB.

void TSQLiteServer::Close(Option_t *)
{
   if (!fSQLite) {
      return;
   }

   if (IsConnected()) {
      sqlite3_close(fSQLite);
      // Mark as disconnected:
      fPort = -1;
      fSQLite = nullptr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// submit "START TRANSACTION" query to database
/// return kTRUE, if successful

Bool_t TSQLiteServer::StartTransaction()
{
   return Exec("BEGIN TRANSACTION");
}

////////////////////////////////////////////////////////////////////////////////
/// Execute SQL command. Result object must be deleted by the user.
/// Returns a pointer to a TSQLResult object if successful, 0 otherwise.
/// The result object must be deleted by the user.

TSQLResult *TSQLiteServer::Query(const char *sql)
{
   if (!IsConnected()) {
      Error("Query", "not connected");
      return 0;
   }

   sqlite3_stmt *preparedStmt = nullptr;

   // -1 as we read until we encounter a \0.
   // NULL because we do not check which char was read last.
#if SQLITE_VERSION_NUMBER >= 3005000
   int retVal = sqlite3_prepare_v2(fSQLite, sql, -1, &preparedStmt, NULL);
#else
   int retVal = sqlite3_prepare(fSQLite, sql, -1, &preparedStmt, NULL);
#endif
   if (retVal != SQLITE_OK) {
      Error("Query", "SQL Error: %d %s", retVal, sqlite3_errmsg(fSQLite));
      return 0;
   }

   return new TSQLiteResult(preparedStmt);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute SQL command which does not produce any result sets.
/// Returns kTRUE if successful.

Bool_t TSQLiteServer::Exec(const char *sql)
{
   if (!IsConnected()) {
      Error("Exec", "not connected");
      return kFALSE;
   }

   char *sqlite_err_msg;
   int ret = sqlite3_exec(fSQLite, sql, NULL, NULL, &sqlite_err_msg);
   if (ret != SQLITE_OK) {
      Error("Exec", "SQL Error: %d %s", ret, sqlite_err_msg);
      sqlite3_free(sqlite_err_msg);
      return kFALSE;
   }
   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Select a database. Always returns non-zero for SQLite,
/// as only one DB exists per file.

Int_t TSQLiteServer::SelectDataBase(const char* /*dbname*/)
{
   Error("SelectDataBase", "SelectDataBase command makes no sense for SQLite!");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// List all available databases. Always returns 0 for SQLite,
/// as only one DB exists per file.

TSQLResult *TSQLiteServer::GetDataBases(const char* /*wild*/)
{
   Error("GetDataBases", "GetDataBases command makes no sense for SQLite!");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// List all tables in the specified database. Wild is for wildcarding
/// "t%" list all tables starting with "t".
/// Returns a pointer to a TSQLResult object if successful, 0 otherwise.
/// The result object must be deleted by the user.

TSQLResult *TSQLiteServer::GetTables(const char* /*dbname*/, const char *wild)
{
   if (!IsConnected()) {
      Error("GetTables", "not connected");
      return 0;
   }

   TString sql = "SELECT name FROM sqlite_master where type='table'";
   if (wild)
      sql += Form(" AND name LIKE '%s'", wild);

   return Query(sql);
}

////////////////////////////////////////////////////////////////////////////////
/// List all columns in specified table (database argument is ignored).
/// Wild is for wildcarding "t%" list all columns starting with "t".
/// Returns a pointer to a TSQLResult object if successful, 0 otherwise.
/// The result object must be deleted by the user.
/// For SQLite, this fails with wildcard, as the column names are not queryable!
/// If no wildcard is used, the result of PRAGMA table_info(table) is returned,
/// which contains the names in field 1.

TSQLResult *TSQLiteServer::GetColumns(const char* /*dbname*/, const char* table,
      const char* wild)
{
   if (!IsConnected()) {
      Error("GetColumns", "not connected");
      return 0;
   }

   if (wild) {
      Error("GetColumns", "Not implementable for SQLite as a query with wildcard, use GetFieldNames() after SELECT instead!");
      return nullptr;
   } else {
      TString sql = Form("PRAGMA table_info('%s')", table);
      return Query(sql);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Produces SQL table info.
/// Object must be deleted by user.

TSQLTableInfo *TSQLiteServer::GetTableInfo(const char* tablename)
{
   if (!IsConnected()) {
      Error("GetTableInfo", "not connected");
      return 0;
   }

   if ((tablename==0) || (*tablename==0)) return nullptr;

   TSQLResult *columnRes = GetColumns("", tablename);

   if (columnRes == nullptr) {
      Error("GetTableInfo", "could not query columns");
      return nullptr;
   }

   TList* lst = nullptr;

   TSQLRow *columnRow;

   while ((columnRow = columnRes->Next()) != nullptr) {
      if (!lst) {
         lst = new TList();
      }

      // Field 3 is 'notnull', i.e. if it is 0, column is nullable
      Bool_t isNullable = (strcmp(columnRow->GetField(3), "0") == 0);

      lst->Add(new TSQLColumnInfo(columnRow->GetField(1), // column name
                                  columnRow->GetField(2), // column type name
                                  isNullable,  // isNullable defined above
                                  -1,   // SQLite is totally free about types
                                  -1,   // SQLite imposes no declarable size-limits
                                  -1,   // Field length only available querying the field
                                  -1,   // no data scale in SQLite
                                  -1)); // SQLite does not enforce any sign(s)
      delete columnRow;
   }
   delete columnRes;

   // lst == NULL is ok as TSQLTableInfo accepts and handles this
   TSQLTableInfo*  info = new TSQLTableInfo(tablename,
                                            lst);

   return info;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a database. Always returns non-zero for SQLite,
/// as it has only one DB per file.

Int_t TSQLiteServer::CreateDataBase(const char* /*dbname*/)
{
   Error("CreateDataBase", "CreateDataBase command makes no sense for SQLite!");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Drop (i.e. delete) a database. Always returns non-zero for SQLite,
/// as it has only one DB per file.

Int_t TSQLiteServer::DropDataBase(const char* /*dbname*/)
{
   Error("DropDataBase", "DropDataBase command makes no sense for SQLite!");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Reload permission tables. Returns 0 if successful, non-zero
/// otherwise. User must have reload permissions.

Int_t TSQLiteServer::Reload()
{
   if (!IsConnected()) {
      Error("Reload", "not connected");
      return -1;
   }

   Error("Reload", "not implemented");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Shutdown the database server. Returns 0 if successful, non-zero
/// otherwise. Makes no sense for SQLite, always returns -1.

Int_t TSQLiteServer::Shutdown()
{
   if (!IsConnected()) {
      Error("Shutdown", "not connected");
      return -1;
   }

   Error("Shutdown", "not implemented");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// We assume prepared statements work for all SQLite-versions.
/// As we actually use the recommended sqlite3_prepare(),
/// or, if possible, sqlite3_prepare_v2(),
/// this already introduces the "compile time check".

Bool_t TSQLiteServer::HasStatement() const
{
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Produce TSQLiteStatement.

TSQLStatement* TSQLiteServer::Statement(const char *sql, Int_t)
{
   if (!sql || !*sql) {
      SetError(-1, "no query string specified", "Statement");
      return nullptr;
   }

   if (!IsConnected()) {
      Error("Statement", "not connected");
      return nullptr;
   }

   sqlite3_stmt *preparedStmt = nullptr;

   // -1 as we read until we encounter a \0.
   // NULL because we do not check which char was read last.
#if SQLITE_VERSION_NUMBER >= 3005000
   int retVal = sqlite3_prepare_v2(fSQLite, sql, -1, &preparedStmt, NULL);
#else
   int retVal = sqlite3_prepare(fSQLite, sql, -1, &preparedStmt, NULL);
#endif
   if (retVal != SQLITE_OK) {
      Error("Statement", "SQL Error: %d %s", retVal, sqlite3_errmsg(fSQLite));
      return nullptr;
   }

   SQLite3_Stmt_t *stmt = new SQLite3_Stmt_t;
   stmt->fConn = fSQLite;
   stmt->fRes  = preparedStmt;

   return new TSQLiteStatement(stmt);
}

////////////////////////////////////////////////////////////////////////////////
/// Return server info, must be deleted by user.

const char *TSQLiteServer::ServerInfo()
{
   if (!IsConnected()) {
      Error("ServerInfo", "not connected");
      return 0;
   }

   return fSrvInfo.Data();
}
