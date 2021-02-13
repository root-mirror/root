// @(#)root/sqlite:$Id$
// Author: o.freyermuth <o.f@cern.ch>, 01/06/2013

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TSQLiteServer
#define ROOT_TSQLiteServer

#include "TSQLServer.h"

class sqlite3;

class TSQLiteServer : public TSQLServer {

 private:
   TString   fSrvInfo;   // Server info string
   sqlite3  *fSQLite{nullptr};    // connection to SQLite DB

 public:
   TSQLiteServer(const char *db, const char *uid = nullptr, const char *pw = nullptr);
   ~TSQLiteServer();

   void           Close(Option_t *opt = "") final;
   Bool_t         StartTransaction() final;
   TSQLResult    *Query(const char *sql) final;
   Bool_t         Exec(const char *sql) final;
   TSQLStatement *Statement(const char *sql, Int_t = 100) final;
   Bool_t         HasStatement() const final;
   Int_t          SelectDataBase(const char *dbname) final;
   TSQLResult    *GetDataBases(const char *wild = nullptr) final;
   TSQLResult    *GetTables(const char *dbname, const char *wild = nullptr) final;
   TSQLResult    *GetColumns(const char *dbname, const char *table, const char *wild = nullptr) final;
   TSQLTableInfo *GetTableInfo(const char* tablename) final;
   Int_t          CreateDataBase(const char *dbname) final;
   Int_t          DropDataBase(const char *dbname) final;
   Int_t          Reload() final;
   Int_t          Shutdown() final;
   const char    *ServerInfo() final;

   ClassDefOverride(TSQLiteServer,0);  // Connection to SQLite DB
};

#endif
