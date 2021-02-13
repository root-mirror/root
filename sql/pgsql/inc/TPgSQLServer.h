// @(#)root/pgsql:$Id$
// Author: g.p.ciceri <gp.ciceri@acm.org> 01/06/2001

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPgSQLServer
#define ROOT_TPgSQLServer

#include "TSQLServer.h"

#include <map>
#include <string>

#include <libpq-fe.h>

class TPgSQLServer : public TSQLServer {

private:
   PGconn  *fPgSQL{nullptr};    // connection to PgSQL server
   TString  fSrvInfo;  // Server info
   std::map<Int_t,std::string> fOidTypNameMap; // Map of oid to typname, used in GetTableInfo()
public:
   TPgSQLServer(const char *db, const char *uid, const char *pw);
   ~TPgSQLServer();

   void           Close(Option_t *opt="") final;
   TSQLResult    *Query(const char *sql) final;
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

   ClassDefOverride(TPgSQLServer,0)  // Connection to PgSQL server
};

#endif
