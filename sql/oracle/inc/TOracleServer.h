// @(#)root/physics:$Id$
// Author: Yan Liu and Shaowen Wang   23/11/04

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TOracleServer
#define ROOT_TOracleServer

#include "TSQLServer.h"

#ifndef R__WIN32
#include <sys/time.h>
#endif

#include <occi.h>

#ifdef CONST
#undef CONST
#endif

class TOracleServer : public TSQLServer {

private:
   oracle::occi::Environment  *fEnv{nullptr};    // environment of Oracle access
   oracle::occi::Connection   *fConn{nullptr};   // connection to Oracle server
   TString       fInfo;  // info string with Oracle version information

   static const char* fgDatimeFormat; //! format for converting date and time stamps into string

public:
   TOracleServer(const char *db, const char *uid, const char *pw);
   ~TOracleServer();

   void        Close(Option_t *opt="") final;
   TSQLResult *Query(const char *sql) final;
   Bool_t      Exec(const char* sql) final;
   TSQLStatement *Statement(const char *sql, Int_t niter = 100) final;
   Bool_t      IsConnected() const final { return (fConn!=0) && (fEnv!=0); }
   Bool_t      HasStatement() const final { return kTRUE; }
   Int_t       SelectDataBase(const char *dbname) final;
   TSQLResult *GetDataBases(const char *wild = nullptr) final;
   TSQLResult *GetTables(const char *dbname, const char *wild = nullptr) final;
   TList      *GetTablesList(const char* wild = nullptr) final;
   TSQLTableInfo *GetTableInfo(const char* tablename) final;
   TSQLResult *GetColumns(const char *dbname, const char *table, const char *wild = nullptr) final;
   Int_t       GetMaxIdentifierLength() final { return 30; }
   Int_t       CreateDataBase(const char *dbname) final;
   Int_t       DropDataBase(const char *dbname) final;
   Int_t       Reload() final;
   Int_t       Shutdown() final;
   const char *ServerInfo() final;

   Bool_t      StartTransaction() final;
   Bool_t      Commit() final;
   Bool_t      Rollback() final;

   static    void     SetDatimeFormat(const char* fmt = "MM/DD/YYYY, HH24:MI:SS");
   static const char* GetDatimeFormat();

   ClassDefOverride(TOracleServer,0)  // Connection to Oracle server
};

#endif
