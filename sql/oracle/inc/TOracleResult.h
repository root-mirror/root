// @(#)root/physics:$Id$
// Author: Yan Liu and Shaowen Wang   23/11/04

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TOracleResult
#define ROOT_TOracleResult

#include "TSQLResult.h"

#include <vector>

#ifndef R__WIN32
#include <sys/time.h>
#endif

#include <occi.h>

#ifdef CONST
#undef CONST
#endif

class TList;

class TOracleResult : public TSQLResult {

private:
   oracle::occi::Connection   *fConn{nullptr};               // connection to Oracle
   oracle::occi::Statement    *fStmt{nullptr};               // executed statement
   oracle::occi::ResultSet    *fResult{nullptr};             // query result (rows)
   std::vector<oracle::occi::MetaData> *fFieldInfo{nullptr}; // info for each field in the row
   Int_t                       fFieldCount{0};               // num of fields in resultset
   UInt_t                      fUpdateCount{0};              // for dml query, mutual exclusive with above
   Int_t                       fResultType{0};               // 0 - nothing; 1 - Select; 2 - table metainfo, 3 - update counter
   TList                      *fPool{nullptr};               // array of results, produced when number of rows are requested
   std::string                 fNameBuffer;                  // buffer for GetFieldName() argument

   Bool_t  IsValid(Int_t field);

   TOracleResult(const TOracleResult&);            // Not implemented;
   TOracleResult &operator=(const TOracleResult&); // Not implemented;

protected:
   void    initResultSet(oracle::occi::Statement *stmt);
   void    ProducePool();

public:
   TOracleResult(oracle::occi::Connection *conn, oracle::occi::Statement *stmt);
   TOracleResult(oracle::occi::Connection *conn, const char *tableName);
   ~TOracleResult();

   void        Close(Option_t *opt="") final;
   Int_t       GetFieldCount() final;
   const char *GetFieldName(Int_t field) final;
   Int_t       GetRowCount() const final;
   TSQLRow    *Next() final;

   Int_t       GetUpdateCount() const { return fUpdateCount; }

   ClassDefOverride(TOracleResult,0)  // Oracle query result
};

#endif
