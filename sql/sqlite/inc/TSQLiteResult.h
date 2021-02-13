// @(#)root/sqlite:$Id$
// Author: o.freyermuth <o.f@cern.ch>, 01/06/2013

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TSQLiteResult
#define ROOT_TSQLiteResult

#include "TSQLResult.h"

class sqlite3_stmt;

class TSQLiteResult : public TSQLResult {

private:
   sqlite3_stmt   *fResult{nullptr};  // query result (rows)

   Bool_t  IsValid(Int_t field);

public:
   TSQLiteResult(void *result);
   ~TSQLiteResult();

   void        Close(Option_t *opt="") final;
   Int_t       GetFieldCount() final;
   const char *GetFieldName(Int_t field) final;
   Int_t       GetRowCount() const final;
   TSQLRow    *Next() final;

   ClassDefOverride(TSQLiteResult, 0)  // SQLite query result
};

#endif
