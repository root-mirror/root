// @(#)root/sqlite:
// Author: o.freyermuth <o.f@cern.ch>, 01/06/2013

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TSQLiteRow
#define ROOT_TSQLiteRow

#include "TSQLRow.h"

class sqlite3_stmt;

class TSQLiteRow : public TSQLRow {

private:
   sqlite3_stmt *fResult{nullptr};       ///<! current result set
   Bool_t        IsValid(Int_t field);

public:
   TSQLiteRow(void *result, ULong_t rowHandle);
   ~TSQLiteRow();

   void        Close(Option_t *opt="") final;
   ULong_t     GetFieldLength(Int_t field) final;
   const char *GetField(Int_t field) final;

   ClassDefOverride(TSQLiteRow,0)  // One row of SQLite query result
};

#endif
