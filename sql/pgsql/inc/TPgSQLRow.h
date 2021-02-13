// @(#)root/pgsql:$Id$
// Author: g.p.ciceri <gp.ciceri@acm.org> 01/06/2001

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPgSQLRow
#define ROOT_TPgSQLRow

#include "TSQLRow.h"

#include <libpq-fe.h>

class TPgSQLRow : public TSQLRow {

private:
   PGresult *fResult{nullptr};       // current result set
   ULong_t   fRowNum{0};       // row number

   Bool_t  IsValid(Int_t field);

public:
   TPgSQLRow(void *result, ULong_t rowHandle);
   ~TPgSQLRow();

   void        Close(Option_t *opt="") final;
   ULong_t     GetFieldLength(Int_t field) final;
   const char *GetField(Int_t field) final;

   ClassDefOverride(TPgSQLRow,0)  // One row of PgSQL query result
};

#endif
