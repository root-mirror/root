// @(#)root/pgsql:$Id$
// Author: g.p.ciceri <gp.ciceri@acm.org> 01/06/2001

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPgSQLResult
#define ROOT_TPgSQLResult

#include "TSQLResult.h"

#include <libpq-fe.h>

class TPgSQLResult : public TSQLResult {

private:
   PGresult   *fResult{nullptr};      // query result (rows)
   ULong_t     fCurrentRow{0};        // info to result row

   Bool_t  IsValid(Int_t field);

public:
   TPgSQLResult(void *result);
   ~TPgSQLResult();

   void        Close(Option_t *opt="") final;
   Int_t       GetFieldCount() final;
   const char *GetFieldName(Int_t field) final;
   TSQLRow    *Next() final;

   ClassDefOverride(TPgSQLResult, 0)  // PgSQL query result
};

#endif
