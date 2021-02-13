// @(#)root/physics:$Id$
// Author: Yan Liu and Shaowen Wang   23/11/04

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TOracleRow
#define ROOT_TOracleRow

#include "TSQLRow.h"

#include <occi.h>

#ifdef CONST
#undef CONST
#endif

class TOracleRow : public TSQLRow {

private:
   oracle::occi::ResultSet *fResult{nullptr};      // current result set
   std::vector<oracle::occi::MetaData> *fFieldInfo{nullptr};   // metadata for columns
   Int_t                    fFieldCount{0};
   char                   **fFieldsBuffer{nullptr};

   Bool_t  IsValid(Int_t field);

   TOracleRow(const TOracleRow&) = delete;            // Not implemented.
   TOracleRow &operator=(const TOracleRow&)= delete; // Not implemented.

protected:
   void        GetRowData();

public:
   TOracleRow(oracle::occi::ResultSet *rs,
              std::vector<oracle::occi::MetaData> *fieldMetaData);
   ~TOracleRow();

   void        Close(Option_t *opt="") final;
   ULong_t     GetFieldLength(Int_t field) final;
   const char *GetField(Int_t field) final;

   ClassDefOverride(TOracleRow,0)  // One row of Oracle query result
};

#endif
