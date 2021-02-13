// @(#)root/tree:$Id$
// Author: Philippe Canal and al. 08/2004

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef TBASKETSQL_CXX
#define TBASKETSQL_CXX

#include "TBasketSQL.h"
#include "TBranch.h"
#include "TFile.h"
#include "TTreeSQL.h"
#include "TBufferSQL.h"

#include <vector>

ClassImp(TBasketSQL);

/** \class TBasketSQL
\ingroup tree

Implement TBasket for a SQL backend.
*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TBasketSQL::TBasketSQL() : TBasket(), fResultPtr(0), fRowPtr(0), fInsertQuery(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Regular constructor.

TBasketSQL::TBasketSQL(const char *name, const char *title, TBranch *branch,
                         TSQLResult ** rs, TString *insert_query,
                         std::vector<Int_t> *vc, TSQLRow **r) :
  fResultPtr(rs),fRowPtr(r)
{
   SetName(name);
   SetTitle(title);
   fClassName   = "TBasketSQL";
   fBufferSize  = branch->GetBasketSize();
   fNevBufSize  = branch->GetEntryOffsetLen();
   fNevBuf      = 0;
   fEntryOffset = 0;  //Must be set to 0 before calling Sizeof
   fDisplacement= 0;  //Must be set to 0 before calling Sizeof
   fBuffer      = 0;  //Must be set to 0 before calling Sizeof
   fInsertQuery = insert_query;

   if (vc==0) {
      fBufferRef = 0;
   } else {
      fBufferRef = new TBufferSQL(TBuffer::kWrite, fBufferSize, vc, fInsertQuery, fRowPtr);
   }
   fHeaderOnly  = kTRUE;
   fLast        = 0; // Must initialize before calling Streamer()
   //Streamer(*fBufferRef);
   fBuffer      = 0;
   fBranch      = branch;
   fHeaderOnly  = kFALSE;
   branch->GetTree()->IncrementTotalBuffers(fBufferSize);

}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TBasketSQL::~TBasketSQL()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TSQLBuffer for this basket.

void TBasketSQL::CreateBuffer(const char *name, TString title,
                              std::vector<Int_t> *vc,
                              TBranch *branch, TSQLResult ** rs)
{
   fResultPtr = rs;
   SetName(name);
   SetTitle(title);
   fClassName   = "TBasketSQL";
   fBufferSize  = branch->GetBasketSize();
   fNevBufSize  = branch->GetEntryOffsetLen();
   fNevBuf      = 0;
   fEntryOffset = 0;  //Must be set to 0 before calling Sizeof
   fDisplacement= 0;  //Must be set to 0 before calling Sizeof
   fBuffer      = 0;  //Must be set to 0 before calling Sizeof

   if(vc==0) {
      fBufferRef = 0;
      Error("CreateBuffer","Need a vector of columns\n");
   } else {
      fBufferRef   = new TBufferSQL(TBuffer::kWrite, fBufferSize, vc, fInsertQuery, fRowPtr);
   }
   fHeaderOnly  = kTRUE;
   fLast        = 0;
   //Streamer(*fBufferRef);
   fBuffer      = 0;
   fBranch      = branch;
   fHeaderOnly  = kFALSE;
   branch->GetTree()->IncrementTotalBuffers(fBufferSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Prepare the basket for the next entry.

void TBasketSQL::PrepareBasket(Long64_t entry)
{
   ((TBufferSQL*)fBufferRef)->ResetOffset();
   ((TTreeSQL*)fBranch->GetTree())->PrepEntry(entry);
   fBufferRef->Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// See TBasket::ReadBasketBytes.  This is not implemented in TBasketSQL.

Int_t TBasketSQL::ReadBasketBytes(Long64_t , TFile *)
{
   Error("ReadBasketBytes","This member function should not be called!");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// See TBasket::ReadBasketBuffers.  This is not implemented in TBasketSQL.

Int_t TBasketSQL::ReadBasketBuffers(Long64_t , Int_t, TFile *)
{
   Error("ReadBasketBuffers","This member function should not be called!");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// See TBasket::Update.

void TBasketSQL::Update(Int_t, Int_t)
{
   ((TBufferSQL*)fBufferRef)->ResetOffset();
   fNevBuf++;
}

#endif
