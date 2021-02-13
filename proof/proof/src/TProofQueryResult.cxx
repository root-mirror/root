// @(#)root/proof:$Id$
// Author: G Ganis Sep 2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TProofQueryResult
\ingroup proofkernel

TQueryResult version adapted to PROOF neeeds

*/

#include "TDSet.h"
#include "TError.h"
#include "TList.h"
#include "TProofQueryResult.h"
#include "TSystem.h"
#include "TMacro.h"

ClassImp(TProofQueryResult);

////////////////////////////////////////////////////////////////////////////////
/// Main constructor.

TProofQueryResult::TProofQueryResult(Int_t sn, const char *opt, TList *inlist,
                                     Long64_t ent, Long64_t fst, TDSet *dset,
                                     const char *sel, TObject *elist)
                  : TQueryResult(sn, opt, inlist, ent, fst, sel)
{
   fNumWrks = -1;
   fStartLog = -1;

   // Add data sets and event lists to the input list
   if (fInputList && dset)
      fInputList->Add(dset);
   if (fInputList && elist)
      fInputList->Add(elist);
}

////////////////////////////////////////////////////////////////////////////////
/// Call when running starts.

void TProofQueryResult::SetRunning(Int_t startlog, const char *par, Int_t nwrks)
{
   Info("SetRunning", "nwrks: %d", nwrks);

   fStatus = kRunning;
   fNumWrks = nwrks;
   fStart.Set();
   fEnd.Set(fStart.Convert()-1);
   fParList = (par && (strlen(par) > 0)) ? par : "-";
   fStartLog = startlog;

   // Add header to log file
   fLogFile->AddLine("+++");
   fLogFile->AddLine(Form("+++ Start processing query # %d (log file offset: %d)",
                     fSeqNum, startlog));
   fLogFile->AddLine("+++");
}
