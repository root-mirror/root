// @(#)root/tree:$Id$
// Author: G Ganis Sep 2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TQueryResult
\ingroup tree

A container class for query results.
*/

#include <cstring>

#include "strlcpy.h"
#include "TBrowser.h"
#include "TEventList.h"
#include "TQueryResult.h"
#include "TRegexp.h"
#include "TROOT.h"
#include "TMacro.h"
#include "TMath.h"
#include "TSelector.h"
#include "TSystem.h"
#include "TTimeStamp.h"

ClassImp(TQueryResult);

////////////////////////////////////////////////////////////////////////////////
/// Main constructor.

TQueryResult::TQueryResult(Int_t seqnum, const char *opt, TList *inlist,
                           Long64_t entries, Long64_t first, const char *selec)
             : fSeqNum(seqnum), fStatus(kSubmitted), fUsedCPU(0.), fOptions(opt),
               fEntries(entries), fFirst(first),
               fBytes(0), fParList("-"), fOutputList(0),
               fFinalized(kFALSE), fArchived(kFALSE), fResultFile("-"),
               fPrepTime(0.), fInitTime(0.), fProcTime(0.), fMergeTime(0.),
               fRecvTime(-1), fTermTime(-1), fNumWrks(-1), fNumMergers(-1)
{
   // Name and unique title
   SetName(TString::Format("q%d", fSeqNum));
   SetTitle(TString::Format("session-localhost-%ld-%d",
                 (Long_t)TTimeStamp().GetSec(), gSystem->GetPid()));

   // Start time
   fStart.Set();
   fEnd.Set(fStart.Convert()-1);

   // Save input list
   fInputList = 0;
   if (inlist) {
      fInputList = (TList *) (inlist->Clone());
      fInputList->SetOwner();
   }

   // Log file
   fLogFile = new TMacro("LogFile");

   // Selector files
   fDraw = selec ? TSelector::IsStandardDraw(selec) : kFALSE;
   if (fDraw) {
      // The input list should contain info about the variables and
      // selection cuts: save them into the macro title
      TString varsel;
      if (fInputList) {
         TIter nxo(fInputList);
         TObject *o = 0;
         while ((o = nxo())) {
            if (!strcmp(o->GetName(),"varexp")) {
               varsel = o->GetTitle();
               Int_t iht = varsel.Index(">>htemp");
               if (iht > -1)
                  varsel.Remove(iht);
               varsel.Form("\"%s\";", varsel.Data());
            }
            if (!strcmp(o->GetName(),"selection"))
               varsel += TString::Format("\"%s\"", o->GetTitle());
         }
         if (gDebug > 0)
            Info("TQueryResult","selec: %s, varsel: %s", selec, varsel.Data());
         // Log notification also in the instance
         fLogFile->AddLine(TString::Format("TQueryResult: selec: %s, varsel: %s",
                                           selec, varsel.Data()));
      }
      // Standard draw action: save only the name
      fSelecImp = new TMacro(selec, varsel);
      fSelecHdr = 0;
   } else {
      // Save selector file
      fSelecHdr = new TMacro;
      fSelecImp = new TMacro;
      SaveSelector(selec);
   }

   // List of libraries loaded at creation
   const char *pl = gSystem->GetLibraries();
   fLibList = (pl && (strlen(pl) > 0)) ? pl : "-";
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TQueryResult::~TQueryResult()
{
   SafeDelete(fInputList);
   SafeDelete(fOutputList);
   SafeDelete(fLogFile);
   SafeDelete(fSelecImp);
   SafeDelete(fSelecHdr);
}

////////////////////////////////////////////////////////////////////////////////
/// Return an instance of TQueryResult containing only the local
/// info fields, i.e. no outputlist, liblist, dset, selectors, etc..
/// Used for fast retrieve of information about existing queries
/// and their status.

TQueryResult *TQueryResult::CloneInfo()
{
   // Create instance
   TQueryResult *qr = new TQueryResult(fSeqNum, fOptions, 0, fEntries,
                                       fFirst, 0);

   // Correct fields
   qr->fStatus = fStatus;
   qr->fStart.Set(fStart.Convert());
   qr->fEnd.Set(fEnd.Convert());
   qr->fUsedCPU = fUsedCPU;
   qr->fEntries = fEntries;
   qr->fFirst = fFirst;
   qr->fBytes = fBytes;
   qr->fParList = fParList;
   qr->fResultFile = fResultFile;
   qr->fArchived = fArchived;
   qr->fPrepTime = fPrepTime;
   qr->fInitTime = fInitTime;
   qr->fProcTime = fProcTime;
   qr->fMergeTime = fMergeTime;
   qr->fRecvTime = fRecvTime;
   qr->fTermTime = fTermTime;
   qr->fNumWrks = fNumWrks;
   qr->fNumMergers = fNumMergers;

   qr->fSelecHdr = 0;
   if (GetSelecHdr()) {
      qr->fSelecHdr = new TMacro();
      qr->fSelecHdr->SetName(GetSelecHdr()->GetName());
      qr->fSelecHdr->SetTitle(GetSelecHdr()->GetTitle());
   }
   qr->fSelecImp = 0;
   if (GetSelecImp()) {
      qr->fSelecImp = new TMacro();
      qr->fSelecImp->SetName(GetSelecImp()->GetName());
      qr->fSelecImp->SetTitle(GetSelecImp()->GetTitle());
   }

   // Name and title
   qr->SetName(GetName());
   qr->SetTitle(GetTitle());

   return qr;
}

////////////////////////////////////////////////////////////////////////////////
/// Save the selector header and implementation into the dedicated
/// TMacro instances. The header is searched for in the same directory
/// of the implementation file.

void TQueryResult::SaveSelector(const char *selector)
{
   if (!selector)
      return;

   // Separate out aclic chars
   TString selec = selector;
   TString aclicMode;
   TString arguments;
   TString io;
   selec = gSystem->SplitAclicMode(selec, aclicMode, arguments, io);

   // Store aclic options, if any
   if (aclicMode.Length() > 0)
      fOptions += TString::Format("#%s", aclicMode.Data());

   // If the selector is in a precompiled shared lib (e.g. in a PAR)
   // we just save the name
   TString selname = gSystem->BaseName(selec);
   fSelecImp->SetName(selname);
   Int_t idx = selname.Index(".");
   if (idx < 0) {
      // Notify
      if (gDebug > 0)
         Info("SaveSelector", "precompiled selector: just save the name");
      fSelecImp->SetTitle(selname);
   } else {
      // We locate the file and save it in compressed form
      if (idx > -1)
         selname.Remove(idx);
      fSelecImp->SetTitle(selname);

      // Locate the implementation file
      char *selc = gSystem->Which(TROOT::GetMacroPath(), selec, kReadPermission);
      if (!selc) {
         if (gDebug > 0)
            Warning("SaveSelector",
                    "could not locate selector implementation file (%s)", selec.Data());
         return;
      }

      // Fill the TMacro instance
      fSelecImp->ReadFile(selc);
      fSelecImp->SetName(gSystem->BaseName(selc));

      // Locate the included header file
      char *p = (char *) strrchr(selc,'.');
      if (p) {
         strlcpy(p+1,"h",strlen(p));
      } else {
         if (gDebug > 0)
            Warning("SaveSelector",
                    "bad formatted name (%s): could not build header file name", selc);
      }
      if (!(gSystem->AccessPathName(selc, kReadPermission))) {
         fSelecHdr->ReadFile(selc);
         fSelecHdr->SetName(gSystem->BaseName(selc));
         fSelecHdr->SetTitle(selname);
      } else {
         if (gDebug > 0)
            Warning("SaveSelector",
                    "could not locate selector header file (%s)", selc);
      }

      delete[] selc;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// End of query settings.

void TQueryResult::RecordEnd(EQueryStatus status, TList *outlist)
{
   // End time
   fEnd.Set();

   // Status
   fStatus = (status < kAborted || status > kCompleted) ? kAborted : status;

   // Clone the results
   if (outlist && fOutputList != outlist) {
      if (fOutputList) {
         fOutputList->Delete();
         SafeDelete(fOutputList);
      }
      if ((fOutputList = (TList *) (outlist->Clone()))) {
         fOutputList->SetOwner();
         Info("RecordEnd", "output list cloned successfully!");
      } else {
         Warning("RecordEnd", "unable to clone output list!!!");
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set processing info.

void TQueryResult::SetProcessInfo(Long64_t ent, Float_t cpu, Long64_t bytes,
                                  Float_t init, Float_t proc)
{
   fEntries = (ent > 0) ? ent : fEntries;
   fUsedCPU = (cpu > 0.) ? cpu : fUsedCPU;
   fBytes = (bytes > 0.) ? bytes : fBytes;
   fInitTime = (init > 0.) ? init : fInitTime;
   fProcTime = (proc > 0.) ? proc : fProcTime;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill log file.

void TQueryResult::AddLogLine(const char *logline)
{
   if (logline)
      fLogFile->AddLine(logline);
}

////////////////////////////////////////////////////////////////////////////////
/// Add obj to the input list

void TQueryResult::AddInput(TObject *obj)
{
   if (fInputList && obj)
      fInputList->Add(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Set (or update) query in archived state.

void TQueryResult::SetArchived(const char *archfile)
{
   if (IsDone()) {
      fArchived = kTRUE;
      if (archfile && (strlen(archfile) > 0))
         fResultFile = archfile;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print query content. Use opt = "F" for a full listing.

void TQueryResult::Print(Option_t *opt) const
{
   // Attention: the list must match EQueryStatus
   const char *qst[] = {
      "aborted  ", "submitted", "running  ", "stopped  ", "completed"
   };

   // Status label
   Int_t st = (fStatus > 0 && fStatus <= kCompleted) ? fStatus : 0;

   // Range label
   Long64_t last = (fEntries > -1) ? fFirst+fEntries-1 : -1;

   // Option
   Bool_t full = ((strchr(opt,'F') || strchr(opt,'f'))) ? kTRUE : kFALSE;

   // Query number to be printed
   Int_t qry = fSeqNum;
   TString qn = opt;
   TRegexp re("N.*N");
   Int_t i1 = qn.Index(re);
   if (i1 != kNPOS) {
      qn.Remove(0, i1+1);
      qn.Remove(qn.Index("N"));
      qry = qn.Atoi();
   }

   // Print separator if full dump
   if (full) Printf("+++");

   TString range;
   if (!full && (last > -1))
      range.Form("evts:%lld-%lld", fFirst, last);

   // Print header
   if (!fDraw) {
      const char *fin = fFinalized ? "finalized" : qst[st];
      const char *arc = fArchived ? "(A)" : "";
      Printf("+++ #:%d ref:\"%s:%s\" sel:%s %9s%s %s",
             qry, GetTitle(), GetName(), fSelecImp->GetTitle(), fin, arc,
             range.Data());
   } else {
      Printf("+++ #:%d ref:\"%s:%s\" varsel:%s %s",
             qry, GetTitle(), GetName(), fSelecImp->GetTitle(),
             range.Data());
   }

   // We are done, if not full dump
   if (!full) return;

   // Time information
   Float_t elapsed = (fProcTime > 0.) ? fProcTime
                                      : (Float_t)(fEnd.Convert() - fStart.Convert());
   Printf("+++        started:   %s", fStart.AsString());
   if (fPrepTime > 0.)
      Printf("+++        prepare:   %.3f sec", fPrepTime);
   Printf("+++        init:      %.3f sec", fInitTime);
   Printf("+++        process:   %.3f sec (CPU time: %.1f sec)", elapsed, fUsedCPU);
   if (fNumMergers > 0) {
      Printf("+++        merge:     %.3f sec (%d mergers)", fMergeTime, fNumMergers);
   } else {
      Printf("+++        merge:     %.3f sec ", fMergeTime);
   }
   if (fRecvTime > 0.)
      Printf("+++        transfer:  %.3f sec", fRecvTime);
   if (fTermTime > 0.)
      Printf("+++        terminate: %.3f sec", fTermTime);

   // Number of events processed, rate, size
   Double_t rate = 0.0;
   if (fEntries > -1 && elapsed > 0)
      rate = fEntries / (Double_t)elapsed ;
   Float_t size = ((Float_t)fBytes) / TMath::Power(2.,20.);
   Printf("+++        processed: %lld events (size: %.3f MBs)", fEntries, size);
   Printf("+++        rate:      %.1f evts/sec", rate);

   Printf("+++        # workers: %d ", fNumWrks);

   // Package information
   if (fParList.Length() > 1)
      Printf("+++        packages:  %s", fParList.Data());

   // Result information
   TString res = fResultFile;
   if (!fArchived) {
      Int_t dq = res.Index("queries");
      if (dq > -1) {
         res.Remove(0,res.Index("queries"));
         res.Insert(0,"<PROOF_SandBox>/");
      }
      if (res.BeginsWith("-")) {
         res = (fStatus == kAborted) ? "not available" : "sent to client";
      }
   }
   if (res.Length() > 1)
      Printf("+++        results:   %s", res.Data());

   if (fOutputList && fOutputList->GetSize() > 0)
      Printf("+++        outlist:   %d objects", fOutputList->GetSize());
}

////////////////////////////////////////////////////////////////////////////////
/// To support browsing of the results.

void TQueryResult::Browse(TBrowser *b)
{
   if (fOutputList)
      b->Add(fOutputList, fOutputList->Class(), "OutputList");
}

////////////////////////////////////////////////////////////////////////////////
/// Set / change the input list.
/// The flag 'adopt' determines whether the list is adopted (default)
/// or cloned. If adopted, object ownership is transferred to this object.
/// The internal fInputList will always be owner of its objects.

void TQueryResult::SetInputList(TList *in, Bool_t adopt)
{
   if (!in || in != fInputList)
      SafeDelete(fInputList);

   if (in && in != fInputList) {
      if (!adopt) {
         fInputList = (TList *) (in->Clone());
      } else {
         fInputList = new TList;
         TIter nxi(in);
         TObject *o = 0;
         while ((o = nxi()))
            fInputList->Add(o);
         in->SetOwner(kFALSE);
      }
      fInputList->SetOwner();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set / change the output list.
/// The flag 'adopt' determines whether the list is adopted (default)
/// or cloned.  If adopted, object ownership is transferred to this object.
/// The internal fOutputList will always be owner of its objects.

void TQueryResult::SetOutputList(TList *out, Bool_t adopt)
{
   if (!out) {
      SafeDelete(fOutputList);
      return;
   }

   if (out && out != fOutputList) {
      TObject *o = 0;
      if (fOutputList) {
         TIter nxoo(fOutputList);
         while ((o = nxoo())) {
            if (out->FindObject(o)) fOutputList->Remove(o);
         }
         SafeDelete(fOutputList);
      }
      if (!adopt) {
         fOutputList = (TList *) (out->Clone());
      } else {
         fOutputList = new TList;
         TIter nxo(out);
         o = 0;
         while ((o = nxo()))
            fOutputList->Add(o);
         out->SetOwner(kFALSE);
      }
      fOutputList->SetOwner();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compare two query result instances for equality.
/// Session name and query number are compared.

Bool_t operator==(const TQueryResult &qr1, const TQueryResult &qr2)
{
   if (!strcmp(qr1.GetTitle(), qr2.GetTitle()))
      if (qr1.GetSeqNum() == qr2.GetSeqNum())
         return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return TRUE if reference ref matches.

Bool_t TQueryResult::Matches(const char *ref)
{
   TString lref; lref.Form("%s:%s", GetTitle(), GetName());

   if (lref == ref)
      return kTRUE;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return first instance of class 'classname' in the input list.
/// Usefull to access TDSet, TEventList, ...

TObject *TQueryResult::GetInputObject(const char *classname) const
{
   TObject *o = 0;
   if (classname && fInputList) {
      TIter nxi(fInputList);
      while ((o = nxi()))
         if (!strncmp(o->ClassName(), classname, strlen(classname)))
            return o;
   }

   // Not found
   return o;
}
