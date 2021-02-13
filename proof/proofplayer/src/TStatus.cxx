// @(#)root/proofplayer:$Id$
// Author: Maarten Ballintijn   7/06/2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TStatus
\ingroup proofkernel

This class holds the status of an ongoing operation and collects
error messages. It provides a Merge() operation allowing it to
be used in PROOF to monitor status in the slaves.
No messages indicates success.

*/

#include "TStatus.h"
#include "Riostream.h"
#include "TBuffer.h"
#include "TClass.h"
#include "TObjString.h"
#include "TProofDebug.h"

ClassImp(TStatus);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TStatus::TStatus() : fIter(&fMsgs), fExitStatus(-1),
                     fVirtMemMax(-1), fResMemMax(-1),
                     fVirtMaxMst(-1), fResMaxMst(-1)
{
   SetName("PROOF_Status");
   fMsgs.SetOwner(kTRUE);
   fInfoMsgs.SetOwner(kTRUE);
   ResetBit(TStatus::kNotOk);
}

////////////////////////////////////////////////////////////////////////////////
/// Add an error message.

void TStatus::Add(const char *mesg)
{
   fMsgs.Add(new TObjString(mesg));
   SetBit(TStatus::kNotOk);
   Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Add an info message.

void TStatus::AddInfo(const char *mesg)
{
   fInfoMsgs.Add(new TObjString(mesg));
}

////////////////////////////////////////////////////////////////////////////////
/// PROOF Merge() function.

Int_t TStatus::Merge(TCollection *li)
{
   TIter stats(li);
   PDB(kOutput,1)
      Info("Merge", "start: max virtual memory: %.2f MB \tmax resident memory: %.2f MB ",
                    GetVirtMemMax()/1024., GetResMemMax()/1024.);
   while (TObject *obj = stats()) {
      TStatus *s = dynamic_cast<TStatus*>(obj);
      if (s == 0) continue;

      TObjString *os = 0;
      // Errors
      TIter nxem(&(s->fMsgs));
      while ((os = (TObjString *) nxem())) {
         Add(os->GetName());
      }

      // Infos (no duplications)
      TIter nxwm(&(s->fInfoMsgs));
      while ((os = (TObjString *) nxwm())) {
         if (!fInfoMsgs.FindObject(os->GetName()))
            AddInfo(os->GetName());
      }

      SetMemValues(s->GetVirtMemMax(), s->GetResMemMax());
      // Check the master values (relevantt if merging submaster info)
      SetMemValues(s->GetVirtMemMax(kTRUE), s->GetResMemMax(kTRUE), kTRUE);
      PDB(kOutput,1)
         Info("Merge", "during: max virtual memory: %.2f MB \t"
                       "max resident memory: %.2f MB ",
                       GetVirtMemMax()/1024., GetResMemMax()/1024.);
      if (GetVirtMemMax(kTRUE) > 0) {
         PDB(kOutput,1)
            Info("Merge", "during: max master virtual memory: %.2f MB \t"
                        "max master resident memory: %.2f MB ",
                        GetVirtMemMax(kTRUE)/1024., GetResMemMax(kTRUE)/1024.);
      }
   }

   return fMsgs.GetSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Standard print function.

void TStatus::Print(Option_t * /*option*/) const
{
   Printf("OBJ: %s\t%s\t%s", IsA()->GetName(), GetName(), (IsOk() ? "OK" : "ERROR"));

   TObjString *os = 0;
   // Errors first
   if (fMsgs.GetSize() > 0) {
      Printf("\n   Errors:");
      TIter nxem(&fMsgs);
      while ((os = (TObjString *) nxem()))
         Printf("\t%s",os->GetName());
      Printf(" ");
   }

   // Infos
   if (fInfoMsgs.GetSize() > 0) {
      Printf("\n   Infos:");
      TIter nxem(&fInfoMsgs);
      while ((os = (TObjString *) nxem()))
         Printf("\t%s",os->GetName());
      Printf(" ");
   }

   Printf(" Max worker virtual memory: %.2f MB \tMax worker resident memory: %.2f MB ",
          GetVirtMemMax()/1024., GetResMemMax()/1024.);
   Printf(" Max master virtual memory: %.2f MB \tMax master resident memory: %.2f MB ",
          GetVirtMemMax(kTRUE)/1024., GetResMemMax(kTRUE)/1024.);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the iterator on the messages.

void TStatus::Reset()
{
   fIter.Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the next message or 0.

const char *TStatus::NextMesg()
{
   TObjString *os = (TObjString *) fIter();
   if (os) return os->GetName();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set max memory values

void TStatus::SetMemValues(Long_t vmem, Long_t rmem, Bool_t master)
{
   if (master) {
      if (vmem > 0. && (fVirtMaxMst < 0. || vmem > fVirtMaxMst)) fVirtMaxMst = vmem;
      if (rmem > 0. && (fResMaxMst < 0. || rmem > fResMaxMst)) fResMaxMst = rmem;
   } else {
      if (vmem > 0. && (fVirtMemMax < 0. || vmem > fVirtMemMax)) fVirtMemMax = vmem;
      if (rmem > 0. && (fResMemMax < 0. || rmem > fResMemMax)) fResMemMax = rmem;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TStatus.

void TStatus::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 4) {
         R__b.ReadClassBuffer(TStatus::Class(), this, R__v, R__s, R__c);
      } else {
         // For version <= 4 masters we need a special streamer
         TNamed::Streamer(R__b);
         std::set<std::string> msgs;
         TClass *cl = TClass::GetClass("set<string>");
         if (cl) {
            UInt_t SS__s = 0, SS__c = 0;
            UInt_t SS__v = cl->GetClassVersion();
            R__b.ReadClassBuffer(cl, &msgs, SS__v, SS__s, SS__c);
         } else {
            Error("Streamer", "no info found for 'set<string>' - skip");
            return;
         }
         std::set<std::string>::const_iterator it;
         for (it = msgs.begin(); it != msgs.end(); ++it) {
            fMsgs.Add(new TObjString((*it).c_str()));
         }
         if (R__v > 2) {
            R__b >> fExitStatus;
         }
         if (R__v > 1) {
            R__b >> fVirtMemMax;
            R__b >> fResMemMax;
         }
         if (R__v > 3) {
            R__b >> fVirtMaxMst;
            R__b >> fResMaxMst;
         }
      }
   } else {
      R__b.WriteClassBuffer(TStatus::Class(),this);
   }
}


