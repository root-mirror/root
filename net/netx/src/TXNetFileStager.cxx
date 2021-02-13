// @(#)root/netx:$Id$
// Author: A. Peters, G. Ganis   7/2/2007

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TXNetFileStager                                                      //
//                                                                      //
// Interface to the 'XRD' staging capabilities.                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TError.h"
#include "TUrl.h"
#include "TXNetFileStager.h"
#include "TXNetSystem.h"
#include "TFileCollection.h"
#include "TStopwatch.h"
#include "TFileInfo.h"

////////////////////////////////////////////////////////////////////////////////
/// Constructor. Init a TXNetSystem instance to the XRD system.

TXNetFileStager::TXNetFileStager(const char *url) : TFileStager("xrd")
{
   fSystem = 0;
   if (url && strlen(url) > 0) {
      GetPrefix(url, fPrefix);

      fSystem = new TXNetSystem(fPrefix);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TXNetFileStager::~TXNetFileStager()
{
   if (fSystem)
      delete fSystem;
   fSystem = 0;
   fPrefix = "";
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the file defined by 'path' is ready to be used.

Bool_t TXNetFileStager::IsStaged(const char *path)
{
   if (!IsValid()) {
      GetPrefix(path, fPrefix);
      fSystem = new TXNetSystem(path);
   }

   if (IsValid()) {
      TString p(path);
      if (!p.BeginsWith("root:"))
         p.Insert(0, fPrefix);
      return (fSystem->IsOnline(p));
   }

   // Failure
   Warning("IsStaged","TXNetSystem not initialized");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Issue a stage request for file defined by 'path'. The string 'opt'
/// defines 'option' and 'priority' for 'Prepare': the format is
///    opt = "option=o priority=p".

Bool_t TXNetFileStager::Stage(TCollection *paths, Option_t *opt)
{
   if (IsValid()) {
      UChar_t o = 8;
      UChar_t p = 0;
      // Parse options, if any
      if (opt && strlen(opt) > 0) {
         TString xo(opt), io;
         Ssiz_t from = 0;
         while (xo.Tokenize(io, from, "[ ,|]")) {
            if (io.Contains("option=")) {
               io.ReplaceAll("option=","");
               if (io.IsDigit()) {
                  Int_t i = io.Atoi();
                  if (i >= 0 && i <= 255)
                     o = (UChar_t) i;
               }
            }
            if (io.Contains("priority=")) {
               io.ReplaceAll("priority=","");
               if (io.IsDigit()) {
                  Int_t i = io.Atoi();
                  if (i >= 0 && i <= 255)
                     p = (UChar_t) i;
               }
            }
         }
      }
      // Run prepare
      return fSystem->Prepare(paths, o, p);
   }

   // Failure
   Warning("Stage","TXNetSystem not initialized");
   return kFALSE;

}

////////////////////////////////////////////////////////////////////////////////
/// Issue a stage request for file defined by 'path'. The string 'opt'
/// defines 'option' and 'priority' for 'Prepare': the format is
///                opt = "option=o priority=p".

Bool_t TXNetFileStager::Stage(const char *path, Option_t *opt)
{
   if (!IsValid()) {
      GetPrefix(path, fPrefix);
      fSystem = new TXNetSystem(path);
   }

   if (IsValid()) {
      UChar_t o = 8;  // XrdProtocol.hh
      UChar_t p = 0;
      // Parse options
      TString xo(opt), io;
      Ssiz_t from = 0;
      while (xo.Tokenize(io, from, "[ ,|]")) {
         if (io.Contains("option=")) {
            io.ReplaceAll("option=","");
            if (io.IsDigit()) {
               Int_t i = io.Atoi();
               if (i >= 0 && i <= 255)
                  o = (UChar_t) i;
            }
         }
         if (io.Contains("priority=")) {
            io.ReplaceAll("priority=","");
            if (io.IsDigit()) {
               Int_t i = io.Atoi();
               if (i >= 0 && i <= 255)
                  p = (UChar_t) i;
            }
         }
      }
      // Make user the full path is used
      TString pp(path);
      if (!pp.BeginsWith("root:"))
         pp.Insert(0, fPrefix);
      return fSystem->Prepare(pp, o, p);
   }

   // Failure
   Warning("Stage","TXNetSystem not initialized");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Isolate prefix in url

void TXNetFileStager::GetPrefix(const char *url, TString &pfx)
{
   if (gDebug > 1)
      ::Info("TXNetFileStager::GetPrefix", "enter: %s", url);

   TUrl u(url);
   pfx = Form("%s://", u.GetProtocol());
   if (strlen(u.GetUser()) > 0)
      pfx += Form("%s@", u.GetUser());
   pfx += u.GetHost();
   if (u.GetPort() != TUrl("root://host").GetPort())
      pfx += Form(":%d", u.GetPort());
   pfx += "/";

   if (gDebug > 1)
      ::Info("TXNetFileStager::GetPrefix", "found prefix: %s", pfx.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Print basic info about this stager

void TXNetFileStager::Print(Option_t *) const
{
   Printf("+++ stager: %s  %s", GetName(), fPrefix.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Get actual end-point url for a path
/// Returns 0 in case of success and 1 if any error occured

Int_t TXNetFileStager::Locate(const char *path, TString &eurl)
{
   if (!IsValid()) {
      GetPrefix(path, fPrefix);
      fSystem = new TXNetSystem(path);
   }

   if (IsValid())
      return fSystem->Locate(path, eurl);

   // Unable to initialize TXNetSystem
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Bulk locate request for a collection of files. A noop prepare command is
/// issued beforehand to fill redirector's cache, then Locate() is issued on
/// each file. Results are saved back to the input collection: when a file is
/// found, the staged bit is set to on, and its endpoint URL is added, if
/// different from the redirector's URL. If a file is not found, the staged
/// bit is set to off.
/// If addDummyUrl is kTRUE, in case file is not staged or redirector is
/// identical to endpoint URL, a dummy URL is prepended, respectively:
/// "noop://redir" and "noop://none".
/// If the collection contains URLs with "anchors" (i.e., #fileName.root),
/// they are ignored by xrootd.
/// The Locate() command preserves anchors, but needs single paths to be full
/// URLs beginning with root://.
/// Returns < 0 in case of errors, and the number of files processed in case
/// of success.

Int_t TXNetFileStager::LocateCollection(TFileCollection *fc,
   Bool_t addDummyUrl)
{
   if (!fc) {
      Error("Locate", "No input collection given!");
      return -1;
   }

   // Fill redirector's cache with an empty prepare request
   //Int_t TXNetSystem::Prepare(TCollection *paths,
   //   UChar_t opt, UChar_t prio, TString *bufout)

   Int_t count = 0;

   TStopwatch ts;
   Double_t timeTaken_s;
   TFileInfo *fi;

   Int_t rv = fSystem->Prepare(fc->GetList(), 0, 0, NULL);
   //                                         o  p

   TIter it(fc->GetList());

   timeTaken_s = ts.RealTime();
   if (gDebug > 0) {
      Info("Locate", "Bulk xprep done in %.1lfs (returned %d)",
        ts.RealTime(), rv);
   }

   ts.Start();
   TString surl, endp;

   while ((fi = dynamic_cast<TFileInfo *>(it.Next())) != NULL) {

      surl = fi->GetCurrentUrl()->GetUrl();

      if (!IsValid()) {
         GetPrefix(surl.Data(), fPrefix);
         if (gDebug > 0) {
            Info("Locate", "Stager non initialized, doing it now for %s",
               fPrefix.Data());
         }
         fSystem = new TXNetSystem(surl.Data());
      }

      // Locating (0=success, 1=error -- 1 includes when file is not staged)
      if (fSystem->Locate(surl.Data(), endp)) {
         // File not staged
         fi->ResetBit(TFileInfo::kStaged);

         if (addDummyUrl)
            fi->AddUrl("noop://none", kTRUE);

         if (gDebug > 1)
            Info("Locate", "Not found: %s", surl.Data());
      }
      else {
         // File staged. Returned endpoint contains the same anchor and options.
         // We just check if it is equal to one of our current URLs.

         fi->SetBit(TFileInfo::kStaged);
         if (surl != endp) {
            fi->AddUrl(endp.Data(), kTRUE);
         }
         else if (addDummyUrl) {
            // Returned URL identical to redirector's URL
            fi->AddUrl("noop://redir", kTRUE);
         }

         if (gDebug > 1)
            Info("Locate", "Found: %s --> %s", surl.Data(), endp.Data());
      }

      count++;
   }

   timeTaken_s += ts.RealTime();
   if (gDebug > 0) {
      Info("Locate", "All locates finished in %.1lfs", ts.RealTime());
      Info("Locate", "Mass prepare and locates took %.1lfs", timeTaken_s);
   }

   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if stager 's' is compatible with current stager.
/// Avoids multiple instantiations of the potentially the same TXNetSystem.

Bool_t TXNetFileStager::Matches(const char *s)
{
   if (IsValid()) {
      TString pfx;
      GetPrefix(s, pfx);
      return ((fPrefix == pfx) ? kTRUE : kFALSE);
   }

   // Not valid
   return kFALSE;
}
