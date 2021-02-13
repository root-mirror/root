// @(#)root/net:$Id$
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
// TFileStager                                                          //
//                                                                      //
// Abstract base class defining an interface to a stager.               //
//                                                                      //
// To open a connection to a stager use the static method               //
// Open("<stager>"), where <stager> contains a keyword allowing to load //
// the relevant plug-in, e.g.                                           //
//           TFileStager::Open("root://lxb6064.cern.ch")                //
// will load TXNetFileStager and initialize it for the redirector at    //
// lxb6046.cern.ch .                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TError.h"
#include "TFileInfo.h"
#include "TFile.h"
#include "TFileStager.h"
#include "TObjString.h"
#include "TPluginManager.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TUrl.h"
#include "TCollection.h"
#include "TFileCollection.h"
#include "THashList.h"

////////////////////////////////////////////////////////////////////////////////
/// Retrieves the staging (online) status for a list of path names. Path names
/// must be of type TUrl, TFileInfo or TObjString. The returned list is the list
/// of staged files as TObjString (we use TObjString, because you can do a FindObject
/// on that list using the file name, which is not possible with TUrl objects.

TList* TFileStager::GetStaged(TCollection *pathlist)
{
   if (!pathlist) {
      Error("GetStaged", "list of pathnames was not specified!");
      return 0;
   }

   TList* stagedlist = new TList();
   TIter nxt(pathlist);
   TObject* o = 0;
   Bool_t local = (!strcmp(GetName(), "local")) ? kTRUE : kFALSE;
   while ((o = nxt())) {
      TString pn = TFileStager::GetPathName(o);
      if (pn == "") {
         Warning("GetStaged", "object is of unexpected type %s - ignoring", o->ClassName());
      } else if (local || IsStaged(pn))
         stagedlist->Add(new TObjString(pn));
   }

   // List of online files
   stagedlist->SetOwner(kTRUE);
   Info("GetStaged", "%d files staged", stagedlist->GetSize());
   return stagedlist;
}

////////////////////////////////////////////////////////////////////////////////
/// Issue a stage request for a list of files.
/// Return the '&' of all single Prepare commands.

Bool_t TFileStager::Stage(TCollection *pathlist, Option_t *opt)
{
   TIter nxt(pathlist);
   TObject *o = 0;
   Bool_t success = kFALSE;

   while ((o = nxt()))  {
      TString pn = TFileStager::GetPathName(o);
      if (pn == "") {
         Warning("Stage", "found object of unexpected type %s - ignoring",
                              o->ClassName());
         continue;
      }

      // Issue to prepare
      success &= Stage(pn, opt);
   }

   // return global flag
   return success;
}

////////////////////////////////////////////////////////////////////////////////
/// Open a stager, after having loaded the relevant plug-in.
/// The format of 'stager' depends on the plug-in.

TFileStager *TFileStager::Open(const char *stager)
{
   TPluginHandler *h;
   TFileStager *s = 0;

   if (!stager) {
      ::Error("TFileStager::Open", "stager name missing: do nothing");
      return 0;
   }

   if (!gSystem->IsPathLocal(stager) &&
      (h = gROOT->GetPluginManager()->FindHandler("TFileStager", stager))) {
      if (h->LoadPlugin() == -1)
         return 0;
      s = (TFileStager *) h->ExecPlugin(1, stager);
   } else
      s = new TFileStager("local");

   return s;
}

////////////////////////////////////////////////////////////////////////////////
/// Just check if the file exists locally

Bool_t TFileStager::IsStaged(const char *f)
{
   // The safest is to open in raw mode
   TUrl u(f);
   u.SetOptions("filetype=raw");
   TFile *ff = TFile::Open(u.GetUrl());
   Bool_t rc = kTRUE;
   if (!ff || ff->IsZombie()) {
      rc = kFALSE;
   }
   if (ff) {
      ff->Close();
      delete ff;
   }
   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Just check if the file exists locally

Int_t TFileStager::Locate(const char *u, TString &f)
{
   if (!IsStaged(u))
      return -1;
   f = u;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Massive location of files. Returns < 0 on error, or number of files
/// processed. Results are returned on the TFileCollection itself

Int_t TFileStager::LocateCollection(TFileCollection *fc, Bool_t)
{
    TFileInfo *fi;
    TString endp;
    TIter it(fc->GetList());
    Int_t count = 0;

    while ((fi = dynamic_cast<TFileInfo *>(it.Next()))) {
       const char *ourl = fi->GetCurrentUrl()->GetUrl();
       if (!ourl) continue;

       if (Locate(ourl, endp) == 0) {
          fi->AddUrl(endp.Data(), kTRUE);
          fi->SetBit(TFileInfo::kStaged);
          fi->ResetUrl();
       }
       else {
          fi->ResetBit(TFileInfo::kStaged);
       }

       count++;
    }

   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the path name contained in object 'o' allowing for
/// TUrl, TObjString or TFileInfo

TString TFileStager::GetPathName(TObject *o)
{
   TString pathname;
   TString cn(o->ClassName());
   if (cn == "TUrl") {
      pathname = ((TUrl*)o)->GetUrl();
   } else if (cn == "TObjString") {
      pathname = ((TObjString*)o)->GetName();
   } else if (cn == "TFileInfo") {
      TFileInfo *fi = (TFileInfo *)o;
      pathname = (fi->GetCurrentUrl()) ? fi->GetCurrentUrl()->GetUrl() : "";
      if (fi->GetCurrentUrl()) {
         if (strlen(fi->GetCurrentUrl()->GetAnchor()) > 0) {
            TUrl url(*(fi->GetCurrentUrl()));
            url.SetAnchor("");
            pathname = url.GetUrl();
         }
      } else {
         pathname = fi->GetCurrentUrl()->GetUrl();
      }
   }

   // Done
   return pathname;
}
