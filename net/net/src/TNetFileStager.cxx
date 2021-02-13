// @(#)root/net:$Id$
// Author: G. Ganis Feb 2011

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TNetFileStager                                                       //
//                                                                      //
// TFileStager implementation for a 'rootd' backend.                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TError.h"
#include "TUrl.h"
#include "TNetFile.h"
#include "TNetFileStager.h"

////////////////////////////////////////////////////////////////////////////////
/// Constructor. Init a TNetSystem instance to the remote rootd.

TNetFileStager::TNetFileStager(const char *url) : TFileStager("net")
{
   fSystem = 0;
   if (url && strlen(url) > 0) {
      GetPrefix(url, fPrefix);

      fSystem = new TNetSystem(fPrefix);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TNetFileStager::~TNetFileStager()
{
   SafeDelete(fSystem);
   fPrefix = "";
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the file defined by 'path' is ready to be used.

Bool_t TNetFileStager::IsStaged(const char *path)
{
   if (!IsValid()) {
      GetPrefix(path, fPrefix);
      fSystem = new TNetSystem(path);
   }

   if (IsValid()) {
      TString p(path);
      if (!p.BeginsWith(fPrefix)) p.Insert(0, fPrefix);
      return (fSystem->AccessPathName(p, kReadPermission) ? kFALSE : kTRUE);
   }

   // Failure
   Warning("IsStaged","TNetSystem not initialized");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Isolate prefix in url

void TNetFileStager::GetPrefix(const char *url, TString &pfx)
{
   if (gDebug > 1)
      ::Info("TNetFileStager::GetPrefix", "enter: %s", url);

   TUrl u(url);
   pfx = TString::Format("%s://", u.GetProtocol());
   if (strlen(u.GetUser()) > 0)
      pfx += TString::Format("%s@", u.GetUser());
   pfx += u.GetHost();
   if (u.GetPort() != TUrl("root://host").GetPort())
      pfx += TString::Format(":%d", u.GetPort());
   pfx += "/";

   if (gDebug > 1)
      ::Info("TNetFileStager::GetPrefix", "found prefix: %s", pfx.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Print basic info about this stager

void TNetFileStager::Print(Option_t *) const
{
   Printf("+++ stager: %s  %s", GetName(), fPrefix.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Get actual end-point url for a path
/// Returns 0 in case of success and 1 if any error occured

Int_t TNetFileStager::Locate(const char *path, TString &eurl)
{
   if (!IsValid()) {
      GetPrefix(path, fPrefix);
      fSystem = new TNetSystem(path);
   }

   if (IsValid()) {
      TString p(path);
      if (!p.BeginsWith(fPrefix)) p.Insert(0, fPrefix);
      if (!fSystem->AccessPathName(p, kReadPermission)) {
         eurl = p;
         return 0;
      }
   }

   // Unable to initialize TNetSystem or file does not exist
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if stager 's' is compatible with current stager.
/// Avoids multiple instantiations of the potentially the same TNetSystem.

Bool_t TNetFileStager::Matches(const char *s)
{
   if (IsValid()) {
      TString pfx;
      GetPrefix(s, pfx);
      return ((fPrefix == pfx) ? kTRUE : kFALSE);
   }

   // Not valid
   return kFALSE;
}
