// @(#)root/net:$Id$
// Author: Adrien Devresse and Fabrizio Furano

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TDavixSystem                                                         //
//                                                                      //
// A TSystem specialization for HTTP and WebDAV                         //
// It supports HTTP and HTTPS in a number of dialects and options       //
//  e.g. S3 is one of them                                              //
// Other caracteristics come from the full support of Davix,            //
//  e.g. full redirection support in any circumstance                   //
//                                                                      //
// Authors:     Adrien Devresse (CERN IT/SDC)                           //
//              Fabrizio Furano (CERN IT/SDC)                           //
//                                                                      //
// September 2013                                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TDavixSystem.h"
#include "TROOT.h"
#include "TSocket.h"
#include "Bytes.h"
#include "TError.h"
#include "TBase64.h"
#include "TVirtualPerfStats.h"
#include "TDavixFileInternal.h"

#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <davix.hpp>
#include <sstream>
#include <string>
#include <cstring>


extern const std::string VERSION;
extern const std::string gUserAgent;

// The prefix that is used to find the variables in the gEnv
#define ENVPFX "Davix."

ClassImp(TDavixSystem);

using namespace Davix;

extern const char* grid_mode_opt;
extern const char* ca_check_opt;
extern const char* s3_seckey_opt;
extern const char* s3_acckey_opt;

////////////////////////////////////////////////////////////////////////////////

TDavixSystem::TDavixSystem(const char *url) :
   TSystem(url),
   d_ptr(new TDavixFileInternal(url, "WEB"))
{
   d_ptr->init();
   SetTitle("WebDAV system administration");
}

////////////////////////////////////////////////////////////////////////////////

TDavixSystem::TDavixSystem() :
   TSystem(),
   d_ptr(new TDavixFileInternal("", "WEB"))
{
   d_ptr->init();
   SetTitle("WebDAV system administration");
}

////////////////////////////////////////////////////////////////////////////////

TDavixSystem::~TDavixSystem()
{
   SafeDelete(d_ptr);
}

////////////////////////////////////////////////////////////////////////////////

void TDavixSystem::FreeDirectory(void *dirp)
{
   d_ptr->davixPosix->closedir(static_cast<DAVIX_DIR *>(dirp), NULL);
   d_ptr->removeDird(dirp);
}

////////////////////////////////////////////////////////////////////////////////

const char *TDavixSystem::GetDirEntry(void *dirp)
{
   struct dirent *dir;
   DavixError *davixErr = NULL;
   if (((dir = d_ptr->davixPosix->readdir(static_cast<DAVIX_DIR *>(dirp), &davixErr)) == NULL)
         && (davixErr != NULL)) {
      Error("DavixReaddir", "failed to readdir the directory: %s (%d)",
            davixErr->getErrMsg().c_str(), davixErr->getStatus());
      DavixError::clearError(&davixErr);
   }
   return (dir) ? (dir->d_name) : NULL;
}

////////////////////////////////////////////////////////////////////////////////

void *TDavixSystem::OpenDirectory(const char *dir)
{
   DavixError *davixErr = NULL;
   DAVIX_DIR *d;
   if ((d = d_ptr->davixPosix->opendir(d_ptr->davixParam, dir, &davixErr)) == NULL) {
      Error("DavixOpendir", "failed to opendir the directory: %s (%d)",
            davixErr->getErrMsg().c_str(), davixErr->getStatus());
      DavixError::clearError(&davixErr);
   } else {
      d_ptr->addDird(d);
   }
   return d;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TDavixSystem::ConsistentWith(const char * /*path*/, void *dirptr)
{
   return (Bool_t) d_ptr->isMyDird(dirptr);
}

////////////////////////////////////////////////////////////////////////////////

Int_t TDavixSystem::GetPathInfo(const char *path, FileStat_t &buf)
{
   struct stat st;

   if (!d_ptr->DavixStat(path, &st)) return 1;
   buf.fDev = 0;
   buf.fIno = 0;
   buf.fMode = st.st_mode; // protection (combination of EFileModeMask bits)

   buf.fUid = st.st_uid; // user id of owner
   buf.fGid = st.st_gid; // group id of owner
   buf.fSize = st.st_size; // total size in bytes
   buf.fMtime = st.st_mtime; // modification date
   buf.fIsLink = kFALSE; // symbolic link
   buf.fUrl = path; // end point url of file

   return 0;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TDavixSystem::IsPathLocal(const char *path)
{
   (void) path;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

Int_t TDavixSystem::Locate(const char *path, TString &endurl)
{
   DavixError *davixErr = NULL;
   ssize_t ret;
   ReplicaVec vecRep;
   DavFile f(*d_ptr->davixContext, Uri(path));
   if ((ret = f.getAllReplicas(d_ptr->davixParam,
                               vecRep,
                               &davixErr)) < 0) {
      Error("DavixLocate", "failed to Locate file: %s (%d)",
            davixErr->getErrMsg().c_str(), davixErr->getStatus());
      DavixError::clearError(&davixErr);
      return 1;
   }
   if (vecRep.size() > 0) {
      endurl = vecRep[0].uri.getString().c_str();
   } else {
      endurl = path;
   }
   if (gDebug > 0)
      Info("DavixLocate", "Davix Locate %s to %s", path, endurl.Data());

   return 0;
}

////////////////////////////////////////////////////////////////////////////////

Int_t TDavixSystem::MakeDirectory(const char *dir)
{
   DavixError *davixErr = NULL;
   int ret;
   if ((ret = d_ptr->davixPosix->mkdir(d_ptr->davixParam, dir, 0755, &davixErr))  < 0) {
      Error("DavixMkdir", "failed to create the directory: %s (%d)",
            davixErr->getErrMsg().c_str(), davixErr->getStatus());
      DavixError::clearError(&davixErr);
   }
   return ret;
}

////////////////////////////////////////////////////////////////////////////////

int TDavixSystem::Unlink(const char *path)
{
   DavixError *davixErr = NULL;
   int ret;
   if ((ret = d_ptr->davixPosix->unlink(d_ptr->davixParam, path, &davixErr))  < 0) {
      Error("DavixUnlink", "failed to unlink the file: %s (%d)",
            davixErr->getErrMsg().c_str(), davixErr->getStatus());
      DavixError::clearError(&davixErr);
   }
   return ret;
}
