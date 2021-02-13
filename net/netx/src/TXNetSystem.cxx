// @(#)root/netx:$Id$
// Author: Frank Winklmeier, Fabrizio Furano

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TXNetSystem                                                          //
//                                                                      //
// Authors: Frank Winklmeier,  Fabrizio Furano                          //
//          INFN Padova, 2005                                           //
//                                                                      //
// TXNetSystem is an extension of TNetSystem able to deal with new      //
// xrootd servers. The class detects the nature of the server and       //
// redirects the calls to TNetSystem in case of a rootd server.         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TEnv.h"
#include "TFileStager.h"
#include "TROOT.h"
#include "TSocket.h"
#include "TString.h"
#include "TUrl.h"
#include "TVirtualMutex.h"
#include "TXNetFile.h"
#include "TXNetSystem.h"

#include "XrdOuc/XrdOucString.hh"
#include "XrdClient/XrdClientVector.hh"
#include "XrdClient/XrdClientAdmin.hh"
#include "XrdClient/XrdClientConn.hh"
#include "XrdClient/XrdClientConst.hh"
#include "XrdClient/XrdClientEnv.hh"
#include "XProtocol/XProtocol.hh"

#include "XrdProofdXrdVers.h"


////////////////////////////////////////////////////////////////////////////////
/// PluginManager creation function

TSystem* ROOT_Plugin_TXNetSystem(const char *url, Bool_t owner) {
   return new TXNetSystem(url, owner);
}


ClassImp(TXNetSystem);

Bool_t TXNetSystem::fgInitDone = kFALSE;
Bool_t TXNetSystem::fgRootdBC = kTRUE;
THashList TXNetSystem::fgAddrFQDN;
THashList TXNetSystem::fgAdminHash;

typedef XrdClientVector<XrdOucString> VecString_t;

////////////////////////////////////////////////////////////////////////////////
/// Create system management class without connecting to server.

TXNetSystem::TXNetSystem(Bool_t owner) : TNetSystem(owner), fDirList(0)
{
   SetTitle("(x)rootd system administration");
   fIsRootd = kFALSE;
   fIsXRootd = kFALSE;
   fDir = "";
   fDirp = 0;
   fUrl = "";
}

////////////////////////////////////////////////////////////////////////////////
/// Create system management class and connect to server specified by url.

TXNetSystem::TXNetSystem(const char *url, Bool_t owner) : TNetSystem(owner), fDirList(0)
{
   SetTitle("(x)rootd system administration");
   fIsRootd = kFALSE;
   fIsXRootd = kFALSE;
   fDir = "";
   fDirp = 0;
   fUrl = url;

   fgAddrFQDN.SetOwner();
   fgAdminHash.SetOwner();

   // Set debug level
   EnvPutInt(NAME_DEBUG, gEnv->GetValue("XNet.Debug", -1));

   // The first time do some global initialization
   if (!fgInitDone)
      InitXrdClient();

   // Fill in user, host, port
   TNetSystem::InitRemoteEntity(url);

   TXNetSystemConnectGuard cguard(this, url);
   if (!cguard.IsValid() && !fIsRootd)
      Error("TXNetSystem","fatal error: connection creation failed.");

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Init a connection to the server.
/// Returns a pointer to the appropriate instance of XrdClientAdmin or 0
/// in case of failure.

XrdClientAdmin *TXNetSystem::Connect(const char *url)
{
   // We need a dummy filename after the server url to connect
   TString dummy = url;
   dummy += "/dummy";

   XrdClientAdmin *cadm = TXNetSystem::GetClientAdmin(dummy);

   if (!cadm) {
      Error("Connect","fatal error: new object creation failed.");
      return cadm;
   }

   // Do not block: restore old value after
   Int_t maxOld = EnvGetLong(NAME_FIRSTCONNECTMAXCNT);

   // Try to connect to the server
   gEnv->SetValue("XNet.FirstConnectMaxCnt", 1);
   EnvPutInt(NAME_FIRSTCONNECTMAXCNT, 1);
   if (cadm->Connect()) {
      fIsXRootd = kTRUE;
      EnvPutInt(NAME_FIRSTCONNECTMAXCNT, maxOld);
   } else {
      EnvPutInt(NAME_FIRSTCONNECTMAXCNT, maxOld);
      if (fgRootdBC) {
         Bool_t isRootd =
            (cadm->GetClientConn()->GetServerType() == kSTRootd);
         Int_t sd = cadm->GetClientConn()->GetOpenSockFD();
         if (isRootd && sd > -1) {
            //
            // Create a TSocket on the open connection
            TSocket *s = new TSocket(sd);

            // We will clean it by ourselves
            R__LOCKGUARD(gROOTMutex);
            gROOT->GetListOfSockets()->Remove(s);

            s->SetOption(kNoBlock, 0);

            // Find out the remote protocol (send the client protocol first)
            Int_t rproto = TXNetFile::GetRootdProtocol(s);
            if (rproto < 0) {
               Error("TXNetSystem", "getting protocol of the rootd server");
               cadm = 0;
               return 0;
            }
            // Finalize TSocket initialization
            s->SetRemoteProtocol(rproto);
            TUrl uut((cadm->GetClientConn()
                             ->GetCurrentUrl()).GetUrl().c_str());
            TString uu;
            TXNetFile::FormUrl(uut,uu);
            if (gDebug > 2)
               Info("Connect"," url: %s",uu.Data());

            s->SetUrl(uu.Data());
            s->SetService("rootd");
            s->SetServType(TSocket::kROOTD);
            //
            // Now we can check if we can create a TNetFile on the
            // open connection
            if (rproto > 13) {
               //
               // Remote support for reuse of open connection
               TNetSystem::Create(uu, s);
            } else {
               //
               // Open connection has been closed because could
               // not be reused; TNetSystem will open a new connection
               TNetSystem::Create(uu);
            }

            // Type of server
            fIsRootd = kTRUE;
            cadm = 0;

         } else {
            Error("Connect", "some severe error occurred while opening"
                  " the connection at %s - exit", url);
            if (cadm && cadm->LastServerError())
               Printf("   '%s'", cadm->LastServerError()->errmsg);
            else
               Printf("   (error message not available)");
            cadm = 0;
            return cadm;
         }
      } else {
         Error("Connect",
               "while opening the connection at %s - exit", url);
         cadm = 0;
         return cadm;
      }
   }

   return cadm;
}

////////////////////////////////////////////////////////////////////////////////
/// One-time initialization of some communication variables for xrootd protocol

void TXNetSystem::InitXrdClient()
{
   // Init vars with default debug level -1, so we do not get warnings
   TXNetFile::SetEnv();

   // Only once
   fgInitDone = kTRUE;

   // Print the tag, if required (only once)
   if (gEnv->GetValue("XNet.PrintTAG",0) == 1)
     Info("TXNetFile","(C) 2005 SLAC TXNetSystem (eXtended TNetSystem) %s",
         gROOT->GetVersion());
}

////////////////////////////////////////////////////////////////////////////////
/// Open a directory. Returns a non-zero pointer (with no special
/// purpose) in case of success, 0 in case of error.

void* TXNetSystem::OpenDirectory(const char* dir)
{
   if (fIsXRootd) {
      // Check if the directory exists
      TXNetSystemConnectGuard cg(this, dir);
      if (cg.IsValid()) {
         fUrl = dir;
         // Extract the directory name
         fDir = TUrl(dir).GetFile();
         fDirp = (void*)&fDir;     // serves as directory pointer

         vecString dirs;
         vecBool existDirs;
         XrdOucString s(fDir.Data());
         dirs.Push_back(s);
         cg.ClientAdmin()->ExistDirs(dirs, existDirs);
         cg.ClientAdmin()->GoBackToRedirector();
         if (existDirs.GetSize()>0 && existDirs[0])
            return fDirp;
         else
            cg.NotifyLastError();
      }
      return 0;
   }

   if (gDebug > 1)
      Info("OpenDirectory", "calling TNetSystem::OpenDirectory");
   return TNetSystem::OpenDirectory(dir);       // for a rootd
}

////////////////////////////////////////////////////////////////////////////////
/// Free(Close) the directory referenced by dirp

void TXNetSystem::FreeDirectory(void *dirp)
{
   if (fIsXRootd) {
      if (dirp != fDirp) {
         Error("FreeDirectory","invalid directory pointer (%p, %p)", dirp, fDirp);
         return;
      }
      fDir = "";
      fDirp = 0;
      fDirEntry = "";
      if (fDirList) {
         ((VecString_t*)fDirList)->Clear();
         delete ((VecString_t*)fDirList);
         fDirList = 0;
      }
      return;
   }

   if (gDebug > 1)
      Info("FreeDirectory","calling TNetSystem::FreeDirectory");
   return TNetSystem::FreeDirectory(dirp);     // for a rootd
}

////////////////////////////////////////////////////////////////////////////////
/// Create a directory. Return 0 on success, -1 otherwise.

Int_t TXNetSystem::MakeDirectory(const char* dir)
{
   if (fIsXRootd) {
      TXNetSystemConnectGuard cg(this, dir);
      if (cg.IsValid()) {
         // use default permissions 755 to create directory
         Bool_t ok = cg.ClientAdmin()->Mkdir(TUrl(dir).GetFile(),7,5,5);
         cg.ClientAdmin()->GoBackToRedirector();
         if (ok) {
            return 0;
         } else {
            cg.NotifyLastError();
            return -1;
         }
      }
   }

   if (gDebug > 1)
      Info("MakeDirectory","Calling TNetSystem::MakeDirectory");
   return TNetSystem::MakeDirectory(dir);     // for a rootd
}

////////////////////////////////////////////////////////////////////////////////
/// Get directory entry for directory referenced by dirp.
/// Returns 0 in case there are no more entries.

const char* TXNetSystem::GetDirEntry(void *dirp)
{
   if (fIsXRootd) {
      if (dirp != fDirp) {
         Error("GetDirEntry","invalid directory pointer");
         return 0;
      }

      // Only request new directory listing the first time called
      if (!fDirList) {
         TXNetSystemConnectGuard cg(this, fUrl);
         if (cg.IsValid()) {
            fDirList = new VecString_t;
            Bool_t ok = cg.ClientAdmin()->DirList(fDir, *(VecString_t*)fDirList);
            cg.ClientAdmin()->GoBackToRedirector();
            if (!ok) {
               cg.NotifyLastError();
               delete (VecString_t*)fDirList;
               fDirList = 0;
               return 0;
            }
         }
      }

      // Return entries one by one with each call of method
      if (fDirList && ((VecString_t*)fDirList)->GetSize() > 0) {
         fDirEntry = ((VecString_t*)fDirList)->Pop_front().c_str();
         return fDirEntry.Data();
      }
      return 0;   // until all of them have been returned
   }

   if (gDebug > 1) Info("GetDirEntry","Calling TNetSystem::GetDirEntry");
   return TNetSystem::GetDirEntry(dirp);      // for a rootd
}

////////////////////////////////////////////////////////////////////////////////
/// Get info about a file. Info is returned in the form of a FileStat_t
/// structure (see TSystem.h).
/// The function returns 0 in case of success and 1 if the file could
/// not be stat'ed.
/// NOTICE: Not all information is available with an xrootd server.

Int_t TXNetSystem::GetPathInfo(const char* path, FileStat_t &buf)
{
   if (fIsXRootd) {
      TXNetSystemConnectGuard cg(this, path);
      if (cg.IsValid()) {

         Long_t id;
         Long64_t size;
         Long_t flags;
         Long_t modtime;

         // Issue the request
         TUrl urlpath(path);
         Bool_t ok = cg.ClientAdmin()->Stat(urlpath.GetFile(),id,size,flags,modtime);
         if (ok) {
            // Save the endpoint path
            urlpath.SetProtocol(cg.ClientAdmin()->GetCurrentUrl().Proto.c_str());
            urlpath.SetHost(cg.ClientAdmin()->GetCurrentUrl().Host.c_str());
            urlpath.SetPort(cg.ClientAdmin()->GetCurrentUrl().Port);
            buf.fUrl = urlpath.GetUrl();
         }
         cg.ClientAdmin()->GoBackToRedirector();

         // Flag offline files
         if (flags & kXR_offline) {
            buf.fMode = kS_IFOFF;
         } else if (ok) {
            buf.fDev = (id >> 24);
            buf.fIno = (id & 0x00FFFFFF);
            buf.fUid = -1;       // not all information available in xrootd
            buf.fGid = -1;       // not available
            buf.fSize = size;
            buf.fMtime = modtime;

            if (flags == 0) buf.fMode = kS_IFREG;
            if (flags & kXR_xset) buf.fMode = (kS_IFREG|kS_IXUSR|kS_IXGRP|kS_IXOTH);
            if (flags & kXR_isDir) buf.fMode = kS_IFDIR;
            if (flags & kXR_other) buf.fMode = kS_IFSOCK;
            if (flags & kXR_readable) buf.fMode |= kS_IRUSR;
            if (flags & kXR_writable) buf.fMode |= kS_IWUSR;

            buf.fIsLink = 0;     // not available
            return 0;
         } else {
            if (gDebug > 0)
               cg.NotifyLastError();
         }
      }
      return 1;
   }

   if (gDebug > 1)
      Info("GetPathInfo","Calling TNetSystem::GetPathInfo");
   return TNetSystem::GetPathInfo(path,buf);       // for a rootd
}

////////////////////////////////////////////////////////////////////////////////
/// Check consistency of this helper with the one required
/// by 'path' or 'dirptr'.

Bool_t TXNetSystem::ConsistentWith(const char *path, void *dirptr)
{
   if (gDebug > 1)
      Info("ConsistentWith",
           "calling for path: %s, dir: %p", path, dirptr);

   return TNetSystem::ConsistentWith(path,dirptr);    // for a rootd
}

////////////////////////////////////////////////////////////////////////////////
/// Returns FALSE if one can access a file using the specified access mode.
/// NB: for the time being mode is ignored for XROOTD (just checks existence
/// of the file or directory).
/// Mode is the same as for the Unix access(2) function.
/// Attention, bizarre convention of return value!!

Bool_t TXNetSystem::AccessPathName(const char *path, EAccessMode mode)
{
   if (fIsXRootd) {
      // Check only if the file or directory exists and
      FileStat_t buf;
      if (GetPathInfo(path, buf) == 0)
         if (buf.fMode != kS_IFSOCK)
            return kFALSE;
      // The file could not be stated
      return kTRUE;
   }

   if (gDebug > 1)
      Info("AccessPathName", "calling TNetSystem::AccessPathName");
   return TNetSystem::AccessPathName(path,mode);    // for a rootd
}

////////////////////////////////////////////////////////////////////////////////
/// Unlink 'path' on the remote server system.
/// Returns 0 on success, -1 otherwise.

int TXNetSystem::Unlink(const char *path)
{
   if (fIsXRootd) {

      TXNetSystemConnectGuard cg(this, path);
      if (cg.IsValid()) {

         Long_t id;
         Long64_t size;
         Long_t flags;
         Long_t modtime;

         // Extract the directory name
         TString edir = TUrl(path).GetFile();
         Bool_t ok = cg.ClientAdmin()->Stat(edir.Data(), id, size, flags, modtime);

         // Flag offline files
         if (ok && !(flags & kXR_offline)) {
            if (flags & kXR_isDir)
               ok = cg.ClientAdmin()->Rmdir(edir.Data());
            else
               ok = cg.ClientAdmin()->Rm(edir.Data());
            cg.ClientAdmin()->GoBackToRedirector();

            // Done
            return ((ok) ? 0 : -1);
         } else if (!ok) {
            cg.ClientAdmin()->GoBackToRedirector();
            cg.NotifyLastError();
         }
      }
   }

   if (gDebug > 1)
      Info("Unlink", "calling TNetSystem::Unlink");
   return -1;    // not implemented for rootd
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the file defined by 'path' is ready to be used

Bool_t TXNetSystem::IsOnline(const char *path)
{
   // This is most efficiently done using GetPathInfo
   FileStat_t st;
   if (GetPathInfo(path, st) != 0) {
      if (gDebug > 0)
         Info("IsOnline", "path '%s' cannot be stat'ed", path);
      return kFALSE;
   }
   if (R_ISOFF(st.fMode)) {
      if (gDebug > 0)
         Info("IsOnline", "path '%s' is offline", path);
      return kFALSE;
   }
   // Done
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Issue a prepare request for file defined by 'path'

Bool_t TXNetSystem::Prepare(const char *path, UChar_t option, UChar_t priority)
{
   TXNetSystemConnectGuard cg(this, path);
   if (cg.IsValid()) {
      XrdOucString pathname = TUrl(path).GetFileAndOptions();
      vecString vs;
      vs.Push_back(pathname);
      cg.ClientAdmin()->Prepare(vs, (kXR_char)option, (kXR_char)priority);
      cg.ClientAdmin()->GoBackToRedirector();
      if (gDebug >0)
         Info("Prepare", "Got Status %d for %s",
              cg.ClientAdmin()->LastServerResp()->status, pathname.c_str());
      if (!(cg.ClientAdmin()->LastServerResp()->status)){
         return kTRUE;
      }
      cg.NotifyLastError();
   }

   // Done
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Issue a prepare request for a list of files defined by 'paths', which must
/// be of one of the following types: TFileInfo, TUrl, TObjString.
/// On output, bufout, if defined, points to a buffer form that can be used
/// with GetPathsInfo.
/// Return the number of paths found or -1 if any error occured.

Int_t TXNetSystem::Prepare(TCollection *paths,
                           UChar_t opt, UChar_t prio, TString *bufout)
{
   if (!paths) {
      Warning("Prepare", "input list is empty!");
      return -1;
   }

   Int_t npaths = 0;

   TXNetSystemConnectGuard cg(this, "");
   if (cg.IsValid()) {

      TString *buf = (bufout) ? bufout : new TString();

      // Prepare the buffer
      TObject *o = 0;
      TUrl u;
      TString path;
      TIter nxt(paths);
      while ((o = nxt()))  {
         // Extract the path name from the allowed object types
         TString pn = TFileStager::GetPathName(o);
         if (pn == "") {
            Warning("Prepare", "object is of unexpected type %s - ignoring", o->ClassName());
            continue;
         }
         u.SetUrl(pn);
         // The path
         path = u.GetFileAndOptions();
         path.ReplaceAll("\n","\r");
         npaths++;
         *buf += Form("%s\n", path.Data());
      }

      Info("Prepare","buffer ready: issuing prepare (opt=%d, prio=%d) ...",
         opt, prio);
      cg.ClientAdmin()->Prepare(buf->Data(), (kXR_char)opt, (kXR_char)prio);
      cg.ClientAdmin()->GoBackToRedirector();
      if (!bufout)
         delete buf;
      if (gDebug >0)
         Info("Prepare", "Got Status %d",
              cg.ClientAdmin()->LastServerResp()->status);
      if (!(cg.ClientAdmin()->LastServerResp()->status)){
         return npaths;
      }
      cg.NotifyLastError();
   }

   // Done
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve status of a '\n'-separated list of files in 'paths'.
/// The information is returned as one UChar_t per file in 'info';
/// 'info' must be allocated by the caller.

Bool_t TXNetSystem::GetPathsInfo(const char *paths, UChar_t *info)
{
   if (!paths) {
      Warning("GetPathsInfo", "input list is empty!");
      return kFALSE;
   }

   TXNetSystemConnectGuard cg(this, "");
   if (cg.IsValid()) {
      cg.ClientAdmin()->SysStatX(paths, info);
      cg.ClientAdmin()->GoBackToRedirector();
      if (gDebug >0)
         Info("GetPathsInfo", "Got Status %d",
              cg.ClientAdmin()->LastServerResp()->status);
      if (!(cg.ClientAdmin()->LastServerResp()->status)){
         return kTRUE;
      }
      cg.NotifyLastError();
   }

   // Done
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns TRUE if the url in 'path' points to the local file system.
/// This is used to avoid going through the NIC card for local operations.

Bool_t TXNetSystem::IsPathLocal(const char *path)
{
   if (fIsXRootd) {
      TXNetSystemConnectGuard cg(this, path);
      if (cg.IsValid()) {
         if (cg.ClientAdmin()->GetClientConn()->GetServerType() != kSTDataXrootd) {
            // Not an end point data server: cannot assert locality
            return kFALSE;
         }
      }
   }
   // Either an end-point data server or 'rootd': check for locality
   return TSystem::IsPathLocal(path);
}

////////////////////////////////////////////////////////////////////////////////
/// Get end-point url of a file. Info is returned in eurl.
/// The function returns 0 in case of success and 1 if the file could
/// not be stat'ed.

Int_t TXNetSystem::Locate(const char *path, TString &eurl)
{
   if (fIsXRootd) {
      TXNetSystemConnectGuard cg(this, path);
      if (cg.IsValid()) {

         // Extract the directory name
         XrdClientLocate_Info li;
         TString edir = TUrl(path).GetFile();

         if (cg.ClientAdmin()->Locate((kXR_char *)edir.Data(), li)) {
            TUrl u(path);
            XrdClientUrlInfo ui((const char *)&li.Location[0]);
            // We got the IP address but we need the FQDN: if we did not resolve
            // it yet do it and cache the result
            TNamed *hn = 0;
            if (fgAddrFQDN.GetSize() <= 0 ||
               !(hn = dynamic_cast<TNamed *>(fgAddrFQDN.FindObject(ui.Host.c_str())))) {
               TInetAddress a(gSystem->GetHostByName(ui.Host.c_str()));
               if (strlen(a.GetHostName()) > 0)
                  hn = new TNamed(ui.Host.c_str(), a.GetHostName());
               else
                  hn = new TNamed(ui.Host.c_str(), ui.Host.c_str());
               fgAddrFQDN.Add(hn);
               if (gDebug > 0)
                  Info("Locate","caching host name: %s", hn->GetTitle());
            }
            if (hn)
               u.SetHost(hn->GetTitle());
            else
               u.SetHost(ui.Host.c_str());
            u.SetPort(ui.Port);
            eurl = u.GetUrl();
            return 0;
         }
         cg.NotifyLastError();
      }
      return 1;
   }

   // Not implemented
   if (gDebug > 0) Info("Locate", "server not Xrootd: method not implemented!");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if an admin for 'url' exists already.
/// Avoid duplications.

XrdClientAdmin *TXNetSystem::GetClientAdmin(const char *url)
{
   XrdClientAdmin *ca = 0;

   // ID key
   TString key = TXNetSystem::GetKey(url);

   // If we have one for 'key', just use it
   TXrdClientAdminWrapper *caw = 0;
   if (fgAdminHash.GetSize() > 0 &&
      (caw = dynamic_cast<TXrdClientAdminWrapper *>(fgAdminHash.FindObject(key.Data()))))
      return caw->fXCA;

   // Create one and save the reference in the hash table
   ca = new XrdClientAdmin(url);
   fgAdminHash.Add(new TXrdClientAdminWrapper(key, ca));

   // Done
   return ca;
}

////////////////////////////////////////////////////////////////////////////////
/// Build from uu a unique ID key used in hash tables

TString TXNetSystem::GetKey(const char *url)
{
   TUrl u(url);
   TString key(u.GetUser());
   if (!key.IsNull())
      key += "@";
   key += u.GetHost();
   if (u.GetPort() > 0) {
      key += ":";
      key += u.GetPort();
   }

   // Done
   return key;
}

//
// Wrapper class
//
////////////////////////////////////////////////////////////////////////////////
/// Destructor: destroy the instance

TXrdClientAdminWrapper::~TXrdClientAdminWrapper()
{
   SafeDelete(fXCA);
}

//
// Guard methods
//
////////////////////////////////////////////////////////////////////////////////
/// Construct a guard object

TXNetSystemConnectGuard::TXNetSystemConnectGuard(TXNetSystem *xn, const char *url)
                        : fClientAdmin(0)
{
    if (xn)
       // Connect
       fClientAdmin = (url && strlen(url) > 0) ? xn->Connect(url)
                                               : xn->Connect(xn->fUrl);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor: close the connection

TXNetSystemConnectGuard::~TXNetSystemConnectGuard()
{
   fClientAdmin = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Print message about last occured error

void TXNetSystemConnectGuard::NotifyLastError()
{
   if (fClientAdmin)
      if (fClientAdmin->GetClientConn())
         Printf("Srv err: %s", fClientAdmin->GetClientConn()->LastServerError.errmsg);
}
