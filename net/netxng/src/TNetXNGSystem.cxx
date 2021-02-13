// @(#)root/netx:$Id$
/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// TNetXNGSystem                                                              //
//                                                                            //
// Authors: Justin Salmon, Lukasz Janyst                                      //
//          CERN, 2013                                                        //
//                                                                            //
// Enables access to XRootD filesystem interface using the new client.        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "TNetXNGSystem.h"
#include "TFileStager.h"
#include "Rtypes.h"
#include "TList.h"
#include "TUrl.h"
#include "TVirtualMutex.h"
#include <XrdCl/XrdClFileSystem.hh>
#include <XrdCl/XrdClXRootDResponses.hh>
#include <XrdVersion.hh>
#if XrdVNUMBER >= 40000
#include <XrdNet/XrdNetAddr.hh>
#else
#include <XrdSys/XrdSysDNS.hh>
#endif


////////////////////////////////////////////////////////////////////////////////
/// PluginManager creation function

TSystem* ROOT_Plugin_TNetXNGSystem(const char *url, Bool_t owner) {
   return new TNetXNGSystem(url, owner);
}


ClassImp(TNetXNGSystem);

THashList TNetXNGSystem::fgAddrFQDN;
TMutex TNetXNGSystem::fgAddrMutex;

namespace
{
   struct DirectoryInfo {
      XrdCl::URL                     *fUrl;         // Path of this directory
      XrdCl::DirectoryList           *fDirList;     // Directory listing
      XrdCl::DirectoryList::Iterator *fDirListIter; // Iterator for this listing

   public:
      DirectoryInfo(const char *dir) : fUrl(new XrdCl::URL(dir)), fDirList(0), fDirListIter(0) {}

      ~DirectoryInfo() {
        delete fUrl;
        delete fDirList;
      }
   };
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor: Create system class without connecting to server
///
/// param owner: (unused)

TNetXNGSystem::TNetXNGSystem(Bool_t /*owner*/) :
   TSystem("-root", "Net file Helper System"), fUrl(0), fFileSystem(0)
{
   // Name must start with '-' to bypass the TSystem singleton check, then
   // be changed to "root"
   SetName("root");
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor: Create system class and connect to server
///
/// param url:   URL of the entry-point server to be contacted
/// param owner: (unused)

TNetXNGSystem::TNetXNGSystem(const char *url, Bool_t /*owner*/) :
   TSystem("-root", "Net file Helper System")
{
   using namespace XrdCl;

   // Name must start with '-' to bypass the TSystem singleton check
   SetName("root");
   fUrl        = new URL(std::string(url));
   fFileSystem = new FileSystem(fUrl->GetURL());
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TNetXNGSystem::~TNetXNGSystem()
{
   delete fFileSystem;
   delete fUrl;
}

////////////////////////////////////////////////////////////////////////////////
/// Open a directory
///
/// param dir: the name of the directory to open
/// returns:   a non-zero pointer (with no special purpose) in case of
///            success, 0 in case of error

void* TNetXNGSystem::OpenDirectory(const char *dir)
{
   using namespace XrdCl;

   DirectoryInfo *dirInfo = new DirectoryInfo(dir);
   fDirPtrs.insert( (void*)dirInfo );
   return (void *) dirInfo;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a directory
///
/// param dir: the directory name
/// returns:   0 on success, -1 otherwise

Int_t TNetXNGSystem::MakeDirectory(const char *dir)
{
   using namespace XrdCl;
   URL url(dir);
   XRootDStatus st = fFileSystem->MkDir(url.GetPath(), MkDirFlags::MakePath,
                                        Access::None);
   if (!st.IsOK()) {
      Error("MakeDirectory", "%s", st.GetErrorMessage().c_str());
      return -1;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Free a directory
///
/// param dirp: the pointer to the directory to be freed

void TNetXNGSystem::FreeDirectory(void *dirp)
{
   fDirPtrs.erase( dirp );
   delete (DirectoryInfo *) dirp;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a directory entry.
///
/// param dirp: the directory pointer
/// returns:    0 in case there are no more entries

const char* TNetXNGSystem::GetDirEntry(void *dirp)
{
   using namespace XrdCl;
   DirectoryInfo *dirInfo = (DirectoryInfo *) dirp;

   if (!dirInfo->fDirList) {
      XRootDStatus st = fFileSystem->DirList(dirInfo->fUrl->GetPath(),
                                             DirListFlags::Locate,
                                             dirInfo->fDirList);
      if (!st.IsOK()) {
         Error("GetDirEntry", "%s", st.GetErrorMessage().c_str());
         return 0;
      }
      dirInfo->fDirListIter = new DirectoryList::Iterator(dirInfo->fDirList->Begin());
   }

   if (*(dirInfo->fDirListIter) != dirInfo->fDirList->End()) {
      const char *filename = (**(dirInfo->fDirListIter))->GetName().c_str();
      (*(dirInfo->fDirListIter))++;
      return filename;

   } else {
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get info about a file (stat)
///
/// param path: the path of the file to stat (in)
/// param buf:  structure that will hold the stat info (out)
/// returns:    0 if success, 1 if the file could not be stat'ed

Int_t TNetXNGSystem::GetPathInfo(const char *path, FileStat_t &buf)
{
   using namespace XrdCl;
   StatInfo *info = 0;
   URL target(path);
   XRootDStatus st = fFileSystem->Stat(target.GetPath(), info);

   if (!st.IsOK()) {

      if (gDebug > 1) {
         Info("GetPathInfo", "Stat error: %s", st.GetErrorMessage().c_str());
      }
      delete info;
      return 1;

   } else {

      // Flag offline files
      if (info->TestFlags(StatInfo::Offline)) {
         buf.fMode = kS_IFOFF;
      } else {
         std::stringstream sstr(info->GetId());
         Long64_t id;
         sstr >> id;

         buf.fDev    = (id >> 32);
         buf.fIno    = (id & 0x00000000FFFFFFFF);
         buf.fUid    = -1;  // not available
         buf.fGid    = -1;  // not available
         buf.fIsLink = 0;   // not available
         buf.fSize   = info->GetSize();
         buf.fMtime  = info->GetModTime();

         if (info->TestFlags(StatInfo::XBitSet))
            buf.fMode = (kS_IFREG | kS_IXUSR | kS_IXGRP | kS_IXOTH);
         if (info->GetFlags() == 0)                 buf.fMode = kS_IFREG;
         if (info->TestFlags(StatInfo::IsDir))      buf.fMode = kS_IFDIR;
         if (info->TestFlags(StatInfo::Other))      buf.fMode = kS_IFSOCK;
         if (info->TestFlags(StatInfo::IsReadable)) buf.fMode |= kS_IRUSR;
         if (info->TestFlags(StatInfo::IsWritable)) buf.fMode |= kS_IWUSR;
      }
   }

   delete info;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check consistency of this helper with the one required by 'path' or
/// 'dirptr'
///
/// param path:   the path to check
/// param dirptr: the directory pointer to check

Bool_t TNetXNGSystem::ConsistentWith(const char *path, void *dirptr)
{
   using namespace XrdCl;

   if( path )
   {
      URL url(path);

      if( gDebug > 1 )
         Info("ConsistentWith", "Protocol: '%s' (%s), Username: '%s' (%s), "
              "Password: '%s' (%s), Hostname: '%s' (%s), Port: %d (%d)",
               fUrl->GetProtocol().c_str(), url.GetProtocol().c_str(),
               fUrl->GetUserName().c_str(), url.GetUserName().c_str(),
               fUrl->GetPassword().c_str(), url.GetPassword().c_str(),
               fUrl->GetHostName().c_str(), url.GetHostName().c_str(),
               fUrl->GetPort(), url.GetPort());

      // Require match of protocol, user, password, host and port
      if( fUrl->GetProtocol() == url.GetProtocol() &&
          fUrl->GetUserName() == url.GetUserName() &&
          fUrl->GetPassword() == url.GetPassword() &&
          fUrl->GetHostName() == url.GetHostName() &&
          fUrl->GetPort() == url.GetPort())
         return kTRUE;
   }

   if( dirptr )
      return fDirPtrs.find( dirptr ) != fDirPtrs.end();

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Unlink a file on the remote server
///
/// param path: the path of the file to unlink
/// returns:    0 on success, -1 otherwise

int TNetXNGSystem::Unlink(const char *path)
{
   using namespace XrdCl;
   StatInfo *info;
   URL url(path);

   // Stat the path to find out if it's a file or a directory
   XRootDStatus st = fFileSystem->Stat(url.GetPath(), info);
   if (!st.IsOK()) {
      Error("Unlink", "%s", st.GetErrorMessage().c_str());
      return -1;
   }

   if (info->TestFlags(StatInfo::IsDir))
      st = fFileSystem->RmDir(url.GetPath());
   else
      st = fFileSystem->Rm(url.GetPath());
   delete info;

   if (!st.IsOK()) {
      Error("Unlink", "%s", st.GetErrorMessage().c_str());
      return -1;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Is this path a local path?
///
/// param path: the URL of the path to check
/// returns:    kTRUE if the path is local, kFALSE otherwise

Bool_t TNetXNGSystem::IsPathLocal(const char *path)
{
   return TSystem::IsPathLocal(path);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the endpoint URL of a file.
///
/// param path:   the entry-point URL of the file (in)
/// param endurl: the endpoint URL of the file (out)
/// returns:      0 in case of success and 1 if the file could not be
///               stat'ed.

Int_t TNetXNGSystem::Locate(const char *path, TString &endurl)
{
   using namespace XrdCl;
   LocationInfo *info = 0;
   URL pathUrl(path);

   // Locate the file
   XRootDStatus st = fFileSystem->Locate(pathUrl.GetPath(), OpenFlags::None,
                                         info);
   if (!st.IsOK()) {
      Error("Locate", "%s", st.GetErrorMessage().c_str());
      delete info;
      return 1;
   }

   // Use the first endpoint address returned by the client
   URL locUrl(info->Begin()->GetAddress());
   TString loc = locUrl.GetHostName();
   delete info;
   info = 0;

   R__LOCKGUARD(&fgAddrMutex);

   // The location returned by the client library is the numeric host address
   // without path. Try to lookup a hostname and replace the path portion of
   // the url before returning the result.
   TNamed *hn = 0;
   if (fgAddrFQDN.GetSize() <= 0 ||
       !(hn = dynamic_cast<TNamed *>(fgAddrFQDN.FindObject(loc)))) {
#if XrdVNUMBER >= 40000
      XrdNetAddr netaddr;
      netaddr.Set(loc.Data());
      const char* name = netaddr.Name();
      if (name) {
         hn = new TNamed(loc.Data(), name);
      } else {
         hn = new TNamed(loc, loc);
      }
#else
      char *addr[1] = {0}, *name[1] = {0};
      int naddr = XrdSysDNS::getAddrName(loc.Data(), 1, addr, name);
      if (naddr == 1) {
         hn = new TNamed(loc.Data(), name[0]);
      } else {
         hn = new TNamed(loc, loc);
      }
      free(addr[0]);
      free(name[0]);
#endif
      fgAddrFQDN.Add(hn);
      if (gDebug > 0)
         Info("Locate","caching host name: %s", hn->GetTitle());
   }

   TUrl res(path);
   res.SetHost(hn->GetTitle());
   res.SetPort(locUrl.GetPort());
   endurl = res.GetUrl();

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Issue a stage request for a single file
///
/// param path: the path of the file to stage
/// param opt:  defines 'option' and 'priority' for 'Prepare': the format is
///             opt = "option=o priority=p"
/// returns:    0 for success, -1 for error

Int_t TNetXNGSystem::Stage(const char* path, UChar_t priority)
{
   TList *files = new TList();
   files->Add((TObject *) new TUrl(path));
   return Stage((TCollection *) files, priority);
}

////////////////////////////////////////////////////////////////////////////////
/// Issue stage requests for multiple files
///
/// param pathlist: list of paths of files to stage
/// param opt:      defines 'option' and 'priority' for 'Prepare': the
///                 format is opt = "option=o priority=p"
/// returns:        0 for success, -1 for error

Int_t TNetXNGSystem::Stage(TCollection *files, UChar_t priority)
{
   using namespace XrdCl;
   std::vector<std::string> fileList;
   TIter it(files);
   TObject *object = 0;

   while ((object = (TObject *) it.Next())) {

      TString path = TFileStager::GetPathName(object);
      if (path == "") {
         Warning("Stage", "object is of unexpected type %s - ignoring",
               object->ClassName());
         continue;
      }

      fileList.push_back(std::string(URL(path.Data()).GetPath()));
   }

   Buffer *response;
   XRootDStatus st = fFileSystem->Prepare(fileList, PrepareFlags::Stage,
                                          (uint8_t) priority, response);
   if (!st.IsOK()) {
      Error("Stage", "%s", st.GetErrorMessage().c_str());
      return -1;
   }

   return 0;
}

