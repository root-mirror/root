// @(#)root/base:$Id$
// Author: Jan Fiete Grosse-Oetringhaus, 04.06.07

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TDataSetManagerFile                                             //
//                                                                      //
// Implementation of TDataSetManager handling datasets from root   //
// files under a specific directory path                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <sys/stat.h>

#include "TDataSetManagerFile.h"

#include "Riostream.h"
#include "TDatime.h"
#include "TEnv.h"
#include "TFileCollection.h"
#include "TFileInfo.h"
#include "TFile.h"
#include "TFileStager.h"
#include "TLockFile.h"
#include "TMap.h"
#include "TRegexp.h"
#include "TMD5.h"
#include "TMacro.h"
#include "TMessage.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TError.h"
#include "TPRegexp.h"
#include "TVirtualMonitoring.h"
#include "THashList.h"
#include "TKey.h"
#include "TTree.h"
#include "TParameter.h"

struct LsTreeEntry_t {
   TObjString  *fGrp;      // Group
   TObjString  *fUsr;      // User
   TObjString  *fMd5;      // Checksum
   TObjString  *fLss;      // 'ls' string
   Long64_t     fMtime;    // modification time
   LsTreeEntry_t(const char *g, const char *u, const char *cs, const char *ls, Long64_t m) :
         fMtime(m) { fGrp = new TObjString(g); fUsr = new TObjString(u);
                     fMd5 = new TObjString(cs); fLss = new TObjString(ls); }
   ~LsTreeEntry_t() { SafeDelete(fGrp); SafeDelete(fUsr); SafeDelete(fMd5); SafeDelete(fLss);}
};

ClassImp(TDataSetManagerFile);

////////////////////////////////////////////////////////////////////////////////
///
/// Main constructor

TDataSetManagerFile::TDataSetManagerFile(const char *group,
                                         const char *user, const char *ins)
                    : TDataSetManager(group, user, ins)
{
   // Parse options
   ParseInitOpts(ins);

   // Init the instance
   Init();
}

////////////////////////////////////////////////////////////////////////////////
///
/// Main constructor

TDataSetManagerFile::TDataSetManagerFile(const char *ins)
                    : TDataSetManager("", "", ins)
{
   // Parse options
   ParseInitOpts(ins);

   // Init the instance
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Do the real inititialization

void TDataSetManagerFile::Init()
{
   fIsRemote = kFALSE;
   if (!fUser.IsNull() && !fGroup.IsNull() && !fDataSetDir.IsNull()) {

      // Make sure that the dataset dir exists
      TString dir;
      dir.Form("%s/%s/%s", fDataSetDir.Data(), fGroup.Data(), fUser.Data());
      if (gSystem->AccessPathName(dir)) {
         if (gSystem->mkdir(dir, kTRUE) != 0) {
            TString emsg = dir;
            // Read only dataset info system: switch to COMMON
            fUser = fCommonUser;
            fGroup = fCommonGroup;
            ResetBit(TDataSetManager::kCheckQuota);
            ResetBit(TDataSetManager::kAllowRegister);
            ResetBit(TDataSetManager::kAllowVerify);
            ResetBit(TDataSetManager::kTrustInfo);
            dir.Form("%s/%s/%s", fDataSetDir.Data(), fGroup.Data(), fUser.Data());
            if (gSystem->AccessPathName(dir)) {
               Error("Init",
                     "could not attach to a valid the dataset dir; paths tried:");
               Error("Init", "    %s", emsg.Data());
               Error("Init", "    %s", dir.Data());
               SetBit(TObject::kInvalidObject);
               return;
            }
         }
         else if (fOpenPerms) {

            // Directory creation was OK: let's open permissions if requested
            TString t;
            Int_t rr = 0;

            t.Form("%s/%s/%s", fDataSetDir.Data(), fGroup.Data(), fUser.Data());
            rr += gSystem->Chmod(t.Data(), 0777);

            t.Form("%s/%s", fDataSetDir.Data(), fGroup.Data());
            rr += gSystem->Chmod(t.Data(), 0777);

            rr += gSystem->Chmod(fDataSetDir.Data(), 0777);

            if (rr < 0) {
               t.Form("%s/%s/%s", fDataSetDir.Data(), fGroup.Data(),
                 fUser.Data());
               Warning("Init",
                 "problems setting perms of dataset directory %s (#%d)",
                 t.Data(), TSystem::GetErrno());
            }

         }
      }

      // If not in sandbox, construct the base URI using session defaults
      // (group, user) (syntax: /group/user/dsname[#[subdir/]objname])
      if (!TestBit(TDataSetManager::kIsSandbox))
         fBase.SetUri(TString(Form("/%s/%s/", fGroup.Data(), fUser.Data())));

      // Local or remote?
      TString locPath;
      TFile::EFileType pathType = TFile::GetType(fDataSetDir, "READ", &locPath);
      if (pathType == TFile::kLocal) {
         // Remote Url pointing to this machine
         fDataSetDir = locPath;
         if (gDebug > 0)
            Info("Init", "repository '%s' is local", fDataSetDir.Data());
      } else if (pathType != TFile::kDefault && pathType != TFile::kFile) {
         fIsRemote = kTRUE;
         if (gDebug > 0)
            Info("Init", "repository '%s' is remote", fDataSetDir.Data());
      }

      // Read locking path from kDataSet_LockLocation
      TString lockloc = TString::Format("%s/%s", fDataSetDir.Data(), kDataSet_LockLocation);
      if (!gSystem->AccessPathName(lockloc, kReadPermission)) {
         // Open the file in RAW mode
         lockloc += "?filetype=raw";
         TFile *f = TFile::Open(lockloc);
         if (f && !(f->IsZombie())) {
            const Int_t blen = 8192;
            char buf[blen];
            Long64_t rest = f->GetSize();
            while (rest > 0) {
               Long64_t len = (rest > blen - 1) ? blen - 1 : rest;
               if (f->ReadBuffer(buf, len)) {
                  fDataSetLockFile = "";
                  break;
               }
               buf[len] = '\0';
               fDataSetLockFile += buf;
               rest -= len;
            }
            f->Close();
            SafeDelete(f);
            fDataSetLockFile.ReplaceAll("\n","");
         } else {
            lockloc.ReplaceAll("?filetype=raw", "");
            Warning("Init", "could not open remore file '%s' with the lock location", lockloc.Data());
         }
      }
      if (fDataSetLockFile.IsNull()) {
         fDataSetLockFile.Form("%s-dataset-lock", fDataSetDir.Data());
         fDataSetLockFile.ReplaceAll("/","%");
         fDataSetLockFile.ReplaceAll(":","%");
         fDataSetLockFile.Insert(0, TString::Format("%s/", gSystem->TempDirectory()));
      }
      if (!fDataSetLockFile.IsNull() && fIsRemote) {
         TUrl lu(fDataSetLockFile, kTRUE);
         if (!strcmp(lu.GetProtocol(), "file")) {
            // Add host and port
            TUrl u(fDataSetDir);
            TString srv(fDataSetDir);
            srv.Remove(srv.Index(u.GetFile()));
            fDataSetLockFile.Insert(0, srv);
         }
      }
   }

   // Limit in seconds after a lock automatically expires
   fLockFileTimeLimit = 120;

   // Default validity of the cache
   fCacheUpdatePeriod = gEnv->GetValue("ProofDataSet.CacheUpdatePeriod", 0);

   // If the MSS url was not given, check if one is defined via env
   if (fMSSUrl.IsNull())
      fMSSUrl = gEnv->GetValue("ProofDataSet.MSSUrl", "");
   // Default highest priority for xrd-backends
   fStageOpts = gEnv->GetValue("DataSet.StageOpts", "p=3");

   // File to check updates and its locking path
   fListFile.Form("%s/%s", fDataSetDir.Data(), kDataSet_DataSetList);

   // Init the local cache directory if the repository is remote
   fUseCache = kFALSE;
   fLocalCacheDir = "";
   InitLocalCache();
}
////////////////////////////////////////////////////////////////////////////////
/// Init the local cache if required

void TDataSetManagerFile::InitLocalCache()
{
   fUseCache = (fIsRemote) ? kTRUE : kFALSE;

   // Check if the caller has given specific instructions
   TString useCache;
   if (TestBit(TDataSetManager::kUseCache)) useCache = "yes";
   if (TestBit(TDataSetManager::kDoNotUseCache)) useCache = "no";
   if (useCache.IsNull()) useCache = gEnv->GetValue("DataSet.UseCache", "");
   if (useCache.IsNull() && gSystem->Getenv("DATASETCACHE"))
      useCache = gSystem->Getenv("DATASETCACHE");
   useCache.ToLower();
   if (!useCache.IsNull())
      fUseCache = (useCache == "no" || useCache == "0") ? kFALSE : kTRUE;

   if (fUseCache) {
      fLocalCacheDir = gSystem->Getenv("DATASETLOCALCACHEDIR");
      if (fLocalCacheDir.IsNull())
         fLocalCacheDir = gEnv->GetValue("DataSet.LocalCacheDir", "");
      if (!fLocalCacheDir.IsNull()) {
         // Make sure that the non-default local cache directory exists and is writable
         if (gSystem->AccessPathName(fLocalCacheDir)) {
            if (gSystem->mkdir(fLocalCacheDir, kTRUE) != 0) {
               // Switch to default
               Warning("InitLocalCache",
                        "non-default local cache directory '%s' could not be created"
                        " - switching to default", fLocalCacheDir.Data());
               fLocalCacheDir = "";
            }
         }
         if (!fLocalCacheDir.IsNull() &&
               gSystem->AccessPathName(fLocalCacheDir, kWritePermission)) {
            Warning("InitLocalCache",
                     "non-default local cache directory '%s' is not writable"
                     " - switching to default",
                     fDataSetDir.Data());
            fLocalCacheDir = "";
         }
      }
      // If not defined yet try the (unique) default
      if (fLocalCacheDir.IsNull()) {
         // Add something related to fDataSetDir
         TString uds(fDataSetDir.Data());
         uds.ReplaceAll("/","%");
         uds.ReplaceAll(":","%");
         if (TString(gSystem->TempDirectory()).EndsWith(fUser.Data())) {
            fLocalCacheDir.Form("%s/%s/%s", gSystem->TempDirectory(),
                                 kDataSet_LocalCache, uds.Data());
         } else {
            fLocalCacheDir.Form("%s/%s/%s/%s", gSystem->TempDirectory(),
                                 fUser.Data(), kDataSet_LocalCache, uds.Data());
         }
         // Make sure that the local cache dir exists and is writable
         if (gSystem->AccessPathName(fLocalCacheDir) && gSystem->mkdir(fLocalCacheDir, kTRUE) != 0) {
            // Disable
            Warning("InitLocalCache",
                     "local cache directory '%s' could not be created"
                     " - disabling cache", fLocalCacheDir.Data());
            fUseCache = kFALSE;
         }
         if (!fLocalCacheDir.IsNull() &&
             gSystem->AccessPathName(fLocalCacheDir, kWritePermission)) {
            Warning("InitLocalCache",
                     "local cache directory '%s' is not writable - disabling cache",
                     fDataSetDir.Data());
            fUseCache = kFALSE;
         }
         if (!fUseCache) fLocalCacheDir = "";
      }
   }
   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Parse the input string and set the init bits accordingly
/// Format is
///    dir:<datasetdir> [mss:<mss-url>] [opt:<base-options>]
/// The <datasetdir> is mandatory.
/// See TDataSetManager::ParseInitOpts for the available
/// base options.
/// The base options are already initialized by the base constructor

void TDataSetManagerFile::ParseInitOpts(const char *ins)
{
   SetBit(TObject::kInvalidObject);
   fOpenPerms = kFALSE;

   // Needs something in
   if (!ins || strlen(ins) <= 0) return;

   // Extract elements
   Int_t from = 0;
   TString s(ins), tok;
   while (s.Tokenize(tok, from, " ")) {
      if (tok.BeginsWith("dir:"))
         fDataSetDir = tok(4, tok.Length());
      if (tok.BeginsWith("mss:"))
         fMSSUrl = tok(4, tok.Length());
      if (tok == "perms:open")
         fOpenPerms = kTRUE;
   }

   // The directory is mandatory
   if (fDataSetDir.IsNull()) return;

   // Object is valid
   ResetBit(TObject::kInvalidObject);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns path of the indicated dataset. The extension is '.root' for all files
/// except for 'dsName==ls' which have extension '.txt'.
/// If 'local' is kTRUE the local cache path is returned instead in the form
/// <cachedir>/<group>.<user>.<dsName>.<ext>.
/// NB: contains a static TString for result, so copy result before using twice.

const char *TDataSetManagerFile::GetDataSetPath(const char *group,
                                                const char *user,
                                                const char *dsName,
                                                TString &md5path, Bool_t local)
{
   if (fgCommonDataSetTag == group)
     group = fCommonGroup;

   if (fgCommonDataSetTag == user)
     user = fCommonUser;

   const char *ext = (!strcmp(dsName, "ls")) ? ".txt" : ".root";
   static TString result;
   if (local) {
      result.Form("%s/%s.%s.%s%s", fLocalCacheDir.Data(), group, user, dsName, ext);
      md5path.Form("%s/%s.%s.%s.md5sum", fLocalCacheDir.Data(), group, user, dsName);
   } else {
      result.Form("%s/%s/%s/%s%s", fDataSetDir.Data(), group, user, dsName, ext);
      md5path.Form("%s/%s/%s/%s.md5sum", fDataSetDir.Data(), group, user, dsName);
   }
   if (gDebug > 0)
      Info("GetDataSetPath","paths: %s, %s ", result.Data(), md5path.Data());
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Save into the <datasetdir>/kDataSet_DataSetList file the name of the updated
/// or created or modified dataset. For still existing datasets, fill the
/// modification date in seconds anf the checksum.
/// Returns 0 on success, -1 on error

Int_t TDataSetManagerFile::NotifyUpdate(const char *group, const char *user,
                                        const char *dsName, Long_t mtime, const char *checksum)
{
   // Update / create list for the owner
   Long_t lsmtime = 0;
   TString lschecksum;
   Int_t lsrc = -1;
   if ((lsrc = CreateLsFile(group, user, lsmtime, lschecksum)) < 0) {
      Warning("NotifyUpdate", "problems (re-)creating the dataset lists for '/%s/%s'",
                              group, user);
   }

   {  TLockFile lock(fDataSetLockFile, fLockFileTimeLimit);
      TString dspath = TString::Format("/%s/%s/%s", group, user, dsName);
      // Check if the global file exists
      Bool_t hasListFile = gSystem->AccessPathName(fListFile) ? kFALSE : kTRUE;
      // Load the info in form of TMacro
      TMD5 *oldMd5 = 0, *newMd5 = 0;
      if (hasListFile && !(oldMd5 = TMD5::FileChecksum(fListFile.Data()))) {
         Error("NotifyUpdate", "problems calculating old checksum of %s", fListFile.Data());
         return -1;
      }
      // Create the TMacro object, filling it with the existing file, if any
      TMacro mac;
      if (hasListFile) mac.ReadFile(fListFile.Data());
      // Locate the line to change or delete
      TObjString *os = mac.GetLineWith(dspath);
      if (os) {
         // Delete the line if required so
         if (!strcmp(checksum, "removed")) {
            mac.GetListOfLines()->Remove(os);
            SafeDelete(os);
         } else {
            // Update the information
            os->SetString(TString::Format("%ld %s %s", mtime, dspath.Data(), checksum));
         }
      } else {
         if (!strcmp(checksum, "removed")) {
            Warning("NotifyUpdate", "entry for removed dataset '%s' not found!", dspath.Data());
         } else {
            // Add new line
            mac.AddLine(TString::Format("%ld %s %s", mtime, dspath.Data(), checksum));
         }
      }
      // Locate the ls line now, is needed
      TString lspath = TString::Format("/%s/%s/ls", group, user);
      os = mac.GetLineWith(lspath);
      if (os) {
         // Delete the line if required so
         if (lsrc == 1) {
            mac.GetListOfLines()->Remove(os);
            SafeDelete(os);
         } else {
            // Update the information
            os->SetString(TString::Format("%ld %s %s", lsmtime, lspath.Data(), lschecksum.Data()));
         }
      } else {
         if (lsrc == 0) {
            // Add new line
            mac.AddLine(TString::Format("%ld %s %s", lsmtime, lspath.Data(), lschecksum.Data()));
         }
      }
      // Write off the new content
      mac.SaveSource(fListFile.Data());
      if (fOpenPerms) {
         if (gSystem->Chmod(fListFile.Data(), 0666) < 0) {
            Warning("NotifyUpdate",
               "can't set permissions of dataset list file %s (#%d)",
               fListFile.Data(), TSystem::GetErrno());
         }
      }
      if (!(newMd5 = TMD5::FileChecksum(fListFile.Data()))) {
         Error("NotifyUpdate", "problems calculating new checksum of %s", fListFile.Data());
         SafeDelete(oldMd5);
         return -1;
      }
      if (oldMd5 && (*newMd5 == *oldMd5))
         Warning("NotifyUpdate", "checksum for %s did not change!", fListFile.Data());
      // Cleanup
      SafeDelete(oldMd5);
      SafeDelete(newMd5);
   }
   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create or recreate the dataset lists for 'uri'.
/// The list are saved in text form in 'uri'/ls.txt for fast browsing and in
/// 'uri'/ls.root in form of TMacro for optimized and portable transfer.
/// Return 0 on success, 1 if the file was empty, -1 on error

Int_t TDataSetManagerFile::CreateLsFile(const char *group, const char *user,
                                        Long_t &mtime, TString &checksum)
{
   mtime = 0;
   checksum = "";
   // Create temporary file; we cannot lock now because we would (dead-)lock
   // during ShowDataSets
   TString tmpfile;
   tmpfile.Form("%s/%s/%s/ls.tmp.txt", fDataSetDir.Data(), group, user);

   // Redirect output to 'tmpfile'
   RedirectHandle_t rh;
   if (gSystem->RedirectOutput(tmpfile.Data(), "w", &rh) != 0) {
      Error("CreateLsFile", "problems redirecting output to %s (errno: %d)",
                            tmpfile.Data(), TSystem::GetErrno());
      return -1;
   }
   // Create the list
   TString uri;
   uri.Form("/%s/%s", group, user);
   ShowDataSets(uri, "forcescan:noheader:");
   // Restore output to standard streams
   if (gSystem->RedirectOutput(0, 0, &rh) != 0) {
      Error("CreateLsFile", "problems restoring output to standard streams (errno: %d)",
                            TSystem::GetErrno());
      return -1;
   }
   // We can lock now
   TLockFile lock(fDataSetLockFile, fLockFileTimeLimit);
   // Rename the temp file
   TString lsfile;
   lsfile.Form("%s/%s/%s/ls.txt", fDataSetDir.Data(), group, user);
   // Remove the file, if existing
   if (!gSystem->AccessPathName(lsfile) && gSystem->Unlink(lsfile) != 0) {
      Error("CreateLsFile", "problems unlinking old file '%s' (errno: %d)",
                            lsfile.Data(), TSystem::GetErrno());
      return -1;
   }
   // Save the new file only if non empty
   FileStat_t st;
   if (gSystem->GetPathInfo(tmpfile, st) == 0 && st.fSize > 0) {
      if (gSystem->Rename(tmpfile, lsfile) != 0) {
         Error("CreateLsFile", "problems renaming '%s' to '%s' (errno: %d)",
                              tmpfile.Data(), lsfile.Data(), TSystem::GetErrno());
         return -1;
      }
#ifndef WIN32
      // Make sure that the ownership and permissions are those expected
      FileStat_t udirst;
      if (!fIsRemote && gSystem->GetPathInfo(gSystem->GetDirName(tmpfile), udirst) == 0) {
         if (chown(lsfile.Data(), udirst.fUid, udirst.fGid) != 0) {
            Warning("CreateLsFile", "problems setting ownership on file '%s' (errno: %d)",
                                    lsfile.Data(), TSystem::GetErrno());
         }
         if (fOpenPerms) {
            if (gSystem->Chmod(lsfile.Data(), 0666) < 0) {
               Warning("NotifyUpdate",
                  "can't set permissions of list file %s (#%d)",
                  lsfile.Data(), TSystem::GetErrno());
            }
         }
         else if (chmod(lsfile.Data(), 0644) != 0) {
            Warning("CreateLsFile", "problems setting permissions on file '%s' (errno: %d)",
                                    lsfile.Data(), TSystem::GetErrno());
         }
      }
#endif
      mtime = st.fMtime;
      TMD5 *md5 = TMD5::FileChecksum(lsfile);
      if (!md5) {
         Error("CreateLsFile", "problems calculating checksum for '%s'", lsfile.Data());
      } else {
         checksum = md5->AsString();
         SafeDelete(md5);
      }
   } else {
      if (!gSystem->AccessPathName(tmpfile) && gSystem->Unlink(tmpfile) != 0) {
         Error("CreateLsFile", "problems unlinking temp file '%s' (errno: %d)",
                               tmpfile.Data(), TSystem::GetErrno());
         return -1;
      }
      // No datasets anymore
      return 1;
   }
   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Adds the dataset in the folder of group, user to the list in target.
/// If dsName is defined, only the information about the specified dataset
/// is processed.
///
/// The unsigned int 'option' is forwarded to GetDataSet and BrowseDataSet.
/// Available options (to be .or.ed):
///    kPrint          print the dataset content
///    kQuotaUpdate    update quotas
///    kExport         use export naming
///    kList           get a list of dataset names
///
/// NB1: options "kPrint", "kQuoatUpdate" and "kExport" are mutually exclusive
/// NB2: for options "kPrint" and "kQuotaUpdate" return is null.

Bool_t TDataSetManagerFile::BrowseDataSets(const char *group, const char *user,
                                           const char *dsName,
                                           UInt_t option, TObject *target)
{
   TString userDirPath;
   userDirPath.Form("%s/%s/%s", fDataSetDir.Data(), group, user);
   void *userDir = gSystem->OpenDirectory(userDirPath);
   if (!userDir)
     return kFALSE;

   // Options
   Bool_t printing = (option & kPrint) ? kTRUE : kFALSE;
   Bool_t exporting = (option & kExport) ? kTRUE : kFALSE;
   Bool_t updating = (option & kQuotaUpdate) ? kTRUE : kFALSE;
   Bool_t printout = (printing && (option & kDebug)) ? kTRUE : kFALSE;
   Bool_t listing = (option & kList) ? kTRUE : kFALSE;

   // If printing is required add kReadShort to the options
   if (printing || updating)
      option |= kReadShort;

   // The last three options are mutually exclusive
   if (((Int_t)printing + (Int_t)exporting + (Int_t)updating + (Int_t)listing) > 1) {
      Error("BrowseDataSets",
            "only one of kPrint, kQuotaUpdate, kExport or kList can be specified at once");
      return kFALSE;
   }
   Bool_t fillmap = (!exporting && !printing && !updating) ? kTRUE : kFALSE;

   // Output object
   TMap *outmap = (fillmap || exporting || listing) ? (TMap *)target : (TMap *)0;
   TList *outlist = (printing) ? (TList *)target : (TList *)0;

   TRegexp rg("^[^./][^/]*.root$");  //check that it is a root file, not starting with "."

   TRegexp *reds = 0;
   if (dsName && strlen(dsName) > 0) reds = new TRegexp(dsName, kTRUE);

   TMap *userMap = 0, *datasetMap = 0;
   // loop over datasets
   const char *dsEnt = 0;
   while ((dsEnt = gSystem->GetDirEntry(userDir))) {
      TString datasetFile(dsEnt);
      if (datasetFile.Index(rg) != kNPOS) {
         TString datasetName(datasetFile(0, datasetFile.Length()-5));

         // Check dataset name, if required
         if (reds && datasetName.Index(*reds) == kNPOS) continue;

         if (gDebug > 0)
            Info("GetDataSets", "found dataset %s of user %s in group %s",
                                 datasetName.Data(), user, group);

         TFileCollection *fileList = GetDataSet(group, user, datasetName, option);
         if (!fileList) {
            Error("GetDataSets", "dataset %s (user %s, group %s) could not be opened",
                                 datasetName.Data(), user, group);
            continue;
         }
         if (gDebug > 0)
            fileList->Print();

         // We found a dataset, now add it to the map

         // COMMON dataset transition
         const char *mapGroup = group;
         if (fCommonGroup == mapGroup)
            mapGroup = fgCommonDataSetTag.Data();
         const char *mapUser = user;
         if (fCommonUser == mapUser)
            mapUser = fgCommonDataSetTag.Data();

         if (fillmap && !listing && outmap) {
            if (!(userMap = dynamic_cast<TMap*> (outmap->GetValue(mapGroup)))) {
               userMap = new TMap;
               userMap->SetOwner();
               outmap->Add(new TObjString(mapGroup), userMap);
            }

            if (!(datasetMap = dynamic_cast<TMap*> (userMap->GetValue(mapUser)))) {
               datasetMap = new TMap;
               datasetMap->SetOwner();
               userMap->Add(new TObjString(mapUser), datasetMap);
            }
         }

         // Action depends on option
         if (exporting) {

            // Just format the dataset name with group and user
            TString dsNameFormatted(Form("/%s/%s/%s", mapGroup,
                                                      mapUser, datasetName.Data()));
            if (outmap)
               outmap->Add(new TObjString(dsNameFormatted), fileList);

         } else if (updating) {

            // Update quotas
            GetQuota(mapGroup, mapUser, datasetName.Data(), fileList);

         } else if (printing) {

            // Prepare the output list
            if (outlist) {
               TString dsNameFormatted(Form("/%s/%s/%s", mapGroup,
                                                         mapUser, datasetName.Data()));
               // Magic number?
               if (dsNameFormatted.Length() < 42)
                  dsNameFormatted.Resize(42);

               // Done
               outlist->Add(fileList->ExportInfo(dsNameFormatted));
               if (printout) {
                  // Recreating the LS file
                  TObjString *os = (TObjString *) outlist->Last();
                  if (os) Printf("%s", os->GetName());
               }
            }
         } else if (listing) {

            // Just a list of available datasets
            if (outmap) {
               outmap->Add(new TObjString(TString::Format("/%s/%s/%s", mapGroup, mapUser, datasetName.Data())),
                           new TObjString(""));
            }
         } else {
            if (fillmap && datasetMap)
               datasetMap->Add(new TObjString(datasetName), fileList);
         }
      }
   }
   gSystem->FreeDirectory(userDir);
   SafeDelete(reds);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// General purpose call to go through the existing datasets.
/// If <user> is 0 or "*", act on all datasets for the given <group>.
/// If <group> is 0 or "*", act on all datasets.
/// If <dsName> is defined, only the information about the specified dataset
/// is processed.
/// Action depends on option; available options:
///
///    kExport         Return a TMap object containing all the information about
///                    datasets in the form:
///                    { <group>, <map of users> }
///                                     |
///                             { <map of datasets>, <dataset>}
///                    (<dataset> are TFileCollection objects)
///    kShowDefault    as kExport with in addition a default selection including
///                    the datasets from the current user, the ones from the group
///                    and the common ones
///
///    kPrint          print the dataset content; no output is returned
///    kList           get a list of available dataset names
///    kForceScan      Re-open files while processing kPrint (do not use the
///                    pre-processed information)
///    kNoHeaderPrint  Labelling header is not printed
///    kQuotaUpdate    update {group, user} quotas; no output is returned
///
/// NB1: options "kPrint", "kQuoatUpdate" and "kExport" are mutually exclusive
/// NB2: for options "kPrint" and "kQuoatUpdate" return is null.

TMap *TDataSetManagerFile::GetDataSets(const char *group, const char *user,
                                       const char *dsName, UInt_t option)
{
   if (group && fgCommonDataSetTag == group)
     group = fCommonGroup;

   if (user && fgCommonDataSetTag == user)
     user = fCommonUser;

   // Special treatment for the COMMON user
   Bool_t notCommonUser = kTRUE;
   if ((user && fCommonUser == user) &&
       (group && fCommonGroup == group)) notCommonUser = kFALSE;

   // convert * to "nothing"
   if (group && (strcmp(group, "*") == 0 || !group[0]))
      group = 0;
   if (user && (strcmp(user, "*") == 0 || !user[0]))
      user = 0;

   Bool_t printing = (option & kPrint) ? kTRUE : kFALSE;
   Bool_t forcescan = (option & kForceScan) ? kTRUE : kFALSE;
   Bool_t printheader = (option & kNoHeaderPrint) ? kFALSE : kTRUE;
   Bool_t exporting = (option & kExport) ? kTRUE : kFALSE;
   Bool_t updating = (option & kQuotaUpdate) ? kTRUE : kFALSE;
   Bool_t refreshingls = (option & kRefreshLs) ? kTRUE : kFALSE;
   Bool_t listing = (option & kList) ? kTRUE : kFALSE;

   // The last three options are mutually exclusive
   if (((Int_t)printing + (Int_t)exporting + (Int_t)updating + (Int_t)listing) > 1) {
      Error("GetDataSets", "only one of '?P', '?Q', '?E' or '?L' can be specified at once");
      return 0;
   }

   TObject *result = 0;
   if (printing) {
      // The output is a list of strings
      TList *ol = new TList();
      ol->SetOwner();
      result = ol;
   } else if (exporting || !updating || listing) {
      TMap *om = new TMap;
      om->SetOwner();
      result = om;
   }

   if (gDebug > 0)
      Info("GetDataSets", "opening dir %s", fDataSetDir.Data());

   Long_t m;
   TString s;
   if (option & kShowDefault) {
      // Add the common ones
      if (refreshingls) {
         if (CreateLsFile(fCommonGroup, fCommonUser, m, s) != 0)
            Warning("GetDataSets", "problems recreating 'ls' info for {%s,%s}",
                                   fCommonGroup.Data(), fCommonUser.Data());
      } else if (!printing || forcescan || (printing &&
                 FillLsDataSet(fCommonGroup, fCommonUser, dsName, (TList *)result, option) != 0)) {
         BrowseDataSets(fCommonGroup, fCommonUser, dsName, option, result);
      }
      user = 0;
   } else {
      // Fill the information at least once
      if (!notCommonUser) notCommonUser = kTRUE;
   }

   // Fill the information only once
   if (notCommonUser) {
      // group, user defined, no looping needed
      if (user && group && strchr(user, '*') && strchr(group, '*')) {
         if (refreshingls) {
            if (CreateLsFile(group, user, m, s) != 0)
               Warning("GetDataSets", "problems recreating 'ls' info for {%s,%s}",
                                      group, user);
         } else if (!printing || forcescan || (printing &&
                    FillLsDataSet(group, user, dsName, (TList *)result, option) != 0)) {
            BrowseDataSets(group, user, dsName, option, result);
         }
         if (!printing) return (TMap *)result;
      } else {
         TRegexp *reg = (group && strlen(group) > 0) ? new TRegexp(group, kTRUE) : 0;
         TRegexp *reu = (user && strlen(user) > 0) ? new TRegexp(user, kTRUE) : 0;
         // Loop needed, either on the local cache or on the real thing
         if (printing && !forcescan &&
             fUseCache && CheckLocalCache(group, user, 0, option) == 0) {
            // Loop on the local cache
            Int_t from = 0;
            TString locupdate, dsn, grp, usr;
            locupdate.Form("%s/%s", fLocalCacheDir.Data(), kDataSet_DataSetList);
            TMacro uptmac(locupdate);
            TIter nxl(uptmac.GetListOfLines());
            TObjString *os = 0;
            while ((os = (TObjString *) nxl())) {
               if (!(os->GetString().Contains("/ls"))) continue;
               from = 0;
               if (!(os->GetString().Tokenize(dsn, from, " "))) continue;
               if (!(os->GetString().Tokenize(dsn, from, " "))) continue;
               from = 0;
               // Get group and apply filter on group
               if (!(dsn.Tokenize(grp, from, "/")) || (reg && (grp.Index(*reg) == kNPOS))) continue;
               // Get user and apply filter on user
               if (!(dsn.Tokenize(usr, from, "/")) || (reu && (usr.Index(*reu) == kNPOS))) continue;
               // Now get the info
               if (FillLsDataSet(grp, usr, dsName, (TList *)result, option) != 0) {
                  // Use the standard way opening all the files
                  BrowseDataSets(grp, usr, dsName, option, result);
               }
            }
         } else {
            // Loop needed on the real thing
            void *dataSetDir = 0;
            if ((dataSetDir = gSystem->OpenDirectory(fDataSetDir))) {
               // loop over groups
               const char *eg = 0;
               while ((eg = gSystem->GetDirEntry(dataSetDir))) {

                  if (strcmp(eg, ".") == 0 || strcmp(eg, "..") == 0)
                     continue;

                  if (reg && (TString(eg).Index(*reg) == kNPOS))
                     continue;

                  TString groupDirPath;
                  groupDirPath.Form("%s/%s", fDataSetDir.Data(), eg);

                  // Make sure it is a directory
                  FileStat_t dirSt;
                  if (gSystem->GetPathInfo(groupDirPath, dirSt) != 0 || !R_ISDIR(dirSt.fMode))
                     continue;

                  void *groupDir = gSystem->OpenDirectory(groupDirPath);
                  if (!groupDir)
                     continue;

                  // loop over users
                  const char *eu = 0;
                  while ((eu = gSystem->GetDirEntry(groupDir))) {

                     if (strcmp(eu, ".") == 0 || strcmp(eu, "..") == 0)
                        continue;

                     if (reu && (TString(eu).Index(*reu) == kNPOS))
                        continue;

                     // If we have the ls.macro use that
                     if (refreshingls) {
                        if (CreateLsFile(eg, eu, m, s) != 0)
                           Warning("GetDataSets", "problems recreating 'ls' info for {%s,%s}",
                                                  eg, eu);
                     } else if (!printing || forcescan || (printing &&
                                 FillLsDataSet(eg, eu, dsName, (TList *)result, option) != 0)) {
                        // Use the standard way opening all the files
                        BrowseDataSets(eg, eu, dsName, option, result);
                     }
                  }
                  gSystem->FreeDirectory(groupDir);
               }
               gSystem->FreeDirectory(dataSetDir);
            }
         }
         SafeDelete(reg);
         SafeDelete(reu);
      }
   }
   // Print the result, if required
   if (printing) {
      TList *output = (TList *)result;
      output->Sort();
      if (printheader) {
         Printf("Dataset repository: %s", fDataSetDir.Data());
         Printf("Dataset URI                               | # Files | Default tree | # Events |   Disk   | Staged");
      }
      TIter iter4(output);
      TObjString *os = 0;
      while ((os = dynamic_cast<TObjString*> (iter4()))) {
         if (os->GetString().BeginsWith("file:")) {
            // Path of the file to be browsed
            TString path(os->GetString()(5, os->GetString().Length()));
            RedirectHandle_t rh(path.Data());
            gSystem->ShowOutput(&rh);
            fflush(stderr);
         } else {
            // Simple line
            Printf("%s", os->String().Data());
         }
      }
      // Cleanup
      SafeDelete(output);
      result = 0;
   }

   return (TMap *)result;
}

////////////////////////////////////////////////////////////////////////////////
/// Check for the 'ls.txt' for 'group' and 'user' and fill the path for the
/// ls file in 'out'.
/// If 'dsname' is defined, open the file and extract the relevant line.
/// Return 0 on success, -1 on failure

Int_t TDataSetManagerFile::FillLsDataSet(const char *group, const char *user,
                                         const char *dsname, TList *out, UInt_t option)
{
   // Check inputs
   if (!group || strlen(group) <= 0 || !user || strlen(user) <= 0 || !out) {
      Error("FillLsDataSet", "at least one of the inputs is invalid (%s,%s,%p)", group, user, out);
      return -1;
   }

   // File with the TMacro
   Int_t crc = -1;
   TString lsfile, lsmd5file;
   if (!fUseCache || (fUseCache && (crc = CheckLocalCache(group, user, "ls", option)) <= 0)) {
      Bool_t local = (crc == 0) ? kTRUE : kFALSE;
      lsfile = GetDataSetPath(group, user, "ls", lsmd5file, local);
   } else {
      // The dataset does not exist anymore
      return 0;
   }

   if (gSystem->AccessPathName(lsfile, kFileExists)) {
      if (gDebug > 0)
         Info("FillLsDataSet", "file '%s' does not exists", lsfile.Data());
      return -1;
   }
   if (gSystem->AccessPathName(lsfile, kReadPermission)) {
      Warning("FillLsDataSet", "file '%s' exists cannot be read (permission denied)", lsfile.Data());
      return -1;
   }

   if (dsname && strlen(dsname) > 0) {
      // Read the macro
      TMacro *mac = new TMacro(lsfile.Data());
      if (!mac) {
         Error("FillLsDataSet", "could not initialize TMacro from '%s'", lsfile.Data());
         return -1;
      }
      // Prepare the string to search for
      TString fullname = TString::Format("/%s/%s/%s", group, user, dsname);
      Bool_t wc = (fullname.Contains("*")) ? kTRUE : kFALSE;
      if (wc) fullname.ReplaceAll("*", ".*");
      TRegexp reds(fullname);
      TIter nxl(mac->GetListOfLines());
      TObjString *o;
      Int_t nf = 0;
      while ((o = (TObjString *) nxl())) {
         if (o->GetString().Index(reds) != kNPOS) {
            out->Add(o->Clone());
            nf++;
            if (!wc) break;
         }
      }
      if (nf > 0 && gDebug > 0)
         Info("FillLsDataSet", "no match for dataset uri '/%s/%s/%s'", group, user, dsname);
      // Delete the macro
      SafeDelete(mac);
   } else {
      // Fill in the file information
      out->Add(new TObjString(TString::Format("file:%s", lsfile.Data())));
   }
   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Returns the dataset <dsName> of user <user> in group <group> .
/// If checksum is non-zero, it will contain the pointer to a TMD5 sum object
/// with the checksum of the file, has to be deleted by the user.
/// If option has the bi kReadShort set, the shortobject is read, that does not
/// contain the list of files. This is much faster.

TFileCollection *TDataSetManagerFile::GetDataSet(const char *group,
                                                 const char *user,
                                                 const char *dsName,
                                                 UInt_t option,
                                                 TMD5 **checksum)
{
   TFileCollection *fileList = 0;
   Bool_t readshort = (option & kReadShort) ? kTRUE : kFALSE;
   // Check is the file is in the cache
   Int_t crc = -1;
   TString path, md5path;
   if (readshort || !fUseCache ||
      (!readshort && fUseCache && (crc = CheckLocalCache(group, user, dsName, option)) <= 0)) {
      Bool_t local = (crc == 0) ? kTRUE : kFALSE;
      path = GetDataSetPath(group, user, dsName, md5path, local);
   } else {
      // The dataset does not exist (anymore?)
      if (gDebug > 0)
         Info("GetDataSet", "dataset %s does not exist", path.Data());
      return fileList;
   }

   // Now we lock because we are going to use the file, if there
   TLockFile lock(fDataSetLockFile, fLockFileTimeLimit);

   // Check if the file can be opened for reading
   if (gSystem->AccessPathName(path, kFileExists)) {
      if (gDebug > 0)
         Info("GetDataSet", "file '%s' does not exists", path.Data());
      return fileList;
   }
   if (gSystem->AccessPathName(path, kReadPermission)) {
      Warning("GetDataSet", "file '%s' exists cannot be read (permission denied)", path.Data());
      return fileList;
   }

   // Get checksum
   if (checksum) {
      // save md5 sum
      *checksum = TMD5::ReadChecksum(md5path);
      if (!(*checksum)) {
         Error("GetDataSet", "could not get checksum of %s from %s", path.Data(), md5path.Data());
         return fileList;
      }
   }

   TFile *f = TFile::Open(path.Data());
   if (!f) {
      Error("GetDataSet", "could not open file %s", path.Data());
      if (checksum) SafeDelete(*checksum);
      return fileList;
   }

   if (option & kReadShort)
     fileList = dynamic_cast<TFileCollection*> (f->Get("dataset_short"));

   if (!fileList)
     fileList = dynamic_cast<TFileCollection*> (f->Get("dataset"));

   f->Close();
   SafeDelete(f);

   return fileList;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the local cache information for group, user, dsName is up-to-date
/// If not, make the relevant updates
/// Return 0 if OK, 1 if the dataset does not exists anymore, -1 on failure

Int_t TDataSetManagerFile::CheckLocalCache(const char *group, const char *user,
                                           const char *dsName, UInt_t option)
{
   // Check first if the global update info is uptodate
   static TMacro *uptmac = 0;
   Bool_t need_last_update = (option & kNoCacheUpdate) ? kFALSE : kTRUE;
   TString locupdtim, locupdate, remupdate;
   locupdtim.Form("%s/%s.update", fLocalCacheDir.Data(), kDataSet_DataSetList);
   locupdate.Form("%s/%s", fLocalCacheDir.Data(), kDataSet_DataSetList);
   remupdate.Form("%s/%s", fDataSetDir.Data(), kDataSet_DataSetList);
   need_last_update = (gSystem->AccessPathName(locupdate)) ? kTRUE : need_last_update;
   TDatime now;
   UInt_t tnow = now.Convert();
   FileStat_t timst, locst, remst;
   if (need_last_update && !gSystem->AccessPathName(locupdtim)) {
      if (gSystem->GetPathInfo(locupdtim, timst) == 0) {
         need_last_update = kFALSE;
         if ((Int_t)tnow > timst.fMtime + fCacheUpdatePeriod) need_last_update = kTRUE;
      }
   }
   if (need_last_update) {
      if (gSystem->GetPathInfo(remupdate, remst) != 0) {
         Error("CheckLocalCache", "cannot get info for remote file '%s' - ignoring", remupdate.Data());
         return -1;
      }
      if (gSystem->GetPathInfo(locupdate, locst) == 0) {
         need_last_update = kFALSE;
         if (remst.fMtime > locst.fMtime) {
            need_last_update = kTRUE;
         } else {
            if (!gSystem->AccessPathName(locupdtim))
               if (gSystem->Utime(locupdtim, tnow, 0) != 0)
                  Warning("CheckLocalCache",
                          "cannot set modification time on file '%s' (errno: %d)",
                          locupdtim.Data(), TSystem::GetErrno());
         }
      }
   }
   // Get the file, if needed
   if (need_last_update) {
      if (!TFile::Cp(remupdate, locupdate, kFALSE)) {
         Error("CheckLocalCache", "cannot get remote file '%s' - ignoring", remupdate.Data());
         return -1;
      }
      // Set the modification time
      if (gSystem->Utime(locupdate, remst.fMtime, 0) != 0) {
         Warning("CheckLocalCache", "cannot set modification time on file '%s' (errno: %d)",
                           locupdate.Data(), TSystem::GetErrno());
      }
      // Touch or create the file to control updates
      if (gSystem->AccessPathName(locupdtim)) {
         FILE *ftim = fopen(locupdtim.Data(), "w");
         if (!ftim) {
            Warning("CheckLocalCache", "problems create file '%s' (errno: %d)",
                                       locupdtim.Data(), TSystem::GetErrno());
         } else {
            if (fclose(ftim) != 0)
               Warning("CheckLocalCache", "problems close file '%s' (errno: %d)",
                                          locupdtim.Data(), TSystem::GetErrno());
            if (gSystem->Utime(locupdtim, now.Convert(), 0) != 0)
               Warning("CheckLocalCache",
                       "cannot set modification time on file '%s' (errno: %d)",
                       locupdtim.Data(), TSystem::GetErrno());
         }
      }
      // Update macro info
      SafeDelete(uptmac);
      uptmac = new TMacro(locupdate);
   } else {
      // Touch or create the file to control updates
      if (gSystem->AccessPathName(locupdtim)) {
         FILE *ftim = fopen(locupdtim.Data(), "w");
         if (!ftim) {
            Warning("CheckLocalCache", "problems create file '%s' (errno: %d)",
                                       locupdtim.Data(), TSystem::GetErrno());
         } else {
            if (fclose(ftim) != 0)
               Warning("CheckLocalCache", "problems close file '%s' (errno: %d)",
                                          locupdtim.Data(), TSystem::GetErrno());
            if (gSystem->GetPathInfo(locupdate, locst) == 0) {
               if (gSystem->Utime(locupdtim, locst.fMtime, 0) != 0)
                  Warning("CheckLocalCache",
                          "cannot set modification time on file '%s' (errno: %d)",
                          locupdtim.Data(), TSystem::GetErrno());
            } else {
               Warning("CheckLocalCache", "cannot get info for file '%s'"
                       " - will not touch '%s'", locupdate.Data(), locupdtim.Data());
            }
         }
      }
      if (!uptmac) uptmac = new TMacro(locupdate);
   }

   // If we are just interested in the global dataset list we are done
   if (!dsName || strlen(dsName) <= 0)
      return 0;

   // Read the information
   TString ds, locpath, path, locmd5path, md5path, remmd5s;
   TMD5 *locmd5 = 0;
   // The paths ...
   path = GetDataSetPath(group, user, dsName, md5path);
   locpath = GetDataSetPath(group, user, dsName, locmd5path, kTRUE);
   ds.Form("/%s/%s/%s", group, user, dsName);
   TObjString *os = uptmac->GetLineWith(ds);
   if (!os) {
      // DataSet does not exist anymore
      if (strcmp(dsName, "ls"))
         Warning("CheckLocalCache", "dataset '%s' does not exists anymore", ds.Data());
      return 1;
   }
   // Extract the relevant information
   TString s;
   Int_t from = 0;
   while (os->GetString().Tokenize(s, from, " ")) {
      if (!s.IsDigit() && s != ds) {
         remmd5s = s;
      }
   }
   if (remmd5s == "---") {
      // DataSet does not exist anymore
      if (strcmp(dsName, "ls"))
         Warning("CheckLocalCache", "dataset '%s' does not exists anymore", ds.Data());
      return 1;
   }
   Bool_t need_update = (option & kNoCacheUpdate) ? kFALSE : kTRUE;
   if (!gSystem->AccessPathName(locpath)) {
      if (need_update) {
         need_update = kFALSE;
         locmd5 = TMD5::ReadChecksum(locmd5path);
         if (!locmd5 && !(locmd5 = TMD5::FileChecksum(locpath))) {
            Warning("CheckLocalCache", "cannot get checksum of '%s' - assuming match failed", ds.Data());
            need_update = kTRUE;
         } else {
            if (remmd5s != locmd5->AsString()) need_update = kTRUE;
         }
      }
   } else {
      need_update = kTRUE;
   }
   // Get the file, if needed
   if (need_update) {
      SafeDelete(locmd5);
      if (!TFile::Cp(path, locpath, kFALSE)) {
         Error("CheckLocalCache", "cannot get remote file '%s' - ignoring", path.Data());
         return -1;
      }
      // Calculate and save the new checksum
      locmd5 = TMD5::FileChecksum(locpath);
      if (locmd5) {
         if (remmd5s != locmd5->AsString())
            Warning("CheckLocalCache", "checksum for freshly downloaded file '%s' does not match the"
                               " one posted in '%s'", locpath.Data(), kDataSet_DataSetList);
         if (TMD5::WriteChecksum(locmd5path, locmd5) != 0)
            Warning("CheckLocalCache", "problems saving checksum to '%s' (errno: %d)",
                               locmd5path.Data(), TSystem::GetErrno());
      } else {
         Warning("CheckLocalCache", "problems calculating checksum for '%s'", locpath.Data());
      }
   }
   SafeDelete(locmd5);
   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Clear cached information matching uri

Int_t TDataSetManagerFile::ClearCache(const char *uri)
{
   // Open the top directory
   void *dirp = gSystem->OpenDirectory(fLocalCacheDir.Data());
   if (!dirp) {
      Error("ClearCache", "cannot open directory '%s' (errno: %d)",
                          fLocalCacheDir.Data(), TSystem::GetErrno());
      return -1;
   }
   TRegexp *re = 0;
   if (uri && strlen(uri) > 0) {
      if (strcmp(uri, "*") && strcmp(uri, "/*") && strcmp(uri, "/*/") &&
          strcmp(uri, "/*/*") && strcmp(uri, "/*/*/") && strcmp(uri, "/*/*/*")) {
         TString u(uri);
         // Remove leading '/'
         if (u(0) == '/') u.Remove(0,1);
         // Change '/' to '%'
         u.ReplaceAll("/", ".");
         // Init the regular expression
         u.ReplaceAll("*", ".*");
         re = new TRegexp(u.Data());
      }
   }

   Printf(" Dataset repository: %s", fDataSetDir.Data());
   Printf(" Local cache directory: %s", fLocalCacheDir.Data());

   Long64_t totsz = 0, nf = 0;
   FileStat_t st;
   TString path;
   const char *e = 0;
   while ((e = gSystem->GetDirEntry(dirp))) {
      // Skip basic entries
      if (!strcmp(e,".") || !strcmp(e,"..")) continue;
      // Apply regular expression, if requested
      if (re && TString(e).Index(*re) == kNPOS) continue;
      // Group directory
      path.Form("%s/%s", fLocalCacheDir.Data(), e);
      // Get file information
      if (gSystem->GetPathInfo(path, st) != 0) {
         Warning("ShowCache", "problems 'stat'-ing '%s' (errno: %d)",
                               path.Data(), TSystem::GetErrno());
         continue;
      }
      // Count
      totsz += st.fSize;
      nf++;
      // Remove the file
      if (gSystem->Unlink(path) != 0) {
         Warning("ClearCache", "problems unlinking '%s' (errno: %d)",
                               path.Data(), TSystem::GetErrno());
      }
   }
   gSystem->FreeDirectory(dirp);
   SafeDelete(re);

   // Notify totals
   Printf(" %lld bytes (%lld files) have been freed", totsz, nf);

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Show cached information matching uri

Int_t TDataSetManagerFile::ShowCache(const char *uri)
{
   // Open the top directory
   void *dirp = gSystem->OpenDirectory(fLocalCacheDir.Data());
   if (!dirp) {
      Error("ShowCache", "cannot open directory '%s' (errno: %d)",
                         fLocalCacheDir.Data(), TSystem::GetErrno());
      return -1;
   }
   TRegexp *re = 0;
   if (uri && strlen(uri) > 0) {
      if (strcmp(uri, "*") && strcmp(uri, "/*") && strcmp(uri, "/*/") &&
          strcmp(uri, "/*/*") && strcmp(uri, "/*/*/") && strcmp(uri, "/*/*/*")) {
         TString u(uri);
         // Remove leading '/'
         if (u(0) == '/') u.Remove(0,1);
         // Change '/' to '%'
         u.ReplaceAll("/", ".");
         // Init the regular expression
         u.ReplaceAll("*", ".*");
         re = new TRegexp(u.Data());
      }
   }

   Printf(" Dataset repository: %s", fDataSetDir.Data());
   Printf(" Local cache directory: %s", fLocalCacheDir.Data());
   Printf(" Last modified        Size(bytes)  File");

   Long64_t totsz = 0, nf = 0;
   FileStat_t st;
   TString path, sz;
   const char *e = 0;
   while ((e = gSystem->GetDirEntry(dirp))) {
      // Skip basic entries
      if (!strcmp(e,".") || !strcmp(e,"..")) continue;
      // Apply regular expression, if requested
      if (re && TString(e).Index(*re) == kNPOS) continue;
      // Group directory
      path.Form("%s/%s", fLocalCacheDir.Data(), e);
      // Get file information
      if (gSystem->GetPathInfo(path, st) != 0) {
         Warning("ShowCache", "problems 'stat'-ing '%s' (errno: %d)",
                               path.Data(), TSystem::GetErrno());
         continue;
      }
      // Count
      totsz += st.fSize;
      nf++;
      // Get modification time in human readable form
      TDatime tmod(st.fMtime);
      sz.Form("%lld", st.fSize);
      sz.Resize(12);
      Printf(" %s  %s %s", tmod.AsSQLString(), sz.Data(), e);
   }
   gSystem->FreeDirectory(dirp);
   SafeDelete(re);

   // Notify totals
   Printf(" %lld files, %lld bytes", nf, totsz);

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Writes indicated dataset.
/// If option has the bit kFileMustExist set, the file must still exist,
/// otherwise the new dataset is not written (returns 3 in this case).
/// If checksum is non-zero the files current checksum is checked against it,
/// if it does not match the file is not written (the function returns 2 in this
/// case, if the file has disappeared it is also not written (i.e. checksum
/// implies the bit kFileMustExist set in option).
/// Returns != 0 for success, 0 for error

Int_t TDataSetManagerFile::WriteDataSet(const char *group, const char *user,
                                        const char *dsName, TFileCollection *dataset,
                                        UInt_t option, TMD5 *checksum)
{
   TString md5path, path, md5sum;
   Long_t mtime = 0;
   {  TLockFile lock(fDataSetLockFile, fLockFileTimeLimit);

      Bool_t checkIfExists = ((option & kFileMustExist) || checksum) ? kTRUE : kFALSE;

      path = GetDataSetPath(group, user, dsName, md5path);

      if (checkIfExists) {
         // check if file still exists, otherwise it was deleted in the meanwhile and is not written here
         Long_t tmp;
         if (gSystem->GetPathInfo(path, 0, (Long_t*) 0, 0, &tmp) != 0) {
            if (gDebug > 0)
               Info("WriteDataSet", "Dataset disappeared. Discarding update.");
            return 3;
         }
      }

      if (checksum) {
         // verify md5 sum, otherwise the file was changed in the meanwhile and is not overwritten here
         TMD5 *checksum2 = TMD5::FileChecksum(path);
         if (!checksum2) {
            Error("WriteDataSet", "Could not get checksum of %s", path.Data());
            return 0;
         }

         Bool_t checksumAgrees = (*checksum == *checksum2);
         delete checksum2;

         if (!checksumAgrees) {
            if (gDebug > 0)
               Info("WriteDataSet", "Dataset changed. Discarding update.");
            return 2;
         }
      }

      // write first in ".<file>" then rename to recover from crash during writing
      TString tempFile(path);
      Int_t index = -1;
      while (tempFile.Index("/", index+1) >= 0)
         index = tempFile.Index("/", index+1);

      tempFile.Insert(index+1, ".");

      TFile *f = TFile::Open(tempFile, "RECREATE");
      if (!f) {
         Error("WriteDataSet", "Could not open dataset for writing %s", tempFile.Data());
         return 0;
      }

      // write full TFileCollection
      dataset->Write("dataset", TObject::kSingleKey | TObject::kOverwrite);

      // write only metadata
      THashList *list = dataset->GetList();
      dataset->SetList(0);
      dataset->Write("dataset_short", TObject::kSingleKey | TObject::kOverwrite);

      f->Close();
      delete f;

      // Restore full list
      dataset->SetList(list);

      // file is written, rename to real filename
      if (gSystem->Rename(tempFile, path) != 0) {
         Error("WriteDataSet", "renaming %s to %s failed; dataset might be corrupted",
                              tempFile.Data(), path.Data());
         // Cleanup any MD5 sum information
         if (!gSystem->AccessPathName(md5path, kWritePermission) && gSystem->Unlink(md5path) != 0)
            Error("WriteDataSet", "unlink of %s failed", md5path.Data());
         return 0;
      }
      else if (fOpenPerms) {
         if (gSystem->Chmod(path.Data(), 0666) < 0) {
            Warning("NotifyUpdate",
               "can't set permissions of dataset file %s (#%d)",
               path.Data(), TSystem::GetErrno());
         }
      }

      // Save md5 sum, otherwise the file was changed in the meanwhile and is not overwritten here
      if (ChecksumDataSet(path, md5path, md5sum) != 0) {
         Error("WriteDataSet", "problems calculating checksum of %s", path.Data());
         return 0;
      }
      else if (fOpenPerms) {
         if (gSystem->Chmod(md5path.Data(), 0666) < 0) {
            Warning("NotifyUpdate",
               "can't set permissions of dataset MD5 checksum file %s (#%d)",
               md5path.Data(), TSystem::GetErrno());
         }
      }

      FileStat_t st;
      if (gSystem->GetPathInfo(path, st) != 0) {
         Error("WriteDataSet", "could not 'stat' the version of '%s'!", path.Data());
         return 0;
      }
      mtime= st.fMtime;
   }

   // The repository was updated
   if (NotifyUpdate(group, user, dsName, mtime, md5sum) != 0)
      Warning("WriteDataSet", "problems notifying update with 'NotifyUpdate'");

   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the checksum of the indicated dataset at 'path' and save it to the
/// appropriate file 'md5path'. The MD5 string is returned in 'md5sum'.
/// Return 0 on success, -1 on error.

Int_t TDataSetManagerFile::ChecksumDataSet(const char *path,
                                           const char *md5path, TString &checksum)
{
   checksum = "";
   // Check inputs
   if (!path || strlen(path) <= 0 || !md5path || strlen(md5path) <= 0) {
      Error("ChecksumDataSet", "one or more inputs are invalid ('%s','%s')",
             path, md5path);
      return -1;
   }
   // Calculate md5 sum
   TMD5 *md5sum = TMD5::FileChecksum(path);
   if (!md5sum) {
      Error("ChecksumDataSet", "problems calculating checksum of '%s'", path);
      return -1;
   }
   // Save it to a file
   if (TMD5::WriteChecksum(md5path, md5sum) != 0) {
      Error("ChecksumDataSet", "problems saving checksum to '%s'", md5path);
      SafeDelete(md5sum);
      return -1;
   }
   // Fill output
   checksum = md5sum->AsString();
   // Done
   SafeDelete(md5sum);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Removes the indicated dataset

Bool_t TDataSetManagerFile::RemoveDataSet(const char *group, const char *user,
                                               const char *dsName)
{
   TString md5path, path;
   {  TLockFile lock(fDataSetLockFile, fLockFileTimeLimit);

      path = GetDataSetPath(group, user, dsName, md5path);
      Int_t rc = 0;
      // Remove the main file
      if ((rc = gSystem->Unlink(path)) != 0)
         Warning("RemoveDataSet", "problems removing main file '%s' (errno: %d)",
                                 path.Data(), TSystem::GetErrno());
      // Remove the checksum file
      if (gSystem->Unlink(md5path) != 0)
         Warning("RemoveDataSet", "problems removing chcksum file '%s' (errno: %d)",
                                 md5path.Data(), TSystem::GetErrno());
   }

   // The repository was updated
   if (gSystem->AccessPathName(path, kFileExists)) {
      if (NotifyUpdate(group, user, dsName, 0, "removed") != 0)
         Warning("RemoveDataSet", "problems notifying update with 'NotifyUpdate'");
      // Success
      return kTRUE;
   }
   // Failure
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the indicated dataset exits

Bool_t TDataSetManagerFile::ExistsDataSet(const char *group, const char *user,
                                               const char *dsName)
{
   TLockFile lock(fDataSetLockFile, fLockFileTimeLimit);

   TString md5path, path(GetDataSetPath(group, user, dsName, md5path));

   return (gSystem->AccessPathName(path) == kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Register a dataset, perfoming quota checkings and verification, if required.
/// If a dataset with the same name already exists the action fails unless 'opts'
/// contains 'O', in which case the old dataset is overwritten, or contains 'U',
/// in which case 'newDataSet' is added to the existing dataset (duplications are
/// ignored, if any).
/// If 'opts' contains 'V' the dataset files are also verified (if the dataset manager
/// is configured to allow so). By default the dataset is not verified.
/// If 'opts' contains 'T' the in the dataset object (status bits, meta,...)
/// is trusted, i.e. not reset (if the dataset manager is configured to allow so).
/// Returns 0 on success, -1 on failure

Int_t TDataSetManagerFile::RegisterDataSet(const char *uri,
                                           TFileCollection *newDataSet,
                                           const char *opts)
{
   if (!TestBit(TDataSetManager::kAllowRegister))
      return -1;

   // Get the dataset name
   TString dsName;
   if (ParseUri(uri, 0, 0, &dsName, 0, kTRUE) == kFALSE) {
      Error("RegisterDataSet", "problem parsing uri: %s", uri);
      return -1;
   }

   // The dataset
   TFileCollection *dataSet = newDataSet;
   // Check option
   TString opt(opts);
   // If in update mode, retrieve the existing dataset, if any
   if (opt.Contains("U", TString::kIgnoreCase)) {
      // Fail if it exists already
      if (ExistsDataSet(fGroup, fUser, dsName)) {
         // Retrieve the dataset
         if (!(dataSet = GetDataSet(fGroup, fUser, dsName))) {
            // Dataset name does exist
            Warning("RegisterDataSet",
                    "dataset '%s' claimed to exists but retrieval failed - ignoring", uri);
            dataSet = newDataSet;
         } else {
            // Add new dataset to existing one
            dataSet->Add(newDataSet);
         }
      }
   } else if (!opt.Contains("O", TString::kIgnoreCase)) {
      // Fail if it exists already
      if (ExistsDataSet(fGroup, fUser, dsName)) {
         //Dataset name does exist
         Error("RegisterDataSet", "dataset '%s' exists already", uri);
         return -1;
      }
   }

   // A temporary list to hold the unique members (i.e. the very set)
   TList *uniqueFileList = new TList();
   TIter nextFile(dataSet->GetList());
   TFileInfo *prevFile = (TFileInfo*)nextFile();
   uniqueFileList->Add(prevFile);
   while (TFileInfo *obj = (TFileInfo*)nextFile()) {
      // Add entities only once to the temporary list
      if (!uniqueFileList->FindObject(obj->GetFirstUrl()->GetUrl()))
         uniqueFileList->Add(obj);
   }

   // Clear dataSet and add contents of uniqueFileList needed otherwise
   // THashList deletes the objects even when nodelete is set
   dataSet->GetList()->SetOwner(0);
   dataSet->GetList()->Clear("nodelete");
   dataSet->GetList()->SetOwner(1);
   dataSet->GetList()->AddAll(uniqueFileList);
   uniqueFileList->SetOwner(kFALSE);
   delete uniqueFileList;

   // Enforce certain settings
   Bool_t reset = kTRUE;
   if (opt.Contains("T", TString::kIgnoreCase)) {
      if (!TestBit(TDataSetManager::kTrustInfo)) {
         Warning("RegisterDataSet", "configured to not trust the information"
                                    " provided by users: ignoring request");
      } else {
         reset = kFALSE;
      }
   }
   if (reset) {
      dataSet->SetName(dsName);
      dataSet->ResetBitAll(TFileInfo::kStaged);
      dataSet->ResetBitAll(TFileInfo::kCorrupted);
      dataSet->RemoveMetaData();
   }

   // Verify the dataset if required
   if (opt.Contains("V", TString::kIgnoreCase)) {
      if (TestBit(TDataSetManager::kAllowVerify)) {
         // Reopen files and notify
         if (TDataSetManager::ScanDataSet(dataSet, 1, 0, 0, kTRUE ) < 0) {
            Error("RegisterDataSet", "problems verifying the dataset");
            return -1;
         }
      } else {
         Warning("RegisterDataSet", "user-driven verification not allowed: ignoring request");
      }
   }

   // Update accumulated information
   dataSet->Update(fAvgFileSize);

   if (TestBit(TDataSetManager::kCheckQuota)) {
      if (dataSet->GetTotalSize() <= 0) {
         Error("RegisterDataSet", "datasets without size information are not accepted:");
         if (fAvgFileSize < 0) {
            Error("RegisterDataSet", "you may want to define an average"
                                     " file size to get an estimated dataset size");
         }
         return -1;
      }
      // now check the quota
      UpdateUsedSpace();
      Long64_t used = GetGroupUsed(fGroup) + dataSet->GetTotalSize();

      Info("RegisterDataSet", "your group %s uses %.1f GB + %.1f GB for the new dataset; "
                              "the available quota is %.1f GB", fGroup.Data(),
                              (Float_t) GetGroupUsed(fGroup)    / 1073741824,
                              (Float_t) dataSet->GetTotalSize() / 1073741824,
                              (Float_t) GetGroupQuota(fGroup)   / 1073741824);
      if (used > GetGroupQuota(fGroup)) {
         Error("RegisterDataSet", "quota exceeded");
         return -1;
      }
   }

   Bool_t success = WriteDataSet(fGroup, fUser, dsName, dataSet);
   if (!success)
      Error("RegisterDataSet", "could not write dataset: %s", dsName.Data());

   // Done
   return ((success) ? 0 : -1);
}
////////////////////////////////////////////////////////////////////////////////
/// Scans the dataset indicated by <uri> and returns the number of missing files.
/// Returns -1 if any failure occurs, >= 0 on success.
/// For more details, see documentation of
/// ScanDataSet(TFileCollection *dataset, const char *option)

Int_t TDataSetManagerFile::ScanDataSet(const char *uri, UInt_t opt)
{
   TString dsName, dsTree;
   if ((opt & kSetDefaultTree)) {
      if (TestBit(TDataSetManager::kAllowRegister)) {
         if (ParseUri(uri, 0, 0, &dsName, &dsTree, kTRUE)) {
            TFileCollection *dataset = GetDataSet(fGroup, fUser, dsName);
            if (!dataset) return -1;
            dataset->SetDefaultTreeName(dsTree.Data());
            Int_t rc = WriteDataSet(fGroup, fUser, dsName, dataset);
            delete dataset;
            return (rc == 0) ? -1 : 0;
         }
      }
   } else {
      if (TestBit(TDataSetManager::kAllowVerify)) {
         if (ParseUri(uri, 0, 0, &dsName, 0, kTRUE, kTRUE)) {
            if (!(dsName.Contains("*"))) {
               if (ScanDataSet(fGroup, fUser, dsName, opt) > 0)
                  return GetNDisapparedFiles();
            } else {
               TString luri = TString::Format("/%s/%s/%s", fGroup.Data(), fUser.Data(), dsName.Data());
               TMap *fcs = GetDataSets(luri, kList);
               if (!fcs) return -1;
               fcs->Print();
               Int_t ndisappeared = 0;
               TIter nxd(fcs);
               TObjString *d = 0;
               while ((d = (TObjString *) nxd())) {
                  if (!(d->GetString().IsNull())) {
                     TString dsn(d->GetName());
                     if (dsn.Contains("/")) dsn.Remove(0, dsn.Last('/') + 1);
                     if (ScanDataSet(fGroup, fUser, dsn, opt) > 0) {
                        ndisappeared += GetNDisapparedFiles();
                     } else {
                        Warning("ScanDataSet", "problems processing dataset: %s", d->GetName());
                     }
                  } else {
                     Warning("ScanDataSet", "empty string found in map while processing: %s", uri);
                  }
               }
               SafeDelete(fcs);
               return ndisappeared;
            }
         }
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// See documentation of ScanDataSet(TFileCollection *dataset, UInt_t option)

Int_t TDataSetManagerFile::ScanDataSet(const char *group, const char *user,
                                       const char *dsName, UInt_t option)
{
   if (!TestBit(TDataSetManager::kAllowVerify))
      return -1;

   TFileCollection *dataset = GetDataSet(group, user, dsName);
   if (!dataset)
      return -1;

   // File selection
   Int_t fopt = ((option & kAllFiles)) ? -1 : 0;
   if (fopt >= 0) {
      if ((option & kStagedFiles)) {
         fopt = 10;
      } else {
         if ((option & kReopen)) fopt++;
         if ((option & kTouch)) fopt++;
      }
      if ((option & kNoStagedCheck)) fopt += 100;
   } else {
      if ((option & kStagedFiles) || (option & kReopen) || (option & kTouch)) {
         Warning("ScanDataSet", "kAllFiles mode: ignoring kStagedFiles or kReopen"
                                " or kTouch requests");
      }
      if ((option & kNoStagedCheck)) fopt -= 100;
   }

   // Type of action
   Int_t sopt = ((option & kNoAction)) ? -1 : 0;
   if (sopt >= 0) {
      if ((option & kLocateOnly) && (option & kStageOnly)) {
         Error("ScanDataSet", "kLocateOnly and kStageOnly cannot be processed concurrently");
         return -1;
      }
      if ((option & kLocateOnly)) sopt = 1;
      if ((option & kStageOnly)) sopt = 2;
   } else if ((option & kLocateOnly) || (option & kStageOnly)) {
      Warning("ScanDataSet", "kNoAction mode: ignoring kLocateOnly or kStageOnly requests");
   }

   Bool_t dbg = ((option & kDebug)) ? kTRUE : kFALSE;
   // Do the scan
   Int_t result = TDataSetManager::ScanDataSet(dataset, fopt, sopt, 0, dbg,
                                   &fNTouchedFiles, &fNOpenedFiles, &fNDisappearedFiles,
                                   (TList *)0, fAvgFileSize, fMSSUrl.Data(), -1, fStageOpts.Data());
   if (result == 2) {
      if (WriteDataSet(group, user, dsName, dataset) == 0) {
         delete dataset;
         return -2;
      }
   }
   delete dataset;

   return result;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Returns all datasets for the <group> and <user> specified by <uri>.
/// If <user> is 0, it returns all datasets for the given <group>.
/// If <group> is 0, it returns all datasets.
/// The returned TMap contains:
///    <group> --> <map of users> --> <map of datasets> --> <dataset> (TFileCollection)
///
/// The unsigned int 'option' is forwarded to GetDataSet and BrowseDataSet.
/// Available options (to be .or.ed):
///    kShowDefault    a default selection is shown that include the ones from
///                    the current user, the ones from the group and the common ones
///    kPrint          print the dataset content
///    kQuotaUpdate    update quotas
///    kExport         use export naming
///
/// NB1: options "kPrint", "kQuoatUpdate" and "kExport" are mutually exclusive
/// NB2: for options "kPrint" and "kQuoatUpdate" return is null.

TMap *TDataSetManagerFile::GetDataSets(const char *uri, UInt_t option)
{
   TString dsUser, dsGroup, dsName;

   if (((option & kPrint) || (option & kExport)) && strlen(uri) <= 0)
      option |= kShowDefault;

   if (ParseUri(uri, &dsGroup, &dsUser, &dsName, 0, kFALSE, kTRUE))
      return GetDataSets(dsGroup, dsUser, dsName, option);
   return (TMap *)0;
}

////////////////////////////////////////////////////////////////////////////////
/// Utility function used in various methods for user dataset upload.

TFileCollection *TDataSetManagerFile::GetDataSet(const char *uri, const char *opts)
{
   TString dsUser, dsGroup, dsName, ss(opts);

   TFileCollection *fc = 0;
   if (!strchr(uri, '*')) {
      if (!ParseUri(uri, &dsGroup, &dsUser, &dsName)) return fc;
      UInt_t opt = (ss.Contains("S:") || ss.Contains("short:")) ? kReadShort : 0;
      ss.ReplaceAll("S:","");
      ss.ReplaceAll("short:","");
      fc = GetDataSet(dsGroup, dsUser, dsName, opt);
   } else {
      TMap *fcs = GetDataSets(uri);
      if (!fcs) return fc;
      TIter nxd(fcs);
      TObject *k = 0;
      TFileCollection *xfc = 0;
      while ((k = nxd()) && (xfc = (TFileCollection *) fcs->GetValue(k))) {
         if (!fc) {
            // The first one
            fc = xfc;
            fcs->Remove(k);
         } else {
            // Add
            fc->Add(xfc);
         }
      }
   }

   if (fc && !ss.IsNull()) {
      // Build up the subset
      TFileCollection *sfc = 0;
      TString s;
      Int_t from = 0;
      while (ss.Tokenize(s, from, ",")) {
         TFileCollection *xfc = fc->GetFilesOnServer(s.Data());
         if (xfc) {
            if (sfc) {
               sfc->Add(xfc);
               delete xfc;
            } else {
               sfc = xfc;
            }
         }
      }
      // Cleanup
      delete fc;
      fc = sfc;
   }
   // Done
   return fc;
}

////////////////////////////////////////////////////////////////////////////////
/// Removes the indicated dataset

Bool_t TDataSetManagerFile::RemoveDataSet(const char *uri)
{
   TString dsName;

   if (TestBit(TDataSetManager::kAllowRegister)) {
      if (ParseUri(uri, 0, 0, &dsName, 0, kTRUE)) {
         Bool_t rc = RemoveDataSet(fGroup, fUser, dsName);
         if (rc) return kTRUE;
         Error("RemoveDataSet", "error removing dataset %s", dsName.Data());
      }
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the indicated dataset exits

Bool_t TDataSetManagerFile::ExistsDataSet(const char *uri)
{
   TString dsUser, dsGroup, dsName;

   if (ParseUri(uri, &dsGroup, &dsUser, &dsName))
      return ExistsDataSet(dsGroup, dsUser, dsName);
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// updates the used space maps

void TDataSetManagerFile::UpdateUsedSpace()
{
   // Clear used space entries
   fGroupUsed.DeleteAll();
   fUserUsed.DeleteAll();

   // Scan the existing datasets
   GetDataSets(0, 0, 0, (UInt_t)kQuotaUpdate);
}

////////////////////////////////////////////////////////////////////////////////
/// Gets last dataset modification time. Returns -1 on error, or number of
/// seconds since epoch on success

Long_t TDataSetManagerFile::GetModTime(const char *uri)
{
   TString group, user, name, md5path;
   if (!ParseUri(uri, &group, &user, &name)) {
      return -1;
   }

   TString path( GetDataSetPath(group, user, name, md5path) );

   Long_t modTime;
   if (gSystem->GetPathInfo(path.Data(),
      (Long_t *)0, (Long_t *)0, (Long_t *)0, &modTime)) {
      return -1;
   }

   return modTime;
}
