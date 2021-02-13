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
// TDataSetManager                                                 //
//                                                                      //
// This class contains functions to handle datasets in PROOF            //
// It is the layer between TProofServ and the file system that stores   //
// the datasets.                                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TDataSetManager.h"

#include "Riostream.h"

#include "TEnv.h"
#include "TError.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TFileInfo.h"
#include "TFileStager.h"
#include "TMD5.h"
#include "THashList.h"
#include "TKey.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TPRegexp.h"
#include "TRegexp.h"
#include "TSystem.h"
#include "TTree.h"
#include "TUrl.h"
#include "TVirtualMonitoring.h"

// One Gigabyte
#define DSM_ONE_GB (1073741824)

// Name for common datasets
TString TDataSetManager::fgCommonDataSetTag = "COMMON";
TList   *TDataSetManager::fgDataSetSrvMaps = 0;

ClassImp(TDataSetManager);

////////////////////////////////////////////////////////////////////////////////
///
/// Main constructor

TDataSetManager::TDataSetManager(const char *group, const char *user,
                                           const char *options)
                     : fGroup(group),
                       fUser(user), fCommonUser(), fCommonGroup(),
                       fGroupQuota(), fGroupUsed(),
                       fUserUsed(), fNTouchedFiles(0), fNOpenedFiles(0),
                       fNDisappearedFiles(0), fMTimeGroupConfig(-1)
{
   // Fill default group and user if none is given
   if (fGroup.IsNull())
      fGroup = "default";
   if (fUser.IsNull()) {
      fUser = "--nouser--";
      // Get user logon name
      UserGroup_t *pw = gSystem->GetUserInfo();
      if (pw) {
         fUser = pw->fUser;
         delete pw;
      }
   }

   fGroupQuota.SetOwner();
   fGroupUsed.SetOwner();
   fUserUsed.SetOwner();

   fCommonUser = "COMMON";
   fCommonGroup = "COMMON";

   fNTouchedFiles = -1;
   fNOpenedFiles = -1;
   fNDisappearedFiles = -1;
   fMTimeGroupConfig = -1;

   fAvgFileSize = 50000000;  // Default 50 MB per file

   // Parse options
   ParseInitOpts(options);

   if (!fUser.IsNull() && !fGroup.IsNull()) {

      // If not in sandbox, construct the base URI using session defaults
      // (group, user) (syntax: /group/user/dsname[#[subdir/]objname])
      if (!TestBit(TDataSetManager::kIsSandbox))
         fBase.SetUri(TString(Form("/%s/%s/", fGroup.Data(), fUser.Data())));

   }

   // List of dataset server mapping instructions
   TString srvmaps(gEnv->GetValue("DataSet.SrvMaps",""));
   TString srvmapsenv(gSystem->Getenv("DATASETSRVMAPS"));
   if (!(srvmapsenv.IsNull())) {
      if (srvmapsenv.BeginsWith("+")) {
         if (!(srvmaps.IsNull())) srvmaps += ",";
         srvmaps += srvmapsenv(1,srvmapsenv.Length());
      } else {
         srvmaps = srvmapsenv;
      }
   }
   if (!(srvmaps.IsNull()) && !(fgDataSetSrvMaps = ParseDataSetSrvMaps(srvmaps)))
      Warning("TDataSetManager", "problems parsing DataSet.SrvMaps input info (%s)"
                                 " - ignoring", srvmaps.Data());

   // Read config file
   ReadGroupConfig(gEnv->GetValue("Proof.GroupFile", ""));
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TDataSetManager::~TDataSetManager()
{
   // Clear used space
   fGroupQuota.DeleteAll();
   fGroupUsed.DeleteAll();
   fUserUsed.DeleteAll();
}

////////////////////////////////////////////////////////////////////////////////
/// Parse the opts string and set the init bits accordingly
/// Available options:
///    Cq:               set kCheckQuota
///    Ar:               set kAllowRegister
///    Av:               set kAllowVerify
///    Ti:               set kTrustInfo
///    Sb:               set kIsSandbox
///    Ca:               set kUseCache or kDoNotUseCache
/// The opts string may also contain additional unrelated info: in such a case
/// the field delimited by the prefix "opt:" is analyzed, e.g. if opts is
/// "/tmp/dataset  opt:Cq:-Ar: root://lxb6046.cern.ch" only the substring
/// "Cq:-Ar:" will be parsed .

void TDataSetManager::ParseInitOpts(const char *opts)
{
   // Default option bits
   ResetBit(TDataSetManager::kCheckQuota);
   SetBit(TDataSetManager::kAllowRegister);
   SetBit(TDataSetManager::kAllowVerify);
   SetBit(TDataSetManager::kTrustInfo);
   ResetBit(TDataSetManager::kIsSandbox);
   ResetBit(TDataSetManager::kUseCache);
   ResetBit(TDataSetManager::kDoNotUseCache);

   if (opts && strlen(opts) > 0) {
      TString opt(opts);
      // If it contains the prefix "opt:", isolate the related field
      Int_t ip = opt.Index("opt:");
      if (ip != kNPOS) opt.Remove(0, ip + 4);
      ip = opt.Index(" ");
      if (ip != kNPOS) opt.Remove(ip);
      // Check the content, now
      if (opt.Contains("Cq:") && !opt.Contains("-Cq:"))
         SetBit(TDataSetManager::kCheckQuota);
      if (opt.Contains("-Ar:"))
         ResetBit(TDataSetManager::kAllowRegister);
      if (opt.Contains("-Av:"))
         ResetBit(TDataSetManager::kAllowVerify);
      if (opt.Contains("-Ti:"))
         ResetBit(TDataSetManager::kTrustInfo);
      if (opt.Contains("Sb:") && !opt.Contains("-Sb:"))
         SetBit(TDataSetManager::kIsSandbox);
      if (opt.Contains("Ca:"))
         SetBit(TDataSetManager::kUseCache);
      if (opt.Contains("-Ca:"))
         SetBit(TDataSetManager::kDoNotUseCache);
   }

   // Check dependencies
   if (TestBit(TDataSetManager::kAllowVerify)) {
      // Dataset verification or requires registration permition
      SetBit(TDataSetManager::kAllowRegister);
   }
   // UseCache has priority
   if (TestBit(TDataSetManager::kUseCache) && TestBit(TDataSetManager::kDoNotUseCache))
      ResetBit(TDataSetManager::kDoNotUseCache);
}

////////////////////////////////////////////////////////////////////////////////
/// Read group config file 'cf'.
/// If cf == 0 re-read, if changed, the file pointed by fGroupConfigFile .
///
/// expects the following directives:
/// Group definition:
///   group <groupname> <user>+
/// disk quota
///   property <groupname> diskquota <quota in GB>
/// average filesize (to be used when the file size is not available)
///   averagefilesize <average size>{G,g,M,m,K,k}

Bool_t TDataSetManager::ReadGroupConfig(const char *cf)
{
   // Validate input
   FileStat_t st;
   if (!cf || (strlen(cf) <= 0) || !strcmp(cf, fGroupConfigFile.Data())) {
      // If this is the first time we cannot do anything
      if (fGroupConfigFile.IsNull()) {
         if (gDebug > 0)
            Info("ReadGroupConfig", "path to config file undefined - nothing to do");
         return kFALSE;
      }
      // Check if fGroupConfigFile has changed
      if (gSystem->GetPathInfo(fGroupConfigFile, st)) {
         Error("ReadGroupConfig", "could not stat %s", fGroupConfigFile.Data());
         return kFALSE;
      }
      if (st.fMtime <= fMTimeGroupConfig) {
         if (gDebug > 0)
            Info("ReadGroupConfig","file has not changed - do nothing");
         return kTRUE;
      }
   }

   // Either new file or the file has changed
   if (cf && (strlen(cf) > 0)) {
      // The file must exist and be readable
      if (gSystem->GetPathInfo(cf, st)) {
         Error("ReadGroupConfig", "could not stat %s", cf);
         return kFALSE;
      }
      if (gSystem->AccessPathName(cf, kReadPermission)) {
         Error("ReadGroupConfig", "cannot read %s", cf);
         return kFALSE;
      }
      // Ok
      fGroupConfigFile = cf;
      fMTimeGroupConfig = st.fMtime;
   }

   if (gDebug > 0)
      Info("ReadGroupConfig","reading group config from %s", cf);

   // Open the config file
   std::ifstream in;
   in.open(cf);
   if (!in.is_open()) {
      Error("ReadGroupConfig", "could not open config file %s", cf);
      return kFALSE;
   }

   // Container for the global common user
   TString tmpCommonUser;

   // Go through
   TString line;
   while (in.good()) {
      // Read new line
      line.ReadLine(in);
      // Explicitely skip comment lines
      if (line[0] == '#') continue;
      // Parse it
      Ssiz_t from = 0;
      TString key;
      if (!line.Tokenize(key, from, " ")) // No token
         continue;
      // Parsing depends on the key
      if (key == "property") {
         // Read group
         TString grp;
         if (!line.Tokenize(grp, from, " ")) {// No token
            if (gDebug > 0)
               Info("ReadGroupConfig","incomplete line: '%s'", line.Data());
            continue;
         }
         // Read type of property
         TString type;
         if (!line.Tokenize(type, from, " ")) // No token
            continue;
         if (type == "diskquota") {
            // Read diskquota
            TString sdq;
            if (!line.Tokenize(sdq, from, " ")) // No token
               continue;
            // Enforce GigaBytes as default
            if (sdq.IsDigit()) sdq += "G";
            Long64_t quota = ToBytes(sdq);
            if (quota > -1) {
               fGroupQuota.Add(new TObjString(grp),
                               new TParameter<Long64_t> ("group quota", quota));
            } else {
               Warning("ReadGroupConfig",
                       "problems parsing string: wrong or unsupported suffix? %s",
                        sdq.Data());
            }
         } else if (type == "commonuser") {
            // Read common user for this group
            TString comusr;
            if (!line.Tokenize(comusr, from, " ")) // No token
               continue;

         }

      } else if (key == "dataset") {
         // Read type
         TString type;
         if (!line.Tokenize(type, from, " ")) {// No token
            if (gDebug > 0)
               Info("ReadGroupConfig","incomplete line: '%s'", line.Data());
            continue;
         }
         if (type == "commonuser") {
            // Read global common user
            TString comusr;
            if (!line.Tokenize(comusr, from, " ")) // No token
               continue;
            fCommonUser = comusr;
         } else if (type == "commongroup") {
            // Read global common group
            TString comgrp;
            if (!line.Tokenize(comgrp, from, " ")) // No token
               continue;
            fCommonGroup = comgrp;
         } else if (type == "diskquota") {
            // Quota check switch
            TString on;
            if (!line.Tokenize(on, from, " ")) // No token
               continue;
            if (on == "on") {
               SetBit(TDataSetManager::kCheckQuota);
            } else if (on == "off") {
               ResetBit(TDataSetManager::kCheckQuota);
            }
         }

      } else if (key == "averagefilesize") {

         // Read average size
         TString avgsize;
         if (!line.Tokenize(avgsize, from, " ")) {// No token
            if (gDebug > 0)
               Info("ReadGroupConfig","incomplete line: '%s'", line.Data());
            continue;
         }
         Long64_t avgsz = ToBytes(avgsize);
         if (avgsz > -1) {
            fAvgFileSize = avgsz;
         } else {
            Warning("ReadGroupConfig",
                    "problems parsing string: wrong or unsupported suffix? %s",
                    avgsize.Data());
         }
      } else if (key == "include") {

         // Read file to include
         TString subfn;
         if (!line.Tokenize(subfn, from, " ")) {// No token
            if (gDebug > 0)
               Info("ReadGroupConfig","incomplete line: '%s'", line.Data());
            continue;
         }
         // The file must be readable
         if (gSystem->AccessPathName(subfn, kReadPermission)) {
            Error("ReadGroupConfig", "request to parse file '%s' which is not readable",
                                     subfn.Data());
            continue;
         }
         if (!ReadGroupConfig(subfn))
            Error("ReadGroupConfig", "problems parsing include file '%s'", subfn.Data());
      }
   }
   in.close();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Static utility function to gt the number of bytes from a string
/// representation in the form "<digit><sfx>" with <sfx> = {"", "k", "M", "G",
/// "T", "P"} (case insensitive).
/// Returns -1 if the format is wrong.

Long64_t TDataSetManager::ToBytes(const char *size)
{
   Long64_t lsize = -1;

   // Check if valid
   if (!size || strlen(size) <= 0) return lsize;

   TString s(size);
   // Determine factor
   Long64_t fact = 1;
   if (!s.IsDigit()) {
      const char *unit[5] = { "k", "M", "G", "T", "P"};
      fact = 1024;
      Int_t jj = 0;
      while (jj <= 4) {
         if (s.EndsWith(unit[jj], TString::kIgnoreCase)) {
            s.Remove(s.Length()-1);
            break;
         }
         fact *= 1024;
         jj++;
      }
   }
   // Apply factor now
   if (s.IsDigit())
      lsize = s.Atoi() * fact;

   // Done
   return lsize;
}

////////////////////////////////////////////////////////////////////////////////
/// Utility function used in various methods for user dataset upload.

TFileCollection *TDataSetManager::GetDataSet(const char *, const char *)
{
   AbstractMethod("GetDataSet");
   return (TFileCollection *)0;
}

////////////////////////////////////////////////////////////////////////////////
/// Removes the indicated dataset

Bool_t TDataSetManager::RemoveDataSet(const char *)
{
   AbstractMethod("RemoveDataSet");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the indicated dataset exits

Bool_t TDataSetManager::ExistsDataSet(const char *)
{
   AbstractMethod("ExistsDataSet");
   return kFALSE;
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

TMap *TDataSetManager::GetDataSets(const char *, UInt_t)
{
   AbstractMethod("GetDataSets");

   return (TMap *)0;
}
////////////////////////////////////////////////////////////////////////////////
/// Scans the dataset indicated by 'uri' following the 'opts' directives
///
/// The 'opts' string contains up to 4 directive fields separated by ':'
///
///  'selection' field :
///    A, allfiles:    process all files
///    D, staged:      process only staged (on Disk) files (if 'allfiles:' is not specified
///                    the default is to process only files marked as non-staged)
///  'pre-action field':
///    O, open:        open the files marked as staged when processing only files
///                    marked as non-staged
///    T, touch:       open and touch the files marked as staged when processing
///                    only files marked as non-staged
///    I, nostagedcheck: do not check the actual stage status on selected files
///
///  'process' field:
///    N, noaction:    do nothing on the selected files
///    P, fullproc:    open the selected files and extract the meta information
///    L, locateonly:  only locate the selected files
///    S, stageonly:   issue a stage request for the selected files not yet staged
///
///  'auxiliary' field
///    V, verbose:     notify the actions
///
/// Returns 0 on success, -1 if any failure occurs.

Int_t TDataSetManager::ScanDataSet(const char *uri, const char *opts)
{
   // Extract the directives
   UInt_t o = 0;
   if (opts && strlen(opts) > 0) {
      // Selection options
      if (strstr(opts, "allfiles:") || strchr(opts, 'A'))
         o |= kAllFiles;
      else if (strstr(opts, "staged:") || strchr(opts, 'D'))
         o |= kStagedFiles;
      // Pre-action options
      if (strstr(opts, "open:") || strchr(opts, 'O'))
         o |= kReopen;
      if (strstr(opts, "touch:") || strchr(opts, 'T'))
         o |= kTouch;
      if (strstr(opts, "nostagedcheck:") || strchr(opts, 'I'))
         o |= kNoStagedCheck;
      // Process options
      if (strstr(opts, "noaction:") || strchr(opts, 'N'))
         o |= kNoAction;
      if (strstr(opts, "locateonly:") || strchr(opts, 'L'))
         o |= kLocateOnly;
      if (strstr(opts, "stageonly:") || strchr(opts, 'S'))
         o |= kStageOnly;
      // Auxilliary options
      if (strstr(opts, "verbose:") || strchr(opts, 'V'))
         o |= kDebug;
   } else {
      // Default
      o = kReopen | kDebug;
   }

   // Run
   return ScanDataSet(uri, o);
}

////////////////////////////////////////////////////////////////////////////////
/// Scans the dataset indicated by <uri> and returns the number of missing files.
/// Returns -1 if any failure occurs.
/// For more details, see documentation of
/// ScanDataSet(TFileCollection *dataset, const char *option)

Int_t TDataSetManager::ScanDataSet(const char *, UInt_t)
{
   AbstractMethod("ScanDataSet");

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Gets quota information from this dataset

void TDataSetManager::GetQuota(const char *group, const char *user,
                                    const char *dsName, TFileCollection *dataset)
{
   if (gDebug > 0)
      Info("GetQuota", "processing dataset %s %s %s", group, user, dsName);

   if (dataset->GetTotalSize() > 0) {
      TParameter<Long64_t> *size =
         dynamic_cast<TParameter<Long64_t>*> (fGroupUsed.GetValue(group));
      if (!size) {
         size = new TParameter<Long64_t> ("group used", 0);
         fGroupUsed.Add(new TObjString(group), size);
      }

      size->SetVal(size->GetVal() + dataset->GetTotalSize());

      TMap *userMap = dynamic_cast<TMap*> (fUserUsed.GetValue(group));
      if (!userMap) {
         userMap = new TMap;
         fUserUsed.Add(new TObjString(group), userMap);
      }

      size = dynamic_cast<TParameter<Long64_t>*> (userMap->GetValue(user));
      if (!size) {
         size = new TParameter<Long64_t> ("user used", 0);
         userMap->Add(new TObjString(user), size);
      }

      size->SetVal(size->GetVal() + dataset->GetTotalSize());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Display quota information

void TDataSetManager::ShowQuota(const char *opt)
{
   UpdateUsedSpace();

   TMap *groupQuotaMap = GetGroupQuotaMap();
   TMap *userUsedMap = GetUserUsedMap();
   if (!groupQuotaMap || !userUsedMap)
      return;

   Bool_t noInfo = kTRUE;
   TIter iter(groupQuotaMap);
   TObjString *group = 0;
   while ((group = dynamic_cast<TObjString*> (iter.Next()))) {
      noInfo = kFALSE;
      Long64_t groupQuota = GetGroupQuota(group->String());
      Long64_t groupUsed = GetGroupUsed(group->String());

      Printf(" +++ Group %s uses %.1f GB out of %.1f GB", group->String().Data(),
                                        (Float_t) groupUsed / DSM_ONE_GB,
                                       (Float_t) groupQuota / DSM_ONE_GB);

      // display also user information
      if (opt && !TString(opt).Contains("U", TString::kIgnoreCase))
         continue;

      TMap *userMap = dynamic_cast<TMap*> (userUsedMap->GetValue(group->String()));
      if (!userMap)
         continue;

      TIter iter2(userMap);
      TObjString *user = 0;
      while ((user = dynamic_cast<TObjString*> (iter2.Next()))) {
         TParameter<Long64_t> *size2 =
            dynamic_cast<TParameter<Long64_t>*> (userMap->GetValue(user->String().Data()));
         if (!size2)
            continue;

         Printf(" +++  User %s uses %.1f GB", user->String().Data(),
                                  (Float_t) size2->GetVal() / DSM_ONE_GB);
      }

      Printf("------------------------------------------------------");
   }
   // Check if something has been printed
   if (noInfo) {
      Printf(" +++ Quota check enabled but no quota info available +++ ");
   }
}

////////////////////////////////////////////////////////////////////////////////
///
/// Prints the quota

void TDataSetManager::PrintUsedSpace()
{
   Info("PrintUsedSpace", "listing used space");

   TIter iter(&fUserUsed);
   TObjString *group = 0;
   while ((group = dynamic_cast<TObjString*> (iter.Next()))) {
      TMap *userMap = dynamic_cast<TMap*> (fUserUsed.GetValue(group->String()));

      TParameter<Long64_t> *size =
         dynamic_cast<TParameter<Long64_t>*> (fGroupUsed.GetValue(group->String()));

      if (userMap && size) {
         Printf("Group %s: %lld B = %.2f GB", group->String().Data(), size->GetVal(),
                                      (Float_t) size->GetVal() / DSM_ONE_GB);

         TIter iter2(userMap);
         TObjString *user = 0;
         while ((user = dynamic_cast<TObjString*> (iter2.Next()))) {
            TParameter<Long64_t> *size2 =
               dynamic_cast<TParameter<Long64_t>*> (userMap->GetValue(user->String().Data()));
            if (size2)
               Printf("  User %s: %lld B = %.2f GB", user->String().Data(), size2->GetVal(),
                                            (Float_t) size2->GetVal() / DSM_ONE_GB);
         }

         Printf("------------------------------------------------------");
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
///
/// Log info to the monitoring server

void TDataSetManager::MonitorUsedSpace(TVirtualMonitoringWriter *monitoring)
{
   Info("MonitorUsedSpace", "sending used space to monitoring server");

   TIter iter(&fUserUsed);
   TObjString *group = 0;
   while ((group = dynamic_cast<TObjString*> (iter.Next()))) {
      TMap *userMap = dynamic_cast<TMap*> (fUserUsed.GetValue(group->String()));
      TParameter<Long64_t> *size =
         dynamic_cast<TParameter<Long64_t>*> (fGroupUsed.GetValue(group->String()));

      if (!userMap || !size)
         continue;

      TList *list = new TList;
      list->SetOwner();
      list->Add(new TParameter<Long64_t>("_TOTAL_", size->GetVal()));
      Long64_t groupQuota = GetGroupQuota(group->String());
      if (groupQuota != -1)
         list->Add(new TParameter<Long64_t>("_QUOTA_", groupQuota));

      TIter iter2(userMap);
      TObjString *user = 0;
      while ((user = dynamic_cast<TObjString*> (iter2.Next()))) {
         TParameter<Long64_t> *size2 =
            dynamic_cast<TParameter<Long64_t>*> (userMap->GetValue(user->String().Data()));
         if (!size2)
            continue;
         list->Add(new TParameter<Long64_t>(user->String().Data(), size2->GetVal()));
      }

      if (!monitoring->SendParameters(list, group->String()))
         Warning("MonitorUsedSpace", "problems sending monitoring parameters");
      delete list;
   }
}

////////////////////////////////////////////////////////////////////////////////
///
/// Returns the used space of that group

Long64_t TDataSetManager::GetGroupUsed(const char *group)
{
   if (fgCommonDataSetTag == group)
      group = fCommonGroup;

   TParameter<Long64_t> *size =
      dynamic_cast<TParameter<Long64_t>*> (fGroupUsed.GetValue(group));
   if (!size) {
      if (gDebug > 0)
         Info("GetGroupUsed", "group %s not found", group);
      return 0;
   }

   return size->GetVal();
}

////////////////////////////////////////////////////////////////////////////////
///
/// returns the quota a group is allowed to have

Long64_t TDataSetManager::GetGroupQuota(const char *group)
{
   if (fgCommonDataSetTag == group)
      group = fCommonGroup;

   TParameter<Long64_t> *value =
      dynamic_cast<TParameter<Long64_t>*> (fGroupQuota.GetValue(group));
   if (!value) {
      if (gDebug > 0)
         Info("GetGroupQuota", "group %s not found", group);
      return 0;
   }
   return value->GetVal();
}

////////////////////////////////////////////////////////////////////////////////
/// updates the used space maps

void TDataSetManager::UpdateUsedSpace()
{
   AbstractMethod("UpdateUsedSpace");
}

////////////////////////////////////////////////////////////////////////////////
/// Register a dataset, perfoming quota checkings, if needed.
/// Returns 0 on success, -1 on failure

Int_t TDataSetManager::RegisterDataSet(const char *,
                                       TFileCollection *, const char *)
{
   AbstractMethod("RegisterDataSet");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Save into the <datasetdir>/dataset.list file the name of the last updated
/// or created or modified dataset
/// Returns 0 on success, -1 on error

Int_t TDataSetManager::NotifyUpdate(const char * /*group*/,
                                    const char * /*user*/,
                                    const char * /*dspath*/,
                                    Long_t /*mtime*/,
                                    const char * /*checksum*/)
{
   AbstractMethod("NotifyUpdate");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Clear cached information matching uri

Int_t TDataSetManager::ClearCache(const char * /*uri*/)
{
   AbstractMethod("ClearCache");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Show cached information matching uri

Int_t TDataSetManager::ShowCache(const char * /*uri*/)
{
   AbstractMethod("ShowCache");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates URI for the dataset manger in the form '[[/dsGroup/]dsUser/]dsName[#dsObjPath]',
/// The optional dsObjPath can be in the form [subdir/]objname]'.

TString TDataSetManager::CreateUri(const char *dsGroup, const char *dsUser,
                                        const char *dsName, const char *dsObjPath)
{
   TString uri;

   if (dsGroup && strlen(dsGroup) > 0) {
      if (dsUser && strlen(dsUser) > 0) {
         uri += Form("/%s/%s/", dsGroup, dsUser);
      } else {
         uri += Form("/%s/*/", dsGroup);
      }
   } else if (dsUser && strlen(dsUser) > 0) {
      uri += Form("%s/", dsUser);
   }
   if (dsName && strlen(dsName) > 0)
      uri += dsName;
   if (dsObjPath && strlen(dsObjPath) > 0)
      uri += Form("#%s", dsObjPath);

   // Done
   return uri;
}

////////////////////////////////////////////////////////////////////////////////
/// Parses a (relative) URI that describes a DataSet on the cluster.
/// The input 'uri' should be in the form '[[/group/]user/]dsname[#[subdir/]objname]',
///  where 'objname' is the name of the object (e.g. the tree name) and the 'subdir'
/// is the directory in the file wher it should be looked for.
/// After resolving against a base URI consisting of proof://masterhost/group/user/
/// - meaning masterhost, group and user of the current session -
/// the path is checked to contain exactly three elements separated by '/':
/// group/user/dsname
/// If wildcards, '*' is allowed in group and user and dsname is allowed to be empty.
/// If onlyCurrent, only group and user of current session are allowed.
/// Only non-null parameters are filled by this function.
/// Returns kTRUE in case of success.

Bool_t TDataSetManager::ParseUri(const char *uri,
                                 TString *dsGroup, TString *dsUser,
                                 TString *dsName, TString *dsTree,
                                 Bool_t onlyCurrent, Bool_t wildcards)
{
   TString uristr(uri);

   // If URI contains fields in the form "Field=Value;" it is a virtual URI and
   // should be treated differently
   if ((uristr.Index('=') >= 0) && (uristr.Index(';') >= 0)) {

      // URI is composed of two parts: a name (dsName), and the tree after the
      // pound sign

      Warning("ParseUri",
        "Dataset URI looks like a virtual URI, treating it as such. "
        "No group and user will be parsed!");

      TPMERegexp reVirtualUri("^([^#]+)(#(.*))?$");
      Int_t nm = reVirtualUri.Match(uristr);

      if (nm >= 2) {
         if (dsGroup) *dsGroup = "";
         if (dsUser) *dsUser = "";
         if (dsName) *dsName = reVirtualUri[1];
         if (dsTree) {
            if (nm == 4) *dsTree = reVirtualUri[3];
            else *dsTree = "";
         }
      }
      else return kFALSE;  // should never happen!

      return kTRUE;
   }

   // Append trailing slash if missing when wildcards are enabled
   Int_t pc = 0;
   if (wildcards && uristr.Length() > 0) {
      pc = uristr.CountChar('/');
      Bool_t endsl = uristr.EndsWith("/") ? kTRUE : kFALSE;
      Bool_t beginsl = uristr.BeginsWith("/") ? kTRUE : kFALSE;
      if (beginsl) {
         if (pc == 1) uristr += "/*/";
         if (pc == 2 && endsl) uristr += "*/";
         if (pc == 2 && !endsl) uristr += "/";
      }
   }

   // Resolve given URI agains the base
   TUri resolved = TUri::Transform(uristr, fBase);
   if (resolved.HasQuery())
      Info ("ParseUri", "URI query part <%s> ignored", resolved.GetQuery().Data());

   TString path(resolved.GetPath());
   // Must be in the form /group/user/dsname
   if ((pc = path.CountChar('/')) != 3) {
      if (!TestBit(TDataSetManager::kIsSandbox)) {
         Error ("ParseUri", "illegal dataset path: '%s'", uri);
         return kFALSE;
      } else if (pc >= 0 && pc < 3) {
         // Add missing slashes
         TString sls("/");
         if (pc == 2) {
            sls = "/";
         } else if (pc == 1) {
            sls.Form("/%s/", fGroup.Data());
         } else if (pc == 0) {
            sls.Form("/%s/%s/", fGroup.Data(), fUser.Data());
         }
         path.Insert(0, sls);
      }
   }
   if (gDebug > 1)
      Info("ParseUri", "path: '%s'", path.Data());

   // Get individual values from tokens
   Int_t from = 1;
   TString group, user, name;
   if (path.Tokenize(group, from, "/")) {
      if (path.Tokenize(user, from, "/")) {
         if (!path.Tokenize(name, from, "/"))
            if (gDebug > 0) Info("ParseUri", "'name' missing");
      } else {
         if (gDebug > 0) Info("ParseUri", "'user' missing");
      }
   } else {
      if (gDebug > 1) Info("ParseUri", "'group' missing");
   }

   // The fragment may contain the subdir and the object name in the form '[subdir/]objname'
   TString tree = resolved.GetFragment();
   if (tree.EndsWith("/"))
      tree.Remove(tree.Length()-1);

   if (gDebug > 1)
      Info("ParseUri", "group: '%s', user: '%s', dsname:'%s', seg: '%s'",
                              group.Data(), user.Data(), name.Data(), tree.Data());

   // Check for unwanted use of wildcards
   if ((user == "*" || group == "*") && !wildcards) {
      Error ("ParseUri", "no wildcards allowed for user/group in this context (uri: '%s')", uri);
      return kFALSE;
   }

   // dsname may only be empty if wildcards expected
   if (name.IsNull() && !wildcards) {
      Error ("ParseUri", "DataSet name is empty");
      return kFALSE;
   }

   // Construct regexp whitelist for checking illegal characters in user/group
   TPRegexp wcExp (wildcards ? "^(?:[A-Za-z0-9-*_.]*|[*])$" : "^[A-Za-z0-9-_.]*$");

   // Check for illegal characters in all components
   if (!wcExp.Match(group)) {
      Error("ParseUri", "illegal characters in group (uri: '%s', group: '%s')", uri, group.Data());
      return kFALSE;
   }

   if (!wcExp.Match(user)) {
      Error("ParseUri", "illegal characters in user (uri: '%s', user: '%s')", uri, user.Data());
      return kFALSE;
   }

   // Construct regexp whitelist for checking illegal characters in name
   if (!wcExp.Match(name)) {
      Error("ParseUri", "illegal characters in name (uri: '%s', name: '%s')", uri, name.Data());
      return kFALSE;
   }

   if (tree.Contains(TRegexp("[^A-Za-z0-9-/_]"))) {
      Error("ParseUri", "Illegal characters in subdir/object name (uri: '%s', obj: '%s')", uri, tree.Data());
      return kFALSE;
   }

   // Check user & group
   if (onlyCurrent && (group.CompareTo(fGroup) || user.CompareTo(fUser))) {
      Error("ParseUri", "only datasets from your group/user allowed");
      return kFALSE;
   }

   // fill parameters passed by reference, if defined
   if (dsGroup)
      *dsGroup = group;
   if (dsUser)
      *dsUser = user;
   if (dsName)
      *dsName = name;
   if (dsTree)
      *dsTree = tree;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Partition dataset 'ds' accordingly to the servers.
/// The returned TMap contains:
///                <server> --> <subdataset> (TFileCollection)
/// where <subdataset> is the subset of 'ds' on <server>
/// The partitioning is done using all the URLs in the TFileInfo's, so the
/// resulting datasets are not mutually exclusive.
/// The string 'exclude' contains a comma-separated list of servers to exclude
/// from the map.

TMap *TDataSetManager::GetSubDataSets(const char *ds, const char *exclude)
{
   TMap *map = (TMap *)0;

   if (!ds || strlen(ds) <= 0) {
      Info("GetDataSets", "dataset name undefined!");
      return map;
   }

   // Get the dataset
   TFileCollection *fc = GetDataSet(ds);
   if (!fc) {
      Info("GetDataSets", "could not retrieve the dataset '%s'", ds);
      return map;
   }

   // Get the subset
   if (!(map = fc->GetFilesPerServer(exclude))) {
      if (gDebug > 0)
         Info("GetDataSets", "could not get map for '%s'", ds);
   }

   // Cleanup
   delete fc;

   // Done
   return map;
}

////////////////////////////////////////////////////////////////////////////////
/// Formatted printout of the content of TFileCollection 'fc'.
/// Options in the form
///           popt = u * 10 + f
///     f    0 => header only, 1 => header + files
///   when printing files
///     u    0 => print file name only, 1 => print full URL

void TDataSetManager::PrintDataSet(TFileCollection *fc, Int_t popt)
{
   if (!fc) return;

   Int_t f = popt%10;
   Int_t u = popt - 10 * f;

   Printf("+++");
   if (fc->GetTitle() && (strlen(fc->GetTitle()) > 0)) {
      Printf("+++ Dumping: %s: ", fc->GetTitle());
   } else {
      Printf("+++ Dumping: %s: ", fc->GetName());
   }
   Printf("%s", fc->ExportInfo("+++ Summary:", 1)->GetName());
   if (f == 1) {
      Printf("+++ Files:");
      Int_t nf = 0;
      TIter nxfi(fc->GetList());
      TFileInfo *fi = 0;
      while ((fi = (TFileInfo *)nxfi())) {
         if (u == 1)
            Printf("+++ %5d. %s", ++nf, fi->GetCurrentUrl()->GetUrl());
         else
            Printf("+++ %5d. %s", ++nf, fi->GetCurrentUrl()->GetFile());
      }
   }
   Printf("+++");
}

////////////////////////////////////////////////////////////////////////////////
/// Prints formatted information about the dataset 'uri'.
/// The type and format of output is driven by 'opt':
///
///   1. opt = "server:srv1[,srv2[,srv3[,...]]]"
///            Print info about the subsets of 'uri' on servers srv1, srv2, ...
///   2. opt = "servers[:exclude:srv1[,srv2[,srv3[,...]]]]"
///            Print info about the subsets of 'uri' on all servers, except
///            the ones in the exclude list srv1, srv2, ...
///   3. opt = <any>
///            Print info about all datasets matching 'uri'
///
///   If 'opt' contains 'full:' the list of files in the datasets are also printed.
///   In case 3. this is enabled only if 'uri' matches a single dataset.
///
///   In case 3, if 'opt' contains
///      'full:'      the list of files in the datasets are also printed.
///      'forcescan:' the dataset are open to get the information; otherwise the
///                   pre-processed information is used.
///      'noheader:'  the labelling header is not printed; usefull when to chain
///                   several printouts
///      'noupdate:'  do not update the cache (which may be slow on very remote
///                   servers)
///      'refresh:'   refresh the information (requires appropriate credentials;
///                   typically it can be done only for owned datasets)

void TDataSetManager::ShowDataSets(const char *uri, const char *opt)
{
   TFileCollection *fc = 0;
   TString o(opt);
   Int_t popt = 0;
   if (o.Contains("full:")) {
      o.ReplaceAll("full:","");
      popt = 1;
   }
   if (o.BeginsWith("server:")) {
      o.ReplaceAll("server:", "");
      TString srv;
      Int_t from = 0;
      while ((o.Tokenize(srv, from, ","))) {
         fc = GetDataSet(uri, srv.Data());
         PrintDataSet(fc, popt);
         delete fc;
      }
   } else if (o.BeginsWith("servers")) {
      o.ReplaceAll("servers", "");
      if (o.BeginsWith(":exclude:"))
         o.ReplaceAll(":exclude:", "");
      else
         o = "";
      TMap *dsmap = GetSubDataSets(uri, o.Data());
      if (dsmap) {
         TIter nxk(dsmap);
         TObject *k = 0;
         while ((k = nxk()) && (fc = (TFileCollection *) dsmap->GetValue(k))) {
            PrintDataSet(fc, popt);
         }
         delete dsmap;
      }
   } else {
      TString u(uri), grp, usr, dsn;
      // Support for "*" or "/*"
      if (u == "" || u == "*" || u == "/*" || u == "/*/" || u == "/*/*") u = "/*/*/";
      if (!ParseUri(u.Data(), &grp, &usr, &dsn, 0, kFALSE, kTRUE))
         Warning("ShowDataSets", "problems parsing URI '%s'", uri);
      // Scan the existing datasets and print the content
      UInt_t xopt = (UInt_t)(TDataSetManager::kPrint);
      if (o.Contains("forcescan:")) xopt |= (UInt_t)(TDataSetManager::kForceScan);
      if (o.Contains("noheader:")) xopt |= (UInt_t)(TDataSetManager::kNoHeaderPrint);
      if (o.Contains("noupdate:")) xopt |= (UInt_t)(TDataSetManager::kNoCacheUpdate);
      if (o.Contains("refresh:")) xopt |= (UInt_t)(TDataSetManager::kRefreshLs);
      if (!u.IsNull() && !u.Contains("*") && !grp.IsNull() && !usr.IsNull() && !dsn.IsNull()) {
         if (ExistsDataSet(uri)) {
            // Single dataset
            if (popt == 0) {
               // Quick listing
               GetDataSets(u.Data(), xopt);
            } else if ((fc = GetDataSet(uri))) {
               // Full print option
               PrintDataSet(fc, 10 + popt);
               delete fc;
            }
            return;
         }
         // Try all the directories
         TRegexp reg(grp, kTRUE), reu(usr, kTRUE);
         if (u.Index(reg) == kNPOS) grp = "*";
         if (u.Index(reu) == kNPOS) usr = "*";
         // Rebuild the uri
         u.Form("/%s/%s/%s", grp.Data(), usr.Data(), dsn.Data());
      }
      GetDataSets(u.Data(), xopt);
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Go through the files in the specified dataset, selecting files according to
/// 'fopt' and doing on these files the actions described by 'sopt'.
/// If required, the information in 'dataset' is updated.
///
/// The int fopt controls which files have to be processed (or added to the list
/// if ropt is 1 - see below); 'fopt' is defined in term of csopt and fsopt:
///                    fopt = sign(fsopt) * csopt * 100 + fsopt
/// where 'fsopt' controls the actual selection
///    -1              all files in the dataset
///     0              process only files marked as 'non-staged'
///   >=1              as 0 but files that are marked 'staged' are open
///   >=2              as 1 but files that are marked 'staged' are touched
///    10              process only files marked as 'staged'; files marked as 'non-staged'
///                    are ignored
/// and 'csopt' controls if an actual check on the staged status (via TFileStager) is done
///     0              check that the file is staged using TFileStager
///     1              do not hard check the staged status
/// (example: use fopt = -101 to check the staged status of all the files, or fopt = 110
///  to re-check the stage status of all the files marked as staged)
///
/// If 'dbg' is true, some information about the ongoing operations is reguraly
/// printed; this can be useful when processing very large datasets, an operation
/// which can take a very long time.
///
/// The int 'sopt' controls what is done on the selected files (this is effective only
/// if ropt is 0 or 2 - see below):
///    -1              no action (fopt = 2 and sopt = -1 touches all staged files)
///     0              do the full process: open the files and fill the meta-information
///                    in the TFileInfo object, including the end-point URL
///     1              only locate the files, by updating the end-point URL (uses TFileStager::Locate
///                    which is must faster of an TFile::Open)
///     2              issue a stage request on the files
///
/// The int 'ropt' controls which actions are performed:
///     0              do the full process: get list of files to process and process them
///     1              get the list of files to be scanned and return it in flist
///     2              process the files in flist (according to sopt)
/// When defined flist is under the responsability the caller.
///
/// If avgsz > 0 it is used for the final update of the dataset global counters.
///
/// If 'mss' is defined use it to initialize the stager (instead of the Url in the
/// TFileInfo objects)
///
/// If maxfiles > 0, select for processing a maximum of 'filesmax' files (but if fopt is 1 or 2
/// all files marked as 'staged' are still open or touched)
///
/// Return code
///     1 dataset was not changed
///     2 dataset was changed
///
/// The number of touched, opened and disappeared files are returned in the respective
/// variables, if these are defined.

Int_t TDataSetManager::ScanDataSet(TFileCollection *dataset,
                                   Int_t fopt, Int_t sopt, Int_t ropt, Bool_t dbg,
                                   Int_t *touched, Int_t *opened, Int_t *disappeared,
                                   TList *flist, Long64_t avgsz, const char *mss,
                                   Int_t maxfiles, const char *stageopts)
{
   // Max number of files
   if (maxfiles > -1 && dbg)
      ::Info("TDataSetManager::ScanDataSet", "processing a maximum of %d files", maxfiles);

   // File selection, Reopen and Touch options
   Bool_t checkstg = (fopt >= 100 || fopt < -1) ? kFALSE : kTRUE;

   // File processing options
   Bool_t noaction   = (sopt == -1) ? kTRUE : kFALSE;
   //Bool_t fullproc   = (sopt == 0)  ? kTRUE : kFALSE;
   Bool_t locateonly = (sopt == 1)  ? kTRUE : kFALSE;
   Bool_t stageonly  = (sopt == 2)  ? kTRUE : kFALSE;

   // Run options
   Bool_t doall       = (ropt == 0) ? kTRUE : kFALSE;
   Bool_t getlistonly = (ropt == 1) ? kTRUE : kFALSE;
   Bool_t scanlist    = (ropt == 2) ? kTRUE : kFALSE;

   if (scanlist && !flist) {
      ::Error("TDataSetManager::ScanDataSet", "input list is mandatory for option 'scan file list'");
      return -1;
   }

   Int_t ftouched = 0;
   Int_t fopened = 0;
   Int_t fdisappeared = 0;

   Bool_t bchanged_ds = kFALSE;

   TList *newStagedFiles = 0;
   TFileInfo *fileInfo = 0;
   TFileStager *stager = 0;
   Bool_t createStager = kFALSE;

   if (doall || getlistonly) {

      // Point to the list
      newStagedFiles = (!doall && getlistonly && flist) ? flist : new TList;
      if (newStagedFiles != flist) newStagedFiles->SetOwner(kFALSE);

      stager = (mss && strlen(mss) > 0) ? TFileStager::Open(mss) : 0;
      createStager = (stager) ? kFALSE : kTRUE;

      Bool_t bchanged_fi = kFALSE;
      Bool_t btouched = kFALSE;
      Bool_t bdisappeared = kFALSE;

      // Check which files have been staged, this can be replaced by a bulk command,
      // once it exists in the xrdclient
      TIter iter(dataset->GetList());
      while ((fileInfo = (TFileInfo *) iter())) {

         // For real time monitoring
         gSystem->DispatchOneEvent(kTRUE);

         bchanged_fi = kFALSE;
         btouched = kFALSE;
         bdisappeared = kFALSE;
         Bool_t newlystaged = CheckStagedStatus(fileInfo, fopt, maxfiles, newStagedFiles->GetEntries(),
                                                stager, createStager, dbg, bchanged_fi, btouched,
                                                bdisappeared);

         if (bchanged_fi) bchanged_ds = kTRUE;
         if (btouched) ftouched++;
         if (bdisappeared) fdisappeared++;

         // Notify
         if (dbg && (ftouched+fdisappeared) % 100 == 0)
            ::Info("TDataSetManager::ScanDataSet", "opening %d: file: %s",
                   ftouched + fdisappeared, fileInfo->GetCurrentUrl()->GetUrl());

         // Register the newly staged file
         if (!noaction && newlystaged) newStagedFiles->Add(fileInfo);
      }
      SafeDelete(stager);

      // If required to only get the list we are done
      if (getlistonly) {
         if (dbg && newStagedFiles->GetEntries() > 0)
            ::Info("TDataSetManager::ScanDataSet", " %d files appear to be newly staged",
                                                   newStagedFiles->GetEntries());
         if (!flist) SafeDelete(newStagedFiles);
         return ((bchanged_ds) ? 2 : 1);
      }
   }

   if (!noaction && (doall || scanlist)) {

      // Point to the list
      newStagedFiles = (!doall && scanlist && flist) ? flist : newStagedFiles;
      if (newStagedFiles != flist) newStagedFiles->SetOwner(kFALSE);

      // loop over now staged files
      if (dbg && newStagedFiles->GetEntries() > 0)
         ::Info("TDataSetManager::ScanDataSet", "opening %d files that appear to be newly staged",
                                                newStagedFiles->GetEntries());

      // If staging files, prepare the stager
      if (locateonly || stageonly) {
         stager = (mss && strlen(mss) > 0) ? TFileStager::Open(mss) : 0;
         createStager = (stager) ? kFALSE : kTRUE;
      }

      // Notify each 'fqnot' files (min 1, max 100)
      Int_t fqnot = (newStagedFiles->GetSize() > 10) ? newStagedFiles->GetSize() / 10 : 1;
      if (fqnot > 100) fqnot = 100;
      Int_t count = 0;
      Bool_t bchanged_fi = kFALSE;
      Bool_t bopened = kFALSE;
      TIter iter(newStagedFiles);
      while ((fileInfo = (TFileInfo *) iter())) {

         if (dbg && (count%fqnot == 0))
            ::Info("TDataSetManager::ScanDataSet", "processing %d.'new' file: %s",
                                                   count, fileInfo->GetCurrentUrl()->GetUrl());
         count++;

         // For real time monitoring
         gSystem->DispatchOneEvent(kTRUE);
         bchanged_fi = kFALSE;
         bopened = kFALSE;

         ProcessFile(fileInfo, sopt, checkstg, doall, stager, createStager,
                      stageopts, dbg, bchanged_fi, bopened);

         bchanged_ds |= bchanged_fi;
         if (bopened) fopened++;
      }
      if (newStagedFiles != flist) SafeDelete(newStagedFiles);

      dataset->RemoveDuplicates();
      dataset->Update(avgsz);
   }

   Int_t result = (bchanged_ds) ? 2 : 1;
   if (result > 0 && dbg)
      ::Info("TDataSetManager::ScanDataSet", "%d files 'new'; %d files touched;"
                                             " %d files disappeared", fopened, ftouched, fdisappeared);

   // Fill outputs, if required
   if (touched) *touched = ftouched;
   if (opened) *opened = fopened;
   if (disappeared) *disappeared = fdisappeared;

   // For real time monitoring
   gSystem->DispatchOneEvent(kTRUE);

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Check stage status of the file described by "fileInfo".
/// fopt is same as "fopt" in TDataSetManager::ScanDataSet, which is repeated below:
/// The int fopt controls which files have to be processed (or added to the list
/// if ropt is 1 - see below); 'fopt' is defined in term of csopt and fsopt:
///                    fopt = sign(fsopt) * csopt * 100 + fsopt
/// where 'fsopt' controls the actual selection
///    -1              all files in the dataset
///     0              process only files marked as 'non-staged'
///   >=1              as 0 but files that are marked 'staged' are open
///   >=2              as 1 but files that are marked 'staged' are touched
///    10              process only files marked as 'staged'; files marked as 'non-staged'
///                    are ignored
/// and 'csopt' controls if an actual check on the staged status (via TFileStager) is done
///     0              check that the file is staged using TFileStager
///     1              do not hard check the staged status
/// (example: use fopt = -101 to check the staged status of all the files, or fopt = 110
///  to re-check the stage status of all the files marked as staged)
///
/// If 'dbg' is true, some information about the ongoing operations is reguraly
/// printed; this can be useful when processing very large datasets, an operation
/// which can take a very long time.
///
/// If maxfiles > 0, select for processing a maximum of 'filesmax' files (but if fopt is 1 or 2
/// all files marked as 'staged' are still open or touched)
///
/// Return code
///     kTRUE the file appears newly staged
///     kFALSE otherwise
///
/// changed is true if the fileinfo is modified
/// touched is true if the file is open and read
/// disappeared is true if the file is marked staged but actually not staged

Bool_t TDataSetManager::CheckStagedStatus(TFileInfo *fileInfo, Int_t fopt, Int_t maxfiles,
                                          Int_t newstagedfiles, TFileStager* stager,
                                          Bool_t createStager, Bool_t dbg, Bool_t& changed,
                                          Bool_t& touched, Bool_t& disappeared)
{
   // File selection, Reopen and Touch options
   Bool_t allf     = (fopt == -1)               ? kTRUE : kFALSE;
   Bool_t checkstg = (fopt >= 100 || fopt < -1) ? kFALSE : kTRUE;
   if (fopt >= 0) fopt %= 100;
   Bool_t nonstgf  = (fopt >= 0 && fopt < 10)   ? kTRUE : kFALSE;
   Bool_t reopen   = (fopt >= 1 && fopt < 10)   ? kTRUE : kFALSE;
   Bool_t touch    = (fopt >= 2 && fopt < 10)   ? kTRUE : kFALSE;
   Bool_t stgf     = (fopt == 10)               ? kTRUE : kFALSE;

   changed = kFALSE;
   touched = kFALSE;
   disappeared = kFALSE;

   // Check which files have been staged, this can be replaced by a bulk command,
   // once it exists in the xrdclient

   if (!allf) {

      fileInfo->ResetUrl();
      if (!fileInfo->GetCurrentUrl()) {
         ::Error("TDataSetManager::CheckStagedStatus", "GetCurrentUrl() returned 0 for %s",
                                                fileInfo->GetFirstUrl()->GetUrl());
         return kFALSE;
      }

      if (nonstgf && fileInfo->TestBit(TFileInfo::kStaged)) {

         // Skip files flagged as corrupted
         if (fileInfo->TestBit(TFileInfo::kCorrupted)) return kFALSE;

         // Skip if we are not asked to re-open the staged files
         if (!reopen) return kFALSE;

         // Set the URL removing the anchor (e.g. #AliESDs.root) because IsStaged()
         // and TFile::Open() with filetype=raw do not accept anchors
         TUrl *curl = fileInfo->GetCurrentUrl();
         const char *furl = curl->GetUrl();
         TString urlmod;
         if (TDataSetManager::CheckDataSetSrvMaps(curl, urlmod) && !(urlmod.IsNull()))
            furl = urlmod.Data();
         TUrl url(furl);
         url.SetAnchor("");

         // Check if file is still available, if touch is set actually read from the file
         TString uopt(url.GetOptions());
         uopt += "filetype=raw&mxredir=2";
         url.SetOptions(uopt.Data());
         TFile *file = TFile::Open(url.GetUrl());
         if (file) {
            if (touch) {
               // Actually access the file
               char tmpChar = 0;
               if (file->ReadBuffer(&tmpChar, 1))
                  ::Warning("TDataSetManager::CheckStagedStatus", "problems reading 1 byte from open file");
               // Count
               touched = kTRUE;
            }
            file->Close();
            delete file;
         } else {
            // File could not be opened, reset staged bit
            if (dbg) ::Info("TDataSetManager::CheckStagedStatus", "file %s disappeared", url.GetUrl());
            fileInfo->ResetBit(TFileInfo::kStaged);
            disappeared = kTRUE;
            changed = kTRUE;

            // Remove invalid URL, if other one left...
            if (fileInfo->GetNUrls() > 1)
               fileInfo->RemoveUrl(curl->GetUrl());
         }
         // Go to next
         return kFALSE;
      } else if (stgf && !(fileInfo->TestBit(TFileInfo::kStaged))) {
         // All staged files are processed: skip non staged
         return kFALSE;
      }
   }

   // Only open maximum number of 'new' files
   if (maxfiles > 0 && newstagedfiles >= maxfiles)
      return kFALSE;

   // Hard check of the staged status, if required
   if (checkstg) {
      // Set the URL removing the anchor (e.g. #AliESDs.root) because IsStaged()
      // and TFile::Open() with filetype=raw do not accept anchors
      TUrl *curl = fileInfo->GetCurrentUrl();
      const char *furl = curl->GetUrl();
      TString urlmod;
      Bool_t mapped = kFALSE;
      if (TDataSetManager::CheckDataSetSrvMaps(curl, urlmod) && !(urlmod.IsNull())) {
         furl = urlmod.Data();
         mapped = kTRUE;
      }
      TUrl url(furl);
      url.SetAnchor("");

      // Get the stager (either the global one or from the URL)
      stager = createStager ? TFileStager::Open(url.GetUrl()) : stager;

      Bool_t result = kFALSE;
      if (stager) {
         result = stager->IsStaged(url.GetUrl());
         if (gDebug > 0)
            ::Info("TDataSetManager::CheckStagedStatus", "IsStaged: %s: %d", url.GetUrl(), result);
         if (createStager)
            SafeDelete(stager);
      } else {
         ::Warning("TDataSetManager::CheckStagedStatus",
                  "could not get stager instance for '%s'", url.GetUrl());
      }

      // Go to next in case of failure
      if (!result) {
         if (fileInfo->TestBit(TFileInfo::kStaged)) {
            // Reset the bit
            fileInfo->ResetBit(TFileInfo::kStaged);
            changed = kTRUE;
         }
         return kFALSE;
      } else {
         if (!(fileInfo->TestBit(TFileInfo::kStaged))) {
            // Set the bit
            fileInfo->SetBit(TFileInfo::kStaged);
            changed = kTRUE;
         }
      }

      // If the url was re-mapped add the new url in front of the list
      if (mapped) {
         url.SetOptions(curl->GetOptions());
         url.SetAnchor(curl->GetAnchor());
         fileInfo->AddUrl(url.GetUrl(), kTRUE);
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Locate, stage, or fully validate file "fileInfo".

void TDataSetManager::ProcessFile(TFileInfo *fileInfo, Int_t sopt, Bool_t checkstg, Bool_t doall,
                                  TFileStager* stager,  Bool_t createStager, const char *stageopts,
                                  Bool_t dbg, Bool_t& changed, Bool_t& opened)
{
   // File processing options
   //Bool_t noaction   = (sopt == -1) ? kTRUE : kFALSE;
   Bool_t fullproc   = (sopt == 0)  ? kTRUE : kFALSE;
   Bool_t locateonly = (sopt == 1)  ? kTRUE : kFALSE;
   Bool_t stageonly  = (sopt == 2)  ? kTRUE : kFALSE;

   changed = kFALSE;
   opened = kFALSE;
   Int_t rc = -1;

   // Set the URL removing the anchor (e.g. #AliESDs.root) because IsStaged()
   // and TFile::Open() with filetype=raw do not accept anchors
   TUrl *curl = fileInfo->GetCurrentUrl();
   const char *furl = curl->GetUrl();
   TString urlmod;
   //Bool_t mapped = kFALSE;
   if (TDataSetManager::CheckDataSetSrvMaps(curl, urlmod) && !(urlmod.IsNull())) {
      furl = urlmod.Data();
      //mapped = kTRUE;
   }
   TUrl url(furl);
   url.SetOptions("");
   url.SetAnchor("");

   if (createStager){
      if (!stager || (stager && !stager->Matches(url.GetUrl()))) {
         SafeDelete(stager);
         if (!(stager = TFileStager::Open(url.GetUrl())) || !(stager->IsValid())) {
            ::Error("TDataSetManager::ProcessFile",
                     "could not get valid stager instance for '%s'", url.GetUrl());
            return;
         }
      }
   }
   // Locate the file, if just requested so
   if (locateonly) {
      TString eurl;
      if (stager && stager->Locate(url.GetUrl(), eurl) == 0) {
         TString opts(curl->GetOptions());
         TString anch(curl->GetAnchor());
         // Get the effective end-point Url
         curl->SetUrl(eurl);
         // Restore original options and anchor, if any
         curl->SetOptions(opts);
         curl->SetAnchor(anch);
         // Flag and count
         changed = kTRUE;
         opened = kTRUE;
      } else {
         // Failure
         ::Error("TDataSetManager::ProcessFile", "could not locate %s", url.GetUrl());
      }

   } else if (stageonly) {
      TString eurl;
      if (stager && !(stager->IsStaged(url.GetUrl()))) {
         if (!(stager->Stage(url.GetUrl(), stageopts))) {
            // Failure
            ::Error("TDataSetManager::ProcessFile",
                     "problems issuing stage request for %s", url.GetUrl());
         }
      }
   } else if (fullproc) {
      TString eurl;
      // Full file validation
      rc = -2;
      Bool_t doscan = kTRUE;
      if (checkstg) {
         doscan = kFALSE;
         if ((doall && fileInfo->TestBit(TFileInfo::kStaged)) ||
             (stager && stager->IsStaged(url.GetUrl()))) doscan = kTRUE;
      }
      if (doscan) {
         if ((rc = TDataSetManager::ScanFile(fileInfo, dbg)) < -1) return;
         changed = kTRUE;
      } else if (stager) {
         ::Warning("TDataSetManager::ProcessFile",
                   "required file '%s' does not look as being online (staged)", url.GetUrl());
      }
      if (rc < 0) return;
      // Count
      opened = kTRUE;
   }
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Open the file described by 'fileinfo' to extract the relevant meta-information.
/// Return 0 if OK, -2 if the file cannot be open, -1 if it is corrupted

Int_t TDataSetManager::ScanFile(TFileInfo *fileinfo, Bool_t dbg)
{
   Int_t rc = -2;
   // We need an input
   if (!fileinfo) {
      ::Error("TDataSetManager::ScanFile", "undefined input (!)");
      return rc;
   }

   TUrl *url = fileinfo->GetCurrentUrl();

   TFile *file = 0;
   Bool_t anchor = kFALSE;

   // Get timeout settings (default none)
   Int_t timeout = gEnv->GetValue("DataSet.ScanFile.OpenTimeout", -1);
   TString fileopt;
   if (timeout > 0) fileopt.Form("TIMEOUT=%d", timeout);

   // To determine the size we have to open the file without the anchor
   // (otherwise we get the size of the contained file - in case of a zip archive)
   // We open in raw mode which makes sure that the opening succeeds, even if
   // the file is corrupted
   const char *furl = url->GetUrl();
   TString urlmod;
   if (TDataSetManager::CheckDataSetSrvMaps(url, urlmod) && !(urlmod.IsNull()))
      furl = urlmod.Data();
   if (strlen(url->GetAnchor()) > 0) {
      anchor = kTRUE;
      // We need a raw open firts to get the real size of the file
      TUrl urlNoAnchor(furl);
      urlNoAnchor.SetAnchor("");
      TString unaopts = urlNoAnchor.GetOptions();
      if (!unaopts.IsNull()) {
         unaopts += "&filetype=raw";
      } else {
         unaopts = "filetype=raw";
      }
      urlNoAnchor.SetOptions(unaopts);
      // Wait max 5 secs per file
      if (!(file = TFile::Open(urlNoAnchor.GetUrl(), fileopt))) return rc;

      // Save some relevant info
      if (file->GetSize() > 0) fileinfo->SetSize(file->GetSize());
      fileinfo->SetBit(TFileInfo::kStaged);

      fileinfo->SetUUID(file->GetUUID().AsString());

      // Add url of the disk server in front of the list
      if (file->GetEndpointUrl()) {
         // add endpoint url if it is not a local file
         TUrl eurl(*(file->GetEndpointUrl()));

         if (strcmp(eurl.GetProtocol(), "file") ||
            !strcmp(eurl.GetProtocol(), url->GetProtocol())) {

            eurl.SetOptions(url->GetOptions());
            eurl.SetAnchor(url->GetAnchor());

            // Fix the hostname
            if (!strcmp(eurl.GetHost(), "localhost") || !strcmp(eurl.GetHost(), "127.0.0.1") ||
               !strcmp(eurl.GetHost(), "localhost.localdomain")) {
               eurl.SetHost(TUrl(gSystem->HostName()).GetHostFQDN());
            }
            // Add only if different
            if (strcmp(eurl.GetUrl(), url->GetUrl()))
               fileinfo->AddUrl(eurl.GetUrl(), kTRUE);

            if (gDebug > 0) ::Info("TDataSetManager::ScanFile", "added URL %s", eurl.GetUrl());
         }
      } else {
         ::Warning("TDataSetManager::ScanFile", "end-point URL undefined for file %s", file->GetName());
      }

      file->Close();
      delete file;
   }

   // OK, set the relevant flags
   rc = -1;

   // Disable warnings when reading a tree without loading the corresponding library
   Int_t oldLevel = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kError+1;

   // Wait max 5 secs per file
   if (!(file = TFile::Open(url->GetUrl(), fileopt))) {
      // If the file could be opened before, but fails now it is corrupt...
      if (dbg) ::Info("TDataSetManager::ScanFile", "marking %s as corrupt", url->GetUrl());
      fileinfo->SetBit(TFileInfo::kCorrupted);
      // Set back old warning level
      gErrorIgnoreLevel = oldLevel;
      return rc;
   } else if (!anchor) {
      // Do the relevant settings
      if (file->GetSize() > 0) fileinfo->SetSize(file->GetSize());
      fileinfo->SetBit(TFileInfo::kStaged);

      // Add url of the disk server in front of the list if it is not a local file
      TUrl eurl(*(file->GetEndpointUrl()));

      if (strcmp(eurl.GetProtocol(), "file") ||
         !strcmp(eurl.GetProtocol(), url->GetProtocol())) {

         eurl.SetOptions(url->GetOptions());
         eurl.SetAnchor(url->GetAnchor());

         // Fix the hostname
         if (!strcmp(eurl.GetHost(), "localhost") || !strcmp(eurl.GetHost(), "127.0.0.1") ||
             !strcmp(eurl.GetHost(), "localhost.localdomain")) {
             eurl.SetHost(TUrl(gSystem->HostName()).GetHostFQDN());
         }
         // Add only if different
         if (strcmp(eurl.GetUrl(), url->GetUrl()))
            fileinfo->AddUrl(eurl.GetUrl(), kTRUE);

         if (gDebug > 0) ::Info("TDataSetManager::ScanFile", "added URL %s", eurl.GetUrl());
      }
      fileinfo->SetUUID(file->GetUUID().AsString());
   }
   rc = 0;

   // Loop over all entries and create/update corresponding metadata.
   // TODO If we cannot read some of the trees, is the file corrupted as well?
   if ((rc = TDataSetManager::FillMetaData(fileinfo, file, "/")) != 0) {
      ::Error("TDataSetManager::ScanFile",
              "problems processing the directory tree in looking for metainfo");
      fileinfo->SetBit(TFileInfo::kCorrupted);
      rc = -1;
   }
   // Set back old warning level
   gErrorIgnoreLevel = oldLevel;

   file->Close();
   delete file;

   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Navigate the directory 'd' (and its subdirectories) looking for TTree objects.
/// Fill in the relevant metadata information in 'fi'. The name of the TFileInfoMeta
/// metadata entry will be "/dir1/dir2/.../tree_name".
/// Return 0 on success, -1 if any problem happens (object found in keys cannot be read,
/// for example)

Int_t TDataSetManager::FillMetaData(TFileInfo *fi, TDirectory *d, const char *rdir)
{
   // Check inputs
   if (!fi || !d || !rdir) {
      ::Error("TDataSetManager::FillMetaData",
              "some inputs are invalid (fi:%p,d:%p,r:%s)", fi, d, rdir);
      return -1;
   }

   if (d->GetListOfKeys()) {
      TIter nxk(d->GetListOfKeys());
      TKey *k = 0;
      while ((k = dynamic_cast<TKey *> (nxk()))) {

         if (TClass::GetClass(k->GetClassName())->InheritsFrom(TDirectory::Class())) {
            // Get the directory
            TDirectory *sd = (TDirectory *) d->Get(k->GetName());
            if (!sd) {
               ::Error("TDataSetManager::FillMetaData", "cannot get sub-directory '%s'", k->GetName());
               return -1;
            }
            if (TDataSetManager::FillMetaData(fi, sd, TString::Format("%s%s/", rdir, k->GetName())) != 0) {
               ::Error("TDataSetManager::FillMetaData", "problems processing sub-directory '%s'", k->GetName());
               return -1;
            }

         } else {
            // We process only trees
            if (!TClass::GetClass(k->GetClassName())->InheritsFrom(TTree::Class())) continue;

            TString ks;
            ks.Form("%s%s", rdir, k->GetName());

            TFileInfoMeta *md = fi->GetMetaData(ks);
            if (!md) {
               // Create it
               md = new TFileInfoMeta(ks, k->GetClassName());
               fi->AddMetaData(md);
               if (gDebug > 0)
                  ::Info("TDataSetManager::FillMetaData", "created meta data for tree %s", ks.Data());
            }
            // Fill values
            TTree *t = dynamic_cast<TTree *> (d->Get(k->GetName()));
            if (t) {
               if (t->GetEntries() >= 0) {
                  md->SetEntries(t->GetEntries());
                  if (t->GetTotBytes() >= 0)
                     md->SetTotBytes(t->GetTotBytes());
                  if (t->GetZipBytes() >= 0)
                     md->SetZipBytes(t->GetZipBytes());
               }
            } else {
               ::Error("TDataSetManager::FillMetaData", "could not get tree '%s'", k->GetName());
               return -1;
            }
         }
      }
   }
   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a server mapping list from the content of 'srvmaps'
/// Return the list (owned by the caller) or 0 if no valid info could be found)

TList *TDataSetManager::ParseDataSetSrvMaps(const TString &srvmaps)
{
   TList *srvmapslist = 0;
   if (srvmaps.IsNull()) {
      ::Warning("TDataSetManager::ParseDataSetSrvMaps",
                "called with an empty string! - nothing to do");
      return srvmapslist;
   }
   TString srvmap, sf, st;
   Int_t from = 0, from1 = 0;
   while (srvmaps.Tokenize(srvmap, from, " ")) {
      sf = ""; st = "";
      if (srvmap.Contains("|")) {
         from1 = 0;
         if (srvmap.Tokenize(sf, from1, "|"))
            if (srvmap.Tokenize(st, from1, "|")) { }
      } else {
         st = srvmap;
      }
      if (st.IsNull()) {
         ::Warning("TDataSetManager::ParseDataSetSrvMaps",
                   "parsing DataSet.SrvMaps: target must be defined"
                   " (token: %s) - ignoring", srvmap.Data());
         continue;
      } else if (!(st.EndsWith("/"))) {
         st += "/";
      }
      // TUrl if wildcards or TObjString
      TString sp;
      TUrl *u = 0;
      if (!(sf.IsNull()) && sf.Contains("*")) {
         u = new TUrl(sf);
         if (!(sf.BeginsWith(u->GetProtocol()))) u->SetProtocol("root");
         sp.Form(":%d", u->GetPort());
         if (!(sf.Contains(sp))) u->SetPort(1094);
         if (!TString(u->GetHost()).Contains("*")) SafeDelete(u);
      }
      if (!srvmapslist) srvmapslist = new TList;
      if (u) {
         srvmapslist->Add(new TPair(u, new TObjString(st)));
      } else {
         srvmapslist->Add(new TPair(new TObjString(sf), new TObjString(st)));
      }
   }
   // Done
   if (srvmapslist) srvmapslist->SetOwner(kTRUE);
   return srvmapslist;
}

////////////////////////////////////////////////////////////////////////////////
/// Static getter for server mapping list

TList *TDataSetManager::GetDataSetSrvMaps()
{
   return fgDataSetSrvMaps;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the dataset server mappings apply to the url defined by 'furl'.
/// Use srvmaplist if defined, else use the default list.
/// If yes, resolve the mapping into file1 and return kTRUE.
/// Otherwise return kFALSE.

Bool_t TDataSetManager::CheckDataSetSrvMaps(TUrl *furl, TString &file1, TList *srvmaplist)
{
   Bool_t replaced = kFALSE;
   if (!furl) return replaced;

   const char *file = furl->GetUrl();
   TList *mlist = (srvmaplist) ? srvmaplist : fgDataSetSrvMaps;
   if (mlist && mlist->GetSize() > 0) {
      TIter nxm(mlist);
      TPair *pr = 0;
      while ((pr = (TPair *) nxm())) {
         Bool_t replace = kFALSE;
         // If TUrl apply reg exp on host
         TUrl *u = dynamic_cast<TUrl *>(pr->Key());
         if (u) {
            if (!strcmp(u->GetProtocol(), furl->GetProtocol())) {
               Ssiz_t len;
               if (!strcmp(u->GetProtocol(), "file")) {
                  TRegexp re(u->GetFileAndOptions(), kTRUE);
                  if (re.Index(furl->GetFileAndOptions(), &len) == 0) replace = kTRUE;
               } else {
                  if (u->GetPort() == furl->GetPort()) {
                     TRegexp re(u->GetHost(), kTRUE);
                     if (re.Index(furl->GetHost(), &len) == 0) replace = kTRUE;
                  }
               }
            }
         } else {
            TObjString *os = dynamic_cast<TObjString *>(pr->Key());
            if (os) {
               if (os->GetString().IsNull() ||
                   !strncmp(file, os->GetName(), os->GetString().Length())) replace = kTRUE;
            }
         }
         if (replace) {
            TObjString *ost = dynamic_cast<TObjString *>(pr->Value());
            if (ost) {
               file1.Form("%s%s", ost->GetName(), furl->GetFileAndOptions());
               replaced = kTRUE;
               break;
            }
         }
      }
   }
   // Done
   return replaced;
}

////////////////////////////////////////////////////////////////////////////////
/// Update scan counters

void TDataSetManager::SetScanCounters(Int_t t, Int_t o, Int_t d)
{
   fNTouchedFiles = (t > -1) ? t : fNTouchedFiles;
   fNOpenedFiles = (o > -1) ? o : fNOpenedFiles;
   fNDisappearedFiles = (d > -1) ? d : fNDisappearedFiles;
}
