// @(#)root/base:$Id$
// Author: Gerhard Erich Bruckner, Jan Fiete Grosse-Oetringhaus  04/06/07

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TFileCollection
\ingroup Base

Class that contains a list of TFileInfo's and accumulated meta
data information about its entries. This class is used to describe
file sets as stored by Grid file catalogs, by PROOF or any other
collection of TFile names.
*/

#include "TFileCollection.h"
#include "THashList.h"
#include "TFileInfo.h"
#include "TIterator.h"
#include "TMap.h"
#include "TObjString.h"
#include "TUri.h"
#include "TUrl.h"
#include "TUUID.h"
#include "TMD5.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TError.h"

#include <iostream>
#include <fstream>


ClassImp(TFileCollection);

////////////////////////////////////////////////////////////////////////////////
/// TFileCollection constructor. Specify a name and title describing
/// the list. If textfile is specified the file is opened and a
/// TFileCollection is created containing the files in the textfile.

TFileCollection::TFileCollection(const char *name, const char *title,
                                 const char *textfile, Int_t nfiles, Int_t firstfile)
   : TNamed(name, title), fList(0), fMetaDataList(0), fDefaultTree(),
     fTotalSize(0), fNFiles(0), fNStagedFiles(0), fNCorruptFiles(0)
{
   fList = new THashList();
   fList->SetOwner();

   fMetaDataList = new TList;
   fMetaDataList->SetOwner();

   AddFromFile(textfile, nfiles, firstfile);
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup.

TFileCollection::~TFileCollection()
{
   delete fList;
   delete fMetaDataList;
}

////////////////////////////////////////////////////////////////////////////////
/// Add TFileInfo to the collection.

Int_t TFileCollection::Add(TFileInfo *info)
{
   if (fList && info) {
      if (!fList->FindObject(info->GetName())) {
         fList->Add(info);
         if (info->GetIndex() < 0) info->SetIndex(fList->GetSize());
         Update();
         return 1;
      } else {
         Warning("Add", "file: '%s' already in the list - ignoring",
                        info->GetCurrentUrl()->GetUrl());
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add content of the TFileCollection to this collection.

Int_t TFileCollection::Add(TFileCollection *coll)
{
   if (fList && coll && coll->GetList()) {
      TIter nxfi(coll->GetList());
      TFileInfo *fi = 0;
      while ((fi = (TFileInfo *) nxfi())) {
         TFileInfo *info = new TFileInfo(*fi);
         fList->Add(info);
         if (fi->GetIndex() < 0) info->SetIndex(fList->GetSize());
      }
      Update();
      return 1;
   } else {
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add file names contained in the specified text file.
/// The file should contain one url per line; empty lines or lines starting with '#'
/// (commented lines) are ignored.
/// If nfiles > 0 only nfiles files are added, starting from file 'firstfile' (>= 1).
/// The method returns the number of added files.

Int_t TFileCollection::AddFromFile(const char *textfile, Int_t nfiles, Int_t firstfile)
{
   if (!fList)
     return 0;

   Int_t nf = 0;
   TString fn(textfile);
   if (!fn.IsNull() && !gSystem->ExpandPathName(fn)) {
      std::ifstream f;
      f.open(fn);
      if (f.is_open()) {
         Bool_t all = (nfiles <= 0) ? kTRUE : kFALSE;
         Int_t ff = (!all && (firstfile < 1)) ? 1 : firstfile;
         Int_t nn = 0;
         while (f.good() && (all || nf < nfiles)) {
            TString line;
            line.ReadToDelim(f);
            // Skip commented or empty lines
            if (!line.IsWhitespace() && !line.BeginsWith("#")) {
               nn++;
               if (all || nn >= ff) {
                  TFileInfo *info = new TFileInfo(line);
                  fList->Add(info);
                  if (info->GetIndex() < 0) info->SetIndex(fList->GetSize());
                  nf++;
               }
            }
         }
         f.close();
         Update();
      } else
         Error("AddFromFile", "unable to open file %s (%s)", textfile, fn.Data());
   }
   return nf;
}

////////////////////////////////////////////////////////////////////////////////
/// Add all files matching the specified pattern to the collection.
/// 'dir' can include wildcards after the last slash, which causes all
/// matching files in that directory to be added.
/// If dir is the full path of a file, only one element is added.
/// Return value is the number of added files.

Int_t TFileCollection::Add(const char *dir)
{
   Int_t nf = 0;

   if (!fList)
      return nf;

   if (!dir || !*dir) {
      Error("Add", "input dir undefined");
      return nf;
   }

   FileStat_t st;
   FileStat_t tmp;
   TString baseDir = gSystem->GetDirName(dir);
   // if the 'dir' or its base dir exist
   if (gSystem->GetPathInfo(dir, st) == 0 ||
       gSystem->GetPathInfo(baseDir, tmp) == 0) {
      // If 'dir' points to a single file, add to the list and exit
      if (R_ISREG(st.fMode)) {
         // regular, single file
         TFileInfo *info = new TFileInfo(dir);
         info->SetBit(TFileInfo::kStaged);
         Add(info);
         nf++;
         Update();
         return nf;
      } else {
         void *dataSetDir = gSystem->OpenDirectory(gSystem->GetDirName(dir).Data());
         if (!dataSetDir) {
            // directory cannot be opened
            Error("Add", "directory %s cannot be opened",
                  gSystem->GetDirName(dir).Data());
         } else {
            const char *ent;
            TString filesExp(TString("^") + gSystem->BaseName(dir) + "$");
            filesExp.ReplaceAll("*",".*");
            TRegexp rg(filesExp);
            while ((ent = gSystem->GetDirEntry(dataSetDir))) {
               TString entryString(ent);
               if (entryString.Index(rg) != kNPOS) {
                  // matching dir entry
                  TString fn = gSystem->GetDirName(dir);
                  fn += "/";
                  fn += ent;
                  gSystem->GetPathInfo(fn, st);
                  if (R_ISREG(st.fMode)) {
                     // regular file
                     TFileInfo *info = new TFileInfo(fn);
                     info->SetBit(TFileInfo::kStaged);
                     Add(info);
                     nf++;
                  }
               }
            }
            // close the directory
            gSystem->FreeDirectory(dataSetDir);
            Update();
         }
      }
   }
   return nf;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove duplicates based on the UUID, typically after a verification.
/// Return the number of entries removed.

Int_t TFileCollection::RemoveDuplicates()
{
   THashList *hl = new THashList;
   hl->SetOwner();

   Int_t n0 = fList->GetSize();
   TIter nxfi(fList);
   TFileInfo *fi = 0;
   while ((fi = (TFileInfo *)nxfi())) {
      if (!(hl->FindObject(fi->GetUUID()->AsString()))) {
         // We hash on the UUID
         fList->Remove(fi);
         fi->SetName(fi->GetUUID()->AsString());
         hl->Add(fi);
      }
   }
   delete fList;
   fList = hl;
   // How many removed?
   Int_t nr = n0 - fList->GetSize();
   if (nr > 0)
      Info("RemoveDuplicates", "%d duplicates found and removed", nr);
   // Done
   return nr;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a subset of the files that have the kStaged & !kCorrupted bit set.

TFileCollection *TFileCollection::GetStagedSubset()
{
   if (!fList)
     return 0;

   TFileCollection *subset = new TFileCollection(GetName(), GetTitle());

   TIter iter(fList);
   TFileInfo *fileInfo = 0;
   while ((fileInfo = dynamic_cast<TFileInfo*>(iter.Next()))) {
      if (fileInfo->TestBit(TFileInfo::kStaged) && !fileInfo->TestBit(TFileInfo::kCorrupted))
         subset->Add(fileInfo);
   }

   subset->Update();

   return subset;
}

////////////////////////////////////////////////////////////////////////////////
/// Merge all TFileCollection objects in li into this TFileCollection object.
/// Updates counters at the end.
/// Returns the number of merged collections or -1 in case of error.

Long64_t TFileCollection::Merge(TCollection *li)
{

   if (!li) return 0;
   if (li->IsEmpty()) return 0;

   Long64_t nentries=0;
   TIter next(li);
   while (TObject *o = next()) {
      TFileCollection* coll = dynamic_cast<TFileCollection*> (o);
      if (!coll) {
         Error("Add", "attempt to add object of class: %s to a %s",
                      o->ClassName(),this->ClassName());
         return -1;
      }
      Add(coll);
      nentries++;
   }
   Update();

   return nentries;
}

////////////////////////////////////////////////////////////////////////////////
/// Update accumulated information about the elements of the collection
/// (e.g. fTotalSize). If 'avgsize' > 0, use an average file size of 'avgsize'
/// bytes when the size info is not available.
/// Also updates the meta data information by summarizing
/// the meta data of the contained objects.
/// Return -1 in case of any failure, 0 if the total size is exact, 1 if
/// incomplete, 2 if complete but (at least partially) estimated.

Int_t TFileCollection::Update(Long64_t avgsize)
{
   if (!fList)
     return -1;

   Int_t rc = 0;

   fTotalSize = 0;
   fNStagedFiles = 0;
   fNCorruptFiles = 0;

   // Clear internal meta information which is going to be rebuilt in this
   // function
   TIter nxm(fMetaDataList);
   TFileInfoMeta *m = 0;
   while ((m = (TFileInfoMeta *)nxm())) {
      if (!(m->TestBit(TFileInfoMeta::kExternal))) {
         fMetaDataList->Remove(m);
         delete m;
      }
   }

   fNFiles = fList->GetEntries();

   TIter iter(fList);
   TFileInfo *fileInfo = 0;
   while ((fileInfo = dynamic_cast<TFileInfo*> (iter.Next()))) {

      if (fileInfo->GetSize() > 0) {
         fTotalSize += fileInfo->GetSize();
      } else {
         rc = 1;
         if (avgsize > 0) {
            rc = 2;
            fTotalSize += avgsize;
         }
      }

      if (fileInfo->TestBit(TFileInfo::kStaged) && !fileInfo->TestBit(TFileInfo::kCorrupted)) {
         fNStagedFiles++;

         if (fileInfo->GetMetaDataList()) {
            TIter metaDataIter(fileInfo->GetMetaDataList());
            // other than TFileInfoMeta is also allowed in list
            TObject *obj = 0;
            while ((obj = metaDataIter.Next())) {
               TFileInfoMeta *metaData = dynamic_cast<TFileInfoMeta*>(obj);
               if (!metaData)
                  continue;
               if (!metaData->IsTree())
                  continue;

               // find corresponding entry in TFileCollection's meta data
               TFileInfoMeta *metaDataSum = dynamic_cast<TFileInfoMeta*>(fMetaDataList->FindObject(metaData->GetName()));
               Bool_t newObj = kFALSE;
               if (!metaDataSum) {
                  // create explicitly, there are some values that do not make sense for the sum
                  metaDataSum = new TFileInfoMeta(metaData->GetName(), metaData->GetTitle());
                  fMetaDataList->Add(metaDataSum);
                  newObj = kTRUE;
               }

               // sum the values
               if (newObj)
                  metaDataSum->SetEntries(metaData->GetEntries());
               else
                  metaDataSum->SetEntries(metaDataSum->GetEntries() + metaData->GetEntries());
            }
         }
      }
      if (fileInfo->TestBit(TFileInfo::kCorrupted))
         fNCorruptFiles++;
   }

   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Prints the contents of the TFileCollection.
/// If option contains:
///
///  - 'M'             print global meta information
///  - 'F'             print all the files in the collection in compact form
///                    (current url, default tree name|class|entries, md5)
///  - 'L'             together with 'F', print all the files in the collection
///                    in long form (uuid, md5, all URLs, all meta objects; on
///                    many lines)
///  - "filter:[SsCc]" invokes PrintDetailed() which prints out dataset
///                    content in a formatted fashion by filtering on files
///                    which are (S)taged or not (s), (C)orrupted or not (c)

void TFileCollection::Print(Option_t *option) const
{
   TString opt(option);
   TPMERegexp re("(^|;| )filter:([SsCc]+)( |;|$)", 4);
   if (re.Match(option) == 4) {
     TString showOnly = re[2];
     PrintDetailed(showOnly);
     return;
   }

   Printf("TFileCollection %s - %s contains: %lld files with a size of"
          " %lld bytes, %.1f %% staged - default tree name: '%s'",
          GetName(), GetTitle(), fNFiles, fTotalSize, GetStagedPercentage(),
          GetDefaultTreeName());

   if (opt.Contains("M", TString::kIgnoreCase)) {
      Printf("The files contain the following trees:");

      TIter metaDataIter(fMetaDataList);
      TFileInfoMeta* metaData = 0;
      while ((metaData = dynamic_cast<TFileInfoMeta*>(metaDataIter.Next()))) {
         if (!metaData->IsTree())
            continue;

         Printf("Tree %s: %lld events", metaData->GetName(), metaData->GetEntries());
      }
   }

   if (fList && opt.Contains("F", TString::kIgnoreCase)) {
      Printf("The collection contains the following files:");
      if (!opt.Contains("L") && !fDefaultTree.IsNull())
         opt += TString::Format(" T:%s", fDefaultTree.Data());
      fList->Print(opt);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print detailed.

void TFileCollection::PrintDetailed(TString &showOnly) const
{
   Bool_t bS, bs, bC, bc;
   bS = bs = bC = bc = kFALSE;

   if (showOnly.Index('S') >= 0) bS = kTRUE;
   if (showOnly.Index('s') >= 0) bs = kTRUE;
   if (showOnly.Index('C') >= 0) bC = kTRUE;
   if (showOnly.Index('c') >= 0) bc = kTRUE;

   // If Ss (or Cc) omitted, show both Ss (or Cc)
   if (!bc && !bC) bc = bC = kTRUE;
   if (!bs && !bS) bs = bS = kTRUE;

   TIter it(fList);
   TFileInfo *info;
   UInt_t countAll = 0;
   UInt_t countMatch = 0;

   Printf("\033[1m   #. SC | Entries | Size       | URL\033[m");

   TString um;
   Double_t sz;

   while ((info = dynamic_cast<TFileInfo *>(it.Next()))) {

      Bool_t s = info->TestBit(TFileInfo::kStaged);
      Bool_t c = info->TestBit(TFileInfo::kCorrupted);

      TUrl *url;

      countAll++;

      if ( ((s && bS) || (!s && bs)) && ((c && bC) || (!c && bc)) ) {

         TFileInfoMeta *meta = info->GetMetaData();  // gets the first one
         Int_t entries = -1;

         if (meta) entries = meta->GetEntries();

         FormatSize(info->GetSize(), um, sz);

         // First line: current URL with all information
         info->ResetUrl();
         TUrl *curUrl = info->GetCurrentUrl();
         const char *curUrlStr = curUrl ? curUrl->GetUrl() : "n.a.";
         Printf("\033[1m%4u.\033[m %c%c | %-7s | %6.1lf %s | %s",
           ++countMatch,
           (s ? 'S' : 's'), (c ? 'C' : 'c'),
           ((entries > 0) ? Form("% 7d", entries) : "n.a."),
           sz, um.Data(), curUrlStr);
         info->NextUrl();

         // Every other URL shown below current one
         while ((url = info->NextUrl())) {
            Printf("         |         |            | %s", url->GetUrl());
         }
         info->ResetUrl();

      } // end match filters

   } // end loop over entries

   if (countAll) {

     Printf(">> There are \033[1m%u\033[m file(s) in dataset: "
        "\033[1m%u (%5.1f%%)\033[m matched your criteria (%s)",
        countAll, countMatch,
        100.*(Float_t)countMatch/(Float_t)countAll, showOnly.Data());

    FormatSize(fTotalSize, um, sz);
    Printf(">> Total size    : \033[1m%.1f %s\033[m", sz, um.Data());
    Printf(">> Staged (S)    : \033[1m%5.1f %%\033[m", GetStagedPercentage());
    Printf(">> Corrupted (C) : \033[1m%5.1f %%\033[m",
      GetCorruptedPercentage());

  }
  else {
    Printf(">> No files in dataset");
  }

  const char *treeName = GetDefaultTreeName();
  Printf(">> Default tree  : \033[1m%s\033[m",
    (treeName ? treeName : "(no default tree)"));

}

////////////////////////////////////////////////////////////////////////////////
/// Format size.

void TFileCollection::FormatSize(Long64_t bytes, TString &um,
   Double_t &size) const
{
   static const char *ums[] = { "byt", "KiB", "MiB", "GiB", "TiB" };
   Int_t maxDiv = sizeof(ums)/sizeof(const char *);
   Int_t nDiv = 0;
   Double_t b = bytes;

   while ((b >= 1024.) && (nDiv+1 < maxDiv)) {
      b /= 1024.;
      nDiv++;
   }

   um = ums[nDiv];
   size = b;
}

////////////////////////////////////////////////////////////////////////////////
/// Calls TUrl::SetAnchor() for all URLs contained in all TFileInfos.

void TFileCollection::SetAnchor(const char *anchor)
{
   if (!fList)
     return;

   TIter iter(fList);
   TFileInfo *fileInfo = 0;
   while ((fileInfo = dynamic_cast<TFileInfo*>(iter.Next()))) {
      fileInfo->ResetUrl();
      TUrl *url = 0;
      while ((url = fileInfo->NextUrl()))
         url->SetAnchor(anchor);
      fileInfo->ResetUrl();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the bit for all TFileInfos

void TFileCollection::SetBitAll(UInt_t f)
{
   if (!fList)
     return;

   TIter iter(fList);
   TFileInfo *fileInfo = 0;
   while ((fileInfo = dynamic_cast<TFileInfo*>(iter.Next())))
      fileInfo->SetBit(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the bit for all TFileInfos

void TFileCollection::ResetBitAll(UInt_t f)
{
   if (!fList)
     return;

   TIter iter(fList);
   TFileInfo *fileInfo = 0;
   while ((fileInfo = dynamic_cast<TFileInfo*>(iter.Next())))
      fileInfo->ResetBit(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the tree set with SetDefaultTreeName if set
/// Returns the name of the first tree in the meta data list.
/// Returns 0 in case no trees are found in the meta data list.

const char *TFileCollection::GetDefaultTreeName() const
{
   if (fDefaultTree.Length() > 0)
     return fDefaultTree;

   TIter metaDataIter(fMetaDataList);
   TFileInfoMeta *metaData = 0;
   while ((metaData = dynamic_cast<TFileInfoMeta*>(metaDataIter.Next()))) {
      if (!metaData->IsTree())
         continue;
      return metaData->GetName();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the number of entries for the specified tree (retrieved from meta data).
/// If tree is not specified, use the default tree name.
/// Returns -1 in case the specified tree is not found.

Long64_t TFileCollection::GetTotalEntries(const char *tree) const
{
   if (!tree || !*tree) {
      tree = GetDefaultTreeName();
      if (!tree)
         return -1;
   }

   TFileInfoMeta *metaData = dynamic_cast<TFileInfoMeta*>(fMetaDataList->FindObject(tree));
   if (!metaData)
      return -1;

   return metaData->GetEntries();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the meta data object with the specified meta name.
/// Returns 0 in case specified meta data is not found.

TFileInfoMeta *TFileCollection::GetMetaData(const char *meta) const
{
   if (!meta || !*meta)
      return 0;

   return dynamic_cast<TFileInfoMeta*>(fMetaDataList->FindObject(meta));
}

////////////////////////////////////////////////////////////////////////////////
/// Moves the indicated meta data in the first position, so that
/// it becomes effectively the default.

void TFileCollection::SetDefaultMetaData(const char *meta)
{
   TFileInfoMeta *fim = GetMetaData(meta);
   if (fim) {
      fMetaDataList->Remove(fim);
      fMetaDataList->AddFirst(fim);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Removes the indicated meta data object in all TFileInfos and this object
/// If no name is given all metadata is removed

void TFileCollection::RemoveMetaData(const char *meta)
{
   if (fList) {
      TIter iter(fList);
      TFileInfo *fileInfo = 0;
      while ((fileInfo = dynamic_cast<TFileInfo*>(iter.Next())))
         fileInfo->RemoveMetaData(meta);
   }

   if (meta) {
      TObject* obj = fMetaDataList->FindObject("meta");
      if (obj) {
         fMetaDataList->Remove(obj);
         delete obj;
      }
   } else
      fMetaDataList->Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Sort the collection.

void TFileCollection::Sort(Bool_t useindex)
{
   if (!fList)
     return;

   // Make sure the relevant bit has the wanted value
   if (useindex) {
      SetBitAll(TFileInfo::kSortWithIndex);
   } else {
      ResetBitAll(TFileInfo::kSortWithIndex);
   }

   fList->Sort();
}

////////////////////////////////////////////////////////////////////////////////
/// Export the relevant info as a string; use 'name' as collection name,
/// if defined, else use GetName().
/// The output object must be destroyed by the caller

TObjString *TFileCollection::ExportInfo(const char *name, Int_t popt)
{
   TString treeInfo;
   if (GetDefaultTreeName()) {
      TFileInfoMeta* meta = GetMetaData(GetDefaultTreeName());
      if (popt == 1) {
         treeInfo = GetDefaultTreeName();
         if (meta)
            treeInfo += TString::Format(", %lld entries", meta->GetEntries());
         TFileInfoMeta *frac = GetMetaData("/FractionOfTotal");
         if (frac)
            treeInfo += TString::Format(", %3.1f %% of total", frac->GetEntries() / 10.);
      } else {
         treeInfo.Form(" %s ", GetDefaultTreeName());
         if (treeInfo.Length() > 14) treeInfo.Replace(13, 1, '>');
         treeInfo.Resize(14);
         if (meta) {
            if (meta->GetEntries() > 99999999) {
               treeInfo += TString::Format("| %8lld ", meta->GetEntries());
            } else {
               treeInfo += TString::Format("| %8.4g ", (Double_t) meta->GetEntries());
            }
         }
      }
   } else {
      treeInfo = "        N/A";
   }
   if (popt == 0) treeInfo.Resize(25);

   // Renormalize the size to kB, MB or GB
   const char *unit[4] = {"kB", "MB", "GB", "TB"};
   Int_t k = 0;
   Long64_t refsz = 1024;
   Long64_t xsz = (Long64_t) (GetTotalSize() / refsz);
   while (xsz > 1024 && k < 3) {
      k++;
      refsz *= 1024;
      xsz = (Long64_t) (GetTotalSize() / refsz);
   }

   // The name
   TString dsname(name);
   if (dsname.IsNull()) dsname = GetName();

   // Create the output string
   TObjString *outs = 0;
   if (popt == 1) {
      outs = new TObjString(Form("%s %lld files, %lld %s, staged %d %%, tree: %s", dsname.Data(),
                                 GetNFiles(), xsz, unit[k],
                                 (Int_t)GetStagedPercentage(), treeInfo.Data()));
   } else {
      outs = new TObjString(Form("%s| %7lld |%s| %5lld %s |  %3d %%", dsname.Data(),
                                 GetNFiles(), treeInfo.Data(), xsz, unit[k],
                                 (Int_t)GetStagedPercentage()));
   }
   // Done
   return outs;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the subset of files served by 'server'. The syntax for 'server' is
/// the standard URI one, i.e. `[<scheme>://]<host>[:port]`

TFileCollection *TFileCollection::GetFilesOnServer(const char *server)
{
   TFileCollection *fc = (TFileCollection *)0;

   // Server specification is mandatory
   if (!server || strlen(server) <= 0) {
      Info("GetFilesOnServer", "server undefined - do nothing");
      return fc;
   }

   // Nothing to do for empty lists
   if (!fList || fList->GetSize() <= 0) {
      Info("GetFilesOnServer", "the list is empty - do nothing");
      return fc;
   }

   // Define the server reference string
   TUri uri(server);
   TString srv, scheme("root"), port;
   if (uri.GetScheme() != "") scheme = uri.GetScheme();
   if (uri.GetPort() != "") port.Form(":%s", uri.GetPort().Data());
   srv.Form("%s://%s%s", scheme.Data(), TUrl(server).GetHostFQDN(), port.Data());
   if (gDebug > 0)
      Info("GetFilesOnServer", "searching for files on server: '%s' (input: '%s')",
                               srv.Data(), server);

   // Prepare the output
   fc = new TFileCollection(GetName());
   TString title;
   if (GetTitle() && strlen(GetTitle()) > 0) {
      title.Form("%s (subset on server %s)", GetTitle(), srv.Data());
   } else {
      title.Form("subset of '%s' on server %s", GetName(), srv.Data());
   }
   fc->SetTitle(title.Data());
   // The default tree name
   fc->SetDefaultTreeName(GetDefaultTreeName());

   // We look for URL starting with srv
   srv.Insert(0, "^");

   // Go through the list
   TIter nxf(fList);
   TFileInfo *fi = 0;
   while ((fi = (TFileInfo *)nxf())) {
      TUrl *xu = 0;
      if ((xu = fi->FindByUrl(srv.Data()))) {
         // Create a new TFileInfo object
         TFileInfo *nfi = new TFileInfo(xu->GetUrl(), fi->GetSize(),
                                        fi->GetUUID() ? fi->GetUUID()->AsString() : 0,
                                        fi->GetMD5() ? fi->GetMD5()->AsString() : 0);
         if (fi->GetMetaDataList()) {
            TIter nxm(fi->GetMetaDataList());
            TFileInfoMeta *md = 0;
            while ((md = (TFileInfoMeta *) nxm())) {
               nfi->AddMetaData(new TFileInfoMeta(*md));
            }
         }
         if (fi->TestBit(TFileInfo::kStaged)) nfi->SetBit(TFileInfo::kStaged);
         if (fi->TestBit(TFileInfo::kCorrupted)) nfi->SetBit(TFileInfo::kCorrupted);
         if (gDebug > 1)
            Info("GetFilesOnServer", "adding: %s", xu->GetUrl());
         fc->Add(nfi);
      }
   }

   // If nothing found, delete the object
   if (fc->GetList()->GetSize() <= 0) {
      delete fc;
      fc = 0;
      Info("GetFilesOnServer", "dataset '%s' has no files on server: '%s' (searched for: '%s')",
                               GetName(), server, srv.Data());
   }

   // Fill up sums on the sub file collection
   if (fc) {
      fc->Update();
      // Fraction of total in permille
      Long64_t xf = (fc->GetTotalSize() * 1000) / GetTotalSize();
      TFileInfoMeta *m = new TFileInfoMeta("FractionOfTotal", "External Info", xf);
      m->SetBit(TFileInfoMeta::kExternal);
      fc->AddMetaData(m);
   }

   // Done
   return fc;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a map of TFileCollections with the files on each data server,
/// excluding servers in the comma-separated list 'exclude'.
/// If curronly is kTRUE, only the URL flagged as current in the TFileInfo
/// are considered.

TMap *TFileCollection::GetFilesPerServer(const char *exclude, Bool_t curronly)
{
   TMap *dsmap = 0;

   // Nothing to do for empty lists
   if (!fList || fList->GetSize() <= 0) {
      Info("GetFilesPerServer", "the list is empty - do nothing");
      return dsmap;
   }

   // List of servers to be ignored
   THashList *excl = 0;
   if (exclude && strlen(exclude) > 0) {
      excl = new THashList;
      excl->SetOwner();
      TUri uri;
      TString srvs(exclude), s, srv, scheme, port;
      Int_t from = 0;
      while (srvs.Tokenize(s, from, ",")) {
         uri.SetUri(s.Data());
         scheme = "root";
         port = "";
         if (uri.GetScheme() != "") scheme = uri.GetScheme();
         if (uri.GetPort() != "") port.Form(":%s", uri.GetPort().Data());
         srv.Form("%s://%s%s", scheme.Data(), TUrl(s.Data()).GetHostFQDN(), port.Data());
         // Add
         excl->Add(new TObjString(srv.Data()));
      }
   }

   // Prepare the output
   dsmap = new TMap();

   // Go through the list
   TIter nxf(fList);
   TFileInfo *fi = 0;
   TUri uri;
   TString key;
   TFileCollection *fc = 0;
   while ((fi = (TFileInfo *)nxf())) {
      // Save current URL
      TUrl *curl = fi->GetCurrentUrl();
      // Loop over URLs
      if (!curronly) fi->ResetUrl();
      TUrl *xurl = 0;
      while ((xurl = (curronly) ? curl : fi->NextUrl())) {
         // Find the key for this server
         key.Form("%s://%s", xurl->GetProtocol(), xurl->GetHostFQDN());
         // Check if this has to be ignored
         if (excl && excl->FindObject(key.Data())) {
            if (curronly) break;
            continue;
         } else if (excl && xurl->GetPort() > 0) {
            // Complete the key, if needed, and recheck
            key += TString::Format(":%d", xurl->GetPort());
            if (excl->FindObject(key.Data())) {
               if (curronly) break;
               continue;
            }
         }
         // Get the map entry for this key
         TPair *ent = 0;
         if (!(ent = (TPair *) dsmap->FindObject(key.Data()))) {
            // Create the TFileCollection
            fc = new TFileCollection(GetName());
            TString title;
            if (GetTitle() && strlen(GetTitle()) > 0) {
               title.Form("%s (subset on server %s)", GetTitle(), key.Data());
            } else {
               title.Form("subset of '%s' on server %s", GetName(), key.Data());
            }
            fc->SetTitle(title.Data());
            // The default tree name
            fc->SetDefaultTreeName(GetDefaultTreeName());
            // Add it to the map
            dsmap->Add(new TObjString(key.Data()), fc);
            // Notify
            if (gDebug > 0)
               Info("GetFilesPerServer", "found server: '%s' (fc: %p)", key.Data(), fc);
         } else {
            // Attach to the TFileCollection
            fc = (TFileCollection *) ent->Value();
         }
         // Create a new TFileInfo object
         TFileInfo *nfi = new TFileInfo(xurl->GetUrl(kTRUE), fi->GetSize(),
                                        fi->GetUUID() ? fi->GetUUID()->AsString() : 0,
                                        fi->GetMD5() ? fi->GetMD5()->AsString() : 0);
         if (fi->GetMetaDataList()) {
            TIter nxm(fi->GetMetaDataList());
            TFileInfoMeta *md = 0;
            while ((md = (TFileInfoMeta *) nxm())) {
               nfi->AddMetaData(new TFileInfoMeta(*md));
            }
         }
         if (fi->TestBit(TFileInfo::kStaged)) nfi->SetBit(TFileInfo::kStaged);
         if (fi->TestBit(TFileInfo::kCorrupted)) nfi->SetBit(TFileInfo::kCorrupted);
         fc->Add(nfi);
         // In current_only mode we are done
         if (curronly) break;
      }
      // Restore current URL
      fi->SetCurrentUrl(curl);
   }

   // Fill up sums on the sub file collections
   TIter nxk(dsmap);
   TObject *k = 0;
   while ((k = nxk()) && (fc = (TFileCollection *) dsmap->GetValue(k))) {
      fc->Update();
      // Fraction of total in permille
      Long64_t xf = (fc->GetTotalSize() * 1000) / GetTotalSize();
      TFileInfoMeta *m = new TFileInfoMeta("FractionOfTotal", "External Info", xf);
      m->SetBit(TFileInfoMeta::kExternal);
      fc->AddMetaData(m);
   }

   // Cleanup
   if (excl) delete excl;

   // Done
   return dsmap;
}

////////////////////////////////////////////////////////////////////////////////
/// Add's a meta data object to the file collection object. The object will be
/// adopted by the TFileCollection and should not be deleted by the user.
/// Typically objects of class TFileInfoMeta or derivatives should be added,
/// but any class is accepted.
///
/// NB : a call to TFileCollection::Update will remove these objects unless the
///      bit TFileInfoMeta::kExternal is set.
/// Returns kTRUE if successful, kFALSE otherwise.

Bool_t TFileCollection::AddMetaData(TObject *meta)
{
   if (meta) {
      if (!fMetaDataList) {
         fMetaDataList = new TList;
         fMetaDataList->SetOwner();
      }
      fMetaDataList->Add(meta);
      return kTRUE;
   }
   return kFALSE;
}
