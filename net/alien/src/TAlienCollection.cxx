// @(#)root/alien:$Id$
// Author: Andreas-Joachim Peters 9/5/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TAlienCollection                                                     //
//                                                                      //
// Class which manages collection of files on AliEn middleware.         //
// The file collection is in the form of an XML file.                   //
//                                                                      //
// The internal list is managed as follows:                             //
// TList* ===> TMap*(file) ===> TMap*(attributes)                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TAlienCollection.h"
#include "TAlienResult.h"
#include "TAlienFile.h"
#include "TDSet.h"
#include "TMap.h"
#include "TFile.h"
#include "TXMLEngine.h"
#include "TObjString.h"
#include "TEntryList.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TError.h"
#include "TFileCollection.h"
#include "TFileInfo.h"
#include "snprintf.h"
#include <cstdlib>

ClassImp(TAlienCollection);

////////////////////////////////////////////////////////////////////////////////
/// Create Alien event collection by reading collection from the specified XML file.
/// You can restrict the number of importet entries using the maxentries value.

TAlienCollection::TAlienCollection(const char *localcollectionfile,
                                   UInt_t maxentries)
{
   fXmlFile = localcollectionfile;
   fFileGroupList = new TList();
   fFileGroupList->SetOwner(kTRUE);
   fFileGroupListIter = new TIter(fFileGroupList);
   fCurrent = 0;
   fNofGroups = 0;
   fNofGroupfiles = 0;
   fHasSUrls = kFALSE;
   fHasSelection = kFALSE;
   fHasOnline = kFALSE;
   fFileStager = 0;
   fExportUrl = "";
   fInfoComment = "";
   fCollectionName = "unnamed";
   fTagFilterList = 0;
   if (localcollectionfile != 0) {
      ParseXML(maxentries);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create Alien event collection using an event list.

TAlienCollection::TAlienCollection(TList *eventlist, UInt_t nofgroups,
                                   UInt_t nofgroupfiles)
{
   fFileGroupList = eventlist;
   fFileGroupList->SetOwner(kTRUE);
   fFileGroupListIter = new TIter(fFileGroupList);
   fCurrent = 0;
   fNofGroups = nofgroups;
   fNofGroupfiles = nofgroupfiles;
   fHasSUrls = kFALSE;
   fHasSelection = kFALSE;
   fHasOnline = kFALSE;
   fFileStager = 0;
   fExportUrl = "";
   fInfoComment = "";
   fCollectionName = "unnamed";
   fTagFilterList = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Clean up event file collection.

TAlienCollection::~TAlienCollection()
{
   if (fFileGroupList)
      delete fFileGroupList;

   if (fFileGroupListIter)
      delete fFileGroupListIter;

   if (fTagFilterList)
      delete fTagFilterList;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method used to create an Alien event collection, by reading
/// an XML collection from the specified url. All ROOT URLs are supported.
/// You can restrict the number of importet entries using the maxentries value

TGridCollection *TAlienCollection::Open(const char *collectionurl,
                                        UInt_t maxentries)
{
   //! stage the url to gSystem->TempDirectory() for remote url

   TString coll(collectionurl);
   Bool_t isRemote = coll.Contains(":/") && !coll.Contains("file:");
   if (isRemote) {
      TUUID uuid;
      coll = gSystem->TempDirectory();
      coll += "/aliencollection.";
      coll += uuid.AsString();
      if (!TFile::Cp(collectionurl, coll.Data())) {
         ::Error("TAlienCollection::Open", "Cannot make a local copy of collection with url %s",
                 collectionurl);
         return 0;
      }
   }

   TAlienCollection *collection =
       new TAlienCollection(coll, maxentries);

   if (isRemote && gSystem->Unlink(coll.Data())) {
      ::Error("TAlienCollection::Open", "Cannot remove the local copy of the collection %s",
              coll.Data());
   }

   return dynamic_cast <TAlienCollection * > (collection);
}

////////////////////////////////////////////////////////////////////////////////
/// Open the file specified by <filename> from the currently active file group in the collection via its TURL.

TFile *TAlienCollection::OpenFile(const char *filename)
{
   const char *turl = GetTURL(filename);
   if (turl) {
      return TFile::Open(turl);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method used to create an Alien event collection, by creating
/// collection from a TGridResult GetCollection result (TAlien::GetCollection)

TAlienCollection *TAlienCollection::OpenAlienCollection(TGridResult * queryresult,
                                             Option_t* /*option*/)
{
   if (!queryresult)
      return 0;

   TList *filelist = new TList();

   TIterator* fileiter = queryresult->MakeIterator();
   TMap *filemap = 0;

   while ((filemap = ((TMap *) fileiter->Next()))) {
      if (!filemap->GetValue("origLFN"))
         continue;

      filemap->Add(new TObjString("lfn"), new TObjString(filemap->GetValue("origLFN")->GetName()));
      filemap->Add(new TObjString("turl"), new TObjString(Form("alien://%s", filemap->GetValue("origLFN")->GetName())));
      TString bname = gSystem->BaseName(filemap->GetValue("origLFN")->GetName());
      filemap->Add(new TObjString("name"), new TObjString(bname.Data()));

      TMap* filegroup = new TMap();
      filegroup->Add(new TObjString(bname.Data()), filemap);
      filegroup->Add(new TObjString(""), filemap);

      // store the filegroup
      filelist->Add(filegroup);
   }
   delete fileiter;

   return new TAlienCollection(filelist, filelist->GetEntries(), 1);;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method used to create an Alien event collection, by creating
/// collection from a TGridResult Query result (TAlien::Query)
/// nogrouping means that files in the same directory are treated as not belonging to a file group

TGridCollection *TAlienCollection::OpenQuery(TGridResult * queryresult,
                                             Bool_t nogrouping)
{
   UInt_t nofgroups = 0;
   UInt_t nofgroupfiles = 0;

   if (!queryresult) {
      return 0;
   }

   TIterator *fileiter = queryresult->MakeIterator();
   TMap *filemap = 0;
   TString prev_bname = "";
   TString prev_dname = "";
   TMap *filegroup = 0;

   TList *filelist = new TList();

   while ((filemap = ((TMap *) fileiter->Next()))) {
      if (!filemap->GetValue("lfn"))
         continue;
      TString dname =
          gSystem->GetDirName(filemap->GetValue("lfn")->GetName());
      TString bname =
          gSystem->BaseName(filemap->GetValue("lfn")->GetName());
      filemap->Add(new TObjString("name"), new TObjString(bname.Data()));
      if ((!nogrouping)
          && ((bname != prev_bname) && (dname == prev_dname)
              && (filegroup))) {
         // add to group
         filegroup->Add(new TObjString(bname.Data()), filemap);
         if (nofgroups == 0)
            nofgroupfiles++;
      } else {
         // new group
         if (filegroup) {
            // store the old filegroup
            filelist->Add(filegroup);
            nofgroups++;
         }
         if (nofgroups == 0)
            nofgroupfiles++;
         filegroup = new TMap();
         filegroup->Add(new TObjString(bname.Data()), filemap);
         filegroup->Add(new TObjString(""), filemap);
      }

      prev_bname = bname;
      prev_dname = dname;
   }
   if (filegroup) {
      nofgroups++;
      filelist->Add(filegroup);
   }
   // set tag filter list
   TList *filterlist = new TList();
   filterlist->Add(new TObjString("type"));
   filterlist->Add(new TObjString("dir"));
   filterlist->Add(new TObjString("perm"));
   filterlist->Add(new TObjString("owner"));
   filterlist->Add(new TObjString("ctime"));
   filterlist->Add(new TObjString("seStringlist"));
   filterlist->Add(new TObjString("aclId"));
   filterlist->Add(new TObjString("expiretime"));
   filterlist->Add(new TObjString("replicated"));
   filterlist->Add(new TObjString("entryId"));
   filterlist->Add(new TObjString("gowner"));
   filterlist->Add(new TObjString("selist"));
   filterlist->Add(new TObjString("select"));
   filterlist->Add(new TObjString("online"));


   TAlienCollection *newcollection =
       new TAlienCollection(filelist, nofgroups, nofgroupfiles);
   if (newcollection) {
      newcollection->SetTagFilterList(filterlist);
   }
   return  dynamic_cast <TAlienCollection *> (newcollection);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset file iterator.

void TAlienCollection::Reset()
{
   fFileGroupListIter->Reset();
   fCurrent = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Parse event file collection XML file.
/// <maxentries> stops the parsing after <maxentries>.

void TAlienCollection::ParseXML(UInt_t maxentries)
{
   TXMLEngine xml;
   UInt_t parsedentries = 0;

   XMLDocPointer_t xdoc = xml.ParseFile(fXmlFile);
   if (!xdoc) {
      Error("ParseXML", "cannot parse the xml file %s", fXmlFile.Data());
      return;
   }

   XMLNodePointer_t xalien = xml.DocGetRootElement(xdoc);
   if (!xalien) {
      Error("ParseXML", "cannot find the <alien> tag in %s",
            fXmlFile.Data());
      return;
   }

   XMLNodePointer_t xcollection = xml.GetChild(xalien);
   if (!xcollection) {
      Error("ParseXML", "cannot find the <collection> tag in %s",
            fXmlFile.Data());
      return;
   }

   if (xml.GetAttr(xcollection,"name")) {
      fCollectionName = TString(xml.GetAttr(xcollection,"name"));
   } else {
      fCollectionName = ("unnamed");
   }

   XMLNodePointer_t xevent = xml.GetChild(xcollection);;
   if (!xevent) {
      Error("ParseXML", "cannot find the <event> tag in %s",
            fXmlFile.Data());
      return;
   }

   fNofGroups = 0;
   fNofGroupfiles = 0;
   fHasSUrls = kFALSE;
   fHasOnline = kFALSE;

   do {
      if (TString(xml.GetNodeName(xevent)) == "event") {
         parsedentries++;
         fNofGroups++;
         TMap *files = new TMap();

         // here is our event
         //      printf("Found event: %s\n",xml.GetAttr(xevent,"name"));

         // files
         XMLNodePointer_t xfile = xml.GetChild(xevent);
         if (!xfile)
            continue;

         Bool_t firstfile = kTRUE;
         do {
            // here we have an event file
            // get the attributes;
            xml.GetAttr(xfile, "lfn");
            xml.GetAttr(xfile, "turl");

            TMap *attributes = new TMap();
            TObjString *oname = new TObjString(xml.GetAttr(xfile, "name"));
            TObjString *oturl = new TObjString(xml.GetAttr(xfile, "turl"));
            TObjString *olfn = new TObjString(xml.GetAttr(xfile, "lfn"));
            TObjString *omd5 = new TObjString(xml.GetAttr(xfile, "md5"));
            TObjString *osize = new TObjString(xml.GetAttr(xfile, "size"));
            TObjString *oguid = new TObjString(xml.GetAttr(xfile, "guid"));
            TObjString *osurl = new TObjString(xml.GetAttr(xfile, "surl"));
            TObjString *osselect =
                new TObjString(xml.GetAttr(xfile, "select"));
            TObjString *ossexporturl =
                new TObjString(xml.GetAttr(xfile, "exporturl"));
            TObjString *osonline =
                new TObjString(xml.GetAttr(xfile, "online"));

            TObjString *oseStringlist =
                new TObjString(xml.GetAttr(xfile, "seStringlist"));
            TObjString *oevlist =
                new TObjString(xml.GetAttr(xfile, "evlist"));
            // if oevlist is defined, we parse it and fill a TEntyList
            if (oevlist && strlen(oevlist->GetName())) {
               TEntryList *xmlentrylist =
                   new TEntryList(oturl->GetName(), oguid->GetName());
               TString stringevlist = oevlist->GetName();
               TObjArray *evlist = stringevlist.Tokenize(",");
               for (Int_t n = 0; n < evlist->GetEntries(); n++) {
                  xmlentrylist->
                      Enter(atol
                            (((TObjString *) evlist->At(n))->GetName()));
               }
               attributes->Add(new TObjString("evlist"), xmlentrylist);
            }
            attributes->Add(new TObjString("name"), oname);
            attributes->Add(new TObjString("turl"), oturl);
            attributes->Add(new TObjString("lfn"), olfn);
            attributes->Add(new TObjString("md5"), omd5);
            attributes->Add(new TObjString("size"), osize);
            attributes->Add(new TObjString("guid"), oguid);
            attributes->Add(new TObjString("seStringlist"), oseStringlist);

            if (osurl && strlen(osurl->GetName())) {
               attributes->Add(new TObjString("surl"), osurl);
               fHasSUrls = kTRUE;
            }
            if (osselect && strlen(osselect->GetName())) {
               attributes->Add(new TObjString("select"), osselect);
               fHasSelection = kTRUE;
            }
            if (osonline && strlen(osonline->GetName())) {
               attributes->Add(new TObjString("online"), osonline);
               fHasOnline = kTRUE;
            }

            if (ossexporturl && strlen(ossexporturl->GetName())) {
               attributes->Add(new TObjString("exporturl"), ossexporturl);
               fExportUrl = ossexporturl->GetName();
            }

            files->Add(new TObjString(xml.GetAttr(xfile, "name")),
                       attributes);

            // we add the first file always as a file without name to the map
            if (firstfile) {
               files->Add(new TObjString(""), attributes);
               firstfile = kFALSE;
            }
            if (fNofGroups == 1)
               fNofGroupfiles++;
         } while ((xfile = xml.GetNext(xfile)));
         fFileGroupList->Add(files);
      }

      if (TString(xml.GetNodeName(xevent)) == "info") {
         if (xml.GetAttr(xevent,"comment")) {
            fInfoComment = TString(xml.GetAttr(xevent,"comment"));
         } else {
            fInfoComment = "";
         }
      }

      if (TString(xml.GetNodeName(xevent)) == "export") {
         if (xml.GetAttr(xevent,"url")) {
            SetExportUrl(xml.GetAttr(xevent,"url"));
         } else {
            fExportUrl = "";
         }
      }
      if (parsedentries >= maxentries)
         return;
   } while ((xevent = xml.GetNext(xevent)));
}

////////////////////////////////////////////////////////////////////////////////
/// Return next event file map.

TMap *TAlienCollection::Next()
{
   fCurrent = (TMap *) fFileGroupListIter->Next();
   return fCurrent;
}

////////////////////////////////////////////////////////////////////////////////
/// Prints statistics, selection and status about the loaded collection.

void TAlienCollection::Status()
{
   TIter *statuslist = new TIter(fFileGroupList);
   statuslist->Reset();
   TMap *oldcurrent = fCurrent;
   TMap *nextgroup = 0;
   UInt_t gc;
   UInt_t fc;
   UInt_t onlinegroups;
   UInt_t offlinegroups;
   UInt_t onlinefiles;
   UInt_t offlinefiles;
   UInt_t totalfiles;
   Long64_t totalfilesize;
   Long64_t onlinefilesize;
   Long64_t offlinefilesize;

   gc = 0;
   fc = 0;
   onlinegroups = offlinegroups = onlinefiles = offlinefiles = totalfiles =
       totalfilesize = onlinefilesize = offlinefilesize = 0;
   while ((nextgroup = (TMap *) statuslist->Next())) {
      gc++;
      fc = 0;
      Bool_t online;
      Bool_t selected;
      online = kTRUE;
      selected = kFALSE;
      TMap *attributes;
      TIterator *nextfile = nextgroup->MakeIterator();
      nextfile->Reset();
      while ((attributes = (TMap *) nextfile->Next())) {
         fCurrent = nextgroup;
         if (TString(attributes->GetName()) != "") {
            totalfiles++;
            totalfilesize += GetSize(attributes->GetName());
         }
         if (IsSelected(attributes->GetName())) {
            selected = kTRUE;
            if (TString(attributes->GetName()) != "") {
               fc++;
               fCurrent = nextgroup;
               if (!IsOnline(attributes->GetName())) {
                  online = kFALSE;
                  offlinefiles++;
                  offlinefilesize += GetSize(attributes->GetName());
               } else {
                  onlinefiles++;
                  onlinefilesize += GetSize(attributes->GetName());
               }
            }
         }
      }
      if (selected) {
         if (online)
            onlinegroups++;
         else
            offlinegroups++;
      }
   }
   fCurrent = oldcurrent;
   Info("Status", "=========================================");
   Info("Status", "     Tot. Number of files: %u", totalfiles);
   Info("Status", "                Tot. Size: %0.2f GB",
        totalfilesize / 1024.0 / 1024.0 / 1024.0);
   Info("Status", "    Number of file groups: %u", gc);
   Info("Status", "Number of files per group: %u", fc);
   Info("Status", "-----------------------------------------");
   Info("Status", "Online (staged [selected]):");
   Info("Status", "    Number of file groups: %u", onlinegroups);
   Info("Status", "          Number of files: %u", onlinefiles);
   Info("Status", "                     Size: %0.2f GB",
        onlinefilesize / 1024.0 / 1024.0 / 1024.0);
   Info("Status", "           Fraction avail: %0.2f %%",
        100.0 * onlinefilesize / (onlinefilesize + offlinefilesize +
                                  0.0000001));
   Info("Status", "-----------------------------------------");
   Info("Status", "Offline (to be staged [selected]):");
   Info("Status", "    Number of file groups: %u", offlinegroups);
   Info("Status", "          Number of files: %u", offlinefiles);
   Info("Status", "                     Size: %0.2f GB",
        offlinefilesize / 1024.0 / 1024.0 / 1024.0);
   Info("Status", "            Fraction miss: %0.2f %%",
        100.0 * offlinefilesize / (onlinefilesize + offlinefilesize +
                                   0.0000001));
   Info("Status", "=========================================\n");

   delete statuslist;
}

////////////////////////////////////////////////////////////////////////////////
/// Set's a key value pair in a tagmap. If it is existing, the existing tag is overwritten. If not, it is created.

void TAlienCollection::SetTag(const char *tag, const char *value,
                              TMap * tagmap)
{
   if ((!tag) || (!value) || (!tagmap)) {
      return;
   }
   TObject *delobject = tagmap->FindObject(tag);
   if (delobject) {
      TObject *keyobject = ((TPair *) delobject)->Key();
      tagmap->Remove(keyobject);;
   }
   tagmap->Add(new TObjString(tag), new TObjString(value));
}


////////////////////////////////////////////////////////////////////////////////
/// Selects all files with name <filename> in the collection
/// All files can be selected using "*" as filename

Bool_t TAlienCollection::SelectFile(const char *filename, Int_t nstart, Int_t nstop)
{
   Int_t cnt=0;
   fHasSelection = kTRUE;
   Reset();
   TMap *nextgroup;
   while ((nextgroup = (TMap *) Next())) {
      cnt++;
      TMap *attributes;
      TIterator *nextfile = nextgroup->MakeIterator();
      nextfile->Reset();
      if ( ((nstart == -1 ) && (nstop == -1)) ||
           ((nstart != -1 ) && (cnt >= nstart) && (nstop == -1)) ||
           ((nstart != -1 ) && (cnt >= nstart) && (nstop != -1) && (cnt <= nstop)) ||
           ((nstop  != -1 ) && (cnt <= nstop)  && (nstart == -1))) {
         while ((attributes = (TMap *) nextfile->Next())) {
            if (TString(attributes->GetName()) != "") {
               if ((TString(attributes->GetName()) == TString(filename)) ||
                   (TString(filename) == TString("*"))) {
                  SetTag("select", "1",
                  ((TMap *) nextgroup->GetValue(attributes->GetName())));
               }
            }
         }
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Deselects the file <filename> from the loaded collection
/// All files can be deselected using "*" as <filename>

Bool_t TAlienCollection::DeselectFile(const char *filename, Int_t nstart, Int_t nstop)
{
   Int_t cnt=0;
   Reset();
   TMap *nextgroup;
   while ((nextgroup = (TMap *) Next())) {
      cnt++;
      TMap *attributes;
      TIterator *nextfile = nextgroup->MakeIterator();
      nextfile->Reset();
      if ( ((nstart == -1 ) && (nstop == -1)) ||
           ((nstart != -1 ) && (cnt >= nstart) && (nstop == -1)) ||
           ((nstart != -1 ) && (cnt >= nstart) && (nstop != -1) && (cnt <= nstop)) ||
           ((nstop  != -1 ) && (cnt <= nstop)  && (nstart == -1))) {
         while ((attributes = (TMap *) nextfile->Next())) {
            if (TString(attributes->GetName()) != "") {
               if ((TString(attributes->GetName()) == TString(filename)) ||
                   (TString(filename) == TString("*"))) {
                  SetTag("select", "0",
                  ((TMap *) nextgroup->GetValue(attributes->GetName())));
               }
            }
         }
      }
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Invert the selection.

Bool_t TAlienCollection::InvertSelection()
{
   Int_t cnt=0;
   fHasSelection = kTRUE;
   Reset();
   TMap *nextgroup;
   while ((nextgroup = (TMap *) Next())) {
      cnt++;
      TMap *attributes;
      TIterator *nextfile = nextgroup->MakeIterator();
      nextfile->Reset();
      while ((attributes = (TMap *) nextfile->Next())) {
         if (IsSelected(attributes->GetName())) {
            SetTag("select", "0",
               ((TMap *) nextgroup->GetValue(attributes->GetName())));
         } else {
            SetTag("select", "1",
               ((TMap *) nextgroup->GetValue(attributes->GetName())));
         }
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// downscales the selection with scaler

Bool_t TAlienCollection::DownscaleSelection(UInt_t scaler)
{
   Int_t cnt = 0;

   Reset();
   TMap *nextgroup;
   while ((nextgroup = (TMap *) Next())) {
      cnt++;
      TMap *attributes;
      TIterator *nextfile = nextgroup->MakeIterator();
      nextfile->Reset();
      if (cnt%scaler) {
         while ((attributes = (TMap *) nextfile->Next())) {
            if (TString(attributes->GetName()) != "") {
               SetTag("select", "0",
                  ((TMap *) nextgroup->GetValue(attributes->GetName())));
            }
         }
      }
   }
   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Return next event file map.

Bool_t TAlienCollection::Remove(TMap * map)
{
   if (fFileGroupList->Remove(map)) {
      return kTRUE;
   } else {
      return kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get a file's transport URL (TURL). Returns 0 in case of error.

const char *TAlienCollection::GetTURL(const char *filename)
{
   if (fCurrent) {
      TMap *obj = (TMap *) fCurrent->GetValue(filename);
      if (obj) {
         if (obj->GetValue("turl")) {
            return (((TObjString *) obj->GetValue("turl"))->GetName());
         }
      }
   }
   Error("GetTURL", "cannot get TURL of file %s", filename);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a file's storage URL (SURL). Returns 0 in case of error.

const char *TAlienCollection::GetSURL(const char *filename)
{
   if (fCurrent) {
      TMap *obj = (TMap *) fCurrent->GetValue(filename);
      if (obj) {
         if (obj->GetValue("surl")) {
            return (((TObjString *) obj->GetValue("surl"))->GetName());
         }
      }
   }
   Error("GetSURL", "cannot get SURL of file %s", filename);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a file's online tag. Returns false if not online or CheckIfOnline was never executed, true if online

Bool_t TAlienCollection::IsOnline(const char *filename)
{
   if (fCurrent) {
      TMap *obj = (TMap *) fCurrent->GetValue(filename);
      if (obj) {
         if (obj->GetValue("online")) {
            TString online =
                (((TObjString *) obj->GetValue("online"))->GetName());
            if (online == "1") {
               return kTRUE;
            } else {
               return kFALSE;
            }
         }
      }
   }
   //   Error("IsOnline", "cannot get online tag of file %s", filename);
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a file's online tag. Returns false if not online or CheckIfOnline was never executed, true if online

Bool_t TAlienCollection::IsSelected(const char *filename)
{
   if (fCurrent) {
      TMap *obj = (TMap *) fCurrent->GetValue(filename);
      if (obj) {
         if (obj->GetValue("select")) {
            TString selected;
            selected =
                (((TObjString *) obj->GetValue("select"))->GetName());
            if (selected == TString("1")) {
               return kTRUE;
            } else {
               return kFALSE;
            }
         }
      }
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a file's event list. Returns 0 in case of error.

TEntryList *TAlienCollection::GetEntryList(const char *filename)
{
   if (fCurrent) {
      TMap *obj = (TMap *) fCurrent->GetValue(filename);
      if (obj) {
         if (obj->GetValue("evlist")) {
            return ((TEntryList *) obj->GetValue("evlist"));
         }
      }
   }
   Error("GetEntryList", "cannot get evelist of file %s", filename);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a file's LFN. Returns 0 in case of error.

const char *TAlienCollection::GetLFN(const char *filename)
{
   if (fCurrent) {
      TMap *obj = (TMap *) fCurrent->GetValue(filename);
      if (obj) {
         if (obj->GetValue("lfn")) {
            return (((TObjString *) obj->GetValue("lfn"))->GetName());
         }
      }
   }
   Error("GetLFN", "cannot get LFN");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a file's LFN. Returns 0 in case of error.

Long64_t TAlienCollection::GetSize(const char *filename)
{
   if (fCurrent) {
      TMap *obj = (TMap *) fCurrent->GetValue(filename);
      if (obj) {
         if (obj->GetValue("size")) {
            TString ssize =
                (((TObjString *) obj->GetValue("size"))->GetName());
            return ssize.Atoll();
         }
      }
   }
   Error("GetSize", "cannot get size of %s", filename);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Print event file collection.

void TAlienCollection::Print(Option_t *) const
{
   Info("Print", "dumping %d elements", fFileGroupList->GetSize());
   TIter next(fFileGroupList);
   TMap *filemap;
   Int_t count = 0;
   while ((filemap = (TMap *) next())) {
      count++;
      Info("Print", "printing element %d", count);
      filemap->Print();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return a TDSet from a collection. Filter, Selections and online tags are not
/// taken into account.

TDSet *TAlienCollection::GetDataset(const char *type, const char *objname,
                                    const char *dir)
{
   Reset();
   TDSet *dset = new TDSet(type, objname, dir);
   if (!dset) {
      return 0;
   }

   while (Next()) {
      if (((TObjString *) fCurrent->GetValue("")))
         dset->Add(((TMap *) (fCurrent->GetValue("")))->GetValue("turl")->
                   GetName());;
   }
   return dset;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a TGridResult.
/// If files have been selected in this collection, only the selected files
/// are returned. If <onlyonline> is set to kTRUE, only files which are
/// 'online' (staged) are included. If no online check was done, TGridResult
/// will be empty. <publicaccess> adds the publicaccess option to the TGridResult entries

TGridResult *TAlienCollection::GetGridResult(const char *filename,
                                             Bool_t onlyonline,
                                             Bool_t publicaccess)
{
   Reset();
   TGridResult *result = new TAlienResult();

   while (Next()) {
      if (((TObjString *) fCurrent->GetValue(filename))) {
         TMap *attributes = (TMap *) fCurrent->GetValue(filename)->Clone();
         if (publicaccess) {
            attributes->Add(new TObjString("options"),
                            new TObjString("&publicaccess=1"));
         }
         if ((!fHasSelection) || (IsSelected(filename))) {
            if ((!onlyonline) || (fHasOnline && IsOnline(filename))) {
               result->Add(attributes);
            }
         }
      }
   }
   return dynamic_cast < TGridResult * >(result);
}

////////////////////////////////////////////////////////////////////////////////
/// return kTRUE if comparator overlaps with this
/// all objects in this collection, which are not defined in the <comparator> collection are removed.

Bool_t TAlienCollection::OverlapCollection(TGridCollection * comparator)
{
   if ((!comparator)) {
      return kFALSE;
   }

loopagain:
   // loop over col1 and try to find it in col2
   this->Reset();
   // loop over all elements in reference (=this)
   TMap *overlapmap;
   while ((overlapmap = this->Next())) {
      comparator->Reset();
      Bool_t found = kFALSE;
      // try to find in the comparator collection
      while ((comparator->Next())) {
         TString s1 = this->GetLFN();
         TString s2 = comparator->GetLFN();
         if (s1 == s2) {
            found = kTRUE;
            break;
         }
      }
      if (!found) {
         this->Remove(overlapmap);
         goto loopagain;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// adds <addcollection> to this collection - equal elements are skipped

void TAlienCollection::Add(TGridCollection * addcollection)
{
   if ((!addcollection)) {
      return;
   }
   // loop over col1 and try to find it in col2
   addcollection->Reset();
   // loop over all elements in reference (=this)
   TMap *addmap;
   while ((addmap = addcollection->Next())) {
      Reset();
      // try to find in the comparator collection
      TString s2 = addcollection->GetLFN();
      while (Next()) {
         TString s1 = GetLFN();
         // printf("%s = %s\n", s1.Data(), s2.Data());
         if (s1 == s2) {
            Error("Add",
                  "File group with lfn %s exists already in this collection - skipping",
                  GetLFN());
            goto leaveloop;
         }
      }
      TMap *clonemap;
      clonemap = (TMap *) addmap->Clone();
      fFileGroupList->Add(clonemap);
leaveloop:
      ;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// adds <addcollection> to this collection - NO check for identical elements

void TAlienCollection::AddFast(TGridCollection * addcollection)
{
   if ((!addcollection)) {
      return;
   }
   addcollection->Reset();
   TMap *addmap, *clonemap;
   while ((addmap = addcollection->Next())) {
      clonemap = (TMap *) addmap->Clone();
      fFileGroupList->Add(clonemap);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// retrieves all the SURLS for the LFNS

Bool_t TAlienCollection::LookupSUrls(Bool_t verbose)
{
   Bool_t ok = kTRUE;
   UInt_t lc = 0;
   Reset();
   TMap *filemap;
   while ((filemap = Next())) {
      // loop over all files in this map
      TIterator *nextfile = filemap->MakeIterator();
      TMap *attributes;
      while ((attributes = (TMap *) nextfile->Next())) {
         if (TString(attributes->GetName()) != "") {
            lc++;
            if (fHasSelection && (!IsSelected(attributes->GetName())))
               continue;
            // there is always an "" entry in the map to point to the first file of a file group
            if (verbose)
               Info("LookupSUrls", "Lookup SURL for %s [%u/%u]",
                    GetTURL(attributes->GetName()), lc,
                    fNofGroups * fNofGroupfiles);
            TString surl =
                TAlienFile::SUrl(GetTURL(attributes->GetName()));
            if (!surl.Data()) {
               ok = kFALSE;
            } else {
               // delete the surl map entry, if it exists already
               TObject *delobject =
                   ((TMap *) filemap->GetValue(attributes->GetName()))->
                   FindObject("surl");
               if (delobject) {
                  TObject *keyobject = ((TPair *) delobject)->Key();
                  ((TMap *) filemap->GetValue(attributes->GetName()))->
                      Remove(keyobject);;
               }

               ((TMap *) filemap->GetValue(attributes->GetName()))->
                   Add(new TObjString("surl"),
                       new TObjString(surl.Data()));
               if (verbose)
                  Info("LookupSUrls", "SURL = %s", surl.Data());

            }

         }
      }
   }
   if (ok) {
      fHasSUrls = kTRUE;
   }
   return ok;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute the 'stage' method for all files in this collection (trigger staging).
/// The <bulk> stage method is currently not working.

Bool_t TAlienCollection::Stage(Bool_t bulk,Option_t* option)
{
   if (!fHasSUrls) {
      Error("Stage",
            "You have to execute LookupSUrls() before you can stage this collection");
      return kFALSE;
   }
   UInt_t fc = 0;
   Reset();
   TMap *filemap;

   if (!bulk) {
      while ((filemap = Next())) {
         // loop over all files in this map
         TIterator *nextfile = filemap->MakeIterator();
         TMap *attributes;
         while ((attributes = (TMap *) nextfile->Next())) {
            if (TString(attributes->GetName()) != "") {
               fc++;
               if (fHasSelection && (!IsSelected(attributes->GetName())))
                  continue;

               if (!fFileStager) {
                  fFileStager =
                      TFileStager::Open(GetSURL(attributes->GetName()));
               }
               if ((fFileStager)->Stage(GetSURL(attributes->GetName()),option)) {
                  // file staged
                  Info("Stage", "[%05u/%05u] <Staged> : %s", fc,
                       GetNofGroups() * GetNofGroupfiles(),
                       GetLFN(attributes->GetName()));
               } else {
                  // file stage failed
                  Error("Stage", "[%05u/%05u] <Failed to stage> : %s",
                        fc, GetNofGroups() * GetNofGroupfiles(),
                        GetLFN(attributes->GetName()));
               }
            }
         }
      }
   } else {
      // bulk request
      TList* stagelist = new TList();
      stagelist->SetOwner(kTRUE);
      Bool_t stageresult=kFALSE;
      Reset();
      while ((filemap = Next())) {
         TIterator *nextfile = filemap->MakeIterator();
         TMap *attributes;
         while ((attributes = (TMap *) nextfile->Next())) {
            if (TString(attributes->GetName()) != "") {
               fc++;
               stagelist->Add( new TUrl((GetSURL(attributes->GetName()))));
            }
         }
      }

      if (fc) {
         if (!fFileStager) {
            fFileStager = TFileStager::Open(stagelist->First()->GetName());
         }

         stageresult = (fFileStager)->Stage(stagelist,option);
      }
      delete stagelist;
      return stageresult;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Run an online check if files are currently accessible (staged) or offline (to be staged).
/// The <bulk> check is currently not working.

Bool_t TAlienCollection::CheckIfOnline(Bool_t bulk)
{
   if (!fHasSUrls) {
      Error("CheckIfOnline",
            "You have to execute LookupSUrls() before you can prepare this collection");
      return kFALSE;
   }

   Reset();
   TMap *filemap;
   UInt_t fc=0;

   if (!bulk) {
      while ((filemap = Next())) {
         // loop over all files in this map
         TIterator *nextfile = filemap->MakeIterator();
         TMap *attributes;
         while ((attributes = (TMap *) nextfile->Next())) {
            if (fHasSelection && (!IsSelected(attributes->GetName())))
               continue;
            if (TString(attributes->GetName()) != "") {
               fc++;
               // check if we have a fFileStager
               if (!fFileStager) {
                  fFileStager =
                      TFileStager::Open(GetSURL(attributes->GetName()));
               }
               Bool_t online = kFALSE;
               if ((fFileStager)->IsStaged(GetSURL(attributes->GetName()))) {
                  // file is online
                  Info("CheckIfOnline", "[%05u/%05u] <Online> : %s", fc,
                       GetNofGroups() * GetNofGroupfiles(),
                       GetLFN(attributes->GetName()));
                  online = kTRUE;
               } else {
                  // file is offline
                  Info("CheckIfOnline", "[%05u/%05u]<Offline> : %s", fc,
                       GetNofGroups() * GetNofGroupfiles(),
                       GetLFN(attributes->GetName()));
                  online = kFALSE;
               }

               {
                  // set the online tag in the collection
                  // delete the online tag entry, if it exists already
                  TObject *delobject =
                      ((TMap *) filemap->GetValue(attributes->GetName()))->
                      FindObject("online");
                  if (delobject) {
                     TObject *keyobject = ((TPair *) delobject)->Key();
                     ((TMap *) filemap->GetValue(attributes->GetName()))->
                         Remove(keyobject);;
                  }
                  if (online)
                     ((TMap *) filemap->GetValue(attributes->GetName()))->
                         Add(new TObjString("online"),
                             new TObjString("1"));
                  else
                     ((TMap *) filemap->GetValue(attributes->GetName()))->
                         Add(new TObjString("online"),
                             new TObjString("0"));
               }
            }
         }
      }
   } else {
      // bulk lookup
      TList *lookuplist = new TList();
      if (lookuplist) {
         lookuplist->SetOwner(kTRUE);
         while ((filemap = Next())) {
            // loop over all files in this map
            TIterator *nextfile = filemap->MakeIterator();
            TMap *attributes;
            while ((attributes = (TMap *) nextfile->Next())) {
               if (TString(attributes->GetName()) != "") {
                  fc++;
                  // check if we have a fFileStager
                  if (!fFileStager) {
                     fFileStager =
                       TFileStager::Open(GetSURL(attributes->GetName()));
                  }
                  lookuplist->
                      Add(new TObjString(GetSURL(attributes->GetName())));
               }
            }
         }

         TList *onlinelist = fFileStager->GetStaged(lookuplist);
         if (!onlinelist) {
            SafeDelete(lookuplist);
            return kFALSE;
         }

         Reset();
         fc=0;
         while ((filemap = Next())) {
            // loop over all files in this map
            TIterator *nextfile = filemap->MakeIterator();
            TMap *attributes;
            while ((attributes = (TMap *) nextfile->Next())) {
               if (TString(attributes->GetName()) != "") {
                  fc++;
                  Bool_t online;
                  // check if it is in the online list
                  if (onlinelist->
                      FindObject(GetSURL(attributes->GetName()))) {
                     // this file is online
                     Info("CheckIfOnline", "[%05u/%05u] <Online> : %s", fc,
                          GetNofGroups() * GetNofGroupfiles(),
                          GetLFN(attributes->GetName()));
                     online = kTRUE;
                  } else {
                     // this file is offline
                     Info("CheckIfOnline", "[%05u/%05u]<Offline> : %s", fc,
                          GetNofGroups() * GetNofGroupfiles(),
                          GetLFN(attributes->GetName()));
                     online = kFALSE;
                  }

                  {
                     // set the online tag in the collection
                     // delete the online tag entry, if it exists already
                     TObject *delobject =
                         ((TMap *) filemap->
                         GetValue(attributes->GetName()))->
                         FindObject("online");
                     if (delobject) {
                        TObject *keyobject = ((TPair *) delobject)->Key();
                        ((TMap *) filemap->
                         GetValue(attributes->GetName()))->
                         Remove(keyobject);
                     }
                     if (online)
                        ((TMap *) filemap->
                         GetValue(attributes->GetName()))->
                         Add(new TObjString("online"), new TObjString("1"));
                     else
                        ((TMap *) filemap->
                         GetValue(attributes->GetName()))->
                         Add(new TObjString("online"), new TObjString("0"));
                  }
               }
            }
         }

         SafeDelete(onlinelist);
         SafeDelete(lookuplist);
      } else {
         fHasOnline = kFALSE;
         return kFALSE;
      }
   }

   fHasOnline = kTRUE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Exports the contens of the TAlienCollection into an XML formatted file.
/// By default exports only selected files. By default exports only accessible (online) files.
/// You can change this behaviour by specifying online=kFALSE or selected=kFALSE.
/// <name> specifies a name you want to assign to this collection.
/// <comment> can be a user comment to this collection.
/// If <exporturl>="" the collection is exported to the URL which was stored inside the collection or
/// was specified by the ExportUrl(const char* url) method.

Bool_t TAlienCollection::ExportXML(const char *exporturl, Bool_t selected, Bool_t online,
                                   const char *name , const char *comment)
{

   TFile *exportfile;

   if ((exporturl == 0) || (TString(exporturl) == "")) {
      if (!(exporturl = GetExportUrl())) {
         Error("ExportXML",
               "You used the option to store to the collection export url, but this is not defined!");
         return kFALSE;
      }
   }
   TUrl expfn(exporturl);
   TString options = expfn.GetOptions();
   if (options.Length()) {
      options += ",filetype=raw";
   } else {
      options = "filetype=raw";
   }
   expfn.SetOptions(options.Data());

   if (!(exportfile = TFile::Open(expfn.GetUrl(), "RECREATE"))) {
      Error("ExportXML", "Cannot open export URL %s", expfn.GetUrl());
      return kFALSE;
   }

   Bool_t expret = ExportXML(exportfile, selected, online, name, comment);
   exportfile->Close();
   return expret;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal Export function to write a collection as an XML file. See above.

Bool_t TAlienCollection::ExportXML(TFile * exportfile, Bool_t selected,
                                   Bool_t online, const char *name,
                                   const char *comment)
{
   char outline[4096];

   // write headers
   snprintf(outline,4096,
           "<?xml version=\"1.0\"?>\n<alien>\n  <collection name=\"%s\">\n",
           name);
   if ((exportfile->WriteBuffer(outline, strlen(outline)))) {
      Error("ExportXML", "Error writing XML to export file");
      exportfile->Close();
      return kFALSE;
   }
   // write events
   TMap *oldcurrent = fCurrent;
   TMap *nextgroup;
   TIter *filegroups = new TIter(fFileGroupList);
   filegroups->Reset();
   UInt_t groupcnt = 0;
   while ((nextgroup = (TMap *) filegroups->Next())) {
      Bool_t isselected;
      isselected = kFALSE;
      TMap *attributes;
      TIterator *nextfile = nextgroup->MakeIterator();
      nextfile->Reset();
      // check if something is selected
      while ((attributes = (TMap *) nextfile->Next())) {
         if (TString(attributes->GetName()) != "") {
            fCurrent = nextgroup;
            if (IsSelected(attributes->GetName())) {
               isselected = kTRUE;
            }

         }
      }
      if ((!selected) || isselected) {
         // loop again and export files
         nextfile->Reset();
         groupcnt++;
         // open new event header
         snprintf(outline,4096, "    <event name=\"%d\">\n", groupcnt);
         if ((exportfile->WriteBuffer(outline, strlen(outline)))) {
            Error("ExportXML", "Error writing XML to export file");
            exportfile->Close();
            return kFALSE;
         }
         while ((attributes = (TMap *) nextfile->Next())) {
            if (TString(attributes->GetName()) != "") {
               fCurrent = nextgroup;
               if ((!selected) || (IsSelected(attributes->GetName()))) {
                  if ((!online) || (IsOnline(attributes->GetName()))) {
                     // export this file
                     /////////////////////////////////////////////////////////////
                     // open file tag
                     snprintf(outline,4096, "      <file ");
                     if ((exportfile->
                          WriteBuffer(outline, strlen(outline)))) {
                        Error("ExportXML",
                              "Error writing XML to export file");
                        exportfile->Close();
                        return kFALSE;
                     }
                     // loop over map keys
                     TIterator *mkeys =
                         ((TMap *) nextgroup->
                          GetValue(attributes->GetName()))->MakeIterator();
                     mkeys->Reset();
                     TObjString *tagname = 0;
                     TObjString *tagval = 0;
                     while ((tagname = (TObjString *) mkeys->Next())) {
                        Bool_t filtered = kFALSE;
                        // check for filtered tags from the filter list
                        if (fTagFilterList) {
                           TIter next(fTagFilterList);
                           TObjString *filtertag;
                           while ((filtertag = (TObjString *) next())) {
                              if (TString(tagname->GetName()) ==
                                  TString(filtertag->GetName())) {
                                 filtered = kTRUE;
                                 break;
                              }
                           }
                        }
                        if (!filtered) {
                           tagval =
                               (TObjString *) ((TMap *) nextgroup->
                                               GetValue(attributes->
                                                        GetName()))->
                               GetValue(tagname);
                           if (TString(tagname->GetName()) != "evlist") {
                              snprintf(outline,4096, "%s=\"%s\" ",
                                      tagname->GetName(),
                                      tagval->GetName());
                           } else {
                              // the eventlist has to be converted from TEventList to a string list with komma separation
                              TEntryList *xmlentrylist =
                                  (TEntryList *) tagval;
                              if (!xmlentrylist)
                                 continue;
                              TString slist = "";
                              for (Int_t i = 0; i < xmlentrylist->GetN(); i++) {
                                 if (i > 0)
                                    slist += ",";
                                 slist += xmlentrylist->GetEntry(i);
                              }
                              snprintf(outline,4096, "%s=\"%s\" ",
                                      tagname->GetName(), slist.Data());
                           }

                           if ((exportfile->
                                WriteBuffer(outline, strlen(outline)))) {
                              Error("ExportXML",
                                    "Error writing XML to export file");
                              exportfile->Close();
                              return kFALSE;
                           }
                        }
                     }

                     // close file tag
                     snprintf(outline,4096, "/>\n");
                     if ((exportfile->
                          WriteBuffer(outline, strlen(outline)))) {
                        Error("ExportXML",
                              "Error writing XML to export file");
                        exportfile->Close();
                        return kFALSE;
                     }
                  }
               }
            }
         }

         // close event
         snprintf(outline,4096, "    </event>\n");
         if ((exportfile->WriteBuffer(outline, strlen(outline)))) {
            Error("ExportXML", "Error writing XML to export file");
            exportfile->Close();
            return kFALSE;
         }
      }
   }

   fCurrent = oldcurrent;

   // write export url if present
   if (GetExportUrl()) {
      snprintf(outline,4096, "    <export url=\"%s\">\n",GetExportUrl());
      if ((exportfile->WriteBuffer(outline, strlen(outline)))) {
         Error("ExportXML", "Error writing XML to export file");
         exportfile->Close();
         return kFALSE;
      }
   }


   // write trailer
   snprintf(outline,4096,
           "    <info comment=\"%s\" />\n</collection>\n</alien>\n",
           comment);
   if ((exportfile->WriteBuffer(outline, strlen(outline)))) {
      Error("ExportXML", "Error writing XML to export file");
      exportfile->Close();
      return kFALSE;
   }

   delete filegroups;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the 'default' export URL for an XML collection. A collection can be restored to the export URL using ExportXML("");

Bool_t TAlienCollection::SetExportUrl(const char *exporturl)
{
   if (exporturl)
      fExportUrl = exporturl;

   if (fExportUrl == "") {
      Info("ExportUrl",
           "There is no remote url defined in this collection");
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Adds to a file given by infile the collection identification , f.e.
/// for collection files sitting in directories like 100/1/AliESD.root
///                                                  ...
///                                                   110/1/AliESD.root
/// "./histo.root" will be converted to "./histo.100_1-110_1.230.root
///
/// The name syntax is <oldname>.<first run>_<first event>-<last run>.<last event>.<nevents>.root
///
/// By default the found files are renamed on the local disk
/// example:
/// - you specify f.e. as arguments GetOutputFileName("*.root",kTRUE);
///  --->> this renames all ROOT files corresponding to the collection contents

const char *TAlienCollection::GetOutputFileName(const char *infile,
                                                Bool_t rename, const char* suffix)
{
   Bool_t first = kTRUE;
   TString firstrun;
   TString firstevent;
   TString lastrun;
   TString lastevent;
   UInt_t nevents = 0;

   Reset();
   while (Next()) {
      TString s1 = gSystem->GetDirName(GetLFN(""));
      TString s2(gSystem->BaseName(s1.Data()));
      TString s3 = gSystem->GetDirName(s1.Data());
      TString s4(gSystem->BaseName(s3.Data()));
      nevents++;
      if (first) {
         first = kFALSE;
         firstevent = s2;
         firstrun = s4;
      }
      lastevent = s2;
      lastrun = s4;
   }

   // list the matching files
   TString pcmd("ls ");
   pcmd += infile;
   //printf("Pipe is %s\n",pcmd.Data());
   FILE *fp = gSystem->OpenPipe(pcmd.Data(), "r");
   if (fp) {
      char rootfile[4096];
      Int_t item;
      while ((item = fscanf(fp, "%4095s", rootfile)) == 1) {
         TString rootdir = gSystem->GetDirName(rootfile);
         TString rootbase(gSystem->BaseName(rootfile));
         TString rootbasenosuffix;
         rootbasenosuffix = rootbase(0, rootbase.First('.'));;
         // build new name like <oldname>.<firstrun>_<firstevent>-<lastrun>_<lastevent>.root
         TString newrootname;
         if (rootdir.Length()) {
            newrootname += rootdir;
            newrootname += "/";
         }
         newrootname += rootbasenosuffix;
         newrootname += ".";
         newrootname += firstrun;
         newrootname += "_";
         newrootname += firstevent;
         newrootname += "-";
         newrootname += lastrun;
         newrootname += "_";
         newrootname += lastevent;
         newrootname += ".";
         newrootname += nevents;
         newrootname += ".";
         newrootname += suffix;
         Info("GetOutputFilename", "Renaming File %s to %s", rootfile,
              newrootname.Data());
         fLastOutFileName = newrootname;
         if (rename) {
            gSystem->Rename(rootfile, newrootname.Data());
         }
      }
      gSystem->ClosePipe(fp);
   }
   return fLastOutFileName.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// creates a TFileCollection objects and fills it with the information from this collection
/// note that TFileCollection has a flat structure and no groups --> all files are filles on a flat level
/// the TFileInfo of each file in the TFileCollection is filled with turl, size, md5, guid
///
/// the object has to be deleted by the user

TFileCollection *TAlienCollection::GetFileCollection(const char* name, const char* title) const
{
   TFileCollection* collection = new TFileCollection(name, title);

   TIter next(fFileGroupList);
   TMap* group = 0;
   while ((group = dynamic_cast<TMap*>(next()))) {
      TIter next2(group);
      TObjString* key = 0;
      while ((key = dynamic_cast<TObjString*> (next2()))) {
         if (key->String().Length() == 0)
            continue;

         TMap* file = dynamic_cast<TMap*> (group->GetValue(key));
         if (!file)
            continue;

         TObjString* turl = dynamic_cast<TObjString*> (file->GetValue("turl"));
         TObjString* size = dynamic_cast<TObjString*> (file->GetValue("size"));
         TObjString* md5 = dynamic_cast<TObjString*> (file->GetValue("md5"));
         TObjString* guid = dynamic_cast<TObjString*> (file->GetValue("guid"));

         if (!turl || turl->String().Length() == 0)
            continue;

         TFileInfo* fileInfo = new TFileInfo(turl->String(), size->String().Atoi(), guid->String(), md5->String());
         collection->Add(fileInfo);
      }
   }

   collection->Update();

   return collection;
}
