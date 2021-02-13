// @(#)root/base:$Id$
// Author: Andreas-Joachim Peters   20/9/2005

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TFileInfo
\ingroup Base

Class describing a generic file including meta information.
*/

#include "TFileInfo.h"
#include "TRegexp.h"
#include "TSystem.h"
#include "TClass.h"
#include "TUrl.h"
#include "TUUID.h"
#include "TMD5.h"

ClassImp(TFileInfo);
ClassImp(TFileInfoMeta);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFileInfo::TFileInfo(const char *in, Long64_t size, const char *uuid,
                     const char *md5, TObject *meta)
   : fCurrentUrl(0), fUrlList(0), fSize(-1), fUUID(0), fMD5(0),
     fMetaDataList(0), fIndex(-1)
{
   // Get initializations form the input string: this will set at least the
   // current URL; but it may set more: see TFileInfo::ParseInput(). Please note
   // that MD5 sum should be provided as a string in md5ascii form.
   ParseInput(in);

   // Now also honour the input arguments: the size
   if (size > -1) fSize = size;
   // The UUID
   if (uuid) {
      SafeDelete(fUUID);
      fUUID = new TUUID(uuid);
   } else if (!fUUID) {
      fUUID = new TUUID;
   }
   // The MD5
   if (md5) {
      SafeDelete(fMD5);
      fMD5 = new TMD5();
      fMD5->SetDigest(md5);  // sets digest from md5ascii representation
   }
   // The meta information
   if (meta) {
      RemoveMetaData(meta->GetName());
      AddMetaData(meta);
   }

   // Now set the name from the UUID
   SetName(fUUID->AsString());
   SetTitle("TFileInfo");

   // By default we ignore the index
   ResetBit(TFileInfo::kSortWithIndex);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TFileInfo::TFileInfo(const TFileInfo &fi) : TNamed(fi.GetName(), fi.GetTitle()),
                                            fCurrentUrl(0), fUrlList(0),
                                            fSize(fi.fSize), fUUID(0), fMD5(0),
                                            fMetaDataList(0), fIndex(fi.fIndex)
{
   if (fi.fUrlList) {
      fUrlList = new TList;
      fUrlList->SetOwner();
      TIter nxu(fi.fUrlList);
      TUrl *u = 0;
      while ((u = (TUrl *)nxu())) {
         fUrlList->Add(new TUrl(u->GetUrl(), kTRUE));
      }
      ResetUrl();
   }
   fSize = fi.fSize;

   if (fi.fUUID)
      fUUID = new TUUID(fi.fUUID->AsString());

   if (fi.fMD5)
      fMD5 = new TMD5(*(fi.fMD5));

   // Staged and corrupted bits
   ResetBit(TFileInfo::kStaged);
   ResetBit(TFileInfo::kCorrupted);
   if (fi.TestBit(TFileInfo::kStaged)) SetBit(TFileInfo::kStaged);
   if (fi.TestBit(TFileInfo::kCorrupted)) SetBit(TFileInfo::kCorrupted);

   if (fi.fMetaDataList) {
      fMetaDataList = new TList;
      fMetaDataList->SetOwner();
      TIter nxm(fi.fMetaDataList);
      TFileInfoMeta *fim = 0;
      while ((fim = (TFileInfoMeta *)nxm())) {
         fMetaDataList->Add(new TFileInfoMeta(*fim));
      }
   }

   // By default we ignore the index
   ResetBit(TFileInfo::kSortWithIndex);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TFileInfo::~TFileInfo()
{
   SafeDelete(fMetaDataList);
   SafeDelete(fUUID);
   SafeDelete(fMD5);
   SafeDelete(fUrlList);
}

////////////////////////////////////////////////////////////////////////////////
/// Parse the input line to extract init information from 'in'; the input
/// string is tokenized on ' '; the tokens can be prefixed by the following
/// keys:
///
///  - `url:<url1>,<url2>,...`     URLs for the file; stored in the order given
///  - `sz:<size>`                 size of the file in bytes
///  - `md5:<md5_ascii>`           MD5 sum of the file in ASCII form
///  - `uuid:<uuid>`               UUID of the file
///
///  - `tree:<name>,<entries>,<first>,<last>`
///                              meta-information about a tree in the file; the
///                              should be in the form "<subdir>/tree-name";'entries' is
///                              the number of entries in the tree; 'first' and 'last'
///                              define the entry range.
///
///  - `obj:<name>,<class>,<entries>`
///                              meta-information about a generic object in the file;
///                              the should be in the form "<subdir>/obj-name"; 'class'
///                              is the object class; 'entries' is the number of occurrences
///                              for this object.
///
///  - `idx:<index>`             Index of this file if sorting with index
///
/// Multiple occurrences of 'tree:' or 'obj:' can be specified.
/// The initializations done via the input string are super-seeded by the ones by other
/// parameters in the constructor, if any.
/// If no key is given, the token is interpreted as URL(s).

void TFileInfo::ParseInput(const char *in)
{
   // Nothing to do if the string is empty
   if (!in || strlen(in) <= 0) return;

   TString sin(in), t;
   Int_t f1 = 0;
   while (sin.Tokenize(t, f1, " ")) {
      if (t.BeginsWith("sz:")) {
         // The size
         t.Replace(0, 3, "");
         if (t.IsDigit()) sscanf(t.Data(), "%lld", &fSize);
      } else if (t.BeginsWith("md5:")) {
         // The MD5
         t.Replace(0, 4, "");
         if (t.Length() >= 32) {
            fMD5 = new TMD5;
            if (fMD5->SetDigest(t) != 0)
               SafeDelete(fMD5);
         }
      } else if (t.BeginsWith("uuid:")) {
         // The UUID
         t.Replace(0, 5, "");
         if (t.Length() > 0) fUUID = new TUUID(t);
      } else if (t.BeginsWith("tree:")) {
         // A tree
         t.Replace(0, 5, "");
         TString nm, se, sf, sl;
         Long64_t ent = -1, fst= -1, lst = -1;
         Int_t f2 = 0;
         if (t.Tokenize(nm, f2, ","))
            if (t.Tokenize(se, f2, ","))
               if (t.Tokenize(sf, f2, ","))
                  t.Tokenize(sl, f2, ",");
         if (!(nm.IsNull())) {
            if (se.IsDigit()) sscanf(se.Data(), "%lld", &ent);
            if (sf.IsDigit()) sscanf(sf.Data(), "%lld", &fst);
            if (sl.IsDigit()) sscanf(sl.Data(), "%lld", &lst);
            TFileInfoMeta *meta = new TFileInfoMeta(nm, "TTree", ent, fst, lst);
            RemoveMetaData(meta->GetName());
            AddMetaData(meta);
         }
      } else if (t.BeginsWith("obj:")) {
         // A generic object
         t.Replace(0, 4, "");
         TString nm, cl, se;
         Long64_t ent = -1;
         Int_t f2 = 0;
         if (t.Tokenize(nm, f2, ","))
            if (t.Tokenize(cl, f2, ","))
               t.Tokenize(se, f2, ",");
         if (cl.IsNull()) cl = "TObject";
         if (!(nm.IsNull())) {
            if (se.IsDigit()) sscanf(se.Data(), "%lld", &ent);
            TFileInfoMeta *meta = new TFileInfoMeta(nm, cl, ent);
            AddMetaData(meta);
         }
      } else if (t.BeginsWith("idx:")) {
         // The size
         t.Replace(0, 4, "");
         if (t.IsDigit()) sscanf(t.Data(), "%d", &fIndex);
      } else {
         // A (set of) URL(s)
         if (t.BeginsWith("url:")) t.Replace(0, 4, "");
         TString u;
         Int_t f2 = 0;
         while (t.Tokenize(u, f2, ",")) {
            if (!(u.IsNull())) AddUrl(u);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the UUID to the value associated to the string 'uuid'. This is
/// useful to set the UUID to the one of the ROOT file during verification.
///
/// NB: we do not change the name in here, because this would screw up lists
///     of these objects hashed on the name. Those lists need to be rebuild.
///     TFileCollection does that in RemoveDuplicates.

void TFileInfo::SetUUID(const char *uuid)
{
   if (uuid) {
      if (fUUID) delete fUUID;
      fUUID = new TUUID(uuid);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current url.

TUrl *TFileInfo::GetCurrentUrl() const
{
   if (!fCurrentUrl)
      const_cast<TFileInfo*>(this)->ResetUrl();
   return fCurrentUrl;
}

////////////////////////////////////////////////////////////////////////////////
/// Iterator function, start iteration by calling ResetUrl().
/// The first call to NextUrl() will return the 1st element,
/// the seconde the 2nd element etc. Returns 0 in case no more urls.

TUrl *TFileInfo::NextUrl()
{
   if (!fUrlList)
      return 0;

   TUrl *returl = fCurrentUrl;

   if (fCurrentUrl)
      fCurrentUrl = (TUrl*)fUrlList->After(fCurrentUrl);

   return returl;
}

////////////////////////////////////////////////////////////////////////////////
/// Find an element from a URL. Returns 0 if not found.

TUrl *TFileInfo::FindByUrl(const char *url, Bool_t withDeflt)
{
   TIter nextUrl(fUrlList);
   TUrl *urlelement;

   TRegexp rg(url);
   while  ((urlelement = (TUrl*) nextUrl())) {
      if (TString(urlelement->GetUrl(withDeflt)).Index(rg) != kNPOS) {
         return urlelement;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new URL. If 'infront' is TRUE the new url is pushed at the beginning
/// of the list; otherwise is pushed back.
/// Returns kTRUE if successful, kFALSE otherwise.

Bool_t TFileInfo::AddUrl(const char *url, Bool_t infront)
{
   if (FindByUrl(url))
      return kFALSE;

   if (!fUrlList) {
      fUrlList = new TList;
      fUrlList->SetOwner();
   }

   TUrl *newurl = new TUrl(url, kTRUE);
   // We set the current Url to the first url added
   if (fUrlList->GetSize() == 0)
      fCurrentUrl = newurl;

   if (infront)
      fUrlList->AddFirst(newurl);
   else
      fUrlList->Add(newurl);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove an URL. Returns kTRUE if successful, kFALSE otherwise.

Bool_t TFileInfo::RemoveUrl(const char *url)
{
   TUrl *lurl;
   if ((lurl = FindByUrl(url))) {
      fUrlList->Remove(lurl);
      if (lurl == fCurrentUrl)
         ResetUrl();
      delete lurl;
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove URL at given position. Returns kTRUE on success, kFALSE on error.

Bool_t TFileInfo::RemoveUrlAt(Int_t i)
{
   TUrl *tUrl;
   if ((tUrl = dynamic_cast<TUrl *>(fUrlList->At(i))) != NULL) {
      fUrlList->Remove(tUrl);
      if (tUrl == fCurrentUrl)
         ResetUrl();
      delete tUrl;
      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set 'url' as current URL, if in the list
/// Return kFALSE if not in the list

Bool_t TFileInfo::SetCurrentUrl(const char *url)
{
   TUrl *lurl;
   if ((lurl = FindByUrl(url))) {
      fCurrentUrl = lurl;
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set 'url' as current URL, if in the list
/// Return kFALSE if not in the list

Bool_t TFileInfo::SetCurrentUrl(TUrl *url)
{
   if (url && fUrlList && fUrlList->FindObject(url)) {
      fCurrentUrl = url;
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add's a meta data object to the file info object. The object will be
/// adopted by the TFileInfo and should not be deleted by the user.
/// Typically objects of class TFileInfoMeta or derivatives should be added,
/// but any class is accepted.
/// Returns kTRUE if successful, kFALSE otherwise.

Bool_t TFileInfo::AddMetaData(TObject *meta)
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

////////////////////////////////////////////////////////////////////////////////
/// Remove the metadata object. If meta is 0 remove all meta data objects.
/// Returns kTRUE if successful, kFALSE otherwise.

Bool_t TFileInfo::RemoveMetaData(const char *meta)
{
   if (fMetaDataList) {
      if (!meta || strlen(meta) <= 0) {
         SafeDelete(fMetaDataList);
         return kTRUE;
      } else {
         TObject *o = fMetaDataList->FindObject(meta);
         if (o) {
            fMetaDataList->Remove(o);
            delete o;
            return kTRUE;
         }
      }
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get meta data object with specified name. If meta is 0
/// get first meta data object. Returns 0 in case no
/// suitable meta data object is found.

TFileInfoMeta *TFileInfo::GetMetaData(const char *meta) const
{
   if (fMetaDataList) {
      TFileInfoMeta *m;
      if (!meta || strlen(meta) <= 0)
         m = (TFileInfoMeta *) fMetaDataList->First();
      else
         m = (TFileInfoMeta *) fMetaDataList->FindObject(meta);
      if (m) {
         TClass *c = m->IsA();
         return (c && c->InheritsFrom(TFileInfoMeta::Class())) ? m : 0;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Compare TFileInfo object by their first urls.

Int_t TFileInfo::Compare(const TObject *obj) const
{
   Int_t rc = 0;
   if (TestBit(TFileInfo::kSortWithIndex)) {
      const TFileInfo *fi = dynamic_cast<const TFileInfo *>(obj);
      if (!fi) {
         rc = -1;
      } else {
         if (fIndex < fi->fIndex) {
            rc = -1;
         } else if (fIndex > fi->fIndex) {
            rc = 1;
         }
      }
   } else {
      if (this == obj) {
         rc = 0;
      } else if (TFileInfo::Class() != obj->IsA()) {
         rc = -1;
      } else {
         rc = (GetFirstUrl()->Compare(((TFileInfo*)obj)->GetFirstUrl()));
      }
   }
   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Print information about this object. If option contains 'L' a long listing
/// will be printed (on multiple lines). Otherwise one line is printed with the
/// following information: current url, default tree name|class|entries, md5;
/// the default tree name is passed via the option ("T:<default_tree>") by the
/// owning TFileCollection.

void TFileInfo::Print(Option_t *option) const
{
   if (GetMD5()) GetMD5()->Final();
   TString opt(option);
   if (opt.Contains("L", TString::kIgnoreCase)) {

      Printf("UUID: %s\nMD5:  %s\nSize: %lld\nIndex: %d",
             GetUUID() ? GetUUID()->AsString() : "undef",
             GetMD5() ? GetMD5()->AsString() : "undef",
             GetSize(), GetIndex());

      TIter next(fUrlList);
      TUrl *u;
      Printf(" === URLs ===");
      while ((u = (TUrl*)next()))
         Printf(" URL:  %s", u->GetUrl());

      TIter nextm(fMetaDataList);
      TObject *m = 0;   // can be any TObject not only TFileInfoMeta
      while ((m = (TObject*) nextm())) {
         Printf(" === Meta Data Object ===");
         m->Print();
      }
   } else {
      TString out("current-url-undef -|-|- md5-undef");
      if (GetCurrentUrl()) out.ReplaceAll("current-url-undef", GetCurrentUrl()->GetUrl());
      // Extract the default tree name, if any
      TString deft;
      if (opt.Contains("T:")) deft = opt(opt.Index("T:")+2, opt.Length());
      TFileInfoMeta *meta = 0;
      if (fMetaDataList && !deft.IsNull()) meta = (TFileInfoMeta *) fMetaDataList->FindObject(deft);
      if (fMetaDataList && !meta) meta = (TFileInfoMeta *) fMetaDataList->First();
      if (meta) out.ReplaceAll("-|-|-", TString::Format("%s|%s|%lld", meta->GetName(),
                                        meta->GetTitle(), meta->GetEntries()));
      if (GetMD5())
         out.ReplaceAll("md5-undef", TString::Format("%s", GetMD5()->AsString()));
      Printf("%s", out.Data());
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Create file meta data object.

TFileInfoMeta::TFileInfoMeta(const char *objPath, const char *objClass,
                             Long64_t entries, Long64_t first, Long64_t last,
                             Long64_t totbytes, Long64_t zipbytes)
              : TNamed(objPath, objClass), fEntries(entries), fFirst(first),
                fLast(last), fTotBytes(totbytes), fZipBytes(zipbytes)
{
   TString p = objPath;
   if (!p.BeginsWith("/")) {
      p.Prepend("/");
      SetName(p);
   }

   TClass *c = TClass::GetClass(objClass);
   fIsTree = (c && c->InheritsFrom("TTree")) ? kTRUE : kFALSE;
   ResetBit(TFileInfoMeta::kExternal);
}

////////////////////////////////////////////////////////////////////////////////
/// Create file meta data object.

TFileInfoMeta::TFileInfoMeta(const char *objPath, const char *objDir,
                             const char *objClass, Long64_t entries,
                             Long64_t first, Long64_t last,
                             Long64_t totbytes, Long64_t zipbytes)
              : TNamed(objPath, objClass), fEntries(entries), fFirst(first),
                fLast(last), fTotBytes(totbytes), fZipBytes(zipbytes)
{
   TString sdir = objDir;
   if (!sdir.BeginsWith("/"))
      sdir.Prepend("/");
   if (!sdir.EndsWith("/"))
      sdir += "/";
   sdir += objPath;
   SetName(sdir);

   TClass *c = TClass::GetClass(objClass);
   fIsTree = (c && c->InheritsFrom("TTree")) ? kTRUE : kFALSE;
   ResetBit(TFileInfoMeta::kExternal);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

TFileInfoMeta::TFileInfoMeta(const TFileInfoMeta &m)
              : TNamed(m.GetName(), m.GetTitle())
{
   fEntries = m.fEntries;
   fFirst = m.fFirst;
   fLast = m.fLast;
   fIsTree = m.fIsTree;
   fTotBytes = m.fTotBytes;
   fZipBytes = m.fZipBytes;
   ResetBit(TFileInfoMeta::kExternal);
   if (m.TestBit(TFileInfoMeta::kExternal)) SetBit(TFileInfoMeta::kExternal);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the object's directory in the ROOT file.

const char *TFileInfoMeta::GetDirectory() const
{
   return gSystem->DirName(GetName());
}

////////////////////////////////////////////////////////////////////////////////
/// Get the object name, with path stripped off. For full path
/// use GetName().

const char *TFileInfoMeta::GetObject() const
{
   return gSystem->BaseName(GetName());
}

////////////////////////////////////////////////////////////////////////////////
/// Print information about this object.

void TFileInfoMeta::Print(Option_t * /* option */) const
{
   Printf(" Name:    %s\n Class:   %s\n Entries: %lld\n"
          " First:   %lld\n Last:    %lld",
          fName.Data(), fTitle.Data(), fEntries, fFirst, fLast);
}
