// @(#)root/proof:$Id$
// Author: Fons Rademakers   11/01/02

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TDSetElement
\ingroup proofkernel
Manages an element of a TDSet.

See TDSet.
*/

#include "TDSet.h"

#include "Riostream.h"
#include "TBuffer.h"
#include "TChain.h"
#include "TClass.h"
#include "TClassTable.h"
#include "TCut.h"
#include "TDataSetManager.h"
#include "TError.h"
#include "TEntryList.h"
#include "TEnv.h"
#include "TEventList.h"
#include "TFile.h"
#include "TFileInfo.h"
#include "TFileStager.h"
#include "TFriendElement.h"
#include "TKey.h"
#include "THashList.h"
#include "TMap.h"
#include "TROOT.h"
#include "TTimeStamp.h"
#include "TTree.h"
#include "TUrl.h"
#include "TRegexp.h"
#include "TVirtualPerfStats.h"
#include "TProof.h"
#include "TProofChain.h"
#include "TProofServ.h"
#include "TPluginManager.h"
#include "TChainElement.h"
#include "TSystem.h"
#include "TSelector.h"
#include "TObjString.h"

#include "TVirtualStreamerInfo.h"
#include "TClassRef.h"

ClassImp(TDSetElement);
ClassImp(TDSet);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TDSetElement::TDSetElement() : TNamed("",""),
                               fDirectory(), fFirst(0), fNum(0), fMsd(),
                               fTDSetOffset(0), fEntryList(0), fValid(kFALSE),
                               fEntries(0), fFriends(0), fDataSet(), fAssocObjList(0),
                               fMaxProcTime(-1)
{
   ResetBit(kWriteV3);
   ResetBit(kHasBeenLookedUp);
   ResetBit(kEmpty);
   ResetBit(kCorrupted);
   ResetBit(kNewRun);
   ResetBit(kNewPacket);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TDSet element.

TDSetElement::TDSetElement(const char *file, const char *objname, const char *dir,
                           Long64_t first, Long64_t num,
                           const char *msd, const char *dataset)
             : TNamed(file, objname)
{
   if (first < 0) {
      Warning("TDSetElement", "first must be >= 0, %lld is not allowed - setting to 0", first);
      fFirst = 0;
   } else {
      fFirst = first;
   }
   if (num < -1) {
      Warning("TDSetElement", "num must be >= -1, %lld is not allowed - setting to -1", num);
      fNum   = -1;
   } else {
      fNum   = num;
   }
   fMsd         = msd;
   fTDSetOffset = 0;
   fEntryList   = 0;
   fFriends     = 0;
   fValid       = kFALSE;
   fEntries     = -1;
   fDataSet     = dataset;
   fAssocObjList = 0;
   if (dir)
      fDirectory = dir;
   fMaxProcTime = -1.;

   ResetBit(kWriteV3);
   ResetBit(kHasBeenLookedUp);
   ResetBit(kEmpty);
   ResetBit(kCorrupted);
   ResetBit(kNewRun);
   ResetBit(kNewPacket);
}

////////////////////////////////////////////////////////////////////////////////
/// copy constructor

TDSetElement::TDSetElement(const TDSetElement& elem)
             : TNamed(elem.GetFileName(), elem.GetObjName())
{
   fDirectory = elem.GetDirectory();
   fFirst = elem.fFirst;
   fNum = elem.fNum;
   fMsd = elem.fMsd;
   fTDSetOffset = elem.fTDSetOffset;
   fEntryList = 0;
   fValid = elem.fValid;
   fEntries = elem.fEntries;
   fFriends = 0;
   fDataSet = elem.fDataSet;
   fAssocObjList = 0;
   fMaxProcTime = elem.fMaxProcTime;
   ResetBit(kWriteV3);
   ResetBit(kHasBeenLookedUp);
   ResetBit(kEmpty);
   ResetBit(kCorrupted);
   ResetBit(kNewRun);
   ResetBit(kNewPacket);
}

////////////////////////////////////////////////////////////////////////////////
/// Clean up the element.

TDSetElement::~TDSetElement()
{
   DeleteFriends();
   if (fAssocObjList) {
      fAssocObjList->SetOwner(kTRUE);
      SafeDelete(fAssocObjList);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset TDSet element.

void TDSetElement::Reset()
{
   fFirst = 0;
   fNum   = -1;
   fTDSetOffset = 0;
   fEntryList   = 0;
   fValid       = kFALSE;
   fEntries     = -1;
   fMaxProcTime = -1.;

   ResetBit(kHasBeenLookedUp);
   ResetBit(kEmpty);
   ResetBit(kNewRun);
   ResetBit(kNewPacket);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if 'elem' is overlapping or subsequent and, if the case, return
/// a merged element.
/// Returns:
///     1    if the elements are overlapping
///     0    if the elements are subsequent
///    -1    if the elements are neither overlapping nor subsequent

Int_t TDSetElement::MergeElement(TDSetElement *elem)
{
   // The element must be defined
   if (!elem) return -1;

   // The file names and object names must be the same
   if (strcmp(GetName(), elem->GetName()) || strcmp(GetTitle(), elem->GetTitle()))
      return -1;

   Int_t rc = -1;
   // Check the overlap or subsequency
   if (fFirst == 0 && fNum == -1) {
      // Overlap, since we cover already the full range
      rc = 1;
   } else if (elem->GetFirst() == 0 && elem->GetNum() == -1) {
      // Overlap, since 'elem' cover already the full range
      fFirst = 0;
      fNum = -1;
      fEntries = elem->GetEntries();
      rc = 1;
   } else if (fFirst >= 0 && fNum > 0 && elem->GetFirst() >= 0 && elem->GetNum() > 0) {
      Long64_t last = fFirst + fNum - 1, lastref = 0;
      Long64_t lastelem = elem->GetFirst() + elem->GetNum() - 1;
      if (elem->GetFirst() == last + 1) {
         lastref = lastelem;
         rc = 0;
      } else if (fFirst == lastelem + 1) {
         fFirst += elem->GetFirst();
         lastref = last;
         rc = 0;
      } else if (elem->GetFirst() < last + 1 && elem->GetFirst() >= fFirst) {
         lastref = (lastelem > last) ? lastelem : last;
         rc = 1;
      } else if (fFirst < lastelem + 1 && fFirst >= elem->GetFirst()) {
         fFirst += elem->GetFirst();
         lastref = (lastelem > last) ? lastelem : last;
         rc = 1;
      }
      fNum = lastref - fFirst + 1;
   }
   if (rc >= 0 && fEntries < 0 && elem->GetEntries() > 0) fEntries = elem->GetEntries();

   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the content of this element in the form of a TFileInfo

TFileInfo *TDSetElement::GetFileInfo(const char *type)
{
   // Create the TFileInfoMeta object
   TFileInfoMeta *meta = 0;
   Long64_t entries = (fEntries < 0 && fNum > 0) ? fNum : fEntries;
   Printf("entries: %lld (%lld)", entries, fNum);
   if (!strcmp(type, "TTree")) {
      meta = new TFileInfoMeta(GetTitle(), "TTree", entries, fFirst,
                                fFirst + entries - 1);
   } else {
      meta = new TFileInfoMeta(GetTitle(), fDirectory, type, entries, fFirst,
                                fFirst + entries - 1);
   }
   TFileInfo *fi = new TFileInfo(GetName(), 0, 0, 0, meta);
   if (!fDataSet.IsNull()) fi->SetTitle(fDataSet.Data());
   if (TestBit(TDSetElement::kCorrupted)) fi->SetBit(TFileInfo::kCorrupted);
   return fi;
}

////////////////////////////////////////////////////////////////////////////////
/// Return directory where to look for object.

const char *TDSetElement::GetDirectory() const
{
   return fDirectory;
}

////////////////////////////////////////////////////////////////////////////////
/// Print a TDSetElement. When option="a" print full data.

void TDSetElement::Print(Option_t *opt) const
{
   if (opt && opt[0] == 'a') {
      Printf("%s file=\"%s\" dir=\"%s\" obj=\"%s\" first=%lld num=%lld msd=\"%s\"",
             IsA()->GetName(), GetName(), fDirectory.Data(), GetTitle(),
             fFirst, fNum, fMsd.Data());
   } else {
      Printf("\tLFN: %s", GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Validate by opening the file.

void TDSetElement::Validate(Bool_t isTree)
{
   Long64_t entries = GetEntries(isTree);
   if (entries < 0) return; // Error should be reported by GetEntries()
   if (fFirst < entries) {
      if (fNum == -1) {
         fNum = entries - fFirst;
         fValid = kTRUE;
      } else {
         if (fNum <= entries - fFirst) {
            fValid = kTRUE;
         } else {
            Error("Validate", "TDSetElement has only %lld entries starting"
                              " with entry %lld, while %lld were requested",
                              entries - fFirst, fFirst, fNum);
         }
      }
   } else {
      Error("Validate", "TDSetElement has only %lld entries with"
            " first entry requested as %lld", entries, fFirst);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Validate by checking against another element.

void TDSetElement::Validate(TDSetElement *elem)
{
   // NOTE: Since this function validates against another TDSetElement,
   //       if the other TDSetElement (elem) did not use -1 to request all
   //       entries, this TDSetElement may get less than all entries if it
   //       requests all (with -1). For the application it was developed for
   //       (TProofSuperMaster::ValidateDSet) it works, since the design was
   //       to send the elements to their mass storage domain and let them
   //       look at the file and send the info back to the supermaster. The
   //       ability to set fValid was also required to be only exist in
   //       TDSetElement through certain function and not be set externally.
   //       TDSetElement may need to be extended for more general applications.

   if (!elem || !elem->GetValid()) {
      Error("Validate", "TDSetElement to validate against is not valid");
      return;
   }

   TString name = TUrl(GetFileName()).GetFileAndOptions();
   TString elemname = TUrl(elem->GetFileName()).GetFileAndOptions();
   if ((name == elemname) &&
       !strcmp(GetDirectory(), elem->GetDirectory()) &&
       !strcmp(GetObjName(), elem->GetObjName())) {
      Long64_t entries = elem->fFirst + elem->fNum;
      if (fFirst < entries) {
         if (fNum == -1) {
            fNum = entries - fFirst;
            fValid = kTRUE;
         } else {
            if (fNum <= entries - fFirst) {
               fValid = kTRUE;
            } else {
               Error("Validate", "TDSetElement requests %lld entries starting"
                                 " with entry %lld, while TDSetElement to validate against"
                                 " has only %lld entries", fNum, fFirst, entries);
            }
         }
      } else {
         Error("Validate", "TDSetElement to validate against has only %lld"
               " entries, but this TDSetElement requested %lld as its first"
               " entry", entries, fFirst);
      }
   } else {
      Error("Validate", "TDSetElements do not refer to same objects");
   }
}

////////////////////////////////////////////////////////////////////////////////
///Compare elements by filename (and the fFirst).

Int_t TDSetElement::Compare(const TObject *obj) const
{
   if (this == obj) return 0;

   const TDSetElement *elem = dynamic_cast<const TDSetElement*>(obj);
   if (!elem) {
      if (obj)
         return (strncmp(GetName(),obj->GetName(),strlen(GetName()))) ? 1 : 0;
      return -1;
   }

   Int_t order = strncmp(GetName(),elem->GetFileName(),strlen(GetName()));
   if (order == 0) {
      if (GetFirst() < elem->GetFirst())
         return -1;
      else if (GetFirst() > elem->GetFirst())
         return 1;
      return 0;
   }
   return order;
}

////////////////////////////////////////////////////////////////////////////////
/// Add friend TDSetElement to this set. The friend element will be copied to this object.

void TDSetElement::AddFriend(TDSetElement *friendElement, const char *alias)
{
   if (!friendElement) {
      Error("AddFriend", "The friend TDSetElement is null!");
      return;
   }
   if (!fFriends) {
      fFriends = new TList();
      fFriends->SetOwner();
   }
   // Add alias (if any) as option 'friend_alias=<alias>|'
   if (alias && strlen(alias) > 0) {
      TUrl uf(friendElement->GetName());
      TString uo(uf.GetOptions());
      uo += TString::Format("friend_alias=%s|", alias);
      uf.SetOptions(uo);
      friendElement->SetName(uf.GetUrl());
   }
   fFriends->Add(new TDSetElement(*friendElement));
}

////////////////////////////////////////////////////////////////////////////////
/// Deletes the list of friends and all the friends on the list.

void TDSetElement::DeleteFriends()
{
   if (!fFriends)
      return;

   fFriends->SetOwner(kTRUE);
   delete fFriends;
   fFriends = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns next TDSetElement.

TDSetElement *TDSet::Next(Long64_t /*totalEntries*/)
{
   if (!fIterator) {
      fIterator = new TIter(fElements);
   }

   fCurrent = (TDSetElement *) fIterator->Next();
   return fCurrent;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns number of entries in tree or objects in file.
/// If not yet defined and 'openfile' is TRUE, get the number from the file
/// (may considerably slow down the application).
/// Returns -1 in case of error.

Long64_t TDSetElement::GetEntries(Bool_t isTree, Bool_t openfile)
{
   if (fEntries > -1 || !openfile)
      return fEntries;

   Double_t start = 0;
   if (gPerfStats) start = TTimeStamp();

   // Take into account possible prefixes
   TFile::EFileType typ = TFile::kDefault;
   TString fname = gEnv->GetValue("Path.Localroot",""), pfx(fname);
   // Get the locality (disable warnings or errors in connection attempts)
   Int_t oldLevel = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kError+1;
   if ((typ = TFile::GetType(GetName(), "", &fname)) != TFile::kLocal) fname = GetName();
   gErrorIgnoreLevel = oldLevel;
   // Open the file
   TFile *file = TFile::Open(fname);

   if (gPerfStats)
      gPerfStats->FileOpenEvent(file, GetName(), start);

   if (file == 0) {
      ::SysError("TDSetElement::GetEntries",
                 "cannot open file %s (type: %d, pfx: %s)", GetName(), typ, pfx.Data());
      return -1;
   }

   // Record end-point Url and mark as looked-up; be careful to change
   // nothing in the file name, otherwise some cross-checks may fail.
   // The lookup is only performed if not yet done
   if (Lookup(kFALSE) != 0) Warning("GetEntries", "lookup problems for %s", GetName());

   TDirectory *dirsave = gDirectory;
   if (!file->cd(fDirectory)) {
      Error("GetEntries", "cannot cd to %s", fDirectory.Data());
      delete file;
      return -1;
   }

   TDirectory *dir = gDirectory;
   dirsave->cd();

   if (isTree) {

      TString on(GetTitle());
      TString sreg(GetTitle());
      // If a wild card we will use the first object of the type
      // requested compatible with the reg expression we got
      if (sreg.Length() <= 0 || sreg == "" || sreg.Contains("*")) {
         if (sreg.Contains("*"))
            sreg.ReplaceAll("*", ".*");
         else
            sreg = ".*";
         TRegexp re(sreg);
         if (dir->GetListOfKeys()) {
            TIter nxk(dir->GetListOfKeys());
            TKey *k = 0;
            Bool_t notfound = kTRUE;
            while ((k = (TKey *) nxk())) {
               if (!strcmp(k->GetClassName(), "TTree")) {
                  TString kn(k->GetName());
                  if (kn.Index(re) != kNPOS) {
                     if (notfound) {
                        on = kn;
                        notfound = kFALSE;
                     } else if (kn != on) {
                       Warning("GetEntries",
                               "additional tree found in the file: %s", kn.Data());
                     }
                  }
               }
            }
         }
      }

      TKey *key = dir->GetKey(on);
      if (key == 0) {
         Error("GetEntries", "cannot find tree \"%s\" in %s",
               GetTitle(), GetName());
         delete file;
         return -1;
      }
      TTree *tree = (TTree *) key->ReadObj();
      if (tree == 0) {
         // Error always reported?
         delete file;
         return -1;
      }
      fEntries = tree->GetEntries();
      delete tree;

   } else {
      TList *keys = dir->GetListOfKeys();
      fEntries = keys->GetSize();
   }

   delete file;
   return fEntries;
}

////////////////////////////////////////////////////////////////////////////////
/// Resolve end-point URL for this element
/// Return 0 on success and -1 otherwise

Int_t TDSetElement::Lookup(Bool_t force)
{
   static Int_t xNetPluginOK = -1;
   static TFileStager *xStager = 0;
   Int_t retVal = 0;

   // Check if required
   if (!force && HasBeenLookedUp())
      return retVal;

   TUrl url(GetName());
   // Save current options and anchor to be set ofthe final end URL
   TString anch = url.GetAnchor();
   TString opts = url.GetOptions();
   // The full path
   TString name(url.GetUrl());

   // Depending on the type of backend, it might not make any sense to lookup
   Bool_t doit = kFALSE;
   TFile::EFileType type = TFile::GetType(name, "");
   if (type == TFile::kNet) {
      TPluginHandler *h = 0;
      // Network files via XROOTD
      if (xNetPluginOK == -1) {
         // Check the plugin the first time
         xNetPluginOK = 0;
         if ((h = gROOT->GetPluginManager()->FindHandler("TFile", name)) &&
            !strcmp(h->GetClass(),"TXNetFile") && h->LoadPlugin() == 0)
            xNetPluginOK = 1;
      }
      doit = (xNetPluginOK == 1) ? kTRUE : kFALSE;
   }

   // Locate the file
   if (doit) {
      if (!xStager || !xStager->Matches(name)) {
         SafeDelete(xStager);
         if (!(xStager = TFileStager::Open(name))) {
            Error("Lookup", "TFileStager instance cannot be instantiated");
            retVal = -1;
         }
      }
      if (xStager && xStager->Locate(name.Data(), name) == 0) {
         // Get the effective end-point Url
         url.SetUrl(name);
         // Restore original options and anchor, if any
         url.SetOptions(opts);
         url.SetAnchor(anch);
         // Save it into the element
         fName = url.GetUrl();
      } else {
         // Failure
         Error("Lookup", "couldn't lookup %s", name.Data());
         retVal = -1;
      }
   }

   // Mark has looked-up
   SetBit(kHasBeenLookedUp);
   return retVal;
}

////////////////////////////////////////////////////////////////////////////////
/// Set entry (or event) list for this element

void TDSetElement::SetEntryList(TObject *aList, Long64_t first, Long64_t num)
{
   if (!aList) {
      // Nothing to do, except making sure to disable any previous setting
      fEntryList = 0;
      return;
   }

   // Link the proper object
   TEventList *evl = 0;
   TEntryList *enl = dynamic_cast<TEntryList*>(aList);
   if (!enl)
      evl = dynamic_cast<TEventList*>(aList);
   if (!enl && !evl) {
      Error("SetEntryList", "type of input object must be either TEntryList "
                            "or TEventList (found: '%s' - do nothing", aList->ClassName());
      return;
   }

   // Action depends on the type
   if (enl) {
      enl->SetEntriesToProcess(num);
   } else {
      for (; num > 0; num--, first++)
         evl->Enter(evl->GetEntry((Int_t)first));
   }
   fEntryList = aList;

   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Add an associated object to the list

void TDSetElement::AddAssocObj(TObject *assocobj)
{
   if (assocobj) {
      if (!fAssocObjList) fAssocObjList = new TList;
      if (fAssocObjList) fAssocObjList->Add(assocobj);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get i-th associated object.
/// If 'isentry' fFirst is subtracted, so that i == fFirst returns the first
/// object in the list.
/// If there are not enough elements in the list, the element i%list_size is
/// returned (if the list has only one element this only one element is always
/// returned.
/// This method is used when packet processing consist in processing the objects
/// in the associated object list.

TObject *TDSetElement::GetAssocObj(Long64_t i, Bool_t isentry)
{
   TObject *o = 0;
   if (!fAssocObjList || fAssocObjList->GetSize() <= 0) return o;

   TString s;
   Int_t pos = -1;
   if (isentry) {
      if (i < fFirst) return o;
      s.Form("%lld", i - fFirst);
   } else {
      if (i < 0) return o;
      s.Form("%lld", i);
   }
   if (!(s.IsDigit())) return o;
   pos = s.Atoi();
   if (pos > fAssocObjList->GetSize() - 1) pos %= fAssocObjList->GetSize();
   return fAssocObjList->At(pos);
}


/** \class TDSet
\ingroup proofkernel

This class implements a data set to be used for PROOF processing.
The TDSet defines the class of which objects will be processed,
the directory in the file where the objects of that type can be
found and the list of files to be processed. The files can be
specified as logical file names (LFN's) or as physical file names
(PFN's). In case of LFN's the resolution to PFN's will be done
according to the currently active GRID interface.
Examples:
  TDSet treeset("TTree", "AOD");
  treeset.Add("lfn:/alien.cern.ch/alice/prod2002/file1");
  ...
  treeset.AddFriend(friendset);

or

  TDSet objset("MyEvent", "*", "/events");
  objset.Add("root://cms.cern.ch/user/prod2002/hprod_1.root");
  ...
  objset.Add(set2003);

Validity of file names will only be checked at processing time
(typically on the PROOF master server), not at creation time.

*/

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TDSet::TDSet()
{
   fElements = new THashList;
   fElements->SetOwner();
   fIsTree    = kFALSE;
   fIterator  = 0;
   fCurrent   = 0;
   fEntryList = 0;
   fProofChain = 0;
   fSrvMaps = 0;
   fSrvMapsIter = 0;
   ResetBit(kWriteV3);
   ResetBit(kEmpty);
   ResetBit(kValidityChecked);
   ResetBit(kSomeInvalid);
   ResetBit(kMultiDSet);

   // Add to the global list
   gROOT->GetListOfDataSets()->Add(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a named TDSet object. The "type" defines the class of which objects
/// will be processed (default 'TTree'). The optional "objname" argument
/// specifies the name of the objects of the specified class.
/// If the "objname" is not given the behaviour depends on the 'type':
/// for 'TTree' the first TTree is analyzed; for other types, all objects of
/// the class found in the specified directory are processed.
/// The "dir" argument specifies in which directory the objects are
/// to be found, the top level directory ("/") is the default.
/// Directories can be specified using wildcards, e.g. "*" or "/*"
/// means to look in all top level directories, "/dir/*" in all
/// directories under "/dir", and "/*/*" to look in all directories
/// two levels deep.
/// For backward compatibility the type can also be passed via 'name',
/// in which case 'type' is ignored.

TDSet::TDSet(const char *name,
             const char *objname, const char *dir, const char *type)
{
   fElements = new THashList;
   fElements->SetOwner();
   fIterator = 0;
   fCurrent  = 0;
   fEntryList = 0;
   fProofChain = 0;
   fSrvMaps = 0;
   fSrvMapsIter = 0;
   ResetBit(kWriteV3);
   ResetBit(kEmpty);
   ResetBit(kValidityChecked);
   ResetBit(kSomeInvalid);
   ResetBit(kMultiDSet);

   fType = "TTree";
   TClass *c = 0;
   // Check name
   if (name && strlen(name) > 0) {
      // In the old constructor signature it was the 'type'
      if (!type) {
         TString cn(name);
         if (cn.Contains(':')) cn.Remove(0, cn.Index(":")+1);
         if (TClass::GetClass(cn))
            fType = cn;
         else
            // Default type is 'TTree'
            fName = name;
      } else {
         // Set name
         fName = name;
         // Check type
         if (strlen(type) > 0)
            if (TClass::GetClass(type))
               fType = type;
      }
   } else if (type && strlen(type) > 0) {
      // Check the type
      if (TClass::GetClass(type))
         fType = type;
   }
   // The correct class type
   c = TClass::GetClass(fType);

   fIsTree = (c->InheritsFrom(TTree::Class())) ? kTRUE : kFALSE;

   if (objname)
      fObjName = objname;

   if (dir)
      fDir = dir;

   // Default name is the object name
   if (fName.Length() <= 0)
      fName = TString::Format("TDSet:%s", fObjName.Data());
   // We set the default title to the 'type'
   fTitle = fType;

   // Add to the global list
   gROOT->GetListOfDataSets()->Add(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a named TDSet object from existing TChain 'chain'.
/// If 'withfriends' is kTRUE add also friends.
/// This constructor substituted the static methods TChain::MakeTDSet
/// removing any residual dependence of 'tree' on 'proof'.

TDSet::TDSet(const TChain &chain, Bool_t withfriends)
{
   fElements = new THashList;
   fElements->SetOwner();
   fIterator = 0;
   fCurrent  = 0;
   fEntryList = 0;
   fProofChain = 0;
   fSrvMaps = 0;
   fSrvMapsIter = 0;
   ResetBit(kWriteV3);
   ResetBit(kEmpty);
   ResetBit(kValidityChecked);
   ResetBit(kSomeInvalid);
   ResetBit(kMultiDSet);

   fType = "TTree";
   fIsTree = kTRUE;
   fObjName = chain.GetName();
   fName = TString::Format("TChain:%s", chain.GetName());

   // First fill elements without friends()
   TIter next(chain.GetListOfFiles());
   TChainElement *elem = 0;
   TString key;
   while ((elem = (TChainElement *)next())) {
      TString file(elem->GetTitle());
      TString tree(elem->GetName());
      Int_t isl = tree.Last('/');
      TString dir = "/";
      if (isl >= 0) {
         // Copy the tree name specification
         TString behindSlash = tree(isl + 1, tree.Length() - isl - 1);
         // and remove it from basename
         tree.Remove(isl);
         dir = tree;
         tree = behindSlash;
      }
      // Find MSD if any
      TString msd(TUrl(file).GetOptions());
      Int_t imsd = kNPOS;
      if ((imsd = msd.Index("msd=")) != kNPOS) {
         msd.Remove(0, imsd+4);
      } else {
         // Not an MSD option
         msd = "";
      }
      Long64_t nent = (elem->GetEntries() > 0 &&
                       elem->GetEntries() != TTree::kMaxEntries) ? elem->GetEntries() : -1;
      if (Add(file, tree, dir, 0, nent, ((msd.IsNull()) ? 0 : msd.Data()))) {
         if (elem->HasBeenLookedUp()) {
            // Save lookup information, if any
            TDSetElement *dse = (TDSetElement *) fElements->Last();
            if (dse) dse->SetLookedUp();
         }
      }
   }
   SetDirectory(0);

   // Add friends now, if requested
   if (withfriends) {
      TList processed;
      TList chainsQueue;
      chainsQueue.Add((TObject *)&chain);
      processed.Add((TObject *)&chain);
      while (chainsQueue.GetSize() > 0) {
         TChain *c = (TChain *) chainsQueue.First();
         chainsQueue.Remove(c);
         TIter friendsIter(c->GetListOfFriends());
         while(TFriendElement *fe = dynamic_cast<TFriendElement*> (friendsIter()) ) {
            if (TChain *fc = dynamic_cast<TChain*>(fe->GetTree())) {
               if (!processed.FindObject(fc)) {    // if not yet processed
                  processed.AddFirst(fc);
                  AddFriend(new TDSet((const TChain &)(*fc), kFALSE), fe->GetName());
                  chainsQueue.Add(fc);                        // for further processing
               }
            } else {
               Reset();
               Error("TDSet", "Only TChains supported. Found illegal tree %s",
                              fe->GetTree()->GetName());
               return;
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup.

TDSet::~TDSet()
{
   SafeDelete(fElements);
   SafeDelete(fIterator);
   SafeDelete(fProofChain);
   fSrvMaps = 0;
   fSrvMapsIter = 0;

   gROOT->GetListOfDataSets()->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Process TDSet on currently active PROOF session.
/// The last argument 'enl' specifies an entry- or event-list to be used as
/// event selection.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.

Long64_t TDSet::Process(TSelector *selector, Option_t *option, Long64_t nentries,
                        Long64_t first, TObject *enl)
{
   if (!IsValid() || !fElements->GetSize()) {
      Error("Process", "not a correctly initialized TDSet");
      return -1;
   }

   // Set entry list
   SetEntryList(enl);

   if (gProof)
      return gProof->Process(this, selector, option, nentries, first);

   Error("Process", "no active PROOF session");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Process TDSet on currently active PROOF session.
/// The last argument 'enl' specifies an entry- or event-list to be used as
/// event selection.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.

Long64_t TDSet::Process(const char *selector, Option_t *option, Long64_t nentries,
                        Long64_t first, TObject *enl)
{
   if (!IsValid() || !fElements->GetSize()) {
      Error("Process", "not a correctly initialized TDSet");
      return -1;
   }

   // Set entry list
   SetEntryList(enl);

   if (gProof)
      return gProof->Process(this, selector, option, nentries, first);

   Error("Process", "no active PROOF session");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Add objects that might be needed during the processing of
/// the selector (see Process()).

void TDSet::AddInput(TObject *obj)
{
   if (gProof) {
      gProof->AddInput(obj);
   } else {
      Error("AddInput","No PROOF session active");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Clear input object list.

void TDSet::ClearInput()
{
   if (gProof)
      gProof->ClearInput();
}

////////////////////////////////////////////////////////////////////////////////
/// Get specified object that has been produced during the processing
/// (see Process()).

TObject *TDSet::GetOutput(const char *name)
{
   if (gProof)
      return gProof->GetOutput(name);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get list with all object created during processing (see Process()).

TList *TDSet::GetOutputList()
{
   if (gProof)
      return gProof->GetOutputList();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Print TDSet basic or full data. When option="a" print full data.

void TDSet::Print(const Option_t *opt) const
{
   const char *clnm = (IsA()) ? IsA()->GetName() : "TDSet";
   Printf("OBJ: %s\ttype %s\t%s\tin %s\telements %d", clnm, GetName(),
          fObjName.Data(), GetTitle(), GetListOfElements()->GetSize());

   if (opt && opt[0] == 'a') {
      TIter next(GetListOfElements());
      TObject *obj;
      while ((obj = next())) {
         obj->Print(opt);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set/change object name.

void TDSet::SetObjName(const char *objname)
{
   if (objname) {
      fObjName = objname;
      TIter next(GetListOfElements());
      TDSetElement *e;
      while ((e = (TDSetElement *) next())) {
         e->SetTitle(objname);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set/change directory.

void TDSet::SetDirectory(const char *dir)
{
   if (dir)
      fDir = dir;
}

////////////////////////////////////////////////////////////////////////////////
/// Add file to list of files to be analyzed. Optionally with the
/// objname and dir arguments the default, TDSet wide, objname and
/// dir can be overridden.

Bool_t TDSet::Add(const char *file, const char *objname, const char *dir,
                  Long64_t first, Long64_t num, const char *msd)
{
   if (!file || !*file) {
      Error("Add", "file name must be specified");
      return kFALSE;
   }

   TString fn = file;
   if (gProof && gProof->IsLite()) {
      TUrl u(file, kTRUE);
      if (!strcmp(u.GetProtocol(), "file")) {
         fn = u.GetFileAndOptions();
         gSystem->ExpandPathName(fn);
         if (!gSystem->IsAbsoluteFileName(fn))
            gSystem->PrependPathName(gSystem->WorkingDirectory(), fn);
      }
   }

   // check, if it already exists in the TDSet
   TDSetElement *el = (TDSetElement *) fElements->FindObject(fn);
   if (!el) {
      if (!objname)
         objname = GetObjName();
      if (!dir)
         dir = GetDirectory();
      fElements->Add(new TDSetElement(fn, objname, dir, first, num, msd));
   } else {
      TString msg;
      msg.Form("duplication detected: %40s is already in dataset - ignored", fn.Data());
      Warning("Add", "%s", msg.Data());
      if (gProofServ) {
         msg.Insert(0, "WARNING: ");
         gProofServ->SendAsynMessage(msg);
      }
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add specified data set to the this set.

Bool_t TDSet::Add(TDSet *dset)
{
   if (!dset)
      return kFALSE;

   if (TestBit(TDSet::kMultiDSet)) {
      fElements->Add(dset);
      return kTRUE;
   }

   if (fType != dset->GetType()) {
      Error("Add", "cannot add a set with a different type");
      return kFALSE;
   }

   TDSetElement *el;
   TIter next(dset->fElements);
   TObject *last = (dset == this) ? fElements->Last() : 0;
   while ((el = (TDSetElement*) next())) {
      Add(el->GetFileName(), el->GetObjName(), el->GetDirectory(),
          el->GetFirst(), el->GetNum(), el->GetMsd());
      if (el == last) break;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add files passed as list of TFileInfo, TUrl or TObjString objects .
/// If TFileInfo, the first entry and the number of entries are also filled.
/// The argument 'meta' can be used to specify one of the subsets in the
/// file as described in the metadata of TFileInfo. By default the first one
/// is taken.
/// If 'availableOnly' is true only files available ('staged' and non corrupted)
/// are taken: those not satisfying this requirement are added to 'badlist', if
/// the latter is defined. By default availableOnly is false.

Bool_t TDSet::Add(TCollection *filelist, const char *meta, Bool_t availableOnly,
                  TCollection *badlist)
{
   if (!filelist)
      return kFALSE;

   TObject *o = 0;
   TIter next(filelist);
   while ((o = next())) {
      TString cn(o->ClassName());
      if (cn == "TFileInfo") {
         TFileInfo *fi = (TFileInfo *)o;
         if (!availableOnly ||
            (fi->TestBit(TFileInfo::kStaged) &&
            !fi->TestBit(TFileInfo::kCorrupted))) {
            Int_t nf = fElements->GetSize();
            if (!Add(fi, meta)) return kFALSE;
            // Duplications count as bad files
            if (fElements->GetSize() <= nf && badlist) badlist->Add(fi);
         } else if (badlist) {
            // Return list of non-usable files
            badlist->Add(fi);
         }
      } else if (cn == "TUrl") {
         Add(((TUrl *)o)->GetUrl());
      } else if (cn == "TObjString") {
         Add(((TObjString *)o)->GetName());
      } else {
         Warning("Add","found object fo unexpected type %s - ignoring", cn.Data());
      }
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set (or unset) the list for mapping servers coordinate for files.
/// Reinitialize the related iterator if needed.
/// Used by TProof.

void TDSet::SetSrvMaps(TList *srvmaps)
{
   fSrvMaps = srvmaps;
   SafeDelete(fSrvMapsIter);
   if (fSrvMaps) fSrvMapsIter = new TIter(fSrvMaps);
}

////////////////////////////////////////////////////////////////////////////////
/// Add file described by 'fi' to list of files to be analyzed.
/// The argument 'meta' can be used to specify a subsets in the
/// file as described in the metadata of TFileInfo. By default the first one
/// is taken.

Bool_t TDSet::Add(TFileInfo *fi, const char *meta)
{
   if (!fi) {
      Error("Add", "TFileInfo object name must be specified");
      return kFALSE;
   }
   TString msg;

   // Check if a remap of the server coordinates is requested
   const char *file = fi->GetFirstUrl()->GetUrl();
   Bool_t setLookedUp = kTRUE;
   TString file1;
   if (TDataSetManager::CheckDataSetSrvMaps(fi->GetFirstUrl(), file1, fSrvMaps) &&
       !(file1.IsNull())) {
      file = file1.Data();
      setLookedUp = kFALSE;
   }
   // Check if it already exists in the TDSet
   if (fElements->FindObject(file)) {
      msg.Form("duplication detected: %40s is already in dataset - ignored", file);
      Warning("Add", "%s", msg.Data());
      if (gProofServ) {
         msg.Insert(0, "WARNING: ");
         gProofServ->SendAsynMessage(msg);
      }
      return kTRUE;
   }

   // If more than one metadata info require the specification of the objpath;
   // the order in which they appear is not guaranteed and the error may be
   // very difficult to find.
   TFileInfoMeta *m = 0;
   if (!meta || strlen(meta) <= 0 || !strcmp(meta, "/")) {
      TList *fil = 0;
      if ((fil = fi->GetMetaDataList()) && fil->GetSize() > 1) {
         msg.Form("\n  Object name unspecified and several objects available.\n");
         msg += "  Please choose one from the list below:\n";
         TIter nx(fil);
         while ((m = (TFileInfoMeta *) nx())) {
            TString nm(m->GetName());
            if (nm.BeginsWith("/")) nm.Remove(0,1);
            msg += Form("  %s  ->   TProof::Process(\"%s#%s\",...)\n",
                        nm.Data(), GetName(), nm.Data());
         }
         if (gProofServ)
            gProofServ->SendAsynMessage(msg);
         else
            Warning("Add", "%s", msg.Data());
         return kFALSE;
      }
   }

   // Get the metadata, if any
   m = fi->GetMetaData(meta);

   // Create the element
   const char *objname = 0;
   const char *dir = 0;
   Long64_t first = 0;
   Long64_t num = -1;
   if (!m) {
      objname = GetObjName();
      dir = GetDirectory();
   } else {
      objname = (m->GetObject() && strlen(m->GetObject())) ? m->GetObject() : GetObjName();
      dir = (m->GetDirectory() && strlen(m->GetDirectory())) ? m->GetDirectory() : GetDirectory();
      first = m->GetFirst();
      num = m->GetEntries();
   }
   const char *dataset = 0;
   if (strcmp(fi->GetTitle(), "TFileInfo")) dataset = fi->GetTitle();
   TDSetElement *el = new TDSetElement(file, objname, dir, first, -1, 0, dataset);
   el->SetEntries(num);

   // Set looked-up bit
   if (fi->TestBit(TFileInfo::kStaged) && setLookedUp)
      el->SetBit(TDSetElement::kHasBeenLookedUp);
   if (fi->TestBit(TFileInfo::kCorrupted))
      el->SetBit(TDSetElement::kCorrupted);

   // Add the element
   fElements->Add(el);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Export TDSetElements files as list of TFileInfo objects in file
/// 'fpath'. If the file exists already the action fails, unless
/// 'opt == "F"'.
/// Return 0 on success, -1 otherwise

Int_t TDSet::ExportFileList(const char *fpath, Option_t *opt)
{
   if (!fElements)
      return -1;
   if (fElements->GetSize() <= 0)
      return 0;

   Bool_t force = (opt[0] == 'F' || opt[0] == 'f');

   if (gSystem->AccessPathName(fpath, kFileExists) == kFALSE) {
      if (force) {
         // Try removing the file
         if (gSystem->Unlink(fpath)) {
            Info("ExportFileList","error removing dataset file: %s", fpath);
            return -1;
         }
      }
   }

   // Create the file list
   TList *fileinfo = new TList;
   fileinfo->SetOwner();

   TDSetElement *dse = 0;
   TIter next(fElements);
   while ((dse = (TDSetElement *) next())) {
      TFileInfoMeta *m = new TFileInfoMeta(dse->GetTitle(), dse->GetDirectory(), GetType(),
                                           dse->GetNum(), dse->GetFirst());
      TFileInfo *fi = new TFileInfo(dse->GetFileName());
      fi->AddMetaData(m);
      fileinfo->Add(fi);
   }

   // Write to file
   TFile *f = TFile::Open(fpath, "RECREATE");
   if (f) {
      f->cd();
      fileinfo->Write("fileList", TObject::kSingleKey);
      f->Close();
   } else {
      Info("ExportFileList","error creating dataset file: %s", fpath);
      SafeDelete(fileinfo);
      return -1;
   }

   // Cleanup
   SafeDelete(f);
   SafeDelete(fileinfo);

   // We are done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add friend dataset to this set. Only possible if the TDSet type is
/// a TTree or derived class. The friendset will be owned by this class
/// and deleted in its destructor.

void TDSet::AddFriend(TDSet *friendset, const char* alias)
{
   if (!friendset) {
      Error("AddFriend", "The friend TDSet is null!");
      return;
   }

   if (!fIsTree) {
      Error("AddFriend", "a friend set can only be added to a TTree TDSet");
      return;
   }
   TList *thisList = GetListOfElements();
   TList *friendsList = friendset->GetListOfElements();
   if (thisList->GetSize() != friendsList->GetSize() && friendsList->GetSize() != 1) {
      Error("AddFriend", "the friend dataset has %d elements while the main one has %d",
            thisList->GetSize(), friendsList->GetSize());
      return;
   }
   TIter next(thisList);
   TIter next2(friendsList);
   TDSetElement *friendElem = 0;
   if (friendsList->GetSize() == 1)
      friendElem = dynamic_cast<TDSetElement*> (friendsList->First());
   while(TDSetElement* e = dynamic_cast<TDSetElement*> (next())) {
      if (friendElem) // just one elem in the friend TDSet
         e->AddFriend(friendElem, alias);
      else
         e->AddFriend(dynamic_cast<TDSetElement*> (next2()), alias);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset or initialize access to the elements.

void TDSet::Reset()
{
   if (!fIterator) {
      fIterator = new TIter(fElements);
   } else {
      fIterator->Reset();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns number of entries in tree or objects in file. Returns -1 in
/// case of error.

Long64_t TDSet::GetEntries(Bool_t isTree, const char *filename, const char *path,
                           TString &objname)
{
   Double_t start = 0;
   if (gPerfStats) start = TTimeStamp();

   // Take into acoount possible prefixes
   TFile::EFileType typ = TFile::kDefault;
   TString fname = gEnv->GetValue("Path.Localroot",""), pfx(fname);
   // Get the locality (disable warnings or errors in connection attempts)
   Int_t oldLevel = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kError+1;
   if ((typ = TFile::GetType(filename, "", &fname)) != TFile::kLocal) fname = filename;
   gErrorIgnoreLevel = oldLevel;
   // Open the file
   TFile *file = TFile::Open(fname);

   if (gPerfStats)
      gPerfStats->FileOpenEvent(file, filename, start);

   if (file == 0) {
      ::SysError("TDSet::GetEntries",
                 "cannot open file %s (type: %d, pfx: %s)", filename, typ, pfx.Data());
      return -1;
   }

   TDirectory *dirsave = gDirectory;
   if (!file->cd(path)) {
      ::Error("TDSet::GetEntries", "cannot cd to %s", path);
      delete file;
      return -1;
   }

   TDirectory *dir = gDirectory;
   dirsave->cd();

   Long64_t entries;
   Bool_t fillname = kFALSE;
   if (isTree) {

      TString on(objname);
      TString sreg(objname);
      // If a wild card we will use the first object of the type
      // requested compatible with the reg expression we got
      if (sreg.Length() <= 0 || sreg == "" || sreg.Contains("*")) {
         fillname = kTRUE;
         if (sreg.Contains("*"))
            sreg.ReplaceAll("*", ".*");
         else
            sreg = ".*";
         TRegexp re(sreg);
         if (dir->GetListOfKeys()) {
            TIter nxk(dir->GetListOfKeys());
            TKey *k = 0;
            Bool_t notfound = kTRUE;
            while ((k = (TKey *) nxk())) {
               if (!strcmp(k->GetClassName(), "TTree")) {
                  TString kn(k->GetName());
                  if (kn.Index(re) != kNPOS) {
                     if (notfound) {
                        on = kn;
                        notfound = kFALSE;
                     } else if (kn != on) {
                       ::Warning("TDSet::GetEntries",
                                 "additional tree found in the file: %s", kn.Data());
                     }
                  }
               }
            }
         }
      }

      TKey *key = dir->GetKey(on);
      if (key == 0) {
         ::Error("TDSet::GetEntries", "cannot find tree \"%s\" in %s",
                 objname.Data(), filename);
         delete file;
         return -1;
      }
      TTree *tree = (TTree *) key->ReadObj();
      if (tree == 0) {
         // Error always reported?
         delete file;
         return -1;
      }
      entries = tree->GetEntries();
      delete tree;

      // Return full name in case of wildcards
      objname = (fillname) ? on : objname;

   } else {
      TList *keys = dir->GetListOfKeys();
      entries = keys->GetSize();
   }

   delete file;
   return entries;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw expression varexp for specified entries.
/// Returns -1 in case of error or number of selected events in case of success.
/// This function accepts a TCut objects as argument.
/// Use the operator+ to concatenate cuts.
/// Example:
///   dset.Draw("x",cut1+cut2+cut3);

Long64_t TDSet::Draw(const char *varexp, const TCut &selection, Option_t *option,
                     Long64_t nentries, Long64_t firstentry)
{
   return Draw(varexp, selection.GetTitle(), option, nentries, firstentry);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw expression varexp for specified entries.
/// Returns -1 in case of error or number of selected events in case of success.
/// For more see TTree::Draw().

Long64_t TDSet::Draw(const char *varexp, const char *selection, Option_t *option,
                     Long64_t nentries, Long64_t firstentry)
{
   if (!IsValid() || !fElements->GetSize()) {
      Error("Draw", "not a correctly initialized TDSet");
      return -1;
   }

   if (gProof)
      return gProof->DrawSelect(this, varexp, selection, option, nentries,
                                firstentry);

   Error("Draw", "no active PROOF session");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Start the TTreeViewer on this TTree.

void TDSet::StartViewer()
{
   if (gROOT->IsBatch()) {
      Warning("StartViewer", "viewer cannot run in batch mode");
      return;
   }

   if (!gProof) {
      Error("StartViewer", "no PROOF found");
      return;
   }
   if (!IsTree()) {
      Error("StartViewer", "TDSet contents should be of type TTree (or subtype)");
      return;
   }
   fProofChain = new TProofChain(this, kTRUE);

   TPluginHandler *h;
   if ((h = gROOT->GetPluginManager()->FindHandler("TVirtualTreeViewer"))) {
      if (h->LoadPlugin() == -1)
         return;
      h->ExecPlugin(1,fProofChain);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a tree header containing the branches' structure of the dataset.

TTree* TDSet::GetTreeHeader(TProof* proof)
{
   return proof->GetTreeHeader(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if all elements are valid.

Bool_t TDSet::ElementsValid()
{
   if (TestBit(TDSet::kValidityChecked))
      return (TestBit(TDSet::kSomeInvalid) ? kFALSE : kTRUE);

   SetBit(TDSet::kValidityChecked);
   ResetBit(TDSet::kSomeInvalid);
   TIter nextElem(GetListOfElements());
   while (TDSetElement *elem = dynamic_cast<TDSetElement*>(nextElem())) {
      if (!elem->GetValid()) {
         SetBit(TDSet::kSomeInvalid);
         return kFALSE;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove TDSetElement 'elem' from the list.
/// Return 0 on success, -1 if the element is not in the list

Int_t TDSet::Remove(TDSetElement *elem, Bool_t deleteElem)
{
   if (!elem || !(((THashList *)(GetListOfElements()))->Remove(elem)))
      return -1;

   if (deleteElem)
      SafeDelete(elem);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Validate the TDSet by opening files.

void TDSet::Validate()
{
   TIter nextElem(GetListOfElements());
   while (TDSetElement *elem = dynamic_cast<TDSetElement*>(nextElem())) {
      if (!elem->GetValid())
         elem->Validate(IsTree());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Resolve the end-point URL for the current elements of this data set
/// If the removeMissing option is set to kTRUE, remove the TDSetElements
/// that can not be located.
/// The method returns the list of removed TDSetElements in *listOfMissingFiles
/// if the latter is defined (the list must be created outside).

void TDSet::Lookup(Bool_t removeMissing, TList **listOfMissingFiles)
{
   // If an entry- or event- list has been given, assign the relevant portions
   // to each element; this allows to look-up only for the elements which have
   // something to be processed, so it is better to do it before the real look-up
   // operations.
   SplitEntryList();

   TString msg("Looking up for exact location of files");
   UInt_t n = 0;
   UInt_t ng = 0;
   UInt_t tot = GetListOfElements()->GetSize();
   UInt_t n2 = (tot > 50) ? (UInt_t) tot / 50 : 1;
   Bool_t st = kTRUE;
   TIter nextElem(GetListOfElements());
   while (TDSetElement *elem = dynamic_cast<TDSetElement*>(nextElem())) {
      if (elem->GetNum() != 0) { // -1 means "all entries"
         ng++;
         if (!elem->GetValid())
            if (elem->Lookup(kFALSE))
               if (removeMissing) {
                  if (Remove(elem, kFALSE))
                     Error("Lookup", "Error removing a missing file");
                  if (listOfMissingFiles)
                   (*listOfMissingFiles)->Add(elem->GetFileInfo(fType));
               }
      }
      n++;
      // Notify the client
      if (gProof && (n > 0 && !(n % n2)))
         gProof->SendDataSetStatus(msg, n, tot, st);
      // Break if we have been asked to stop
      if (gProof && gProof->GetRunStatus() != TProof::kRunning)
         break;
   }
   // Notify the client if not all the files have entries to be processed
   // (which may happen if an entry-list is used)
   if (ng < tot && gProofServ) {
      msg = Form("Files with entries to be processed: %d (out of %d)\n", ng, tot);
      gProofServ->SendAsynMessage(msg);
   } else {
      // Final notification to the client
      if (gProof) gProof->SendDataSetStatus(msg, n, tot, st);
   }
   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Flag all the elements as looked-up, so to avoid opening the files
/// if the functionality is not supported

void TDSet::SetLookedUp()
{
   TIter nextElem(GetListOfElements());
   while (TDSetElement *elem = dynamic_cast<TDSetElement*>(nextElem()))
      elem->SetLookedUp();
}

////////////////////////////////////////////////////////////////////////////////
/// Validate the TDSet against another TDSet.
/// Only validates elements in common from input TDSet.

void TDSet::Validate(TDSet* dset)
{
   THashList bestElements;
   bestElements.SetOwner();
   TList namedHolder;
   namedHolder.SetOwner();
   TIter nextOtherElem(dset->GetListOfElements());
   while (TDSetElement *elem = dynamic_cast<TDSetElement*>(nextOtherElem())) {
      if (!elem->GetValid()) continue;
      TString dir_file_obj = elem->GetDirectory();
      dir_file_obj += "_";
      dir_file_obj += TUrl(elem->GetFileName()).GetFileAndOptions();
      dir_file_obj += "_";
      dir_file_obj += elem->GetObjName();
      TPair *p = dynamic_cast<TPair*>(bestElements.FindObject(dir_file_obj));
      if (p) {
         TDSetElement *prevelem = dynamic_cast<TDSetElement*>(p->Value());
         if (prevelem) {
            Long64_t entries = prevelem->GetFirst()+prevelem->GetNum();
            if (entries<elem->GetFirst()+elem->GetNum()) {
               bestElements.Remove(p);
               bestElements.Add(new TPair(p->Key(), elem));
               delete p;
            }
         }
      } else {
         TNamed* named = new TNamed(dir_file_obj, dir_file_obj);
         namedHolder.Add(named);
         bestElements.Add(new TPair(named, elem));
      }
   }

   TIter nextElem(GetListOfElements());
   while (TDSetElement *elem = dynamic_cast<TDSetElement*>(nextElem())) {
      if (!elem->GetValid()) {
         TString dir_file_obj = elem->GetDirectory();
         dir_file_obj += "_";
         dir_file_obj += TUrl(elem->GetFileName()).GetFileAndOptions();
         dir_file_obj += "_";
         dir_file_obj += elem->GetObjName();
         if (TPair *p = dynamic_cast<TPair*>(bestElements.FindObject(dir_file_obj))) {
            TDSetElement* validelem = dynamic_cast<TDSetElement*>(p->Value());
            elem->Validate(validelem);
         }
      }
   }
}

//
// To handle requests coming from version 3 client / masters we need
// a special streamer
////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TDSetElement.

void TDSetElement::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      ResetBit(kWriteV3);
      if (R__v > 4) {
         R__b.ReadClassBuffer(TDSetElement::Class(), this, R__v, R__s, R__c);
      } else {
         // For version 3 client / masters we need a special streamer
         SetBit(kWriteV3);
         if (R__v > 3) {
            TNamed::Streamer(R__b);
         } else {
            // Old versions were not deriving from TNamed and had the
            // file name and the object type name in the first two members
            TObject::Streamer(R__b);
            TString name, title;
            R__b >> name >> title;
            SetNameTitle(name, title);
         }
         // Now we read the standard part
         R__b >> fDirectory;
         R__b >> fFirst;
         R__b >> fNum;
         R__b >> fMsd;
         R__b >> fTDSetOffset;
         TEventList *evl;
         R__b >> evl;
         R__b >> fValid;
         R__b >> fEntries;

         // Special treatment waiting for proper retrieving of stl containers
         FriendsList_t *friends = new FriendsList_t;
         static TClassRef classFriendsList = TClass::GetClass(typeid(FriendsList_t));
         R__b.ReadClassBuffer( classFriendsList, friends, classFriendsList->GetClassVersion(), 0, 0);
         if (friends) {
            // Convert friends to a TList (to be written)
            fFriends = new TList();
            fFriends->SetOwner();
            for (FriendsList_t::iterator i = friends->begin();
                 i != friends->end(); ++i) {
               TDSetElement *dse = (TDSetElement *) i->first->Clone();
               fFriends->Add(new TPair(dse, new TObjString(i->second.Data())));
            }
         }
         // the value for fIsTree (only older versions are sending it)
         Bool_t tmpIsTree;
         R__b >> tmpIsTree;
         R__b.CheckByteCount(R__s, R__c, TDSetElement::IsA());
      }
   } else {
      if (TestBit(kWriteV3)) {
         // For version 3 client / masters we need a special streamer
         R__b << Version_t(3);
         TObject::Streamer(R__b);
         R__b << TString(GetName());
         R__b << TString(GetTitle());
         R__b << fDirectory;
         R__b << fFirst;
         R__b << fNum;
         R__b << fMsd;
         R__b << fTDSetOffset;
         R__b << (TEventList *)0;
         R__b << fValid;
         R__b << fEntries;

         // Special treatment waiting for proper retrieving of stl containers
         FriendsList_t *friends = new FriendsList_t;
         if (fFriends) {
            TIter nxf(fFriends);
            TPair *p = 0;
            while ((p = (TPair *)nxf()))
               friends->push_back(std::make_pair((TDSetElement *)p->Key(),
                                   TString(((TObjString *)p->Value())->GetName())));
         }
         static TClassRef classFriendsList = TClass::GetClass(typeid(FriendsList_t));
         R__b.WriteClassBuffer( classFriendsList, &friends );

         // Older versions had an unused boolean called fIsTree: we fill it
         // with its default value
         R__b << kFALSE;
      } else {
         R__b.WriteClassBuffer(TDSetElement::Class(),this);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TDSet.

void TDSet::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      ResetBit(kWriteV3);
      if (R__v > 3) {
         R__b.ReadClassBuffer(TDSet::Class(), this, R__v, R__s, R__c);
      } else {
         // For version 3 client / masters we need a special streamer
         SetBit(kWriteV3);
         TNamed::Streamer(R__b);
         R__b >> fDir;
         R__b >> fType;
         R__b >> fObjName;
         TList elems;
         elems.Streamer(R__b);
         elems.SetOwner(kFALSE);
         if (elems.GetSize() > 0) {
            fElements = new THashList;
            fElements->SetOwner();
            TDSetElement *e = 0;
            TIter nxe(&elems);
            while ((e = (TDSetElement *)nxe())) {
               fElements->Add(e);
            }
         } else {
            fElements = 0;
         }
         R__b >> fIsTree;
      }
   } else {
      if (TestBit(kWriteV3)) {
         // For version 3 client / masters we need a special streamer
         R__b << Version_t(3);
         TNamed::Streamer(R__b);
         R__b << fDir;
         R__b << fType;
         R__b << fObjName;
         TList elems;
         if (fElements) {
            elems.SetOwner(kFALSE);
            if (fElements->GetSize() > 0) {
               TDSetElement *e = 0;
               TIter nxe(fElements);
               while ((e = (TDSetElement *)nxe()))
                  elems.Add(e);
            }
         }
         elems.Streamer(R__b);
         R__b << fIsTree;
      } else {
         R__b.WriteClassBuffer(TDSet::Class(),this);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set/Reset the 'OldStreamer' bit in this instance and its elements.
/// Needed for backward compatibility in talking to old client / masters.

void TDSet::SetWriteV3(Bool_t on)
{
   if (on)
      SetBit(TDSet::kWriteV3);
   else
      ResetBit(TDSet::kWriteV3);
   // Loop over dataset elements
   TIter nxe(GetListOfElements());
   TObject *o = 0;
   while ((o = nxe()))
      if (on)
         o->SetBit(TDSetElement::kWriteV3);
      else
         o->ResetBit(TDSetElement::kWriteV3);
}

////////////////////////////////////////////////////////////////////////////////
/// Set entry (or event) list for this data set

void TDSet::SetEntryList(TObject *aList)
{
   if (!aList) {
      // Nothing to do, except making sure to disable any previous setting
      fEntryList = 0;
      // Reset the element lists
      // TEntryList
      TIter next(fElements);
      TDSetElement *el=0;
      while ((el=(TDSetElement*)next())){
         el->SetEntryList(aList);
         el->Reset();
      }
      return;
   }

   if (TestBit(TDSet::kMultiDSet)) {

      // Global entry list for all the datasets
      TIter nxds(fElements);
      TDSet *ds = 0;
      while ((ds = (TDSet *) nxds()))
         ds->SetEntryList(aList);

   } else {

      // Link the proper object
      TEventList *evl = 0;
      TEntryList *enl = dynamic_cast<TEntryList*>(aList);
      if (!enl)
         evl = dynamic_cast<TEventList*>(aList);
      if (!enl && !evl) {
         Error("SetEntryList", "type of input object must be either TEntryList "
                              "or TEventList (found: '%s' - do nothing", aList->ClassName());
         return;
      }

      // Action depends on the type
      fEntryList = (enl) ? enl : (TEntryList *)evl;
   }
   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Splits the main entry (or event) list into sub-lists for the elements of
/// thet data set

void TDSet::SplitEntryList()
{
   if (TestBit(TDSet::kMultiDSet)) {
      // Global entry list for all the datasets
      TIter nxds(fElements);
      TDSet *ds = 0;
      while ((ds = (TDSet *) nxds()))
         ds->SplitEntryList();
      // Done
      return;
   }

   if (!fEntryList) {
      if (gDebug > 0)
         Info("SplitEntryList", "no entry- (or event-) list to split - do nothing");
      return;
   }

   // Action depend on type of list
   TEntryList *enl = dynamic_cast<TEntryList *>(fEntryList);
   if (enl) {
      // TEntryList
      TIter next(fElements);
      TDSetElement *el=0;
      TEntryList *sublist = 0;
      while ((el=(TDSetElement*)next())){
         sublist = enl->GetEntryList(el->GetObjName(), el->GetFileName());
         if (sublist){
            el->SetEntryList(sublist);
            el->SetNum(sublist->GetN());
         } else {
            sublist = new TEntryList("", "");
            el->SetEntryList(sublist);
            el->SetNum(0);
         }
      }
   } else {
      TEventList *evl = dynamic_cast<TEventList *>(fEntryList);
      if (evl) {
         // TEventList
         TIter next(fElements);
         TDSetElement *el, *prev;

         prev = dynamic_cast<TDSetElement*> (next());
         if (!prev)
            return;
         Long64_t low = prev->GetTDSetOffset();
         Long64_t high = low;
         Long64_t currPos = 0;
         do {
            el = dynamic_cast<TDSetElement*> (next());
            // kMaxLong64 means infinity
            high = (el == 0) ? kMaxLong64 : el->GetTDSetOffset();
#ifdef DEBUG
            while (currPos < evl->GetN() && evl->GetEntry(currPos) < low) {
               Error("SplitEntryList",
                     "TEventList: event outside of the range of any of the TDSetElements");
               currPos++;        // unnecessary check
            }
#endif
            TEventList* nevl = new TEventList();
            while (currPos < evl->GetN() && evl->GetEntry((Int_t)currPos) < high) {
               nevl->Enter(evl->GetEntry((Int_t)currPos) - low);
               currPos++;
            }
            prev->SetEntryList(nevl);
            prev->SetNum(nevl->GetN());
            low = high;
            prev = el;
         } while (el);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of files in the dataset

Int_t TDSet::GetNumOfFiles()
{
   Int_t nf = -1;
   if (fElements) {
      nf = 0;
      if (TestBit(TDSet::kMultiDSet)) {
         TIter nxds(fElements);
         TDSet *ds = 0;
         while ((ds = (TDSet *) nxds()))
            if (ds->GetListOfElements()) nf += ds->GetListOfElements()->GetSize();
      } else {
         nf = fElements->GetSize();
      }
   }
   // Done
   return nf;
}
