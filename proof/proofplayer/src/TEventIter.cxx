// @(#)root/proofplayer:$Id$
// Author: Maarten Ballintijn   07/01/02
// Modified: Long Tran-Thanh    04/09/07  (Addition of TEventIterUnit)

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TEventIter
\ingroup proofkernel

Special iterator class used in TProofPlayer to iterate over events
or objects in the packets.

*/

#include "TEnv.h"
#include "TEventIter.h"
#include "TFriendElement.h"
#include "TCollection.h"
#include "TDSet.h"
#include "TFile.h"
#include "TKey.h"
#include "TProofDebug.h"
#include "TSelector.h"
#include "TTimeStamp.h"
#include "TTree.h"
#include "TTreeCache.h"
#include "TTreeCacheUnzip.h"
#include "TVirtualPerfStats.h"
#include "TEventList.h"
#include "TEntryList.h"
#include "TList.h"
#include "TMap.h"
#include "TRegexp.h"
#include "TProofServ.h"
#include "TSystem.h"

#include "TError.h"

#if defined(R__MACOSX)
#include "fcntl.h"
#endif

ClassImp(TEventIter);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TEventIter::TEventIter()
{
   fDSet  = 0;
   fElem  = 0;
   fFile  = 0;
   fDir   = 0;
   fSel   = 0;
   fFirst = 0;
   fCur   = -1;
   fNum   = 0;
   fStop  = kFALSE;
   fOldBytesRead = 0;
   fEventList = 0;
   fEventListPos = 0;
   fEntryList = 0;
   fEntryListPos = 0;
   fElemFirst = 0;
   fElemNum = 0;
   fElemCur = -1;
   ResetBit(TEventIter::kData);

   if ((fPackets = new TList)) {
      TString n("ProcessedPackets_");
      if (gProofServ) n += gProofServ->GetOrdinal();
      fPackets->SetName(n);
      Info("TEventIter", "fPackets list '%s' created", n.Data());
   } else {
      Warning("TEventIter", "fPackets list could not be created");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TEventIter::TEventIter(TDSet *dset, TSelector *sel, Long64_t first, Long64_t num)
   : fDSet(dset), fSel(sel)
{
   fElem  = 0;
   fFile  = 0;
   fDir   = 0;
   fFirst = first;
   fCur   = -1;
   fNum   = num;
   fStop  = kFALSE;
   fEventList = 0;
   fEventListPos = 0;
   fEntryList = 0;
   fEntryListPos = 0;
   fOldBytesRead = 0;
   fElemFirst = 0;
   fElemNum = 0;
   fElemCur = -1;
   ResetBit(TEventIter::kData);

   if ((fPackets = new TList)) {
      TString n("ProcessedPackets_");
      if (gProofServ) n += gProofServ->GetOrdinal();
      fPackets->SetName(n);
      Info("TEventIter", "fPackets list '%s' created", n.Data());
   } else {
      Warning("TEventIter", "fPackets list could not be created");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TEventIter::~TEventIter()
{
   if (fPackets) {
      fPackets->SetOwner(kTRUE);
      SafeDelete(fPackets);
   }
   delete fFile;
}

////////////////////////////////////////////////////////////////////////////////
/// Invalidated the current packet (if any) by setting the TDSetElement::kCorrupted bit

void TEventIter::InvalidatePacket()
{
   if (fElem) fElem->SetBit(TDSetElement::kCorrupted);
}

////////////////////////////////////////////////////////////////////////////////
/// Set flag to stop the process

void TEventIter::StopProcess(Bool_t /*abort*/)
{
   fStop = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create and instance of the appropriate iterator

TEventIter *TEventIter::Create(TDSet *dset, TSelector *sel, Long64_t first, Long64_t num)
{
   if (dset->TestBit(TDSet::kEmpty)) {
      return new TEventIterUnit(dset, sel, num);
   }  else if (dset->IsTree()) {
      return new TEventIterTree(dset, sel, first, num);
   } else {
      return new TEventIterObj(dset, sel, first, num);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Load directory

Int_t TEventIter::LoadDir()
{
   Int_t ret = 0;

   // Check Filename
   if ( fFile == 0 || fFilename != fElem->GetFileName() ) {
      fDir = 0;
      delete fFile; fFile = 0;

      fFilename = fElem->GetFileName();

      TDirectory *dirsave = gDirectory;

      Double_t start = 0;
      if (gPerfStats) start = TTimeStamp();

      // Take into acoount possible prefixes
      TFile::EFileType typ = TFile::kDefault;
      TString fname = gEnv->GetValue("Path.Localroot","");
      if (!fname.IsNull())
         typ = TFile::GetType(fFilename, "", &fname);
      if (typ != TFile::kLocal)
         fname = fFilename;
      fFile = TFile::Open(fname);

      if (gPerfStats) {
         gPerfStats->FileOpenEvent(fFile, fFilename, start);
         fOldBytesRead = 0;
      }

      if (dirsave) dirsave->cd();

      if (!fFile || fFile->IsZombie() ) {
         if (fFile)
            Error("Process","Cannot open file: %s (%s)",
               fFilename.Data(), strerror(fFile->GetErrno()) );
         else
            Error("Process","Cannot open file: %s (errno unavailable)",
               fFilename.Data());
         // cleanup ?
         return -1;
      }
      PDB(kLoop,2) Info("LoadDir","Opening file: %s", fFilename.Data() );
      ret = 1;
   }

   // Check Directory
   if ( fDir == 0 || fPath != fElem->GetDirectory() ) {
      TDirectory *dirsave = gDirectory;

      fPath = fElem->GetDirectory();
      if ( !fFile->cd(fPath) ) {
         Error("Process","Cannot cd to: %s",
            fPath.Data() );
         return -1;
      }
      PDB(kLoop,2) Info("Process","Cd to: %s", fPath.Data() );
      fDir = gDirectory;
      if (dirsave) dirsave->cd();
      ret = 1;
   }

   return ret;
}

//______________________________________________________________________________
Long64_t TEventIter::GetEntryNumber(Long64_t next)
{
   // Get the entry number, taking into account event/entry lists

   Long64_t entry = next;
   // Set entry number; if data iteration we may need to test the entry or event lists
   if (TestBit(TEventIter::kData)) {
      if (fEntryList){
         entry = fEntryList->GetEntry(next);
      } else if (fEventList) {
         entry = fEventList->GetEntry(next);
      }
   }
   // Pre-event processing
   PreProcessEvent(entry);
   // Done
   return entry;
}

//------------------------------------------------------------------------


ClassImp(TEventIterUnit);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TEventIterUnit::TEventIterUnit()
{
   fDSet = 0;
   fElem = 0;
   fSel = 0;
   fNum = 0;
   fCurrent = 0;
   fStop = kFALSE;
   fOldBytesRead = 0; // Measures the bytes written
}

////////////////////////////////////////////////////////////////////////////////
/// Main constructor

TEventIterUnit::TEventIterUnit(TDSet* dset, TSelector *sel, Long64_t num)
{
   fDSet = dset;
   fElem = 0;
   fSel = sel;
   fNum = num;
   fCurrent = 0;
   fStop = kFALSE;
   fOldBytesRead = 0; // Measures the bytes written
}

////////////////////////////////////////////////////////////////////////////////
/// Get loop range

Int_t TEventIterUnit::GetNextPacket(Long64_t &fst, Long64_t &num)
{
   if (gPerfStats) {
      Long64_t totBytesWritten = TFile::GetFileBytesWritten();
      Long64_t bytesWritten = totBytesWritten - fOldBytesRead;
      PDB(kLoop, 2) Info("GetNextPacket", "bytes written: %lld", bytesWritten);
      gPerfStats->SetBytesRead(bytesWritten);
      fOldBytesRead = totBytesWritten;
   }

   if (fDSet->TestBit(TDSet::kIsLocal)) {
      if (fElem) {
         if (fPackets) {
            fPackets->Add(fElem);
            PDB(kLoop, 2)
               Info("GetNextEvent", "packet added to list (sz: %d)", fPackets->GetSize());
            fElem = 0;
         } else {
            SafeDelete(fElem);
         }
         return -1;
      } else {
         fElem = new TDSetElement("", "", "", 0, fNum);
         fElem->SetBit(TDSetElement::kEmpty);
      }
   } else {
      if (fPackets && fElem) {
         fPackets->Add(fElem);
         PDB(kLoop, 2)
            Info("GetNextEvent", "packet added to list (sz: %d)", fPackets->GetSize());
         fElem = 0;
      } else {
         SafeDelete(fElem);
      }
      if (!(fElem = fDSet->Next()))
         return -1;
   }
   fElem->SetBit(TDSetElement::kNewPacket);

   if (!fElem->TestBit(TDSetElement::kEmpty)) {
      Error("GetNextPacket", "data element must be set to kEmtpy");
      return -1;
   }

   // Set output
   num = fElem->GetNum();
   if (num == 0) return -1;
   fst = fElem->GetFirst();

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get next event

Long64_t TEventIterUnit::GetNextEvent()
{
   if (fStop || fNum == 0)
      return -1;

   if (fElem) fElem->ResetBit(TDSetElement::kNewPacket);

   while (fElem == 0 || fCurrent == 0) {

      if (gPerfStats) {
         Long64_t totBytesWritten = TFile::GetFileBytesWritten();
         Long64_t bytesWritten = totBytesWritten - fOldBytesRead;
         PDB(kLoop, 2) Info("GetNextEvent", "bytes written: %lld", bytesWritten);
         gPerfStats->SetBytesRead(bytesWritten);
         fOldBytesRead = totBytesWritten;
      }

      if (fDSet->TestBit(TDSet::kIsLocal)) {
         if (fElem) {
            if (fPackets) {
               fPackets->Add(fElem);
               PDB(kLoop, 2)
                  Info("GetNextEvent", "packet added to list (sz: %d)", fPackets->GetSize());
               fElem = 0;
            } else {
               SafeDelete(fElem);
            }
            return -1;
         } else {
            fElem = new TDSetElement("", "", "", 0, fNum);
            fElem->SetBit(TDSetElement::kEmpty);
         }
      } else {
         if (fPackets && fElem) {
            fPackets->Add(fElem);
            PDB(kLoop, 2)
               Info("GetNextEvent", "packet added to list (sz: %d)", fPackets->GetSize());
            fElem = 0;
         } else {
            SafeDelete(fElem);
         }
         if (!(fElem = fDSet->Next()))
            return -1;
      }
      fElem->SetBit(TDSetElement::kNewPacket);

      if (!fElem->TestBit(TDSetElement::kEmpty)) {
         Error("GetNextEvent", "data element must be set to kEmtpy");
         return -1;
      }

      fNum = fElem->GetNum();
      if (!(fCurrent = fNum)) {
         fNum = 0;
         return -1;
      }
      fFirst = fElem->GetFirst();
   }
   Long64_t event = fNum - fCurrent + fFirst ;
   --fCurrent;
   return event;
}

//------------------------------------------------------------------------


ClassImp(TEventIterObj);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TEventIterObj::TEventIterObj()
{
   fKeys     = 0;
   fNextKey  = 0;
   fObj      = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TEventIterObj::TEventIterObj(TDSet *dset, TSelector *sel, Long64_t first, Long64_t num)
   : TEventIter(dset,sel,first,num)
{
   fClassName = dset->GetType();
   fKeys     = 0;
   fNextKey  = 0;
   fObj      = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor

TEventIterObj::~TEventIterObj()
{
   // delete fKeys ?
   delete fNextKey;
   delete fObj;
}

////////////////////////////////////////////////////////////////////////////////
/// Get loop range

Int_t TEventIterObj::GetNextPacket(Long64_t &first, Long64_t &num)
{
   SafeDelete(fElem);

   if (fStop || fNum == 0) return -1;

   while (fElem == 0 || fCur < fFirst-1) {

      if (gPerfStats && fFile) {
         Long64_t bytesRead = fFile->GetBytesRead();
         gPerfStats->SetBytesRead(bytesRead - fOldBytesRead);
         fOldBytesRead = bytesRead;
      }

      if (fElem) {
         // Save it to the list of processed packets
         if (fPackets) {
            fPackets->Add(fElem);
            fElem = 0;
         } else {
            SafeDelete(fElem);
         }
      }
      fElem = fDSet->Next(fKeys->GetSize());
      if (fElem && fElem->GetEntryList()) {
         Error("GetNextPacket", "entry- or event-list not available");
         return -1;
      }

      if ( fElem == 0 ) {
         fNum = 0;
         return -1;
      }
      fElem->SetBit(TDSetElement::kNewPacket);

      Int_t r = LoadDir();

      if ( r == -1 ) {

         // Error has been reported
         fNum = 0;
         return -1;

      } else if ( r == 1 ) {

         // New file and/or directory
         fKeys = fDir->GetListOfKeys();
         fNextKey = new TIter(fKeys);
      }

      // Validate values for this element
      fElemFirst = fElem->GetFirst();
      fElemNum = fElem->GetNum();
      if (fElem->GetEntryList()) {
         if (!(fEntryList = dynamic_cast<TEntryList *>(fElem->GetEntryList())))
            fEventList = dynamic_cast<TEventList *>(fElem->GetEntryList());
      }
      fEventListPos = 0;
      if (fEntryList)
         fElemNum = fEntryList->GetEntriesToProcess();
      else if (fEventList)
         fElemNum = fEventList->GetN();

      Long64_t tnum = fKeys->GetSize();

      if ( fElemFirst > tnum ) {
         Error("GetNextPacket","First (%lld) higher then number of keys (%lld) in %s",
               fElemFirst, tnum, fElem->GetName());
         fNum = 0;
         return -1;
      }

      if ( fElemNum == -1 ) {
         fElemNum = tnum - fElemFirst;
      } else if ( fElemFirst+fElemNum  > tnum ) {
         Error("GetNextPacket","Num (%lld) + First (%lld) larger then number of keys (%lld) in %s",
            fElemNum, fElemFirst, tnum, fElem->GetDirectory());
         fElemNum = tnum - fElemFirst;
      }

      // Skip this element completely?
      if ( fCur + fElemNum < fFirst ) {
         fCur += fElemNum;
         continue;
      }

      // Position within this element
      fNextKey->Reset();
      for(fElemCur = -1; fElemCur < fElemFirst-1 ; fElemCur++, fNextKey->Next()) { }
   }

   first = ++fElemCur;
   num = fElemNum;

   // Done
   return 0;
}
////////////////////////////////////////////////////////////////////////////////
/// To be executed before by TProofPlayer calling TSelector::Process

void TEventIterObj::PreProcessEvent(Long64_t)
{
   --fNum;
   ++fCur;
   TKey *key = (TKey*) fNextKey->Next();
   TDirectory *dirsave = gDirectory;
   fDir->cd();
   fObj = key->ReadObj();
   if (dirsave) dirsave->cd();
   fSel->SetObject(fObj);
}

////////////////////////////////////////////////////////////////////////////////
/// Get next event

Long64_t TEventIterObj::GetNextEvent()
{
   if (fStop || fNum == 0) return -1;

   if (fElem) fElem->ResetBit(TDSetElement::kNewPacket);

   while ( fElem == 0 || fElemNum == 0 || fCur < fFirst-1 ) {

      if (gPerfStats && fFile) {
         Long64_t bytesRead = fFile->GetBytesRead();
         gPerfStats->SetBytesRead(bytesRead - fOldBytesRead);
         fOldBytesRead = bytesRead;
      }

      if (fElem) {
         // Save it to the list of processed packets
         if (fPackets) {
            fPackets->Add(fElem);
            fElem = 0;
         } else {
            SafeDelete(fElem);
         }
      }
      fElem = fDSet->Next(fKeys->GetSize());
      if (fElem && fElem->GetEntryList()) {
         Error("GetNextEvent", "Entry- or event-list not available");
         return -1;
      }

      if ( fElem == 0 ) {
         fNum = 0;
         return -1;
      }
      fElem->SetBit(TDSetElement::kNewPacket);

      Int_t r = LoadDir();

      if ( r == -1 ) {

         // Error has been reported
         fNum = 0;
         return -1;

      } else if ( r == 1 ) {

         // New file and/or directory
         fKeys = fDir->GetListOfKeys();
         fNextKey = new TIter(fKeys);
      }

      // Validate values for this element
      fElemFirst = fElem->GetFirst();
      fElemNum = fElem->GetNum();
      if (fElem->GetEntryList()) {
         if (!(fEntryList = dynamic_cast<TEntryList *>(fElem->GetEntryList())))
            fEventList = dynamic_cast<TEventList *>(fElem->GetEntryList());
      }
      fEventListPos = 0;
      if (fEntryList)
         fElemNum = fEntryList->GetEntriesToProcess();
      else if (fEventList)
         fElemNum = fEventList->GetN();

      Long64_t num = fKeys->GetSize();

      if ( fElemFirst > num ) {
         Error("GetNextEvent","First (%lld) higher then number of keys (%lld) in %s",
               fElemFirst, num, fElem->GetName());
         fNum = 0;
         return -1;
      }

      if ( fElemNum == -1 ) {
         fElemNum = num - fElemFirst;
      } else if ( fElemFirst+fElemNum  > num ) {
         Error("GetNextEvent","Num (%lld) + First (%lld) larger then number of keys (%lld) in %s",
            fElemNum, fElemFirst, num, fElem->GetDirectory());
         fElemNum = num - fElemFirst;
      }

      // Skip this element completely?
      if ( fCur + fElemNum < fFirst ) {
         fCur += fElemNum;
         continue;
      }

      // Position within this element. TODO: more efficient?
      fNextKey->Reset();
      for(fElemCur = -1; fElemCur < fElemFirst-1 ; fElemCur++, fNextKey->Next()) { }
   }

   --fElemNum;
   ++fElemCur;

   // Pre-event processing
   PreProcessEvent(fElemCur);

   return fElemCur;
}

//------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TEventIterTree::TFileTree::TFileTree(const char *name, TFile *f, Bool_t islocal)
               : TNamed(name, ""), fUsed(kFALSE), fIsLocal(islocal), fFile(f)
{
   fTrees = new TList;
   fTrees->SetOwner();
}
////////////////////////////////////////////////////////////////////////////////
/// Default dtor.

TEventIterTree::TFileTree::~TFileTree()
{
   // Avoid destroying the cache; must be placed before deleting the trees
   TTree *tree = (TTree *)fTrees->First();
   while (tree) {
      fFile->SetCacheRead(0, tree);
      tree = (TTree *)fTrees->After(tree);
   }
   SafeDelete(fTrees);
   SafeDelete(fFile);
}

ClassImp(TEventIterTree);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TEventIterTree::TEventIterTree()
{
   fTree = 0;
   fTreeCache = 0;
   fUseTreeCache = 1;
   fCacheSize = -1;
   fTreeCacheIsLearning = kTRUE;
   fFileTrees = 0;
   fUseParallelUnzip = 0;
   fDontCacheFiles = kFALSE;
   SetBit(TEventIter::kData);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TEventIterTree::TEventIterTree(TDSet *dset, TSelector *sel, Long64_t first, Long64_t num)
   : TEventIter(dset,sel,first,num)
{
   fTreeName = dset->GetObjName();
   fTree = 0;
   fTreeCache = 0;
   fTreeCacheIsLearning = kTRUE;
   fFileTrees = new TList;
   fFileTrees->SetOwner();
   fUseTreeCache = gEnv->GetValue("ProofPlayer.UseTreeCache", 1);
   fCacheSize = gEnv->GetValue("ProofPlayer.CacheSize", -1);
   fUseParallelUnzip = gEnv->GetValue("ProofPlayer.UseParallelUnzip", 0);
   if (fUseParallelUnzip) {
      TTreeCacheUnzip::SetParallelUnzip(TTreeCacheUnzip::kEnable);
   } else {
      TTreeCacheUnzip::SetParallelUnzip(TTreeCacheUnzip::kDisable);
   }
   fDontCacheFiles = gEnv->GetValue("ProofPlayer.DontCacheFiles", 0);
   SetBit(TEventIter::kData);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TEventIterTree::~TEventIterTree()
{
   // Delete the tree cache ...
   SafeDelete(fTreeCache);
   // ... and the remaining open files
   SafeDelete(fFileTrees);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the size in bytes of the cache, if any
/// Return -1 if not used

Long64_t TEventIterTree::GetCacheSize()
{
   if (fUseTreeCache) return fCacheSize;
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of entries in the learning phase

Int_t TEventIterTree::GetLearnEntries()
{
   return TTreeCache::GetLearnEntries();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a Tree for the main TDSetElement and for all the friends.
/// Returns the main tree or 0 in case of an error.

TTree* TEventIterTree::GetTrees(TDSetElement *elem)
{
   // Reset used flags
   TIter nxft(fFileTrees);
   TFileTree *ft = 0;
   while ((ft = (TFileTree *)nxft()))
      ft->fUsed = kFALSE;

   Bool_t localfile = kFALSE;
   TTree* main = Load(elem, localfile);

   if (main && main != fTree) {
      // Set the file cache
      if (fUseTreeCache) {
         TFile *curfile = main->GetCurrentFile();
         if (curfile) {
            if (!fTreeCache) {
               main->SetCacheSize(fCacheSize);
               fTreeCache = (TTreeCache *)curfile->GetCacheRead(main);
               if (fCacheSize < 0) fCacheSize = main->GetCacheSize();
            } else {
               fTreeCache->ResetCache();
               curfile->SetCacheRead(fTreeCache, main);
               fTreeCache->UpdateBranches(main);
            }
            if (fTreeCache) {
               fTreeCacheIsLearning = fTreeCache->IsLearning();
               if (fTreeCacheIsLearning)
                  Info("GetTrees","the tree cache is in learning phase");
            }
         } else {
            Warning("GetTrees", "default tree does nto have a file attached: corruption? Tree cache untouched");
         }
      } else {
         // Disable the cache
         main->SetCacheSize(0);
      }
   }
   Bool_t loc = kFALSE;
   // Also the friends
   TList *friends = elem->GetListOfFriends();
   if (friends) {
      TIter nxf(friends);
      TDSetElement *dse = 0;
      while ((dse = (TDSetElement *) nxf())) {
         // The alias, if any, is in the element name options ('friend_alias=<alias>|')
         TUrl uf(dse->GetName());
         TString uo(uf.GetOptions()), alias;
         Int_t from = kNPOS;
         if ((from = uo.Index("friend_alias=")) != kNPOS) {
            from += strlen("friend_alias=");
            if (!uo.Tokenize(alias, from, "|"))
               Warning("GetTrees", "empty 'friend_alias' found for tree friend");
            // The options may be used for other things, so remove the internal strings once decoded
            uo.ReplaceAll(TString::Format("friend_alias=%s|", alias.Data()), "");
            uf.SetOptions(uo);
            dse->SetName(uf.GetUrl());
         }
         TTree *friendTree = Load(dse, loc, dse->GetObjName());
         if (friendTree && main) {
            // Make sure it has not yet been added
            Bool_t addfriend = kTRUE;
            TList *frnds = main->GetListOfFriends();
            if (frnds) {
               TIter xnxf(frnds);
               TFriendElement *fe = 0;
               while ((fe = (TFriendElement *) xnxf())) {
                  if (fe->GetTree() == friendTree) {
                     addfriend = kFALSE;
                     break;
                  }
               }
            }
            if (addfriend) {
               if (alias.IsNull())
                  main->AddFriend(friendTree);
               else
                  main->AddFriend(friendTree, alias);
            }
         } else {
            return 0;
         }
      }
   }

   // Remove instances not used
   nxft.Reset();
   while ((ft = (TFileTree *)nxft())) {
      if (!(ft->fUsed)) {
         fFileTrees->Remove(ft);
         delete ft;
      }
   }

   // Done, successfully or not
   return main;
}

////////////////////////////////////////////////////////////////////////////////
/// Load a tree from s TDSetElement

TTree* TEventIterTree::Load(TDSetElement *e, Bool_t &localfile, const char *objname)
{
   if (!e) {
      Error("Load", "undefined element");
      return (TTree *)0;
   }

   const char *fn = e->GetFileName();
   const char *dn = e->GetDirectory();
   const char *tn = 0;
   if (objname && strlen(objname) > 0) {
      tn = objname;
   } else {
      tn = (fDSet->GetObjName() && strlen(fDSet->GetObjName()) > 0)
         ? fDSet->GetObjName() : e->GetObjName();
      if (!tn || (tn && strlen(tn) <= 0)) tn = "*";
   }
   PDB(kLoop,2)
      Info("Load","loading: fn:'%s' dn:'%s' tn:'%s'", fn, dn, tn);

   TFile *f = 0;

   // Check if the file is already open
   TString names(fn);
   TString name;
   Ssiz_t from = 0;
   TFileTree *ft = 0;
   while (names.Tokenize(name,from,"|")) {
      TString key(TUrl(name).GetFileAndOptions());
      if ((ft = (TFileTree *) fFileTrees->FindObject(key.Data()))) {
         f = ft->fFile;
         break;
      }
   }

   // Open the file, if needed
   if (!f) {
      TFile::EFileType typ = TFile::kDefault;
      TString fname = gEnv->GetValue("Path.Localroot","");
      if (!fname.IsNull())
         typ = TFile::GetType(fn, "", &fname);
      if (typ != TFile::kLocal) {
         fname = fn;
      } else {
         localfile = kTRUE;
      }

      // Open the file
      f = TFile::Open(fname);
      if (!f) {
         Error("Load","file '%s' ('%s') could not be open", fn, fname.Data());
         return (TTree *)0;
      }

#if defined(R__MACOSX)
      // If requested set the no cache mode
      if (fDontCacheFiles && localfile) {
         fcntl(f->GetFd(), F_NOCACHE, 1);
      }
#endif

      // Create TFileTree instance in the list
      ft = new TFileTree(TUrl(f->GetName()).GetFileAndOptions(), f, localfile);
      fFileTrees->Add(ft);
   } else {
      // Fill locality boolean
      localfile = ft->fIsLocal;
      PDB(kLoop,2)
         Info("Load","file '%s' already open (local:%d)", fn, localfile);
   }

   // Check if the tree is already loaded
   if (ft && ft->fTrees->GetSize() > 0) {
      TTree *t = 0;
      if (!strcmp(tn, "*"))
         t = (TTree *) ft->fTrees->First();
      else
         t = (TTree *) ft->fTrees->FindObject(tn);
      if (t) {
         ft->fUsed = kTRUE;
         return t;
      }
   }

   TDirectory *dd = f;
   // Change dir, if required
   if (dn && !(dd = f->GetDirectory(dn))) {
      Error("Load","Cannot get to: %s", dn);
      return (TTree *)0;
   }
   PDB(kLoop,2)
      Info("Load","got directory: %s", dn);

   // If a wild card we will use the first object of the type
   // requested compatible with the reg expression we got
   TString on(tn);
   TString sreg(tn);
   if (sreg.Length() <= 0 || sreg == "" || sreg.Contains("*")) {
      if (sreg.Contains("*"))
         sreg.ReplaceAll("*", ".*");
      else
         sreg = ".*";
      TRegexp re(sreg);
      if (dd->GetListOfKeys()) {
         TIter nxk(dd->GetListOfKeys());
         TKey *k = 0;
         while ((k = (TKey *) nxk())) {
            if (!strcmp(k->GetClassName(), "TTree")) {
               TString kn(k->GetName());
               if (kn.Index(re) != kNPOS) {
                  on = kn;
                  break;
               }
            }
         }
      }
   }

   // Point to the key
   TKey *key = dd->GetKey(gSystem->BaseName(on));
   if (key == 0) {
      Error("Load", "Cannot find tree \"%s\" in %s", tn, fn);
      return (TTree*)0;
   }

   PDB(kLoop,2) Info("Load", "Reading: %s", tn);

   TTree *tree = dynamic_cast<TTree*> (key->ReadObj());
   dd->cd();

   if (tree == 0) {
      Error("Load", "Cannot <dynamic_cast> obj to tree \"%s\"", tn);
      return (TTree*)0;
   }

   // Add track in the cache
   ft->fTrees->Add(tree);
   ft->fUsed = kTRUE;
   PDB(kLoop,2)
      Info("Load","TFileTree for '%s' flagged as 'in-use' ...", ft->GetName());

   // Done
   return tree;
}

////////////////////////////////////////////////////////////////////////////////
/// Get loop range

Int_t TEventIterTree::GetNextPacket(Long64_t &first, Long64_t &num)
{
   if (first > -1) fEntryListPos = first;

   if (fStop || fNum == 0) return -1;

   Bool_t attach = kFALSE;

   // When files are aborted during processing (via TSelector::kAbortFile) the player
   // invalidates the element by settign this bit. We need to ask for a new packet
   Bool_t corrupted = kFALSE;
   Long64_t rest = -1;
   if (fElem) {
      corrupted = (fElem->TestBit(TDSetElement::kCorrupted)) ? kTRUE : kFALSE;
      rest = fElem->GetNum();
      if (fElemCur >= 0) rest -= (fElemCur + 1 - fElemFirst);
   }

   SafeDelete(fElem);

   while (fElem == 0 || fElemNum == 0 || fCur < fFirst-1) {

      if (gPerfStats && fTree) {
         Long64_t totBytesRead = fTree->GetCurrentFile()->GetBytesRead();
         Long64_t bytesRead = totBytesRead - fOldBytesRead;
         gPerfStats->SetBytesRead(bytesRead);
         fOldBytesRead = totBytesRead;
      }

      if (fElem) {
         // Save it to the list of processed packets
         if (fPackets) {
            fPackets->Add(fElem);
            fElem = 0;
         } else {
            SafeDelete(fElem);
         }
      }
      while (!fElem) {
         // For a corrupted/invalid file the request for a new packet is with totalEntries = -1
         // (the default) so that the packetizer invalidates the element
         if (corrupted) {
            fElem = fDSet->Next(rest);
         } else if (fTree) {
            fElem = fDSet->Next(fTree->GetEntries());
         } else {
            fElem = fDSet->Next();
         }

         if (!fElem) {
            // End of processing
            fNum = 0;
            return -1;
         }
         corrupted = kFALSE;
         fElem->SetBit(TDSetElement::kNewPacket);
         fElem->ResetBit(TDSetElement::kCorrupted);

         TTree *newTree = GetTrees(fElem);
         if (newTree) {
            if (newTree != fTree) {
               // The old tree is owned by TFileTree and will be deleted there
               fTree = newTree;
               attach = kTRUE;
               fOldBytesRead = (fTree->GetCurrentFile()) ? fTree->GetCurrentFile()->GetBytesRead() : 0;
            }
            // Set range to be analysed
            if (fTreeCache)
               fTreeCache->SetEntryRange(fElem->GetFirst(),
                                         fElem->GetFirst() + fElem->GetNum() - 1);
         } else {
            // Could not open this element: ask for another one
            SafeDelete(fElem);
            // The current tree, if any, is not valid anymore
            fTree = 0;
         }
      }

      // Validate values for this element
      fElemFirst = fElem->GetFirst();
      fElemNum = fElem->GetNum();
      fEntryList = 0;
      fEventList = 0;
      if (fElem->GetEntryList()) {
         if (!(fEntryList = dynamic_cast<TEntryList *>(fElem->GetEntryList())))
            fEventList = dynamic_cast<TEventList *>(fElem->GetEntryList());
      }
      fEntryListPos = fElemFirst;
      fEventListPos = 0;
      if (fEntryList)
         fElemNum = fEntryList->GetEntriesToProcess();
      else if (fEventList)
         fElemNum = fEventList->GetN();

      Long64_t tnum = (Long64_t) fTree->GetEntries();

      if (!fEntryList && !fEventList) {
         if ( fElemFirst > tnum ) {
            Error("GetNextPacket", "first (%lld) higher then number of entries (%lld) in %s",
                                  fElemFirst, tnum, fElem->GetObjName());
            fNum = 0;
            return -1;
         }
         if ( fElemNum == -1 ) {
            fElemNum = tnum - fElemFirst;
         } else if ( fElemFirst+fElemNum  > tnum ) {
            Error("GetNextPacket", "num (%lld) + first (%lld) larger then number of entries (%lld) in %s",
                                  fElemNum, fElemFirst, tnum, fElem->GetName());
            fElemNum = tnum - fElemFirst;
         }

         // Skip this element completely?
         if ( fCur + fElemNum < fFirst ) {
            fCur += fElemNum;
            continue;
         }
         // Position within this element
         fElemCur = fElemFirst-1;
      }
   }

   if (attach) {
      PDB(kLoop,1) Info("GetNextPacket", "call Init(%p) and Notify()",fTree);
      fSel->Init(fTree);
      fSel->Notify();
      TIter next(fSel->GetOutputList());
      TEntryList *elist=0;
      while ((elist=(TEntryList*)next())){
         if (elist->InheritsFrom(TEntryList::Class()))
            elist->SetTree(fTree->GetName(), fElem->GetFileName());
      }
      if (fSel->GetAbort() == TSelector::kAbortProcess) {
         // the error has been reported
         return -1;
      }
      attach = kFALSE;
   }

   // Fill the output now
   num = fElemNum;
   if (fEntryList) {
      first = fEntryListPos;
   } else if (fEventList){
      first = fEventListPos;
   } else {
      first = fElemFirst;
   }

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Actions to be done just before processing entry 'entry'.
/// Called by TProofPlayer.

void TEventIterTree::PreProcessEvent(Long64_t entry)
{
   if (!(fEntryList || fEventList)) {
      --fNum;
      ++fCur;
   }

   // Signal ending of learning phase
   if (fTreeCache && fTreeCacheIsLearning) {
      if (!(fTreeCache->IsLearning())) {
         fTreeCacheIsLearning = kFALSE;
         if (gProofServ) gProofServ->RestartComputeTime();
      }
   }

   // For prefetching
   if (fTree->LoadTree(entry) < 0) {
      Warning("PreEventProcess", "problems setting entry in TTree");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get next event

Long64_t TEventIterTree::GetNextEvent()
{
   if (fStop || fNum == 0) return -1;

   Bool_t attach = kFALSE;

   // When files are aborted during processing (via TSelector::kAbortFile) the player
   // invalidates the element by settign this bit. We need to ask for a new packet
   Bool_t corrupted = (fElem && fElem->TestBit(TDSetElement::kCorrupted)) ? kTRUE : kFALSE;

   if (fElem) fElem->ResetBit(TDSetElement::kNewPacket);

   while ( fElem == 0 || fElemNum == 0 || fCur < fFirst-1 || corrupted) {

      if (gPerfStats && fTree) {
         Long64_t totBytesRead = fTree->GetCurrentFile()->GetBytesRead();
         Long64_t bytesRead = totBytesRead - fOldBytesRead;
         gPerfStats->SetBytesRead(bytesRead);
         fOldBytesRead = totBytesRead;
      }

      Long64_t rest = -1;
      if (fElem) {
         rest = fElem->GetNum();
         if (fElemCur >= 0) rest -= (fElemCur + 1 - fElemFirst);
         // Save it to the list of processed packets
         if (fPackets) {
            fPackets->Add(fElem);
            fElem = 0;
         } else {
            SafeDelete(fElem);
         }
      }

      while (!fElem) {
         // For a corrupted/invalid file the request for a new packet is with totalEntries = -1
         // (the default) so that the packetizer invalidates the element
         if (corrupted) {
            fElem = fDSet->Next(rest);
         } else if (fTree) {
            fElem = fDSet->Next(fTree->GetEntries());
         } else {
            fElem = fDSet->Next();
         }

         if (!fElem) {
            // End of processing
            fNum = 0;
            return -1;
         }
         corrupted = kFALSE;
         fElem->SetBit(TDSetElement::kNewPacket);
         fElem->ResetBit(TDSetElement::kCorrupted);

         TTree *newTree = GetTrees(fElem);
         if (newTree) {
            if (newTree != fTree) {
               // The old tree is owned by TFileTree and will be deleted there
               fTree = newTree;
               attach = kTRUE;
               fOldBytesRead = (fTree->GetCurrentFile()) ? fTree->GetCurrentFile()->GetBytesRead() : 0;
            }
            // Set range to be analysed
            if (fTreeCache)
               fTreeCache->SetEntryRange(fElem->GetFirst(),
                                         fElem->GetFirst() + fElem->GetNum() - 1);
         } else {
            // Could not open this element: ask for another one
            SafeDelete(fElem);
            // The current tree, if any, is not valid anymore
            fTree = 0;
         }
      }

      // Validate values for this element
      fElemFirst = fElem->GetFirst();
      fElemNum = fElem->GetNum();
      fEntryList = dynamic_cast<TEntryList *>(fElem->GetEntryList());
      fEventList = (fEntryList) ? (TEventList *)0
                                : dynamic_cast<TEventList *>(fElem->GetEntryList());
      fEntryListPos = fElemFirst;
      fEventListPos = 0;
      if (fEntryList)
         fElemNum = fEntryList->GetEntriesToProcess();
      else if (fEventList)
         fElemNum = fEventList->GetN();

      Long64_t num = (Long64_t) fTree->GetEntries();

      if (!fEntryList && !fEventList) {
         if ( fElemFirst > num ) {
            Error("GetNextEvent", "first (%lld) higher then number of entries (%lld) in %s",
                                  fElemFirst, num, fElem->GetObjName());
            fNum = 0;
            return -1;
         }
         if ( fElemNum == -1 ) {
            fElemNum = num - fElemFirst;
         } else if ( fElemFirst+fElemNum  > num ) {
            Error("GetNextEvent", "num (%lld) + first (%lld) larger then number of entries (%lld) in %s",
                                  fElemNum, fElemFirst, num, fElem->GetName());
            fElemNum = num - fElemFirst;
         }

         // Skip this element completely?
         if ( fCur + fElemNum < fFirst ) {
            fCur += fElemNum;
            continue;
         }
         // Position within this element. TODO: more efficient?
         fElemCur = fElemFirst-1;
      }
   }

   if ( attach ) {
      PDB(kLoop,1) Info("GetNextEvent", "call Init(%p) and Notify()",fTree);
      fSel->Init(fTree);
      fSel->Notify();
      TIter next(fSel->GetOutputList());
      TEntryList *elist=0;
      while ((elist=(TEntryList*)next())){
         if (elist->InheritsFrom(TEntryList::Class()))
            elist->SetTree(fTree->GetName(), fElem->GetFileName());
      }
      if (fSel->GetAbort() == TSelector::kAbortProcess) {
         // the error has been reported
         return -1;
      }
      attach = kFALSE;
   }
   Long64_t rv;

   if (fEntryList){
      --fElemNum;
      rv = fEntryList->GetEntry(fEntryListPos);
      fEntryListPos++;
   } else if (fEventList) {
      --fElemNum;
      rv = fEventList->GetEntry(fEventListPos);
      fEventListPos++;
   } else {
      --fElemNum;
      ++fElemCur;
      rv = fElemCur;
   }

   // Pre-event processing
   PreProcessEvent(rv);

   return rv;
}
