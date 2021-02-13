// @(#)root/memstat:$Id$
// Author: Anar Manafov (A.Manafov@gsi.de) 2008-03-02

/*************************************************************************
* Copyright (C) 1995-2010, Rene Brun and Fons Rademakers.               *
* All rights reserved.                                                  *
*                                                                       *
* For the licensing terms see $ROOTSYS/LICENSE.                         *
* For the list of contributors see $ROOTSYS/README/CREDITS.             *
*************************************************************************/
// STD
#include <cstdlib>
// ROOT
#include "TSystem.h"
#include "TError.h"
#include "Riostream.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TArrayL64.h"
#include "TH1.h"
#include "TMD5.h"
#include "TMath.h"
// Memstat
#include "TMemStatBacktrace.h"
#include "TMemStatMng.h"

using namespace Memstat;

ClassImp(TMemStatMng);

TMemStatMng* TMemStatMng::fgInstance = NULL;

//****************************************************************************//
//
//****************************************************************************//

TMemStatMng::TMemStatMng():
   TObject(),
#if !defined(__APPLE__)
   fPreviousMallocHook(TMemStatHook::GetMallocHook()),
   fPreviousFreeHook(TMemStatHook::GetFreeHook()),
#endif
   fDumpFile(NULL),
   fDumpTree(NULL),
   fUseGNUBuiltinBacktrace(kFALSE),
   fBeginTime(0),
   fPos(0),
   fTimems(0),
   fNBytes(0),
   fBtID(0),
   fMaxCalls(5000000),
   fBufferSize(10000),
   fBufN(0),
   fBufPos(0),
   fBufTimems(0),
   fBufNBytes(0),
   fBufBtID(0),
   fIndex(0),
   fMustWrite(0),
   fFAddrsList(0),
   fHbtids(0),
   fBTCount(0),
   fBTIDCount(0),
   fSysInfo(0)
{
   // Default constructor
}

////////////////////////////////////////////////////////////////////////////////
///Initialize MemStat manager - used only by instance method

void TMemStatMng::Init()
{
   fBeginTime = fTimeStamp.AsDouble();

   fDumpFile = new TFile(Form("memstat_%d.root", gSystem->GetPid()), "recreate");
   Int_t opt = 200000;
   if(!fDumpTree) {
      fDumpTree = new TTree("T", "Memory Statistics");
      fDumpTree->Branch("pos",   &fPos,   "pos/l", opt);
      fDumpTree->Branch("time",  &fTimems, "time/I", opt);
      fDumpTree->Branch("nbytes", &fNBytes, "nbytes/I", opt);
      fDumpTree->Branch("btid",  &fBtID,  "btid/I", opt);
   }

   fBTCount = 0;

   fBTIDCount = 0;

   fFAddrsList = new TObjArray();
   fFAddrsList->SetOwner(kTRUE);
   fFAddrsList->SetName("FAddrsList");

   fHbtids  = new TH1I("btids", "table of btids", 10000, 0, 1);   //where fHbtids is a member of the manager class
   fHbtids->SetDirectory(0);
   // save the histogram and the TObjArray to the tree header
   fDumpTree->GetUserInfo()->Add(fHbtids);
   fDumpTree->GetUserInfo()->Add(fFAddrsList);
   // save the system info to a tree header
   std::string sSysInfo(gSystem->GetBuildNode());
   sSysInfo += " | ";
   sSysInfo += gSystem->GetBuildCompilerVersion();
   sSysInfo += " | ";
   sSysInfo += gSystem->GetFlagsDebug();
   sSysInfo += " ";
   sSysInfo += gSystem->GetFlagsOpt();
   fSysInfo = new TNamed("SysInfo", sSysInfo.c_str());

   fDumpTree->GetUserInfo()->Add(fSysInfo);
   fDumpTree->SetAutoSave(10000000);
}

////////////////////////////////////////////////////////////////////////////////
/// GetInstance - a static function
/// Initialize a singleton of MemStat manager

TMemStatMng* TMemStatMng::GetInstance()
{
   if(!fgInstance) {
      fgInstance = new TMemStatMng;
      fgInstance->Init();
   }
   return fgInstance;
}

////////////////////////////////////////////////////////////////////////////////
/// Close - a static function
/// This method stops the manager,
/// flashes all the buffered data and closes the output tree.

void TMemStatMng::Close()
{
   // TODO: This is a temporary solution until we find a properalgorithm for SaveData
   //fgInstance->fDumpFile->WriteObject(fgInstance->fFAddrsList, "FAddrsList");

/*  std::ofstream f("mem_stat_debug.txt");
  int *btids = fgInstance->fHbtids->GetArray();
   if( !btids )
      return;
   int btid(1);
   int count(0);
   bool bStop(false);
   int count_empty(0);
   for (int i = 0; i < fgInstance->fBTChecksums.size(); ++i)
   {
     if (bStop)
         break;
     count = btids[btid-1];
     f << "++++++++++++++++++++++++\n";
     f << "BTID: " << btid << "\n";
     if ( count <= 0 )
        ++count_empty;
     for (int j = btid+1; j <= (btid+count); ++j )
         {
           TNamed *nm = (TNamed*)fgInstance->fFAddrsList->At(btids[j]);
           if( !nm )
               {
                   f << "Bad ID" << std::endl;
                   bStop = true;
                }
           f << "-------> " << nm->GetTitle() << "\n";
         }
     btid = btid + count + 1;
   }
   f.close();
   ::Info("TMemStatMng::Close", "btids without a stack %d\n", count_empty);
*/

   // to be documented
   fgInstance->FillTree();
   fgInstance->Disable();
   fgInstance->fDumpTree->AutoSave();
   fgInstance->fDumpTree->GetUserInfo()->Delete();

   ::Info("TMemStatMng::Close", "Tree saved to file %s\n", fgInstance->fDumpFile->GetName());
   ::Info("TMemStatMng::Close", "Tree entries = %d, file size = %g MBytes\n", (Int_t)fgInstance->fDumpTree->GetEntries(),1e-6*Double_t(fgInstance->fDumpFile->GetEND()));

   delete fgInstance->fDumpFile;
   //fgInstance->fDumpFile->Close();
   //delete fgInstance->fFAddrsList;
   //delete fgInstance->fSysInfo;

   delete fgInstance;
   fgInstance = NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// if an instance is destructed - the hooks are reseted to old hooks

TMemStatMng::~TMemStatMng()
{
   if(this != TMemStatMng::GetInstance())
      return;

   Info("~TMemStatMng", ">>> All free/malloc calls count: %d", fBTIDCount);
   Info("~TMemStatMng", ">>> Unique BTIDs count: %zu", fBTChecksums.size());

   Disable();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the maximum number of alloc/free calls to be buffered.
///if the alloc and free are in the buffer, the corresponding entries
///are not saved tio the Tree, reducing considerably the Tree output size

void TMemStatMng::SetBufferSize(Int_t buffersize)
{
   fBufferSize = buffersize;
   if (fBufferSize < 1) fBufferSize = 1;
   fBufN = 0;
   fBufPos    = new ULong64_t[fBufferSize];
   fBufTimems = new Int_t[fBufferSize];
   fBufNBytes = new Int_t[fBufferSize];
   fBufBtID   = new Int_t[fBufferSize];
   fIndex     = new Int_t[fBufferSize];
   fMustWrite = new Bool_t[fBufferSize];
}

////////////////////////////////////////////////////////////////////////////////
/// Set the maximum number of new/delete registered in the output Tree.

void TMemStatMng::SetMaxCalls(Int_t maxcalls)
{
   fMaxCalls = maxcalls;
}

////////////////////////////////////////////////////////////////////////////////
/// Enable memory hooks

void TMemStatMng::Enable()
{
   if(this != GetInstance())
      return;
#if defined(__APPLE__)
   TMemStatHook::trackZoneMalloc(MacAllocHook, MacFreeHook);
#else
   // set hook to our functions
   TMemStatHook::SetMallocHook(AllocHook);
   TMemStatHook::SetFreeHook(FreeHook);
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Disble memory hooks

void TMemStatMng::Disable()
{
   //FillTree();
   if(this != GetInstance())
      return;
#if defined(__APPLE__)
   TMemStatHook::untrackZoneMalloc();
#else
   // set hook to our functions
   TMemStatHook::SetMallocHook(fPreviousMallocHook);
   TMemStatHook::SetFreeHook(fPreviousFreeHook);
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// AllocHook - a static function
/// a special memory hook for Mac OS X memory zones.
/// Triggered when memory is allocated.

void TMemStatMng::MacAllocHook(void *ptr, size_t size)
{
   TMemStatMng* instance = TMemStatMng::GetInstance();
   // Restore all old hooks
   instance->Disable();

   // Call our routine
   instance->AddPointer(ptr, Int_t(size));

   // Restore our own hooks
   instance->Enable();
}

////////////////////////////////////////////////////////////////////////////////
/// AllocHook - a static function
/// a special memory hook for Mac OS X memory zones.
/// Triggered when memory is deallocated.

void TMemStatMng::MacFreeHook(void *ptr)
{
   TMemStatMng* instance = TMemStatMng::GetInstance();
   // Restore all old hooks
   instance->Disable();

   // Call our routine
   instance->AddPointer(ptr, -1);

   // Restore our own hooks
   instance->Enable();
}

////////////////////////////////////////////////////////////////////////////////
/// AllocHook - a static function
/// A glibc memory allocation hook.

void *TMemStatMng::AllocHook(size_t size, const void* /*caller*/)
{
   TMemStatMng* instance = TMemStatMng::GetInstance();
   // Restore all old hooks
   instance->Disable();

   // Call recursively
   void *result = malloc(size);
   // Call our routine
   instance->AddPointer(result, Int_t(size));
   //  TTimer::SingleShot(0, "TYamsMemMng", instance, "SaveData()");

   // Restore our own hooks
   instance->Enable();

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// FreeHook - a static function
/// A glibc memory deallocation hook.

void TMemStatMng::FreeHook(void* ptr, const void* /*caller*/)
{
   TMemStatMng* instance = TMemStatMng::GetInstance();
   // Restore all old hooks
   instance->Disable();

   // Call recursively
   free(ptr);

   // Call our routine
   instance->AddPointer(ptr, -1);

   // Restore our own hooks
   instance->Enable();
}

////////////////////////////////////////////////////////////////////////////////
/// An internal function, which returns a bitid for a corresponding CRC digest
/// cache variables

Int_t TMemStatMng::generateBTID(UChar_t *CRCdigest, Int_t stackEntries,
                                void **stackPointers)
{
   static Int_t old_btid = -1;
   static SCustomDigest old_digest;

   Int_t ret_val = -1;
   bool startCheck(false);
   if(old_btid >= 0) {
      for(int i = 0; i < g_digestSize; ++i) {
         if(old_digest.fValue[i] != CRCdigest[i]) {
            startCheck = true;
            break;
         }
      }
      ret_val = old_btid;
   } else {
      startCheck = true;
   }

   // return cached value
   if(!startCheck)
      return ret_val;

   old_digest = SCustomDigest(CRCdigest);
   CRCSet_t::const_iterator found = fBTChecksums.find(CRCdigest);

   if(fBTChecksums.end() == found) {
      // check the size of the BT array container
      const int nbins = fHbtids->GetNbinsX();
      //check that the current allocation in fHbtids is enough, otherwise expend it with
      if(fBTCount + stackEntries + 1 >= nbins) {
         fHbtids->SetBins(nbins * 2, 0, 1);
      }

      int *btids = fHbtids->GetArray();
      // the first value is a number of entries in a given stack
      btids[fBTCount++] = stackEntries;
      ret_val = fBTCount;
      if(stackEntries <= 0) {
         Warning("AddPointer",
                 "A number of stack entries is equal or less than zero. For btid %d", ret_val);
      }

      // add new BT's CRC value
      std::pair<CRCSet_t::iterator, bool> res = fBTChecksums.insert(CRCSet_t::value_type(CRCdigest, ret_val));
      if(!res.second)
         Error("AddPointer", "Can't added a new BTID to the container.");

      // save all symbols of this BT
      for(int i = 0; i < stackEntries; ++i) {
         ULong_t func_addr = (ULong_t)(stackPointers[i]);
         Int_t idx = fFAddrs.find(func_addr);
         // check, whether it's a new symbol
         if(idx < 0) {
            TString strFuncAddr;
            strFuncAddr += func_addr;
            TString strSymbolInfo;
            getSymbolFullInfo(stackPointers[i], &strSymbolInfo);

            TNamed *nm = new TNamed(strFuncAddr, strSymbolInfo);
            fFAddrsList->Add(nm);
            idx = fFAddrsList->GetEntriesFast() - 1;
            // TODO: more detailed error message...
            if(!fFAddrs.add(func_addr, idx))
               Error("AddPointer", "Can't add a function return address to the container");
         }

         // even if we have -1 as an index we add it to the container
         btids[fBTCount++] = idx;
      }

   } else {
      // reuse an existing BT
      ret_val = found->second;
   }

   old_btid = ret_val;

   return ret_val;
}

////////////////////////////////////////////////////////////////////////////////
/// Add pointer to table.
/// This method is called every time when any of the hooks are triggered.
/// The memory de-/allocation information will is recorded.

void TMemStatMng::AddPointer(void *ptr, Int_t size)
{
   void *stptr[g_BTStackLevel + 1];
   const int stackentries = getBacktrace(stptr, g_BTStackLevel, fUseGNUBuiltinBacktrace);

   // save only unique BTs
   TMD5 md5;
   md5.Update(reinterpret_cast<UChar_t*>(stptr), sizeof(void*) * stackentries);
   UChar_t digest[g_digestSize];
   md5.Final(digest);

   // for Debug. A counter of all (de)allacations.
   ++fBTIDCount;

   Int_t btid(generateBTID(digest, stackentries, stptr));

   if(btid <= 0)
      Error("AddPointer", "bad BT id");

   fTimeStamp.Set();
   Double_t CurTime = fTimeStamp.AsDouble();
   fBufTimems[fBufN] = Int_t(10000.*(CurTime - fBeginTime));
   ULong_t ul = (ULong_t)(ptr);
   fBufPos[fBufN]    = (ULong64_t)(ul);
   fBufNBytes[fBufN] = size;
   fBufBtID[fBufN]   = btid;
   fBufN++;
   if (fBufN >= fBufferSize) {
      FillTree();
   }
}

////////////////////////////////////////////////////////////////////////////////
///loop on all entries in the buffer and fill the output Tree
///entries with alloc and free in the buffer are eliminated

void TMemStatMng::FillTree()
{

   //eliminate alloc/free pointing to the same location in the current buffer
   TMath::Sort(fBufN,fBufPos,fIndex,kFALSE);
   memset(fMustWrite,0,fBufN*sizeof(Bool_t));
   Int_t i=0,j;
   while (i<fBufN) {
      Int_t indi = fIndex[i];
      Int_t indmin = indi;
      Int_t indmax = indi;
      j = i+1;;
      ULong64_t pos = fBufPos[indi];
      while (j < fBufN) {
         Int_t indj = fIndex[j];
         ULong64_t posj = fBufPos[indj];
         if (posj != pos) break;
         if (indmin > indj) indmin = indj;
         if (indmax < indj) indmax = indj;
         j++;
      }
      if (indmin == indmax) fMustWrite[indmin] = kTRUE;
      if (fBufNBytes[indmin] == -1) fMustWrite[indmin] = kTRUE;
      if (fBufNBytes[indmax] > 0)   fMustWrite[indmax] = kTRUE;
      i = j;
   }

   // now fill the Tree with the remaining allocs/frees
   for (i=0;i<fBufN;i++) {
      if (!fMustWrite[i]) continue;
      fPos    = fBufPos[i];
      fTimems = fBufTimems[i];
      fNBytes = fBufNBytes[i];
      fBtID   = fBufBtID[i];
      fDumpTree->Fill();
   }

   fBufN = 0;
   if (fDumpTree->GetEntries() >= fMaxCalls) TMemStatMng::GetInstance()->Disable();
}
