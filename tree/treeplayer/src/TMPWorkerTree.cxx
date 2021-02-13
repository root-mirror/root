/* @(#)root/multiproc:$Id$ */
// Author: Enrico Guiraud July 2015
// Modified: G Ganis Jan 2017

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "MPCode.h"
#include "MPSendRecv.h"
#include "TError.h"
#include "TMPWorkerTree.h"
#include "TEnv.h"
#include <string>

//////////////////////////////////////////////////////////////////////////
///
/// \class TMPWorkerTree
///
/// This class works in conjuction with TTreeProcessorMP, reacting to messages
/// received from it as specified by the Notify and HandleInput methods.
///
/// \class TMPWorkerTreeFunc
///
/// Templated derivation of TMPWorkerTree handlign generic function tree processing. 
///
/// \class TMPWorkerTreeSel
///
/// Templated derivation of TMPWorkerTree handlign selector tree processing. 
///
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
/// Class constructors.
/// Note that this does not set variables like fPid or fS (worker's socket).\n
/// These operations are handled by the Init method, which is called after
/// forking.\n
/// This separation is in place because the instantiation of a worker
/// must be done once _before_ forking, while the initialization of the
/// members must be done _after_ forking by each of the children processes.
TMPWorkerTree::TMPWorkerTree()
   : TMPWorker(), fFileNames(), fTreeName(), fTree(nullptr), fFile(nullptr), fEntryList(nullptr), fFirstEntry(0),
     fTreeCache(0), fTreeCacheIsLearning(kFALSE), fUseTreeCache(kTRUE), fCacheSize(-1)
{
   Setup();
}

TMPWorkerTree::TMPWorkerTree(const std::vector<std::string> &fileNames, TEntryList *entries,
                             const std::string &treeName, UInt_t nWorkers, ULong64_t maxEntries, ULong64_t firstEntry)
   : TMPWorker(nWorkers, maxEntries), fFileNames(fileNames), fTreeName(treeName), fTree(nullptr), fFile(nullptr),
     fEntryList(entries), fFirstEntry(firstEntry), fTreeCache(0), fTreeCacheIsLearning(kFALSE), fUseTreeCache(kTRUE),
     fCacheSize(-1)
{
   Setup();
}

TMPWorkerTree::TMPWorkerTree(TTree *tree, TEntryList *entries, UInt_t nWorkers, ULong64_t maxEntries,
                             ULong64_t firstEntry)
   : TMPWorker(nWorkers, maxEntries), fTree(tree), fFile(nullptr), fEntryList(entries), fFirstEntry(firstEntry),
     fTreeCache(0), fTreeCacheIsLearning(kFALSE), fUseTreeCache(kTRUE), fCacheSize(-1)
{
   Setup();
}

TMPWorkerTree::~TMPWorkerTree()
{
   // Properly close the open file, if any
   CloseFile();
}

//////////////////////////////////////////////////////////////////////////
/// Auxilliary method for common initializations
void TMPWorkerTree::Setup()
{
   Int_t uc = gEnv->GetValue("MultiProc.UseTreeCache", 1);
   if (uc != 1) fUseTreeCache = kFALSE;
   fCacheSize = gEnv->GetValue("MultiProc.CacheSize", -1);
}

//////////////////////////////////////////////////////////////////////////
/// Handle file closing.

void TMPWorkerTree::CloseFile()
{
   // Avoid destroying the cache; must be placed before deleting the trees
   if (fFile) {
      if (fTree) fFile->SetCacheRead(0, fTree);
      delete fFile ;
      fFile = 0;
   }
}

//////////////////////////////////////////////////////////////////////////
/// Handle file opening.

TFile *TMPWorkerTree::OpenFile(const std::string& fileName)
{

   TFile *fp = TFile::Open(fileName.c_str());
   if (fp == nullptr || fp->IsZombie()) {
      std::stringstream ss;
      ss << "could not open file " << fileName;
      std::string errmsg = ss.str();
      SendError(errmsg, MPCode::kProcError);
      return nullptr;
   }

   return fp;
}

//////////////////////////////////////////////////////////////////////////
/// Retrieve a tree from an open file.

TTree *TMPWorkerTree::RetrieveTree(TFile *fp)
{
   //retrieve the TTree with the specified name from file
   //we are not the owner of the TTree object, the file is!
   TTree *tree = nullptr;
   if(fTreeName == "") {
      // retrieve the first TTree
      // (re-adapted from TEventIter.cxx)
      if (fp->GetListOfKeys()) {
         for(auto k : *fp->GetListOfKeys()) {
            TKey *key = static_cast<TKey*>(k);
            if (!strcmp(key->GetClassName(), "TTree") || !strcmp(key->GetClassName(), "TNtuple"))
               tree = static_cast<TTree*>(fp->Get(key->GetName()));
         }
      }
   } else {
      tree = static_cast<TTree*>(fp->Get(fTreeName.c_str()));
   }
   if (tree == nullptr) {
      std::stringstream ss;
      ss << "cannot find tree with name " << fTreeName << " in file " << fp->GetName();
      std::string errmsg = ss.str();
      SendError(errmsg, MPCode::kProcError);
      return nullptr;
   }

   return tree;
}

//////////////////////////////////////////////////////////////////////////
/// Tree cache handling

void TMPWorkerTree::SetupTreeCache(TTree *tree)
{
   if (fUseTreeCache) {
      TFile *curfile = tree->GetCurrentFile();
      if (curfile) {
         if (!fTreeCache) {
            tree->SetCacheSize(fCacheSize);
            fTreeCache = (TTreeCache *)curfile->GetCacheRead(tree);
            if (fCacheSize < 0) fCacheSize = tree->GetCacheSize();
         } else {
            fTreeCache->UpdateBranches(tree);
            fTreeCache->ResetCache();
            curfile->SetCacheRead(fTreeCache, tree);
         }
         if (fTreeCache) {
            fTreeCacheIsLearning = fTreeCache->IsLearning();
         }
      } else {
         Warning("SetupTreeCache", "default tree does not have a file attached: corruption? Tree cache untouched");
      }
   } else {
      // Disable the cache
      tree->SetCacheSize(0);
   }
}

//////////////////////////////////////////////////////////////////////////
/// Init overload definign max entries

void TMPWorkerTree::Init(Int_t fd, UInt_t workerN)
{

   TMPWorker::Init(fd, workerN);
   fMaxNEntries = EvalMaxEntries(fMaxNEntries);
}

//////////////////////////////////////////////////////////////////////////
/// Max entries evaluation

ULong64_t TMPWorkerTree::EvalMaxEntries(ULong64_t maxEntries)
{
   // E.g.: when dividing 10 entries between 3 workers, the first
   //       two will process 10/3 == 3 entries, the last one will process
   //       10 - 2*(10/3) == 4 entries.
   if(GetNWorker() < fNWorkers-1)
      return maxEntries/fNWorkers;
   else
      return maxEntries - (fNWorkers-1)*(maxEntries/fNWorkers);
}

//////////////////////////////////////////////////////////////////////////
/// Generic input handling

void TMPWorkerTree::HandleInput(MPCodeBufPair& msg)
{
   UInt_t code = msg.first;

   if (code == MPCode::kProcRange
         || code == MPCode::kProcFile
         || code == MPCode::kProcTree) {
      //execute fProcFunc on a file or a range of entries in a file
      Process(code, msg);
   } else if (code == MPCode::kSendResult) {
      //send back result
      SendResult();
   } else {
      //unknown code received
      std::string reply = "S" + std::to_string(GetNWorker());
      reply += ": unknown code received: " + std::to_string(code);
      MPSend(GetSocket(), MPCode::kError, reply.c_str());
   }
}



//////////////////////////////////////////////////////////////////////////
/// Selector processing SendResult and Process overload

void TMPWorkerTreeSel::SendResult()
{
   //send back result
   fSelector.SlaveTerminate();
   MPSend(GetSocket(), MPCode::kProcResult, fSelector.GetOutputList());
}

/// Selector specialization
void TMPWorkerTreeSel::Process(UInt_t code, MPCodeBufPair &msg)
{
   //evaluate the index of the file to process in fFileNames
   //(we actually don't need the parameter if code == kProcTree)

   Long64_t start = 0;
   Long64_t finish = 0;
   TEntryList *enl = 0;
   std::string errmsg;
   if (LoadTree(code, msg, start, finish, &enl, errmsg) != 0) {
      SendError(errmsg);
      return;
   }

   if (fCallBegin) {
      fSelector.SlaveBegin(nullptr);
      fCallBegin = false;
   }

   fSelector.Init(fTree);
   fSelector.Notify();
   for (Long64_t entry = start; entry < finish; ++entry) {
      Long64_t e = (enl) ? enl->GetEntry(entry) : entry;
      fSelector.Process(e);
   }

   // update the number of processed entries
   fProcessedEntries += finish - start;

   MPSend(GetSocket(), MPCode::kIdling);

   return;
}

/// Load the requierd tree and evaluate the processing range

Int_t TMPWorkerTree::LoadTree(UInt_t code, MPCodeBufPair &msg, Long64_t &start, Long64_t &finish, TEntryList **enl,
                              std::string &errmsg)
{
   // evaluate the index of the file to process in fFileNames
   //(we actually don't need the parameter if code == kProcTree)

   start = 0;
   finish = 0;
   errmsg = "";

   UInt_t fileN = 0;
   UInt_t nProcessed = 0;
   Bool_t setupcache = true;

   std::string mgroot = "[S" + std::to_string(GetNWorker()) + "]: ";

   TTree *tree = 0;
   if (code ==  MPCode::kProcTree) {

      mgroot += "MPCode::kProcTree: ";

      // The tree must be defined at this level
      if(fTree == nullptr) {
         errmsg = mgroot + std::string("tree undefined!");
         return -1;
      }

      //retrieve the total number of entries ranges processed so far by TPool
      nProcessed = ReadBuffer<UInt_t>(msg.second.get());

      //create entries range
      //example: for 21 entries, 4 workers we want ranges 0-5, 5-10, 10-15, 15-21
      //and this worker must take the rangeN-th range
      Long64_t nEntries = fTree->GetEntries();
      UInt_t nBunch = nEntries / fNWorkers;
      UInt_t rangeN = nProcessed % fNWorkers;
      start = rangeN * nBunch;
      if (rangeN < (fNWorkers - 1)) {
         finish = (rangeN+1)*nBunch;
      } else {
         finish = nEntries;
      }

      //process tree
      tree = fTree;
      CloseFile(); // May not be needed
      if (fTree->GetCurrentFile()) {
         // We need to reopen the file locally (TODO: to understand and fix this)
         if ((fFile = TFile::Open(fTree->GetCurrentFile()->GetName())) && !fFile->IsZombie()) {
            if (!(tree = (TTree *) fFile->Get(fTree->GetName()))) {
               errmsg = mgroot + std::string("unable to retrieve tree from open file ") +
                        std::string(fTree->GetCurrentFile()->GetName());
               delete fFile;
               return -1;
            }
            fTree = tree;
         } else {
            //errors are handled inside OpenFile
            errmsg = mgroot + std::string("unable to open file ") + std::string(fTree->GetCurrentFile()->GetName());
            if (fFile && fFile->IsZombie()) delete fFile;
            return -1;
         }
      }

   } else {

      if (code == MPCode::kProcRange) {
         mgroot += "MPCode::kProcRange: ";
         //retrieve the total number of entries ranges processed so far by TPool
         nProcessed = ReadBuffer<UInt_t>(msg.second.get());
         //evaluate the file and the entries range to process
         fileN = nProcessed / fNWorkers;
      } else if (code == MPCode::kProcFile) {
         mgroot += "MPCode::kProcFile: ";
         //evaluate the file and the entries range to process
         fileN = ReadBuffer<UInt_t>(msg.second.get());
      } else {
         errmsg += "MPCode undefined!";
         return -1;
      }

      // Open the file if required
      if (fFile && strcmp(fFileNames[fileN].c_str(), fFile->GetName())) CloseFile();
      if (!fFile) {
         fFile = OpenFile(fFileNames[fileN]);
         if (fFile == nullptr) {
            // errors are handled inside OpenFile
            errmsg = mgroot + std::string("unable to open file ") + fFileNames[fileN];
            return -1;
         }
      }

      //retrieve the TTree with the specified name from file
      //we are not the owner of the TTree object, the file is!
      tree = RetrieveTree(fFile);
      if (tree == nullptr) {
         //errors are handled inside RetrieveTree
         errmsg = mgroot + std::string("unable to retrieve tree from open file ") + fFileNames[fileN];
         return -1;
      }

      // Prepare to setup the cache, if required
      setupcache = (tree != fTree) ? true : false;

      // Store as reference
      fTree = tree;

      //create entries range
      if (code == MPCode::kProcRange) {
         //example: for 21 entries, 4 workers we want ranges 0-5, 5-10, 10-15, 15-21
         //and this worker must take the rangeN-th range
         Long64_t nEntries = tree->GetEntries();
         UInt_t nBunch = nEntries / fNWorkers;
         if(nEntries % fNWorkers) nBunch++;
         UInt_t rangeN = nProcessed % fNWorkers;
         start = rangeN * nBunch;
         if(rangeN < (fNWorkers-1))
            finish = (rangeN+1)*nBunch;
         else
            finish = nEntries;
      } else {
         start = 0;
         finish = tree->GetEntries();
      }
   }

   // Setup the cache, if required
   if (setupcache) SetupTreeCache(fTree);

   // Get the entrylist, if required
   if (fEntryList && enl) {
      if ((*enl = fEntryList->GetEntryList(fTree->GetName(), TUrl(fFile->GetName()).GetFile()))) {
         // create entries range
         if (code == MPCode::kProcRange) {
            // example: for 21 entries, 4 workers we want ranges 0-5, 5-10, 10-15, 15-21
            // and this worker must take the rangeN-th range
            ULong64_t nEntries = (*enl)->GetN();
            UInt_t nBunch = nEntries / fNWorkers;
            if (nEntries % fNWorkers) nBunch++;
            UInt_t rangeN = nProcessed % fNWorkers;
            start = rangeN * nBunch;
            if (rangeN < (fNWorkers - 1))
               finish = (rangeN + 1) * nBunch;
            else
               finish = nEntries;
         } else {
            start = 0;
            finish = (*enl)->GetN();
         }
      } else {
         Warning("LoadTree", "failed to get entry list for: %s %s", fTree->GetName(), TUrl(fFile->GetName()).GetFile());
      }
   }

   //check if we are going to reach the max of entries
   //change finish accordingly
   if (fMaxNEntries)
      if (fProcessedEntries + finish - start > fMaxNEntries)
         finish = start + fMaxNEntries - fProcessedEntries;

   if (gDebug > 0 && fFile)
      Info("LoadTree", "%s %d %d file: %s %lld %lld", mgroot.c_str(), nProcessed, fileN, fFile->GetName(), start,
           finish);

   return 0;
}
