/* @(#)root/multiproc:$Id$ */
// Author: G Ganis Jan 2017

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMPWorkerTree
#define ROOT_TMPWorkerTree

#include "TMPWorker.h"
#include "TFile.h"
#include "TEntryList.h"
#include "TEventList.h"
#include "TH1.h"
#include "TKey.h"
#include "TSelector.h"
#include "TTree.h"
#include "TTreeCache.h"
#include "TTreeReader.h"

#include <memory> //unique_ptr
#include <string>
#include <sstream>
#include <type_traits> //std::result_of
#include <unistd.h> //pid_t
#include <vector>

class TMPWorkerTree : public TMPWorker {
   /// \cond
//   ClassDef(TMPWorkerTree, 0);
   /// \endcond
public:
   TMPWorkerTree();
   TMPWorkerTree(const std::vector<std::string> &fileNames, TEntryList *entries, const std::string &treeName,
                 UInt_t nWorkers, ULong64_t maxEntries, ULong64_t firstEntry);
   TMPWorkerTree(TTree *tree, TEntryList *entries, UInt_t nWorkers, ULong64_t maxEntries, ULong64_t firstEntry);
   virtual ~TMPWorkerTree();

   // It doesn't make sense to copy a TMPWorker (each one has a uniq_ptr to its socket)
   TMPWorkerTree(const TMPWorkerTree &) = delete;
   TMPWorkerTree &operator=(const TMPWorkerTree &) = delete;

protected:

   void         CloseFile();
   ULong64_t    EvalMaxEntries(ULong64_t maxEntries);
   void         HandleInput(MPCodeBufPair& msg); ///< Execute instructions received from a MP client
   void Init(int fd, UInt_t workerN);
   Int_t LoadTree(UInt_t code, MPCodeBufPair &msg, Long64_t &start, Long64_t &finish, TEntryList **enl,
                  std::string &errmsg);
   TFile       *OpenFile(const std::string& fileName);
   virtual void Process(UInt_t, MPCodeBufPair &) {}
   TTree       *RetrieveTree(TFile *fp);
   virtual void SendResult() { }
   void         Setup();
   void         SetupTreeCache(TTree *tree);

   std::vector<std::string> fFileNames; ///< the files to be processed by all workers
   std::string fTreeName;               ///< the name of the tree to be processed
   TTree *fTree;                        ///< pointer to the tree to be processed. It is only used if the tree is directly passed to TProcessExecutor::Process as argument
   TFile *fFile;                        ///< last open file
   TEntryList *fEntryList;              ///< entrylist
   ULong64_t fFirstEntry;               ///< first entry to br processed

private:

   // TTree cache handling
   TTreeCache *fTreeCache;              ///< instance of the tree cache for the tree
   Bool_t      fTreeCacheIsLearning;    ///< Whether cache is in learning phase
   Bool_t      fUseTreeCache;           ///< Control usage of the tree cache
   Long64_t    fCacheSize;              ///< Cache size
};

template<class F>
class TMPWorkerTreeFunc : public TMPWorkerTree {
public:
   TMPWorkerTreeFunc(F procFunc, const std::vector<std::string> &fileNames, TEntryList *entries,
                     const std::string &treeName, UInt_t nWorkers, ULong64_t maxEntries, ULong64_t firstEntry)
      : TMPWorkerTree(fileNames, entries, treeName, nWorkers, maxEntries, firstEntry), fProcFunc(procFunc),
        fReducedResult(), fCanReduce(false)
   {
   }
   TMPWorkerTreeFunc(F procFunc, TTree *tree, TEntryList *entries, UInt_t nWorkers, ULong64_t maxEntries,
                     ULong64_t firstEntry)
      : TMPWorkerTree(tree, entries, nWorkers, maxEntries, firstEntry), fProcFunc(procFunc), fReducedResult(),
        fCanReduce(false)
   {
   }
   virtual ~TMPWorkerTreeFunc() {}

private:
   void Process(UInt_t code, MPCodeBufPair &msg);
   void SendResult();

   F  fProcFunc; ///< copy the function to be executed
   typename std::result_of<F(std::reference_wrapper<TTreeReader>)>::type fReducedResult; ///< the results of the executions of fProcFunc merged together
   bool fCanReduce; ///< true if fReducedResult can be reduced with a new result, false until we have produced one result
};

class TMPWorkerTreeSel : public TMPWorkerTree {
public:
   TMPWorkerTreeSel(TSelector &selector, const std::vector<std::string> &fileNames, TEntryList *entries,
                    const std::string &treeName, UInt_t nWorkers, ULong64_t maxEntries, ULong64_t firstEntry)
      : TMPWorkerTree(fileNames, entries, treeName, nWorkers, maxEntries, firstEntry), fSelector(selector),
        fCallBegin(true)
   {
   }
   TMPWorkerTreeSel(TSelector &selector, TTree *tree, TEntryList *entries, UInt_t nWorkers, ULong64_t maxEntries,
                    ULong64_t firstEntry)
      : TMPWorkerTree(tree, entries, nWorkers, maxEntries, firstEntry), fSelector(selector), fCallBegin(true)
   {
   }
   virtual ~TMPWorkerTreeSel() {}

private:
   void Process(UInt_t code, MPCodeBufPair &msg);
   void SendResult();

   TSelector &fSelector; ///< pointer to the selector to be used to process the tree. It is null if we are not using a TSelector.
   bool fCallBegin = true;
};

//////////////////////////////////////////////////////////////////////////
/// Auxilliary templated functions
/// If the user lambda returns a TH1F*, TTree*, TEventList*, we incur in the
/// problem of that object being automatically owned by the current open file.
/// For these three types, we call SetDirectory(nullptr) to detach the returned
/// object from the file we are reading the TTree from.
/// Note: the only sane case in which this should happen is when a TH1F* is
/// returned.
template <class T, typename std::enable_if<std::is_pointer<T>::value && std::is_constructible<TObject *, T>::value &&
                                           !std::is_constructible<TCollection *, T>::value>::type * = nullptr>
void DetachRes(T res)
{
   auto th1p = dynamic_cast<TH1*>(res);
   if(th1p != nullptr) {
      th1p->SetDirectory(nullptr);
      return;
   }
   auto ttreep = dynamic_cast<TTree*>(res);
   if(ttreep != nullptr) {
      ttreep->SetDirectory(nullptr);
      return;
   }
   auto tentrylist = dynamic_cast<TEntryList*>(res);
   if(tentrylist != nullptr) {
      tentrylist->SetDirectory(nullptr);
      return;
   }
   auto teventlist = dynamic_cast<TEventList*>(res);
   if(teventlist != nullptr) {
      teventlist->SetDirectory(nullptr);
      return;
   }
   return;
}

// Specialization for TCollections
template <class T, typename std::enable_if<std::is_pointer<T>::value &&
                                           std::is_constructible<TCollection *, T>::value>::type * = nullptr>
void DetachRes(T res)
{
   if (res) {
      TIter nxo(res);
      TObject *obj = 0;
      while ((obj = nxo())) {
         DetachRes(obj);
      }
   }
}

//////////////////////////////////////////////////////////////////////////
/// Generic function processing SendResult and Process overload

template<class F>
void TMPWorkerTreeFunc<F>::SendResult()
{
   //send back result
   MPSend(GetSocket(), MPCode::kProcResult, fReducedResult);
}

template <class F>
void TMPWorkerTreeFunc<F>::Process(UInt_t code, MPCodeBufPair &msg)
{

   Long64_t start = 0;
   Long64_t finish = 0;
   TEntryList *enl = 0;
   std::string reply, errmsg, sn = "[S" + std::to_string(GetNWorker()) + "]: ";
   if (LoadTree(code, msg, start, finish, &enl, errmsg) != 0) {
      reply = sn + errmsg;
      MPSend(GetSocket(), MPCode::kProcError, reply.c_str());
      return;
   }

   // create a TTreeReader that reads this range of entries
   TTreeReader reader(fTree, enl);

   TTreeReader::EEntryStatus status = reader.SetEntriesRange(start, finish);
   if(status != TTreeReader::kEntryValid) {
      reply = sn + "could not set TTreeReader to range " + std::to_string(start) + " " + std::to_string(finish - 1);
      MPSend(GetSocket(), MPCode::kProcError, reply.c_str());
      return;
   }

   //execute function
   auto res = fProcFunc(reader);

   //detach result from file if needed (currently needed for TH1, TTree, TEventList)
   DetachRes(res);

   //update the number of processed entries
   fProcessedEntries += finish - start;

   if(fCanReduce) {
      PoolUtils::ReduceObjects<TObject *> redfunc;
      fReducedResult = static_cast<decltype(fReducedResult)>(redfunc({res, fReducedResult})); //TODO try not to copy these into a vector, do everything by ref. std::vector<T&>?
   } else {
      fCanReduce = true;
      fReducedResult = res;
   }

   if(fMaxNEntries == fProcessedEntries)
      //we are done forever
      MPSend(GetSocket(), MPCode::kProcResult, fReducedResult);
   else
      //we are done for now
      MPSend(GetSocket(), MPCode::kIdling);
}

#endif
