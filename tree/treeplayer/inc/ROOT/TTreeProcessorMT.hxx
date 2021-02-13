// @(#)root/thread:$Id$
// Authors: Enric Tejedor, Enrico Guiraud CERN 05/06/2018

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTreeProcessorMT
#define ROOT_TTreeProcessorMT

#include "TKey.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TEntryList.h"
#include "TTreeReader.h"
#include "TError.h"
#include "TEntryList.h"
#include "TFriendElement.h"
#include "ROOT/RMakeUnique.hxx"
#include "ROOT/TThreadedObject.hxx"
#include "ROOT/TThreadExecutor.hxx"

#include <functional>
#include <utility> // std::pair
#include <vector>

/** \class TTreeView
    \brief A helper class that encapsulates a file and a tree.

A helper class that encapsulates a TFile and a TTree, along with their names.
It is used together with TTProcessor and ROOT::TThreadedObject, so that
in the TTProcessor::Process method each thread can work on its own
<TFile,TTree> pair.

This class can also be used with a collection of file names or a TChain, in case
the tree is stored in more than one file. A view will always contain only the
current (active) tree and file objects.

A copy constructor is defined for TTreeView to work with ROOT::TThreadedObject.
The latter makes a copy of a model object every time a new thread accesses
the threaded object.
*/

namespace ROOT {
namespace Internal {
/// Names, aliases, and file names of a TTree's or TChain's friends
using NameAlias = std::pair<std::string, std::string>;
struct FriendInfo {
   /// Pairs of names and aliases of friend trees/chains
   std::vector<Internal::NameAlias> fFriendNames;
   /// Names of the files where each friend is stored. fFriendFileNames[i] is the list of files for friend with
   /// name fFriendNames[i]
   std::vector<std::vector<std::string>> fFriendFileNames;
};

class TTreeView {
private:
   std::vector<std::unique_ptr<TChain>> fFriends; ///< Friends of the tree/chain, if present
   std::unique_ptr<TEntryList> fEntryList;        ///< TEntryList for fChain, if present
   // NOTE: fFriends and fEntryList MUST come before fChain to be deleted after it, because neither friend trees nor
   // entrylists are deregistered from the main tree at destruction (ROOT-9283 tracks the issue for friends).
   std::unique_ptr<TChain> fChain;                ///< Chain on which to operate

   void MakeChain(const std::vector<std::string> &treeName, const std::vector<std::string> &fileNames,
                  const FriendInfo &friendInfo, const std::vector<Long64_t> &nEntries,
                  const std::vector<std::vector<Long64_t>> &friendEntries);

public:
   TTreeView() = default;
   // no-op, we don't want to copy the local TChains
   TTreeView(const TTreeView &) {}
   std::unique_ptr<TTreeReader> GetTreeReader(Long64_t start, Long64_t end, const std::vector<std::string> &treeName,
                                              const std::vector<std::string> &fileNames, const FriendInfo &friendInfo,
                                              const TEntryList &entryList, const std::vector<Long64_t> &nEntries,
                                              const std::vector<std::vector<Long64_t>> &friendEntries);
};
} // End of namespace Internal

class TTreeProcessorMT {
private:
   const std::vector<std::string> fFileNames; ///< Names of the files
   const std::vector<std::string> fTreeNames; ///< TTree names (always same size and ordering as fFileNames)
   /// User-defined selection of entry numbers to be processed, empty if none was provided
   TEntryList fEntryList;
   const Internal::FriendInfo fFriendInfo;
   ROOT::TThreadExecutor fPool; ///<! Thread pool for processing.

   /// Thread-local TreeViews
   // Must be declared after fPool, for IMT to be initialized first!
   ROOT::TThreadedObject<ROOT::Internal::TTreeView> fTreeView{TNumSlots{ROOT::GetThreadPoolSize()}};

   Internal::FriendInfo GetFriendInfo(TTree &tree);
   std::vector<std::string> FindTreeNames();
   static unsigned int fgMaxTasksPerFilePerWorker;

public:
   TTreeProcessorMT(std::string_view filename, std::string_view treename = "", UInt_t nThreads = 0u);
   TTreeProcessorMT(const std::vector<std::string_view> &filenames, std::string_view treename = "",
                    UInt_t nThreads = 0u);
   TTreeProcessorMT(TTree &tree, const TEntryList &entries, UInt_t nThreads = 0u);
   TTreeProcessorMT(TTree &tree, UInt_t nThreads = 0u);

   void Process(std::function<void(TTreeReader &)> func);
   static void SetMaxTasksPerFilePerWorker(unsigned int m);
   static unsigned int GetMaxTasksPerFilePerWorker();
};

} // End of namespace ROOT

#endif // defined TTreeProcessorMT
