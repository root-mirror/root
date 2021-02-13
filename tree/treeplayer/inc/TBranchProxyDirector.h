// @(#)root/base:$Id$
// Author: Philippe Canal 13/05/2003

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun, Fons Rademakers and al.           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBranchProxyDirector
#define ROOT_TBranchProxyDirector

#include "RtypesCore.h"
#include <vector>
#include <list>
#include <algorithm>

class TH1F;
class TTree;

namespace ROOT {
namespace Detail {
   class TBranchProxy;
   class TFriendProxy;
}

namespace Internal{
   class TFriendProxy;

   // Helper function to call SetReadEntry on all TFriendProxy
   void ResetReadEntry(TFriendProxy *fp);

   class TBranchProxyDirector {

      //This class could actually be the selector itself.
      TTree   *fTree;  // TTree we are currently looking at.
      Long64_t fEntry; // Entry currently being read (in the local TTree rather than the TChain)

      std::list<Detail::TBranchProxy*> fDirected;
      std::vector<TFriendProxy*> fFriends;

      TBranchProxyDirector(const TBranchProxyDirector &) : fTree(nullptr), fEntry(-1) {;}
      TBranchProxyDirector& operator=(const TBranchProxyDirector&) {return *this;}

   public:

      TBranchProxyDirector(TTree* tree, Long64_t i);
      TBranchProxyDirector(TTree* tree, Int_t i);     // cint has (had?) a problem casting int to long long

      void     Attach(Detail::TBranchProxy* p);
      void     Attach(TFriendProxy* f);
      TH1F*    CreateHistogram(const char *options);

      /// Return the current 'local' entry number; i.e. in the 'local' TTree rather than the TChain.
      /// This value will be passed directly to TBranch::GetEntry.
      Long64_t GetReadEntry() const { return fEntry; }

      TTree*   GetTree() const { return fTree; };
      // void   Print();

      /// Move to a new entry to read
      /// entry is the 'local' entry number; i.e. in the 'local' TTree rather than the TChain.
      /// This value will be passed directly to TBranch::GetEntry.
      void     SetReadEntry(Long64_t entry) {
         fEntry = entry;
         if (!fFriends.empty()) {
            std::for_each(fFriends.begin(), fFriends.end(), ResetReadEntry);
         }
      }
      TTree*   SetTree(TTree *newtree);
      Bool_t   Notify();

   };

} // namespace Internal
} // namespace ROOT

#endif
