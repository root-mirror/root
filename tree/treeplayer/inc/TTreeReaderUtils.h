// @(#)root/tree:$Id$
// Author: Axel Naumann, 2010-10-12

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTreeReaderUtils
#define ROOT_TTreeReaderUtils


////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TTreeReaderUtils                                                       //
//                                                                        //
// TTreeReader's helpers.                                                 //
//                                                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TBranchProxyDirector.h"
#include "TBranchProxy.h"
#include "TTreeReaderValue.h"

#include <string>

class TDictionary;
class TTree;

namespace ROOT {
   namespace Detail {
      class TBranchProxy;
   }

namespace Internal {
   class TBranchProxyDirector;
   class TTreeReaderArrayBase;

   class TNamedBranchProxy {
   public:
      TNamedBranchProxy(): fDict(0), fContentDict(0) {}
      TNamedBranchProxy(TBranchProxyDirector* boss, TBranch* branch, const char* fullname, const char* membername):
         fProxy(boss, fullname, branch, membername), fDict(0), fContentDict(0), fFullName(fullname) {}

      // Constructor for friend case, the fullname (containing the name of the friend tree) may be different
      // from the lookup name (without the name of the friend)
      TNamedBranchProxy(TBranchProxyDirector* boss, TBranch* branch, const char* fullname, const char* proxyname, const char* membername):
         fProxy(boss, proxyname, branch, membername), fDict(0), fContentDict(0), fFullName(fullname) {}

      const char* GetName() const { return fFullName.c_str(); }
      const Detail::TBranchProxy* GetProxy() const { return &fProxy; }
      Detail::TBranchProxy* GetProxy() { return &fProxy; }
      TDictionary* GetDict() const { return fDict; }
      void SetDict(TDictionary* dict) { fDict = dict; }
      TDictionary* GetContentDict() const { return fContentDict; }
      void SetContentDict(TDictionary* dict) { fContentDict = dict; }

   private:
      Detail::TBranchProxy fProxy;
      TDictionary*         fDict;
      TDictionary*         fContentDict; // type of content, if a collection
      std::string          fFullName;
   };

   // Used by TTreeReaderArray
   class TVirtualCollectionReader {
   public:
      TTreeReaderValueBase::EReadStatus fReadStatus;

      TVirtualCollectionReader() : fReadStatus(TTreeReaderValueBase::kReadNothingYet) {}

      virtual ~TVirtualCollectionReader();
      virtual size_t GetSize(Detail::TBranchProxy*) = 0;
      virtual void* At(Detail::TBranchProxy*, size_t /*idx*/) = 0;
   };

}
}

#endif // defined TTreeReaderUtils
