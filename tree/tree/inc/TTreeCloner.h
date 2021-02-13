// @(#)root/tree:$Id$
// Author: Philippe Canal 07/11/2005

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTreeCloner
#define ROOT_TTreeCloner

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TTreeCloner                                                          //
//                                                                      //
// Class implementing or helping  the various TTree cloning method      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObjArray.h"

class TBranch;
class TTree;
class TFileCacheRead;

class TTreeCloner {
   TString    fWarningMsg;       ///< Text of the error message lead to an 'invalid' state

   Bool_t     fIsValid;
   Bool_t     fNeedConversion;   ///< True if the fast merge is not possible but a slow merge might possible.
   UInt_t     fOptions;
   TTree     *fFromTree;
   TTree     *fToTree;
   Option_t  *fMethod;
   TObjArray  fFromBranches;
   TObjArray  fToBranches;

   UInt_t     fMaxBaskets;
   UInt_t    *fBasketBranchNum;  ///<[fMaxBaskets] Index of the branch(es) of the basket.
   UInt_t    *fBasketNum;        ///<[fMaxBaskets] index of the basket within the branch.

   Long64_t  *fBasketSeek;       ///<[fMaxBaskets] list of basket position to be read.
   Long64_t  *fBasketEntry;      ///<[fMaxBaskets] list of basket start entries.
   UInt_t    *fBasketIndex;      ///<[fMaxBaskets] ordered list of basket indices to be written.

   UShort_t   fPidOffset;        ///< Offset to be added to the copied key/basket.

   UInt_t     fCloneMethod;      ///< Indicates which cloning method was selected.
   Long64_t   fToStartEntries;   ///< Number of entries in the target tree before any addition.

   Int_t           fCacheSize;   ///< Requested size of the file cache
   TFileCacheRead *fFileCache;   ///< File Cache used to reduce the number of individual reads
   TFileCacheRead *fPrevCache;   ///< Cache that set before the TTreeCloner ctor for the 'from' TTree if any.

   enum ECloneMethod {
      kDefault             = 0,
      kSortBasketsByBranch = 1,
      kSortBasketsByOffset = 2,
      kSortBasketsByEntry  = 3
   };

   class CompareSeek {
      TTreeCloner *fObject;
   public:
      CompareSeek(TTreeCloner *obj) : fObject(obj) {}
      Bool_t operator()(UInt_t i1, UInt_t i2);
   };

   class CompareEntry {
      TTreeCloner *fObject;
   public:
      CompareEntry(TTreeCloner *obj) : fObject(obj) {}
      Bool_t operator()(UInt_t i1, UInt_t i2);
   };

   friend class CompareSeek;
   friend class CompareEntry;

   void ImportClusterRanges();
   void CreateCache();
   UInt_t FillCache(UInt_t from);
   void RestoreCache();

private:
   TTreeCloner(const TTreeCloner&) = delete;
   TTreeCloner &operator=(const TTreeCloner&) = delete;

public:
   enum EClonerOptions {
      kNone       = 0,
      kNoWarnings = BIT(1),
      kIgnoreMissingTopLevel = BIT(2),
      kNoFileCache = BIT(3)
   };

   TTreeCloner(TTree *from, TTree *to, Option_t *method, UInt_t options = kNone);
   virtual ~TTreeCloner();

   void   CloseOutWriteBaskets();
   UInt_t CollectBranches(TBranch *from, TBranch *to);
   UInt_t CollectBranches(TObjArray *from, TObjArray *to);
   UInt_t CollectBranches();
   void   CollectBaskets();
   void   CopyMemoryBaskets();
   void   CopyStreamerInfos();
   void   CopyProcessIds();
   const char *GetWarning() const { return fWarningMsg; }
   Bool_t Exec();
   Bool_t IsValid() { return fIsValid; }
   Bool_t NeedConversion() { return fNeedConversion; }
   void   SetCacheSize(Int_t size);
   void   SortBaskets();
   void   WriteBaskets();

   ClassDef(TTreeCloner,0); // helper used for the fast cloning of TTrees.
};

#endif
