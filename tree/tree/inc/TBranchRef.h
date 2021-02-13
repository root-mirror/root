// @(#)root/tree:$Id$
// Author: Rene Brun   19/08/2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBranchRef
#define ROOT_TBranchRef


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TBranchRef                                                           //
//                                                                      //
// A Branch to support referenced objects on other branches             //
//////////////////////////////////////////////////////////////////////////


#include "TBranch.h"

#ifdef R__LESS_INCLUDES
class TRefTable;
#else
#include "TRefTable.h"
#endif

class TTree;

class TBranchRef : public TBranch {
private:
   Long64_t   fRequestedEntry;  ///<! Cursor indicating which entry is being requested.

protected:
   TRefTable *fRefTable;        ///< pointer to the TRefTable

   void    ReadLeavesImpl(TBuffer &b);
   void    FillLeavesImpl(TBuffer &b);

public:
   TBranchRef();
   TBranchRef(TTree *tree);
   virtual ~TBranchRef();
   virtual void    Clear(Option_t *option="");
   TRefTable      *GetRefTable() const {return fRefTable;}
   virtual Bool_t  Notify();
   virtual void    Print(Option_t *option="") const;
   virtual void    Reset(Option_t *option="");
   virtual void    ResetAfterMerge(TFileMergeInfo *);
   virtual Int_t   SetParent(const TObject* obj, Int_t branchID);
   virtual void    SetRequestedEntry(Long64_t entry) {fRequestedEntry = entry;}

private:
   virtual Int_t   FillImpl(ROOT::Internal::TBranchIMTHelper *);

   ClassDef(TBranchRef,1);  //to support referenced objects on other branches
};

#endif
