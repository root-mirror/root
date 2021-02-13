// @(#)root/base:$Id$
// Author: Philippe Canal  13/05/2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun, Fons Rademakers and al.           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** TBranchProxyDirector
This class is used to 'drive' and hold a serie of TBranchProxy objects
which represent and give access to the content of TTree object.
This is intended to be used as part of a generate Selector class
which will hold the directory and its associate
*/

#include "TBranchProxyDirector.h"
#include "TBranchProxy.h"
#include "TFriendProxy.h"
#include "TTree.h"
#include "TEnv.h"
#include "TH1F.h"
#include "TVirtualPad.h"
#include "TList.h"

#include <algorithm>

ClassImp(ROOT::Internal::TBranchProxyDirector);

namespace ROOT {
namespace Internal {

   // Helper function to call Reset on each TBranchProxy
   void NotifyDirected(Detail::TBranchProxy *x) { x->Notify(); }

   // Helper function to call SetReadEntry on all TFriendProxy
   void ResetReadEntry(TFriendProxy *fp) { fp->ResetReadEntry(); }

   // Helper class to call Update on all TFriendProxy
   struct Update {
      Update(TTree *newtree) : fNewTree(newtree) {}
      TTree *fNewTree;
      void operator()(TFriendProxy *x) { x->Update(fNewTree); }
   };


   TBranchProxyDirector::TBranchProxyDirector(TTree* tree, Long64_t i) :
      fTree(tree),
      fEntry(i)
   {
      // Simple constructor
   }

   TBranchProxyDirector::TBranchProxyDirector(TTree* tree, Int_t i) :
      // cint has a problem casting int to long long
      fTree(tree),
      fEntry(i)
   {
      // Simple constructor
   }

   void TBranchProxyDirector::Attach(Detail::TBranchProxy* p) {

      // Attach a TBranchProxy object to this director.  The director just
      // 'remembers' this BranchProxy and does not own it.  It will be use
      // to apply Tree wide operation (like reseting).
      fDirected.push_back(p);
   }

   void TBranchProxyDirector::Attach(TFriendProxy* p) {

      // Attach a TFriendProxy object to this director.  The director just
      // 'remembers' this BranchProxy and does not own it.  It will be use
      // to apply Tree wide operation (like reseting).
      fFriends.push_back(p);
   }

   TH1F* TBranchProxyDirector::CreateHistogram(const char *options) {
      // Create a temporary 1D histogram.

      Int_t nbins = gEnv->GetValue("Hist.Binning.1D.x",100);
      Double_t vmin=0, vmax=0;
      Double_t xmin=0, xmax=0;
      Bool_t canExtend = kTRUE;
      TString opt( options );
      Bool_t optSame = opt.Contains("same");
      if (optSame) canExtend = kFALSE;

      if (gPad && optSame) {
         TListIter np(gPad->GetListOfPrimitives());
         TObject *op;
         TH1 *oldhtemp = 0;
         while ((op = np()) && !oldhtemp) {
            if (op->InheritsFrom(TH1::Class())) oldhtemp = (TH1 *)op;
         }
         if (oldhtemp) {
            nbins = oldhtemp->GetXaxis()->GetNbins();
            vmin = oldhtemp->GetXaxis()->GetXmin();
            vmax = oldhtemp->GetXaxis()->GetXmax();
         } else {
            vmin = gPad->GetUxmin();
            vmax = gPad->GetUxmax();
         }
      } else {
         vmin = xmin;
         vmax = xmax;
         if (xmin < xmax) canExtend = kFALSE;
      }
      TH1F *hist = new TH1F("htemp","htemp",nbins,vmin,vmax);
      hist->SetLineColor(fTree->GetLineColor());
      hist->SetLineWidth(fTree->GetLineWidth());
      hist->SetLineStyle(fTree->GetLineStyle());
      hist->SetFillColor(fTree->GetFillColor());
      hist->SetFillStyle(fTree->GetFillStyle());
      hist->SetMarkerStyle(fTree->GetMarkerStyle());
      hist->SetMarkerColor(fTree->GetMarkerColor());
      hist->SetMarkerSize(fTree->GetMarkerSize());
      if (canExtend) hist->SetCanExtend(TH1::kAllAxes);
      hist->GetXaxis()->SetTitle("var");
      hist->SetBit(kCanDelete);
      hist->SetDirectory(0);

      if (opt.Length() && opt.Contains("e")) hist->Sumw2();
      return hist;
   }

   TTree* TBranchProxyDirector::SetTree(TTree *newtree) {

      // Set the BranchProxy to be looking at a new tree.
      // Reset all.
      // Return the old tree.

      TTree* oldtree = fTree;
      fTree = newtree;
      if(!Notify()) return nullptr;
      return oldtree;
   }

   Bool_t TBranchProxyDirector::Notify() {
      fEntry = -1;
      bool retVal = true;
      for_each(fDirected.begin(),fDirected.end(),NotifyDirected);
      for (auto brProxy : fDirected) {
         retVal = retVal && brProxy->Notify();
      }
      Update update(fTree);
      for_each(fFriends.begin(),fFriends.end(),update);
      return retVal;
   }

} // namespace Internal
} // namespace ROOT
