// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//______________________________________________________________________________
// TTreeTools
//
// Collection of classes for TTree interaction.

#include <ROOT/REveTreeTools.hxx>

#include "TTree.h"
#include "TTreeFormula.h"

using namespace ROOT::Experimental;

/** \class REveSelectorToEventList
\ingroup REve
TSelector that stores entry numbers of matching TTree entries into
an event-list.
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveSelectorToEventList::REveSelectorToEventList(TEventList* evl, const char* sel) :
   TSelectorDraw(), fEvList(evl)
{
   fInput.Add(new TNamed("varexp", ""));
   fInput.Add(new TNamed("selection", sel));
   SetInputList(&fInput);
}

////////////////////////////////////////////////////////////////////////////////
/// Process entry.

Bool_t REveSelectorToEventList::Process(Long64_t entry)
{
   if(GetSelect()->EvalInstance(0) != 0)
      fEvList->Enter(entry);
   return kTRUE;
}

/** \class REvePointSelectorConsumer
\ingroup REve
REvePointSelectorConsumer is a virtual base for classes that can be
filled from TTree data via the REvePointSelector class.
*/

/** \class REvePointSelector
\ingroup REve
REvePointSelector is a sub-class of TSelectorDraw for direct
extraction of point-like data from a Tree.
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REvePointSelector::REvePointSelector(TTree* t,
                                     REvePointSelectorConsumer* c,
                                     const char* vexp, const char* sel) :
   TSelectorDraw(),

   fTree      (t),
   fConsumer  (c),
   fVarexp    (vexp),
   fSelection (sel),
   fSubIdExp  (),
   fSubIdNum  (0)
{
   SetInputList(&fInput);
}

////////////////////////////////////////////////////////////////////////////////
/// Process the tree, select points matching 'selection'.

Long64_t REvePointSelector::Select(const char* selection)
{
   TString var(fVarexp);
   if (fSubIdExp.IsNull()) {
      fSubIdNum = 0;
   } else {
      fSubIdNum = fSubIdExp.CountChar(':') + 1;
      var += ":" + fSubIdExp;
   }

   TString sel;
   if (selection)
      sel = selection;
   else
      sel = fSelection;

   fInput.Delete();
   fInput.Add(new TNamed("varexp",    var.Data()));
   fInput.Add(new TNamed("selection", sel.Data()));

   if (fConsumer)
      fConsumer->InitFill(fSubIdNum);

   if (fTree)
      fTree->Process(this, "goff");

   return fSelectedRows;
}

////////////////////////////////////////////////////////////////////////////////
/// Process tree 't', select points matching 'selection'.

Long64_t REvePointSelector::Select(TTree *t, const char *selection)
{
   fTree = t;
   return Select(selection);
}

////////////////////////////////////////////////////////////////////////////////
/// Callback from tree-player after a chunk of data has been processed.
/// This is forwarded to the current point-consumer.

void REvePointSelector::TakeAction()
{
   fSelectedRows += fNfill;
   // printf("REvePointSelector::TakeAction nfill=%d, nall=%lld\n", fNfill, fSelectedRows);
   if (fConsumer) {
      fConsumer->TakeAction(this);
   }
}
