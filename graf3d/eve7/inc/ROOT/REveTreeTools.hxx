// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveTreeTools
#define ROOT7_REveTreeTools

#include "TSelectorDraw.h"
#include "TEventList.h"

namespace ROOT {
namespace Experimental {

/////////////////////////////////////////////////////////////////////////////////
/// REveSelectorToEventList
/// TSelector that stores entry numbers of matching TTree entries into an event-list.
/////////////////////////////////////////////////////////////////////////////////

class REveSelectorToEventList : public TSelectorDraw
{
   REveSelectorToEventList(const REveSelectorToEventList &) = delete;
   REveSelectorToEventList &operator=(const REveSelectorToEventList &) = delete;

protected:
   TEventList *fEvList{nullptr};
   TList fInput;

public:
   REveSelectorToEventList(TEventList *evl, const char *sel);
   virtual ~REveSelectorToEventList() {}

   virtual Int_t Version() const { return 1; }
   virtual Bool_t Process(Long64_t entry);

   ClassDef(REveSelectorToEventList, 1); // TSelector that stores entry numbers of matching TTree entries into an event-list.
};

/////////////////////////////////////////////////////////////////////////////////
/// REvePointSelectorConsumer
/// Virtual base for classes that can be filled from TTree data via the REvePointSelector class.
/////////////////////////////////////////////////////////////////////////////////

class REvePointSelector;

class REvePointSelectorConsumer
{
public:
   enum ETreeVarType_e { kTVT_XYZ, kTVT_RPhiZ };

protected:
   ETreeVarType_e fSourceCS; // Coordinate-System of the source tree variables

public:
   REvePointSelectorConsumer(ETreeVarType_e cs = kTVT_XYZ) : fSourceCS(cs) {}
   virtual ~REvePointSelectorConsumer() {}

   virtual void InitFill(Int_t /*subIdNum*/) {}
   virtual void TakeAction(REvePointSelector *) = 0;

   ETreeVarType_e GetSourceCS() const { return fSourceCS; }
   void SetSourceCS(ETreeVarType_e cs) { fSourceCS = cs; }
};

/////////////////////////////////////////////////////////////////////////////////
/// REvePointSelector
/// TSelector for direct extraction of point-like data from a Tree.
/////////////////////////////////////////////////////////////////////////////////

class REvePointSelector : public TSelectorDraw
{
   REvePointSelector(const REvePointSelector &) = delete;
   REvePointSelector &operator=(const REvePointSelector &) = delete;

protected:
   TTree *fTree{nullptr};
   REvePointSelectorConsumer *fConsumer{nullptr};

   TString fVarexp;
   TString fSelection;

   TString fSubIdExp;
   Int_t fSubIdNum;

   TList fInput;

public:
   REvePointSelector(TTree *t = nullptr, REvePointSelectorConsumer *c = nullptr, const char *vexp = "", const char *sel = "");
   virtual ~REvePointSelector() {}

   virtual Long64_t Select(const char *selection = nullptr);
   virtual Long64_t Select(TTree *t, const char *selection = nullptr);
   virtual void TakeAction();

   TTree *GetTree() const { return fTree; }
   void SetTree(TTree *t) { fTree = t; }

   REvePointSelectorConsumer *GetConsumer() const { return fConsumer; }
   void SetConsumer(REvePointSelectorConsumer *c) { fConsumer = c; }

   const char *GetVarexp() const { return fVarexp; }
   void SetVarexp(const char *v) { fVarexp = v; }

   const char *GetSelection() const { return fSelection; }
   void SetSelection(const char *s) { fSelection = s; }

   const char *GetSubIdExp() const { return fSubIdExp; }
   void SetSubIdExp(const char *s) { fSubIdExp = s; }

   Int_t GetSubIdNum() const { return fSubIdNum; }

   ClassDef(REvePointSelector, 1); // TSelector for direct extraction of point-like data from a Tree.
};

} // namespace Experimental
} // namespace ROOT

#endif
