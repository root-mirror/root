// @(#)root/tree:$Id$
// Author: Rene Brun   30/08/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TVirtualTreePlayer
#define ROOT_TVirtualTreePlayer


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TVirtualTreePlayer                                                   //
//                                                                      //
// Abstract base class for Tree players.                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"

#include <limits>

class TClass;
class TTree;
class TH1;
class TTreeFormula;
class TSQLResult;
class TSelector;
class TPrincipal;
class TVirtualIndex;


class TVirtualTreePlayer : public TObject {

private:
   static TClass              *fgPlayer;  ///< Pointer to class of Tree player
   static TVirtualTreePlayer  *fgCurrent; ///< Pointer to current Tree player

public:
   static constexpr Long64_t kMaxEntries = std::numeric_limits<Long64_t>::max();

   TVirtualTreePlayer() { }
   virtual ~TVirtualTreePlayer();
   virtual TVirtualIndex *BuildIndex(const TTree *T, const char *majorname, const char *minorname) = 0;
   virtual TTree         *CopyTree(const char *selection, Option_t *option=""
                                   ,Long64_t nentries=kMaxEntries, Long64_t firstentry=0) = 0;
   virtual Long64_t       DrawScript(const char *wrapperPrefix,
                                     const char *macrofilename, const char *cutfilename,
                                     Option_t *option, Long64_t nentries, Long64_t firstentry) = 0;
   virtual Long64_t       DrawSelect(const char *varexp, const char *selection, Option_t *option
                                     ,Long64_t nentries, Long64_t firstentry) = 0;
   virtual Int_t          Fit(const char *formula ,const char *varexp, const char *selection,Option_t *option ,Option_t *goption
                              ,Long64_t nentries, Long64_t firstentry) = 0;
   virtual Int_t          GetDimension() const = 0;
   virtual TH1           *GetHistogram() const = 0;
   virtual Int_t          GetNfill() const = 0;
   virtual Long64_t       GetEntries(const char *) = 0;
   virtual Long64_t       GetSelectedRows() const = 0;
   virtual TSelector     *GetSelector() const = 0;
   virtual TSelector     *GetSelectorFromFile() const = 0;
   virtual TTreeFormula  *GetSelect() const = 0;
   virtual TTreeFormula  *GetVar(Int_t) const = 0;
   virtual TTreeFormula  *GetVar1() const = 0;
   virtual TTreeFormula  *GetVar2() const = 0;
   virtual TTreeFormula  *GetVar3() const = 0;
   virtual TTreeFormula  *GetVar4() const = 0;
   virtual Double_t      *GetVal(Int_t) const = 0;
   virtual Double_t      *GetV1() const = 0;
   virtual Double_t      *GetV2() const = 0;
   virtual Double_t      *GetV3() const = 0;
   virtual Double_t      *GetV4() const = 0;
   virtual Double_t      *GetW() const = 0;
   virtual Int_t          MakeClass(const char *classname, const char *option) = 0;
   virtual Int_t          MakeCode(const char *filename) = 0;
   virtual Int_t          MakeProxy(const char *classname,
                                    const char *macrofilename = 0, const char *cutfilename = 0,
                                    const char *option = 0, Int_t maxUnrolling = 3) = 0;
   virtual Int_t          MakeReader(const char *classname, Option_t *option) = 0;
   virtual TPrincipal    *Principal(const char *varexp="", const char *selection="", Option_t *option="np"
                                    ,Long64_t nentries=kMaxEntries, Long64_t firstentry=0) = 0;
   virtual Long64_t       Process(const char *filename,Option_t *option="", Long64_t nentries=kMaxEntries, Long64_t firstentry=0) = 0;
   virtual Long64_t       Process(TSelector *selector,Option_t *option="",  Long64_t nentries=kMaxEntries, Long64_t firstentry=0) = 0;
   virtual Long64_t       Scan(const char *varexp, const char *selection, Option_t *option
                               ,Long64_t nentries, Long64_t firstentry) = 0;
   virtual TSQLResult    *Query(const char *varexp, const char *selection, Option_t *option
                                ,Long64_t nentries, Long64_t firstentry) = 0;
   virtual void           SetEstimate(Long64_t n) = 0;
   virtual void           SetTree(TTree *t) = 0;
   virtual void           StartViewer(Int_t ww, Int_t wh) = 0;
   virtual Int_t          UnbinnedFit(const char *formula ,const char *varexp, const char *selection,Option_t *option
                                      ,Long64_t nentries, Long64_t firstentry) = 0;
   virtual void           UpdateFormulaLeaves() = 0;

   static  TVirtualTreePlayer *GetCurrentPlayer();
   static  TVirtualTreePlayer *TreePlayer(TTree *obj);
   static void        SetPlayer(const char *player);

   ClassDef(TVirtualTreePlayer,0);  //Abstract interface for Tree players
};

#endif
