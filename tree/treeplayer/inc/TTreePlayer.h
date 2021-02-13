// @(#)root/treeplayer:$Id$
// Author: Rene Brun   12/01/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTreePlayer
#define ROOT_TTreePlayer


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TTreePlayer                                                          //
//                                                                      //
// A TTree object is a list of TBranch.                                 //
//   To Create a TTree object one must:                                 //
//    - Create the TTree header via the TTree constructor               //
//    - Call the TBranch constructor for every branch.                  //
//                                                                      //
//   To Fill this object, use member function Fill with no parameters.  //
//     The Fill function loops on all defined TBranch.                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TVirtualTreePlayer.h"

#include "TSelectorDraw.h"
#include "TTree.h"

class TVirtualIndex;

class TTreePlayer : public TVirtualTreePlayer {

private:
   TTreePlayer(const TTreePlayer &) = delete;
   TTreePlayer& operator=(const TTreePlayer &) = delete;

protected:
   TTree         *fTree;            //!  Pointer to current Tree
   Bool_t         fScanRedirect;    //  Switch to redirect TTree::Scan output to a file
   const char    *fScanFileName;    //  Name of the file where Scan is redirected
   Int_t          fDimension;       //  Dimension of the current expression
   Long64_t       fSelectedRows;    //  Number of selected entries
   TH1           *fHistogram;       //! Pointer to histogram used for the projection
   TSelectorDraw *fSelector;        //! Pointer to current selector
   TSelector     *fSelectorFromFile;//! Pointer to a user defined selector created by this TTreePlayer object
   TClass        *fSelectorClass;   //! Pointer to the actual class of the TSelectorFromFile
   TList         *fInput;           //! input list to the selector
   TList         *fFormulaList;     //! Pointer to a list of coordinated list TTreeFormula (used by Scan and Query)
   TSelector     *fSelectorUpdate;  //! Set to the selector address when it's entry list needs to be updated by the UpdateFormulaLeaves function

protected:
   const   char  *GetNameByIndex(TString &varexp, Int_t *index,Int_t colindex);
   void           DeleteSelectorFromFile();

public:
   TTreePlayer();
   virtual ~TTreePlayer();
   virtual TVirtualIndex *BuildIndex(const TTree *T, const char *majorname, const char *minorname);
   virtual TTree    *CopyTree(const char *selection, Option_t *option
                              ,Long64_t nentries, Long64_t firstentry);
   virtual Long64_t  DrawScript(const char* wrapperPrefix,
                                const char *macrofilename, const char *cutfilename,
                                Option_t *option, Long64_t nentries, Long64_t firstentry);
   virtual Long64_t  DrawSelect(const char *varexp, const char *selection, Option_t *option
                                ,Long64_t nentries, Long64_t firstentry);
   virtual Int_t     Fit(const char *formula ,const char *varexp, const char *selection,Option_t *option ,
                         Option_t *goption ,Long64_t nentries, Long64_t firstentry);
   virtual Int_t     GetDimension() const {return fDimension;}
   TH1              *GetHistogram() const {return fHistogram;}
   virtual Long64_t  GetEntries(const char *selection);
   virtual Long64_t  GetEntriesToProcess(Long64_t firstentry, Long64_t nentries) const;
   virtual Int_t     GetNfill() const {return fSelector->GetNfill();}
   const char       *GetScanFileName() const {return fScanFileName;}
   TTreeFormula     *GetSelect() const    {return fSelector->GetSelect();}
   virtual Long64_t  GetSelectedRows() const {return fSelectedRows;}
   virtual TSelector *GetSelector() const {return fSelector;}
   virtual TSelector *GetSelectorFromFile() const {return fSelectorFromFile;}
   // See TSelectorDraw::GetVar
   TTreeFormula     *GetVar(Int_t i) const {return fSelector->GetVar(i);};
   // See TSelectorDraw::GetVar
   TTreeFormula     *GetVar1() const {return fSelector->GetVar1();}
   // See TSelectorDraw::GetVar
   TTreeFormula     *GetVar2() const {return fSelector->GetVar2();}
   // See TSelectorDraw::GetVar
   TTreeFormula     *GetVar3() const {return fSelector->GetVar3();}
   // See TSelectorDraw::GetVar
   TTreeFormula     *GetVar4() const {return fSelector->GetVar4();}
   // See TSelectorDraw::GetVal
   virtual Double_t *GetVal(Int_t i) const {return fSelector->GetVal(i);};
   // See TSelectorDraw::GetVal
   virtual Double_t *GetV1() const   {return fSelector->GetV1();}
   // See TSelectorDraw::GetVal
   virtual Double_t *GetV2() const   {return fSelector->GetV2();}
   // See TSelectorDraw::GetVal
   virtual Double_t *GetV3() const   {return fSelector->GetV3();}
   // See TSelectorDraw::GetVal
   virtual Double_t *GetV4() const   {return fSelector->GetV4();}
   virtual Double_t *GetW() const    {return fSelector->GetW();}
   virtual Int_t     MakeClass(const char *classname, Option_t *option);
   virtual Int_t     MakeCode(const char *filename);
   virtual Int_t     MakeProxy(const char *classname,
                               const char *macrofilename = 0, const char *cutfilename = 0,
                               const char *option = 0, Int_t maxUnrolling = 3);
   virtual Int_t     MakeReader(const char *classname, Option_t *option);
   TPrincipal       *Principal(const char *varexp, const char *selection, Option_t *option
                               ,Long64_t nentries, Long64_t firstentry);
   virtual Long64_t  Process(const char *filename,Option_t *option, Long64_t nentries, Long64_t firstentry);
   virtual Long64_t  Process(TSelector *selector,Option_t *option,  Long64_t nentries, Long64_t firstentry);
   virtual void      RecursiveRemove(TObject *obj);
   virtual Long64_t  Scan(const char *varexp, const char *selection, Option_t *option
                          ,Long64_t nentries, Long64_t firstentry);
   Bool_t            ScanRedirected() {return fScanRedirect;}
   virtual TSQLResult *Query(const char *varexp, const char *selection, Option_t *option
                             ,Long64_t nentries, Long64_t firstentry);
   virtual void      SetEstimate(Long64_t n);
   void              SetScanRedirect(Bool_t on=kFALSE) {fScanRedirect = on;}
   void              SetScanFileName(const char *name) {fScanFileName=name;}
   virtual void      SetTree(TTree *t) {fTree = t;}
   virtual void      StartViewer(Int_t ww, Int_t wh);
   virtual Int_t     UnbinnedFit(const char *formula ,const char *varexp, const char *selection,Option_t *option
                                 ,Long64_t nentries, Long64_t firstentry);
   virtual void      UpdateFormulaLeaves();

   ClassDef(TTreePlayer,3);  //Manager class to play with TTrees
};

#endif
