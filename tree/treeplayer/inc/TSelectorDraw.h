// @(#)root/treeplayer:$Id$
// Author: Rene Brun   08/01/2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TSelectorDraw
#define ROOT_TSelectorDraw


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TSelectorDraw                                                        //
//                                                                      //
// A specialized TSelector for TTree::Draw.                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TSelector.h"

#include <vector>

class TTreeFormula;
class TTreeFormulaManager;
class TH1;
class TEntryListArray;

class TSelectorDraw : public TSelector {

protected:
   enum EStatusBits { kWarn = BIT(12) };

   TTree         *fTree;           //  Pointer to current Tree
   TTreeFormula **fVar;            //![fDimension] Array of pointers to variables formula
   TTreeFormula  *fSelect;         //  Pointer to selection formula
   TTreeFormulaManager *fManager;  //  Pointer to the formula manager
   TObject       *fTreeElist;      //  pointer to Tree Event list
   TEntryListArray *fTreeElistArray;   //!  pointer to Tree Event list array
   TH1           *fOldHistogram;   //! Pointer to previously used histogram
   Int_t          fAction;         //! Action type
   Long64_t       fDraw;           //! Last entry loop number when object was drawn
   Int_t          fNfill;          //! Total number of histogram fills
   Int_t          fMultiplicity;   //  Indicator of the variability of the size of entries
   Int_t          fDimension;      //  Dimension of the current expression
   Long64_t       fSelectedRows;   //  Number of selected entries
   Long64_t       fOldEstimate;    //  value of Tree fEstimate when selector is called
   Int_t          fForceRead;      //  Force Read flag
   Int_t         *fNbins;          //![fDimension] Number of bins per dimension
   Double_t      *fVmin;           //![fDimension] Minima of varexp columns
   Double_t      *fVmax;           //![fDimension] Maxima of varexp columns
   Double_t       fWeight;         //  Tree weight (see TTree::SetWeight)
   Double_t     **fVal;            //![fSelectedRows][fDimension] Local buffer for the variables
   Int_t          fValSize;
   Double_t      *fW;              //![fSelectedRows]Local buffer for weights
   Bool_t        *fVarMultiple;    //![fDimension] true if fVar[i] has a variable index
   Bool_t         fSelectMultiple; //  true if selection has a variable index
   Bool_t         fCleanElist;     //  true if original Tree elist must be saved
   Bool_t         fObjEval;        //  true if fVar1 returns an object (or pointer to).
   Long64_t       fCurrentSubEntry; // Current subentry when fSelectMultiple is true. Used to fill TEntryListArray

protected:
   virtual void      ClearFormula();
   virtual Bool_t    CompileVariables(const char *varexp="", const char *selection="");
   virtual void      InitArrays(Int_t newsize);

private:
   TSelectorDraw(const TSelectorDraw&);             // not implemented
   TSelectorDraw& operator=(const TSelectorDraw&);  // not implemented

public:
   TSelectorDraw();
   virtual ~TSelectorDraw();

   virtual void      Begin(TTree *tree);
   virtual Int_t     GetAction() const {return fAction;}
   virtual Bool_t    GetCleanElist() const {return fCleanElist;}
   virtual Int_t     GetDimension() const {return fDimension;}
   virtual Long64_t  GetDrawFlag() const {return fDraw;}
   TObject          *GetObject() const {return fObject;}
   Int_t             GetMultiplicity() const   {return fMultiplicity;}
   virtual Int_t     GetNfill() const {return fNfill;}
   TH1              *GetOldHistogram() const {return fOldHistogram;}
   TTreeFormula     *GetSelect() const    {return fSelect;}
   virtual Long64_t  GetSelectedRows() const {return fSelectedRows;}
   TTree            *GetTree() const {return fTree;}
   TTreeFormula     *GetVar(Int_t i) const;
   // See TSelectorDraw::GetVar
   TTreeFormula     *GetVar1() const {return GetVar(0);}
   // See TSelectorDraw::GetVar
   TTreeFormula     *GetVar2() const {return GetVar(1);}
   // See TSelectorDraw::GetVar
   TTreeFormula     *GetVar3() const {return GetVar(2);}
   // See TSelectorDraw::GetVar
   TTreeFormula     *GetVar4() const {return GetVar(3);}
   virtual Double_t *GetVal(Int_t i) const;
   // See TSelectorDraw::GetVal
   virtual Double_t *GetV1() const   {return GetVal(0);}
   // See TSelectorDraw::GetVal
   virtual Double_t *GetV2() const   {return GetVal(1);}
   // See TSelectorDraw::GetVal
   virtual Double_t *GetV3() const   {return GetVal(2);}
   // See TSelectorDraw::GetVal
   virtual Double_t *GetV4() const   {return GetVal(3);}
   virtual Double_t *GetW() const    {return fW;}
   virtual Bool_t    Notify();
   virtual Bool_t    Process(Long64_t /*entry*/) { return kFALSE; }
   virtual void      ProcessFill(Long64_t entry);
   virtual void      ProcessFillMultiple(Long64_t entry);
   virtual void      ProcessFillObject(Long64_t entry);
   virtual void      SetEstimate(Long64_t n);
   virtual UInt_t    SplitNames(const TString &varexp, std::vector<TString> &names);
   virtual void      TakeAction();
   virtual void      TakeEstimate();
   virtual void      Terminate();

   ClassDef(TSelectorDraw,1);  //A specialized TSelector for TTree::Draw
};

#endif
