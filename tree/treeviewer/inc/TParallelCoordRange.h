// @(#)root/treeviewer:$Id$
// Author: Bastien Dalla Piazza  02/08/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TParallelCoordRange
#define ROOT_TParallelCoordRange

#include "TNamed.h"
#include "TAttLine.h"
#include "TList.h"

class TParallelCoordVar;
class TParallelCoord;
class TParallelCoordSelect;
class TPoint;
class TString;

class TParallelCoordRange : public TNamed, public TAttLine {
public:
   enum EStatusBits {
      kShowOnPad = BIT(15),
      kLiveUpdate = BIT(16)
   };

private:
   Double_t          fMin;    // Min value for the range.
   Double_t          fMax;    // Max value for the range.
   const Double_t    fSize;   // Size of the painted range.
   TParallelCoordVar *fVar;   // Variable owning the range.
   TParallelCoordSelect* fSelect; // Selection owning the range.

   void              PaintSlider(Double_t value,Bool_t fill=kFALSE);
   TPoint*           GetBindingLinePoints(Int_t pos,Int_t mindragged);
   TPoint*           GetSliderPoints(Double_t value);
   TPoint*           GetSliderPoints(Int_t pos);

public:
   TParallelCoordRange();
   TParallelCoordRange(TParallelCoordVar *var, Double_t min=0, Double_t max=0, TParallelCoordSelect* sel=NULL);
   virtual ~TParallelCoordRange();

   virtual void BringOnTop() ;// *MENU*
   virtual void Delete(const Option_t* options=""); // *MENU*
   virtual Int_t DistancetoPrimitive(Int_t px, Int_t py);
   virtual void Draw(Option_t *options="");
   virtual void ExecuteEvent(Int_t entry, Int_t px, Int_t py);
   virtual Double_t GetMin() {return fMin;}
   virtual Double_t GetMax() {return fMax;}
   TParallelCoordVar* GetVar() {return fVar;}
   TParallelCoordSelect* GetSelection() {return fSelect;}
   Bool_t IsIn(Double_t evtval);
   virtual void Paint(Option_t *options);
   virtual void Print(Option_t *options) const; // *MENU*
   virtual void SendToBack(); // *MENU*
   virtual void SetLineColor(Color_t col);
   virtual void SetLineWidth(Width_t wid);

   ClassDef(TParallelCoordRange,1); // A TParallelCoordRange is a range used for parallel coordinates plots.
};


class TParallelCoordSelect : public TList, public TAttLine {
public:
   enum {
      kActivated = BIT(18),
      kShowRanges = BIT(19)
   };

private:
   TString fTitle;            // Title of the selection.

public:
   TParallelCoordSelect();    // Default constructor.
   TParallelCoordSelect(const char* title); // Normal constructor.
   ~TParallelCoordSelect();   // Destructor.

   const char* GetTitle() const {return fTitle.Data();}
   void        SetActivated(Bool_t on);
   void        SetShowRanges(Bool_t s);
   void        SetTitle(const char* title) {fTitle = title;}

   ClassDef(TParallelCoordSelect,1); // A TParallelCoordSelect is a specialised TList to hold TParallelCoordRanges used by TParallelCoord.
};

#endif

