// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEvePointSetArrayEditor.h"
#include "TEvePointSet.h"
#include "TEveGValuators.h"

#include "TGDoubleSlider.h"

/** \class TEvePointSetArrayEditor
\ingroup TEve
Editor for TEvePointSetArray class.
*/

ClassImp(TEvePointSetArrayEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEvePointSetArrayEditor::TEvePointSetArrayEditor(const TGWindow *p,
                                                 Int_t width, Int_t height,
                                                 UInt_t options, Pixel_t back) :
   TGedFrame(p,width, height, options | kVerticalFrame, back),
   fM(0),
   fRange(0)
{
   fM = 0;
   MakeTitle("TEvePointSetArray");

   fRange = new TEveGDoubleValuator(this,"Range", 200, 0);
   fRange->SetNELength(6);
   fRange->Build();
   fRange->GetSlider()->SetWidth(224);
   fRange->Connect("ValueSet()",
                   "TEvePointSetArrayEditor", this, "DoRange()");
   AddFrame(fRange, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TEvePointSetArrayEditor::~TEvePointSetArrayEditor()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEvePointSetArrayEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TEvePointSetArray*>(obj);

   // printf("FullRange(%f, %f) Selected(%f,%f)\n",
   //        fM->GetMin(), fM->GetMax(), fM->GetCurMin(), fM->GetCurMax());

   fRange->SetLimits(fM->fMin, fM->fMax, TGNumberFormat::kNESRealTwo);
   fRange->SetValues(fM->fCurMin, fM->fCurMax);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting the range of the separating quantity.

void TEvePointSetArrayEditor::DoRange()
{
   fM->SetRange(fRange->GetMin(), fRange->GetMax());
   Update();
}
