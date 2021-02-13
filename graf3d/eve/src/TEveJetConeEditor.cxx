// @(#)root/eve:$Id$
// Author: Matevz Tadel, Jochen Thaeder 2009

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveJetConeEditor.h"
#include "TEveJetCone.h"

/** \class TEveJetConeEditor
\ingroup TEve
GUI editor for TEveJetCone.
*/

ClassImp(TEveJetConeEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveJetConeEditor::TEveJetConeEditor(const TGWindow *p, Int_t width, Int_t height,
                                     UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM(0)
   // Initialize widget pointers to 0
{
   MakeTitle("TEveJetCone");

   // Create widgets
   // fXYZZ = new TGSomeWidget(this, ...);
   // AddFrame(fXYZZ, new TGLayoutHints(...));
   // fXYZZ->Connect("SignalName()", "Reve::TEveJetConeEditor", this, "DoXYZZ()");
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveJetConeEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TEveJetCone*>(obj);

   // Set values of widgets
   // fXYZZ->SetValue(fM->GetXYZZ());
}

