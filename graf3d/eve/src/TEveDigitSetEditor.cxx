// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveDigitSetEditor.h"
#include "TEveDigitSet.h"

#include "TEveGValuators.h"
#include "TEveRGBAPaletteEditor.h"
#include "TEveGedEditor.h"

#include "TVirtualPad.h"
#include "TH1F.h"
#include "TStyle.h"

#include "TGLabel.h"
#include "TG3DLine.h"

/** \class TEveDigitSetEditor
\ingroup TEve
Editor for TEveDigitSet class.
*/

ClassImp(TEveDigitSetEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveDigitSetEditor::TEveDigitSetEditor(const TGWindow *p, Int_t width, Int_t height,
                                       UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM       (0),
   fPalette (0),

   fHistoButtFrame(0),
   fInfoFrame(0)
{
   MakeTitle("Palette controls");

   fPalette = new TEveRGBAPaletteSubEditor(this);
   AddFrame(fPalette, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 0, 0, 0));
   fPalette->Connect("Changed()", "TEveDigitSetEditor", this, "Update()");

   CreateInfoTab();
}

////////////////////////////////////////////////////////////////////////////////
/// Create information tab.

void TEveDigitSetEditor::CreateInfoTab()
{
   fInfoFrame = CreateEditorTabSubFrame("Info");

   TGCompositeFrame *title1 = new TGCompositeFrame(fInfoFrame, 180, 10,
                                                   kHorizontalFrame |
                                                   kLHintsExpandX   |
                                                   kFixedWidth      |
                                                   kOwnBackground);

   title1->AddFrame(new TGLabel(title1, "TEveDigitSet Info"),
                    new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   title1->AddFrame(new TGHorizontal3DLine(title1),
                    new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   fInfoFrame->AddFrame(title1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));


   fHistoButtFrame = new TGHorizontalFrame(fInfoFrame);
   TGTextButton* b = 0;
   b = new TGTextButton(fHistoButtFrame, "Histo");
   b->SetToolTipText("Show histogram over full range.");
   fHistoButtFrame->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
   b->Connect("Clicked()", "TEveDigitSetEditor", this, "DoHisto()");

   b = new TGTextButton(fHistoButtFrame, "Range Histo");
   b->SetToolTipText("Show histogram over selected range.");
   fHistoButtFrame->AddFrame(b, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
   b->Connect("Clicked()", "TEveDigitSetEditor", this, "DoRangeHisto()");
   fInfoFrame->AddFrame(fHistoButtFrame, new TGLayoutHints(kLHintsExpandX, 2, 0, 0, 0));
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveDigitSetEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TEveDigitSet*>(obj);

   if (fM->fValueIsColor || fM->fPalette == 0) {
      fPalette->UnmapWindow();
   } else {
      fPalette->SetModel(fM->fPalette);
      fPalette->MapWindow();
   }

   if (fM->fHistoButtons)
      fHistoButtFrame->MapWindow();
   else
      fHistoButtFrame->UnmapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Show histogram slot.

void TEveDigitSetEditor::DoHisto()
{
   Int_t min, max;
   if (fM->fPalette) {
      min = fM->fPalette->GetLowLimit();
      max = fM->fPalette->GetHighLimit();
   } else {
      fM->ScanMinMaxValues(min, max);
   }
   PlotHisto(min, max);
}

////////////////////////////////////////////////////////////////////////////////
/// Show ranged histogram slot.

void TEveDigitSetEditor::DoRangeHisto()
{
   Int_t min, max;
   if (fM->fPalette) {
      min = fM->fPalette->GetMinVal();
      max = fM->fPalette->GetMaxVal();
   } else {
      fM->ScanMinMaxValues(min, max);
   }
   PlotHisto(min, max);
}

////////////////////////////////////////////////////////////////////////////////
/// Plots a histogram from digit vales with given range.

void TEveDigitSetEditor::PlotHisto(Int_t min, Int_t max)
{
   Int_t nbins = max-min+1;
   while (nbins > 200)
      nbins /= 2;

   TH1F* h = new TH1F(fM->GetName(), fM->GetTitle(), nbins, min-0.5, max+0.5);
   h->SetDirectory(0);
   h->SetBit(kCanDelete);
   TEveChunkManager::iterator qi(fM->fPlex);
   while (qi.next())
      h->Fill(((TEveDigitSet::DigitBase_t*)qi())->fValue);

   gStyle->SetOptStat(1111111);
   h->Draw();
   gPad->Modified();
   gPad->Update();
}
