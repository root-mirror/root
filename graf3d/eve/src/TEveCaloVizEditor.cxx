// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveCaloVizEditor.h"
#include "TEveCalo.h"
#include "TEveGValuators.h"
#include "TEveRGBAPaletteEditor.h"
#include "TEveCaloData.h"

#include "TGClient.h"
#include "TGFont.h"
#include "TGedEditor.h"

#include "TGLabel.h"
#include "TGNumberEntry.h"
#include "TGDoubleSlider.h"
#include "TG3DLine.h"
#include "TGButtonGroup.h"
#include "TColor.h"
#include "TGColorSelect.h"

#include "TMathBase.h"
#include "TMath.h"

/** \class TEveCaloVizEditor
\ingroup TEve
GUI editor for TEveCaloEditor.
*/

ClassImp(TEveCaloVizEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveCaloVizEditor::TEveCaloVizEditor(const TGWindow *p, Int_t width, Int_t height,
                                     UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM(0),

   fPlotE(0),
   fPlotEt(0),

   fScaleAbs(0),
   fMaxValAbs(0),
   fMaxTowerH(0),

   fEtaRng(0),
   fPhi(0),
   fPhiOffset(0),
   fDataFrame(0),
   fSliceFrame(0)
{
   MakeTitle("TEveCaloViz");

   TGLabel* label  = 0;
   Int_t    labelW = 45;

   // scaling
   TGHorizontalFrame* scf = new TGHorizontalFrame(this);

   label = new TGLabel(scf, "ScaleAbsolute:");
   scf->AddFrame(label, new TGLayoutHints(kLHintsLeft  | kLHintsBottom));

   fScaleAbs  = new TGCheckButton(scf);
   scf->AddFrame(fScaleAbs, new TGLayoutHints(kLHintsLeft, 3, 5, 3, 0));
   fScaleAbs->Connect("Toggled(Bool_t)", "TEveCaloVizEditor", this, "DoScaleAbs()");


   fMaxValAbs = new TEveGValuator(scf, "MaxEVal:", 70, 0);
   fMaxValAbs->SetLabelWidth(56);
   fMaxValAbs->SetNELength(5);
   fMaxValAbs->SetShowSlider(kFALSE);
   fMaxValAbs->Build();
   fMaxValAbs->SetLimits(0, 1000);
   fMaxValAbs->Connect("ValueSet(Double_t)", "TEveCaloVizEditor", this, "DoMaxValAbs()");
   scf->AddFrame(fMaxValAbs, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   AddFrame(scf, new TGLayoutHints(kLHintsTop, 4, 1, 1, 0));


   // tower height
   fMaxTowerH = new TEveGValuator(this, "MaxTowerH:", 96, 0);
   fMaxTowerH->SetLabelWidth(71);
   fMaxTowerH->SetNELength(5);
   fMaxTowerH->SetShowSlider(kFALSE);
   fMaxTowerH->Build();
   fMaxTowerH->SetLimits(0.1, 500, 501, TGNumberFormat::kNESRealOne);
   fMaxTowerH->Connect("ValueSet(Double_t)", "TEveCaloVizEditor", this, "DoMaxTowerH()");
   AddFrame(fMaxTowerH, new TGLayoutHints(kLHintsTop, 4, 1, 1, 1));


   //______________________________________________________________________________

   fDataFrame = CreateEditorTabSubFrame("Data");

   // E/Et Plot
   {
      TGHorizontalFrame *group = new TGHorizontalFrame(fDataFrame);
      TGCompositeFrame  *labfr = new TGHorizontalFrame(group, 28, 20, kFixedSize);

      TGFont *myfont = gClient->GetFont("-adobe-times-bold-r-*-*-12-*-*-*-*-*-iso8859-1");
      label = new TGLabel(labfr, "Plot:");
      label->SetTextFont(myfont);
      labfr->AddFrame(label, new TGLayoutHints(kLHintsLeft  | kLHintsBottom));
      group->AddFrame(labfr, new TGLayoutHints(kLHintsLeft));

      fPlotE = new TGRadioButton(group, new TGHotString("E"), 11);
      fPlotE->Connect("Clicked()", "TEveCaloVizEditor", this, "DoPlot()");
      group->AddFrame(fPlotE, new TGLayoutHints(kLHintsLeft  | kLHintsBottom, 2, 2, 0, 0));

      fPlotEt = new TGRadioButton(group, new TGHotString("Et"), 22);
      fPlotEt->Connect("Clicked()", "TEveCaloVizEditor", this, "DoPlot()");
      group->AddFrame(fPlotEt, new TGLayoutHints(kLHintsLeft  | kLHintsBottom,  2, 2, 0, 0));

      fDataFrame->AddFrame(group, new TGLayoutHints(kLHintsTop, 4, 1, 1, 0));
   }

   // eta
   fEtaRng = new TEveGDoubleValuator(fDataFrame,"Eta rng:", 40, 0);
   fEtaRng->SetNELength(6);
   fEtaRng->SetLabelWidth(labelW);
   fEtaRng->Build();
   fEtaRng->GetSlider()->SetWidth(195);
   fEtaRng->SetLimits(-5.5, 5.5, TGNumberFormat::kNESRealTwo);
   fEtaRng->Connect("ValueSet()", "TEveCaloVizEditor", this, "DoEtaRange()");
   fDataFrame->AddFrame(fEtaRng, new TGLayoutHints(kLHintsTop, 1, 1, 4, 5));

   // phi
   fPhi = new TEveGValuator(fDataFrame, "Phi:", 90, 0);
   fPhi->SetLabelWidth(labelW);
   fPhi->SetNELength(6);
   fPhi->Build();
   fPhi->SetLimits(-TMath::Pi(), TMath::Pi(), TGNumberFormat::kNESRealTwo);
   fPhi->Connect("ValueSet(Double_t)", "TEveCaloVizEditor", this, "DoPhi()");
   fDataFrame->AddFrame(fPhi, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   fPhiOffset = new TEveGValuator(fDataFrame, "PhiOff:", 90, 0);
   fPhiOffset->SetLabelWidth(labelW);
   fPhiOffset->SetNELength(6);
   fPhiOffset->Build();
   fPhiOffset->SetLimits(0, TMath::Pi(), TGNumberFormat::kNESRealTwo);
   fPhiOffset->Connect("ValueSet(Double_t)", "TEveCaloVizEditor", this, "DoPhi()");
   fDataFrame->AddFrame(fPhiOffset, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   fSliceFrame = new TGVerticalFrame(fDataFrame);
   fDataFrame->AddFrame(fSliceFrame);
}

////////////////////////////////////////////////////////////////////////////////
/// Create slice info gui.

void TEveCaloVizEditor::MakeSliceInfo()
{
   Int_t ns = fM->GetData()->GetNSlices();
   Int_t nf = fSliceFrame->GetList()->GetSize();

   if (ns > nf)
   {
      for (Int_t i=nf; i<ns; ++i)
      {
         TGHorizontalFrame* f = new TGHorizontalFrame(fSliceFrame);

         TEveGValuator* threshold = new TEveGValuator(f,"", 90, 0, i);
         threshold->SetLabelWidth(50);
         threshold->SetNELength(6);
         threshold->SetShowSlider(kFALSE);
         threshold->Build();
         threshold->SetLimits(0, 1000, TGNumberFormat::kNESRealTwo);
         threshold->Connect("ValueSet(Double_t)", "TEveCaloVizEditor", this, "DoSliceThreshold()");
         f->AddFrame(threshold, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

         TGColorSelect* color = new TGColorSelect(f, 0, i);
         f->AddFrame(color, new TGLayoutHints(kLHintsLeft|kLHintsTop, 3, 1, 0, 1));
         color->Connect("ColorSelected(Pixel_t)", "TEveCaloVizEditor", this, "DoSliceColor(Pixel_t)");

         TGNumberEntry* transparency = new TGNumberEntry(f, 0., 2, i,
                                           TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
                                           TGNumberFormat::kNELLimitMinMax, 0, 100);
         transparency->SetHeight(18);
         transparency->GetNumberEntry()->SetToolTipText("Transparency: 0 is opaque, 100 fully transparent.");
         f->AddFrame(transparency, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
         transparency->Connect("ValueSet(Long_t)", "TEveCaloVizEditor", this, "DoSliceTransparency(Long_t)");


         fSliceFrame->AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));
      }
      nf = ns;
   }

   TIter frame_iterator(fSliceFrame->GetList());
   for (Int_t i=0; i<nf; ++i)
   {
      TGFrameElement    *el = (TGFrameElement*)    frame_iterator();
      TGHorizontalFrame *fr = (TGHorizontalFrame*) el->fFrame;
      if (i < ns)
      {
         TEveCaloData::SliceInfo_t &si = fM->GetData()->RefSliceInfo(i);

         TEveGValuator *threshold = (TEveGValuator*) ((TGFrameElement*) fr->GetList()->At(0))->fFrame;
         TGColorSelect *color     = (TGColorSelect*) ((TGFrameElement*) fr->GetList()->At(1) )->fFrame;
         TGNumberEntry *transp    = (TGNumberEntry*) ((TGFrameElement*) fr->GetList()->At(2))->fFrame;

         threshold->GetLabel()->SetText(si.fName);
         threshold->SetValue(si.fThreshold);
         color->SetColor(TColor::Number2Pixel(si.fColor), kFALSE);
         transp->SetNumber(si.fTransparency);

         if (! fr->IsMapped()) {
            fr->MapSubwindows();
            fr->MapWindow();
         }
      }
      else
      {
         if (fr->IsMapped()) {
            fr->UnmapWindow();
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveCaloVizEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TEveCaloViz*>(obj);
   if (fM->GetPlotEt())
   {
      fPlotEt->SetState(kButtonDown, kFALSE);
      fPlotE->SetState(kButtonUp, kFALSE);
   }
   else
   {
      fPlotE->SetState(kButtonDown, kFALSE);
      fPlotEt->SetState(kButtonUp, kFALSE);
   }

   if (fM->fData)
   {
      TGCompositeFrame* p = GetGedEditor()->GetEditorTab("Data");
      if (p->GetList()->IsEmpty())
      {
         p->MapWindow();
         p->MapSubwindows();
      }

      fScaleAbs->SetState(fM->GetScaleAbs() ? kButtonDown : kButtonUp);
      fMaxValAbs->SetValue(fM->GetMaxValAbs());
      fMaxTowerH->SetValue(fM->GetMaxTowerH());

      Double_t min, max;
      fM->GetData()->GetEtaLimits(min, max);
      fEtaRng->SetLimits((Float_t)min, (Float_t)max);
      fEtaRng->SetValues(fM->fEtaMin, fM->fEtaMax);

      fM->GetData()->GetPhiLimits(min, max);
      fPhi->SetLimits(min, max, 101, TGNumberFormat::kNESRealTwo);
      fPhi->SetValue(fM->fPhi);
      if ( fM->GetData()->GetWrapTwoPi())
      {
         fPhi->SetToolTip("Center angle in radians");
         fPhiOffset->SetLimits(1e-3, TMath::Pi(), 101, TGNumberFormat::kNESRealTwo);
      }
      else
      {
         Float_t d = (max-min)*0.5;
         fPhiOffset->SetLimits(1e-3, d, 101, TGNumberFormat::kNESRealTwo);
      }
      fPhiOffset->SetValue(fM->fPhiOffset);
      fPhiOffset->SetToolTip("Phi range in radians");

      MakeSliceInfo();
   }
   else
   {

      fDataFrame->UnmapWindow();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting max tower height.

void TEveCaloVizEditor::DoMaxTowerH()
{
   fM->SetMaxTowerH(fMaxTowerH->GetValue());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for enabling/disabling absolute scale.

void TEveCaloVizEditor::DoScaleAbs()
{
   fM->SetScaleAbs(fScaleAbs->IsOn());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting max E in for absolute scale.

void TEveCaloVizEditor::DoMaxValAbs()
{
   fM->SetMaxValAbs(fMaxValAbs->GetValue());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting E/Et plot.

void TEveCaloVizEditor::DoPlot()
{
   TGButton *btn = (TGButton *) gTQSender;
   Int_t id = btn->WidgetId();

   if (id == fPlotE->WidgetId())
      fPlotEt->SetState(kButtonUp);
   else
      fPlotE->SetState(kButtonUp);

   fM->SetPlotEt(fPlotEt->IsDown());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting eta range.

void TEveCaloVizEditor::DoEtaRange()
{
   fM->SetEta(fEtaRng->GetMin(), fEtaRng->GetMax());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting phi range.

void TEveCaloVizEditor::DoPhi()
{
   fM->SetPhiWithRng(fPhi->GetValue(), fPhiOffset->GetValue());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for SliceThreshold.

void TEveCaloVizEditor::DoSliceThreshold()
{
   TEveGValuator *st = (TEveGValuator *) gTQSender;
   fM->SetDataSliceThreshold(st->WidgetId(), st->GetValue());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for slice info Color.

void TEveCaloVizEditor::DoSliceColor(Pixel_t pixel)
{
   TGColorSelect *cs = (TGColorSelect *) gTQSender;
   fM->SetDataSliceColor(cs->WidgetId(), Color_t(TColor::GetColor(pixel)));
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for slice transparency.

void TEveCaloVizEditor::DoSliceTransparency(Long_t t)
{
   TGNumberEntry *cs = (TGNumberEntry*) gTQSender;
   fM->GetData()->SetSliceTransparency(cs->WidgetId(), t);
   Update();
}

/** \class TEveCalo3DEditor
\ingroup TEve
GUI editor for TEveCalo3D.
*/

ClassImp(TEveCalo3DEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveCalo3DEditor::TEveCalo3DEditor(const TGWindow *p, Int_t width, Int_t height,
                                       UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM(0),
   fFrameTransparency(0)
{
   MakeTitle("TEveCalo3D");

   TGHorizontalFrame* f = new TGHorizontalFrame(this);
   TGLabel* lab = new TGLabel(f, "Frame transparency: ");
   f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 1, 1, 1));

   fFrameTransparency = new TGNumberEntry(f, 0., 2, -1,
                                     TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
                                     TGNumberFormat::kNELLimitMinMax, 0, 100);

   fFrameTransparency->SetHeight(18);
   fFrameTransparency->GetNumberEntry()->SetToolTipText("Transparency: 0 is opaque, 100 fully transparent.");
   f->AddFrame(fFrameTransparency, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
   fFrameTransparency->Connect("ValueSet(Long_t)","TEveCalo3DEditor", this, "DoFrameTransparency()");

   AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveCalo3DEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TEveCalo3D*>(obj);
   fFrameTransparency->SetNumber(fM->GetFrameTransparency());
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for frame transparency.

void TEveCalo3DEditor::DoFrameTransparency()
{
   fM->SetFrameTransparency((Char_t)(fFrameTransparency->GetNumber()));
   Update();
}
