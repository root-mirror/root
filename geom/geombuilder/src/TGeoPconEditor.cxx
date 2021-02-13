// @(#):$Id$
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoPconEditor
\ingroup Geometry_builder

Editor for a TGeoPcon.

\image html geom_pcon_pic.png

\image html geom_pcon_ed.png

*/

#include "TGeoPconEditor.h"
#include "TGeoTabManager.h"
#include "TGeoPcon.h"
#include "TGeoManager.h"
#include "TVirtualGeoPainter.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TGCanvas.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"

ClassImp(TGeoPconEditor);

enum ETGeoPconWid {
   kPCON_NAME, kPCON_NZ, kPCON_PHI1, kPCON_DPHI, kPCON_APPLY, kPCON_UNDO
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for polycone editor

TGeoPconEditor::TGeoPconEditor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fShape   = 0;
   fNsections = 0;
   fSections = 0;
   fNsecti = 0;
   fPhi1i = 0;
   fDPhii = 0;
   fZi = 0;
   fRmini = 0;
   fRmaxi = 0;
   fIsModified = kFALSE;
   fIsShapeEditable = kFALSE;

   fLHsect = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,0,2,2);

   // TextEntry for shape name
   MakeTitle("Name");
   fShapeName = new TGTextEntry(this, new TGTextBuffer(50), kPCON_NAME);
   fShapeName->Resize(135, fShapeName->GetDefaultHeight());
   fShapeName->SetToolTipText("Enter the polycone name");
   fShapeName->Associate(this);
   AddFrame(fShapeName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   MakeTitle("Parameters");
   // Number entry for Nsections
   TGTextEntry *nef;
   TGCompositeFrame *f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Nz"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fENz = new TGNumberEntry(f1, 0., 5, kPCON_NZ);
   fENz->SetNumAttr(TGNumberFormat::kNEAPositive);
   fENz->SetNumStyle(TGNumberFormat::kNESInteger);
   fENz->Resize(100,fENz->GetDefaultHeight());
   nef = (TGTextEntry*)fENz->GetNumberEntry();
   nef->SetToolTipText("Enter the  number of Z sections");
   fENz->Associate(this);
   f1->AddFrame(fENz, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for Phi1
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Phi1"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEPhi1 = new TGNumberEntry(f1, 0., 5, kPCON_PHI1);
   fEPhi1->Resize(100,fEPhi1->GetDefaultHeight());
   nef = (TGTextEntry*)fEPhi1->GetNumberEntry();
   nef->SetToolTipText("Enter the starting phi angle [deg]");
   fEPhi1->Associate(this);
   f1->AddFrame(fEPhi1, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for Dphi
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Dphi"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDPhi = new TGNumberEntry(f1, 0., 5, kPCON_DPHI);
   fEDPhi->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEDPhi->Resize(100,fEDPhi->GetDefaultHeight());
   nef = (TGTextEntry*)fEDPhi->GetNumberEntry();
   nef->SetToolTipText("Enter the phi range [deg]");
   fEDPhi->Associate(this);
   f1->AddFrame(fEDPhi, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // TGCanvas containing sections
   MakeTitle("Pcon sections");
   fCan = new TGCanvas(this, 160, 200,  kSunkenFrame | kDoubleBorder);
   TGCompositeFrame *cont = new TGCompositeFrame(fCan->GetViewPort(), 155, 20, kVerticalFrame | kFixedWidth);
   fCan->SetContainer(cont);
   // labels for #i, Z, Rmin, Rmax
   f1 = new TGCompositeFrame(cont, 160, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "#"), new TGLayoutHints(kLHintsLeft, 2, 20, 6, 0));
   f1->AddFrame(new TGLabel(f1, "Z"), new TGLayoutHints(kLHintsLeft, 2, 20, 6, 0));
   f1->AddFrame(new TGLabel(f1, "Rmin"), new TGLayoutHints(kLHintsLeft, 2, 20, 6, 0));
   f1->AddFrame(new TGLabel(f1, "Rmax"), new TGLayoutHints(kLHintsLeft, 2, 10, 6, 0));
   cont->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 0,0,0,0));
   CreateSections(2);

   AddFrame(fCan, new TGLayoutHints(kLHintsLeft, 0, 0, 4, 4));

   // Delayed draw
   fDFrame = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth | kSunkenFrame);
   fDelayed = new TGCheckButton(fDFrame, "Delayed draw");
   fDFrame->AddFrame(fDelayed, new TGLayoutHints(kLHintsLeft , 2, 2, 4, 4));
   AddFrame(fDFrame,  new TGLayoutHints(kLHintsLeft, 6, 6, 4, 4));

   // Buttons
   fBFrame = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   fApply = new TGTextButton(fBFrame, "Apply");
   fBFrame->AddFrame(fApply, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
   fApply->Associate(this);
   fUndo = new TGTextButton(fBFrame, "Undo");
   fBFrame->AddFrame(fUndo, new TGLayoutHints(kLHintsRight , 2, 2, 4, 4));
   fUndo->Associate(this);
   AddFrame(fBFrame,  new TGLayoutHints(kLHintsLeft, 6, 6, 4, 4));
   fUndo->SetSize(fApply->GetSize());
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPconEditor::~TGeoPconEditor()
{
   if (fSections) delete fSections;
   if (fZi) delete [] fZi;
   if (fRmini) delete [] fRmini;
   if (fRmaxi) delete [] fRmaxi;
   TGFrameElement *el;
   TIter next(GetList());
   while ((el = (TGFrameElement *)next())) {
      if (el->fFrame->IsComposite())
         TGeoTabManager::Cleanup((TGCompositeFrame*)el->fFrame);
   }
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TGeoPconEditor::ConnectSignals2Slots()
{
   fENz->Connect("ValueSet(Long_t)", "TGeoPconEditor", this, "DoNz()");
   fEPhi1->Connect("ValueSet(Long_t)", "TGeoPconEditor", this, "DoPhi()");
   fEDPhi->Connect("ValueSet(Long_t)", "TGeoPconEditor", this, "DoPhi()");
   fApply->Connect("Clicked()", "TGeoPconEditor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoPconEditor", this, "DoUndo()");
   fShapeName->Connect("TextChanged(const char *)", "TGeoPconEditor", this, "DoModified()");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Connect to a given pcon.

void TGeoPconEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA() != TGeoPcon::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoPcon*)obj;
   const char *sname = fShape->GetName();
   if (!strcmp(sname, fShape->ClassName())) fShapeName->SetText("-no_name");
   else fShapeName->SetText(sname);

   Int_t nsections = fShape->GetNz();
   fNsecti = nsections;
   fENz->SetNumber(nsections);
   fEPhi1->SetNumber(fShape->GetPhi1());
   fPhi1i = fShape->GetPhi1();
   fEDPhi->SetNumber(fShape->GetDphi());
   fDPhii = fShape->GetDphi();
   CreateSections(nsections);
   UpdateSections();

   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Change dynamically the number of sections.

void TGeoPconEditor::CreateSections(Int_t inew)
{
   if (inew == fNsections) return;
   if (!fSections) fSections = new TObjArray(8);
   TGCompositeFrame *cont = (TGCompositeFrame*)fCan->GetContainer();
   TGeoPconSection *sect;
   Int_t isect;
   // new sections requested
   if (inew>fNsections) {
      for (isect=fNsections; isect<inew; isect++) {
         sect = new TGeoPconSection(cont, 150, 10, isect);
         fSections->Add(sect);
         cont->AddFrame(sect, fLHsect);
         sect->Connect("Changed(Int_t)", "TGeoPconEditor", this, "DoSectionChange(Int_t)");
      }
   } else {
   // some sections need to be removed
      for (isect=inew; isect<fNsections; isect++) {
         sect = (TGeoPconSection*)fSections->At(isect);
         sect->HideDaughters();
         cont->HideFrame(sect);
         cont->RemoveFrame(sect);
//         sect->Disconnect("Changed(Int_t)", "TGeoPconEditor", this, "DoSectionChange(Int_t)");
         fSections->RemoveAt(isect);
         delete sect;
      }
   }
   fNsections = inew;
   fCan->MapSubwindows();
   cont->Layout();
   cont->MapWindow();
   fCan->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Check validity of sections

Bool_t TGeoPconEditor::CheckSections(Bool_t change)
{
   TGeoPconSection *sect;
   Double_t zmin = 0;
   Double_t rmin = 0, rmax = 1.;
   for (Int_t isect=0; isect<fNsections; isect++) {
      sect = (TGeoPconSection*)fSections->At(isect);
      if (isect && (sect->GetZ()<zmin)) {
         if (!change) return kFALSE;
         sect->SetZ(zmin+1.);
      }
      zmin = sect->GetZ();
      if (sect->GetRmin()<0 ||
          (sect->GetRmax()<0) || ((sect->GetRmin()==0) && (sect->GetRmax()==0))) {
         if (!change) return kFALSE;
         sect->SetRmin(rmin);
         sect->SetRmax(rmax);
      }
      rmin = sect->GetRmin();
      rmax = sect->GetRmax();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Update sections according fShape.

void TGeoPconEditor::UpdateSections()
{
   if (fZi) delete [] fZi;
   if (fRmini) delete [] fRmini;
   if (fRmaxi) delete [] fRmaxi;
   fZi = new Double_t[fNsections];
   fRmini = new Double_t[fNsections];
   fRmaxi = new Double_t[fNsections];
   TGeoPconSection *sect;
   for (Int_t isect=0; isect<fNsections; isect++) {
      sect = (TGeoPconSection*)fSections->At(isect);
      sect->SetZ(fShape->GetZ(isect));
      fZi[isect] = fShape->GetZ(isect);
      sect->SetRmin(fShape->GetRmin(isect));
      fRmini[isect] = fShape->GetRmin(isect);
      sect->SetRmax(fShape->GetRmax(isect));
      fRmaxi[isect] = fShape->GetRmax(isect);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check if shape drawing is delayed.

Bool_t TGeoPconEditor::IsDelayed() const
{
   return (fDelayed->GetState() == kButtonDown);
}

////////////////////////////////////////////////////////////////////////////////
/// Perform name change

void TGeoPconEditor::DoName()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying modifications.

void TGeoPconEditor::DoApply()
{
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled();
   if (!CheckSections()) return;
   // check if number of sections changed
   Bool_t recreate = kFALSE;
   Int_t nz = fENz->GetIntNumber();
   Double_t phi1 = fEPhi1->GetNumber();
   Double_t dphi = fEDPhi->GetNumber();
   if (nz != fShape->GetNz()) recreate = kTRUE;
   TGeoPconSection *sect;
   Int_t isect;
   if (recreate) {
      Double_t *array = new Double_t[3*(nz+1)];
      array[0] = phi1;
      array[1] = dphi;
      array[2] = nz;
      for (isect=0; isect<nz; isect++) {
         sect = (TGeoPconSection*)fSections->At(isect);
         array[3+3*isect] = sect->GetZ();
         array[4+3*isect] = sect->GetRmin();
         array[5+3*isect] = sect->GetRmax();
      }
      fShape->SetDimensions(array);
      delete [] array;
      if (fPad) {
         if (gGeoManager && gGeoManager->GetPainter() && gGeoManager->GetPainter()->IsPaintingShape()) {
            TView *view = fPad->GetView();
            if (!view) {
               fShape->Draw();
               fPad->GetView()->ShowAxis();
            } else {
               const Double_t *orig = fShape->GetOrigin();
               view->SetRange(orig[0]-fShape->GetDX(), orig[1]-fShape->GetDY(), orig[2]-fShape->GetDZ(),
                              orig[0]+fShape->GetDX(), orig[1]+fShape->GetDY(), orig[2]+fShape->GetDZ());
               Update();
            }
         } else Update();
      }
      return;
   }
   // No need to call SetDimensions
   if (TMath::Abs(phi1-fShape->GetPhi1())>1.e-6) fShape->Phi1() = phi1;
   if (TMath::Abs(dphi-fShape->GetDphi())>1.e-6)  fShape->Dphi() = dphi;
   for (isect=0; isect<fNsections; isect++) {
      sect = (TGeoPconSection*)fSections->At(isect);
      fShape->Z(isect) = sect->GetZ();
      fShape->Rmin(isect) = sect->GetRmin();
      fShape->Rmax(isect) = sect->GetRmax();
   }
   fShape->ComputeBBox();
   if (fPad) {
      if (gGeoManager && gGeoManager->GetPainter() && gGeoManager->GetPainter()->IsPaintingShape()) {
         TView *view = fPad->GetView();
         if (!view) {
            fShape->Draw();
            fPad->GetView()->ShowAxis();
         } else {
            const Double_t *orig = fShape->GetOrigin();
            view->SetRange(orig[0]-fShape->GetDX(), orig[1]-fShape->GetDY(), orig[2]-fShape->GetDZ(),
                           orig[0]+fShape->GetDX(), orig[1]+fShape->GetDY(), orig[2]+fShape->GetDZ());
            Update();
         }
      } else Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change parameters of section isect;

void TGeoPconEditor::DoSectionChange(Int_t isect)
{
   TGeoPconSection *sect, *sectlo=0, *secthi=0;
   sect = (TGeoPconSection*)fSections->At(isect);
   if (isect) sectlo = (TGeoPconSection*)fSections->At(isect-1);
   if (isect<fNsections-1) secthi = (TGeoPconSection*)fSections->At(isect+1);
   Double_t z = sect->GetZ();
   if (sectlo && z<sectlo->GetZ()) {
      z = sectlo->GetZ();
      sect->SetZ(z);
   }
   if (secthi && z>secthi->GetZ()) {
      z = secthi->GetZ();
      sect->SetZ(z);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Change number of sections.

void TGeoPconEditor::DoNz()
{
   Int_t nz = fENz->GetIntNumber();
   if (nz < 2) {
      nz = 2;
      fENz->SetNumber(nz);
   }
   CreateSections(nz);
   CheckSections(kTRUE);
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Change phi range.

void TGeoPconEditor::DoPhi()
{
   Double_t phi1 = fEPhi1->GetNumber();
   Double_t dphi = fEDPhi->GetNumber();
   if (TMath::Abs(phi1)>360) fEPhi1->SetNumber(0);
   if (dphi>360) fEDPhi->SetNumber(360);
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for signaling modifications.

void TGeoPconEditor::DoModified()
{
   fApply->SetEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoPconEditor::DoUndo()
{
   fENz->SetNumber(fNsecti);
   CreateSections(fNsecti);
   fEPhi1->SetNumber(fPhi1i);
   fEDPhi->SetNumber(fDPhii);
   TGeoPconSection *sect;
   for (Int_t isect=0; isect<fNsections; isect++) {
      sect = (TGeoPconSection*)fSections->At(isect);
      sect->SetZ(fZi[isect]);
      sect->SetRmin(fRmini[isect]);
      sect->SetRmax(fRmaxi[isect]);
   }
   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

/** \class TGeoPconSection
\ingroup Geometry_builder

Utility frame used by TGeoPcon editor.

*/

ClassImp(TGeoPconSection);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGeoPconSection::TGeoPconSection(const TGWindow *p, UInt_t w, UInt_t h, Int_t id)
                :TGCompositeFrame(p,w,h,kHorizontalFrame | kFixedWidth)
{
   fNumber = id;
   TGTextEntry *nef;
   // Label with number
   AddFrame(new TGLabel(this, TString::Format("#%i",id)), new TGLayoutHints(kLHintsLeft, 2, 4, 6, 0));

   // Z entry
   fEZ = new TGNumberEntry(this, 0., 5);
   fEZ->Resize(40,fEZ->GetDefaultHeight());
   nef = (TGTextEntry*)fEZ->GetNumberEntry();
   nef->SetToolTipText("Enter the Z position");
   fEZ->Associate(this);
   AddFrame(fEZ, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
   // Rmin entry
   fERmin = new TGNumberEntry(this, 0., 5);
   fERmin->SetNumAttr(TGNumberFormat::kNEAPositive);
   fERmin->Resize(40,fERmin->GetDefaultHeight());
   nef = (TGTextEntry*)fERmin->GetNumberEntry();
   nef->SetToolTipText("Enter the minimum radius");
   fERmin->Associate(this);
   AddFrame(fERmin, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
   // Rmax entry
   fERmax = new TGNumberEntry(this, 0., 5);
   fERmax->SetNumAttr(TGNumberFormat::kNEAPositive);
   fERmax->Resize(40,fERmax->GetDefaultHeight());
   nef = (TGTextEntry*)fERmax->GetNumberEntry();
   nef->SetToolTipText("Enter the maximum radius");
   fERmax->Associate(this);
   AddFrame(fERmax, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));

   ConnectSignals2Slots();
   MapSubwindows();
   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoPconSection::~TGeoPconSection()
{
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Hide daughter frames

void TGeoPconSection::HideDaughters()
{
   fEZ->UnmapWindow();
   fERmin->UnmapWindow();
   fERmax->UnmapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Emit Changed(Int_t) signal.

void TGeoPconSection::Changed(Int_t i)
{
   Emit("Changed(Int_t)", i);
}

////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TGeoPconSection::ConnectSignals2Slots()
{
   fEZ->Connect("ValueSet(Long_t)", "TGeoPconSection", this, "DoZ()");
   fERmin->Connect("ValueSet(Long_t)", "TGeoPconSection", this, "DoRmin()");
   fERmax->Connect("ValueSet(Long_t)", "TGeoPconSection", this, "DoRmax()");
}

////////////////////////////////////////////////////////////////////////////////
/// Z value getter

Double_t TGeoPconSection::GetZ() const
{
   return fEZ->GetNumber();
}

////////////////////////////////////////////////////////////////////////////////
/// Rmin value getter

Double_t TGeoPconSection::GetRmin() const
{
   return fERmin->GetNumber();
}

////////////////////////////////////////////////////////////////////////////////
/// Rmax value getter

Double_t TGeoPconSection::GetRmax() const
{
   return fERmax->GetNumber();
}

////////////////////////////////////////////////////////////////////////////////
/// Z value setter

void TGeoPconSection::SetZ(Double_t z)
{
   fEZ->SetNumber(z);
}

////////////////////////////////////////////////////////////////////////////////
/// Rmin value setter

void TGeoPconSection::SetRmin(Double_t rmin)
{
   fERmin->SetNumber(rmin);
}

////////////////////////////////////////////////////////////////////////////////
/// Rmax value setter

void TGeoPconSection::SetRmax(Double_t rmax)
{
   fERmax->SetNumber(rmax);
}

////////////////////////////////////////////////////////////////////////////////
/// Z slot.

void TGeoPconSection::DoZ()
{
   Changed(fNumber);
}

////////////////////////////////////////////////////////////////////////////////
/// Rmin slot.

void TGeoPconSection::DoRmin()
{
   Double_t rmin = fERmin->GetNumber();
   Double_t rmax = fERmax->GetNumber();
   if (rmin>rmax-1.e-8) fERmin->SetNumber(rmax);
   Changed(fNumber);
}

////////////////////////////////////////////////////////////////////////////////
/// Rmax slot.

void TGeoPconSection::DoRmax()
{
   Double_t rmin = fERmin->GetNumber();
   Double_t rmax = fERmax->GetNumber();
   if (rmax<rmin+1.e-8) fERmax->SetNumber(rmin);
   Changed(fNumber);
}
