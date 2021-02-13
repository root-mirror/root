// @(#)root/asimage:$Id$
// Author: Reiner Rohlfs   24/03/2002

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun, Fons Rademakers and Reiner Rohlfs *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TASPaletteEditor
\ingroup asimage

This is a GUI window to edit a color palette.
It is called by a pull down menu item of TASImage.
*/

#include "TASImage.h"
#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TASPaletteEditor.h"
#include "TGXYLayout.h"
#include "TGButton.h"
#include "TGComboBox.h"
#include "TGFileDialog.h"
#include "TLine.h"
#include "TROOT.h"
#include "TMath.h"
#include "RConfigure.h"
#include "TVirtualX.h"
#include "strlcpy.h"
#include "snprintf.h"

#ifdef R__HAS_COCOA
#   define X_DISPLAY_MISSING 1
#endif

#ifdef WIN32
#include "Windows4root.h"
#endif

#ifndef WIN32
#   include <afterbase.h>
#else
#   include <win32/config.h>
#   include <win32/afterbase.h>
#endif
#   include <afterimage.h>
extern "C" {
#   include <bmp.h>
}


static const char *gFileTypes[] = {
   "ROOT palette file",  "*.pal.root",
   "ASCII palette file", "*.pal.txt",
   0,                    0
};

static UShort_t gRedRainbow[12] = {
   0x0000, 0x7000, 0x0000, 0x0000, 0x0000,
   0xffff, 0xffff, 0x7000, 0x8000, 0xffff
};
static UShort_t gGreenRainbow[12] = {
   0x0000, 0x0000, 0x0000, 0xffff, 0xffff,
   0xffff, 0x0000, 0x0000, 0x8000, 0xffff
};
static UShort_t gBlueRainbow[12] = {
   0x0000, 0x7000, 0xffff, 0xffff, 0x0000,
   0x0000, 0x0000, 0x0000, 0xa000, 0xffff
};


ClassImp(TASPaletteEditor);

////////////////////////////////////////////////////////////////////////////////
/// Palette editor constructor.
/// The palette editor allows the editing of the color palette of the image.

TASPaletteEditor::TASPaletteEditor(TAttImage *attImage, UInt_t w, UInt_t h)
   : TPaletteEditor(attImage, w, h), TGMainFrame(0, w, h)
{
   SetLayoutManager(new TGXYLayout(this));
   fHisto        = 0;
   fLimitLine[0] = 0;
   fLimitLine[1] = 0;
   fRampFactor   = 0;
   fImagePad     = gPad;

   fPaletteList = new TList;
   fPaletteList->SetOwner();

   fPalette = new TImagePalette(attImage->GetPalette());
   fPaletteList->Add(fPalette);

   // buttons
   TGTextButton *button;

   button = new TGTextButton(this, "&Apply", 1);
   button->SetToolTipText("Apply the palette to the image");
   AddFrame(button, new TGXYLayoutHints(70, 1, 8, 1.8));

   button = new TGTextButton(this, "&Ok", 2);
   button->SetToolTipText("Same as Apply and Cancel button");
   AddFrame(button, new TGXYLayoutHints(70, 3, 8, 1.8));

   button = new TGTextButton(this, "&Cancel", 3);
   button->SetToolTipText("Close this window");
   AddFrame(button, new TGXYLayoutHints(70, 5, 8, 1.8));

   button = new TGTextButton(this, "&Save", 4);
   button->SetToolTipText("Save the palette in a ROOT or an ASCII file");
   AddFrame(button, new TGXYLayoutHints(70, 7.5, 8, 1.8));

   button = new TGTextButton(this, "O&pen", 5);
   button->SetToolTipText("Read a palette from a ROOT or an ASCII file");
   AddFrame(button, new TGXYLayoutHints(70, 9.5, 8, 1.8));

   button = new TGTextButton(this, "&New", 6);
   button->SetToolTipText("Create a new palette (not yet implemented)");
   button->SetState(kButtonDisabled);
   AddFrame(button, new TGXYLayoutHints(70, 12, 8, 1.8));

   button = new TGTextButton(this, "&Edit", 7);
   button->SetToolTipText("Edit a palette (not yet implemented)");
   button->SetState(kButtonDisabled);
   AddFrame(button, new TGXYLayoutHints(70, 14, 8, 1.8));

   fAutoUpdate = new TGCheckButton(this, "Auto Update", 13);
   fAutoUpdate->SetToolTipText("Automatic update of the image (without Apply button)");
   AddFrame(fAutoUpdate, new TGXYLayoutHints(50, 1, 20, 1.8));

   fUnDoButton = new TGTextButton(this, "&Undo", 20);
   fUnDoButton->SetToolTipText("Undo the last modification (repeatable)");
   AddFrame(fUnDoButton, new TGXYLayoutHints(50, 3, 8, 1.8));

   fReDoButton = new TGTextButton(this, "&Redo", 21);
   fReDoButton->SetToolTipText("Undo the last undo operation (repeatable)");
   AddFrame(fReDoButton, new TGXYLayoutHints(60, 3, 8, 1.8));

   button = new TGTextButton(this, "&Log", 8);
   button->SetToolTipText("Apply a log operation to the anchor points of the palette");
   AddFrame(button, new TGXYLayoutHints(50, 15, 8, 1.8));

   button = new TGTextButton(this, "E&xp", 9);
   button->SetToolTipText("Apply a exp operation to the anchor points of the palette");
   AddFrame(button, new TGXYLayoutHints(50, 17, 8, 1.8));

   button = new TGTextButton(this, "L&in", 10);
   button->SetToolTipText("Make the distance of all anchor points constant");
   AddFrame(button, new TGXYLayoutHints(50, 19, 8, 1.8));

   button = new TGTextButton(this, "In&vert", 11);
   button->SetToolTipText("Invert the order of the colors");
   AddFrame(button, new TGXYLayoutHints(60, 17, 8, 1.8));

   fStepButton = new TGCheckButton(this, "Step", 12);
   fStepButton->SetToolTipText("Apply a step function to the palette");
   AddFrame(fStepButton, new TGXYLayoutHints(60, 19, 8, 1.8));

   // ramp: 1, 2 or 4
   TGGroupFrame *rampFrame = new TGGroupFrame(this, "Ramps");
   rampFrame->SetLayoutManager(new TGXYLayout(rampFrame));
   AddFrame(rampFrame, new TGXYLayoutHints(50, 8.5, 14, 6,
            TGXYLayoutHints::kLRubberX | TGXYLayoutHints::kLRubberY |
            TGXYLayoutHints::kLRubberH | TGXYLayoutHints::kLRubberW));

   fRamps[0] = new TGRadioButton(rampFrame, "1", 1);
   fRamps[0]->SetToolTipText("Repeat the palette once");
   rampFrame->AddFrame(fRamps[0], new TGXYLayoutHints(2, 1.4, 5, 1.8));

   fRamps[1] = new TGRadioButton(rampFrame, "2", 2);
   fRamps[1]->SetToolTipText("Repeat the palette twice");
   rampFrame->AddFrame(fRamps[1], new TGXYLayoutHints(2, 3.3, 5, 1.8));

   fRamps[2] = new TGRadioButton(rampFrame, "4", 4);
   fRamps[2]->SetToolTipText("Repeat the palette four times");
   rampFrame->AddFrame(fRamps[2], new TGXYLayoutHints(8, 3.3, 5, 1.8));

   fRamps[0]->Associate(this);
   fRamps[1]->Associate(this);
   fRamps[2]->Associate(this);

   // the histogram of the data
   fHistCanvas = new TRootEmbeddedCanvas("data hist", this, 300, 50);
   AddFrame(fHistCanvas, new TGXYLayoutHints(1, 1, 48, 20,
            TGXYLayoutHints::kLRubberW | TGXYLayoutHints::kLRubberH));

   const ASImage *image = ((TASImage*)attImage)->GetImage();
   if (image && image->alt.vector) {
      Int_t pixel;
      Double_t *data = image->alt.vector;
      Int_t numPixel = image->width * image->height;
      Int_t numBins = numPixel / 20;
      numBins = (numBins < 10) ? 10 : (numBins > 200) ? 200 : numBins;

      // get min and max value of image
      fMinValue = fMaxValue = *image->alt.vector;
      for (pixel = 1; pixel < numPixel; pixel++) {
         if (fMinValue > *(data + pixel)) fMinValue = *(data + pixel);
         if (fMaxValue < *(data + pixel)) fMaxValue = *(data + pixel);
      }

      fHisto = new TH1D("Statistics", "Pixel histogram of unzoomed image    ",
                        numBins, fMinValue, fMaxValue);
      for (pixel = 0; pixel < numPixel; pixel++)
         fHisto->Fill(*(data + pixel));

      fHisto->Draw("HIST");
      fHisto->GetXaxis()->SetLabelFont(63);
      fHisto->GetXaxis()->SetLabelSize(10);
      fHisto->GetYaxis()->SetLabelFont(63);
      fHisto->GetYaxis()->SetLabelSize(10);

      fLimitLine[0] = new LimitLine(fMinValue + fPalette->fPoints[1] * (fMaxValue - fMinValue),
                                    0, fHisto->GetMaximum(), this);
      fLimitLine[0]->Draw();
      fLimitLine[1] = new LimitLine(fMinValue + fPalette->fPoints[fPalette->fNumPoints - 2] *
                                    (fMaxValue - fMinValue), 0, fHisto->GetMaximum(), this);
      fLimitLine[1]->Draw();
   }

   // the combobox of different palettes
   fComboBox = new TGComboBox(this, 100);
   AddFrame(fComboBox, new TGXYLayoutHints(50, 6, 14, 2));

   fComboBox->AddEntry("Rainbow", 0);
   fComboBox->AddEntry("Grey", 1);
   fComboBox->AddEntry("Hot",  2);
   fComboBox->AddEntry("Cold", 3);
   fComboBox->AddEntry("Bowlerhat", 4);
   fComboBox->AddEntry("", 5);


   // the palette
   fPaletteCanvas = new TRootEmbeddedCanvas("palette", this, 300, 50);
   AddFrame(fPaletteCanvas, new TGXYLayoutHints(1, 22, 78, 2.5,
            TGXYLayoutHints::kLRubberW | TGXYLayoutHints::kLRubberY));

   fPaintPalette = new PaintPalette(&fPalette, attImage);
   fPaintPalette->Draw();

   MapSubwindows();
   Layout();

   SetWindowName("Palette Editor");
   SetIconName("Palette Editor");

   MapWindow();

   UpdateScreen(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Palette editor destructor. Deletes all frames and their layout hints.

TASPaletteEditor::~TASPaletteEditor()
{
   TGFrameElement *ptr;

   // delete all frames and layout hints
   if (fList) {
      TIter next(fList);
      while ((ptr = (TGFrameElement *) next())) {
         if (ptr->fLayout)
            delete ptr->fLayout;
         if (ptr->fFrame)
            delete ptr->fFrame;
      }
   }

   delete fHisto;
   delete fPaintPalette;
   delete fLimitLine[0];
   delete fLimitLine[1];
   delete fPaletteList;
}

////////////////////////////////////////////////////////////////////////////////
/// Close editor.

void TASPaletteEditor::CloseWindow()
{
   TPaletteEditor::CloseWindow();
   delete this;
}

////////////////////////////////////////////////////////////////////////////////
/// Process all editor mouse events

Bool_t TASPaletteEditor::ProcessMessage(Long_t msg, Long_t param1, Long_t param2)
{
   switch (GET_MSG(msg)) {

      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {

            case kCM_COMBOBOX:
               NewPalette(param2);
               break;

            case kCM_RADIOBUTTON:
               SetRamp(param1);
               break;

            case kCM_CHECKBUTTON:
               if (param1 == 12)
                  SetStep();
               break;

            case kCM_BUTTON:
               switch (param1) {

                  case 1 :  // Apply
                     fAttImage->SetPalette(fPalette);
                     fImagePad->Modified();
                     fImagePad->Update();
                     break;

                  case 2 :  // OK
                     fAttImage->SetPalette(fPalette);
                     fImagePad->Modified();
                     fImagePad->Update();
                     CloseWindow();
                     break;

                  case 3 : // Cancel
                     CloseWindow();
                     break;

                  case 4 : // Save
                     Save();
                     break;

                  case 5 : // Open
                     Open();
                     break;

                  case 8: // log
                     LogPalette();
                     break;

                  case 9: // exp
                     ExpPalette();
                     break;

                  case 10: // lin
                     LinPalette();
                     break;

                  case 11: // invert
                     InvertPalette();
                     break;


                  case 20: // undo
                     fPalette = (TImagePalette*)(fPaletteList->Before(fPalette));
                     if (fAutoUpdate->GetState() == kButtonDown) {
                        fAttImage->SetPalette(fPalette);
                        fImagePad->Modified();
                        fImagePad->Update();
                     }
                     if (fPalette) UpdateScreen(kTRUE);
                     break;

                  case 21: // redo
                     fPalette = (TImagePalette*)(fPaletteList->After(fPalette));
                     if (fAutoUpdate->GetState() == kButtonDown) {
                        fAttImage->SetPalette(fPalette);
                        fImagePad->Modified();
                        fImagePad->Update();
                     }
                     if (fPalette) UpdateScreen(kTRUE);
                     break;

                  default: ;
               }
               break;

            default: ;
         }
         break;

      default: ;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// The newPalette is inserted in the list of palettes (fPaletteList) and
/// fPalette is set to the newPalette. Protected method,

void TASPaletteEditor::InsertNewPalette(TImagePalette *newPalette)
{
   // first remove all palettes in the list which are behind the
   // current palette
   TObject *obj;
   while ((obj = fPaletteList->After(fPalette)) != 0)
      delete fPaletteList->Remove(obj);

   // add new palette and make it to the current palette
   fPaletteList->Add(newPalette);
   fPalette = newPalette;

   // update the image
   if (fAutoUpdate->GetState() == kButtonDown) {
      fAttImage->SetPalette(fPalette);
      fImagePad->Modified();
      fImagePad->Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Saves the current palette either into a ROOT file or in an ASCII file.
/// It is called by the Save - button. Protected method.

void TASPaletteEditor::Save()
{
   TGFileInfo fi;
   fi.fFileTypes = gFileTypes;
   static Bool_t overwr = kFALSE;
   fi.fOverwrite = overwr;

   new TGFileDialog(gClient->GetRoot(), this, kFDSave, &fi);
   overwr = fi.fOverwrite;
   if (fi.fFilename == 0)
      return;

   if (strcmp(".pal.txt", fi.fFilename + strlen(fi.fFilename) - 8) == 0) {
      // write into an ASCII file
      FILE *fl = fopen(fi.fFilename, "w");
      if (!fl) return;
      fprintf(fl, "%u\n", fPalette->fNumPoints);
      for (Int_t pt = 0; pt < Int_t(fPalette->fNumPoints); pt++)
         fprintf(fl, "%10.9f %04hx %04hx %04hx %04hx\n",
                 fPalette->fPoints[pt],
                 fPalette->fColorRed[pt],
                 fPalette->fColorGreen[pt],
                 fPalette->fColorBlue[pt],
                 fPalette->fColorAlpha[pt] );
      fclose(fl);
   } else {
      // write into a ROOT file
      char fn[512];
      if (strcmp(".pal.root", fi.fFilename + strlen(fi.fFilename) - 9) != 0)
         snprintf(fn,512, "%s%s", fi.fFilename, ".pal.root");
      else
         strlcpy(fn, fi.fFilename,512);

      gROOT->ProcessLine(Form("gROOT->SaveObjectAs((TASPaletteEditor*)0x%lx,\"%s\",\"%s\");",(ULong_t)this,fn,"q"));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Opens either a ROOT file or an ASCII file and reads a palette.
/// It is called by the Open - button. Protected method.

void TASPaletteEditor::Open()
{
   TGFileInfo fi;
   fi.fFileTypes = gFileTypes;

   new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
   if (fi.fFilename == 0)
      return;

   TImagePalette *newPalette;

   if (strcmp(".pal.txt", fi.fFilename + strlen(fi.fFilename) - 8) == 0) {
      FILE *fl = fopen(fi.fFilename, "r");
      if (!fl) return;
      UInt_t numPoints;
      // coverity [Calling risky function : FALSE]
      if (fscanf(fl, "%u\n", &numPoints)) {;}
      newPalette = new TImagePalette(numPoints);
      for (Int_t pt = 0; pt < Int_t(numPoints); pt++)
         // coverity [Calling risky function : FALSE]
         if (fscanf(fl, "%lf %hx %hx %hx %hx\n",
                newPalette->fPoints + pt,
                newPalette->fColorRed + pt,
                newPalette->fColorGreen + pt,
                newPalette->fColorBlue + pt,
                    newPalette->fColorAlpha + pt )) {;}
             fclose(fl);
   } else {
      // read from a ROOT file
      char fn[512];
      if (strcmp(".pal.root", fi.fFilename + strlen(fi.fFilename) - 9) != 0)
         snprintf(fn,512, "%s%s", fi.fFilename, ".pal.root");
      else
         strlcpy(fn, fi.fFilename,512);
      TDirectory *dirsav = gDirectory;

      TFile *fsave = new TFile(fn, "READ");
      if (!fsave->IsOpen()) {
         delete fsave;
         return;
      }

      newPalette = (TImagePalette*)fsave->Get("TImagePalette");
      delete fsave;
      if (dirsav) dirsav->cd();
      if (!newPalette)
         return;
   }

   InsertNewPalette(newPalette);
   UpdateScreen(kTRUE);

   fComboBox->Select(5);  // empty entry
}

////////////////////////////////////////////////////////////////////////////////
/// All widgets of the screen are updated with the current palette.
/// Protected method.

void TASPaletteEditor::UpdateScreen(Bool_t histoUpdate)
{
   // update the color palette
   fPaletteCanvas->GetCanvas()->Modified();
   fPaletteCanvas->GetCanvas()->Update();

   if (histoUpdate) {
      // update the limit lines
      Double_t xPos = fMinValue + fPalette->fPoints[1] * (fMaxValue - fMinValue);
      fLimitLine[0]->SetX1(xPos);
      fLimitLine[0]->SetX2(xPos);

      xPos = fMinValue + fPalette->fPoints[fPalette->fNumPoints - 2] * (fMaxValue - fMinValue);
      fLimitLine[1]->SetX1(xPos);
      fLimitLine[1]->SetX2(xPos);

      fHistCanvas->GetCanvas()->Modified();
      fHistCanvas->GetCanvas()->Update();
   }

   // update undo / redo button
   fUnDoButton->SetState(fPalette == fPaletteList->First() ? kButtonDisabled : kButtonUp);
   fReDoButton->SetState(fPalette == fPaletteList->Last()  ? kButtonDisabled : kButtonUp);

   // test if it is a step palette
   EButtonState step = kButtonDown;

   Int_t pt;
   for (pt = 2; pt < Int_t(fPalette->fNumPoints - 2); pt += 2)
      if (TMath::Abs(fPalette->fPoints[pt] - fPalette->fPoints[pt + 1])  > 0.0001 ||
          fPalette->fColorRed[pt]   != fPalette->fColorRed[pt-1]   ||
          fPalette->fColorGreen[pt] != fPalette->fColorGreen[pt-1] ||
          fPalette->fColorBlue[pt]  != fPalette->fColorBlue[pt-1])
         step = kButtonUp;
   fStepButton->SetState(step);

   // find the ramp factor
   fRampFactor = 4;
   Int_t off = (fPalette->fNumPoints - 2) / 4;
   for (pt = 0; pt < Int_t(fPalette->fNumPoints - 2) / 4 * 3; pt++)
      if (fPalette->fColorRed[pt + 1]   != fPalette->fColorRed[pt + 1 + off]   ||
          fPalette->fColorGreen[pt + 1] != fPalette->fColorGreen[pt + 1 + off] ||
          fPalette->fColorBlue[pt + 1]  != fPalette->fColorBlue[pt + 1 + off]  ||
          fPalette->fColorAlpha[pt + 1] != fPalette->fColorAlpha[pt + 1 + off]) {
         fRampFactor = 2;
         break;
      }

   off = (fPalette->fNumPoints - 2) / 2;
   for (pt = 0; pt < Int_t(fPalette->fNumPoints - 2) / 2; pt++)
      if (fPalette->fColorRed[pt + 1]   != fPalette->fColorRed[pt + 1 + off]   ||
          fPalette->fColorGreen[pt + 1] != fPalette->fColorGreen[pt + 1 + off] ||
          fPalette->fColorBlue[pt + 1]  != fPalette->fColorBlue[pt + 1 + off]  ||
          fPalette->fColorAlpha[pt + 1] != fPalette->fColorAlpha[pt + 1 + off]) {
         fRampFactor = 1;
         break;
      }

   fRamps[0]->SetState(fRampFactor == 1 ? kButtonDown : kButtonUp);
   fRamps[1]->SetState(fRampFactor == 2 ? kButtonDown : kButtonUp);
   fRamps[2]->SetState(fRampFactor == 4 ? kButtonDown : kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// The anchor points are rescaled by a log operation.
/// It is called by the log - button. Protected method.

void TASPaletteEditor::LogPalette()
{
   TImagePalette *newPalette = new TImagePalette(*fPalette);

   Double_t delta = fPalette->fPoints[fPalette->fNumPoints-2] - fPalette->fPoints[1];

   for (Int_t pt = 2; pt < Int_t(fPalette->fNumPoints - 2); pt++)
      newPalette->fPoints[pt] = fPalette->fPoints[1] +
         TMath::Log(fPalette->fPoints[pt] - fPalette->fPoints[1] + 1) /
         TMath::Log(delta + 1) * delta;

   InsertNewPalette(newPalette);
   UpdateScreen(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// The anchor points are rescaled by a exp operation.
/// It is called by the exp - button. Protected method.

void TASPaletteEditor::ExpPalette()
{
   TImagePalette *newPalette = new TImagePalette(*fPalette);

   Double_t delta = fPalette->fPoints[fPalette->fNumPoints-2] - fPalette->fPoints[1];

   for (Int_t pt = 2; pt < Int_t(fPalette->fNumPoints - 2); pt++)
      newPalette->fPoints[pt] = fPalette->fPoints[1] +
         TMath::Exp((fPalette->fPoints[pt] - fPalette->fPoints[1]) *
         TMath::Log(delta + 1) / delta) - 1;

   InsertNewPalette(newPalette);
   UpdateScreen(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// The anchor points are rescaled to be linar.
/// It is called by the lin - button. Protected method.

void TASPaletteEditor::LinPalette()
{
   TImagePalette *newPalette = new TImagePalette(*fPalette);

   Double_t delta = fPalette->fPoints[fPalette->fNumPoints-2] - fPalette->fPoints[1];
   if (fStepButton->GetState() == kButtonUp) {
      for (Int_t pt = 2; pt < Int_t(fPalette->fNumPoints - 2); pt++)
         newPalette->fPoints[pt] = fPalette->fPoints[1] +
            delta * (pt - 1) / (fPalette->fNumPoints - 3);
   } else {
      for (Int_t pt = 0; pt < Int_t(fPalette->fNumPoints - 4); pt += 2) {
         newPalette->fPoints[pt + 3] = fPalette->fPoints[1] + delta * (pt + 2) /
                                       (fPalette->fNumPoints - 2) ;
         newPalette->fPoints[pt + 2] = newPalette->fPoints[pt + 3];
      }
   }

   InsertNewPalette(newPalette);
   UpdateScreen(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// The palette is inverted.
/// It is called by the invert - button. Protected method.

void TASPaletteEditor::InvertPalette()
{
   TImagePalette *newPalette = new TImagePalette(*fPalette);

   Int_t pt;
   for (pt = 0; pt < Int_t(fPalette->fNumPoints); pt++)  {
      newPalette->fColorRed[pt]   = fPalette->fColorRed[fPalette->fNumPoints - 1 - pt];
      newPalette->fColorGreen[pt] = fPalette->fColorGreen[fPalette->fNumPoints - 1 - pt];
      newPalette->fColorBlue[pt]  = fPalette->fColorBlue[fPalette->fNumPoints - 1 - pt];
      newPalette->fColorAlpha[pt] = fPalette->fColorAlpha[fPalette->fNumPoints - 1 - pt];
   }

   for (pt = 2; pt < Int_t(fPalette->fNumPoints - 2); pt++)
      newPalette->fPoints[pt] = fPalette->fPoints[1] +
         fPalette->fPoints[fPalette->fNumPoints - 2] -
         fPalette->fPoints[fPalette->fNumPoints - 1 - pt];

   InsertNewPalette(newPalette);
   UpdateScreen(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// A new palette is created, depending on the id.
/// It is called by the combo box. Protected method.

void TASPaletteEditor::NewPalette(Long_t id)
{
   if (id == 5) // empty entry
      return;

   TImagePalette *newPalette;

   Double_t delta = fPalette->fPoints[fPalette->fNumPoints-2] - fPalette->fPoints[1];
   UInt_t   numPt;

   numPt = id == 0 ? 12 : 13;
   newPalette = new TImagePalette(numPt);
   Int_t pt;
   for (pt = 1; pt < Int_t(numPt - 1); pt++) {
      newPalette->fPoints[pt] = fPalette->fPoints[1] + (pt - 1) * delta / (numPt - 3);
      newPalette->fColorAlpha[pt] = 0xffff;
   }

   switch (id) {
      case 0:  // rainbow
         memcpy(newPalette->fColorRed + 1,   gRedRainbow,   12 * sizeof(UShort_t));
         memcpy(newPalette->fColorGreen + 1, gGreenRainbow, 12 * sizeof(UShort_t));
         memcpy(newPalette->fColorBlue + 1,  gBlueRainbow,  12 * sizeof(UShort_t));
         break;

      case 1:  // gray
         for (pt = 1; pt < Int_t(numPt - 1); pt++) {
            newPalette->fColorRed[pt]   = 0xffff * (pt - 1) / (numPt - 3);
            newPalette->fColorGreen[pt] = 0xffff * (pt - 1) / (numPt - 3);
            newPalette->fColorBlue[pt]  = 0xffff * (pt - 1) / (numPt - 3);
         }
         break;

      case 2:  // hot (red)
         for (pt = 1; pt < Int_t(numPt - 1) / 2; pt++) {
            newPalette->fColorRed[pt]   = 0xffff * (pt - 1) / ((numPt - 3) / 2);
            newPalette->fColorGreen[pt] = 0;
            newPalette->fColorBlue[pt]  = 0;
         }

         for (; pt < Int_t(numPt - 1); pt++) {
            newPalette->fColorRed[pt]   = 0xffff;
            newPalette->fColorGreen[pt] = 0xffff * (pt - (numPt - 1) / 2) / ((numPt - 3) / 2);
            newPalette->fColorBlue[pt]  = 0xffff * (pt - (numPt - 1) / 2) / ((numPt - 3) / 2);
         }
         break;

      case 3:  // cold (blue)
         for (pt = 1; pt < Int_t(numPt - 1) / 2; pt++) {
            newPalette->fColorRed[pt]   = 0;
            newPalette->fColorGreen[pt] = 0;
            newPalette->fColorBlue[pt]  = 0xffff * (pt - 1) / ((numPt - 3) / 2);
         }

         for (; pt < Int_t(numPt - 1); pt++) {
            newPalette->fColorRed[pt]   = 0xffff * (pt - (numPt - 1) / 2) / ((numPt - 3) / 2);
            newPalette->fColorGreen[pt] = 0xffff * (pt - (numPt - 1) / 2) / ((numPt - 3) / 2);
            newPalette->fColorBlue[pt]  = 0xffff;
         }
         break;

      case 4:  // bolwerhat
         for (pt = 1; pt < Int_t(numPt + 1) / 2; pt++) {
            newPalette->fColorRed[pt]   = newPalette->fColorRed[numPt - pt - 1]
                                        = 0xffff * (pt - 1) / ((numPt - 3) / 2);
            newPalette->fColorGreen[pt] = newPalette->fColorGreen[numPt - pt - 1]
                                        = 0xffff * (pt - 1) / ((numPt - 3) / 2);
            newPalette->fColorBlue[pt]  = newPalette->fColorBlue[numPt - pt - 1]
                                        = 0xffff * (pt - 1) / ((numPt - 3) / 2);
         }
         break;
   }

   newPalette->fPoints[0]     = 0;
   newPalette->fColorRed[0]   = newPalette->fColorRed[1];
   newPalette->fColorGreen[0] = newPalette->fColorGreen[1];
   newPalette->fColorBlue[0]  = newPalette->fColorBlue[1];
   newPalette->fColorAlpha[0] = newPalette->fColorAlpha[1];

   newPalette->fPoints[newPalette->fNumPoints-1]     = 1.0;
   newPalette->fColorRed[newPalette->fNumPoints-1]   = newPalette->fColorRed[newPalette->fNumPoints-2];
   newPalette->fColorGreen[newPalette->fNumPoints-1] = newPalette->fColorGreen[newPalette->fNumPoints-2];
   newPalette->fColorBlue[newPalette->fNumPoints-1]  = newPalette->fColorBlue[newPalette->fNumPoints-2];
   newPalette->fColorAlpha[newPalette->fNumPoints-1] = newPalette->fColorAlpha[newPalette->fNumPoints-2];

   InsertNewPalette(newPalette);
   UpdateScreen(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a step palette. This is called by the step - check button.
/// Protected method.

void TASPaletteEditor::SetStep()
{
   TImagePalette *newPalette;

   if (fStepButton->GetState() == kButtonDown) {
      // change colors in steps
      newPalette = new TImagePalette(fPalette->fNumPoints * 2 - 2);
      Double_t fkt = (Double_t)(fPalette->fNumPoints - 3) / (fPalette->fNumPoints - 2);
      for (Int_t pt = 1; pt < Int_t(fPalette->fNumPoints - 1); pt++) {
         newPalette->fPoints[pt * 2 - 1] = fPalette->fPoints[1] + (fPalette->fPoints[pt] - fPalette->fPoints[1]) * fkt;
         newPalette->fPoints[pt * 2] = fPalette->fPoints[1] + (fPalette->fPoints[pt + 1] - fPalette->fPoints[1]) * fkt;
         newPalette->fColorRed[pt * 2 - 1]   = newPalette->fColorRed[pt * 2]   = fPalette->fColorRed[pt];
         newPalette->fColorGreen[pt * 2 - 1] = newPalette->fColorGreen[pt * 2] = fPalette->fColorGreen[pt];
         newPalette->fColorBlue[pt * 2 - 1]  = newPalette->fColorBlue[pt * 2]  = fPalette->fColorBlue[pt];
         newPalette->fColorAlpha[pt * 2 - 1] = newPalette->fColorAlpha[pt * 2] = fPalette->fColorAlpha[pt];
      }
   } else {
      // continuous change of colors
      newPalette = new TImagePalette(fPalette->fNumPoints / 2 + 1);
      Double_t fkt = (Double_t) (fPalette->fPoints[fPalette->fNumPoints - 2] - fPalette->fPoints[1]) /
                                (fPalette->fPoints[fPalette->fNumPoints - 3] - fPalette->fPoints[1]);
      for (Int_t pt = 1; pt < Int_t(newPalette->fNumPoints - 1); pt++) {
         newPalette->fPoints[pt] = fPalette->fPoints[pt * 2 -1] * fkt;
         newPalette->fColorRed[pt]   = fPalette->fColorRed[pt * 2 - 1];
         newPalette->fColorGreen[pt] = fPalette->fColorGreen[pt * 2 - 1];
         newPalette->fColorBlue[pt]  = fPalette->fColorBlue[pt * 2 - 1];
         newPalette->fColorAlpha[pt] = fPalette->fColorAlpha[pt * 2 - 1];
      }
   }

   newPalette->fPoints[0]     = fPalette->fPoints[0];
   newPalette->fColorRed[0]   = fPalette->fColorRed[0];
   newPalette->fColorGreen[0] = fPalette->fColorGreen[0];
   newPalette->fColorBlue[0]  = fPalette->fColorBlue[0];
   newPalette->fColorAlpha[0] = fPalette->fColorAlpha[0];

   newPalette->fPoints[newPalette->fNumPoints-2]     = fPalette->fPoints[fPalette->fNumPoints-2];
   newPalette->fPoints[newPalette->fNumPoints-1]     = fPalette->fPoints[fPalette->fNumPoints-1];
   newPalette->fColorRed[newPalette->fNumPoints-1]   = fPalette->fColorRed[fPalette->fNumPoints-1];
   newPalette->fColorGreen[newPalette->fNumPoints-1] = fPalette->fColorGreen[fPalette->fNumPoints-1];
   newPalette->fColorBlue[newPalette->fNumPoints-1]  = fPalette->fColorBlue[fPalette->fNumPoints-1];
   newPalette->fColorAlpha[newPalette->fNumPoints-1] = fPalette->fColorAlpha[fPalette->fNumPoints-1];

   InsertNewPalette(newPalette);
   UpdateScreen(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// The palette is repeated up to 4 times.
/// This is called by one of the ramp radio buttons. Protected method.

void TASPaletteEditor::SetRamp(Long_t ramp)
{
   if (ramp == fRampFactor)
      return;

   Int_t ptPerRamp = (fPalette->fNumPoints - 2) / fRampFactor;
   TImagePalette *newPalette = new TImagePalette(ptPerRamp * ramp + 2);

   Double_t delta = fPalette->fPoints[fPalette->fNumPoints-2] - fPalette->fPoints[1];
   for (Int_t rp = 0; rp < ramp; rp++) {
      for (Int_t pt = 0; pt < ptPerRamp; pt++) {
         newPalette->fPoints[1 + pt + rp * ptPerRamp] = fPalette->fPoints[1] +
              delta / ramp * rp +
              (fPalette->fPoints[1+pt] - fPalette->fPoints[1]) * fRampFactor / ramp;
         newPalette->fColorRed  [1 + pt + rp * ptPerRamp] = fPalette->fColorRed  [1 + pt];
         newPalette->fColorGreen[1 + pt + rp * ptPerRamp] = fPalette->fColorGreen[1 + pt];
         newPalette->fColorBlue [1 + pt + rp * ptPerRamp] = fPalette->fColorBlue [1 + pt];
         newPalette->fColorAlpha[1 + pt + rp * ptPerRamp] = fPalette->fColorAlpha[1 + pt];
      }
   }

   newPalette->fPoints[0]     = fPalette->fPoints[0];
   newPalette->fColorRed[0]   = fPalette->fColorRed[0];
   newPalette->fColorGreen[0] = fPalette->fColorGreen[0];
   newPalette->fColorBlue[0]  = fPalette->fColorBlue[0];
   newPalette->fColorAlpha[0] = fPalette->fColorAlpha[0];

   newPalette->fPoints[newPalette->fNumPoints-2]     = fPalette->fPoints[fPalette->fNumPoints-2];
   newPalette->fPoints[newPalette->fNumPoints-1]     = fPalette->fPoints[fPalette->fNumPoints-1];
   newPalette->fColorRed[newPalette->fNumPoints-1]   = fPalette->fColorRed[fPalette->fNumPoints-1];
   newPalette->fColorGreen[newPalette->fNumPoints-1] = fPalette->fColorGreen[fPalette->fNumPoints-1];
   newPalette->fColorBlue[newPalette->fNumPoints-1]  = fPalette->fColorBlue[fPalette->fNumPoints-1];
   newPalette->fColorAlpha[newPalette->fNumPoints-1] = fPalette->fColorAlpha[fPalette->fNumPoints-1];

   InsertNewPalette(newPalette);
   UpdateScreen(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Updates the range of the palette.
/// This is called after the blue limit lines were moved to define
/// a new range.

void TASPaletteEditor::UpdateRange()
{
   if (fMaxValue == fMinValue)
      return;

   TImagePalette *newPalette = new TImagePalette(*fPalette);

   Double_t l0 = fLimitLine[0]->GetX1();
   Double_t l1 = fLimitLine[1]->GetX1();
   l0 = (l0 < fMinValue) ? fMinValue : ((l0 > fMaxValue) ?  fMaxValue : l0);
   l1 = (l1 < fMinValue) ? fMinValue : ((l1 > fMaxValue) ?  fMaxValue : l1);
   if (l0 > l1) {
      Double_t tmp = l0;
      l0 = l1;
      l1 = tmp;
   }

   Double_t oldDelta = fPalette->fPoints[fPalette->fNumPoints - 2] - fPalette->fPoints[1];
   Double_t newDelta = (l1 - l0) / (fMaxValue - fMinValue);
   Double_t newOff = (l0 - fMinValue) / (fMaxValue - fMinValue);

   if (newDelta < 0.001 || oldDelta < 0.001)
      return;

   for (Int_t pt = 1; pt < Int_t(fPalette->fNumPoints - 1); pt++)
      newPalette->fPoints[pt] = newOff +
            (fPalette->fPoints[pt] - fPalette->fPoints[1]) * newDelta / oldDelta;

   InsertNewPalette(newPalette);
   UpdateScreen(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Actually paint the paletter.

void TASPaletteEditor::PaintPalette::Paint(Option_t *)
{
   // get geometry of pad
   Int_t to_w = TMath::Abs(gPad->XtoPixel(gPad->GetX2()) -
                           gPad->XtoPixel(gPad->GetX1()));
   Int_t to_h = TMath::Abs(gPad->YtoPixel(gPad->GetY2()) -
                           gPad->YtoPixel(gPad->GetY1()));

   ASGradient grad;

   grad.npoints = (*fPalette)->fNumPoints - 2;
   grad.type = GRADIENT_Left2Right;
   grad.color = new ARGB32[grad.npoints];
   grad.offset = new double[grad.npoints];
   for (Int_t pt = 0; pt < grad.npoints; pt++) {
      grad.offset[pt] = ((*fPalette)->fPoints[pt + 1] - (*fPalette)->fPoints[1]) /
                        ((*fPalette)->fPoints[(*fPalette)->fNumPoints - 2] - (*fPalette)->fPoints[1]);
      grad.color[pt] = (((ARGB32)((*fPalette)->fColorBlue[pt + 1]   & 0xff00)) >>  8) |
                        (((ARGB32)((*fPalette)->fColorGreen[pt + 1] & 0xff00))      ) |
                        (((ARGB32)((*fPalette)->fColorRed[pt + 1]   & 0xff00)) <<  8) |
                        (((ARGB32)((*fPalette)->fColorAlpha[pt + 1] & 0xff00)) << 16);
   }

   ASImage * grad_im = make_gradient((ASVisual*)TASImage::GetVisual(), &grad , to_w, to_h,
                                     SCL_DO_COLOR, ASA_ARGB32, 0,
                                     fAttImage->GetImageQuality());
   delete [] grad.color;
   delete [] grad.offset;

   Window_t wid = (Window_t)gVirtualX->GetWindowID(gPad->GetPixmapID());
   TASImage::Image2Drawable(grad_im, wid, 0, 0);
   destroy_asimage(&grad_im);
}

////////////////////////////////////////////////////////////////////////////////
/// The blue limit line in the pixel value histogram.

TASPaletteEditor::LimitLine::LimitLine(Coord_t x, Coord_t y1, Coord_t y2,
                                       TASPaletteEditor *gui)
   : TLine(x, y1, x, y2)
{
   fGui = gui;
   SetLineColor(4);
   SetLineWidth(2);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the limit lines.

void TASPaletteEditor::LimitLine::Paint(Option_t *option)
{
   fY1 = gPad->GetUymin();
   fY2 = gPad->GetUymax();

   TLine::Paint(option);
}

////////////////////////////////////////////////////////////////////////////////

void TASPaletteEditor::LimitLine::ExecuteEvent(Int_t event,
                                               Int_t px, Int_t /*py*/)
{
   static Int_t oldX;

   if (!gPad) return;

   switch(event) {
      case kMouseMotion:
         gPad->SetCursor(kMove);
         break;

      case kButton1Down:
         gVirtualX->SetLineColor(-1);
         TAttLine::Modify();  //Change line attributes only if necessary
         oldX = gPad->XtoAbsPixel(fX1);
         break;

      case kButton1Motion:
         gVirtualX->DrawLine(oldX, gPad->YtoPixel(fY1), oldX, gPad->YtoPixel(fY2));
         oldX = px;
         gVirtualX->DrawLine(oldX, gPad->YtoPixel(fY1), oldX, gPad->YtoPixel(fY2));
         gVirtualX->Update();
         break;

      case kButton1Up:
         gVirtualX->SetLineColor(-1);
         TAttLine::Modify();  //Change line attributes only if necessary
         fX1 = fX2 = gPad->AbsPixeltoX(oldX);
         fGui->UpdateRange();
         gPad->Modified(kTRUE);
         gPad->Update();
         break;

      default:
         break;
   }
}
