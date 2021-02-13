// @(#)root/gpad:$Id$
// Author: Rene Brun   03/07/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TDialogCanvas.h"
#include "TGroupButton.h"
#include "TStyle.h"

ClassImp(TDialogCanvas);

/** \class TDialogCanvas
\ingroup gpad

A canvas specialized to set attributes.

It contains, in general, TGroupButton objects.
When the APPLY button is executed, the actions corresponding
to the active buttons are executed via the Interpreter.

See examples in TAttLineCanvas, TAttFillCanvas, TAttTextCanvas, TAttMarkerCanvas
*/

////////////////////////////////////////////////////////////////////////////////
/// DialogCanvas default constructor

TDialogCanvas::TDialogCanvas() : TCanvas()
{
   fRefObject = 0;
   fRefPad    = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// DialogCanvas constructor

TDialogCanvas::TDialogCanvas(const char *name, const char *title, Int_t ww, Int_t wh)
             : TCanvas(name,title,-ww,wh)
{
   SetFillColor(36);
   fRefObject = 0;
   fRefPad    = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// DialogCanvas constructor

TDialogCanvas::TDialogCanvas(const char *name, const char *title, Int_t wtopx, Int_t wtopy, UInt_t ww, UInt_t wh)
             : TCanvas(name,title,-wtopx,wtopy,ww,wh)
{
   SetFillColor(36);
   fRefObject = 0;
   fRefPad    = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// DialogCanvas default destructor

TDialogCanvas::~TDialogCanvas()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Called when the APPLY button is executed

void TDialogCanvas::Apply(const char *action)
{
   if (!fRefPad) return;
   SetCursor(kWatch);

   TIter next(fPrimitives);
   TObject *refobj = fRefObject;
   TObject *obj;
   TGroupButton *button;
   if (!strcmp(action,"gStyle")) fRefObject = gStyle;

   while ((obj = next())) {
      if (obj->InheritsFrom(TGroupButton::Class())) {
         button = (TGroupButton*)obj;
         if (button->GetBorderMode() < 0) button->ExecuteAction();
      }
   }
   fRefObject = refobj;
   if (!gROOT->GetSelectedPad()) return;
   gROOT->GetSelectedPad()->Modified();
   gROOT->GetSelectedPad()->Update();
}


////////////////////////////////////////////////////////////////////////////////
/// Create APPLY, gStyle and CLOSE buttons

void TDialogCanvas::BuildStandardButtons()
{
   TGroupButton *apply = new TGroupButton("APPLY","Apply","",.05,.01,.3,.09);
   apply->SetTextSize(0.55);
   apply->SetBorderSize(3);
   apply->SetFillColor(44);
   apply->Draw();

   apply = new TGroupButton("APPLY","gStyle","",.375,.01,.625,.09);
   apply->SetTextSize(0.55);
   apply->SetBorderSize(3);
   apply->SetFillColor(44);
   apply->Draw();

   apply = new TGroupButton("APPLY","Close","",.70,.01,.95,.09);
   apply->SetTextSize(0.55);
   apply->SetBorderSize(3);
   apply->SetFillColor(44);
   apply->Draw();
}


////////////////////////////////////////////////////////////////////////////////
/// Set world coordinate system for the pad

void TDialogCanvas::Range(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
   TPad::Range(x1,y1,x2,y2);
}


////////////////////////////////////////////////////////////////////////////////
/// Recursively remove object from a pad and its sub-pads

void TDialogCanvas::RecursiveRemove(TObject *obj)
{
   TPad::RecursiveRemove(obj);
   if (fRefObject == obj) fRefObject = 0;
   if (fRefPad    == obj) fRefPad    = 0;
}
