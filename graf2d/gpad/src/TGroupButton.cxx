// @(#)root/gpad:$Id: faa839c41a9482cbcb403b991070b15d606c9137 $
// Author: Rene Brun   01/07/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TBox.h"
#include "TGroupButton.h"
#include "TDialogCanvas.h"
#include "TCanvas.h"
#include "TText.h"
#include "TInterpreter.h"
#include "strlcpy.h"
#include "snprintf.h"

#include <cstring>
#include <iostream>

ClassImp(TGroupButton);


/** \class TGroupButton
\ingroup gpad

A specialized TButton used in a group of Buttons.
When a button from a group of TGroupButtons is selected, all other buttons
from the group with the same name are disabled.

For examples of use of TGroupButton objects, see:
AttFillCanvas, TAttLineCanvas, TAttTextCanvas and TAttMarkerCanvas.

A TGroupButton object is a specialized TPad including possible list
of primitives used to build selections and options menus in a canvas.
*/


////////////////////////////////////////////////////////////////////////////////
/// GroupButton default constructor.

TGroupButton::TGroupButton(): TButton()
{
   SetFraming();
}

////////////////////////////////////////////////////////////////////////////////
/// GroupButton normal constructor.

TGroupButton::TGroupButton(const char *groupname, const char *title, const char *method, Double_t x1, Double_t y1,Double_t x2, Double_t  y2)
           :TButton(title,method,x1,y1,x2,y2)
{
   SetName((char*)groupname);
   SetFraming();
}

////////////////////////////////////////////////////////////////////////////////
/// GroupButton default destructor.

TGroupButton::~TGroupButton()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Display Color Table in an attribute canvas.

void TGroupButton::DisplayColorTable(const char *action, Double_t x0, Double_t y0, Double_t wc, Double_t hc)
{
   TGroupButton *colorpad;
   Int_t i, j;
   Int_t color;
   Double_t xlow, ylow, hs, ws;

   // draw color table buttons
   hs = hc/5;
   ws = wc/10;
   char command[32];
   for (i=0;i<10;i++) {
      xlow = x0 + ws*i;
      for (j=0;j<5;j++) {
         ylow = y0 + hs*j;
         color = 10*j + i + 1;
         snprintf(command,32,"%s(%d)",action,10*j+i+1);
         colorpad = new TGroupButton("Color","",command,xlow, ylow, xlow+0.9*ws, ylow+0.9*hs);
         colorpad->SetFillColor(color);
         colorpad->SetBorderSize(1);
         if (i == 0 && j == 0) colorpad->SetBorderMode(-1);
         colorpad->Draw();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action of this button.
///
/// If an object has been selected before executing the APPLY button
/// in the control canvas, The member function and its parameters
/// for this object is executed via the interpreter.

void TGroupButton::ExecuteAction()
{
   TVirtualPad *pad;
   char line[128];
   strlcpy(line,GetMethod(),128);
   char *method = line;
   if(!line[0]) return;
   char *params = strchr(method,'(');
   if (params) {
      *params = 0;
      params++;
      char *end = strrchr(params,')');
      if (end) *end = 0;
   }
   TDialogCanvas *canvas = (TDialogCanvas*)GetMother();
   TObject *obj = canvas->GetRefObject();
   if (!obj) return;
   if (strcmp(method,"PIXELS")) {
      obj->Execute(method,params);
   } else {
      TText *text = (TText*)GetListOfPrimitives()->First();
      Int_t npixels = Int_t((YtoPixel(0) - YtoPixel(1))*text->GetTextSize());
      Double_t dy;
      pad = gROOT->GetSelectedPad();
      if (!params) return;
      Int_t nmax = (Int_t)(params-method);
      if (obj->InheritsFrom("TPaveLabel")) {
         TBox *pl = (TBox*)obj;
         dy = pad->AbsPixeltoY(0) - pad->AbsPixeltoY(npixels);
         snprintf(params,nmax,"%f",dy/(pl->GetY2() - pl->GetY1()));
         obj->Execute("SetTextSize",params);
      } else {
         if (obj->InheritsFrom("TPave")) {
            dy = pad->AbsPixeltoY(0) - pad->AbsPixeltoY(npixels);
            snprintf(params,nmax,"%f",dy/(pad->GetY2() - pad->GetY1()));
            obj->Execute("SetTextSize",params);
         } else {
            snprintf(params,nmax,"%d",npixels);
            obj->Execute("SetTextSizePixels",params);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.
///
/// This member function is called when a Button object is clicked.

void TGroupButton::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (fMother->IsEditable()) {
      TPad::ExecuteEvent(event,px,py);
      return;
   }

   TCanvas *c = gPad->GetCanvas();
   if (!c) return;
   TIter next(c->GetListOfPrimitives());
   TObject *obj;
   TGroupButton *button;
   TPad *pad;
   TDialogCanvas *canvas;

   switch (event) {

   case kButton1Down:

   case kMouseMotion:

      break;

   case kButton1Motion:

      break;

   case kButton1Up:
      //Clicked on APPLY button?
      if (!strcasecmp(GetName(),"APPLY")) {
         canvas = (TDialogCanvas*)GetMother();
         if (!strcasecmp(GetTitle(),"CLOSE")) {
            canvas->Close();
            return;
         }
         pad = canvas->GetRefPad();
         if (pad) pad->GetCanvas()->FeedbackMode(kFALSE);
         canvas->Apply(GetTitle());   //just in case the apply button executes some code
         if (pad) {
            pad->Modified(kTRUE);
            pad->Update();
         }
         break;
      }
      //Unset other buttons with same name
      while ((obj = next())) {
         if (obj == this) continue;
         if (obj->InheritsFrom(TGroupButton::Class())) {
            button = (TGroupButton*)obj;
            if (!strcmp(button->GetName(),GetName())) {
               if (button->GetBorderMode() < 0) {
                  button->SetBorderMode(1);
                  button->Modified();
               }
            }
         }
      }
      //Set button on
      SetBorderMode(-1);
      Modified();
      c->Modified();
      gPad->Update();
      break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TGroupButton::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   TPad *padsav = (TPad*)gPad;
   char quote = '"';
   if (gROOT->ClassSaved(TGroupButton::Class())) {
      out<<"   ";
   } else {
      out<<"   TGroupButton *";
   }
   out<<"button = new TGroupButton("<<quote<<GetName()<<quote<<", "<<quote<<GetTitle()
      <<quote<<","<<quote<<GetMethod()<<quote
      <<","<<fXlowNDC
      <<","<<fYlowNDC
      <<","<<fXlowNDC+fWNDC
      <<","<<fYlowNDC+fHNDC
      <<");"<<std::endl;

   SaveFillAttributes(out,"button",0,1001);
   SaveLineAttributes(out,"button",1,1,1);
   SaveTextAttributes(out,"button",22,0,1,62,.75);

   if (GetBorderSize() != 2) {
      out<<"   button->SetBorderSize("<<GetBorderSize()<<");"<<std::endl;
   }
   if (GetBorderMode() != 1) {
      out<<"   button->SetBorderMode("<<GetBorderMode()<<");"<<std::endl;
   }

   out<<"   button->Draw();"<<std::endl;
   out<<"   button->cd();"<<std::endl;

   TIter next(GetListOfPrimitives());
   TObject *obj = next();  //do not save first primitive

   while ((obj = next()))
         obj->SavePrimitive(out, (Option_t *)next.GetOption());

   out<<"   "<<padsav->GetName()<<"->cd();"<<std::endl;
   padsav->cd();
}
