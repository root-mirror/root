// @(#)root/gpad:$Id$
// Author: Rene Brun   23/11/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TSlider.h"
#include "TSliderBox.h"

#include <cstring>
#include <iostream>

ClassImp(TSlider);

/** \class TSlider
\ingroup gpad

A specialized TPad including a TSliderBox object.
The TSliderBox can be moved in the pad.

Slider drawing options include the possibility to change the slider
starting and ending positions or only one of them.

The current slider position can be retrieved via the functions
TSlider::GetMinimum and TSlider::GetMaximum
 These two functions return numbers in the range [0,1].

If a method has been set (via TSlider::SetMethod), the expression  is
executed via the interpreter when the button 1 is released.

  if no method has been set, and an object is referenced (TSlider::SetObject
  has been called), while the slider is being moved/resized,
  the object ExecuteEvent function is called.

### Example 1 using TSlider::SetMethod

#### macro xyslider.C

~~~ {.cpp}
void xyslider()
{
   // Example of macro featuring two sliders
   TFile *f = new TFile("hsimple.root");
   TH2F *hpxpy = (TH2F*)f->Get("hpxpy");
   TCanvas *c1 = new TCanvas("c1");
   TPad *pad = new TPad("pad","Color plot pad",0.1,0.1,0.98,0.98);
   pad->Draw();
   pad->cd();
   hpxpy->Draw("COLZ");
   c1->cd();

   // Create two sliders in main canvas. When button1 will be released
   // the macro action.C will be called.
   TSlider *xslider = new TSlider("xslider","x",0.1,0.02,0.98,0.08);
   xslider->SetMethod(".x action.C");
   TSlider *yslider = new TSlider("yslider","y",0.02,0.1,0.06,0.98);
   yslider->SetMethod(".x action.C");
}
~~~

#### macro action.C

~~~ {.cpp}
void action()
{
   Int_t nx = hpxpy->GetXaxis()->GetNbins();
   Int_t ny = hpxpy->GetYaxis()->GetNbins();
   Int_t binxmin = nx*xslider->GetMinimum();
   Int_t binxmax = nx*xslider->GetMaximum();
   hpxpy->GetXaxis()->SetRange(binxmin,binxmax);
   Int_t binymin = ny*yslider->GetMinimum();
   Int_t binymax = ny*yslider->GetMaximum();
   hpxpy->GetYaxis()->SetRange(binymin,binymax);
   pad->cd();
   hpxpy->Draw("COLZ");
   c1->Update();
}
~~~

The canvas and the sliders created in the above macro are shown in the picture
below.

\image html gpad_slider.png

### Example 2 using TSlider::SetObject    macro xyslider.C

  Same example as above. Instead of TSlider::SetMethod:
~~~ {.cpp}
    Myclass *obj = new Myclass(); // Myclass derived from TObject
    xslider->SetObject(obj);
    yslider->SetObject(obj);
~~~

When the slider will be changed, MyClass::ExecuteEvent will be called with px=0
and py = 0
*/

////////////////////////////////////////////////////////////////////////////////
/// slider default constructor.

TSlider::TSlider(): TPad()
{
   fObject  = 0;
   fMethod  = "";
   fMinimum = 0;
   fMaximum = 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Slider normal constructor.
///
///   x1,y1,x2,y2 are in pad user coordinates

TSlider::TSlider(const char *name, const char *title, Double_t x1, Double_t y1,Double_t x2, Double_t  y2, Color_t color, Short_t bordersize, Short_t bordermode)
           :TPad(name,title,0.1,0.1,0.9,0.9,color,bordersize,bordermode)
{
   if (!gPad) return;

   Double_t x1pad = gPad->GetX1();
   Double_t x2pad = gPad->GetX2();
   Double_t y1pad = gPad->GetY1();
   Double_t y2pad = gPad->GetY2();
   Double_t xmin  = (x1-x1pad)/(x2pad-x1pad);
   Double_t ymin  = (y1-y1pad)/(y2pad-y1pad);
   Double_t xmax  = (x2-x1pad)/(x2pad-x1pad);
   Double_t ymax  = (y2-y1pad)/(y2pad-y1pad);
   SetPad(xmin,ymin,xmax,ymax);
   Range(0,0,1,1);

   SetBit(kCanDelete);
   Modified(kTRUE);

   fMinimum = 0;
   fMaximum = 1;
   fObject  = 0;
   fMethod  = "";
   Double_t dx = PixeltoX(bordersize);
   Double_t dy = PixeltoY(-bordersize);
   TSliderBox *sbox = new TSliderBox(dx,dy,1-dx,1-dy,color,bordersize,-bordermode);
   sbox->SetSlider(this);
   fPrimitives->Add(sbox);
   AppendPad();
}

////////////////////////////////////////////////////////////////////////////////
/// slider default destructor.

TSlider::~TSlider()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this slider with its current attributes.

void TSlider::Paint(Option_t *option)
{
   TPad::Paint(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TSlider::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   TPad *padsav = (TPad*)gPad;
   char quote = '"';
   if (gROOT->ClassSaved(TSlider::Class())) {
      out<<"   ";
   } else {
      out<<"   TSlider *";
   }
   out<<"slider = new TSlider("<<quote<<GetName()<<quote<<", "<<quote<<GetTitle()
      <<quote
      <<","<<fXlowNDC
      <<","<<fYlowNDC
      <<","<<fXlowNDC+fWNDC
      <<","<<fYlowNDC+fHNDC
      <<");"<<std::endl;

   SaveFillAttributes(out,"slider",0,1001);
   SaveLineAttributes(out,"slider",1,1,1);

   if (GetBorderSize() != 2) {
      out<<"   slider->SetBorderSize("<<GetBorderSize()<<");"<<std::endl;
   }
   if (GetBorderMode() != -1) {
      out<<"   slider->SetBorderMode("<<GetBorderMode()<<");"<<std::endl;
   }
   Int_t lenMethod = strlen(GetMethod());
   if (lenMethod > 0) {
      out<<"   slider->SetMethod("<<quote<<GetMethod()<<quote<<");"<<std::endl;
   }

   out<<"   "<<padsav->GetName()<<"->cd();"<<std::endl;
   padsav->cd();
}

////////////////////////////////////////////////////////////////////////////////
/// Set Slider range in [0,1]

void TSlider::SetRange(Double_t xmin, Double_t xmax)
{
   TSliderBox *sbox = (TSliderBox*)fPrimitives->FindObject("TSliderBox");
   if (sbox) {
      if (fAbsWNDC > fAbsHNDC) {
         sbox->SetX1(xmin);
         sbox->SetX2(xmax);
      } else {
         sbox->SetY1(xmin);
         sbox->SetY2(xmax);
      }
   }
   fMinimum = xmin;
   fMaximum = xmax;
   Modified();
}
