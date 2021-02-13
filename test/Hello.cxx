// @(#)root/test:$Id$
// Author: Valeriy Onuchin & Fons Rademakers   04/10/98

///////////////////////////////////////////////////////////////////
//  Animated Text with cool wave effect.
//
//  ROOT implementation of the hello world example borrowed
//  from the Qt hello world example.
//
//  To run this example do the following:
//  $ root
//  root [0] gSystem.Load("Hello")
//  root [1] Hello h
//  <enjoy>
//  root [2] .q
//
//  Other ROOT fun examples: Tetris, Aclock ...
//
//  Begin_Html
// <img src="http://emcal06.rhic.bnl.gov/~onuchin/root/gif/hello_clock.gif">
//  End_Html
//
///////////////////////////////////////////////////////////////////

#include <TSystem.h>
#include "Hello.h"
#include "TList.h"

ClassImp(Hello);


TChar::TChar(char ch, Coord_t x, Coord_t y) : TText(x, y, "")
{
   // Create single character text.

   SetTitle(Form("%c",ch));
}

Float_t TChar::GetWidth()
{
   // Return character width

   UInt_t w,h;

   if (!TVirtualPad::Pad()) return 0;

   Float_t wh = (Float_t)gPad->XtoAbsPixel(gPad->GetX2());
   GetTextExtent(w, h, (char*)GetTitle());
   return  w/wh;
}


Hello::Hello(const char *text) : TTimer(40, kTRUE)
{
   // Hello constructor

   TChar *ch;
   fI = 0;
   fList = new TList();

   if (!TVirtualPad::Pad())
      fPad = new TCanvas("Hello:canvas","ROOT says Hello!",200,200,400,200);
   else
      fPad = (TPad*)gPad;

   Connect(fPad, "Closed()", "TTimer", this, "TurnOff()");

   while(text[fI]) {               // create  list of characters
      ch = new TChar(text[fI]);
      ch->SetTextFont(72);         // times-bold-r-normal
      ch->SetTextSize(0.3);
      ch->Modify();
      fI++;
      fList->Add(ch);
   }

   Draw();                         // append this to current pad
   Paint();                        // calculate new coordinates of chars and paint it
   gSystem->AddTimer(this);        // start timer= start animation
}

Hello::~Hello()
{
   // Clean up hello.

   fList->Delete();
   delete fList;
}

void Hello::ExecuteEvent(Int_t event, Int_t, Int_t)
{
   // Actions when mouse clicked.

   if (event == kButton1Up) {
      TTimer::Remove();
   }
}

Float_t Hello::GetWidth()
{
   // Return text width

   TChar *ch;
   TIter nextin(fList);
   Float_t width = 0;

   while ((ch = (TChar*)nextin()))
      width = width + ch->GetWidth();

   return width;
}

void Hello::Paint(Option_t *)
{
   // Paint text

   static int sin_tbl[16] = {
        0, 38, 71, 92, 100, 92, 71, 38, 0, -38, -71, -92, -100, -92, -71, -38
   };

   TChar *ch;
   TIter nextin(fList);
   Float_t width = GetWidth();

   Coord_t xnext = (1-width)/2.;    // draw text in the center of the pad
   Coord_t y;
   Coord_t y0 = 0.5;

   fI = (fI+1) & 15;
   int i = 0;

   while ((ch = (TChar*)nextin())) {
      int i16 = (fI+i) & 15;
      y = y0 - (Coord_t)sin_tbl[i16]/1500.;   // sin-wave coordinates
      ch->SetY(y);
      ch->SetX(xnext);
      ch->SetTextColor(i16+200);              // let's use x3d colors
      i++;
      xnext = xnext + ch->GetWidth();         // x position of next character
      ch->Paint("");                          // paint character
   }
}

void Hello::Print(Option_t * opt) const
{
   fList->Print(opt);
}

void Hello::ls(Option_t *opt) const
{
   fList->ls(opt);
}

Bool_t Hello::Notify()
{
   // Actions after timer time-out

   fPad->Modified();
   fPad->Update();
   TTimer::Reset();
   return kFALSE;
}
