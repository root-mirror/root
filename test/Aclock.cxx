// @(#)root/test:$Id$
// Author: Valeriy Onuchin & Fons Rademakers   04/10/98

///////////////////////////////////////////////////////////////////
//  ROOT implementation of the X11 xclock.
//
//  To run this example do the following:
//
//  $ root
//  root [0] gSystem->Load("libGpad")
//  root [1] gSystem->Load("Aclock")
//  root [1] Aclock a
//  <enjoy>
//  root [2] .q
//
//  Other ROOT fun examples: Tetris, Hello ...
//
///////////////////////////////////////////////////////////////////

#include <TSystem.h>
#include "Aclock.h"


ClassImp(Aclock);


Float_t MinuteHand::fgMinuteHandX[] = { -0.05, 0, 0.05 };
Float_t MinuteHand::fgMinuteHandY[] = { -0.04, 0.625, -0.04 };

Float_t HourHand::fgHourHandX[] = { -0.05, 0, 0.05 };
Float_t HourHand::fgHourHandY[] = { -0.04, 0.4, -0.04 };

Float_t SecondHand::fgSecondHandX[] = { 0., 0.035, 0., -0.035 };
Float_t SecondHand::fgSecondHandY[] = { 0.78, 0.73, 0.68, 0.73 };

TDatime *ClockHand::fgTime = new TDatime;


TPolygon::TPolygon(Int_t n, Float_t *x, Float_t *y) : TPolyLine(n+1,x,y), fPad(0)
{
   // Create filled polygon. Polygon will be added to current pad.

   if (!TVirtualPad::Pad()) { Error("Constructor","Create TPad first!"); return; }

   fPad  = (TPad*)TVirtualPad::Pad();
   fX[n] = x[0];             // create an extra point to connect polyline ends
   fY[n] = y[0];
   Draw();                   // append to fPad
}

void TPolygon::Paint(Option_t *)
{
   // Paint filled polygon.

   TAttLine::Modify();
   TAttFill::Modify();
   fPad->PaintFillArea(fN-1, fX, fY);
   fPad->PaintPolyLine(fN, fX, fY, "");
}



ClockHand::ClockHand(Int_t n, Float_t *x, Float_t *y) : TPolygon(n,x,y)
{
   // Create clockhand

   fPrevTimeValue = 0;
   fX0 = new Float_t[GetN()];   // initial shape and position
   fY0 = new Float_t[GetN()];   // initial shape and position

   for (int i = 0; i < GetN(); i++) {
     fX0[i] = fX[i];
     fY0[i] = fY[i];
   }
   SetFillColor(9);            // default fill color
}

void ClockHand::Update()
{
   // Update hand position

   if (!IsModified()) return;

   fPrevTimeValue = GetTimeValue();
   Move(GetHandAngle());
}

void ClockHand::Move(Float_t clock_angle)
{
   // Move clockhand.

  Int_t n = GetN();

  // ClockPoints used to rotate,scale and shift initial points
  static ClockPoints *points = new ClockPoints();

  Float_t wh = (Float_t)fPad->XtoPixel(fPad->GetX2());
  Float_t hh = (Float_t)fPad->YtoPixel(fPad->GetY1());

  for (int i = 0; i < n; i++) {
     points->SetXY(fX0[i],fY0[i]);
     points->Rotate(clock_angle);

     if (wh < hh) points->Scale(0.5,0.5*wh/hh);    //  scale to have a circular clock
     else         points->Scale(0.5*hh/wh,0.5);

     points->Shift(0.5,0.5);                      // move to the center of pad
     fX[i] = points->GetX();
     fY[i] = points->GetY();
   }
}


Aclock::Aclock(Int_t csize) : TTimer(500, kTRUE)
{
   // Create a clock in a new canvas.

   fPad = new TCanvas("Aclock:canvas","xclock",-csize,csize);

   fPad->SetFillColor(14);     // grey

   fMinuteHand = new MinuteHand();
   fSecondHand = new SecondHand();
   fHourHand   = new HourHand();

   SetBit(kCanDelete);
   Draw();                         // append this Aclock to fPad
   Animate();
   gSystem->AddTimer(this);        // start timer = start animation
}

Aclock::~Aclock()
{
   // Clean up the clock.

   delete fMinuteHand;
   delete fSecondHand;
   delete fHourHand;
}

void Aclock::Paint(Option_t *)
{
   // Just draw clock scale (time and minutes ticks)

   static ClockPoints *point1 = new ClockPoints();
   static ClockPoints *point2 = new ClockPoints();

   Float_t wh = (Float_t)fPad->XtoPixel(fPad->GetX2());
   Float_t hh = (Float_t)fPad->YtoPixel(fPad->GetY1());

   for (int i = 0; i < 60; i++ ) {             // draw minute/hour ticks
      point1->SetXY(0.,0.9);

      if (!(i%5)) point2->SetXY(0.,0.8);       // hour  ticks  are longer
      else        point2->SetXY(0.,0.87);

      Float_t angle = 6.*i;
      point1->Rotate(angle);
      point2->Rotate(angle);

      if (wh < hh) {                   // scale in oder to draw circle scale
         point1->Scale(0.5,0.5*wh/hh);
         point2->Scale(0.5,0.5*wh/hh);
      } else {
         point1->Scale(0.5*hh/wh,0.5);
         point2->Scale(0.5*hh/wh,0.5);
      }

      point1->Shift(0.5,0.5);              // move to center of pad
      point2->Shift(0.5,0.5);

      fPad->PaintLine(point1->GetX(),point1->GetY(),point2->GetX(),point2->GetY());
   }
}

void Aclock::Animate()
{
   // Update clock hand positions and redraw the fPad

   if (!fSecondHand->IsModified()) return;

   fHourHand->Update();            //  update position every minute
   fMinuteHand->Update();          //  update position every minute
   fSecondHand->Update();          //  update position every second
   fPad->Modified();               //  drawing ...
   fPad->Update();
}

Bool_t Aclock::Notify()
{
   // Actions after timer's time-out

   Animate();
   TTimer::Reset();
   return kFALSE;
}
