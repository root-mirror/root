// @(#)root/guibuilder:$Id$
// Author: Valeriy Onuchin   12/09/04

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGuiBldHintsButton.h"
#include "TGResourcePool.h"
#include "TVirtualX.h"

//_____________________________________________________________________________
//
// TGuiBldHintsButton
//
// Special button class used for editing layout hints in the ROOT GUI Builder.
//_____________________________________________________________________________

ClassImp(TGuiBldHintsButton);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGuiBldHintsButton::TGuiBldHintsButton(const TGWindow* p, Int_t id) :
                     TGButton(p, id)
{
   fStayDown = kTRUE;

   switch ((ELayoutHints)fWidgetId) {
      case kLHintsCenterX:
      case kLHintsExpandX:
         Resize(40, 15);
         break;
      case kLHintsCenterY:
      case kLHintsExpandY:
         Resize(15, 40);
         break;
      default:
         Resize(15, 15);
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw button.

void TGuiBldHintsButton::DoRedraw()
{
   TGButton::DoRedraw();

   switch (fWidgetId) {
      case kLHintsCenterX:
         DrawCenterX();
         break;
      case kLHintsCenterY:
         DrawCenterY();
         break;
      case kLHintsExpandX:
         DrawExpandX();
         break;
      case kLHintsExpandY:
         DrawExpandY();
         break;
      case (kLHintsTop | kLHintsLeft):
         DrawTopLeft();
         break;
      case (kLHintsTop | kLHintsRight):
         DrawTopRight();
         break;
      case (kLHintsBottom | kLHintsLeft):
         DrawBottomLeft();
         break;
      case (kLHintsBottom | kLHintsRight):
         DrawBottomRight();
         break;
      default:
         DrawExpandX();
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw expand X button.

void TGuiBldHintsButton::DrawExpandX()
{
   const int dist = 3;
   const int amplitude = TMath::Min(3, (int)fHeight/3);
   int base = fHeight/2;
   int i = 0;
   const TGResourcePool *pool = fClient->GetResourcePool();
   const TGGC* gc = pool->GetWhiteGC();

   if ((fState == kButtonDown) || (fState == kButtonEngaged)) {
      base++;
   }

   for ( i = 1; i < (int)fWidth/3 - 2; ++i ) {
      gVirtualX->DrawLine(fId, gc->GetGC(), i * dist, base - amplitude,
                           i * dist + dist/2, base + amplitude);
   }
   gc = IsEnabled() ? pool->GetSelectedBckgndGC() : pool->GetFrameShadowGC();

   for ( i = 1; i < (int)fWidth/3 - 2; ++i ) {
      gVirtualX->DrawLine(fId, gc->GetGC(), i * dist + dist/2, base + amplitude,
                           i * dist + dist, base - amplitude);
   }
   gVirtualX->DrawLine(fId, gc->GetGC(), 3, 6, 3, fHeight - 6);
   gVirtualX->DrawLine(fId, gc->GetGC(), fWidth - 6, 6, fWidth - 6, fHeight - 6);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw expand Y button.

void TGuiBldHintsButton::DrawExpandY()
{
   const int dist = 3;
   const int amplitude = TMath::Min(3, (int)fWidth/3);
   int base = fWidth/2;
   int i = 0;

   if ((fState == kButtonDown) || (fState == kButtonEngaged)) {
      base++;
   }
   const TGResourcePool *pool = fClient->GetResourcePool();
   const TGGC* gc = pool->GetWhiteGC();

   for ( i = 1; i < (int)fHeight/3 - 2; ++i ) {
      gVirtualX->DrawLine(fId, gc->GetGC(), base - amplitude, i * dist,
                           base + amplitude,i * dist + dist/2);
   }

   gc = IsEnabled() ? pool->GetSelectedBckgndGC() : pool->GetFrameShadowGC();

   for ( i = 1; i < (int)fHeight/3 - 2; ++i ) {
      gVirtualX->DrawLine(fId, gc->GetGC(), base + amplitude, i * dist + dist/2,
                           base - amplitude, i * dist + dist );
   }
   gVirtualX->DrawLine(fId, gc->GetGC(), 6, 3, fWidth - 6, 3);
   gVirtualX->DrawLine(fId, gc->GetGC(), 6, fHeight - 6, fWidth - 6, fHeight - 6);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw center X buton.

void TGuiBldHintsButton::DrawCenterX()
{
   int base = fHeight/2;
   int x = 6;
   int y = 6;

   const TGResourcePool *pool = fClient->GetResourcePool();
   const TGGC* gc = pool->GetWhiteGC();

   if ((fState == kButtonDown) || (fState == kButtonEngaged)) {
      base++;
      x++;
      y++;
   }

   gVirtualX->DrawLine(fId, gc->GetGC(), x, base, x + fWidth - 12, base);

   gc = IsEnabled() ? pool->GetSelectedBckgndGC() : pool->GetFrameShadowGC();

   gVirtualX->DrawLine(fId, gc->GetGC(), x, base - 1, x + fWidth/2 - 12, base - 1);
   gVirtualX->DrawLine(fId, gc->GetGC(), x + fWidth/2, base - 1, x + fWidth - 12, base - 1);
   gVirtualX->DrawLine(fId, gc->GetGC(), x, base + 1, x + fWidth/2 - 12, base +  1);
   gVirtualX->DrawLine(fId, gc->GetGC(), x + fWidth/2, base + 1, x + fWidth - 12, base + 1);

   Point_t arrow[3];
   arrow[0].fX = arrow[01].fX = x + fWidth/2 - 12;
   arrow[2].fX = x + fWidth/2 - 6;
   arrow[2].fY = y + fHeight/2 - 6;
   arrow[0].fY = arrow[2].fY - 4;
   arrow[1].fY = arrow[2].fY + 4;
   gVirtualX->FillPolygon(fId, gc->GetGC(), (Point_t*)&arrow, 3);

   arrow[0].fX = arrow[01].fX = x + fWidth/2;
   gVirtualX->FillPolygon(fId, gc->GetGC(), (Point_t*)&arrow, 3);

   gVirtualX->DrawLine(fId, gc->GetGC(), x, y, x, y + fHeight - 12);
   gVirtualX->DrawLine(fId, gc->GetGC(), x + fWidth - 12, y, x + fWidth - 12, y + fHeight - 12);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw center Y button.

void TGuiBldHintsButton::DrawCenterY()
{
   int base = fWidth/2;
   int x = 6;
   int y = 6;

   const TGResourcePool *pool = fClient->GetResourcePool();
   const TGGC* gc = pool->GetWhiteGC();

   if ((fState == kButtonDown) || (fState == kButtonEngaged)) {
      base++;
      x++;
      y++;
   }

   gVirtualX->DrawLine(fId, gc->GetGC(), base, y, base, y + fHeight - 12);

   gc = IsEnabled() ? pool->GetSelectedBckgndGC() : pool->GetFrameShadowGC();

   gVirtualX->DrawLine(fId, gc->GetGC(), base - 1, y,  base - 1, y + fHeight/2 - 12);
   gVirtualX->DrawLine(fId, gc->GetGC(),  base - 1, y + fHeight/2, base - 1, y + fHeight - 12);
   gVirtualX->DrawLine(fId, gc->GetGC(), base + 1, y,  base + 1, y + fHeight/2 - 12);
   gVirtualX->DrawLine(fId, gc->GetGC(),  base + 1, y + fHeight/2, base + 1, y + fHeight - 12);

   Point_t arrow[3];
   arrow[0].fY = arrow[01].fY = y + fHeight/2 - 12;
   arrow[2].fY = y + fHeight/2 - 6;
   arrow[2].fX = x + fWidth/2 - 6;
   arrow[0].fX = arrow[2].fX - 4;
   arrow[1].fX = arrow[2].fX + 4;
   gVirtualX->FillPolygon(fId, gc->GetGC(), (Point_t*)&arrow, 3);

   arrow[0].fY = arrow[01].fY = y + fHeight/2;
   gVirtualX->FillPolygon(fId, gc->GetGC(), (Point_t*)&arrow, 3);
   gVirtualX->DrawLine(fId, gc->GetGC(), x, y, x + fWidth - 12, y);
   gVirtualX->DrawLine(fId, gc->GetGC(), x, y + fHeight - 12,  x + fWidth - 12, y + fHeight - 12);
}

////////////////////////////////////////////////////////////////////////////////
/// DrawTopLeft.

void TGuiBldHintsButton::DrawTopLeft()
{
}

////////////////////////////////////////////////////////////////////////////////
/// DrawTopRight.

void TGuiBldHintsButton::DrawTopRight()
{
}

////////////////////////////////////////////////////////////////////////////////
/// DrawBottomLeft.

void TGuiBldHintsButton::DrawBottomLeft()
{
}

////////////////////////////////////////////////////////////////////////////////
/// DrawBottomRight.

void TGuiBldHintsButton::DrawBottomRight()
{
}

