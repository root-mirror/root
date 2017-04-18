// @(#)root/ged:$Id$
// Author: Marek Biskup, Ilka Antcheva   22/07/03

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGedMarkerSelect, TGedMarkerPopup                                    //
//                                                                      //
// The TGedMarkerPopup is a popup containing buttons to                 //
// select marker style.                                                 //
//                                                                      //
// The TGedMarkerSelect widget is a button showing selected marker      //
// and a little down arrow. When clicked on the arrow the               //
// TGedMarkerPopup pops up.                                             //
//                                                                      //
// Selecting a marker in this widget will generate the event:           //
// kC_MARKERSEL, kMAR_SELCHANGED, widget id, style.                     //
//                                                                      //
// and the signal:                                                      //
// MarkerSelected(Style_t marker)                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGedMarkerSelect.h"
#include "TGPicture.h"
#include "TGToolTip.h"
#include "TGButton.h"
#include "Riostream.h"

ClassImp(TGedMarkerSelect)
ClassImp(TGedMarkerPopup)

struct MarkerDescription_t {
   const char* fFilename;  // xpm file name
   const char* fName;      // type number for tooltip
   Int_t       fNumber;    // marker type number
};

static MarkerDescription_t  gMarkers[] = {

   {"marker1.xpm", "1", 1},
   {"marker6.xpm", "6", 6},
   {"marker7.xpm", "7", 7},
   {"marker2.xpm", "2", 2},
   {"marker3.xpm", "3", 3},
   {"marker4.xpm", "4", 4},
   {"marker5.xpm", "5", 5},
   {"marker20.xpm", "20", 20},
   {"marker21.xpm", "21", 21},
   {"marker22.xpm", "22", 22},
   {"marker23.xpm", "23", 23},
   {"marker24.xpm", "24", 24},
   {"marker25.xpm", "25", 25},
   {"marker26.xpm", "26", 26},
   {"marker27.xpm", "27", 27},
   {"marker28.xpm", "28", 28},
   {"marker29.xpm", "29", 29},
   {"marker30.xpm", "30", 30},
   {0, 0, 0},
};

static MarkerDescription_t* GetMarkerByNumber(Int_t number)
{
   for (Int_t i = 0; gMarkers[i].fFilename != 0; i++) {
      if (gMarkers[i].fNumber == number)
         return &gMarkers[i];
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create marker popup window.

TGedMarkerPopup::TGedMarkerPopup(const TGWindow *p, const TGWindow *m, Style_t markerStyle)
   : TGedPopup(p, m, 30, 30, kDoubleBorder | kRaisedFrame | kOwnBackground,
               GetDefaultFrameBackground())
{
   TGButton *b;
   fCurrentStyle = markerStyle;

   Pixel_t white;
   gClient->GetColorByName("white", white); // white background
   SetBackgroundColor(white);

   SetLayoutManager(new TGTileLayout(this, 1));

   for (int i = 0; gMarkers[i].fFilename != 0; i++) {
      AddFrame(b = new TGPictureButton(this, gMarkers[i].fFilename,
               gMarkers[i].fNumber, TGButton::GetDefaultGC()(), kSunkenFrame),
               new TGLayoutHints(kLHintsLeft, 14, 14, 14, 14));
      b->SetToolTipText(gMarkers[i].fName);
   }

   Resize(65, 94);
   MapSubwindows();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGedMarkerPopup::~TGedMarkerPopup()
{
   TGFrameElement *el;
   TIter next(GetList());
   while ((el = (TGFrameElement *)next())) {
      if (el->fFrame->InheritsFrom(TGPictureButton::Class()))
         fClient->FreePicture(((TGPictureButton *)el->fFrame)->GetPicture());
   }
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Process messages generated by the marker popup window.

Bool_t TGedMarkerPopup::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
   if (GET_MSG(msg) == kC_COMMAND && GET_SUBMSG(msg) == kCM_BUTTON) {
      SendMessage(fMsgWindow, MK_MSG(kC_MARKERSEL, kMAR_SELCHANGED), 0, parm1);
      EndPopup();
   }

   if (parm2) {
      // no warning
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create and show marker popup window.

TGedMarkerSelect::TGedMarkerSelect(const TGWindow *p, Style_t markerStyle, Int_t id)
   : TGedSelect(p, id)
{
   fPicture = 0;
   SetPopup(new TGedMarkerPopup(gClient->GetDefaultRoot(), this, markerStyle));
   SetMarkerStyle(markerStyle);

}

////////////////////////////////////////////////////////////////////////////////
/// Process messages according to the user input.

Bool_t TGedMarkerSelect::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
   if (GET_MSG(msg) == kC_MARKERSEL && GET_SUBMSG(msg) == kMAR_SELCHANGED) {
      SetMarkerStyle(parm2);
      parm1 = (Long_t)fWidgetId;
      SendMessage(fMsgWindow, MK_MSG(kC_MARKERSEL, kMAR_SELCHANGED),
                  parm1, parm2);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw selected marker type as current one.

void TGedMarkerSelect::DoRedraw()
{
   TGedSelect::DoRedraw();

   Int_t  x, y;
   UInt_t w, h;

   if (IsEnabled()) {
      // pattern rectangle

      x = fBorderWidth + 2;
      y = fBorderWidth + 2;  // 1;
      h = fHeight - (fBorderWidth * 2) - 4;  // -3;  // 14
      w = h;
      if (fState == kButtonDown) {
         ++x; ++y;
      }
      gVirtualX->DrawRectangle(fId, GetShadowGC()(), x, y, w - 1, h - 1);

      if (fPicture != 0) fPicture->Draw(fId, fDrawGC->GetGC(), x + 1, y + 1);
   } else { // sunken rectangle
      x = fBorderWidth + 2;
      y = fBorderWidth + 2;  // 1;
      w = 42;
      h = fHeight - (fBorderWidth * 2) - 4;  // 3;
      Draw3dRectangle(kSunkenFrame, x, y, w, h);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker.

void TGedMarkerSelect::SetMarkerStyle(Style_t markerStyle)
{
   fMarkerStyle = markerStyle;
   gClient->NeedRedraw(this);

   if (fPicture) {
      gClient->FreePicture(fPicture);
      fPicture = 0;
   }

   MarkerDescription_t *md = GetMarkerByNumber(fMarkerStyle);

   if (md)  fPicture = gClient->GetPicture(md->fFilename);

   MarkerSelected(fMarkerStyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Save the pattern select widget as a C++ statement(s) on output stream out

void TGedMarkerSelect::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   out <<"   TGedMarkerSelect *";
   out << GetName() << " = new TGedMarkerSelect(" << fParent->GetName()
       << "," << fMarkerStyle << "," << WidgetId() << ");" << std::endl;
}
