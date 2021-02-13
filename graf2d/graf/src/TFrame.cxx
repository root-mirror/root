// @(#)root/graf:$Id$
// Author: Rene Brun   31/10/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TVirtualPad.h"
#include "TFrame.h"
#include "TStyle.h"

ClassImp(TFrame);

/** \class TFrame
\ingroup BasicGraphics

Define a Frame.

A `TFrame` is a `TWbox` for drawing histogram frames.
*/

////////////////////////////////////////////////////////////////////////////////
/// Frame default constructor.

TFrame::TFrame(): TWbox()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Frame normal constructor.

TFrame::TFrame(Double_t x1, Double_t y1,Double_t x2, Double_t  y2)
       :TWbox(x1,y1,x2,y2)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Frame copy constructor.

TFrame::TFrame(const TFrame &frame) : TWbox(frame)
{
   ((TFrame&)frame).Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Frame default destructor.

TFrame::~TFrame()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this frame to frame.

void TFrame::Copy(TObject &frame) const
{
   TWbox::Copy(frame);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this frame with its current attributes.

void TFrame::Draw(Option_t *option)
{
   AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.
///
///  This member function is called when a TFrame object is clicked.

void TFrame::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (!gPad) return;

   if (!gPad->IsEditable()) return;

   TWbox::ExecuteEvent(event, px, py);

   Bool_t opaque  = gPad->OpaqueMoving();

   if ((event == kButton1Up) || ((opaque)&&(event == kButton1Motion))) {
      // update pad margins
      Double_t xmin         = gPad->GetUxmin();
      Double_t xmax         = gPad->GetUxmax();
      Double_t ymin         = gPad->GetUymin();
      Double_t ymax         = gPad->GetUymax();
      Double_t dx           = xmax-xmin;
      Double_t dy           = ymax-ymin;
      Double_t leftMargin   = (fX1-gPad->GetX1())/(gPad->GetX2()-gPad->GetX1());
      Double_t topMargin    = (gPad->GetY2()-fY2)/(gPad->GetY2()-gPad->GetY1());
      Double_t rightMargin  = (gPad->GetX2()-fX2)/(gPad->GetX2()-gPad->GetX1());
      Double_t bottomMargin = (fY1-gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());
      // margin may get very small negative values
      if (leftMargin   < 0) leftMargin   = 0;
      if (topMargin    < 0) topMargin    = 0;
      if (rightMargin  < 0) rightMargin  = 0;
      if (bottomMargin < 0) bottomMargin = 0;
      gPad->SetLeftMargin(leftMargin);
      gPad->SetRightMargin(rightMargin);
      gPad->SetBottomMargin(bottomMargin);
      gPad->SetTopMargin(topMargin);
      Double_t dxr  = dx/(1 - gPad->GetLeftMargin() - gPad->GetRightMargin());
      Double_t dyr  = dy/(1 - gPad->GetBottomMargin() - gPad->GetTopMargin());

      // Range() could change the size of the pad pixmap and therefore should
      // be called before the other paint routines
      gPad->Range(xmin - dxr*gPad->GetLeftMargin(),
                         ymin - dyr*gPad->GetBottomMargin(),
                         xmax + dxr*gPad->GetRightMargin(),
                         ymax + dyr*gPad->GetTopMargin());
      gPad->RangeAxis(xmin, ymin, xmax, ymax);
      fX1 = xmin;
      fY1 = ymin;
      fX2 = xmax;
      fY2 = ymax;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this wbox with its current attributes.

void TFrame::Paint(Option_t *option)
{
   const TPickerStackGuard stackGuard(this);

   if (!gPad->PadInHighlightMode() || (gPad->PadInHighlightMode() && this == gPad->GetSelected())) {
      TWbox::Paint(option);

      gPad->PaintBox(fX1,fY1,fX2,fY2,"s");
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Do not pop frame's, if allowed they would cover the picture they frame.

void TFrame::Pop()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TFrame::SavePrimitive(std::ostream &, Option_t * /*= ""*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Replace current frame attributes by current style.

void TFrame::UseCurrentStyle()
{
   if (gStyle->IsReading()) {
      SetFillColor(gStyle->GetFrameFillColor());
      SetLineColor(gStyle->GetFrameLineColor());
      SetFillStyle(gStyle->GetFrameFillStyle());
      SetLineStyle(gStyle->GetFrameLineStyle());
      SetLineWidth(gStyle->GetFrameLineWidth());
      SetBorderSize(gStyle->GetFrameBorderSize());
      SetBorderMode(gStyle->GetFrameBorderMode());
   } else {
      gStyle->SetFrameFillColor(GetFillColor());
      gStyle->SetFrameLineColor(GetLineColor());
      gStyle->SetFrameFillStyle(GetFillStyle());
      gStyle->SetFrameLineStyle(GetLineStyle());
      gStyle->SetFrameLineWidth(GetLineWidth());
      gStyle->SetFrameBorderSize(GetBorderSize());
      gStyle->SetFrameBorderMode(GetBorderMode());
   }
}
