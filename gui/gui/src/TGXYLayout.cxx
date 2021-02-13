// @(#)root/gui:$Id$
// Author: Reiner Rohlfs   24/03/2002

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun, Fons Rademakers and Reiner Rohlfs *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGXYLayout                                                           //
//                                                                      //
// Is a layout manager where the position and the size of each widget   //
// in the frame are defined by X / Y - coordinates. The coordinates     //
// for each widget are defined by the TGXYLayoutHints. Therefore it     //
// is not possible to share a layout hint for several widgets.          //
//                                                                      //
// The coordinates (X, Y) and the size (W, H) are defined in units      //
// of the size of a typical character. Also the size of the             //
// TGCompositeFrame for which a TGXYLayout manager is used has to be    //
// defined in its constructor in units of the size of a character!      //
//                                                                      //
// It is not possible to use any other layout hint than the             //
// TGXYLayoutHints for this layout manager!                             //
//                                                                      //
// The rubberFlag in the constructor of the TGLXYLayoutHins defines     //
// how the position and the size of a widget is recalculated if the     //
// size of the frame is increased:                                      //
// - kLRubberX: The X - position (left edge) is increased by the same   //
//              factor as the width of the frame increases.             //
// - kLRubberY: The Y - position (upper edge) is increased by the same  //
//              factor as the height of the frame increases.            //
// - kLRubberW: The width of the widget is increased by the same        //
//              factor as the width of the frame increases.             //
// - kLRubberH: The height of the widget is increased by the same       //
//              factor as the height of the frame increases.            //
// But the size never becomes smaller than defined by the               //
// TGXYLayoutHints and the X and Y coordinates becomes never smaller    //
// than defined by the layout hints.                                    //
//                                                                      //
// TGXYLayoutHints                                                      //
//                                                                      //
// This layout hint must be used for the TGXYLouyout manager!           //
//                                                                      //
//                                                                      //
// Example how to use this layout manager:                              //
//                                                                      //
// TGMyFrame::TGMyFrame()                                               //
//    : TGMainFrame(gClient->GetRoot(), 30, 12)                         //
//    // frame is 30 character long and 12 character heigh              //
// {                                                                    //
//    SetLayoutManager(new TGXYLayout(this));                           //
//                                                                      //
//    // create a button of size 8 X 1.8 at position 20 / 1             //
//    TGTextButton * button;                                            //
//    button = new TGTextButton(this, "&Apply", 1);                     //
//    AddFrame(button, new TGXYLayoutHints(20, 1, 8, 1.8));             //
//                                                                      //
//    // create a listbox of size 18 X 10 at position 1 / 1.            //
//    // The height will increase if the frame height increases         //
//    TGListBox * listBox;                                              //
//    listBox = new TGListBox(this, 2);                                 //
//    AddFrame(listBox, new TGXYLayoutHints(1, 1, 18, 10,               //
//             TGXYLayoutHints::kLRubberX |                             //
//             TGXYLayoutHints::kLRubberY |                             //
//             TGXYLayoutHints::kLRubberH ));                           //
//    .                                                                 //
//    .                                                                 //
//    .                                                                 //
// }                                                                    //
//                                                                      //
// Normally there is one layout hint per widget. Therefore these        //
// can be deleted like in the following example in the desctuctor       //
// of the frame:                                                        //
//                                                                      //
// TGMyFrame::~TGMyFrame()                                              //
// {                                                                    //
//    // Destructor, deletes all frames and their layout hints.         //
//                                                                      //
//    TGFrameElement *ptr;                                              //
//                                                                      //
//    // delete all frames and layout hints                             //
//    if (fList) {                                                      //
//       TIter next(fList);                                             //
//       while ((ptr = (TGFrameElement *) next())) {                    //
//          if (ptr->fLayout)                                           //
//             delete ptr->fLayout;                                     //
//          if (ptr->fFrame)                                            //
//             delete ptr->fFrame;                                      //
//       }                                                              //
//    }                                                                 //
// }                                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGXYLayout.h"
#include "TGFrame.h"
#include "TGLabel.h"
#include "TVirtualX.h"

#include <iostream>


ClassImp(TGXYLayout);
ClassImp(TGXYLayoutHints);

////////////////////////////////////////////////////////////////////////////////
/// Constructor. The x, y, w and h define the position of the widget in
/// its frame and the size of the widget. The unit is the size of a
/// character. The rubberFlag defines how to move and to resize the
/// widget when the frame is resized. Default is moving the X and Y
/// position but keep the size of the widget.

TGXYLayoutHints::TGXYLayoutHints(Double_t x, Double_t y, Double_t w, Double_t h,
                                 UInt_t rubberFlag)
   : TGLayoutHints(kLHintsNormal, 0,0,0,0)
{
   fX    = x;
   fY    = y;
   fW    = w;
   fH    = h;
   fFlag = rubberFlag;
}

////////////////////////////////////////////////////////////////////////////////
/// Save XY layout hints as a C++ statement(s) on output stream.

void TGXYLayoutHints::SavePrimitive(std::ostream &out, Option_t * /*option = ""*/)
{
   TString flag = "";
   if (fFlag & kLRubberX) {
      if (flag.Length() == 0)  flag  = "TGXYLayoutHints::kLRubberX";
      else                     flag += " | TGXYLayoutHints::kLRubberX";
   }
   if (fFlag & kLRubberY) {
      if  (flag.Length() == 0) flag  = "TGXYLayoutHints::kLRubberY";
      else                      flag += " | TGXYLayoutHints::kLRubberY";
   }
   if (fFlag & kLRubberW) {
      if (flag.Length() == 0) flag  = "TGXYLayoutHints::kLRubberW";
      else                     flag += " | TGXYLayoutHints::kLRubberW";
   }
   if (fFlag & kLRubberH) {
      if (flag.Length() == 0) flag  = "TGXYLayoutHints::kLRubberH";
      else                     flag += " | TGXYLayoutHints::kLRubberH";
   }

   out << ", new TGXYLayoutHints(" << GetX() << ", " << GetY() << ", "
       << GetW() << ", " << GetH();

   if (!flag.Length())
      out << ")";
   else
      out << ", " << flag << ")";

}

////////////////////////////////////////////////////////////////////////////////
/// Constructor. The main is the frame for which this layout manager works.

TGXYLayout::TGXYLayout(TGCompositeFrame *main)
{
   UInt_t  width, height;
   Int_t   dummy;

   fMain = main;
   fList = main->GetList();
   fFirst = kTRUE;
   fFirstWidth = fFirstHeight = 0;

   FontStruct_t fs = TGLabel::GetDefaultFontStruct();

   // get standard width an height of a character
   fTWidth = gVirtualX->TextWidth(fs, "1234567890", 10) / 10;
   gVirtualX->GetFontProperties(fs, fTHeight, dummy);

   // the size of the main window are defined in units of a character
   // but the system does not understand this. We have to recalculate
   // the size into pixels.
   width  = main->GetWidth() * fTWidth;
   height = main->GetHeight() * fTHeight;

   main->Resize(width, height);
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGXYLayout::TGXYLayout(const TGXYLayout& xyl) :
  TGLayoutManager(xyl),
  fList(xyl.fList),
  fMain(xyl.fMain),
  fFirst(xyl.fFirst),
  fFirstWidth(xyl.fFirstWidth),
  fFirstHeight(xyl.fFirstHeight),
  fTWidth(xyl.fTWidth),
  fTHeight(xyl.fTHeight)
{
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGXYLayout& TGXYLayout::operator=(const TGXYLayout& xyl)
{
   if(this!=&xyl) {
      TGLayoutManager::operator=(xyl);
      fList=xyl.fList;
      fMain=xyl.fMain;
      fFirst=xyl.fFirst;
      fFirstWidth=xyl.fFirstWidth;
      fFirstHeight=xyl.fFirstHeight;
      fTWidth=xyl.fTWidth;
      fTHeight=xyl.fTHeight;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Recalculates the postion and the size of all widgets.

void TGXYLayout::Layout()
{
   TGFrameElement   *ptr;
   TGXYLayoutHints  *layout;
   Double_t          xFactor;
   Double_t          yFactor;
   Int_t             newX, newY;
   UInt_t            newW, newH;
   Double_t          temp;

   if (!fList) return;

   if (fFirst) {
      // save the original size of the frame. It is used to determin
      // if the user has changed the window
      fFirstWidth   = fMain->GetWidth();
      fFirstHeight  = fMain->GetHeight();
      fFirst        = kFALSE;
   }

   // get the factor of the increacement of the window
   xFactor = (Double_t)fMain->GetWidth() / (Double_t)fFirstWidth;
   if (xFactor < 1.0) xFactor = 1.0;
   yFactor = (Double_t)fMain->GetHeight() / (Double_t)fFirstHeight;
   if (yFactor < 1.0) yFactor = 1.0;

   // set the position an size for each widget and call the layout
   // function for each widget
   TIter next(fList);
   while ((ptr = (TGFrameElement *) next()))  {
      if (ptr->fState & kIsVisible) {
         layout = (TGXYLayoutHints*)ptr->fLayout;
         if (layout == 0)
            continue;

         temp = layout->GetX() * fTWidth ;
         if (layout->GetFlag() & TGXYLayoutHints::kLRubberX)
            temp *= xFactor;
         newX = (Int_t)(temp + 0.5);

         temp = layout->GetY() * fTHeight;
         if (layout->GetFlag() & TGXYLayoutHints::kLRubberY)
            temp *= yFactor;
         newY = (Int_t)(temp + 0.5);

         temp = layout->GetW() * fTWidth;
         if (layout->GetFlag() & TGXYLayoutHints::kLRubberW)
            temp *= xFactor;
         newW = (UInt_t)(temp + 0.5);

         temp = layout->GetH() * fTHeight;
         if (layout->GetFlag() & TGXYLayoutHints::kLRubberH)
            temp *= yFactor;
         newH = (UInt_t)(temp + 0.5);
         ptr->fFrame->MoveResize(newX, newY, newW, newH);
         ptr->fFrame->Layout();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the original size of the frame.

TGDimension TGXYLayout::GetDefaultSize() const
{
   TGDimension size(fFirstWidth, fFirstHeight);

   return size;
}

////////////////////////////////////////////////////////////////////////////////
/// Save XY layout manager as a C++ statement(s) on output stream.

void TGXYLayout::SavePrimitive(std::ostream &out, Option_t * /*option = ""*/)
{
   out << "new TGXYLayout(" << fMain->GetName() << ")";

}
