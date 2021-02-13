// @(#)root/gui:$Id$
// Author: Fons Rademakers   6/09/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGHorizontal3DLine and TGVertical3DLine                              //
//                                                                      //
// A horizontal 3D line is a line that typically separates a toolbar    //
// from the menubar.                                                    //
// A vertical 3D line is a line that can be used to separate groups of  //
// widgets.                                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TG3DLine.h"
#include "TVirtualX.h"

#include <iostream>


ClassImp(TGHorizontal3DLine);
ClassImp(TGVertical3DLine);

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGHorizontal3DLine::TGHorizontal3DLine(const TGWindow *p, UInt_t w, UInt_t h,
                                       UInt_t options, Pixel_t back) :
                    TGFrame(p, w, h, options, back)
{
   SetWindowName();
   fEditDisabled = kEditDisableHeight;
}

////////////////////////////////////////////////////////////////////////////////
/// draw border

void TGHorizontal3DLine::DrawBorder()
{
   gVirtualX->DrawLine(fId, GetShadowGC()(),  0, 0, fWidth-2, 0);
   gVirtualX->DrawLine(fId, GetHilightGC()(), 0, 1, fWidth-1, 1);
   gVirtualX->DrawLine(fId, GetHilightGC()(), fWidth-1, 0, fWidth-1, 1);
}


////////////////////////////////////////////////////////////////////////////////
/// Save an vertical 3D line as a C++ statement(s) on output stream out.

void TGHorizontal3DLine::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   out << "   TGHorizontal3DLine *";
   out << GetName() << " = new TGHorizontal3DLine(" << fParent->GetName()
       << "," << GetWidth() << "," << GetHeight();

   if (fBackground == GetDefaultFrameBackground()) {
      if (!GetOptions()) {
         out << ");" << std::endl;
      } else {
         out << "," << GetOptionString() << ");" << std::endl;
      }
   } else {
      out << "," << GetOptionString() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGVertical3DLine::TGVertical3DLine(const TGWindow *p, UInt_t w, UInt_t h,
                                   UInt_t options, Pixel_t back) :
                  TGFrame(p, w, h, options, back)
{
   SetWindowName();
   fEditDisabled = kEditDisableWidth;
}

////////////////////////////////////////////////////////////////////////////////
/// draw border

void TGVertical3DLine::DrawBorder()
{
   gVirtualX->DrawLine(fId, GetShadowGC()(),  0, 0, 0, fHeight-2);
   gVirtualX->DrawLine(fId, GetHilightGC()(), 1, 0, 1, fHeight-1);
   gVirtualX->DrawLine(fId, GetHilightGC()(), 0, fHeight-1, 1, fHeight-1);
}

////////////////////////////////////////////////////////////////////////////////
/// Save an vertical 3D line as a C++ statement(s) on output stream out.

void TGVertical3DLine::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   out << "   TGVertical3DLine *";
   out << GetName() << " = new TGVertical3DLine(" << fParent->GetName()
       << "," << GetWidth() << "," << GetHeight();

   if (fBackground == GetDefaultFrameBackground()) {
      if (!GetOptions()) {
         out << ");" << std::endl;
      } else {
         out << "," << GetOptionString() <<");" << std::endl;
      }
   } else {
      out << "," << GetOptionString() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;
}
