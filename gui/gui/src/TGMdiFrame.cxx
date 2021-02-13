// @(#)root/gui:$Id$
// Author: Bertrand Bellenot   20/08/2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**************************************************************************

    This file is part of TGMdi an extension to the xclass toolkit.
    Copyright (C) 1998-2002 by Harald Radke, Hector Peraza.

    This application is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This application is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library; if not, write to the Free
    Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

**************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGMdiFrame.                                                          //
//                                                                      //
// This file contains the TGMdiFrame class.                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGMdiFrame.h"
#include "TGMdiMainFrame.h"
#include "TGMdiDecorFrame.h"

#include <iostream>

ClassImp(TGMdiFrame);

////////////////////////////////////////////////////////////////////////////////
/// TGMdiFrame constructor.

TGMdiFrame::TGMdiFrame(TGMdiMainFrame *main, Int_t w, Int_t h, UInt_t options,
                       Pixel_t back) :
   TGCompositeFrame(main->GetContainer(), w, h,
                    options | kOwnBackground | kMdiFrame, back)
{
   fMain = main;
   fMain->AddMdiFrame(this);  // this reparents the window
   fMdiHints = kMdiDefaultHints;
}

////////////////////////////////////////////////////////////////////////////////
/// TGMdiFrame destructor.

TGMdiFrame::~TGMdiFrame()
{
   Cleanup();
   fMain->RemoveMdiFrame(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Close MDI frame window.

Bool_t TGMdiFrame::CloseWindow()
{
   DeleteWindow();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Typically call this method in the slot connected to the CloseWindow()
/// signal to prevent the calling of the default or any derived CloseWindow()
/// methods to prevent premature or double deletion of this window.

void TGMdiFrame::DontCallClose()
{
   SetBit(kDontCallClose);
}

////////////////////////////////////////////////////////////////////////////////
/// Set MDI hints, also used to identify titlebar buttons.

void TGMdiFrame::SetMdiHints(ULong_t mdihints)
{
   fMdiHints = mdihints;
   ((TGMdiDecorFrame *)fParent)->SetMdiButtons(mdihints);
}

////////////////////////////////////////////////////////////////////////////////
/// Set MDI window name (set titlebar title).

void TGMdiFrame::SetWindowName(const char *name)
{
   ((TGMdiDecorFrame *)fParent)->SetWindowName(name);
   fMain->UpdateWinListMenu();
}

////////////////////////////////////////////////////////////////////////////////
/// Set MDI window icon (titlebar icon).

void TGMdiFrame::SetWindowIcon(const TGPicture *pic)
{
   ((TGMdiDecorFrame *)fParent)->SetWindowIcon(pic);
   fMain->UpdateWinListMenu();
}

////////////////////////////////////////////////////////////////////////////////
/// Return MDI window name.

const char *TGMdiFrame::GetWindowName()
{
   return ((TGMdiDecorFrame *)fParent)->GetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to picture used as MDI window icon (on titlebar).

const TGPicture *TGMdiFrame::GetWindowIcon()
{
   return ((TGMdiDecorFrame *)fParent)->GetWindowIcon();
}

////////////////////////////////////////////////////////////////////////////////
/// Move MDI window at position x, y.

void TGMdiFrame::Move(Int_t x, Int_t y)
{
   ((TGMdiDecorFrame *)fParent)->Move(x, y);
   fX = x; fY = y;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a MDI option string - used in SavePrimitive().

TString TGMdiFrame::GetMdiHintsString() const
{
   TString hints;
   if (fMdiHints == kMdiDefaultHints)
      hints = "kMdiDefaultHints";
   else {
      if (fMdiHints & kMdiClose) {
         if (hints.Length() == 0) hints = "kMdiClose";
         else                     hints += " | kMdiClose";
      }
      if (fMdiHints & kMdiRestore) {
         if (hints.Length() == 0) hints = "kMdiRestore";
         else                     hints += " | kMdiRestore";
      }
      if (fMdiHints & kMdiMove) {
         if (hints.Length() == 0) hints = "kMdiMove";
         else                     hints += " | kMdiMove";
      }
      if (fMdiHints & kMdiSize) {
         if (hints.Length() == 0) hints = "kMdiSize";
         else                     hints += " | kMdiSize";
      }
      if (fMdiHints & kMdiMinimize) {
         if (hints.Length() == 0) hints = "kMdiMinimize";
         else                     hints += " | kMdiMinimize";
      }
      if (fMdiHints & kMdiMaximize) {
         if (hints.Length() == 0) hints = "kMdiMaximize";
         else                     hints += " | kMdiMaximize";
      }
      if (fMdiHints & kMdiHelp) {
         if (hints.Length() == 0) hints = "kMdiHelp";
         else                     hints += " | kMdiHelp";
      }
      if (fMdiHints & kMdiMenu) {
         if (hints.Length() == 0) hints = "kMdiMenu";
         else                     hints += " | kMdiMenu";
      }
   }
   return hints;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a MDIframe as a C++ statement(s) on output stream out

void TGMdiFrame::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';

   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   TGMdiTitleBar *tb = fMain->GetWindowList()->GetDecorFrame()->GetTitleBar();

   out << std::endl <<"   // MDI frame "<< quote << GetWindowName() << quote << std::endl;
   out << "   TGMdiFrame *";
   out << GetName() << " = new TGMdiFrame(" << fMain->GetName()
       << "," << GetWidth() + GetBorderWidth()*2
       << "," << GetHeight() + tb->GetHeight() + GetBorderWidth()*2;

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

   SavePrimitiveSubframes(out, option);

   out << "   " << GetName() << "->SetWindowName(" << quote << GetWindowName()
       << quote << ");" << std::endl;
   out << "   " << GetName() << "->SetMdiHints(" << GetMdiHintsString()
       << ");" << std::endl;
   if ((GetX() != 5) && (GetY() != 23))
      out << "   " << GetName() << "->Move(" << GetX() << "," << GetY()
          << ");" << std::endl;

   out << "   " << GetName() << "->MapSubwindows();" << std::endl;
   out << "   " << GetName() << "->Layout();" << std::endl;
}
