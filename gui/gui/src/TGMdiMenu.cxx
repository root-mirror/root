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
// TGMdiMenu.                                                           //
//                                                                      //
// This file contains the TGMdiMenuBar class.                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGMdi.h"
#include "TGMdiMenu.h"
#include "TList.h"

#include <iostream>


ClassImp(TGMdiMenuBar);

////////////////////////////////////////////////////////////////////////////////
/// TGMdiMenuBar constructor.

TGMdiMenuBar::TGMdiMenuBar(const TGWindow *p, int w, int h) :
   TGCompositeFrame(p, w, h, kHorizontalFrame)
{
   fLHint = new TGLayoutHints(kLHintsNormal);
   fLeftHint = new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 1, 1, 1, 1);
   fBarHint = new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 1, 1, 1, 1);
   fRightHint = new TGLayoutHints(kLHintsRight | kLHintsCenterY, 1, 2, 1, 1);

   fLeft = new TGCompositeFrame(this, 10, 10, kHorizontalFrame);
   fBar = new TGMenuBar(this, 1, 20, kHorizontalFrame);
   fRight = new TGCompositeFrame(this, 10, 10, kHorizontalFrame);

   AddFrame(fLeft,  fLeftHint);
   AddFrame(fBar,   fBarHint);
   AddFrame(fRight, fRightHint);
}

////////////////////////////////////////////////////////////////////////////////
/// TGMdiMenuBar destructor.

TGMdiMenuBar::~TGMdiMenuBar()
{
   if (!MustCleanup()) {
      delete fLHint;
      delete fLeftHint;
      delete fRightHint;
      delete fBarHint;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add popup menu to the MDI menu bar with layout hints l.

void TGMdiMenuBar::AddPopup(TGHotString *s, TGPopupMenu *menu, TGLayoutHints *l)
{
   fBar->AddPopup(s, menu, l);
   // Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// This is called from TGMdiMainFrame on Maximize().

void TGMdiMenuBar::AddFrames(TGMdiTitleIcon *icon, TGMdiButtons *buttons)
{
   // Hide all frames first
   TGFrameElement *el;
   TIter nextl(fLeft->GetList());
   while ((el = (TGFrameElement *) nextl())) {
      fLeft->HideFrame(el->fFrame);
   }
   TIter nextr(fRight->GetList());
   while ((el = (TGFrameElement *) nextr())) {
      fRight->HideFrame(el->fFrame);
   }
   // Then add specified frames
   icon->ReparentWindow(fLeft);
   buttons->ReparentWindow(fRight);
   fLeft->AddFrame(icon, fLHint);
   fRight->AddFrame(buttons, fLHint);
}

////////////////////////////////////////////////////////////////////////////////
/// This is called from TGMdiMainFrame on Restore()

void TGMdiMenuBar::RemoveFrames(TGMdiTitleIcon *icon, TGMdiButtons *buttons)
{
   // Remove specified frames
   fLeft->RemoveFrame(icon);
   fRight->RemoveFrame(buttons);
   // Then show (restore) last frames
   TGFrameElement *el;
   el = (TGFrameElement *)fLeft->GetList()->Last();
   if (el)
      fLeft->ShowFrame(el->fFrame);
   el = (TGFrameElement *)fRight->GetList()->Last();
   if (el)
      fRight->ShowFrame(el->fFrame);
}

////////////////////////////////////////////////////////////////////////////////
/// This is called from TGMdiMainFrame on Maximize().

void TGMdiMenuBar::ShowFrames(TGMdiTitleIcon *icon, TGMdiButtons *buttons)
{
   // Hide all frames first
   TGFrameElement *el;
   TIter nextl(fLeft->GetList());
   while ((el = (TGFrameElement *) nextl())) {
      fLeft->HideFrame(el->fFrame);
   }
   TIter nextr(fRight->GetList());
   while ((el = (TGFrameElement *) nextr())) {
      fRight->HideFrame(el->fFrame);
   }
   // Then show specified frames
   fLeft->ShowFrame(icon);
   fRight->ShowFrame(buttons);
}

////////////////////////////////////////////////////////////////////////////////
/// Used to hide specific frames from menu bar

void TGMdiMenuBar::HideFrames(TGMdiTitleIcon *icon, TGMdiButtons *buttons)
{
   // Hide specified frames
   fLeft->HideFrame(icon);
   fRight->HideFrame(buttons);

   // Then show (restore) last frames
   TGFrameElement *el;
   el = (TGFrameElement *)fLeft->GetList()->Last();
   if (el)
      fLeft->ShowFrame(el->fFrame);
   el = (TGFrameElement *)fRight->GetList()->Last();
   if (el)
      fRight->ShowFrame(el->fFrame);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a MDI menu as a C++ statement(s) on output stream out

void TGMdiMenuBar::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   out << std::endl;
   out << "   // MDI menu bar" << std::endl;

   out << "   TGMdiMenuBar *";
   out << GetName() << " = new TGMdiMenuBar(" << fParent->GetName()
       << "," << GetWidth() << "," << GetHeight() << ");" << std::endl;
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   if (!fList) return;

   out << "   TGMenuBar *" << fBar->GetName() << " = " << GetName()
       << "->GetMenuBar();" << std::endl;

   TGFrameElement *el;
   TIter next(fBar->GetList());

   while ((el = (TGFrameElement *)next())) {
      el->fFrame->SavePrimitive(out, option);
      el->fLayout->SavePrimitive(out, option);
      out << ");" << std::endl;
   }
}
