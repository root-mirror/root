// @(#)root/gui:$Id$
// Author: Fons Rademakers   22/02/98

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TRootControlBar
#define ROOT_TRootControlBar


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRootControlBar                                                      //
//                                                                      //
// This class provides an interface to the GUI dependent functions of   //
// the TControlBar class. A control bar is a horizontal or vertical bar //
// with a number of buttons (text or picture buttons).                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TControlBarImp.h"
#include "TGFrame.h"

class TControlBar;
class TList;


class TRootControlBar : public TGMainFrame, public TControlBarImp {

private:
   TList          *fWidgets; // list of TGTextButton or TGPictureButtons
   TGLayoutHints  *fL1;      // button layout hints
   UInt_t          fBwidth;  // button width in pixels

public:
   TRootControlBar(TControlBar *c = 0, const char *title = "ROOT Control Bar",
                   Int_t x = -999, Int_t y = -999);
   ~TRootControlBar() override;

   void Create() override;
   void Hide() override;
   void Show() override;

   TList *GetWidgets() const { return fWidgets; }

   // overridden from TGMainFrame
   void   CloseWindow() override;
   Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2) override;
   void   ReallyDelete() override;
   void   SetButtonState(const char *label, Int_t state = 0) override;
   void   SetButtonWidth(UInt_t width) override;
   void   SetFont(const char *fontName) override;
   void   SetTextColor(const char *colorName) override;

   ClassDef(TRootControlBar,0)  //ROOT native GUI implementation of TControlBar
};

#endif
