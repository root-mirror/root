// @(#)root/gui:$Id$
// Author: David Gonzalez Maline  19/07/2006

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGInputDialog
#define ROOT_TGInputDialog

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// Input Dialog Widget                                                   //
//                                                                       //
// An Input dialog box                                                   //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "TGFrame.h"

class TGLabel;
class TGTextEntry;
class TGTextButton;

class TGInputDialog : public TGTransientFrame {

private:
   TGLabel          *fLabel{nullptr};   // text entry label
   TGTextEntry      *fTE{nullptr};      // text entry widget
   TGTextButton     *fOk{nullptr};      // ok button
   TGTextButton     *fCancel{nullptr};  // cancel button
   char             *fRetStr{nullptr};  // address to store return string
   char             *fOwnBuf{nullptr};  // internal buffer when return string not specified

   TGInputDialog(const TGInputDialog&) = delete;
   TGInputDialog &operator= (const TGInputDialog&) = delete;

public:
   TGInputDialog(const TGWindow *p = nullptr, const TGWindow *main = nullptr,
                 const char *prompt = nullptr, const char *defval = nullptr,
                 char *retstr = nullptr, UInt_t options = kVerticalFrame);
   ~TGInputDialog();
   virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t);

   ClassDef(TGInputDialog, 0)  // Simple input dialog

};

#endif
