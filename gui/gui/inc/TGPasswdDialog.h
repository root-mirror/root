// @(#)root/gui:$Id$
// Author: G. Ganis  10/10/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGPasswdDialog
#define ROOT_TGPasswdDialog


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGPasswdDialog                                                       //
//                                                                      //
// Graphic dialog to enter passwords                                    //
//                                                                      //
// Usage:                                                               //
//                                                                      //
// {                                                                    //
//   // Buffer for the passwd                                           //
//   char pwdbuf[128]                                                   //
//                                                                      //
//   Open the dialog box                                                //
//   TGPasswdDialog dialog("My prompt", pwdbuf, 128);                   //
//                                                                      //
//   // Wait until the user is done                                     //
//   while (gROOT->IsInterrupted())                                     //
//      gSystem->DispatchOneEvent(kFALSE);                              //
//                                                                      //
//   // Password is now in pwdbuf                                       //
//   ...                                                                //
//                                                                      //
// }                                                                    //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"

class TGTransientFrame;
class TGTextButton;
class TGTextEntry;
class TGTextBuffer;


class TGPasswdDialog {

private:
   char             *fPwdBuf;     // buffer where to store the passwd
   Int_t             fPwdLenMax;  // passwd buffer length
   TGTransientFrame *fDialog;     // main frame of this widget
   TGTextButton     *fOk;         // Ok button
   TGTextEntry      *fPasswd;     // Password TextEntry
   TGTextBuffer     *fPasswdText; // Passwd Buffer

public:
   TGPasswdDialog(const char *prompt, char *pwdbuf, Int_t pwdlenmax,
                  UInt_t w = 400, UInt_t h = 400);
   virtual ~TGPasswdDialog();

   void   ReturnPressed();

   // slots
   void   CloseWindow();
   void   DoClose();

   ClassDef(TGPasswdDialog,0)  // Dialog for entering passwords
};

#endif
