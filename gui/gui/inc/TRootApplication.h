// @(#)root/gui:$Id$
// Author: Fons Rademakers   15/01/98

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#ifndef ROOT_TRootApplication
#define ROOT_TRootApplication

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRootApplication                                                     //
//                                                                      //
// This class create the ROOT native GUI version of the ROOT            //
// application environment. This in contrast the Win32 version.         //
// Once the native widgets work on Win32 this class can be folded into  //
// the TApplication class (since all graphic will go via TVirtualX).    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TApplicationImp.h"


class TGClient;


class TRootApplication : public TApplicationImp {

private:
   TGClient    *fClient;        // pointer to the client environment
   char        *fDisplay;       // display server to connect to

   TRootApplication() { fClient = 0; fDisplay = 0; }
   void GetOptions(Int_t *argc, char **argv);

public:
   TRootApplication(const char *appClassName, Int_t *argc, char **argv);
   ~TRootApplication() override;

   TGClient     *Client() const { return fClient; }

   void    Show() override { }
   void    Hide() override { }
   void    Iconify() override { }
   Bool_t  IsCmdThread() override;
   void    Init() override { }
   void    Open() override { }
   void    Raise() override { }
   void    Lower() override { }

   ClassDef(TRootApplication,0)  // ROOT native GUI application environment
};

#endif
