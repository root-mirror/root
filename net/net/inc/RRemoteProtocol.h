// @(#)root/net:$Id$
// Author: G. Ganis  10/5/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RRemoteProtocol
#define ROOT_RRemoteProtocol

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RRemoteProtocol                                                      //
//                                                                      //
// Protocol and parameters for remote running                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RtypesCore.h"

// Protocol version we run
// 1              Initial version
const Int_t       kRRemote_Protocol = 1;

// Message types
enum ERootRemMsgTypes {
   kRRT_Undef           = -1,
   kRRT_Fatal           = 0,
   kRRT_Reset           = 1,
   kRRT_CheckFile       = 2,
   kRRT_File            = 3,
   kRRT_LogFile         = 4,
   kRRT_LogDone         = 5,
   kRRT_Protocol        = 6,
   kRRT_GetObject       = 7,
   kRRT_Message         = 8,
   kRRT_Terminate       = 9,
   kRRT_SendFile        = 10
};

// Interrupts
enum ERootRemInterrupt {
   kRRI_Hard          = 1,
   kRRI_Soft          = 2,
   kRRI_Shutdown      = 3
};

#endif
