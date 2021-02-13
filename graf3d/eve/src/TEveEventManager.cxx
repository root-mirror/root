// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveEventManager.h"

#include "TInterpreter.h"

/** \class TEveEventManager
\ingroup TEve
Base class for event management and navigation.
*/

ClassImp(TEveEventManager);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveEventManager::TEveEventManager(const char* n, const char* t) :
   TEveElementList(n, t),
   fNewEventCommands()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual function to be called after a new event is loaded.
/// It iterates over the list of registered commands
/// (fNewEventCommands) and executes them in given order.

void TEveEventManager::AfterNewEventLoaded()
{
   for (std::vector<TString>::iterator i = fNewEventCommands.begin(); i != fNewEventCommands.end(); ++i)
   {
      gInterpreter->ProcessLine(*i);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Register a command to be executed on each new event.

void TEveEventManager::AddNewEventCommand(const TString& cmd)
{
   fNewEventCommands.push_back(cmd);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the first command equal to cmd.

void TEveEventManager::RemoveNewEventCommand(const TString& cmd)
{
   for (std::vector<TString>::iterator i = fNewEventCommands.begin(); i != fNewEventCommands.end(); ++i)
   {
      if (cmd == *i) {
         fNewEventCommands.erase(i);
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the list of commands to be executed on each new event.

void TEveEventManager::ClearNewEventCommands()
{
   fNewEventCommands.clear();
}
