// @(#)root/utils:$Id$
// Author: Axel Naumann, 2014-04-07

/*************************************************************************
 * Copyright (C) 1995-2014, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// Provides bindings to TCling (compiled with rtti) from rootcling (compiled
// without rtti).


#ifndef ROOT_ROOTCLINGIO_H_H
#define ROOT_ROOTCLINGIO_H_H

extern "C" {
   void InitializeStreamerInfoROOTFile(const char *filename);
   void AddStreamerInfoToROOTFile(const char *normName);
   void AddTypedefToROOTFile(const char *tdname);
   void AddEnumToROOTFile(const char *tdname);
   bool CloseStreamerInfoROOTFile(bool writeEmptyRootPCM);
}

#endif //ROOT_ROOTCLINGIO_H_H
