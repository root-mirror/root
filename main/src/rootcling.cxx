// Authors: Axel Naumann, Philippe Canal, Danilo Piparo

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "rootclingTCling.h"
#include "rootclingIO.h"
#include "rootcling_impl.h"
#include "RConfigure.h"
#include <ROOT/RConfig.hxx>
#include "TSystem.h"

extern "C" {
   R__DLLEXPORT void usedToIdentifyRootClingByDlSym() {}
}

// force compiler to emit symbol for function above
static void (*dlsymaddr)() = &usedToIdentifyRootClingByDlSym;

int main(int argc, char **argv)
{
   (void) dlsymaddr; // avoid unused variable warning

   ROOT::Internal::RootCling::DriverConfig config{};

   config.fBuildingROOTStage1 = false;
   config.fPRootDir = &gRootDir;
   config.fTROOT__GetExtraInterpreterArgs = &TROOT__GetExtraInterpreterArgs;
   config.fTROOT__GetIncludeDir = &TROOT__GetIncludeDir;
   config.fTROOT__GetEtcDir = &TROOT__GetEtcDir;
   config.fTCling__GetInterpreter = &TCling__GetInterpreter;
   config.fInitializeStreamerInfoROOTFile = &InitializeStreamerInfoROOTFile;
   config.fAddStreamerInfoToROOTFile = &AddStreamerInfoToROOTFile;
   config.fAddTypedefToROOTFile = &AddTypedefToROOTFile;
   config.fAddEnumToROOTFile = &AddEnumToROOTFile;
   config.fCloseStreamerInfoROOTFile = &CloseStreamerInfoROOTFile;

   return ROOT_rootcling_Driver(argc, argv, config);
}
