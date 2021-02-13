// Authors: Axel Naumann, Philippe Canal, Danilo Piparo

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

namespace cling {
   class Interpreter;
}

namespace ROOT {
namespace Internal {
namespace RootCling {
   struct DriverConfig {
      bool fBuildingROOTStage1 = false;
      const char** fPRootDir = nullptr;

      // Function that might (rootcling) or might not (rootcling_stage1) be there.
      const char ** * (*fTROOT__GetExtraInterpreterArgs)() = nullptr;
      const char *(*fTROOT__GetIncludeDir)() = nullptr;
      const char *(*fTROOT__GetEtcDir)() = nullptr;
      cling::Interpreter *(*fTCling__GetInterpreter)() = nullptr;
      void (*fInitializeStreamerInfoROOTFile)(const char *filename) = nullptr;
      void (*fAddStreamerInfoToROOTFile)(const char *normName) = nullptr;
      void (*fAddTypedefToROOTFile)(const char *tdname) = nullptr;
      void (*fAddEnumToROOTFile)(const char *tdname) = nullptr;
      bool (*fCloseStreamerInfoROOTFile)(bool writeEmptyRootPCM) = nullptr;
   };

   struct TROOTSYSSetter {
     TROOTSYSSetter();
   };
} // namespace RootCling
} // namespace Internal
} // namespace ROOT

#ifndef R__DLLEXPORT
#ifdef _MSC_VER
#define R__DLLEXPORT __declspec(dllexport)
#else
#define R__DLLEXPORT __attribute__ ((visibility ("default")))
#endif
#endif

extern "C" R__DLLEXPORT
int ROOT_rootcling_Driver(int argc, char **argv, const ROOT::Internal::RootCling::DriverConfig& config);
