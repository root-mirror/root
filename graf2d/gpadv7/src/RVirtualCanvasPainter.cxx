/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RVirtualCanvasPainter.hxx>

#include "ROOT/RAttrBase.hxx" // for GPadLog()
#include <ROOT/RLogger.hxx>
#include <TSystem.h> // TSystem::Load

#include <exception>

namespace {
static int LoadCanvasPainterLibraryOnce() {
  static int loadResult = gSystem->Load("libROOTCanvasPainter");
  if (loadResult != 0)
     R__LOG_ERROR(ROOT::Experimental::GPadLog()) << "Loading of libROOTCanvasPainter failed!";
  return loadResult;
}
static void LoadCanvasPainterLibrary() {
  static int loadResult = LoadCanvasPainterLibraryOnce();
  (void) loadResult;
}
} // unnamed namespace

using namespace ROOT::Experimental::Internal;

/// The implementation is here to pin the vtable.
RVirtualCanvasPainter::~RVirtualCanvasPainter() = default;

std::unique_ptr<RVirtualCanvasPainter::Generator> &RVirtualCanvasPainter::GetGenerator()
{
   /// The generator for implementations.
   static std::unique_ptr<Generator> generator;
   return generator;
}

std::unique_ptr<RVirtualCanvasPainter> RVirtualCanvasPainter::Create(ROOT::Experimental::RCanvas &canv)
{
   if (!GetGenerator()) {
      LoadCanvasPainterLibrary();
      if (!GetGenerator()) {
         R__LOG_ERROR(GPadLog()) << "RVirtualCanvasPainter::Generator failed to register!";
         throw std::runtime_error("RVirtualCanvasPainter::Generator failed to initialize");
      }
   }
   return GetGenerator()->Create(canv);
}

/// The implementation is here to pin the vtable.
RVirtualCanvasPainter::Generator::~Generator() = default;

