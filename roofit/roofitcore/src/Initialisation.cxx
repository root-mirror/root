#include "RConfigure.h"
#include "TEnv.h"

#include <iostream>

/**
\file Initialisation.cxx
\ingroup Roofitcore

Run static initialisers on first load of RooFitCore.
**/

namespace {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Print RooFit banner.
void doBanner() {
#ifndef __ROOFIT_NOBANNER
  if (gEnv->GetValue("RooFit.Banner", 1) == 0)
    return;

  /// RooFit version tag.
  constexpr char VTAG[] = "3.60";

  std::cout << '\n'
      << "\033[1mRooFit v" << VTAG << " -- Developed by Wouter Verkerke and David Kirkby\033[0m " << '\n'
      << "                Copyright (C) 2000-2013 NIKHEF, University of California & Stanford University" << '\n'
      << "                All rights reserved, please read http://roofit.sourceforge.net/license.txt" << '\n'
      << std::endl;
#endif
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// A RAII that performs RooFit's static initialisation.
static struct RooFitInitialiser {
  RooFitInitialiser() {
    doBanner();
  }
} __rooFitInitialiser;

}
