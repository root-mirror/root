/// \file
/// \ingroup tutorial_v7
///
/// \macro_code
///
/// \date 2019-04-11
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
/// \author Sergey Linev <S.Linev@gsi.de>
/// \author Iliana Betsou <Iliana.Betsou@cern.ch>

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RDirectory.hxx>
#include <ROOT/RFitPanel.hxx>
#include "TH1.h"
#include "TFile.h"

using namespace ROOT::Experimental;

void fitpanel6()
{
   TFile::Open("hsimple.root");
   if (gFile) {
      gFile->Get("hpx");
      gFile->Get("hpxpy");
      gFile->Get("hprof");
   }

   // create panel
   auto panel = std::make_shared<RFitPanel>("FitPanel");
   RDirectory::Heap().Add("fitpanel", panel);

   auto test = RDirectory::Heap().Create<TH1F>("test", "test","This is test histogram",100,-4,4);
   test->FillRandom("gaus", 10000);

   panel->AssignHistogram(test);

   panel->Show();
}

