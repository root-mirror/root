/// \file
/// \ingroup tutorial_v7
///
/// \macro_code
///
/// \date 2015-03-22
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!
/// \author Axel Naumann <axel@cern.ch>

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RCanvas.hxx"
#include "ROOT/RPad.hxx"
#include "ROOT/RLine.hxx"

void pad()
{
   using namespace ROOT::Experimental;

   auto canvas = RCanvas::Create("what to do with a pad!");
   auto pads   = canvas->Divide(3, 3);

   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
         pads[i][j]->Draw<RLine>()->SetP1({0.1_normal, 0.1_normal}).SetP2({0.9_normal, 0.9_normal});
         pads[i][j]->Draw<RLine>()->SetP1({0.1_normal, 0.9_normal}).SetP2({0.9_normal, 0.1_normal});
      }

   canvas->Show();
}
