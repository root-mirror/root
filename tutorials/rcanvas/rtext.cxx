/// \file
/// \ingroup tutorial_rcanvas
///
/// This macro demonstrate the text attributes for RText. Angle, size and color are
/// changed in a loop. The text alignment and the text font are fixed.
///
/// \macro_image (rcanvas_js)
/// \macro_code
///
/// \date 2017-10-17
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!
/// \author Olivier Couet <Olivier.Couet@cern.ch>

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RCanvas.hxx"
#include "ROOT/RColor.hxx"
#include "ROOT/RText.hxx"
#include "ROOT/RPadPos.hxx"

void rtext()
{
   using namespace ROOT::Experimental;

   // Create a canvas to be displayed.
   auto canvas = RCanvas::Create("RText example");

   for (int i=0; i<=360; i+=10) {
      auto draw = canvas->Draw<RText>(RPadPos(0.5_normal, 0.6_normal), "____  Hello World");

      draw->text.color = RColor((int) (0.38*i), (int) (0.64*i), (int) (0.76*i));
      draw->text.size = 10+i/10;
      draw->text.angle = i;
      draw->text.align = 13;
      draw->text.font = 4;
   }

   canvas->Show();
}
