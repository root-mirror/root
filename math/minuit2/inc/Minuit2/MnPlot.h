// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnPlot
#define ROOT_Minuit2_MnPlot

#include "Minuit2/MnConfig.h"
#include <vector>
#include <utility>

namespace ROOT {

namespace Minuit2 {

/** MnPlot produces a text-screen graphical output of (x,y) points, e.g.
    from Scan or Contours.
*/

class MnPlot {

public:
   MnPlot() : fPageWidth(80), fPageLength(30) {}

   MnPlot(unsigned int width, unsigned int length) : fPageWidth(width), fPageLength(length)
   {
      if (fPageWidth > 120)
         fPageWidth = 120;
      if (fPageLength > 56)
         fPageLength = 56;
   }

   ~MnPlot() {}

   void operator()(const std::vector<std::pair<double, double>> &) const;
   void operator()(double, double, const std::vector<std::pair<double, double>> &) const;

   unsigned int Width() const { return fPageWidth; }
   unsigned int Length() const { return fPageLength; }

private:
   unsigned int fPageWidth;
   unsigned int fPageLength;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnPlot
