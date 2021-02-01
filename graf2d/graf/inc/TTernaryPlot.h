// @(#)root/graf:$Id$
// Author: Olivier Couet   27/10/20

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTernaryPlot
#define ROOT_TTernaryPlot

#include "TObject.h"


class TTernaryPlot: public TObject {

protected:
   Int_t              fNpoints;   // Number of points
   Int_t              fMaxSize;   // Current dimension of arrays fX and fY

   Double_t          *fX;         // [fNpoints] array of X points
   Double_t          *fY;         // [fNpoints] array of Y points
   Double_t yC;

public:

   TTernaryPlot(Int_t n);

   virtual ~TTernaryPlot();

   void SetPoint(Double_t u, Double_t v, Option_t *option);
   void Draw(Option_t *option = "");

   ClassDef(TTernaryPlot,1)  //A ternary plot
};

#endif

