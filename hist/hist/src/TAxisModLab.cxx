// @(#)root/graf:$Id$
// Author: Olivier Couet

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <cstdlib>

#include "TAxisModLab.h"

ClassImp(TAxisModLab);

/** \class TAxisModLab
\ingroup BasicGraphics

TAxis helper class used to store the modified labels.
*/

////////////////////////////////////////////////////////////////////////////////
/// TAxisModLab default constructor.

TAxisModLab::TAxisModLab() {
   fLabNum    = 0;
   fTextAngle = -1.;
   fTextSize  = -1.;
   fTextAlign = -1;
   fTextColor = -1;
   fTextFont  = -1;
   fLabText   = "";
}

////////////////////////////////////////////////////////////////////////////////
/// Set modified label number.

void TAxisModLab::SetLabNum(Int_t l) {
   if (l!=0) fLabNum = l;
}

////////////////////////////////////////////////////////////////////////////////
/// Set modified label angle.

void TAxisModLab::SetAngle(Double_t a) {
   if (a>=0.) fTextAngle = a;
}

////////////////////////////////////////////////////////////////////////////////
/// Set modified label size.

void TAxisModLab::SetSize(Double_t s) {
   if (s>=0.) fTextSize  = s;
}

////////////////////////////////////////////////////////////////////////////////
/// Set modified label alignment.

void TAxisModLab::SetAlign(Int_t a) {
   if (a>0) fTextAlign = a;
}

////////////////////////////////////////////////////////////////////////////////
/// Set modified label color.

void TAxisModLab::SetColor(Int_t c) {
   if (c>0) fTextColor = c;
}

////////////////////////////////////////////////////////////////////////////////
/// Set modified label font.

void TAxisModLab::SetFont(Int_t f) {
   if (f>0) fTextFont  = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set modified label text.

void TAxisModLab::SetText(TString s) {
   fLabText   = s;
}
