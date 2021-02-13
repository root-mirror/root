// @(#)root/gpad:$Id: TCreatePrimitives.h,v 1.0

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TCreatePrimitives
#define ROOT_TCreatePrimitives


#include "TLatex.h"
#include "TCurlyArc.h"
#include "TArrow.h"
#include "TArc.h"
#include "TPavesText.h"
#include "TPaveLabel.h"
#include "TDiamond.h"
#include "TGraph.h"

class TCreatePrimitives {

private:

   static TLine *fgLine;
   static TLatex *fgText;
   static TCurlyLine *fgCLine;
   static TArrow *fgArrow;
   static TCurlyArc *fgCArc;
   static TArc *fgArc;
   static TEllipse *fgEllipse;
   static TPave *fgPave;
   static TPaveText *fgPaveText;
   static TPavesText *fgPavesText;
   static TDiamond *fgDiamond;
   static TPaveLabel *fgPaveLabel;
   static TGraph *fgPolyLine;
   static TBox *fgPadBBox;

public:

   TCreatePrimitives();
   virtual ~TCreatePrimitives();
   static void Ellipse(Int_t event, Int_t px, Int_t py,Int_t mode);
   static void Line(Int_t event, Int_t px, Int_t py, Int_t mode);
   static void Pad(Int_t event, Int_t px, Int_t py, Int_t);
   static void Pave(Int_t event, Int_t px, Int_t py, Int_t mode);
   static void PolyLine(Int_t event, Int_t px, Int_t py, Int_t mode);
   static void Text(Int_t event, Int_t px, Int_t py, Int_t mode);
};

#endif
