// @(#)root/graf:$Id$
// Author: Rene Brun   19/11/95

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPavesText
#define ROOT_TPavesText


#include "TPaveText.h"

class TPavesText : public TPaveText {

protected:
   Int_t      fNpaves;        ///< Number of stacked paves

public:
   TPavesText();
   TPavesText(Double_t x1, Double_t y1,Double_t x2 ,Double_t y2, Int_t npaves=5, Option_t *option="br");
   TPavesText(const TPavesText &pavestext);
   virtual ~TPavesText();

   virtual void  Draw(Option_t *option="");
   virtual Int_t GetNpaves() {return fNpaves;}
   virtual void  Paint(Option_t *option="");
   virtual void  SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void  SetNpaves(Int_t npaves=5) {fNpaves=npaves;} // *MENU*

   ClassDef(TPavesText,1)  //Stacked Paves with text strings
};

#endif

