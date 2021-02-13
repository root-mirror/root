// @(#)root/graf:$Id$
// Author: Sebastian Boser, 02/02/06

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGraphPolar
#define ROOT_TGraphPolar

#include "TGraphErrors.h"

#ifdef R__LESS_INCLUDES
class TGraphPolargram;
#else
#include "TGraphPolargram.h"
#endif

class TGraphPolar: public TGraphErrors {

private:
   Bool_t fOptionAxis;          ///< Force drawing of new coord system

protected:
   TGraphPolargram* fPolargram; ///< The polar coordinates system
   Double_t* fXpol;             ///< [fNpoints] points in polar coordinates
   Double_t* fYpol;             ///< [fNpoints] points in polar coordinates


public:
   TGraphPolar();
   TGraphPolar(Int_t n, const Double_t* theta=0, const Double_t* r=0,
                        const Double_t* etheta=0, const Double_t* er=0);
   virtual ~TGraphPolar();

   TGraphPolargram *GetPolargram() {return fPolargram;};

   void             Draw(Option_t* options = "");
   Bool_t           GetOptionAxis() {return fOptionAxis;};
   void             SetMaxRadial(Double_t maximum = 1); //*MENU*
   void             SetMinRadial(Double_t minimum = 0); //*MENU*
   void             SetMaximum(Double_t maximum = 1) {SetMaxRadial(maximum);}
   void             SetMinimum(Double_t minimum = 0) {SetMinRadial(minimum);}
   void             SetMaxPolar(Double_t maximum = 6.28318530717958623); //*MENU*
   void             SetMinPolar(Double_t minimum = 0); //*MENU*
   void             SetOptionAxis(Bool_t opt) {fOptionAxis = opt;};
   void             SetPolargram(TGraphPolargram *p) {fPolargram = p;};
   Double_t        *GetXpol();
   Double_t        *GetYpol();

   ClassDef(TGraphPolar,1); // Polar graph
};

#endif
