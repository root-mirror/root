// @(#)root/hist:$Id: TGraph2DErrors.h,v 1.00
// Author: Olivier Couet 26/11/2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGraph2DErrors
#define ROOT_TGraph2DErrors


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGraph2DErrors                                                       //
//                                                                      //
// a 2D Graph with error bars                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGraph2D.h"

class TGraph2DErrors : public TGraph2D {

private:


protected:
   Double_t *fEX; ///<[fNpoints] array of X errors
   Double_t *fEY; ///<[fNpoints] array of Y errors
   Double_t *fEZ; ///<[fNpoints] array of Z errors

public:
   TGraph2DErrors();
   TGraph2DErrors(Int_t n);
   TGraph2DErrors(Int_t n, Double_t *x, Double_t *y, Double_t *z,
                  Double_t *ex=0, Double_t *ey=0, Double_t *ez=0, Option_t *option="");
   TGraph2DErrors(const TGraph2DErrors&);
   TGraph2DErrors& operator=(const TGraph2DErrors&);
   virtual ~TGraph2DErrors();
   Double_t        GetErrorX(Int_t bin) const;
   Double_t        GetErrorY(Int_t bin) const;
   Double_t        GetErrorZ(Int_t bin) const;
   Double_t       *GetEX() const {return fEX;}
   Double_t       *GetEY() const {return fEY;}
   Double_t       *GetEZ() const {return fEZ;}
   Double_t        GetXmaxE() const;
   Double_t        GetXminE() const;
   Double_t        GetYmaxE() const;
   Double_t        GetYminE() const;
   Double_t        GetZmaxE() const;
   Double_t        GetZminE() const;
   virtual void    Print(Option_t *chopt="") const;
   Int_t           RemovePoint(Int_t ipoint); // *MENU*
   virtual void    Scale(Double_t c1=1., Option_t *option="z"); // *MENU*
   virtual void    Set(Int_t n);
   virtual void    SetPoint(Int_t i, Double_t x, Double_t y, Double_t z);
   virtual void    SetPointError(Int_t i, Double_t ex, Double_t ey, Double_t ez);

   ClassDef(TGraph2DErrors,1)  //A 2D graph with error bars
};

#endif


