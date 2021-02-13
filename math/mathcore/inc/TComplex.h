// @(#)root/mathcore:$Id$
// Author: Federico Carminati   22/04/2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TComplex
#define ROOT_TComplex

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TComplex                                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMath.h"

#include "Rtypes.h"
#include <string>


class TComplex {

protected:
   Double_t fRe;    // real part
   Double_t fIm;    // imaginary part

public:
   // ctors and dtors
   TComplex(): fRe(0), fIm(0) {}
   TComplex(Double_t re, Double_t im=0, Bool_t polar=kFALSE);
   virtual ~TComplex() {}

   // constants
   static TComplex I() {return TComplex(0,1);}
   static TComplex One() {return TComplex(1,0);}

   // getters and setters
   Double_t Re() const {return fRe;}
   Double_t Im() const {return fIm;}
   Double_t Rho() const {return TMath::Sqrt(fRe*fRe+fIm*fIm);}
   Double_t Rho2() const {return fRe*fRe+fIm*fIm;}
   Double_t Theta() const {return (fIm||fRe)?TMath::ATan2(fIm,fRe):0;}
   TComplex operator()(Double_t x, Double_t y, Bool_t polar=kFALSE)
      { if (polar) { fRe = x*TMath::Cos(y); fIm = x*TMath::Sin(y); }
        else { fRe = x; fIm = y; } return *this; }

   // Simple operators complex - complex
   TComplex operator *(const TComplex & c) const
      {return TComplex(fRe*c.fRe-fIm*c.fIm,fRe*c.fIm+fIm*c.fRe);}
   TComplex operator +(const TComplex & c) const
      {return TComplex(fRe+c.fRe, fIm+c.fIm);}
   TComplex operator /(const TComplex & c) const
      {return TComplex(fRe*c.fRe+fIm*c.fIm,-fRe*c.fIm+fIm*c.fRe)/c.Rho2();}
   TComplex operator -(const TComplex & c) const
      {return TComplex(fRe-c.fRe, fIm-c.fIm);}

   TComplex operator *=(const TComplex & c)
      {return ((*this) = (*this) * c);}
   TComplex operator +=(const TComplex & c)
      {return ((*this) = (*this) + c);}
   TComplex operator /=(const TComplex & c)
      {return ((*this) = (*this) / c);}
   TComplex operator -=(const TComplex & c)
      {return ((*this) = (*this) - c);}

   TComplex operator -()
      {return TComplex(-fRe,-fIm);}
   TComplex operator +()
      {return *this;}

   // Simple operators complex - double
   TComplex operator *(Double_t c) const
      {return TComplex(fRe*c,fIm*c);}
   TComplex operator +(Double_t c) const
      {return TComplex(fRe+c, fIm);}
   TComplex operator /(Double_t c) const
      {return TComplex(fRe/c,fIm/c);}
   TComplex operator -(Double_t c) const
      {return TComplex(fRe-c, fIm);}

   // Simple operators double - complex
   friend TComplex operator *(Double_t d, const TComplex & c)
      {return TComplex(d*c.fRe,d*c.fIm);}
   friend TComplex operator +(Double_t d, const TComplex & c)
      {return TComplex(d+c.fRe, c.fIm);}
   friend TComplex operator /(Double_t d, const TComplex & c)
      {return TComplex(d*c.fRe,-d*c.fIm)/c.Rho2();}
   friend TComplex operator -(Double_t d, const TComplex & c)
      {return TComplex(d-c.fRe, -c.fIm);}

   // Convertors
   operator Double_t () const {return fRe;}
   operator Float_t  () const {return static_cast<Float_t>(fRe);}
   operator Int_t    () const {return static_cast<Int_t>(fRe);}

   // TMath:: extensions
   static TComplex Sqrt(const TComplex &c)
      {return TComplex(TMath::Sqrt(c.Rho()),0.5*c.Theta(),kTRUE);}

   static TComplex Exp(const TComplex &c)
      {return TComplex(TMath::Exp(c.fRe),c.fIm,kTRUE);}
   static TComplex Log(const TComplex &c)
      {return TComplex(0.5*TMath::Log(c.Rho2()),c.Theta());}
   static TComplex Log2(const TComplex &c)
      {return Log(c)/TMath::Log(2);}
   static TComplex Log10(const TComplex &c)
      {return Log(c)/TMath::Log(10);}

   static TComplex Sin(const TComplex &c)
      {return TComplex(TMath::Sin(c.fRe)*TMath::CosH(c.fIm),
                       TMath::Cos(c.fRe)*TMath::SinH(c.fIm));}
   static TComplex Cos(const TComplex &c)
      {return TComplex(TMath::Cos(c.fRe)*TMath::CosH(c.fIm),
                       -TMath::Sin(c.fRe)*TMath::SinH(c.fIm));}
   static TComplex Tan(const TComplex &c)
      {TComplex cc=Cos(c); return Sin(c)*Conjugate(cc)/cc.Rho2();}

   static TComplex ASin(const TComplex &c)
      {return -I()*Log(I()*c+TMath::Sign(1.,c.Im())*Sqrt(1.-c*c));}
   static TComplex ACos(const TComplex &c)
      {return -I()*Log(c+TMath::Sign(1.,c.Im())*Sqrt(c*c-1.));}
   static TComplex ATan(const TComplex &c)
      {return -0.5*I()*Log((1.+I()*c)/(1.-I()*c));}

   static TComplex SinH(const TComplex &c)
      {return TComplex(TMath::SinH(c.fRe)*TMath::Cos(c.fIm),
                       TMath::CosH(c.fRe)*TMath::Sin(c.fIm));}
   static TComplex CosH(const TComplex &c)
      {return TComplex(TMath::CosH(c.fRe)*TMath::Cos(c.fIm),
                       TMath::SinH(c.fRe)*TMath::Sin(c.fIm));}
   static TComplex TanH(const TComplex &c)
      {TComplex cc=CosH(c); return SinH(c)*Conjugate(cc)/cc.Rho2();}

   static TComplex ASinH(const TComplex &c)
      {return Log(c+TMath::Sign(1.,c.Im())*Sqrt(c*c+1.));}
   static TComplex ACosH(const TComplex &c)
      {return Log(c+TMath::Sign(1.,c.Im())*Sqrt(c*c-1.));}
   static TComplex ATanH(const TComplex &c)
      {return 0.5*Log((1.+c)/(1.-c));}

   static Double_t Abs(const TComplex &c)
      {return c.Rho();}

   static TComplex Power(const TComplex& x, const TComplex& y)
      {Double_t lrho=TMath::Log(x.Rho());
       Double_t theta=x.Theta();
       return TComplex(TMath::Exp(lrho*y.Re()-theta*y.Im()),
                       lrho*y.Im()+theta*y.Re(),kTRUE);}
   static TComplex Power(const TComplex& x, Double_t y)
      {return TComplex(TMath::Power(x.Rho(),y),x.Theta()*y,kTRUE);}
   static TComplex Power(Double_t x, const TComplex& y)
      {Double_t lrho=TMath::Log(TMath::Abs(x));
       Double_t theta=(x>0)?0:TMath::Pi();
       return TComplex(TMath::Exp(lrho*y.Re()-theta*y.Im()),
                       lrho*y.Im()+theta*y.Re(),kTRUE);}
   static TComplex Power(const TComplex& x, Int_t y)
      {return TComplex(TMath::Power(x.Rho(),y),x.Theta()*y,kTRUE);}

   static Int_t Finite(const TComplex& c)
      {return TMath::Min(TMath::Finite(c.Re()),TMath::Finite(c.Im()));}
   static Int_t IsNaN(const TComplex& c)
      {return TMath::IsNaN(c.Re()) || TMath::IsNaN(c.Im());}

   static TComplex Min(const TComplex &a, const TComplex &b)
      {return a.Rho()<=b.Rho()?a:b;}
   static TComplex Max(const TComplex &a, const TComplex &b)
      {return a.Rho()>=b.Rho()?a:b;}
   static TComplex Normalize(const TComplex &c)
      {return TComplex(1.,c.Theta(),kTRUE);}
   static TComplex Conjugate(const TComplex &c)
      {return TComplex(c.Re(),-c.Im());}
   static TComplex Range(const TComplex &lb, const TComplex &ub, const TComplex &c)
     {return Max(lb,Min(c,ub));}

   // I/O
   friend std::ostream& operator<<(std::ostream& out, const TComplex& c);
   friend std::istream& operator>>(std::istream& in, TComplex& c);

   ClassDef(TComplex,1)  //Complex Class
};

namespace cling {
std::string printValue(TComplex *c);
}

#endif
