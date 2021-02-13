// @(#)root/hist:$Id$
// Author: Federico Carminati   28/02/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TSpline
#define ROOT_TSpline

#ifdef R__LESS_INCLUDES
#include "TNamed.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TAttMarker.h"
#else
#include "TGraph.h"
#endif

class TH1;
class TH1F;
class TF1;
class TGraph;

class TSpline : public TNamed, public TAttLine,
                public TAttFill, public TAttMarker
{
protected:
   Double_t  fDelta;     // Distance between equidistant knots
   Double_t  fXmin;      // Minimum value of abscissa
   Double_t  fXmax;      // Maximum value of abscissa
   Int_t     fNp;        // Number of knots
   Bool_t    fKstep;     // True of equidistant knots
   TH1F     *fHistogram; // Temporary histogram
   TGraph   *fGraph;     // Graph for drawing the knots
   Int_t     fNpx;       // Number of points used for graphical representation

   TSpline(const TSpline&);
   TSpline& operator=(const TSpline&);
   virtual void     BuildCoeff()=0;

public:
   TSpline() : fDelta(-1), fXmin(0), fXmax(0),
      fNp(0), fKstep(kFALSE), fHistogram(0), fGraph(0), fNpx(100) {}
   TSpline(const char *title, Double_t delta, Double_t xmin,
      Double_t xmax, Int_t np, Bool_t step) :
      TNamed("Spline",title), TAttFill(0,1),
      fDelta(delta), fXmin(xmin),
      fXmax(xmax), fNp(np), fKstep(step),
      fHistogram(0), fGraph(0), fNpx(100) {}
   virtual ~TSpline();

   virtual void     GetKnot(Int_t i, Double_t &x, Double_t &y) const =0;
   virtual Int_t    DistancetoPrimitive(Int_t px, Int_t py);
   virtual void     Draw(Option_t *option="");
   virtual void     ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual Double_t GetDelta() const {return fDelta;}
   TH1F            *GetHistogram() const {return fHistogram;}
   virtual Int_t    GetNp()    const {return fNp;}
   virtual Int_t    GetNpx()   const {return fNpx;}
   virtual Double_t GetXmin()  const {return fXmin;}
   virtual Double_t GetXmax()  const {return fXmax;}
   virtual void     Paint(Option_t *option="");
   virtual Double_t Eval(Double_t x) const=0;
   virtual void     SaveAs(const char * /*filename*/,Option_t * /*option*/) const {;}
   void             SetNpx(Int_t n) {fNpx=n;}

   ClassDef (TSpline,2) // Spline base class
};


//______________________________________________________________________________
class TSplinePoly : public TObject
{
protected:
   Double_t fX;     // abscissa
   Double_t fY;     // constant term

public:
   TSplinePoly() :
      fX(0), fY(0) {}
   TSplinePoly(Double_t x, Double_t y) :
      fX(x), fY(y) {}
   TSplinePoly(TSplinePoly const &other);
   TSplinePoly &operator=(TSplinePoly const &other);

   Double_t &X() {return fX;}
   Double_t &Y() {return fY;}
   void GetKnot(Double_t &x, Double_t &y) const {x=fX; y=fY;}

   virtual Double_t Eval(Double_t) const {return fY;}

   private:
   void CopyPoly(TSplinePoly const &other);

   ClassDef(TSplinePoly,2) // Spline polynomial terms
};

inline TSplinePoly::TSplinePoly(TSplinePoly const &other)
:
   TObject(other), fX(0), fY(0)
{
  CopyPoly(other);
}


//______________________________________________________________________________
class TSplinePoly3 : public TSplinePoly
{
private:
   Double_t fB; // first order expansion coefficient :  fB*1! is the first derivative at x
   Double_t fC; // second order expansion coefficient : fC*2! is the second derivative at x
   Double_t fD; // third order expansion coefficient :  fD*3! is the third derivative at x

public:
   TSplinePoly3() :
      fB(0), fC(0), fD(0) {}
   TSplinePoly3(Double_t x, Double_t y, Double_t b, Double_t c, Double_t d) :
      TSplinePoly(x,y), fB(b), fC(c), fD(d) {}
   TSplinePoly3(TSplinePoly3 const &other);
   TSplinePoly3 &operator=(TSplinePoly3 const &other);

   Double_t &B() {return fB;}
   Double_t &C() {return fC;}
   Double_t &D() {return fD;}
   Double_t Eval(Double_t x) const {
      Double_t dx=x-fX;
      return (fY+dx*(fB+dx*(fC+dx*fD)));
   }
   Double_t Derivative(Double_t x) const {
      Double_t dx=x-fX;
      return (fB+dx*(2*fC+3*fD*dx));
   }

private:
   void CopyPoly(TSplinePoly3 const &other);

   ClassDef(TSplinePoly3,1)  // Third spline polynomial terms
};

inline TSplinePoly3::TSplinePoly3(TSplinePoly3 const &other)
   :
  TSplinePoly(other), fB(0), fC(0), fD(0)
{
   CopyPoly(other);
}

//______________________________________________________________________________
class TSplinePoly5 : public TSplinePoly
{
private:
   Double_t fB; // first order expansion coefficient :  fB*1! is the first derivative at x
   Double_t fC; // second order expansion coefficient : fC*2! is the second derivative at x
   Double_t fD; // third order expansion coefficient :  fD*3! is the third derivative at x
   Double_t fE; // fourth order expansion coefficient : fE*4! is the fourth derivative at x
   Double_t fF; // fifth order expansion coefficient :  fF*5! is the fifth derivative at x

public:
   TSplinePoly5() :
      fB(0), fC(0), fD(0), fE(0), fF(0) {}
   TSplinePoly5(Double_t x, Double_t y, Double_t b, Double_t c,
      Double_t d, Double_t e, Double_t f) :
      TSplinePoly(x,y), fB(b), fC(c), fD(d), fE(e), fF(f) {}
   TSplinePoly5(TSplinePoly5 const &other);
   TSplinePoly5 &operator=(TSplinePoly5 const &other);

   Double_t &B() {return fB;}
   Double_t &C() {return fC;}
   Double_t &D() {return fD;}
   Double_t &E() {return fE;}
   Double_t &F() {return fF;}
   Double_t Eval(Double_t x) const {
      Double_t dx=x-fX;
      return (fY+dx*(fB+dx*(fC+dx*(fD+dx*(fE+dx*fF)))));
   }
   Double_t Derivative(Double_t x) const{
      Double_t dx=x-fX;
      return (fB+dx*(2*fC+dx*(3*fD+dx*(4*fE+dx*(5*fF)))));
   }

private:
   void CopyPoly(TSplinePoly5 const &other);

   ClassDef(TSplinePoly5,1)  // Quintic spline polynomial terms
};

inline TSplinePoly5::TSplinePoly5(TSplinePoly5 const &other)
:
    TSplinePoly(other), fB(0), fC(0), fD(0), fE(0), fF(0)
{
  CopyPoly(other);
}


//______________________________________________________________________________
class TSpline3 : public TSpline
{
protected:
   TSplinePoly3  *fPoly;       //[fNp] Array of polynomial terms
   Double_t       fValBeg;     // Initial value of first or second derivative
   Double_t       fValEnd;     // End value of first or second derivative
   Int_t          fBegCond;    // 0=no beg cond, 1=first derivative, 2=second derivative
   Int_t          fEndCond;    // 0=no end cond, 1=first derivative, 2=second derivative

   void   BuildCoeff();
   void   SetCond(const char *opt);

public:
   TSpline3() : TSpline() , fPoly(0), fValBeg(0), fValEnd(0),
      fBegCond(-1), fEndCond(-1) {}
   TSpline3(const char *title,
            Double_t x[], Double_t y[], Int_t n, const char *opt=0,
            Double_t valbeg=0, Double_t valend=0);
   TSpline3(const char *title,
            Double_t xmin, Double_t xmax,
            Double_t y[], Int_t n, const char *opt=0,
            Double_t valbeg=0, Double_t valend=0);
   TSpline3(const char *title,
            Double_t x[], const TF1 *func, Int_t n, const char *opt=0,
            Double_t valbeg=0, Double_t valend=0);
   TSpline3(const char *title,
            Double_t xmin, Double_t xmax,
            const TF1 *func, Int_t n, const char *opt=0,
            Double_t valbeg=0, Double_t valend=0);
   TSpline3(const char *title,
            const TGraph *g, const char *opt=0,
            Double_t valbeg=0, Double_t valend=0);
   TSpline3(const TH1 *h, const char *opt=0,
            Double_t valbeg=0, Double_t valend=0);
   TSpline3(const TSpline3&);
   TSpline3& operator=(const TSpline3&);
   Int_t    FindX(Double_t x) const;
   Double_t Eval(Double_t x) const;
   Double_t Derivative(Double_t x) const;
   virtual ~TSpline3() {if (fPoly) delete [] fPoly;}
   void GetCoeff(Int_t i, Double_t &x, Double_t &y, Double_t &b,
                 Double_t &c, Double_t &d) {x=fPoly[i].X();y=fPoly[i].Y();
                  b=fPoly[i].B();c=fPoly[i].C();d=fPoly[i].D();}
   void GetKnot(Int_t i, Double_t &x, Double_t &y) const
      {x=fPoly[i].X(); y=fPoly[i].Y();}
   virtual  void     SaveAs(const char *filename,Option_t *option="") const;
   virtual  void     SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual  void     SetPoint(Int_t i, Double_t x, Double_t y);
   virtual  void     SetPointCoeff(Int_t i, Double_t b, Double_t c, Double_t d);
   static void Test();

   ClassDef (TSpline3,2)  // Class to create third natural splines
};


//______________________________________________________________________________
class TSpline5 : public TSpline
{
protected:
   TSplinePoly5  *fPoly;     //[fNp] Array of polynomial terms

   void BuildCoeff();
   void BoundaryConditions(const char *opt, Int_t &beg, Int_t &end,
                           const char *&cb1, const char *&ce1, const char *&cb2,
                           const char *&ce2);
   void SetBoundaries(Double_t b1, Double_t e1, Double_t b2, Double_t e2,
                      const char *cb1, const char *ce1, const char *cb2,
                      const char *ce2);
public:
   TSpline5() : TSpline() , fPoly(0) {}
   TSpline5(const char *title,
            Double_t x[], Double_t y[], Int_t n,
            const char *opt=0, Double_t b1=0, Double_t e1=0,
            Double_t b2=0, Double_t e2=0);
   TSpline5(const char *title,
            Double_t xmin, Double_t xmax,
            Double_t y[], Int_t n,
            const char *opt=0, Double_t b1=0, Double_t e1=0,
            Double_t b2=0, Double_t e2=0);
   TSpline5(const char *title,
            Double_t x[], const TF1 *func, Int_t n,
            const char *opt=0, Double_t b1=0, Double_t e1=0,
            Double_t b2=0, Double_t e2=0);
   TSpline5(const char *title,
            Double_t xmin, Double_t xmax,
            const TF1 *func, Int_t n,
            const char *opt=0, Double_t b1=0, Double_t e1=0,
            Double_t b2=0, Double_t e2=0);
   TSpline5(const char *title,
            const TGraph *g,
            const char *opt=0, Double_t b1=0, Double_t e1=0,
            Double_t b2=0, Double_t e2=0);
   TSpline5(const TH1 *h,
            const char *opt=0, Double_t b1=0, Double_t e1=0,
            Double_t b2=0, Double_t e2=0);
   TSpline5(const TSpline5&);
   TSpline5& operator=(const TSpline5&);
   Int_t    FindX(Double_t x) const;
   Double_t Eval(Double_t x) const;
   Double_t Derivative(Double_t x) const;
   virtual ~TSpline5() {if (fPoly) delete [] fPoly;}
   void GetCoeff(Int_t i, Double_t &x, Double_t &y, Double_t &b,
                 Double_t &c, Double_t &d, Double_t &e, Double_t &f)
      {x=fPoly[i].X();y=fPoly[i].Y();b=fPoly[i].B();
      c=fPoly[i].C();d=fPoly[i].D();
      e=fPoly[i].E();f=fPoly[i].F();}
   void GetKnot(Int_t i, Double_t &x, Double_t &y) const
      {x=fPoly[i].X(); y=fPoly[i].Y();}
   virtual  void     SaveAs(const char *filename,Option_t *option="") const;
   virtual  void     SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual  void     SetPoint(Int_t i, Double_t x, Double_t y);
   virtual  void     SetPointCoeff(Int_t i, Double_t b, Double_t c, Double_t d,
                                   Double_t e, Double_t f);
   static void Test();

   ClassDef (TSpline5,2) // Class to create quintic natural splines
};

#endif
