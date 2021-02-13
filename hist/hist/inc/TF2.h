// @(#)root/hist:$Id$
// Author: Rene Brun   23/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
// ---------------------------------- F2.h

#ifndef ROOT_TF2
#define ROOT_TF2



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TF2                                                                  //
//                                                                      //
// The Parametric 2-D function                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TF1.h"
#include "TArrayD.h"

class TF2 : public TF1 {

protected:
   Double_t  fYmin;        //Lower bound for the range in y
   Double_t  fYmax;        //Upper bound for the range in y
   Int_t     fNpy;         //Number of points along y used for the graphical representation
   TArrayD   fContour;     //Array to display contour levels

public:
   TF2();
   TF2(const char *name, const char *formula, Double_t xmin=0, Double_t xmax=1, Double_t ymin=0, Double_t ymax=1, Option_t * opt = nullptr);
#ifndef __CINT__
   TF2(const char *name, Double_t (*fcn)(Double_t *, Double_t *), Double_t xmin=0, Double_t xmax=1, Double_t ymin=0, Double_t ymax=1, Int_t npar=0,Int_t ndim = 2);
   TF2(const char *name, Double_t (*fcn)(const Double_t *, const Double_t *), Double_t xmin=0, Double_t xmax=1, Double_t ymin=0, Double_t ymax=1, Int_t npar=0, Int_t ndim = 2);
#endif

   // constructor using a functor

   TF2(const char *name, ROOT::Math::ParamFunctor f, Double_t xmin = 0, Double_t xmax = 1, Double_t ymin = 0, Double_t ymax = 1, Int_t npar = 0, Int_t ndim = 2);


   // Template constructors from a pointer to any C++ class of type PtrObj with a specific member function of type
   // MemFn.
   template <class PtrObj, typename MemFn>
   TF2(const char *name, const  PtrObj& p, MemFn memFn, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Int_t npar, Int_t ndim = 2) :
      TF1(name,p,memFn,xmin,xmax,npar,ndim),
   fYmin(ymin), fYmax(ymax), fNpy(30), fContour(0)
   {
      fNpx = 30;
   }
   /// backward compatible ctor
   template <class PtrObj, typename MemFn>
   TF2(const char *name, const  PtrObj& p, MemFn memFn, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Int_t npar, const char * , const char *) :
      TF1(name,p,memFn,xmin,xmax,npar,2),
   fYmin(ymin), fYmax(ymax), fNpy(30), fContour(0)
   {
      fNpx = 30;
   }

   // Template constructors from any  C++ callable object,  defining  the operator() (double * , double *)
   // and returning a double.
   template <typename Func>
   TF2(const char *name, Func f, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Int_t npar,Int_t ndim = 2) :
      TF1(name,f,xmin,xmax,npar,ndim),
   fYmin(ymin), fYmax(ymax), fNpy(30), fContour(0)
   {
      fNpx = 30;
   }
   /// backward compatible ctor
   template <typename Func>
   TF2(const char *name, Func f, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Int_t npar,const char *) :
      TF1(name,f,xmin,xmax,npar,2),
   fYmin(ymin), fYmax(ymax), fNpy(30), fContour(0)
   {
      fNpx = 30;
   }


   TF2(const TF2 &f2);
   TF2 &operator=(const TF2& rhs);
   virtual   ~TF2();
   virtual void     Copy(TObject &f2) const;
   virtual Int_t    DistancetoPrimitive(Int_t px, Int_t py);
   virtual void     Draw(Option_t *option="");
   virtual TF1     *DrawCopy(Option_t *option="") const;
   virtual TObject *DrawDerivative(Option_t * ="al") {return 0;}
   virtual TObject *DrawIntegral(Option_t * ="al")   {return 0;}
   //virtual void     DrawF2(const char *formula, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Option_t *option="");
   virtual void     ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual Int_t    GetContour(Double_t *levels=0);
   virtual Double_t GetContourLevel(Int_t level) const;
          Int_t     GetNpy() const {return fNpy;}
   virtual char    *GetObjectInfo(Int_t px, Int_t py) const;
   Double_t GetRandom(TRandom * rng = nullptr, Option_t * opt = nullptr);
   Double_t GetRandom(Double_t xmin, Double_t xmax, TRandom * rng = nullptr, Option_t * opt = nullptr);
   virtual void     GetRandom2(Double_t &xrandom, Double_t &yrandom, TRandom * rng = nullptr);
   using TF1::GetRange;
   virtual void     GetRange(Double_t &xmin, Double_t &ymin, Double_t &xmax, Double_t &ymax) const;
   virtual void     GetRange(Double_t &xmin, Double_t &ymin, Double_t &zmin, Double_t &xmax, Double_t &ymax, Double_t &zmax) const;
   virtual Double_t GetSave(const Double_t *x);
   virtual Double_t GetMinimumXY(Double_t &x, Double_t &y) const;
   virtual Double_t GetMaximumXY(Double_t &x, Double_t &y) const;
   using TF1::GetMinimum;
   using TF1::GetMaximum;
   virtual Double_t GetMinimum(Double_t *x ) const;
   virtual Double_t GetMaximum(Double_t *x ) const;
   virtual Double_t GetYmin() const {return fYmin;}
   virtual Double_t GetYmax() const {return fYmax;}
   using TF1::Integral;
   virtual Double_t Integral(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t epsrel=1.e-6);
   virtual Bool_t   IsInside(const Double_t *x) const;
   virtual TH1     *CreateHistogram();
   virtual void     Paint(Option_t *option="");
   virtual void     Save(Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Double_t zmin, Double_t zmax);
   virtual void     SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void     SetNpy(Int_t npy=100); // *MENU*
   virtual void     SetContour(Int_t nlevels=20, const Double_t *levels=0);
   virtual void     SetContourLevel(Int_t level, Double_t value);
   virtual void     SetRange(Double_t xmin, Double_t xmax);
   virtual void     SetRange(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax); // *MENU*
   virtual void     SetRange(Double_t xmin, Double_t ymin, Double_t zmin, Double_t xmax, Double_t ymax, Double_t zmax);

   //Moments
   virtual Double_t Moment2(Double_t nx, Double_t ax, Double_t bx, Double_t ny, Double_t ay, Double_t by, Double_t epsilon=0.000001);
   virtual Double_t CentralMoment2(Double_t nx, Double_t ax, Double_t bx, Double_t ny, Double_t ay, Double_t by, Double_t epsilon=0.000001);

   virtual Double_t Mean2X(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t epsilon=0.000001) {return Moment2(1,ax,bx,0,ay,by,epsilon);}
   virtual Double_t Mean2Y(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t epsilon=0.000001) {return Moment2(0,ax,bx,1,ay,by,epsilon);}

   virtual Double_t Variance2X(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t epsilon=0.000001) {return CentralMoment2(2,ax,bx,0,ay,by,epsilon);}
   virtual Double_t Variance2Y(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t epsilon=0.000001) {return CentralMoment2(0,ax,bx,2,ay,by,epsilon);}

   virtual Double_t Covariance2XY(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t epsilon=0.000001) {return CentralMoment2(1,ax,bx,1,ay,by,epsilon);}

protected:

   virtual Double_t FindMinMax(Double_t* x, bool findmax) const;

   ClassDef(TF2,4)  //The Parametric 2-D function
};

inline void TF2::SetRange(Double_t xmin, Double_t xmax)
   { TF1::SetRange(xmin, xmax); }
inline void TF2::SetRange(Double_t xmin, Double_t ymin, Double_t, Double_t xmax, Double_t ymax, Double_t)
   { SetRange(xmin, ymin, xmax, ymax); }

#endif
