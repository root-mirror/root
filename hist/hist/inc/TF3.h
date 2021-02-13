// @(#)root/hist:$Id$
// Author: Rene Brun   27/10/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
// ---------------------------------- F3.h

#ifndef ROOT_TF3
#define ROOT_TF3



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TF3                                                                  //
//                                                                      //
// The Parametric 3-D function                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TF2.h"

class TF3 : public TF2 {

protected:
   Double_t  fZmin;        //Lower bound for the range in z
   Double_t  fZmax;        //Upper bound for the range in z
   Int_t     fNpz;         //Number of points along z used for the graphical representation
   Bool_t    fClipBoxOn{kFALSE}; //! is clip box on
   Double_t  fClipBox[3];        //! coordinates of clipbox
public:
   TF3();
   TF3(const char *name, const char *formula, Double_t xmin=0, Double_t xmax=1, Double_t ymin=0,
       Double_t ymax=1, Double_t zmin=0, Double_t zmax=1, Option_t * opt = nullptr);
#ifndef __CINT__
   TF3(const char *name, Double_t (*fcn)(Double_t *, Double_t *), Double_t xmin=0, Double_t xmax=1, Double_t ymin=0,
       Double_t ymax=1, Double_t zmin=0, Double_t zmax=1, Int_t npar=0, Int_t ndim = 3);
   TF3(const char *name, Double_t (*fcn)(const Double_t *, const Double_t *), Double_t xmin=0, Double_t xmax=1, Double_t ymin=0,
       Double_t ymax=1, Double_t zmin=0, Double_t zmax=1, Int_t npar=0, Int_t ndim = 3);
#endif

   // constructor using a functor

   TF3(const char *name, ROOT::Math::ParamFunctor f, Double_t xmin = 0, Double_t xmax = 1, Double_t ymin = 0, Double_t ymax = 1, Double_t zmin=0, Double_t zmax=1, Int_t npar = 0, Int_t ndim = 3);


   // Template constructors from a pointer to any C++ class of type PtrObj with a specific member function of type
   // MemFn.
   template <class PtrObj, typename MemFn>
   TF3(const char *name, const  PtrObj& p, MemFn memFn, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Double_t zmin, Double_t zmax, Int_t npar,
       Int_t ndim = 3) :
      TF2(name,p,memFn,xmin,xmax,ymin,ymax,npar,ndim),
      fZmin(zmin), fZmax(zmax), fNpz(30)
   {   }
   /// Backward compatible ctor
   template <class PtrObj, typename MemFn>
   TF3(const char *name, const  PtrObj& p, MemFn memFn, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Double_t zmin, Double_t zmax, Int_t npar,
       const char * , const char *  ) :
      TF2(name,p,memFn,xmin,xmax,ymin,ymax,npar,3),
      fZmin(zmin), fZmax(zmax), fNpz(30)
   {   }
   // Template constructors from any  C++ callable object,  defining  the operator() (double * , double *)
   // and returning a double.
   template <typename Func>
   TF3(const char *name, Func f, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Double_t zmin, Double_t zmax, Int_t npar,
       Int_t ndim = 3 ) :
      TF2(name,f,xmin,xmax,ymin,ymax,npar,ndim),
      fZmin(zmin), fZmax(zmax), fNpz(30)
   { }
   /// backward compatible ctor
   template <typename Func>
   TF3(const char *name, Func f, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Double_t zmin, Double_t zmax, Int_t npar,
       const char *  ) :
      TF2(name,f,xmin,xmax,ymin,ymax,npar,3),
      fZmin(zmin), fZmax(zmax), fNpz(30)
   { }

   TF3(const TF3 &f3);
   TF3& operator=(const TF3 &rhs);
   virtual   ~TF3();
   virtual void     Copy(TObject &f3) const;
   virtual Int_t    DistancetoPrimitive(Int_t px, Int_t py);
   virtual void     Draw(Option_t *option="");
   virtual TObject *DrawDerivative(Option_t * ="al") {return 0;}
   virtual TObject *DrawIntegral(Option_t * ="al")   {return 0;}
   virtual void     ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual Double_t GetMinimumXYZ(Double_t &x, Double_t &y, Double_t &z);
   virtual Double_t GetMaximumXYZ(Double_t &x, Double_t &y, Double_t &z);
          Int_t     GetNpz() const {return fNpz;}
   virtual void     GetRandom3(Double_t &xrandom, Double_t &yrandom, Double_t &zrandom, TRandom * rng = nullptr);
   using TF1::GetRange;
   virtual void     GetRange(Double_t &xmin, Double_t &xmax) const;
   virtual void     GetRange(Double_t &xmin, Double_t &ymin, Double_t &xmax, Double_t &ymax) const ;
   virtual void     GetRange(Double_t &xmin, Double_t &ymin, Double_t &zmin, Double_t &xmax, Double_t &ymax, Double_t &zmax) const;
   virtual Double_t GetSave(const Double_t *x);
   virtual Double_t GetZmin() const {return fZmin;}
   virtual Double_t GetZmax() const {return fZmax;}
   using TF2::Integral;
   virtual Double_t Integral(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t az, Double_t bz, Double_t epsrel=1.e-6);
   virtual Bool_t   IsInside(const Double_t *x) const;
   virtual TH1     *CreateHistogram();
   virtual void     Paint(Option_t *option="");
   virtual void     Save(Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Double_t zmin, Double_t zmax);
   virtual void     SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void     SetClippingBoxOff(); // *MENU*
   virtual Bool_t   GetClippingBoxOn() const { return fClipBoxOn; }
   virtual void     SetClippingBoxOn(Double_t xclip=0, Double_t yclip=0, Double_t zclip=0); // *MENU*
   virtual const Double_t *GetClippingBox() const { return fClipBoxOn ? fClipBox : nullptr; }
   virtual void     SetNpz(Int_t npz=30);
   virtual void     SetRange(Double_t xmin, Double_t xmax);
   virtual void     SetRange(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax);
   virtual void     SetRange(Double_t xmin, Double_t ymin, Double_t zmin, Double_t xmax, Double_t ymax, Double_t zmax); // *MENU*

   //Moments
   virtual Double_t Moment3(Double_t nx, Double_t ax, Double_t bx, Double_t ny, Double_t ay, Double_t by, Double_t nz, Double_t az, Double_t bz, Double_t epsilon=0.000001);
   virtual Double_t CentralMoment3(Double_t nx, Double_t ax, Double_t bx, Double_t ny, Double_t ay, Double_t by, Double_t nz, Double_t az, Double_t bz, Double_t epsilon=0.000001);

   virtual Double_t Mean3X(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t az, Double_t bz, Double_t epsilon=0.000001) {return Moment3(1,ax,bx,0,ay,by,0,az,bz,epsilon);}
   virtual Double_t Mean3Y(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t az, Double_t bz, Double_t epsilon=0.000001) {return Moment3(0,ax,bx,1,ay,by,0,az,bz,epsilon);}
   virtual Double_t Mean3Z(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t az, Double_t bz, Double_t epsilon=0.000001) {return Moment3(0,ax,bx,0,ay,by,1,az,bz,epsilon);}

   virtual Double_t Variance3X(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t az, Double_t bz, Double_t epsilon=0.000001) {return CentralMoment3(2,ax,bx,0,ay,by,0,az,bz,epsilon);}
   virtual Double_t Variance3Y(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t az, Double_t bz, Double_t epsilon=0.000001) {return CentralMoment3(0,ax,bx,2,ay,by,0,az,bz,epsilon);}
   virtual Double_t Variance3Z(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t az, Double_t bz, Double_t epsilon=0.000001) {return CentralMoment3(0,ax,bx,0,ay,by,2,az,bz,epsilon);}

   virtual Double_t Covariance3XY(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t az, Double_t bz, Double_t epsilon=0.000001) {return CentralMoment3(1,ax,bx,1,ay,by,0,az,bz,epsilon);}
   virtual Double_t Covariance3XZ(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t az, Double_t bz, Double_t epsilon=0.000001) {return CentralMoment3(1,ax,bx,0,ay,by,1,az,bz,epsilon);}
   virtual Double_t Covariance3YZ(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t az, Double_t bz, Double_t epsilon=0.000001) {return CentralMoment3(0,ax,bx,1,ay,by,1,az,bz,epsilon);}

protected:

   virtual Double_t FindMinMax(Double_t* x, bool findmax) const;

   ClassDef(TF3,3)  //The Parametric 3-D function
};

inline void TF3::GetRange(Double_t &xmin, Double_t &xmax) const
   { TF2::GetRange(xmin, xmax); }
inline void TF3::GetRange(Double_t &xmin, Double_t &ymin, Double_t &xmax, Double_t &ymax) const
   { TF2::GetRange(xmin, ymin, xmax, ymax); }
inline void TF3::SetRange(Double_t xmin, Double_t xmax)
   { TF2::SetRange(xmin, xmax); }
inline void TF3::SetRange(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax)
   { TF2::SetRange(xmin, ymin, xmax, ymax); }

#endif
