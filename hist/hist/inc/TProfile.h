// @(#)root/hist:$Id$
// Author: Rene Brun   29/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TProfile
#define ROOT_TProfile


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TProfile                                                             //
//                                                                      //
// Profile histogram class.                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TH1.h"

class TProfileHelper;

enum EErrorType { kERRORMEAN = 0, kERRORSPREAD, kERRORSPREADI, kERRORSPREADG };

class TF1;

class TProfile : public TH1D {

public:
   friend class TProfileHelper;
   friend class TH1Merger;

protected:
   TArrayD fBinEntries;   // number of entries per bin
   EErrorType fErrorMode; // Option to compute errors
   Double_t fYmin;        // Lower limit in Y (if set)
   Double_t fYmax;        // Upper limit in Y (if set)
   Bool_t fScaling;       //! True when TProfile::Scale is called
   Double_t fTsumwy;      // Total Sum of weight*Y
   Double_t fTsumwy2;     // Total Sum of weight*Y*Y
   TArrayD fBinSumw2;     // Array of sum of squares of weights per bin

   static Bool_t fgApproximate; // bin error approximation option

   virtual Int_t    BufferFill(Double_t, Double_t) {return -2;} //may not use
   virtual Int_t    BufferFill(Double_t x, Double_t y, Double_t w);

   // helper methods for the Merge unification in TProfileHelper
   void SetBins(const Int_t* nbins, const Double_t* range) { SetBins(nbins[0], range[0], range[1]); };
   Int_t Fill(const Double_t* v) { return Fill(v[0], v[1], v[2]); };

   virtual Double_t RetrieveBinContent(Int_t bin) const { return (fBinEntries.fArray[bin] > 0) ? fArray[bin]/fBinEntries.fArray[bin] : 0; }
   //virtual void     UpdateBinContent(Int_t bin, Double_t content);
   virtual Double_t GetBinErrorSqUnchecked(Int_t bin) const { Double_t err = GetBinError(bin); return err*err; }

private:
   Int_t Fill(Double_t) { MayNotUse("Fill(Double_t)"); return -1;}
   void FillN(Int_t, const Double_t *, const Double_t *, Int_t) { MayNotUse("FillN(Int_t, Double_t*, Double_t*, Int_t)"); }
   Double_t *GetB()  {return &fBinEntries.fArray[0];}
   Double_t *GetB2() {return (fBinSumw2.fN ? &fBinSumw2.fArray[0] : 0 ); }
   Double_t *GetW()  {return &fArray[0];}
   Double_t *GetW2() {return &fSumw2.fArray[0];}
   void SetBins(Int_t, Double_t, Double_t, Int_t, Double_t, Double_t)
      { MayNotUse("SetBins(Int_t, Double_t, Double_t, Int_t, Double_t, Double_t"); }
   void SetBins(Int_t, const Double_t*, Int_t, const Double_t*)
      { MayNotUse("SetBins(Int_t, const Double_t*, Int_t, const Double_t*"); }
   void SetBins(Int_t, Double_t, Double_t, Int_t, Double_t, Double_t, Int_t, Double_t, Double_t)
      { MayNotUse("SetBins(Int_t, Double_t, Double_t, Int_t, Double_t, Double_t, Int_t, Double_t, Double_t"); }
   void SetBins(Int_t, const Double_t *, Int_t, const Double_t *, Int_t, const Double_t *)
      { MayNotUse("SetBins(Int_t, const Double_t*, Int_t, const Double_t*, Int_t, const Double_t*"); }

public:
   TProfile();
   TProfile(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup, Option_t *option="");
   TProfile(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup,Double_t ylow,Double_t yup,Option_t *option="");
   TProfile(const char *name,const char *title,Int_t nbinsx,const Float_t  *xbins, Option_t *option="");
   TProfile(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins, Option_t *option="");
   TProfile(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins, Double_t ylow,Double_t yup, Option_t *option="");
   TProfile(const TProfile &profile);
   TProfile &operator=(const TProfile &profile);
   virtual ~TProfile();
   virtual Bool_t   Add(TF1 *h1, Double_t c1=1, Option_t *option="");
   virtual Bool_t   Add(const TH1 *h1, Double_t c1=1);
   virtual Bool_t   Add(const TH1 *h1, const TH1 *h2, Double_t c1=1, Double_t c2=1); // *MENU*
   static  void     Approximate(Bool_t approx=kTRUE);
   virtual Int_t    BufferEmpty(Int_t action=0);
           void     BuildOptions(Double_t ymin, Double_t ymax, Option_t *option);
   virtual void     Copy(TObject &hnew) const;
   virtual Bool_t   Divide(TF1 *h1, Double_t c1=1);
   virtual Bool_t   Divide(const TH1 *h1);
   virtual Bool_t   Divide(const TH1 *h1, const TH1 *h2, Double_t c1=1, Double_t c2=1, Option_t *option=""); // *MENU*
   virtual void     ExtendAxis(Double_t x, TAxis *axis);
   virtual Int_t    Fill(Double_t x, Double_t y);
   virtual Int_t    Fill(const char *namex, Double_t y);
   virtual Int_t    Fill(Double_t x, Double_t y, Double_t w);
   virtual Int_t    Fill(const char *namex, Double_t y, Double_t w);
   virtual void     FillN(Int_t ntimes, const Double_t *x, const Double_t *y, const Double_t *w, Int_t stride=1);
   virtual Double_t GetBinContent(Int_t bin) const;
   virtual Double_t GetBinContent(Int_t bin, Int_t) const {return GetBinContent(bin);}
   virtual Double_t GetBinContent(Int_t bin, Int_t, Int_t) const {return GetBinContent(bin);}
   virtual Double_t GetBinError(Int_t bin) const;
   virtual Double_t GetBinError(Int_t bin, Int_t) const {return GetBinError(bin);}
   virtual Double_t GetBinError(Int_t bin, Int_t, Int_t) const {return GetBinError(bin);}
   virtual Double_t GetBinEntries(Int_t bin) const;
   virtual Double_t GetBinEffectiveEntries(Int_t bin) const;
   virtual TArrayD *GetBinSumw2() {return &fBinSumw2;}
   virtual const TArrayD *GetBinSumw2() const {return &fBinSumw2;}
   Option_t        *GetErrorOption() const;
   virtual void     GetStats(Double_t *stats) const;
   virtual Double_t GetYmin() const {return fYmin;}
   virtual Double_t GetYmax() const {return fYmax;}
   virtual void     LabelsDeflate(Option_t *axis="X");
   virtual void     LabelsInflate(Option_t *axis="X");
   virtual void     LabelsOption(Option_t *option="h", Option_t *axis="X");
   virtual Long64_t Merge(TCollection *list);
   virtual Bool_t   Multiply(TF1 *h1, Double_t c1=1);
   virtual Bool_t   Multiply(const TH1 *h1);
   virtual Bool_t   Multiply(const TH1 *h1, const TH1 *h2, Double_t c1=1, Double_t c2=1, Option_t *option=""); // *MENU*
           TH1D    *ProjectionX(const char *name="_px", Option_t *option="e") const;
   virtual void     PutStats(Double_t *stats);
           TH1     *Rebin(Int_t ngroup=2, const char*newname="", const Double_t *xbins=0);
   virtual void     Reset(Option_t *option="");
   virtual void     SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void     Scale(Double_t c1=1, Option_t *option="");
   virtual void     SetBinEntries(Int_t bin, Double_t w);
   virtual void     SetBins(Int_t nbins, Double_t xmin, Double_t xmax);
   virtual void     SetBins(Int_t nx, const Double_t *xbins);
   virtual void     SetBinsLength(Int_t n=-1);
   virtual void     SetBuffer(Int_t buffersize, Option_t *option="");
   virtual void     SetErrorOption(Option_t *option=""); // *MENU*
   virtual void     Sumw2(Bool_t flag = kTRUE);

   ClassDef(TProfile,7)  //Profile histogram class
};

#endif
