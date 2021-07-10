// @(#)root/hist:$Id$
// Author: Rene Brun   03/03/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGraphAsymmErrors
#define ROOT_TGraphAsymmErrors


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGraphAsymmErrors                                                    //
//                                                                      //
// a Graph with asymmetric error bars                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGraph.h"

class TGraphAsymmErrors : public TGraph {

protected:
   Double_t    *fEXlow;        ///<[fNpoints] array of X low errors
   Double_t    *fEXhigh;       ///<[fNpoints] array of X high errors
   Double_t    *fEYlow;        ///<[fNpoints] array of Y low errors
   Double_t    *fEYhigh;       ///<[fNpoints] array of Y high errors

   virtual void    SwapPoints(Int_t pos1, Int_t pos2);

   virtual Double_t** Allocate(Int_t size);
   virtual void       CopyAndRelease(Double_t **newarrays,
                                     Int_t ibegin, Int_t iend, Int_t obegin);
   virtual Bool_t     CopyPoints(Double_t **arrays, Int_t ibegin, Int_t iend,
                                 Int_t obegin);
   Bool_t             CtorAllocate();
   virtual void       FillZero(Int_t begin, Int_t end,
                               Bool_t from_ctor = kTRUE);
   virtual Bool_t     DoMerge(const TGraph * g);

public:
   TGraphAsymmErrors();
   TGraphAsymmErrors(Int_t n);
   TGraphAsymmErrors(Int_t n, const Float_t *x, const Float_t *y, const Float_t *exl=0, const Float_t *exh=0, const Float_t *eyl=0, const Float_t *eyh=0);
   TGraphAsymmErrors(Int_t n, const Double_t *x, const Double_t *y, const Double_t *exl=0, const Double_t *exh=0, const Double_t *eyl=0, const Double_t *eyh=0);
   TGraphAsymmErrors(const TVectorF &vx, const TVectorF &vy, const TVectorF &vexl, const TVectorF &vexh, const TVectorF &veyl, const TVectorF &veyh);
   TGraphAsymmErrors(const TVectorD &vx, const TVectorD &vy, const TVectorD &vexl, const TVectorD &vexh, const TVectorD &veyl, const TVectorD &veyh);
   TGraphAsymmErrors(const TGraphAsymmErrors &gr);
   TGraphAsymmErrors& operator=(const TGraphAsymmErrors &gr);
   TGraphAsymmErrors(const TH1 *h);
   TGraphAsymmErrors(const TH1* pass, const TH1* total, Option_t *option="");
   TGraphAsymmErrors(const char *filename, const char *format="%lg %lg %lg %lg %lg %lg", Option_t *option="");

   virtual ~TGraphAsymmErrors();

   virtual void    Apply(TF1 *f);
   virtual void    BayesDivide(const TH1* pass, const TH1* total, Option_t *opt="");
   virtual void    Divide(const TH1* pass, const TH1* total, Option_t *opt="cp");
   virtual void    ComputeRange(Double_t &xmin, Double_t &ymin, Double_t &xmax, Double_t &ymax) const;
   Double_t        GetErrorX(Int_t bin)   const;
   Double_t        GetErrorY(Int_t bin)   const;
   Double_t        GetErrorXlow(Int_t i)  const;
   Double_t        GetErrorXhigh(Int_t i) const;
   Double_t        GetErrorYlow(Int_t i)  const;
   Double_t        GetErrorYhigh(Int_t i) const;
   Double_t       *GetEXlow()  const {return fEXlow;}
   Double_t       *GetEXhigh() const {return fEXhigh;}
   Double_t       *GetEYlow()  const {return fEYlow;}
   Double_t       *GetEYhigh() const {return fEYhigh;}
   virtual Int_t   Merge(TCollection* list);
   virtual void    Print(Option_t *chopt="") const;
   virtual void    SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void    Scale(Double_t c1=1., Option_t *option="y"); // *MENU*
   virtual void    SetPointError(Double_t exl, Double_t exh, Double_t eyl, Double_t eyh); // *MENU*
   virtual void    SetPointError(Int_t i, Double_t exl, Double_t exh, Double_t eyl, Double_t eyh);
   virtual void    SetPointEXlow(Int_t i, Double_t exl);
   virtual void    SetPointEXhigh(Int_t i, Double_t exh);
   virtual void    SetPointEYlow(Int_t i, Double_t eyl);
   virtual void    SetPointEYhigh(Int_t i, Double_t eyh);

   ClassDef(TGraphAsymmErrors,3)  //A graph with asymmetric error bars
};

#endif
