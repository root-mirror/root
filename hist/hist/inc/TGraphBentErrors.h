// @(#)root/hist:$Id$
// Author: Dave Morrison  30/06/2003

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGraphBentErrors
#define ROOT_TGraphBentErrors

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGraphBentErrors                                                     //
//                                                                      //
// a Graph with bent, asymmetric error bars                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGraph.h"

class TGraphBentErrors : public TGraph {

protected:
   Double_t    *fEXlow;        ///<[fNpoints] array of X low errors
   Double_t    *fEXhigh;       ///<[fNpoints] array of X high errors
   Double_t    *fEYlow;        ///<[fNpoints] array of Y low errors
   Double_t    *fEYhigh;       ///<[fNpoints] array of Y high errors

   Double_t    *fEXlowd;       ///<[fNpoints] array of X low displacements
   Double_t    *fEXhighd;      ///<[fNpoints] array of X high displacements
   Double_t    *fEYlowd;       ///<[fNpoints] array of Y low displacements
   Double_t    *fEYhighd;      ///<[fNpoints] array of Y high displacements

   virtual void       SwapPoints(Int_t pos1, Int_t pos2);

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
   TGraphBentErrors();
   TGraphBentErrors(Int_t n);
   TGraphBentErrors(Int_t n,
                    const Float_t *x, const Float_t *y,
                    const Float_t *exl=0, const Float_t *exh=0,
                    const Float_t *eyl=0, const Float_t *eyh=0,
                    const Float_t *exld=0, const Float_t *exhd=0,
                    const Float_t *eyld=0, const Float_t *eyhd=0);
   TGraphBentErrors(Int_t n,
                    const Double_t *x, const Double_t *y,
                    const Double_t *exl=0, const Double_t *exh=0,
                    const Double_t *eyl=0, const Double_t *eyh=0,
                    const Double_t *exld=0, const Double_t *exhd=0,
                    const Double_t *eyld=0, const Double_t *eyhd=0);
   TGraphBentErrors(const TGraphBentErrors &gr);
   virtual ~TGraphBentErrors();
   virtual void    Apply(TF1 *f);
   virtual void    ComputeRange(Double_t &xmin, Double_t &ymin,
                                Double_t &xmax, Double_t &ymax) const;
   Double_t        GetErrorX(Int_t bin)     const;
   Double_t        GetErrorY(Int_t bin)     const;
   Double_t        GetErrorXlow(Int_t bin)  const;
   Double_t        GetErrorXhigh(Int_t bin) const;
   Double_t        GetErrorYlow(Int_t bin)  const;
   Double_t        GetErrorYhigh(Int_t bin) const;
   Double_t       *GetEXlow()   const {return fEXlow;}
   Double_t       *GetEXhigh()  const {return fEXhigh;}
   Double_t       *GetEYlow()   const {return fEYlow;}
   Double_t       *GetEYhigh()  const {return fEYhigh;}
   Double_t       *GetEXlowd()  const {return fEXlowd;}
   Double_t       *GetEXhighd() const {return fEXhighd;}
   Double_t       *GetEYlowd()  const {return fEYlowd;}
   Double_t       *GetEYhighd() const {return fEYhighd;}
   virtual void    Print(Option_t *chopt="") const;
   virtual void    SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void    Scale(Double_t c1=1., Option_t *option="y"); // *MENU*
   virtual void    SetPointError(Double_t exl, Double_t exh,
                                 Double_t eyl, Double_t eyh,
                                 Double_t exld=0, Double_t exhd=0,
                                 Double_t eyld=0, Double_t eyhd=0); // *MENU*
   virtual void    SetPointError(Int_t i,
                                 Double_t exl, Double_t exh,
                                 Double_t eyl, Double_t eyh,
                                 Double_t exld=0, Double_t exhd=0,
                                 Double_t eyld=0, Double_t eyhd=0);

   ClassDef(TGraphBentErrors,1)  //A graph with bent, asymmetric error bars
};

inline Double_t **TGraphBentErrors::Allocate(Int_t size) {
   return AllocateArrays(10, size);
}

#endif
