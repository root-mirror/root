// @(#)root/hist:$Id$
// Author: Simon Spies 18/02/19

/*************************************************************************
 * Copyright (C) 2018-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGraphMultiErrors
#define ROOT_TGraphMultiErrors

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGraphMultiErrors                                                    //
//                                                                      //
// a Graph with asymmetric error bars and multiple y errors             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGraph.h"
#include <vector>

class TArrayF;
class TArrayD;

class TGraphMultiErrors : public TGraph {

protected:
   Int_t fNYErrors;                     ///<  The amount of different y-errors
   Int_t fSumErrorsMode;                ///<  How y errors are summed: kOnlyFirst = Only First; kSquareSum = Squared Sum; kSum =
                                        ///<  Absolute Addition
   Double_t *fExL;                      ///<[fNpoints] array of X low errors
   Double_t *fExH;                      ///<[fNpoints] array of X high errors
   std::vector<TArrayD> fEyL;           ///<  Two dimensional array of Y low errors
   std::vector<TArrayD> fEyH;           ///<  Two dimensional array of Y high errors
   mutable Double_t *fEyLSum = nullptr; ///<! Array of summed Y low errors for fitting
   mutable Double_t *fEyHSum = nullptr; ///<! Array of summed Y high errors for fitting
   std::vector<TAttFill> fAttFill;      ///<  The AttFill attributes of the different errors
   std::vector<TAttLine> fAttLine;      ///<  The AttLine attributes of the different errors

   virtual Double_t **Allocate(Int_t size);
   Bool_t CtorAllocate();

   virtual void CopyAndRelease(Double_t **newarrays, Int_t ibegin, Int_t iend, Int_t obegin);
   virtual Bool_t CopyPoints(Double_t **arrays, Int_t ibegin, Int_t iend, Int_t obegin);
   virtual void FillZero(Int_t begin, Int_t end, Bool_t from_ctor = kTRUE);

   void CalcYErrorsSum() const;
   virtual Bool_t DoMerge(const TGraph *tg);
   virtual void SwapPoints(Int_t pos1, Int_t pos2);

public:
   enum ESummationModes {
      kOnlyFirst = 0, ///< Only take errors from first dimension
      kSquareSum = 1, ///< Calculate the square sum of all errors
      kAbsSum = 2     ///< Calculate the absolute sum of all errors
   };

   TGraphMultiErrors();
   TGraphMultiErrors(const Char_t *name, const Char_t *title);
   TGraphMultiErrors(Int_t np, Int_t ne = 1);
   TGraphMultiErrors(const Char_t *name, const Char_t *title, Int_t np, Int_t ne = 1);

   TGraphMultiErrors(Int_t np, const Float_t *x, const Float_t *y, const Float_t *exL = nullptr,
                     const Float_t *exH = nullptr, const Float_t *eyL = nullptr, const Float_t *eyH = nullptr,
                     Int_t m = kOnlyFirst);
   TGraphMultiErrors(const Char_t *name, const Char_t *title, Int_t np, const Float_t *x, const Float_t *y,
                     const Float_t *exL = nullptr, const Float_t *exH = nullptr, const Float_t *eyL = nullptr,
                     const Float_t *eyH = nullptr, Int_t m = kOnlyFirst);
   TGraphMultiErrors(Int_t np, const Double_t *x, const Double_t *y, const Double_t *exL = nullptr,
                     const Double_t *exH = nullptr, const Double_t *eyL = nullptr, const Double_t *eyH = nullptr,
                     Int_t m = kOnlyFirst);
   TGraphMultiErrors(const Char_t *name, const Char_t *title, Int_t np, const Double_t *x, const Double_t *y,
                     const Double_t *exL = nullptr, const Double_t *exH = nullptr, const Double_t *eyL = nullptr,
                     const Double_t *eyH = nullptr, Int_t m = kOnlyFirst);

   TGraphMultiErrors(Int_t np, Int_t ne, const Float_t *x, const Float_t *y, const Float_t *exL, const Float_t *exH,
                     std::vector<std::vector<Float_t>> eyL, std::vector<std::vector<Float_t>> eyH,
                     Int_t m = kOnlyFirst);
   TGraphMultiErrors(const Char_t *name, const Char_t *title, Int_t np, Int_t ne, const Float_t *x, const Float_t *y,
                     const Float_t *exL, const Float_t *exH, std::vector<std::vector<Float_t>> eyL,
                     std::vector<std::vector<Float_t>> eyH, Int_t m = kOnlyFirst);
   TGraphMultiErrors(Int_t np, Int_t ne, const Double_t *x, const Double_t *y, const Double_t *exL, const Double_t *exH,
                     std::vector<std::vector<Double_t>> eyL, std::vector<std::vector<Double_t>> eyH,
                     Int_t m = kOnlyFirst);
   TGraphMultiErrors(const Char_t *name, const Char_t *title, Int_t np, Int_t ne, const Double_t *x, const Double_t *y,
                     const Double_t *exL, const Double_t *exH, std::vector<std::vector<Double_t>> eyL,
                     std::vector<std::vector<Double_t>> eyH, Int_t m = kOnlyFirst);

   TGraphMultiErrors(Int_t np, Int_t ne, const Float_t *x, const Float_t *y, const Float_t *exL, const Float_t *exH,
                     std::vector<TArrayF> eyL, std::vector<TArrayF> eyH, Int_t m = kOnlyFirst);
   TGraphMultiErrors(const Char_t *name, const Char_t *title, Int_t np, Int_t ne, const Float_t *x, const Float_t *y,
                     const Float_t *exL, const Float_t *exH, std::vector<TArrayF> eyL, std::vector<TArrayF> eyH,
                     Int_t m = kOnlyFirst);
   TGraphMultiErrors(Int_t np, Int_t ne, const Double_t *x, const Double_t *y, const Double_t *exL, const Double_t *exH,
                     std::vector<TArrayD> eyL, std::vector<TArrayD> eyH, Int_t m = kOnlyFirst);
   TGraphMultiErrors(const Char_t *name, const Char_t *title, Int_t np, Int_t ne, const Double_t *x, const Double_t *y,
                     const Double_t *exL, const Double_t *exH, std::vector<TArrayD> eyL, std::vector<TArrayD> eyH,
                     Int_t m = kOnlyFirst);

   TGraphMultiErrors(const TVectorF &tvX, const TVectorF &tvY, const TVectorF &tvExL, const TVectorF &tvExH,
                     const TVectorF &tvEyL, const TVectorF &tvEyH, Int_t m = kOnlyFirst);
   TGraphMultiErrors(const TVectorD &tvX, const TVectorD &tvY, const TVectorD &tvExL, const TVectorD &tvExH,
                     const TVectorD &tvEyL, const TVectorD &tvEyH, Int_t m = kOnlyFirst);

   TGraphMultiErrors(Int_t ne, const TVectorF &tvX, const TVectorF &tvY, const TVectorF &tvExL, const TVectorF &tvExH,
                     const TVectorF *tvEyL, const TVectorF *tvEyH, Int_t m = kOnlyFirst);
   TGraphMultiErrors(Int_t ne, const TVectorD &tvX, const TVectorD &tvY, const TVectorD &tvExL, const TVectorD &tvExH,
                     const TVectorD *tvEyL, const TVectorD *tvEyH, Int_t m = kOnlyFirst);

   TGraphMultiErrors(const TGraphMultiErrors &tgme);
   TGraphMultiErrors &operator=(const TGraphMultiErrors &tgme);

   TGraphMultiErrors(const TH1 *th, Int_t ne = 1);
   TGraphMultiErrors(const TH1 *pass, const TH1 *total, Int_t ne = 1, Option_t *option = "");

   virtual ~TGraphMultiErrors();

   virtual void AddYError(Int_t np, const Double_t *eyL = nullptr, const Double_t *eyH = nullptr);
   virtual void Apply(TF1 *f);
   virtual void BayesDivide(const TH1 *pass, const TH1 *total, Option_t *opt = "");
   void Divide(const TH1 *pass, const TH1 *total, Option_t *opt = "cp");
   virtual void ComputeRange(Double_t &xmin, Double_t &ymin, Double_t &xmax, Double_t &ymax) const;
   virtual void DeleteYError(Int_t e);

   virtual Double_t GetErrorX(Int_t i) const;
   virtual Double_t GetErrorY(Int_t i) const;
   virtual Double_t GetErrorY(Int_t i, Int_t e) const;

   virtual Double_t GetErrorXlow(Int_t i) const;
   virtual Double_t GetErrorXhigh(Int_t i) const;
   virtual Double_t GetErrorYlow(Int_t i) const;
   virtual Double_t GetErrorYhigh(Int_t i) const;
   virtual Double_t GetErrorYlow(Int_t i, Int_t e) const;
   virtual Double_t GetErrorYhigh(Int_t i, Int_t e) const;

   virtual Double_t *GetEXlow() const { return fExL; }
   virtual Double_t *GetEXhigh() const { return fExH; }
   virtual Double_t *GetEYlow() const;
   virtual Double_t *GetEYhigh() const;
   virtual Double_t *GetEYlow(Int_t e);
   virtual Double_t *GetEYhigh(Int_t e);

   virtual TAttFill *GetAttFill(Int_t e);
   virtual TAttLine *GetAttLine(Int_t e);

   using TAttFill::GetFillColor;
   using TAttFill::GetFillStyle;

   virtual Color_t GetFillColor(Int_t e) const;
   virtual Style_t GetFillStyle(Int_t e) const;

   using TAttLine::GetLineColor;
   using TAttLine::GetLineStyle;
   using TAttLine::GetLineWidth;

   virtual Color_t GetLineColor(Int_t e) const;
   virtual Style_t GetLineStyle(Int_t e) const;
   virtual Width_t GetLineWidth(Int_t e) const;

   Int_t GetSumErrorsMode() const { return fSumErrorsMode; }
   Int_t GetNYErrors() const { return fNYErrors; }

   virtual void Print(Option_t *chopt = "") const;
   virtual void SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void Scale(Double_t c1=1., Option_t *option="y"); // *MENU*

   virtual void SetPointError(Double_t exL, Double_t exH, Double_t eyL1, Double_t eyH1, Double_t eyL2 = 0.,
                              Double_t eyH2 = 0., Double_t eyL3 = 0., Double_t eyH3 = 0.); // *MENU*
   virtual void SetPointError(Int_t i, Int_t ne, Double_t exL, Double_t exH, const Double_t *eyL, const Double_t *eyH);

   virtual void SetPointEX(Int_t i, Double_t exL, Double_t exH);
   virtual void SetPointEXlow(Int_t i, Double_t exL);
   virtual void SetPointEXhigh(Int_t i, Double_t exH);
   virtual void SetPointEY(Int_t i, Int_t ne, const Double_t *eyL, const Double_t *eyH);
   virtual void SetPointEYlow(Int_t i, Int_t ne, const Double_t *eyL);
   virtual void SetPointEYhigh(Int_t i, Int_t ne, const Double_t *eyH);
   virtual void SetPointEY(Int_t i, Int_t e, Double_t eyL, Double_t eyH);
   virtual void SetPointEYlow(Int_t i, Int_t e, Double_t eyL);
   virtual void SetPointEYhigh(Int_t i, Int_t e, Double_t eyH);

   virtual void SetEY(Int_t e, Int_t np, const Double_t *eyL, const Double_t *eyH);
   virtual void SetEYlow(Int_t e, Int_t np, const Double_t *eyL);
   virtual void SetEYhigh(Int_t e, Int_t np, const Double_t *eyH);

   virtual void SetSumErrorsMode(Int_t m);

   virtual void SetAttFill(Int_t e, TAttFill *taf);
   virtual void SetAttLine(Int_t e, TAttLine *tal);

   using TAttFill::SetFillColor;
   using TAttFill::SetFillColorAlpha;
   using TAttFill::SetFillStyle;

   virtual void SetFillColor(Int_t e, Color_t fcolor);
   virtual void SetFillColorAlpha(Int_t e, Color_t fcolor, Float_t falpha);
   virtual void SetFillStyle(Int_t e, Style_t fstyle);

   using TAttLine::SetLineColor;
   using TAttLine::SetLineColorAlpha;
   using TAttLine::SetLineStyle;
   using TAttLine::SetLineWidth;

   virtual void SetLineColor(Int_t e, Color_t lcolor);
   virtual void SetLineColorAlpha(Int_t e, Color_t lcolor, Float_t lalpha);
   virtual void SetLineStyle(Int_t e, Style_t lstyle);
   virtual void SetLineWidth(Int_t e, Width_t lwidth);

   ClassDef(TGraphMultiErrors, 1) // A Graph with asymmetric error bars and multiple y error dimensions
};

#endif // ROOT_TGraphMultiErrors
