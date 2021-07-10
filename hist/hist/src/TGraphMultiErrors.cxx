// @(#)root/hist:$Id$
// Author: Simon Spies 18/02/19

/*************************************************************************
 * Copyright (C) 2018-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TStyle.h"
#include "TVirtualPad.h"
#include "TEfficiency.h"
#include "Riostream.h"

#include "TArrayD.h"
#include "TVector.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "Math/QuantFuncMathCore.h"

#include "TGraphMultiErrors.h"

ClassImp(TGraphMultiErrors)

/** \class TGraphMultiErrors
    \ingroup Hist
TGraph with asymmetric error bars and multiple y error dimensions.

The TGraphMultiErrors painting is performed thanks to the TGraphPainter
class. All details about the various painting options are given in this class.

The picture below gives an example:

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1", "A Simple Graph with multiple y-errors", 200, 10, 700, 500);
   c1->SetGrid();
   c1->GetFrame()->SetBorderSize(12);
   const Int_t np = 5;
   Double_t x[np]       = {0, 1, 2, 3, 4};
   Double_t y[np]       = {0, 2, 4, 1, 3};
   Double_t exl[np]     = {0.3, 0.3, 0.3, 0.3, 0.3};
   Double_t exh[np]     = {0.3, 0.3, 0.3, 0.3, 0.3};
   Double_t eylstat[np] = {1, 0.5, 1, 0.5, 1};
   Double_t eyhstat[np] = {0.5, 1, 0.5, 1, 2};
   Double_t eylsys[np]  = {0.5, 0.4, 0.8, 0.3, 1.2};
   Double_t eyhsys[np]  = {0.6, 0.7, 0.6, 0.4, 0.8};
   auto gme = new TGraphMultiErrors("gme", "TGraphMultiErrors Example", np, x, y, exl, exh, eylstat, eyhstat);
   gme->AddYError(np, eylsys, eyhsys);
   gme->SetMarkerStyle(20);
   gme->SetLineColor(kRed);
   gme->GetAttLine(0)->SetLineColor(kRed);
   gme->GetAttLine(1)->SetLineColor(kBlue);
   gme->GetAttFill(1)->SetFillStyle(0);
   gme->Draw("APS ; Z ; 5 s=0.5");
}
End_Macro
*/
////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors default constructor.

TGraphMultiErrors::TGraphMultiErrors()
   : TGraph(), fNYErrors(0), fSumErrorsMode(TGraphMultiErrors::kOnlyFirst), fExL(nullptr), fExH(nullptr)
{
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors default constructor with name and title.

TGraphMultiErrors::TGraphMultiErrors(const Char_t *name, const Char_t *title) : TGraphMultiErrors()
{
   SetNameTitle(name, title);
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors normal constructor with np points and ne y-errors.
///
/// All values are initialized to 0.

TGraphMultiErrors::TGraphMultiErrors(Int_t np, Int_t ne)
   : TGraph(np), fNYErrors(ne), fSumErrorsMode(TGraphMultiErrors::kOnlyFirst)
{
   CtorAllocate();
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors normal constructor with `name`, `title`, `np` points and `ne` y-errors.
///
/// All values are initialized to 0.

TGraphMultiErrors::TGraphMultiErrors(const Char_t *name, const Char_t *title, Int_t np, Int_t ne)
   : TGraphMultiErrors(np, ne)
{
   SetNameTitle(name, title);
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors normal constructor with `np` points and a single y-error.
///
/// The signature of this constructor is equal to the corresponding constructor of TGraphAsymmErrors.
/// If `exL`,`exH` or `eyL`,`exH` are NULL, the corresponding values are preset to zero.

TGraphMultiErrors::TGraphMultiErrors(Int_t np, const Float_t *x, const Float_t *y, const Float_t *exL,
                                     const Float_t *exH, const Float_t *eyL, const Float_t *eyH, Int_t m)
   : TGraph(np, x, y), fNYErrors(1), fSumErrorsMode(m)
{
   if (!CtorAllocate())
      return;

   for (Int_t i = 0; i < fNpoints; i++) {
      if (exL)
         fExL[i] = exL[i];
      else
         fExL[i] = 0.;
      if (exH)
         fExH[i] = exH[i];
      else
         fExH[i] = 0.;
      if (eyL)
         fEyL[0][i] = eyL[i];
      else
         fEyL[0][i] = 0.;
      if (eyH)
         fEyH[0][i] = eyH[i];
      else
         fEyH[0][i] = 0.;
   }

   CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors normal constructor with `name`, `title`, `np` points and a single y-error.
///
/// If `exL`,`exH` or `eyL`,`eyH` are NULL, the corresponding values are preset to zero.

TGraphMultiErrors::TGraphMultiErrors(const Char_t *name, const Char_t *title, Int_t np, const Float_t *x,
                                     const Float_t *y, const Float_t *exL, const Float_t *exH, const Float_t *eyL,
                                     const Float_t *eyH, Int_t m)
   : TGraphMultiErrors(np, x, y, exL, exH, eyL, eyH, m)
{
   SetNameTitle(name, title);
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors normal constructor with `np` points and a single y-error.
///
/// The signature of this constructor is equal to the corresponding constructor of TGraphAsymmErrors.
/// If `exL`,`exH` or `eyL`,`exH` are NULL, the corresponding values are preset to zero.

TGraphMultiErrors::TGraphMultiErrors(Int_t np, const Double_t *x, const Double_t *y, const Double_t *exL,
                                     const Double_t *exH, const Double_t *eyL, const Double_t *eyH, Int_t m)
   : TGraph(np, x, y), fNYErrors(1), fSumErrorsMode(m)
{
   if (!CtorAllocate())
      return;

   Int_t n = fNpoints * sizeof(Double_t);

   if (exL)
      memcpy(fExL, exL, n);
   else
      memset(fExL, 0, n);
   if (exH)
      memcpy(fExH, exH, n);
   else
      memset(fExH, 0, n);

   if (eyL)
      fEyL[0].Set(fNpoints, eyL);
   else
      fEyL[0].Reset(0.);

   if (eyH)
      fEyH[0].Set(fNpoints, eyH);
   else
      fEyH[0].Reset(0.);

   CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors normal constructor with name, title, `np` points and a single y-error.
///
/// If `exL`,`exH` or `eyL`,`exH` are NULL, the corresponding values are preset to zero.

TGraphMultiErrors::TGraphMultiErrors(const Char_t *name, const Char_t *title, Int_t np, const Double_t *x,
                                     const Double_t *y, const Double_t *exL, const Double_t *exH, const Double_t *eyL,
                                     const Double_t *eyH, Int_t m)
   : TGraphMultiErrors(np, x, y, exL, exH, eyL, eyH, m)
{
   SetNameTitle(name, title);
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors normal constructor with `np` points and `ne` y-errors.
///
/// If `exL`,`exH` are NULL, the corresponding values are preset to zero.
/// The multiple y-errors are passed as std::vectors of std::vectors.

TGraphMultiErrors::TGraphMultiErrors(Int_t np, Int_t ne, const Float_t *x, const Float_t *y, const Float_t *exL,
                                     const Float_t *exH, std::vector<std::vector<Float_t>> eyL,
                                     std::vector<std::vector<Float_t>> eyH, Int_t m)
   : TGraph(np, x, y), fNYErrors(ne), fSumErrorsMode(m)
{
   if (!CtorAllocate())
      return;

   for (Int_t i = 0; i < fNpoints; i++) {
      if (exL)
         fExL[i] = exL[i];
      else
         fExL[i] = 0.;
      if (exH)
         fExH[i] = exH[i];
      else
         fExH[i] = 0.;

      for (Int_t j = 0; j < fNYErrors; j++) {
         if (Int_t(eyL.size()) > j && Int_t(eyL[j].size()) > i)
            fEyL[j][i] = eyL[j][i];
         else
            fEyL[j][i] = 0.;
         if (Int_t(eyH.size()) > j && Int_t(eyH[j].size()) > i)
            fEyH[j][i] = eyH[j][i];
         else
            fEyH[j][i] = 0.;
      }
   }

   CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors normal constructor with name, title, `np` points and `ne` y-errors.
///
/// If `exL`,`exH` are NULL, the corresponding values are preset to zero.
/// The multiple y-errors are passed as std::vectors of std::vectors.

TGraphMultiErrors::TGraphMultiErrors(const Char_t *name, const Char_t *title, Int_t np, Int_t ne, const Float_t *x,
                                     const Float_t *y, const Float_t *exL, const Float_t *exH,
                                     std::vector<std::vector<Float_t>> eyL, std::vector<std::vector<Float_t>> eyH,
                                     Int_t m)
   : TGraphMultiErrors(np, ne, x, y, exL, exH, eyL, eyH, m)
{
   SetNameTitle(name, title);
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors normal constructor with `np` points and `ne` y-errors.
///
/// If `exL`,`exH` are NULL, the corresponding values are preset to zero.
/// The multiple y-errors are passed as std::vectors of std::vectors.

TGraphMultiErrors::TGraphMultiErrors(Int_t np, Int_t ne, const Double_t *x, const Double_t *y, const Double_t *exL,
                                     const Double_t *exH, std::vector<std::vector<Double_t>> eyL,
                                     std::vector<std::vector<Double_t>> eyH, Int_t m)
   : TGraph(np, x, y), fNYErrors(ne), fSumErrorsMode(m)
{
   if (!CtorAllocate())
      return;

   Int_t n = fNpoints * sizeof(Double_t);

   if (exL)
      memcpy(fExL, exL, n);
   else
      memset(fExL, 0, n);
   if (exH)
      memcpy(fExH, exH, n);
   else
      memset(fExH, 0, n);

   for (Int_t i = 0; i < fNpoints; i++) {
      for (Int_t j = 0; j < fNYErrors; j++) {
         if (Int_t(eyL.size()) > j && Int_t(eyL[j].size()) > i)
            fEyL[j][i] = eyL[j][i];
         else
            fEyL[j][i] = 0.;
         if (Int_t(eyH.size()) > j && Int_t(eyH[j].size()) > i)
            fEyH[j][i] = eyH[j][i];
         else
            fEyH[j][i] = 0.;
      }
   }

   CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors normal constructor with `name`, `title`, `np` points and `ne` y-errors.
///
/// If `exL`,`exH` are NULL, the corresponding values are preset to zero.
/// The multiple y-errors are passed as std::vectors of std::vectors.

TGraphMultiErrors::TGraphMultiErrors(const Char_t *name, const Char_t *title, Int_t np, Int_t ne, const Double_t *x,
                                     const Double_t *y, const Double_t *exL, const Double_t *exH,
                                     std::vector<std::vector<Double_t>> eyL, std::vector<std::vector<Double_t>> eyH,
                                     Int_t m)
   : TGraphMultiErrors(np, ne, x, y, exL, exH, eyL, eyH, m)
{
   SetNameTitle(name, title);
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors normal constructor with `np` points and `ne` y-errors.
///
/// If `exL`,`exH` are NULL, the corresponding values are preset to zero.
/// The multiple y-errors are passed as std::vectors of TArrayF objects.

TGraphMultiErrors::TGraphMultiErrors(Int_t np, Int_t ne, const Float_t *x, const Float_t *y, const Float_t *exL,
                                     const Float_t *exH, std::vector<TArrayF> eyL, std::vector<TArrayF> eyH, Int_t m)
   : TGraph(np, x, y), fNYErrors(ne), fSumErrorsMode(m)
{
   if (!CtorAllocate())
      return;

   for (Int_t i = 0; i < fNpoints; i++) {
      if (exL)
         fExL[i] = exL[i];
      else
         fExL[i] = 0.;
      if (exH)
         fExH[i] = exH[i];
      else
         fExH[i] = 0.;

      for (Int_t j = 0; j < fNYErrors; j++) {
         if (Int_t(eyL.size()) > j && eyL[j].GetSize() > i)
            fEyL[j][i] = eyL[j][i];
         else
            fEyL[j][i] = 0.;
         if (Int_t(eyH.size()) > j && eyH[j].GetSize() > i)
            fEyH[j][i] = eyH[j][i];
         else
            fEyH[j][i] = 0.;
      }
   }

   CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors normal constructor with name, title, `np` points and `ne` y-errors.
///
/// If `exL`,`exH` are NULL, the corresponding values are preset to zero.
/// The multiple y-errors are passed as std::vectors of TArrayF objects.

TGraphMultiErrors::TGraphMultiErrors(const Char_t *name, const Char_t *title, Int_t np, Int_t ne, const Float_t *x,
                                     const Float_t *y, const Float_t *exL, const Float_t *exH, std::vector<TArrayF> eyL,
                                     std::vector<TArrayF> eyH, Int_t m)
   : TGraphMultiErrors(np, ne, x, y, exL, exH, eyL, eyH, m)
{
   SetNameTitle(name, title);
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors normal constructor with `np` points and `ne` y-errors.
///
/// If `exL`,`exH` are NULL, the corresponding values are preset to zero.
/// The multiple y-errors are passed as std::vectors of TArrayD objects.

TGraphMultiErrors::TGraphMultiErrors(Int_t np, Int_t ne, const Double_t *x, const Double_t *y, const Double_t *exL,
                                     const Double_t *exH, std::vector<TArrayD> eyL, std::vector<TArrayD> eyH, Int_t m)
   : TGraph(np, x, y), fNYErrors(ne), fSumErrorsMode(m)
{
   if (!CtorAllocate())
      return;

   Int_t n = fNpoints * sizeof(Double_t);

   if (exL)
      memcpy(fExL, exL, n);
   else
      memset(fExL, 0, n);
   if (exH)
      memcpy(fExH, exH, n);
   else
      memset(fExH, 0, n);

   for (Int_t i = 0; i < fNpoints; i++) {
      for (Int_t j = 0; j < fNYErrors; j++) {
         if (Int_t(eyL.size()) > j && eyL[j].GetSize() > i)
            fEyL[j][i] = eyL[j][i];
         else
            fEyL[j][i] = 0.;
         if (Int_t(eyH.size()) > j && eyH[j].GetSize() > i)
            fEyH[j][i] = eyH[j][i];
         else
            fEyH[j][i] = 0.;
      }
   }

   CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors normal constructor with `name`, `title`, `np` points and `ne` y-errors.
///
/// If `exL`,`exH` are NULL, the corresponding values are preset to zero.
/// The multiple y-errors are passed as std::vectors of TArrayD objects.

TGraphMultiErrors::TGraphMultiErrors(const Char_t *name, const Char_t *title, Int_t np, Int_t ne, const Double_t *x,
                                     const Double_t *y, const Double_t *exL, const Double_t *exH,
                                     std::vector<TArrayD> eyL, std::vector<TArrayD> eyH, Int_t m)
   : TGraphMultiErrors(np, ne, x, y, exL, exH, eyL, eyH, m)
{
   SetNameTitle(name, title);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with six vectors of floats in input and a single y error dimension.
/// The signature of this constructor is equal to the corresponding constructor of TGraphAsymmErrors.
/// A grapherrors is built with the X coordinates taken from `tvX` the Y coordinates from `tvY`
/// and the errors from vectors `tvExL`, `tvExH` and `tvEyL`, `tvEyH`.
/// The number of points in the graph is the minimum of number of points
/// in `tvX` and `tvY`.

TGraphMultiErrors::TGraphMultiErrors(const TVectorF &tvX, const TVectorF &tvY, const TVectorF &tvExL,
                                     const TVectorF &tvExH, const TVectorF &tvEyL, const TVectorF &tvEyH, Int_t m)
   : TGraph(), fNYErrors(1), fSumErrorsMode(m)
{
   fNpoints = TMath::Min(tvX.GetNrows(), tvY.GetNrows());

   if (!TGraph::CtorAllocate())
      return;

   if (!CtorAllocate())
      return;

   Int_t itvXL = tvX.GetLwb();
   Int_t itvYL = tvY.GetLwb();
   Int_t itvExLL = tvExL.GetLwb();
   Int_t itvExHL = tvExH.GetLwb();
   Int_t itvEyLL = tvEyL.GetLwb();
   Int_t itvEyHL = tvEyH.GetLwb();

   for (Int_t i = 0; i < fNpoints; i++) {
      fX[i] = tvX(itvXL + i);
      fY[i] = tvY(itvYL + i);
      fExL[i] = tvExL(itvExLL + i);
      fExH[i] = tvExH(itvExHL + i);
      fEyL[0][i] = tvEyL(itvEyLL + i);
      fEyH[0][i] = tvEyH(itvEyHL + i);
   }

   CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with six vectors of doubles in input and a single y error dimension.
/// The signature of this constructor is equal to the corresponding constructor of TGraphAsymmErrors.
/// A grapherrors is built with the X coordinates taken from `tvX` the Y coordinates from `tvY`
/// and the errors from vectors `tvExL`, `tvExH` and `tvEyL`, `tvEyH`.
/// The number of points in the graph is the minimum of number of points
/// in `tvX` and `tvY`.

TGraphMultiErrors::TGraphMultiErrors(const TVectorD &tvX, const TVectorD &tvY, const TVectorD &tvExL,
                                     const TVectorD &tvExH, const TVectorD &tvEyL, const TVectorD &tvEyH, Int_t m)
   : TGraph(), fNYErrors(1), fSumErrorsMode(m)
{
   fNpoints = TMath::Min(tvX.GetNrows(), tvY.GetNrows());

   if (!TGraph::CtorAllocate())
      return;

   if (!CtorAllocate())
      return;

   Int_t itvXL = tvX.GetLwb();
   Int_t itvYL = tvY.GetLwb();
   Int_t itvExLL = tvExL.GetLwb();
   Int_t itvExHL = tvExH.GetLwb();
   Int_t itvEyLL = tvEyL.GetLwb();
   Int_t itvEyHL = tvEyH.GetLwb();

   for (Int_t i = 0; i < fNpoints; i++) {
      fX[i] = tvX(i + itvXL);
      fY[i] = tvY(i + itvYL);
      fExL[i] = tvExL(i + itvExLL);
      fExH[i] = tvExH(i + itvExHL);
      fEyL[0][i] = tvEyL(i + itvEyLL);
      fEyH[0][i] = tvEyH(i + itvEyHL);
   }

   CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with multiple vectors of floats in input and multiple y error dimension.
/// A grapherrors is built with the X coordinates taken from `tvX` the Y coordinates from `tvY`
/// and the errors from vectors `tvExL`, `tvExH` and `tvEyL/H[yErrorDimension]`.
/// The number of points in the graph is the minimum of number of points
/// in `tvX` and `tvY`.

TGraphMultiErrors::TGraphMultiErrors(Int_t ne, const TVectorF &tvX, const TVectorF &tvY, const TVectorF &tvExL,
                                     const TVectorF &tvExH, const TVectorF *tvEyL, const TVectorF *tvEyH, Int_t m)
   : TGraph(), fNYErrors(ne), fSumErrorsMode(m)
{
   fNpoints = TMath::Min(tvX.GetNrows(), tvY.GetNrows());

   if (!TGraph::CtorAllocate())
      return;

   if (!CtorAllocate())
      return;

   Int_t itvXL = tvX.GetLwb();
   Int_t itvYL = tvY.GetLwb();
   Int_t itvExLL = tvExL.GetLwb();
   Int_t itvExHL = tvExH.GetLwb();

   for (Int_t i = 0; i < fNpoints; i++) {
      fX[i] = tvX(i + itvXL);
      fY[i] = tvY(i + itvYL);
      fExL[i] = tvExL(i + itvExLL);
      fExH[i] = tvExH(i + itvExHL);

      for (Int_t j = 0; j < ne; j++) {
         fEyL[j][i] = tvEyL[j](i + tvEyL[j].GetLwb());
         fEyH[j][i] = tvEyH[j](i + tvEyH[j].GetLwb());
      }
   }

   CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with multiple vectors of doubles in input and multiple y error dimensions
/// A grapherrors is built with the X coordinates taken from `tvX` the Y coordinates from `tvY`
/// and the errors from vectors `tvExL`, `tvExH` and `tvEyL/H[yErrorDimension]`.
/// The number of points in the graph is the minimum of number of points
/// in `tvX` and `tvY`.

TGraphMultiErrors::TGraphMultiErrors(Int_t ne, const TVectorD &tvX, const TVectorD &tvY, const TVectorD &tvExL,
                                     const TVectorD &tvExH, const TVectorD *tvEyL, const TVectorD *tvEyH, Int_t m)
   : TGraph(), fNYErrors(ne), fSumErrorsMode(m)
{
   fNpoints = TMath::Min(tvX.GetNrows(), tvY.GetNrows());

   if (!TGraph::CtorAllocate())
      return;

   if (!CtorAllocate())
      return;

   Int_t itvXL = tvX.GetLwb();
   Int_t itvYL = tvY.GetLwb();
   Int_t itvExLL = tvExL.GetLwb();
   Int_t itvExHL = tvExH.GetLwb();

   for (Int_t i = 0; i < fNpoints; i++) {
      fX[i] = tvX(i + itvXL);
      fY[i] = tvY(i + itvYL);
      fExL[i] = tvExL(i + itvExLL);
      fExH[i] = tvExH(i + itvExHL);

      for (Int_t j = 0; j < ne; j++) {
         fEyL[j][i] = tvEyL[j](i + tvEyL[j].GetLwb());
         fEyH[j][i] = tvEyH[j](i + tvEyH[j].GetLwb());
      }
   }

   CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors copy constructor.

TGraphMultiErrors::TGraphMultiErrors(const TGraphMultiErrors &tgme) : TGraph(tgme)
{
   fNYErrors = tgme.fNYErrors;
   fSumErrorsMode = tgme.fSumErrorsMode;

   if (!CtorAllocate())
      return;

   Int_t n = fNpoints * sizeof(Double_t);
   memcpy(fExL, tgme.fExL, n);
   memcpy(fExH, tgme.fExH, n);

   for (Int_t j = 0; j < fNYErrors; j++) {
      fEyL[j] = tgme.fEyL[j];
      fEyH[j] = tgme.fEyH[j];
      tgme.fAttFill[j].Copy(fAttFill[j]);
      tgme.fAttLine[j].Copy(fAttLine[j]);
   }

   CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors assignment operator.

TGraphMultiErrors &TGraphMultiErrors::operator=(const TGraphMultiErrors &tgme)
{
   if (this != &tgme) {
      TGraph::operator=(tgme);
      // delete arrays
      if (fExL)
         delete[] fExL;
      if (fExH)
         delete[] fExH;
      if (fEyLSum)
         delete[] fEyLSum;
      if (fEyHSum)
         delete[] fEyHSum;

      fNYErrors = tgme.fNYErrors;
      fSumErrorsMode = tgme.fSumErrorsMode;

      if (!CtorAllocate())
         return *this;

      Int_t n = fNpoints * sizeof(Double_t);
      memcpy(fExL, tgme.fExL, n);
      memcpy(fExH, tgme.fExH, n);
      memcpy(fEyLSum, tgme.fEyLSum, n);
      memcpy(fEyHSum, tgme.fEyHSum, n);

      for (Int_t j = 0; j < fNYErrors; j++) {
         fEyL[j] = tgme.fEyL[j];
         fEyH[j] = tgme.fEyH[j];
         tgme.fAttFill[j].Copy(fAttFill[j]);
         tgme.fAttLine[j].Copy(fAttLine[j]);
      }
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors constructor importing its parameters from the TH1 object passed as argument.
/// The low and high errors are set to the bin error of the histogram.

TGraphMultiErrors::TGraphMultiErrors(const TH1 *h, Int_t ne)
   : TGraph(h), fNYErrors(ne), fSumErrorsMode(TGraphMultiErrors::kOnlyFirst)
{
   if (!CtorAllocate())
      return;

   for (Int_t i = 0; i < fNpoints; i++) {
      fExL[i] = h->GetBinWidth(i + 1) * gStyle->GetErrorX();
      fExH[i] = h->GetBinWidth(i + 1) * gStyle->GetErrorX();
      fEyL[0][i] = h->GetBinError(i + 1);
      fEyH[0][i] = h->GetBinError(i + 1);

      for (Int_t j = 1; j < fNYErrors; j++) {
         fEyL[j][i] = 0.;
         fEyH[j][i] = 0.;
      }
   }

   CalcYErrorsSum();

   TAttFill::Copy(fAttFill[0]);
   TAttLine::Copy(fAttLine[0]);
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TGraphMultiErrors by dividing two input TH1 histograms:
/// pass/total. (see TGraphMultiErrors::Divide)

TGraphMultiErrors::TGraphMultiErrors(const TH1 *pass, const TH1 *total, Int_t ne, Option_t *option)
   : TGraph(pass ? pass->GetNbinsX() : 0), fNYErrors(ne), fSumErrorsMode(TGraphMultiErrors::kOnlyFirst)
{
   if (!pass || !total) {
      Error("TGraphMultiErrors", "Invalid histogram pointers");
      return;
   }

   if (!CtorAllocate())
      return;

   std::string sname = "divide_" + std::string(pass->GetName()) + "_by_" + std::string(total->GetName());
   SetName(sname.c_str());
   SetTitle(pass->GetTitle());

   // copy style from pass
   pass->TAttLine::Copy(*this);
   pass->TAttFill::Copy(*this);
   pass->TAttMarker::Copy(*this);

   Divide(pass, total, option);
   CalcYErrorsSum();

   TAttFill::Copy(fAttFill[0]);
   TAttLine::Copy(fAttLine[0]);
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphMultiErrors default destructor.

TGraphMultiErrors::~TGraphMultiErrors()
{
   if (fExL)
      delete[] fExL;
   if (fExH)
      delete[] fExH;
   fEyL.resize(0);
   fEyH.resize(0);
   if (fEyLSum)
      delete[] fEyLSum;
   if (fEyHSum)
      delete[] fEyHSum;
   fAttFill.resize(0);
   fAttLine.resize(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Should be called from ctors after `fNpoints` has been set
/// Note: This function should be called only from the constructor
/// since it does not delete previously existing arrays

Bool_t TGraphMultiErrors::CtorAllocate()
{
   if (!fNpoints || !fNYErrors) {
      fExL = fExH = nullptr;
      fEyL.resize(0);
      fEyH.resize(0);
      return kFALSE;
   }

   fExL = new Double_t[fMaxSize];
   fExH = new Double_t[fMaxSize];
   fEyL.resize(fNYErrors, TArrayD(fMaxSize));
   fEyH.resize(fNYErrors, TArrayD(fMaxSize));
   fEyLSum = new Double_t[fMaxSize];
   fEyHSum = new Double_t[fMaxSize];
   fAttFill.resize(fNYErrors);
   fAttLine.resize(fNYErrors);

   Int_t n = fMaxSize * sizeof(Double_t);
   memset(fExL, 0, n);
   memset(fExH, 0, n);
   memset(fEyLSum, 0, n);
   memset(fEyHSum, 0, n);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy and release.

void TGraphMultiErrors::CopyAndRelease(Double_t **newarrays, Int_t ibegin, Int_t iend, Int_t obegin)
{
   CopyPoints(newarrays, ibegin, iend, obegin);
   if (newarrays) {
      delete[] fX;
      fX = newarrays[0];
      delete[] fY;
      fY = newarrays[1];

      delete[] fExL;
      fExL = newarrays[2];
      delete[] fExH;
      fExH = newarrays[3];

      if (fEyLSum)
         delete[] fEyLSum;
      fEyLSum = newarrays[4];
      if (fEyHSum)
         delete[] fEyHSum;
      fEyHSum = newarrays[5];

      delete[] newarrays;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy errors from `fE***` to `arrays[***]`
/// or to `f***` Copy points.

Bool_t TGraphMultiErrors::CopyPoints(Double_t **arrays, Int_t ibegin, Int_t iend, Int_t obegin)
{
   if (TGraph::CopyPoints(arrays, ibegin, iend, obegin)) {
      Int_t n = (iend - ibegin) * sizeof(Double_t);

      if (arrays) {
         memmove(&arrays[2][obegin], &fExL[ibegin], n);
         memmove(&arrays[3][obegin], &fExH[ibegin], n);
         memmove(&arrays[4][obegin], &fEyLSum[ibegin], n);
         memmove(&arrays[5][obegin], &fEyHSum[ibegin], n);
      } else {
         memmove(&fExL[obegin], &fExL[ibegin], n);
         memmove(&fExH[obegin], &fExH[ibegin], n);
         memmove(&fEyLSum[obegin], &fEyLSum[ibegin], n);
         memmove(&fEyHSum[obegin], &fEyHSum[ibegin], n);
      }

      return kTRUE;
   } else
      return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set zero values for point arrays in the range `[begin, end]`.

void TGraphMultiErrors::FillZero(Int_t begin, Int_t end, Bool_t from_ctor)
{
   if (!from_ctor)
      TGraph::FillZero(begin, end, from_ctor);

   Int_t n = (end - begin) * sizeof(Double_t);
   memset(fExL + begin, 0, n);
   memset(fExH + begin, 0, n);
   memset(fEyLSum + begin, 0, n);
   memset(fEyHSum + begin, 0, n);

   for (Int_t j = 0; j < fNYErrors; j++) {
      memset(fEyL[j].GetArray() + begin, 0, n);
      memset(fEyH[j].GetArray() + begin, 0, n);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Recalculates the summed y error arrays.

void TGraphMultiErrors::CalcYErrorsSum() const
{
   if (!fEyLSum)
      fEyLSum = new Double_t[fNpoints];
   if (!fEyHSum)
      fEyHSum = new Double_t[fNpoints];

   for (Int_t i = 0; i < fNpoints; i++) {
      fEyLSum[i] = GetErrorYlow(i);
      fEyHSum[i] = GetErrorYhigh(i);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Protected function to perform the merge operation of a graph with multiple asymmetric errors.

Bool_t TGraphMultiErrors::DoMerge(const TGraph *tg)
{
   if (tg->GetN() == 0)
      return kFALSE;

   if (tg->IsA() == TGraphMultiErrors::Class()) {
      auto tgme = (TGraphMultiErrors *)tg;

      for (Int_t i = 0; i < tgme->GetN(); i++) {
         Int_t ipoint = GetN();
         Double_t x, y;
         tgme->GetPoint(i, x, y);
         SetPoint(ipoint, x, y);
         SetPointEX(ipoint, tgme->GetErrorXlow(i), tgme->GetErrorXhigh(i));
         for (Int_t j = 0; j < tgme->GetNYErrors(); j++)
            SetPointEY(ipoint, j, tgme->GetErrorYlow(i, j), tgme->GetErrorYhigh(i, j));
      }

      return kTRUE;
   } else {
      Warning("DoMerge", "Merging a %s is not compatible with a TGraphMultiErrors - Errors will be ignored",
              tg->IsA()->GetName());
      return TGraph::DoMerge(tg);
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Swap points.

void TGraphMultiErrors::SwapPoints(Int_t pos1, Int_t pos2)
{
   SwapValues(fExL, pos1, pos2);
   SwapValues(fExH, pos1, pos2);

   for (Int_t j = 0; j <= fNYErrors; j++) {
      SwapValues(fEyL[j].GetArray(), pos1, pos2);
      SwapValues(fEyH[j].GetArray(), pos1, pos2);
   }

   TGraph::SwapPoints(pos1, pos2);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new y error to the graph and fill it with the values from `eyL` and `eyH`

void TGraphMultiErrors::AddYError(Int_t np, const Double_t *eyL, const Double_t *eyH)
{
   fEyL.emplace_back(np, eyL);
   fEyH.emplace_back(np, eyH);
   fEyL.back().Set(fNpoints);
   fEyH.back().Set(fNpoints);
   fAttFill.emplace_back();
   fAttLine.emplace_back();

   fNYErrors += 1;

   CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// Allocate internal data structures for `size` points.
Double_t **TGraphMultiErrors::Allocate(Int_t size)
{
   return AllocateArrays(6, size);
}

////////////////////////////////////////////////////////////////////////////////
/// Apply a function to all data points \f$ y = f(x,y) \f$.
///
/// Errors are calculated as \f$ eyh = f(x,y+eyh)-f(x,y) \f$ and
/// \f$ eyl = f(x,y)-f(x,y-eyl) \f$
///
/// Only the first error dimension is affected.
///
/// Special treatment has to be applied for the functions where the
/// role of "up" and "down" is reversed.
///
/// Function suggested/implemented by Miroslav Helbich <helbich@mail.desy.de>

void TGraphMultiErrors::Apply(TF1 *f)
{
   Double_t x, y, eyL, eyH, eyLNew, eyHNew, fxy;

   if (fHistogram) {
      delete fHistogram;
      fHistogram = nullptr;
   }

   for (Int_t i = 0; i < fNpoints; i++) {
      GetPoint(i, x, y);
      eyL = GetErrorYlow(i, 0);
      eyH = GetErrorYhigh(i, 0);

      fxy = f->Eval(x, y);
      SetPoint(i, x, fxy);

      if (f->Eval(x, y - eyL) < f->Eval(x, y + eyH)) {
         eyLNew = TMath::Abs(fxy - f->Eval(x, y - eyL));
         eyHNew = TMath::Abs(f->Eval(x, y + eyH) - fxy);
      } else {
         eyHNew = TMath::Abs(fxy - f->Eval(x, y - eyL));
         eyLNew = TMath::Abs(f->Eval(x, y + eyH) - fxy);
      }

      // systematic errors and error on x doesn't change
      SetPointEY(i, 0, eyLNew, eyHNew);
   }

   if (gPad)
      gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// This function is only kept for backward compatibility.
/// You should rather use the Divide method.
/// It calls `Divide(pass,total,"cl=0.683 b(1,1) mode")` which is equivalent to the
/// former BayesDivide method.

void TGraphMultiErrors::BayesDivide(const TH1 *pass, const TH1 *total, Option_t *)
{
   Divide(pass, total, "cl=0.683 b(1,1) mode");
}

////////////////////////////////////////////////////////////////////////////////
/// This function was adapted from the TGraphAsymmErrors class.
/// See TGraphAsymmErrors::Divide for the documentation
///
/// Only the first error dimension is affected.

void TGraphMultiErrors::Divide(const TH1 *pass, const TH1 *total, Option_t *opt)
{
   // check pointers
   if (!pass || !total) {
      Error("Divide", "one of the passed pointers is zero");
      return;
   }

   // check dimension of histograms; only 1-dimensional ones are accepted
   if ((pass->GetDimension() > 1) || (total->GetDimension() > 1)) {
      Error("Divide", "passed histograms are not one-dimensional");
      return;
   }

   // check whether histograms are filled with weights -> use number of effective
   // entries
   Bool_t bEffective = false;
   // compare sum of weights with sum of squares of weights
   // re-compute here to be sure to get the right values
   Double_t psumw = 0;
   Double_t psumw2 = 0;
   if (pass->GetSumw2()->fN > 0) {
      for (int i = 0; i < pass->GetNbinsX(); ++i) {
         psumw += pass->GetBinContent(i);
         psumw2 += pass->GetSumw2()->At(i);
      }
   } else {
      psumw = pass->GetSumOfWeights();
      psumw2 = psumw;
   }
   if (TMath::Abs(psumw - psumw2) > 1e-6)
      bEffective = true;

   Double_t tsumw = 0;
   Double_t tsumw2 = 0;
   if (total->GetSumw2()->fN > 0) {
      for (int i = 0; i < total->GetNbinsX(); ++i) {
         tsumw += total->GetBinContent(i);
         tsumw2 += total->GetSumw2()->At(i);
      }
   } else {
      tsumw = total->GetSumOfWeights();
      tsumw2 = tsumw;
   }
   if (TMath::Abs(tsumw - tsumw2) > 1e-6)
      bEffective = true;

   // we do not want to ignore the weights
   // if (bEffective && (pass->GetSumw2()->fN == 0 || total->GetSumw2()->fN == 0) ) {
   //    Warning("Divide","histogram have been computed with weights but the sum of weight squares are not stored in the
   //    histogram. Error calculation is performed ignoring the weights"); bEffective = false;
   // }

   // parse option
   TString option = opt;
   option.ToLower();

   Bool_t bVerbose = false;
   // pointer to function returning the boundaries of the confidence interval
   //(is only used in the frequentist cases.)
   // Double_t (*pBound)(Int_t,Int_t,Double_t,Bool_t) = &TEfficiency::ClopperPearson; // default method
   Double_t (*pBound)(Double_t, Double_t, Double_t, Bool_t) = &TEfficiency::ClopperPearson; // default method
   // confidence level
   Double_t conf = 0.682689492137;
   // values for bayesian statistics
   Bool_t bIsBayesian = false;
   Double_t alpha = 1;
   Double_t beta = 1;

   // verbose mode
   if (option.Contains("v")) {
      option.ReplaceAll("v", "");
      bVerbose = true;
      if (bEffective)
         Info("Divide", "weight will be considered in the Histogram Ratio");
   }

   // confidence level
   if (option.Contains("cl=")) {
      Double_t level = -1;
      // coverity [secure_coding : FALSE]
      sscanf(strstr(option.Data(), "cl="), "cl=%lf", &level);
      if ((level > 0) && (level < 1))
         conf = level;
      else
         Warning("Divide", "given confidence level %.3lf is invalid", level);
      option.ReplaceAll("cl=", "");
   }

   // normal approximation
   if (option.Contains("n")) {
      option.ReplaceAll("n", "");
      pBound = &TEfficiency::Normal;
   }

   // clopper pearson interval
   if (option.Contains("cp")) {
      option.ReplaceAll("cp", "");
      pBound = &TEfficiency::ClopperPearson;
   }

   // wilson interval
   if (option.Contains("w")) {
      option.ReplaceAll("w", "");
      pBound = &TEfficiency::Wilson;
   }

   // agresti coull interval
   if (option.Contains("ac")) {
      option.ReplaceAll("ac", "");
      pBound = &TEfficiency::AgrestiCoull;
   }
   // Feldman-Cousins interval
   if (option.Contains("fc")) {
      option.ReplaceAll("fc", "");
      pBound = &TEfficiency::FeldmanCousins;
   }
   // mid-P Lancaster interval (In a later ROOT Version!)
   if (option.Contains("midp")) {
      option.ReplaceAll("midp", "");
      // pBound = &TEfficiency::MidPInterval;
   }

   // bayesian with prior
   if (option.Contains("b(")) {
      Double_t a = 0;
      Double_t b = 0;
      sscanf(strstr(option.Data(), "b("), "b(%lf,%lf)", &a, &b);
      if (a > 0)
         alpha = a;
      else
         Warning("Divide", "given shape parameter for alpha %.2lf is invalid", a);
      if (b > 0)
         beta = b;
      else
         Warning("Divide", "given shape parameter for beta %.2lf is invalid", b);
      option.ReplaceAll("b(", "");
      bIsBayesian = true;
   }

   // use posterior mode
   Bool_t usePosteriorMode = false;
   if (bIsBayesian && option.Contains("mode")) {
      usePosteriorMode = true;
      option.ReplaceAll("mode", "");
   }

   Bool_t plot0Bins = false;
   if (option.Contains("e0")) {
      plot0Bins = true;
      option.ReplaceAll("e0", "");
   }

   Bool_t useShortestInterval = false;
   if (bIsBayesian && (option.Contains("sh") || (usePosteriorMode && !option.Contains("cen")))) {
      useShortestInterval = true;
   }

   // interpret as Poisson ratio
   Bool_t bPoissonRatio = false;
   if (option.Contains("pois")) {
      bPoissonRatio = true;
      option.ReplaceAll("pois", "");
   }

   // weights works only in case of Normal approximation or Bayesian for binomial interval
   // in case of Poisson ratio we can use weights by rescaling the obtained results using the effective entries
   if ((bEffective && !bPoissonRatio) && !bIsBayesian && pBound != &TEfficiency::Normal) {
      Warning("Divide", "Histograms have weights: only Normal or Bayesian error calculation is supported");
      Info("Divide", "Using now the Normal approximation for weighted histograms");
   }

   if (bPoissonRatio) {
      if (pass->GetDimension() != total->GetDimension()) {
         Error("Divide", "passed histograms are not of the same dimension");
         return;
      }

      if (!TEfficiency::CheckBinning(*pass, *total)) {
         Error("Divide", "passed histograms are not consistent");
         return;
      }
   } else {
      // check consistency of histograms, allowing weights
      if (!TEfficiency::CheckConsistency(*pass, *total, "w")) {
         Error("Divide", "passed histograms are not consistent");
         return;
      }
   }

   // Set the graph to have a number of points equal to the number of histogram
   // bins
   Int_t nbins = pass->GetNbinsX();
   Set(nbins);

   // Ok, now set the points for each bin
   // (Note: the TH1 bin content is shifted to the right by one:
   //  bin=0 is underflow, bin=nbins+1 is overflow.)

   // this keeps track of the number of points added to the graph
   Int_t npoint = 0;
   // number of total and passed events
   Double_t t = 0, p = 0;
   Double_t tw = 0, tw2 = 0, pw = 0, pw2 = 0, wratio = 1; // for the case of weights
   // loop over all bins and fill the graph
   for (Int_t b = 1; b <= nbins; ++b) {
      // efficiency with lower and upper boundary of confidence interval default value when total =0;
      Double_t eff = 0., low = 0., upper = 0.;

      // special case in case of weights we have to consider the sum of weights and the sum of weight squares
      if (bEffective) {
         tw = total->GetBinContent(b);
         tw2 = (total->GetSumw2()->fN > 0) ? total->GetSumw2()->At(b) : tw;
         pw = pass->GetBinContent(b);
         pw2 = (pass->GetSumw2()->fN > 0) ? pass->GetSumw2()->At(b) : pw;

         if (bPoissonRatio) {
            // tw += pw;
            // tw2 += pw2;
            // compute ratio on the effective entries ( p and t)
            // special case is when (pw=0, pw2=0) in this case we cannot get the bin weight.
            // we use then the overall weight of the full histogram
            if (pw == 0 && pw2 == 0)
               p = 0;
            else
               p = (pw * pw) / pw2;

            if (tw == 0 && tw2 == 0)
               t = 0;
            else
               t = (tw * tw) / tw2;

            if (pw > 0 && tw > 0)
               // this is the ratio of the two bin weights ( pw/p  / t/tw )
               wratio = (pw * t) / (p * tw);
            else if (pw == 0 && tw > 0)
               // case p histogram has zero  compute the weights from all the histogram
               // weight of histogram - sumw2/sumw
               wratio = (psumw2 * t) / (psumw * tw);
            else if (tw == 0 && pw > 0)
               // case t histogram has zero  compute the weights from all the histogram
               // weight of histogram - sumw2/sumw
               wratio = (pw * tsumw) / (p * tsumw2);
            else if (p > 0)
               wratio = pw / p; // not sure if needed
            else {
               // case both pw and tw are zero - we skip these bins
               if (!plot0Bins)
                  continue; // skip bins with total <= 0
            }

            t += p;
            // std::cout << p << "   " << t << "  " << wratio << std::endl;
         } else if (tw <= 0 && !plot0Bins)
            continue; // skip bins with total <= 0

         // in the case of weights have the formula only for
         // the normal and  bayesian statistics (see below)

      }

      // use bin contents
      else {
         t = TMath::Nint(total->GetBinContent(b));
         p = TMath::Nint(pass->GetBinContent(b));

         if (bPoissonRatio)
            t += p;

         if (t == 0. && !plot0Bins)
            continue; // skip bins with total = 0
      }

      // using bayesian statistics
      if (bIsBayesian) {
         if ((bEffective && !bPoissonRatio) && tw2 <= 0) {
            // case of bins with zero errors
            eff = pw / tw;
            low = eff;
            upper = eff;
         } else {
            Double_t aa, bb;

            if (bEffective && !bPoissonRatio) {
               // tw/tw2 re-normalize the weights
               double norm = tw / tw2; // case of tw2 = 0 is treated above
               aa = pw * norm + alpha;
               bb = (tw - pw) * norm + beta;
            } else {
               aa = double(p) + alpha;
               bb = double(t - p) + beta;
            }
            if (usePosteriorMode)
               eff = TEfficiency::BetaMode(aa, bb);
            else
               eff = TEfficiency::BetaMean(aa, bb);

            if (useShortestInterval) {
               TEfficiency::BetaShortestInterval(conf, aa, bb, low, upper);
            } else {
               low = TEfficiency::BetaCentralInterval(conf, aa, bb, false);
               upper = TEfficiency::BetaCentralInterval(conf, aa, bb, true);
            }
         }
      }
      // case of non-bayesian statistics
      else {
         if (bEffective && !bPoissonRatio) {

            if (tw > 0) {

               eff = pw / tw;

               // use normal error calculation using variance of MLE with weights (F.James 8.5.2)
               // this is the same formula used in ROOT for TH1::Divide("B")

               double variance = (pw2 * (1. - 2 * eff) + tw2 * eff * eff) / (tw * tw);
               double sigma = sqrt(variance);

               double prob = 0.5 * (1. - conf);
               double delta = ROOT::Math::normal_quantile_c(prob, sigma);
               low = eff - delta;
               upper = eff + delta;
               if (low < 0)
                  low = 0;
               if (upper > 1)
                  upper = 1.;
            }
         } else {
            // when not using weights (all cases) or in case of  Poisson ratio with weights
            if (t != 0.)
               eff = ((Double_t)p) / t;

            low = pBound(t, p, conf, false);
            upper = pBound(t, p, conf, true);
         }
      }
      // treat as Poisson ratio
      if (bPoissonRatio) {
         Double_t ratio = eff / (1 - eff);
         // take the intervals in eff as intervals in the Poisson ratio
         low = low / (1. - low);
         upper = upper / (1. - upper);
         eff = ratio;
         if (bEffective) {
            // scale result by the ratio of the weight
            eff *= wratio;
            low *= wratio;
            upper *= wratio;
         }
      }
      // Set the point center and its errors
      if (TMath::Finite(eff)) {
         SetPoint(npoint, pass->GetBinCenter(b), eff);
         SetPointEX(npoint, pass->GetBinCenter(b) - pass->GetBinLowEdge(b),
                    pass->GetBinLowEdge(b) - pass->GetBinCenter(b) + pass->GetBinWidth(b));
         SetPointEY(npoint, 0, eff - low, upper - eff);
         npoint++; // we have added a point to the graph
      }
   }

   Set(npoint); // tell the graph how many points we've really added
   if (npoint < nbins)
      Warning("Divide", "Number of graph points is different than histogram bins - %d points have been skipped",
              nbins - npoint);

   if (bVerbose) {
      Info("Divide", "made a graph with %d points from %d bins", npoint, nbins);
      Info("Divide", "used confidence level: %.2lf\n", conf);
      if (bIsBayesian)
         Info("Divide", "used prior probability ~ beta(%.2lf,%.2lf)", alpha, beta);
      Print();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute Range.

void TGraphMultiErrors::ComputeRange(Double_t &xmin, Double_t &ymin, Double_t &xmax, Double_t &ymax) const
{
   TGraph::ComputeRange(xmin, ymin, xmax, ymax);

   for (Int_t i = 0; i < fNpoints; i++) {
      if (fX[i] - fExL[i] < xmin) {
         if (gPad && gPad->GetLogx()) {
            if (fExL[i] < fX[i])
               xmin = fX[i] - fExL[i];
            else
               xmin = TMath::Min(xmin, fX[i] / 3.);
         } else
            xmin = fX[i] - fExL[i];
      }

      if (fX[i] + fExH[i] > xmax)
         xmax = fX[i] + fExH[i];

      Double_t eyLMax = 0., eyHMax = 0.;
      for (Int_t j = 0; j < fNYErrors; j++) {
         eyLMax = TMath::Max(eyLMax, fEyL[j][i]);
         eyHMax = TMath::Max(eyHMax, fEyH[j][i]);
      }

      if (fY[i] - eyLMax < ymin) {
         if (gPad && gPad->GetLogy()) {
            if (eyLMax < fY[i])
               ymin = fY[i] - eyLMax;
            else
               ymin = TMath::Min(ymin, fY[i] / 3.);
         } else
            ymin = fY[i] - eyLMax;
      }

      if (fY[i] + eyHMax > ymax)
         ymax = fY[i] + eyHMax;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Deletes the y error with the index `e`.
/// Note that you must keep at least 1 error

void TGraphMultiErrors::DeleteYError(Int_t e)
{
   if (fNYErrors == 1 || e >= fNYErrors)
      return;

   fEyL.erase(fEyL.begin() + e);
   fEyH.erase(fEyH.begin() + e);
   fAttFill.erase(fAttFill.begin() + e);
   fAttLine.erase(fAttLine.begin() + e);

   fNYErrors -= 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get error on x coordinate for point `i`.
/// In case of asymmetric errors the mean of the square sum is returned

Double_t TGraphMultiErrors::GetErrorX(Int_t i) const
{
   if (i < 0 || i >= fNpoints || (!fExL && !fExH))
      return -1.;

   Double_t exL = fExL ? fExL[i] : 0.;
   Double_t exH = fExH ? fExH[i] : 0.;
   return TMath::Sqrt((exL * exL + exH * exH) / 2.);
}

////////////////////////////////////////////////////////////////////////////////
/// Get error on y coordinate for point `i`.
/// The multiple errors of the dimensions are summed according to `fSumErrorsMode`.
/// In case of asymmetric errors the mean of the square sum is returned

Double_t TGraphMultiErrors::GetErrorY(Int_t i) const
{
   if (i < 0 || i >= fNpoints || (fEyL.empty() && fEyH.empty()))
      return -1.;

   Double_t eyL = GetErrorYlow(i);
   Double_t eyH = GetErrorYhigh(i);
   return TMath::Sqrt((eyL * eyL + eyH * eyH) / 2.);
}

////////////////////////////////////////////////////////////////////////////////
/// Get error e on y coordinate for point `i`.
/// In case of asymmetric errors the mean of the square sum is returned

Double_t TGraphMultiErrors::GetErrorY(Int_t i, Int_t e) const
{
   if (i < 0 || i >= fNpoints || e >= fNYErrors || (fEyL.empty() && fEyH.empty()))
      return -1.;

   Double_t eyL = fEyL.empty() ? 0. : fEyL[e][i];
   Double_t eyH = fEyH.empty() ? 0. : fEyH[e][i];
   return TMath::Sqrt((eyL * eyL + eyH * eyH) / 2.);
}

////////////////////////////////////////////////////////////////////////////////
/// Get low error on x coordinate for point `i`.

Double_t TGraphMultiErrors::GetErrorXlow(Int_t i) const
{
   if (i < 0 || i >= fNpoints || !fExL)
      return -1.;
   else
      return fExL[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Get high error on x coordinate for point `i`.

Double_t TGraphMultiErrors::GetErrorXhigh(Int_t i) const
{
   if (i < 0 || i >= fNpoints || !fExH)
      return -1.;
   else
      return fExH[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Get low error on y coordinate for point `i`.
/// The multiple errors of the dimensions are summed according to `fSumErrorsMode`.

Double_t TGraphMultiErrors::GetErrorYlow(Int_t i) const
{
   if (i < 0 || i >= fNpoints || fEyL.empty())
      return -1.;

   if (fSumErrorsMode == TGraphMultiErrors::kOnlyFirst)
      return fEyL[0][i];
   else if (fSumErrorsMode == TGraphMultiErrors::kSquareSum) {
      Double_t sum = 0.;
      for (Int_t j = 0; j < fNYErrors; j++)
         sum += fEyL[j][i] * fEyL[j][i];
      return TMath::Sqrt(sum);
   } else if (fSumErrorsMode == TGraphMultiErrors::kAbsSum) {
      Double_t sum = 0.;
      for (Int_t j = 0; j < fNYErrors; j++)
         sum += fEyL[j][i];
      return sum;
   }

   return -1.;
}

////////////////////////////////////////////////////////////////////////////////
/// Get high error on y coordinate for point `i`.
/// The multiple errors of the dimensions are summed according to `fSumErrorsMode`.

Double_t TGraphMultiErrors::GetErrorYhigh(Int_t i) const
{
   if (i < 0 || i >= fNpoints || fEyH.empty())
      return -1.;

   if (fSumErrorsMode == TGraphMultiErrors::kOnlyFirst)
      return fEyH[0][i];
   else if (fSumErrorsMode == TGraphMultiErrors::kSquareSum) {
      Double_t sum = 0.;
      for (Int_t j = 0; j < fNYErrors; j++)
         sum += fEyH[j][i] * fEyH[j][i];
      return TMath::Sqrt(sum);
   } else if (fSumErrorsMode == TGraphMultiErrors::kAbsSum) {
      Double_t sum = 0.;
      for (Int_t j = 0; j < fNYErrors; j++)
         sum += fEyH[j][i];
      return sum;
   }

   return -1.;
}

////////////////////////////////////////////////////////////////////////////////
/// Get low error e on y coordinate for point `i`.

Double_t TGraphMultiErrors::GetErrorYlow(Int_t i, Int_t e) const
{
   if (i < 0 || i >= fNpoints || e >= fNYErrors || fEyL.empty())
      return -1.;

   return fEyL[e][i];
}

////////////////////////////////////////////////////////////////////////////////
/// Get high error e on y coordinate for point `i`.

Double_t TGraphMultiErrors::GetErrorYhigh(Int_t i, Int_t e) const
{
   if (i < 0 || i >= fNpoints || e >= fNYErrors || fEyH.empty())
      return -1.;

   return fEyH[e][i];
}

////////////////////////////////////////////////////////////////////////////////
/// Get all low errors on y coordinates as an array summed according to `fSumErrorsMode`.

Double_t *TGraphMultiErrors::GetEYlow() const
{
   if (!fEyLSum)
      CalcYErrorsSum();

   return fEyLSum;
}

////////////////////////////////////////////////////////////////////////////////
/// Get all high errors on y coordinates as an array summed according to `fSumErrorsMode`.

Double_t *TGraphMultiErrors::GetEYhigh() const
{
   if (!fEyHSum)
      CalcYErrorsSum();

   return fEyHSum;
}

////////////////////////////////////////////////////////////////////////////////
/// Get all low errors `e` on y coordinates as an array.

Double_t *TGraphMultiErrors::GetEYlow(Int_t e)
{
   if (e >= fNYErrors || fEyL.empty())
      return nullptr;
   else
      return fEyL[e].GetArray();
}

////////////////////////////////////////////////////////////////////////////////
/// Get all high errors `e` on y coordinates as an array.

Double_t *TGraphMultiErrors::GetEYhigh(Int_t e)
{
   if (e >= fNYErrors || fEyH.empty())
      return nullptr;
   else
      return fEyH[e].GetArray();
}

////////////////////////////////////////////////////////////////////////////////
/// Get AttFill pointer for specified error dimension.

TAttFill *TGraphMultiErrors::GetAttFill(Int_t e)
{
   if (e >= 0 && e < fNYErrors)
      return &fAttFill.at(e);
   else
      return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Get AttLine pointer for specified error dimension.

TAttLine *TGraphMultiErrors::GetAttLine(Int_t e)
{
   if (e >= 0 && e < fNYErrors)
      return &fAttLine.at(e);
   else
      return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Get Fill Color for specified error e (-1 = Global and x errors).

Color_t TGraphMultiErrors::GetFillColor(Int_t e) const
{
   if (e == -1)
      return GetFillColor();
   else if (e >= 0 && e < fNYErrors)
      return fAttFill[e].GetFillColor();
   else
      return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get Fill Style for specified error e (-1 = Global and x errors).

Style_t TGraphMultiErrors::GetFillStyle(Int_t e) const
{
   if (e == -1)
      return GetFillStyle();
   else if (e >= 0 && e < fNYErrors)
      return fAttFill[e].GetFillStyle();
   else
      return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get Line Color for specified error e (-1 = Global and x errors).

Color_t TGraphMultiErrors::GetLineColor(Int_t e) const
{
   if (e == -1)
      return GetLineColor();
   else if (e >= 0 && e < fNYErrors)
      return fAttLine[e].GetLineColor();
   else
      return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get Line Style for specified error e (-1 = Global and x errors).

Style_t TGraphMultiErrors::GetLineStyle(Int_t e) const
{
   if (e == -1)
      return GetLineStyle();
   else if (e >= 0 && e < fNYErrors)
      return fAttLine[e].GetLineStyle();
   else
      return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get Line Width for specified error e (-1 = Global and x errors).

Width_t TGraphMultiErrors::GetLineWidth(Int_t e) const
{
   if (e == -1)
      return GetLineWidth();
   else if (e >= 0 && e < fNYErrors)
      return fAttLine[e].GetLineWidth();
   else
      return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Print graph and errors values.

void TGraphMultiErrors::Print(Option_t *) const
{
   for (Int_t i = 0; i < fNpoints; i++) {
      printf("x[%d]=%g, y[%d]=%g", i, fX[i], i, fY[i]);
      if (fExL)
         printf(", exl[%d]=%g", i, fExL[i]);
      if (fExH)
         printf(", exh[%d]=%g", i, fExH[i]);
      if (!fEyL.empty())
         for (Int_t j = 0; j < fNYErrors; j++)
            printf(", eyl[%d][%d]=%g", j, i, fEyL[j][i]);
      if (!fEyH.empty())
         for (Int_t j = 0; j < fNYErrors; j++)
            printf(", eyh[%d][%d]=%g", j, i, fEyH[j][i]);
      printf("\n");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TGraphMultiErrors::SavePrimitive(std::ostream &out, Option_t *option)
{
   char quote = '"';
   out << "   " << std::endl;

   if (gROOT->ClassSaved(TGraphMultiErrors::Class()))
      out << "   ";
   else
      out << "   TGraphMultiErrors* ";

   out << "tgme = new TGraphMultiErrors(" << fNpoints << ", " << fNYErrors << ");" << std::endl;
   out << "   tgme->SetName(" << quote << GetName() << quote << ");" << std::endl;
   out << "   tgme->SetTitle(" << quote << GetTitle() << quote << ");" << std::endl;

   SaveFillAttributes(out, "tgme", 0, 1001);
   SaveLineAttributes(out, "tgme", 1, 1, 1);
   SaveMarkerAttributes(out, "tgme", 1, 1, 1);

   for (Int_t j = 0; j < fNYErrors; j++) {
      fAttFill[j].SaveFillAttributes(out, Form("tgme->GetAttFill(%d)", j), 0, 1001);
      fAttLine[j].SaveLineAttributes(out, Form("tgme->GetAttLine(%d)", j), 1, 1, 1);
   }

   for (Int_t i = 0; i < fNpoints; i++) {
      out << "   tgme->SetPoint(" << i << ", " << fX[i] << ", " << fY[i] << ");" << std::endl;
      out << "   tgme->SetPointEX(" << i << ", " << fExL[i] << ", " << fExH[i] << ");" << std::endl;

      for (Int_t j = 0; j < fNYErrors; j++)
         out << "   tgme->SetPointEY(" << i << ", " << j << ", " << fEyL[j][i] << ", " << fEyH[j][i] << ");"
             << std::endl;
   }

   static Int_t frameNumber = 0;
   if (fHistogram) {
      frameNumber++;
      TString hname = fHistogram->GetName();
      hname += frameNumber;
      fHistogram->SetName(Form("Graph_%s", hname.Data()));
      fHistogram->SavePrimitive(out, "nodraw");
      out << "   tgme->SetHistogram(" << fHistogram->GetName() << ");" << std::endl;
      out << "   " << std::endl;
   }

   // save list of functions
   TIter next(fFunctions);
   TObject *obj;
   while ((obj = next())) {
      obj->SavePrimitive(out, "nodraw");
      if (obj->InheritsFrom("TPaveStats")) {
         out << "   tgme->GetListOfFunctions()->Add(ptstats);" << std::endl;
         out << "   ptstats->SetParent(tgme->GetListOfFunctions());" << std::endl;
      } else
         out << "   tgme->GetListOfFunctions()->Add(" << obj->GetName() << ");" << std::endl;
   }

   const char *l = strstr(option, "multigraph");
   if (l)
      out << "   multigraph->Add(tgme, " << quote << l + 10 << quote << ");" << std::endl;
   else
      out << "   tgme->Draw(" << quote << option << quote << ");" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply the values and errors of a TGraphMultiErrors by a constant c1.
///
/// If option contains "x" the x values and errors are scaled
/// If option contains "y" the y values and (multiple) errors are scaled
/// If option contains "xy" both x and y values and (multiple) errors are scaled

void TGraphMultiErrors::Scale(Double_t c1, Option_t *option)
{
   TGraph::Scale(c1, option);
   TString opt = option; opt.ToLower();
   if (opt.Contains("x") && GetEXlow()) {
      for (Int_t i=0; i<GetN(); i++)
         GetEXlow()[i] *= c1;
   }
   if (opt.Contains("x") && GetEXhigh()) {
      for (Int_t i=0; i<GetN(); i++)
         GetEXhigh()[i] *= c1;
   }
   if (opt.Contains("y")) {
      for (size_t d=0; d<fEyL.size(); d++)
         for (Int_t i=0; i<fEyL[d].GetSize(); i++)
            fEyL[d][i] *= c1;
      for (size_t d=0; d<fEyH.size(); d++)
         for (Int_t i=0; i<fEyH[d].GetSize(); i++)
            fEyH[d][i] *= c1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set ex and ey values for point pointed by the mouse.
///
/// Up to 3 y error dimensions possible.

void TGraphMultiErrors::SetPointError(Double_t exL, Double_t exH, Double_t eyL1, Double_t eyH1, Double_t eyL2,
                                      Double_t eyH2, Double_t eyL3, Double_t eyH3)
{
   Int_t px = gPad->GetEventX();
   Int_t py = gPad->GetEventY();

   // localize point to be deleted
   Int_t ipoint = -2;
   Int_t i;
   // start with a small window (in case the mouse is very close to one point)
   for (i = 0; i < fNpoints; i++) {
      Int_t dpx = px - gPad->XtoAbsPixel(gPad->XtoPad(fX[i]));
      Int_t dpy = py - gPad->YtoAbsPixel(gPad->YtoPad(fY[i]));

      if (dpx * dpx + dpy * dpy < 25) {
         ipoint = i;
         break;
      }
   }

   if (ipoint == -2)
      return;

   SetPointEX(ipoint, exL, exH);

   if (fNYErrors > 0)
      SetPointEY(ipoint, 0, eyL1, eyH1);
   if (fNYErrors > 1)
      SetPointEY(ipoint, 1, eyL2, eyH2);
   if (fNYErrors > 2)
      SetPointEY(ipoint, 2, eyL3, eyH3);
   gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Set ex and ey values for point `i`.

void TGraphMultiErrors::SetPointError(Int_t i, Int_t ne, Double_t exL, Double_t exH, const Double_t *eyL,
                                      const Double_t *eyH)
{
   SetPointEX(i, exL, exH);
   SetPointEY(i, ne, eyL, eyH);
}

////////////////////////////////////////////////////////////////////////////////
/// Set ex values for point `i`.

void TGraphMultiErrors::SetPointEX(Int_t i, Double_t exL, Double_t exH)
{
   SetPointEXlow(i, exL);
   SetPointEXhigh(i, exH);
}

////////////////////////////////////////////////////////////////////////////////
/// Set exL value for point `i`.

void TGraphMultiErrors::SetPointEXlow(Int_t i, Double_t exL)
{
   if (i < 0)
      return;

   if (i >= fNpoints) {
      // re-allocate the object
      TGraphMultiErrors::SetPoint(i, 0., 0.);
   }

   fExL[i] = exL;
}

////////////////////////////////////////////////////////////////////////////////
/// Set exH value for point `i`.

void TGraphMultiErrors::SetPointEXhigh(Int_t i, Double_t exH)
{
   if (i < 0)
      return;

   if (i >= fNpoints) {
      // re-allocate the object
      TGraphMultiErrors::SetPoint(i, 0., 0.);
   }

   fExH[i] = exH;
}

////////////////////////////////////////////////////////////////////////////////
/// Set ey values for point `i`.

void TGraphMultiErrors::SetPointEY(Int_t i, Int_t ne, const Double_t *eyL, const Double_t *eyH)
{
   SetPointEYlow(i, ne, eyL);
   SetPointEYhigh(i, ne, eyH);
}

////////////////////////////////////////////////////////////////////////////////
/// Set eyL values for point `i`.

void TGraphMultiErrors::SetPointEYlow(Int_t i, Int_t ne, const Double_t *eyL)
{
   for (Int_t j = 0; j < fNYErrors; j++) {
      if (j < ne)
         SetPointEYlow(i, j, eyL[j]);
      else
         SetPointEYlow(i, j, 0.);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set eyH values for point `i`.

void TGraphMultiErrors::SetPointEYhigh(Int_t i, Int_t ne, const Double_t *eyH)
{
   for (Int_t j = 0; j < fNYErrors; j++) {
      if (j < ne)
         SetPointEYhigh(i, j, eyH[j]);
      else
         SetPointEYhigh(i, j, 0.);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set error e ey values for point `i`.

void TGraphMultiErrors::SetPointEY(Int_t i, Int_t e, Double_t eyL, Double_t eyH)
{
   SetPointEYlow(i, e, eyL);
   SetPointEYhigh(i, e, eyH);
}

////////////////////////////////////////////////////////////////////////////////
/// Set error e eyL value for point `i`.

void TGraphMultiErrors::SetPointEYlow(Int_t i, Int_t e, Double_t eyL)
{
   if (i < 0 || e < 0)
      return;

   if (i >= fNpoints)
      // re-allocate the object
      TGraphMultiErrors::SetPoint(i, 0., 0.);

   while (e >= fNYErrors)
      AddYError(fNpoints);

   fEyL[e][i] = eyL;
   if (fEyLSum)
      fEyLSum[i] = GetErrorYlow(i);
   else
      CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// Set error e eyH value for point `i`.

void TGraphMultiErrors::SetPointEYhigh(Int_t i, Int_t e, Double_t eyH)
{
   if (i < 0 || e < 0)
      return;

   if (i >= fNpoints)
      // re-allocate the object
      TGraphMultiErrors::SetPoint(i, 0., 0.);

   while (e >= fNYErrors)
      AddYError(fNpoints);

   fEyH[e][i] = eyH;
   if (fEyHSum)
      fEyHSum[i] = GetErrorYhigh(i);
   else
      CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// Set error e ey values.

void TGraphMultiErrors::SetEY(Int_t e, Int_t np, const Double_t *eyL, const Double_t *eyH)
{
   SetEYlow(e, np, eyL);
   SetEYhigh(e, np, eyH);
}

////////////////////////////////////////////////////////////////////////////////
/// Set error e eyL values.

void TGraphMultiErrors::SetEYlow(Int_t e, Int_t np, const Double_t *eyL)
{
   for (Int_t i = 0; i < fNpoints; i++) {
      if (i < np)
         SetPointEYlow(i, e, eyL[i]);
      else
         SetPointEYlow(i, e, 0.);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set error e eyH values.

void TGraphMultiErrors::SetEYhigh(Int_t e, Int_t np, const Double_t *eyH)
{
   for (Int_t i = 0; i < fNpoints; i++) {
      if (i < np)
         SetPointEYhigh(i, e, eyH[i]);
      else
         SetPointEYhigh(i, e, 0.);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the sum errors mode and recalculate summed errors.
void TGraphMultiErrors::SetSumErrorsMode(Int_t m)
{
   if (fSumErrorsMode == m)
      return;
   fSumErrorsMode = m;
   CalcYErrorsSum();
}

////////////////////////////////////////////////////////////////////////////////
/// Set TAttFill parameters of error e by copying from another TAttFill (-1 = Global and x errors).

void TGraphMultiErrors::SetAttFill(Int_t e, TAttFill *taf)
{
   if (e == -1)
      taf->TAttFill::Copy(*this);
   else if (e >= 0 && e < fNYErrors)
      taf->TAttFill::Copy(fAttFill[e]);
}

////////////////////////////////////////////////////////////////////////////////
/// Set TAttLine parameters of error e by copying from another TAttLine (-1 = Global and x errors).

void TGraphMultiErrors::SetAttLine(Int_t e, TAttLine *taf)
{
   if (e == -1)
      taf->TAttLine::Copy(*this);
   else if (e >= 0 && e < fNYErrors)
      taf->TAttLine::Copy(fAttLine[e]);
}

////////////////////////////////////////////////////////////////////////////////
/// Set Fill Color of error e (-1 = Global and x errors).

void TGraphMultiErrors::SetFillColor(Int_t e, Color_t fcolor)
{
   if (e == -1)
      SetFillColor(fcolor);
   else if (e >= 0 && e < fNYErrors)
      fAttFill[e].SetFillColor(fcolor);
}

////////////////////////////////////////////////////////////////////////////////
/// Set Fill Color and Alpha of error e (-1 = Global and x errors).

void TGraphMultiErrors::SetFillColorAlpha(Int_t e, Color_t fcolor, Float_t falpha)
{
   if (e == -1)
      SetFillColorAlpha(fcolor, falpha);
   else if (e >= 0 && e < fNYErrors)
      fAttFill[e].SetFillColorAlpha(fcolor, falpha);
}

////////////////////////////////////////////////////////////////////////////////
/// Set Fill Style of error e (-1 = Global and x errors).

void TGraphMultiErrors::SetFillStyle(Int_t e, Style_t fstyle)
{
   if (e == -1)
      SetFillStyle(fstyle);
   else if (e >= 0 && e < fNYErrors)
      fAttFill[e].SetFillStyle(fstyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Set Line Color of error e (-1 = Global and x errors).

void TGraphMultiErrors::SetLineColor(Int_t e, Color_t lcolor)
{
   if (e == -1)
      SetLineColor(lcolor);
   else if (e >= 0 && e < fNYErrors)
      fAttLine[e].SetLineColor(lcolor);
}

////////////////////////////////////////////////////////////////////////////////
/// Set Line Color and Alpha of error e (-1 = Global and x errors).

void TGraphMultiErrors::SetLineColorAlpha(Int_t e, Color_t lcolor, Float_t lalpha)
{
   if (e == -1)
      SetLineColorAlpha(lcolor, lalpha);
   else if (e >= 0 && e < fNYErrors)
      fAttLine[e].SetLineColorAlpha(lcolor, lalpha);
}

////////////////////////////////////////////////////////////////////////////////
/// Set Line Style of error e (-1 = Global and x errors).

void TGraphMultiErrors::SetLineStyle(Int_t e, Style_t lstyle)
{
   if (e == -1)
      SetLineStyle(lstyle);
   else if (e >= 0 && e < fNYErrors)
      fAttLine[e].SetLineStyle(lstyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Set Line Width of error e (-1 = Global and x errors).

void TGraphMultiErrors::SetLineWidth(Int_t e, Width_t lwidth)
{
   if (e == -1)
      SetLineWidth(lwidth);
   else if (e >= 0 && e < fNYErrors)
      fAttLine[e].SetLineWidth(lwidth);
}
