// @(#)root/hist:$Id$
// Author: Dave Morrison  30/06/2003

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <cstring>
#include <iostream>

#include "TROOT.h"
#include "TGraphBentErrors.h"
#include "TMath.h"
#include "TVirtualPad.h"
#include "TH1.h"
#include "TF1.h"

ClassImp(TGraphBentErrors);


////////////////////////////////////////////////////////////////////////////////

/** \class  TGraphBentErrors
    \ingroup Hist
A TGraphBentErrors is a TGraph with bent, asymmetric error bars.

The TGraphBentErrors painting is performed thanks to the TGraphPainter
class. All details about the various painting options are given in this class.

The picture below gives an example:
Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","A Simple Graph with bent error bars",200,10,700,500);
   const Int_t n = 10;
   Double_t x[n]    = {-0.22, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
   Double_t y[n]    = {1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
   Double_t exl[n]  = {.05,.1,.07,.07,.04,.05,.06,.07,.08,.05};
   Double_t eyl[n]  = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};
   Double_t exh[n]  = {.02,.08,.05,.05,.03,.03,.04,.05,.06,.03};
   Double_t eyh[n]  = {.6,.5,.4,.3,.2,.2,.3,.4,.5,.6};
   Double_t exld[n] = {.0,.0,.0,.0,.0,.0,.0,.0,.0,.0};
   Double_t eyld[n] = {.0,.0,.05,.0,.0,.0,.0,.0,.0,.0};
   Double_t exhd[n] = {.0,.0,.0,.0,.0,.0,.0,.0,.0,.0};
   Double_t eyhd[n] = {.0,.0,.0,.0,.0,.0,.0,.0,.05,.0};
   auto gr = new TGraphBentErrors(n,x,y,exl,exh,eyl,eyh,exld,exhd,eyld,eyhd);
   gr->SetTitle("TGraphBentErrors Example");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("ALP");
}
End_Macro
*/


////////////////////////////////////////////////////////////////////////////////
/// TGraphBentErrors default constructor.

TGraphBentErrors::TGraphBentErrors(): TGraph()
{
   if (!CtorAllocate()) return;
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphBentErrors copy constructor

TGraphBentErrors::TGraphBentErrors(const TGraphBentErrors &gr)
       : TGraph(gr)
{
   if (!CtorAllocate()) return;
   Int_t n = fNpoints*sizeof(Double_t);
   memcpy(fEXlow, gr.fEXlow, n);
   memcpy(fEYlow, gr.fEYlow, n);
   memcpy(fEXhigh, gr.fEXhigh, n);
   memcpy(fEYhigh, gr.fEYhigh, n);
   memcpy(fEXlowd, gr.fEXlowd, n);
   memcpy(fEYlowd, gr.fEYlowd, n);
   memcpy(fEXhighd, gr.fEXhighd, n);
   memcpy(fEYhighd, gr.fEYhighd, n);
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphBentErrors normal constructor.
///
///  the arrays are preset to zero

TGraphBentErrors::TGraphBentErrors(Int_t n)
       : TGraph(n)
{
   if (!CtorAllocate()) return;
   FillZero(0, fNpoints);
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphBentErrors normal constructor.
///
/// if exl,h or eyl,h are null, the corresponding arrays are preset to zero

TGraphBentErrors::TGraphBentErrors(Int_t n,
                                   const Float_t *x, const Float_t *y,
                                   const Float_t *exl, const Float_t *exh,
                                   const Float_t *eyl, const Float_t *eyh,
                                   const Float_t *exld, const Float_t *exhd,
                                   const Float_t *eyld, const Float_t *eyhd)
  : TGraph(n,x,y)
{
   if (!CtorAllocate()) return;

   for (Int_t i=0;i<n;i++) {
      if (exl) fEXlow[i]  = exl[i];
      else     fEXlow[i]  = 0;
      if (exh) fEXhigh[i] = exh[i];
      else     fEXhigh[i] = 0;
      if (eyl) fEYlow[i]  = eyl[i];
      else     fEYlow[i]  = 0;
      if (eyh) fEYhigh[i] = eyh[i];
      else     fEYhigh[i] = 0;

      if (exld) fEXlowd[i]  = exld[i];
      else     fEXlowd[i]  = 0;
      if (exhd) fEXhighd[i] = exhd[i];
      else     fEXhighd[i] = 0;
      if (eyld) fEYlowd[i]  = eyld[i];
      else     fEYlowd[i]  = 0;
      if (eyhd) fEYhighd[i] = eyhd[i];
      else     fEYhighd[i] = 0;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphBentErrors normal constructor.
///
/// if exl,h or eyl,h are null, the corresponding arrays are preset to zero

TGraphBentErrors::TGraphBentErrors(Int_t n,
                                   const Double_t *x, const Double_t *y,
                                   const Double_t *exl, const Double_t *exh,
                                   const Double_t *eyl, const Double_t *eyh,
                                   const Double_t *exld, const Double_t *exhd,
                                   const Double_t *eyld, const Double_t *eyhd)
  : TGraph(n,x,y)
{
   if (!CtorAllocate()) return;
   n = sizeof(Double_t)*fNpoints;

      if (exl) memcpy(fEXlow, exl, n);
      else     memset(fEXlow, 0, n);
      if (exh) memcpy(fEXhigh, exh, n);
      else     memset(fEXhigh, 0, n);
      if (eyl) memcpy(fEYlow, eyl, n);
      else     memset(fEYlow, 0, n);
      if (eyh) memcpy(fEYhigh, eyh, n);
      else     memset(fEYhigh, 0, n);

      if (exld) memcpy(fEXlowd, exld, n);
      else      memset(fEXlowd, 0, n);
      if (exhd) memcpy(fEXhighd, exhd, n);
      else      memset(fEXhighd, 0, n);
      if (eyld) memcpy(fEYlowd,  eyld, n);
      else      memset(fEYlowd, 0, n);
      if (eyhd) memcpy(fEYhighd, eyhd, n);
      else      memset(fEYhighd, 0, n);
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphBentErrors default destructor.

TGraphBentErrors::~TGraphBentErrors()
{
   delete [] fEXlow;
   delete [] fEXhigh;
   delete [] fEYlow;
   delete [] fEYhigh;

   delete [] fEXlowd;
   delete [] fEXhighd;
   delete [] fEYlowd;
   delete [] fEYhighd;
}


////////////////////////////////////////////////////////////////////////////////
/// Apply a function to all data points \f$ y = f(x,y) \f$.
///
/// Errors are calculated as \f$ eyh = f(x,y+eyh)-f(x,y) \f$ and
/// \f$ eyl = f(x,y)-f(x,y-eyl) \f$.
///
/// Special treatment has to be applied for the functions where the
/// role of "up" and "down" is reversed.
///
/// Function suggested/implemented by Miroslav Helbich <helbich@mail.desy.de>

void TGraphBentErrors::Apply(TF1 *f)
{
   Double_t x,y,exl,exh,eyl,eyh,eyl_new,eyh_new,fxy;

   if (fHistogram) {
      delete fHistogram;
      fHistogram = 0;
   }
   for (Int_t i=0;i<GetN();i++) {
      GetPoint(i,x,y);
      exl=GetErrorXlow(i);
      exh=GetErrorXhigh(i);
      eyl=GetErrorYlow(i);
      eyh=GetErrorYhigh(i);

      fxy = f->Eval(x,y);
      SetPoint(i,x,fxy);

      // in the case of the functions like y-> -1*y the roles of the
      // upper and lower error bars is reversed
      if (f->Eval(x,y-eyl)<f->Eval(x,y+eyh)) {
         eyl_new = TMath::Abs(fxy - f->Eval(x,y-eyl));
         eyh_new = TMath::Abs(f->Eval(x,y+eyh) - fxy);
      }
      else {
         eyh_new = TMath::Abs(fxy - f->Eval(x,y-eyl));
         eyl_new = TMath::Abs(f->Eval(x,y+eyh) - fxy);
      }

      //error on x doesn't change
      SetPointError(i,exl,exh,eyl_new,eyh_new);
   }
   if (gPad) gPad->Modified();
}


////////////////////////////////////////////////////////////////////////////////
/// Compute range.

void TGraphBentErrors::ComputeRange(Double_t &xmin, Double_t &ymin, Double_t &xmax, Double_t &ymax) const
{
   TGraph::ComputeRange(xmin,ymin,xmax,ymax);

   for (Int_t i=0;i<fNpoints;i++) {
      if (fX[i] -fEXlow[i] < xmin) {
         if (gPad && gPad->GetLogx()) {
            if (fEXlow[i] < fX[i]) xmin = fX[i]-fEXlow[i];
            else                   xmin = TMath::Min(xmin,fX[i]/3);
         } else {
            xmin = fX[i]-fEXlow[i];
         }
      }
      if (fX[i] +fEXhigh[i] > xmax) xmax = fX[i]+fEXhigh[i];
      if (fY[i] -fEYlow[i] < ymin) {
         if (gPad && gPad->GetLogy()) {
            if (fEYlow[i] < fY[i]) ymin = fY[i]-fEYlow[i];
            else                   ymin = TMath::Min(ymin,fY[i]/3);
         } else {
            ymin = fY[i]-fEYlow[i];
         }
      }
      if (fY[i] +fEYhigh[i] > ymax) ymax = fY[i]+fEYhigh[i];
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Copy and release.

void TGraphBentErrors::CopyAndRelease(Double_t **newarrays,
                                      Int_t ibegin, Int_t iend, Int_t obegin)
{
   CopyPoints(newarrays, ibegin, iend, obegin);
   if (newarrays) {
      delete[] fEXlow;
      fEXlow = newarrays[0];
      delete[] fEXhigh;
      fEXhigh = newarrays[1];
      delete[] fEYlow;
      fEYlow = newarrays[2];
      delete[] fEYhigh;
      fEYhigh = newarrays[3];
      delete[] fEXlowd;
      fEXlowd = newarrays[4];
      delete[] fEXhighd;
      fEXhighd = newarrays[5];
      delete[] fEYlowd;
      fEYlowd = newarrays[6];
      delete[] fEYhighd;
      fEYhighd = newarrays[7];
      delete[] fX;
      fX = newarrays[8];
      delete[] fY;
      fY = newarrays[9];
      delete[] newarrays;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Copy errors from `fE*** `to `arrays[***]`
/// or to `f***` Copy points.

Bool_t TGraphBentErrors::CopyPoints(Double_t **arrays,
                                    Int_t ibegin, Int_t iend, Int_t obegin)
{
   if (TGraph::CopyPoints(arrays ? arrays+8 : 0, ibegin, iend, obegin)) {
      Int_t n = (iend - ibegin)*sizeof(Double_t);
      if (arrays) {
         memmove(&arrays[0][obegin], &fEXlow[ibegin], n);
         memmove(&arrays[1][obegin], &fEXhigh[ibegin], n);
         memmove(&arrays[2][obegin], &fEYlow[ibegin], n);
         memmove(&arrays[3][obegin], &fEYhigh[ibegin], n);
         memmove(&arrays[4][obegin], &fEXlowd[ibegin], n);
         memmove(&arrays[5][obegin], &fEXhighd[ibegin], n);
         memmove(&arrays[6][obegin], &fEYlowd[ibegin], n);
         memmove(&arrays[7][obegin], &fEYhighd[ibegin], n);
      } else {
         memmove(&fEXlow[obegin], &fEXlow[ibegin], n);
         memmove(&fEXhigh[obegin], &fEXhigh[ibegin], n);
         memmove(&fEYlow[obegin], &fEYlow[ibegin], n);
         memmove(&fEYhigh[obegin], &fEYhigh[ibegin], n);
         memmove(&fEXlowd[obegin], &fEXlowd[ibegin], n);
         memmove(&fEXhighd[obegin], &fEXhighd[ibegin], n);
         memmove(&fEYlowd[obegin], &fEYlowd[ibegin], n);
         memmove(&fEYhighd[obegin], &fEYhighd[ibegin], n);
      }
      return kTRUE;
   } else {
      return kFALSE;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Should be called from ctors after `fNpoints` has been set.

Bool_t TGraphBentErrors::CtorAllocate(void)
{
   if (!fNpoints) {
      fEXlow = fEYlow = fEXhigh = fEYhigh = 0;
      fEXlowd = fEYlowd = fEXhighd = fEYhighd = 0;
      return kFALSE;
   }
   fEXlow = new Double_t[fMaxSize];
   fEYlow = new Double_t[fMaxSize];
   fEXhigh = new Double_t[fMaxSize];
   fEYhigh = new Double_t[fMaxSize];
   fEXlowd = new Double_t[fMaxSize];
   fEYlowd = new Double_t[fMaxSize];
   fEXhighd = new Double_t[fMaxSize];
   fEYhighd = new Double_t[fMaxSize];
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Protected function to perform the merge operation of a graph with asymmetric errors.

Bool_t TGraphBentErrors::DoMerge(const TGraph *g)
{
   if (g->GetN() == 0) return kFALSE;

   Double_t * exl = g->GetEXlow();
   Double_t * exh = g->GetEXhigh();
   Double_t * eyl = g->GetEYlow();
   Double_t * eyh = g->GetEYhigh();

   Double_t * exld = g->GetEXlowd();
   Double_t * exhd = g->GetEXhighd();
   Double_t * eyld = g->GetEYlowd();
   Double_t * eyhd = g->GetEYhighd();

   if (exl == 0 || exh == 0 || eyl == 0 || eyh == 0 ||
       exld == 0 || exhd == 0 || eyld == 0 || eyhd == 0) {
      if (g->IsA() != TGraph::Class() )
         Warning("DoMerge","Merging a %s is not compatible with a TGraphBentErrors - errors will be ignored",g->IsA()->GetName());
      return TGraph::DoMerge(g);
   }
   for (Int_t i = 0 ; i < g->GetN(); i++) {
      Int_t ipoint = GetN();
      Double_t x = g->GetX()[i];
      Double_t y = g->GetY()[i];
      SetPoint(ipoint, x, y);
      SetPointError(ipoint, exl[i],  exh[i],  eyl[i],  eyh[i],
                            exld[i], exhd[i], eyld[i], eyhd[i] );
   }

   return kTRUE;

}
////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It returns the error along X at point `i`.

Double_t TGraphBentErrors::GetErrorX(Int_t i) const
{
   if (i < 0 || i >= fNpoints) return -1;
   if (!fEXlow && !fEXhigh) return -1;
   Double_t elow=0, ehigh=0;
   if (fEXlow)  elow  = fEXlow[i];
   if (fEXhigh) ehigh = fEXhigh[i];
   return TMath::Sqrt(0.5*(elow*elow + ehigh*ehigh));
}


////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It returns the error along Y at point `i`.

Double_t TGraphBentErrors::GetErrorY(Int_t i) const
{
   if (i < 0 || i >= fNpoints) return -1;
   if (!fEYlow && !fEYhigh) return -1;
   Double_t elow=0, ehigh=0;
   if (fEYlow)  elow  = fEYlow[i];
   if (fEYhigh) ehigh = fEYhigh[i];
   return TMath::Sqrt(0.5*(elow*elow + ehigh*ehigh));
}


////////////////////////////////////////////////////////////////////////////////
/// Get high error on X[i].

Double_t TGraphBentErrors::GetErrorXhigh(Int_t i) const
{
   if (i<0 || i>fNpoints) return -1;
   if (fEXhigh) return fEXhigh[i];
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Get low error on X[i].

Double_t TGraphBentErrors::GetErrorXlow(Int_t i) const
{
   if (i<0 || i>fNpoints) return -1;
   if (fEXlow) return fEXlow[i];
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Get high error on Y[i].

Double_t TGraphBentErrors::GetErrorYhigh(Int_t i) const
{
   if (i<0 || i>fNpoints) return -1;
   if (fEYhigh) return fEYhigh[i];
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Get low error on Y[i].

Double_t TGraphBentErrors::GetErrorYlow(Int_t i) const
{
   if (i<0 || i>fNpoints) return -1;
   if (fEYlow) return fEYlow[i];
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Set zero values for point arrays in the range `[begin, end]`

void TGraphBentErrors::FillZero(Int_t begin, Int_t end,
                                 Bool_t from_ctor)
{
   if (!from_ctor) {
      TGraph::FillZero(begin, end, from_ctor);
   }
   Int_t n = (end - begin)*sizeof(Double_t);
   memset(fEXlow + begin, 0, n);
   memset(fEXhigh + begin, 0, n);
   memset(fEYlow + begin, 0, n);
   memset(fEYhigh + begin, 0, n);
   memset(fEXlowd + begin, 0, n);
   memset(fEXhighd + begin, 0, n);
   memset(fEYlowd + begin, 0, n);
   memset(fEYhighd + begin, 0, n);
}


////////////////////////////////////////////////////////////////////////////////
/// Print graph and errors values.

void TGraphBentErrors::Print(Option_t *) const
{
   for (Int_t i=0;i<fNpoints;i++) {
      printf("x[%d]=%g, y[%d]=%g, exl[%d]=%g, exh[%d]=%g, eyl[%d]=%g, eyh[%d]=%g\n"
         ,i,fX[i],i,fY[i],i,fEXlow[i],i,fEXhigh[i],i,fEYlow[i],i,fEYhigh[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply the values and errors of a TGraphBentErrors by a constant c1.
///
/// If option contains "x" the x values and errors are scaled
/// If option contains "y" the y values and errors are scaled
/// If option contains "xy" both x and y values and errors are scaled

void TGraphBentErrors::Scale(Double_t c1, Option_t *option)
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
   if (opt.Contains("y") && GetEYlow()) {
      for (Int_t i=0; i<GetN(); i++)
         GetEYlow()[i] *= c1;
   }
   if (opt.Contains("y") && GetEYhigh()) {
      for (Int_t i=0; i<GetN(); i++)
         GetEYhigh()[i] *= c1;
   }
   if (opt.Contains("x") && GetEXlowd()) {
      for (Int_t i=0; i<GetN(); i++)
         GetEXlowd()[i] *= c1;
   }
   if (opt.Contains("x") && GetEXhighd()) {
      for (Int_t i=0; i<GetN(); i++)
         GetEXhighd()[i] *= c1;
   }
   if (opt.Contains("y") && GetEYlowd()) {
      for (Int_t i=0; i<GetN(); i++)
         GetEYlowd()[i] *= c1;
   }
   if (opt.Contains("y") && GetEYhighd()) {
      for (Int_t i=0; i<GetN(); i++)
         GetEYhighd()[i] *= c1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out.

void TGraphBentErrors::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';
   out << "   " << std::endl;
   static Int_t frameNumber = 2000;
   frameNumber++;

   Int_t i;
   TString fXName    = TString(GetName()) + Form("_fx%d",frameNumber);
   TString fYName    = TString(GetName()) + Form("_fy%d",frameNumber);
   TString fElXName  = TString(GetName()) + Form("_felx%d",frameNumber);
   TString fElYName  = TString(GetName()) + Form("_fely%d",frameNumber);
   TString fEhXName  = TString(GetName()) + Form("_fehx%d",frameNumber);
   TString fEhYName  = TString(GetName()) + Form("_fehy%d",frameNumber);
   TString fEldXName = TString(GetName()) + Form("_feldx%d",frameNumber);
   TString fEldYName = TString(GetName()) + Form("_feldy%d",frameNumber);
   TString fEhdXName = TString(GetName()) + Form("_fehdx%d",frameNumber);
   TString fEhdYName = TString(GetName()) + Form("_fehdy%d",frameNumber);
   out << "   Double_t " << fXName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fX[i] << "," << std::endl;
   out << "   " << fX[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fYName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fY[i] << "," << std::endl;
   out << "   " << fY[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fElXName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fEXlow[i] << "," << std::endl;
   out << "   " << fEXlow[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fElYName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fEYlow[i] << "," << std::endl;
   out << "   " << fEYlow[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fEhXName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fEXhigh[i] << "," << std::endl;
   out << "   " << fEXhigh[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fEhYName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fEYhigh[i] << "," << std::endl;
   out << "   " << fEYhigh[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fEldXName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fEXlowd[i] << "," << std::endl;
   out << "   " << fEXlowd[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fEldYName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fEYlowd[i] << "," << std::endl;
   out << "   " << fEYlowd[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fEhdXName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fEXhighd[i] << "," << std::endl;
   out << "   " << fEXhighd[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fEhdYName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fEYhighd[i] << "," << std::endl;
   out << "   " << fEYhighd[fNpoints-1] << "};" << std::endl;

   if (gROOT->ClassSaved(TGraphBentErrors::Class())) out << "   ";
   else out << "   TGraphBentErrors *";
   out << "grbe = new TGraphBentErrors("<< fNpoints << ","
                                    << fXName     << ","  << fYName  << ","
                                    << fElXName   << ","  << fEhXName << ","
                                    << fElYName   << ","  << fEhYName << ","
                                    << fEldXName  << ","  << fEhdXName << ","
                                    << fEldYName  << ","  << fEhdYName << ");"
                                    << std::endl;

   out << "   grbe->SetName(" << quote << GetName() << quote << ");" << std::endl;
   out << "   grbe->SetTitle(" << quote << GetTitle() << quote << ");" << std::endl;

   SaveFillAttributes(out,"grbe",0,1001);
   SaveLineAttributes(out,"grbe",1,1,1);
   SaveMarkerAttributes(out,"grbe",1,1,1);

   if (fHistogram) {
      TString hname = fHistogram->GetName();
      hname += frameNumber;
      fHistogram->SetName(Form("Graph_%s",hname.Data()));
      fHistogram->SavePrimitive(out,"nodraw");
      out<<"   grbe->SetHistogram("<<fHistogram->GetName()<<");"<<std::endl;
      out<<"   "<<std::endl;
   }

   // save list of functions
   TIter next(fFunctions);
   TObject *obj;
   while ((obj = next())) {
      obj->SavePrimitive(out, Form("nodraw #%d\n",++frameNumber));
      if (obj->InheritsFrom("TPaveStats")) {
         out << "   grbe->GetListOfFunctions()->Add(ptstats);" << std::endl;
         out << "   ptstats->SetParent(grbe->GetListOfFunctions());" << std::endl;
      } else {
         TString objname;
         objname.Form("%s%d",obj->GetName(),frameNumber);
         if (obj->InheritsFrom("TF1")) {
            out << "   " << objname << "->SetParent(grbe);\n";
         }
         out << "   grbe->GetListOfFunctions()->Add("
             << objname << ");" << std::endl;
      }
   }

   const char *l = strstr(option,"multigraph");
   if (l) {
      out<<"   multigraph->Add(grbe,"<<quote<<l+10<<quote<<");"<<std::endl;
   } else {
      out<<"   grbe->Draw("<<quote<<option<<quote<<");"<<std::endl;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Set ex and ey values for point pointed by the mouse.

void TGraphBentErrors::SetPointError(Double_t exl, Double_t exh, Double_t eyl, Double_t eyh,
                                     Double_t exld, Double_t exhd, Double_t eyld, Double_t eyhd)
{
   Int_t px = gPad->GetEventX();
   Int_t py = gPad->GetEventY();

   //localize point to be deleted
   Int_t ipoint = -2;
   Int_t i;
   // start with a small window (in case the mouse is very close to one point)
   for (i=0;i<fNpoints;i++) {
      Int_t dpx = px - gPad->XtoAbsPixel(gPad->XtoPad(fX[i]));
      Int_t dpy = py - gPad->YtoAbsPixel(gPad->YtoPad(fY[i]));
      if (dpx*dpx+dpy*dpy < 25) {ipoint = i; break;}
   }
   if (ipoint == -2) return;

   fEXlow[ipoint]   = exl;
   fEYlow[ipoint]   = eyl;
   fEXhigh[ipoint]  = exh;
   fEYhigh[ipoint]  = eyh;
   fEXlowd[ipoint]  = exld;
   fEXhighd[ipoint] = exhd;
   fEYlowd[ipoint]  = eyld;
   fEYhighd[ipoint] = eyhd;
   gPad->Modified();
}


////////////////////////////////////////////////////////////////////////////////
/// Set ex and ey values for point number `i`.

void TGraphBentErrors::SetPointError(Int_t i, Double_t exl, Double_t exh, Double_t eyl, Double_t eyh,
                                     Double_t exld, Double_t exhd, Double_t eyld, Double_t eyhd)
{
   if (i < 0) return;
   if (i >= fNpoints) {
   // re-allocate the object
      TGraphBentErrors::SetPoint(i,0,0);
   }
   fEXlow[i]   = exl;
   fEYlow[i]   = eyl;
   fEXhigh[i]  = exh;
   fEYhigh[i]  = eyh;
   fEXlowd[i]  = exld;
   fEXhighd[i] = exhd;
   fEYlowd[i]  = eyld;
   fEYhighd[i] = eyhd;
}


////////////////////////////////////////////////////////////////////////////////
/// Swap points.

void TGraphBentErrors::SwapPoints(Int_t pos1, Int_t pos2)
{
   SwapValues(fEXlow,  pos1, pos2);
   SwapValues(fEXhigh, pos1, pos2);
   SwapValues(fEYlow,  pos1, pos2);
   SwapValues(fEYhigh, pos1, pos2);

   SwapValues(fEXlowd,  pos1, pos2);
   SwapValues(fEXhighd, pos1, pos2);
   SwapValues(fEYlowd,  pos1, pos2);
   SwapValues(fEYhighd, pos1, pos2);

   TGraph::SwapPoints(pos1, pos2);
}
