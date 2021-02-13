// @(#)root/hist:$Id$
// Author: Rene Brun   14/01/2002
/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TH1.h"
#include "TMath.h"
#include "THLimitsFinder.h"

THLimitsFinder *THLimitsFinder::fgLimitsFinder = 0;

ClassImp(THLimitsFinder);

/** \class THLimitsFinder
    \ingroup Hist
Class to compute nice axis limits.

This class is called by default by the histogramming system
and also by TTree::Draw, TTreePlayer::DrawSelect. TGaxis use it also to
optimize axis labels.

A different finder may be specified via THLimitsFinder::SetFinder.
*/

////////////////////////////////////////////////////////////////////////////////

THLimitsFinder::THLimitsFinder()
{
}

////////////////////////////////////////////////////////////////////////////////

THLimitsFinder::~THLimitsFinder()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Compute the best axis limits for the X axis.
///
/// If the bit kIsInteger is set, the number of channels is also recomputed.
/// The axis parameters are replaced by the optimized parameters.
///
/// Example:
/// With the input parameters xmin=-1.467 and xmax=2.344, the function
/// will compute better limits -1.8 and 2.7 and store them in the axis.

Int_t THLimitsFinder::FindGoodLimits(TH1 *h, Double_t xmin, Double_t xmax)
{
   Int_t newbins;
   TAxis *xaxis = h->GetXaxis();

   if (xmin >= xmax) {
      if (xaxis->GetLabels()) {xmin  = 0; xmax  = xmin +xaxis->GetNbins();}
      else                    {xmin -= 1; xmax += 1;}
   }

   THLimitsFinder::OptimizeLimits(xaxis->GetNbins(),
                                  newbins,xmin,xmax,
                                  xaxis->TestBit(TAxis::kIsInteger));

   h->SetBins(newbins,xmin,xmax);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the best axis limits for the X and Y axis.
///
/// If the bit kIsInteger is set, the number of channels is also recomputed.
/// The axis parameters are replaced by the optimized parameters

Int_t THLimitsFinder::FindGoodLimits(TH1 *h, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax)
{
   Int_t newbinsx,newbinsy;
   TAxis *xaxis = h->GetXaxis();
   TAxis *yaxis = h->GetYaxis();

   if (xmin >= xmax) {
      if (xaxis->GetLabels()) {xmin  = 0; xmax  = xmin +xaxis->GetNbins();}
      else                    {xmin -= 1; xmax += 1;}
   }
   if (ymin >= ymax) {
      if (yaxis->GetLabels()) {ymin  = 0; ymax  = ymin +yaxis->GetNbins();}
      else                    {ymin -= 1; ymax += 1;}
   }

   THLimitsFinder::OptimizeLimits(xaxis->GetNbins(),
                                  newbinsx,xmin,xmax,
                                  xaxis->TestBit(TAxis::kIsInteger));

   THLimitsFinder::OptimizeLimits(yaxis->GetNbins(),
                                  newbinsy,ymin,ymax,
                                  yaxis->TestBit(TAxis::kIsInteger));

   h->SetBins(newbinsx,xmin,xmax,newbinsy,ymin,ymax);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the best axis limits for the X, Y and Z axis.
///
/// If the bit kIsInteger is set, the number of channels is also recomputed.
/// The axis parameters are replaced by the optimized parameters

Int_t THLimitsFinder::FindGoodLimits(TH1 *h, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Double_t zmin, Double_t zmax)
{
   Int_t newbinsx,newbinsy,newbinsz;
   TAxis *xaxis = h->GetXaxis();
   TAxis *yaxis = h->GetYaxis();
   TAxis *zaxis = h->GetZaxis();

   if (xmin >= xmax) {
      if (xaxis->GetLabels()) {xmin  = 0; xmax  = xmin +xaxis->GetNbins();}
      else                    {xmin -= 1; xmax += 1;}
   }
   if (ymin >= ymax) {
      if (yaxis->GetLabels()) {ymin  = 0; ymax  = ymin +yaxis->GetNbins();}
      else                    {ymin -= 1; ymax += 1;}
   }
   if (zmin >= zmax) {
      if (zaxis->GetLabels()) {zmin  = 0; zmax  = zmin +zaxis->GetNbins();}
      else                    {zmin -= 1; zmax += 1;}
   }

   THLimitsFinder::OptimizeLimits(xaxis->GetNbins(),
                                  newbinsx,xmin,xmax,
                                  xaxis->TestBit(TAxis::kIsInteger));

   THLimitsFinder::OptimizeLimits(yaxis->GetNbins(),
                                  newbinsy,ymin,ymax,
                                  yaxis->TestBit(TAxis::kIsInteger));

   THLimitsFinder::OptimizeLimits(zaxis->GetNbins(),
                                  newbinsz,zmin,zmax,
                                  zaxis->TestBit(TAxis::kIsInteger));

   h->SetBins(newbinsx,xmin,xmax,newbinsy,ymin,ymax,newbinsz,zmin,zmax);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to the current finder.
///
/// Create one if none exists
/// Use SetLimitsFinder to set a user defined finder.

THLimitsFinder *THLimitsFinder::GetLimitsFinder()
{
   if (!fgLimitsFinder) fgLimitsFinder = new THLimitsFinder();
   return fgLimitsFinder;
}

////////////////////////////////////////////////////////////////////////////////
/// This static function can be used to specify a finder derived from THLimitsFinder.
///
/// The finder may redefine the functions FindGoodLimits.
/// Note that the redefined functions may call THLimitsFinder::FindGoodLimits.

void THLimitsFinder::SetLimitsFinder(THLimitsFinder *finder)
{
   fgLimitsFinder = finder;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function to compute reasonable axis limits
///
/// Input parameters:
///
///  \param[in]   A1,A2          : Original axis limits
///  \param[in]   BinLow,BinHigh : Optimized axis limits. They should be initialized by the
///                                calling method for instance to 0.
///  \param[out]  nold           : Original number of divisions.
///  \param[out]  nbins          : Optimized number of divisions.
///  \param[out]  BinWidth       : Optimized bin width. It should be initialized by the
///                                calling method for instance to 0.
///  \param[in]  option          : "T" means Time axis.

void THLimitsFinder::Optimize(Double_t A1,  Double_t A2,  Int_t nold ,
                              Double_t &BinLow, Double_t &BinHigh,
                              Int_t &nbins, Double_t &BinWidth,
                              Option_t *option)
{
   Int_t lwid, kwid;
   Int_t ntemp = 0;
   Int_t jlog  = 0;
   Double_t siground = 0;
   Double_t alb, awidth, sigfig;
   Double_t timemulti = 1;
   Int_t roundmode =0;

   Int_t optionTime;
   if(strchr(option,'t')) optionTime = 1;  else optionTime = 0;

   nbins = nold;

   Double_t al = TMath::Min(A1,A2);
   Double_t ah = TMath::Max(A1,A2);
   if (al == ah) ah = al+1;
   // if nold  ==  -1 , program uses binwidth input from calling routine
   if (nold == -1 && BinWidth > 0 ) goto L90;
   ntemp = TMath::Max(nold,2);
   if (ntemp < 1) ntemp = 1;

L20:
   awidth = (ah-al)/Double_t(ntemp);
   timemulti = 1;
   if (awidth >= FLT_MAX) goto LOK;  //in float.h
   if (awidth <= 0)       goto LOK;

//      If time representation, bin width should be rounded to seconds
//      minutes, hours or days

   if (optionTime && awidth>=60) { // if width in seconds, treat it as normal
      //   width in minutes
      awidth /= 60; timemulti *=60;
      roundmode = 1; // round minutes (60)
      //   width in hours ?
      if (awidth>=60) {
         awidth /= 60; timemulti *= 60;
         roundmode = 2; // round hours (24)
         //   width in days ?
         if (awidth>=24) {
            awidth /= 24; timemulti *= 24;
            roundmode = 3; // round days (30)
            //   width in months ?
            if (awidth>=30.43685) { // Mean month length in 1900.
               awidth /= 30.43685; timemulti *= 30.43685;
               roundmode = 2; // round months (12)
               //   width in years ?
               if (awidth>=12) {
                  awidth /= 12; timemulti *= 12;
                  roundmode = 0; // round years (10)
               }
            }
         }
      }
   }
//      Get nominal bin width in exponential form

   jlog   = Int_t(TMath::Log10(awidth));
   if (jlog <-200 || jlog > 200) {
      BinLow   = 0;
      BinHigh  = 1;
      BinWidth = 0.01;
      nbins    = 100;
      return;
   }
   if (awidth <= 1 && (!optionTime || timemulti==1) ) jlog--;
   sigfig = awidth*TMath::Power(10,-jlog) -1e-10;
   //in the above statement, it is important to subtract 1e-10
   //to avoid precision problems if the tests below

//      Round mantissa

   switch (roundmode) {

//      Round mantissa up to 1, 1.5, 2, 3, or 6 in case of minutes
      case 1: // case 60
         if      (sigfig <= 1)    siground = 1;
         else if (sigfig <= 1.5 && jlog==1)    siground = 1.5;
         else if (sigfig <= 2)    siground = 2;
         else if (sigfig <= 3 && jlog ==1)    siground = 3;
         else if (sigfig <= 5 && sigfig>3 && jlog ==0) siground = 5; //added (Damir in 3.10/02)
         else if (jlog==0)        {siground = 1; jlog++;}
         else                     siground = 6;
         break;
      case 2: // case 12 and 24

//      Round mantissa up to 1, 1.2, 2, 2.4, 3 or 6 in case of hours or months
         if      (sigfig <= 1 && jlog==0)    siground = 1;
         else if (sigfig <= 1.2 && jlog==1)    siground = 1.2;
         else if (sigfig <= 2 && jlog==0)    siground = 2;
         else if (sigfig <= 2.4 && jlog==1)    siground = 2.4;
         else if (sigfig <= 3)    siground = 3;
         else if (sigfig <= 6)    siground = 6;
         else if (jlog==0)        siground = 12;
         else                     siground = 2.4;
         break;

//-      Round mantissa up to 1, 1.4, 2, or 7 in case of days (weeks)
      case 3: // case 30
         if      (sigfig <= 1 && jlog==0)    siground = 1;
         else if (sigfig <= 1.4 && jlog==1)    siground = 1.4;
         else if (sigfig <= 3 && jlog ==1)    siground = 3;
         else                     siground = 7;
         break;
      default :

//      Round mantissa up to 1, 2, 2.5, 5, or 10 in case of decimal number
         if      (sigfig <= 1)    siground = 1;
         else if (sigfig <= 2)    siground = 2;
         else if (sigfig <= 5 && (!optionTime || jlog<1))  siground = 5;
         else if (sigfig <= 6 && optionTime && jlog==1)    siground = 6;
         else                    {siground = 1;   jlog++; }
         break;
   }

   BinWidth = siground*TMath::Power(10,jlog);
   if (optionTime) BinWidth *= timemulti;

//      Get new bounds from new width BinWidth

L90:
   alb  = al/BinWidth;
   if (TMath::Abs(alb) > 1e9) {
      BinLow  = al;
      BinHigh = ah;
      if (nbins > 10*nold && nbins > 10000) nbins = nold;
      return;
   }
   lwid   = Int_t(alb);
   if (alb < 0) lwid--;
   BinLow     = BinWidth*Double_t(lwid);
   alb        = ah/BinWidth + 1.00001;
   kwid = Int_t(alb);
   if (alb < 0) kwid--;
   BinHigh = BinWidth*Double_t(kwid);
   nbins = kwid - lwid;
   if (nold == -1) goto LOK;
   if (nold <= 5) {          //    Request for one bin is difficult case
      if (nold > 1 || nbins == 1)goto LOK;
      BinWidth = BinWidth*2;
      nbins    = 1;
      goto LOK;
   }
   if (2*nbins == nold && !optionTime) {ntemp++; goto L20; }

LOK:
   Double_t oldBinLow = BinLow;
   Double_t oldBinHigh = BinHigh;
   Int_t oldnbins = nbins;

   Double_t atest = BinWidth*0.0001;
   //if (TMath::Abs(BinLow-A1)  >= atest) { BinLow  += BinWidth;  nbins--; } //replaced by Damir in 3.10/02
   //if (TMath::Abs(BinHigh-A2) >= atest) { BinHigh -= BinWidth;  nbins--; } //by the next two lines
   if (al-BinLow  >= atest) { BinLow  += BinWidth;  nbins--; }
   if (BinHigh-ah >= atest) { BinHigh -= BinWidth;  nbins--; }
   if (!optionTime && BinLow >= BinHigh) {
      //this case may happen when nbins <=5
      BinLow = oldBinLow;
      BinHigh = oldBinHigh;
      nbins = oldnbins;
   }
   else if (optionTime && BinLow>=BinHigh) {
      nbins = 2*oldnbins;
      BinHigh = oldBinHigh;
      BinLow = oldBinLow;
      BinWidth = (oldBinHigh - oldBinLow)/nbins;
      atest = BinWidth*0.0001;
      if (al-BinLow  >= atest) { BinLow  += BinWidth;  nbins--; }
      if (BinHigh-ah >= atest) { BinHigh -= BinWidth;  nbins--; }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Optimize axis limits.
///
/// When isInter=kTRUE, the function makes an integer binwidth
/// and recompute the number of bins accordingly.

void THLimitsFinder::OptimizeLimits(Int_t nbins, Int_t &newbins, Double_t &xmin, Double_t &xmax, Bool_t isInteger)
{
   Double_t binlow = 0,binhigh = 0,binwidth=0;
   Int_t n=0;
   Double_t dx = 0.1*(xmax-xmin);
   if (isInteger) dx = 5*(xmax-xmin)/nbins;
   Double_t umin = xmin - dx;
   Double_t umax = xmax + dx;
   if (umin < 0 && xmin >= 0) umin = 0;
   if (umax > 0 && xmax <= 0) umax = 0;

   THLimitsFinder::Optimize(umin,umax,nbins,binlow,binhigh,n,binwidth,"");

   if (binwidth <= 0 || binwidth > 1.e+39) {
      xmin = -1;
      xmax = 1;
   } else {
      xmin    = binlow;
      xmax    = binhigh;
   }
   if (isInteger) {
      Long64_t ixmin = Long64_t(xmin);
      Long64_t ixmax = Long64_t(xmax);
      Double_t dxmin = Double_t(ixmin);
      Double_t dxmax = Double_t(ixmax);
      if (xmin < 0 && xmin != dxmin) xmin = dxmin - 1;
      else                           xmin = dxmin;
      if (xmax > 0 && xmax != dxmax)      xmax = dxmax + 1;
      else if (xmax ==0 && xmax == dxmax) xmax = 1;
      else                                xmax = dxmax;
      if (xmin >= xmax) xmax = xmin+1;
      Long64_t bw = Long64_t((xmax-xmin)/nbins);
      if (bw == 0) bw = 1;
      nbins = Int_t((xmax-xmin)/bw);
      if (xmin +nbins*bw < umax) {nbins++; xmax = xmin +nbins*bw;}
      if (xmin > umin)           {nbins++; xmin = xmax -nbins*bw;}
   }
   newbins = nbins;
}
