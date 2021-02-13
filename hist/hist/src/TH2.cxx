// @(#)root/hist:$Id$
// Author: Rene Brun   26/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TBuffer.h"
#include "TClass.h"
#include "THashList.h"
#include "TH2.h"
#include "TVirtualPad.h"
#include "TF2.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TMatrixFBase.h"
#include "TMatrixDBase.h"
#include "THLimitsFinder.h"
#include "TError.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TVirtualHistPainter.h"
#include "snprintf.h"

ClassImp(TH2);

/** \addtogroup Hist
@{
\class TH2C
\brief 2-D histogram with a byte per channel (see TH1 documentation)
\class TH2S
\brief 2-D histogram with a short per channel (see TH1 documentation)
\class TH2I
\brief 2-D histogram with an int per channel (see TH1 documentation)}
\class TH2F
\brief 2-D histogram with a float per channel (see TH1 documentation)}
\class TH2D
\brief 2-D histogram with a double per channel (see TH1 documentation)}
@}
*/

/** \class TH2
 Service class for 2-Dim histogram classes

- TH2C a 2-D histogram with one byte per cell (char)
- TH2S a 2-D histogram with two bytes per cell (short integer)
- TH2I a 2-D histogram with four bytes per cell (32 bits integer)
- TH2F a 2-D histogram with four bytes per cell (float)
- TH2D a 2-D histogram with eight bytes per cell (double)
*/


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2::TH2()
{
   fDimension   = 2;
   fScalefactor = 1;
   fTsumwy      = fTsumwy2 = fTsumwxy = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// See comments in the TH1 base class constructors.

TH2::TH2(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
                                     ,Int_t nbinsy,Double_t ylow,Double_t yup)
     :TH1(name,title,nbinsx,xlow,xup)
{
   fDimension   = 2;
   fScalefactor = 1;
   fTsumwy      = fTsumwy2 = fTsumwxy = 0;
   if (nbinsy <= 0) {Warning("TH2","nbinsy is <=0 - set to nbinsy = 1"); nbinsy = 1; }
   fYaxis.Set(nbinsy,ylow,yup);
   fNcells      = fNcells*(nbinsy+2); // fNCells is set in the TH1 constructor
}


////////////////////////////////////////////////////////////////////////////////
/// See comments in the TH1 base class constructors.

TH2::TH2(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
                                     ,Int_t nbinsy,Double_t ylow,Double_t yup)
     :TH1(name,title,nbinsx,xbins)
{
   fDimension   = 2;
   fScalefactor = 1;
   fTsumwy      = fTsumwy2 = fTsumwxy = 0;
   if (nbinsy <= 0) {Warning("TH2","nbinsy is <=0 - set to nbinsy = 1"); nbinsy = 1; }
   fYaxis.Set(nbinsy,ylow,yup);
   fNcells      = fNcells*(nbinsy+2); // fNCells is set in the TH1 constructor
}


////////////////////////////////////////////////////////////////////////////////
/// See comments in the TH1 base class constructors.

TH2::TH2(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
                                     ,Int_t nbinsy,const Double_t *ybins)
     :TH1(name,title,nbinsx,xlow,xup)
{
   fDimension   = 2;
   fScalefactor = 1;
   fTsumwy      = fTsumwy2 = fTsumwxy = 0;
   if (nbinsy <= 0) {Warning("TH2","nbinsy is <=0 - set to nbinsy = 1"); nbinsy = 1; }
   if (ybins) fYaxis.Set(nbinsy,ybins);
   else       fYaxis.Set(nbinsy,0,1);
   fNcells      = fNcells*(nbinsy+2); // fNCells is set in the TH1 constructor
}


////////////////////////////////////////////////////////////////////////////////
/// See comments in the TH1 base class constructors.

TH2::TH2(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
                                           ,Int_t nbinsy,const Double_t *ybins)
     :TH1(name,title,nbinsx,xbins)
{
   fDimension   = 2;
   fScalefactor = 1;
   fTsumwy      = fTsumwy2 = fTsumwxy = 0;
   if (nbinsy <= 0) {Warning("TH2","nbinsy is <=0 - set to nbinsy = 1"); nbinsy = 1; }
   if (ybins) fYaxis.Set(nbinsy,ybins);
   else       fYaxis.Set(nbinsy,0,1);
   fNcells      = fNcells*(nbinsy+2); // fNCells is set in the TH1 constructor
}


////////////////////////////////////////////////////////////////////////////////
/// See comments in the TH1 base class constructors.

TH2::TH2(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins
                                           ,Int_t nbinsy,const Float_t *ybins)
     :TH1(name,title,nbinsx,xbins)
{
   fDimension   = 2;
   fScalefactor = 1;
   fTsumwy      = fTsumwy2 = fTsumwxy = 0;
   if (nbinsy <= 0) {Warning("TH2","nbinsy is <=0 - set to nbinsy = 1"); nbinsy = 1; }
   if (ybins) fYaxis.Set(nbinsy,ybins);
   else       fYaxis.Set(nbinsy,0,1);
   fNcells      = fNcells*(nbinsy+2); // fNCells is set in the TH1 constructor.
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.
/// The list of functions is not copied. (Use Clone if needed)

TH2::TH2(const TH2 &h) : TH1()
{
   ((TH2&)h).Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH2::~TH2()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Fill histogram with all entries in the buffer.
///  - action = -1 histogram is reset and refilled from the buffer (called by THistPainter::Paint)
///  - action =  0 histogram is filled from the buffer
///  - action =  1 histogram is filled and buffer is deleted
///                The buffer is automatically deleted when the number of entries
///                in the buffer is greater than the number of entries in the histogram

Int_t TH2::BufferEmpty(Int_t action)
{
   // do we need to compute the bin size?
   if (!fBuffer) return 0;
   Int_t nbentries = (Int_t)fBuffer[0];

   // nbentries correspond to the number of entries of histogram

   if (nbentries == 0) return 0;
   if (nbentries < 0 && action == 0) return 0;    // case histogram has been already filled from the buffer

   Double_t *buffer = fBuffer;
   if (nbentries < 0) {
      nbentries  = -nbentries;
      //  a reset might call BufferEmpty() giving an infinite loop
      // Protect it by setting fBuffer = 0
      fBuffer=0;
       //do not reset the list of functions
      Reset("ICES");
      fBuffer = buffer;
   }

   if (CanExtendAllAxes() || fXaxis.GetXmax() <= fXaxis.GetXmin() || fYaxis.GetXmax() <= fYaxis.GetXmin()) {
      //find min, max of entries in buffer
      Double_t xmin = fBuffer[2];
      Double_t xmax = xmin;
      Double_t ymin = fBuffer[3];
      Double_t ymax = ymin;
      for (Int_t i=1;i<nbentries;i++) {
         Double_t x = fBuffer[3*i+2];
         if (x < xmin) xmin = x;
         if (x > xmax) xmax = x;
         Double_t y = fBuffer[3*i+3];
         if (y < ymin) ymin = y;
         if (y > ymax) ymax = y;
      }
      if (fXaxis.GetXmax() <= fXaxis.GetXmin() || fYaxis.GetXmax() <= fYaxis.GetXmin()) {
         THLimitsFinder::GetLimitsFinder()->FindGoodLimits(this,xmin,xmax,ymin,ymax);
      } else {
         fBuffer = 0;
         Int_t keep = fBufferSize; fBufferSize = 0;
         if (xmin <  fXaxis.GetXmin()) ExtendAxis(xmin,&fXaxis);
         if (xmax >= fXaxis.GetXmax()) ExtendAxis(xmax,&fXaxis);
         if (ymin <  fYaxis.GetXmin()) ExtendAxis(ymin,&fYaxis);
         if (ymax >= fYaxis.GetXmax()) ExtendAxis(ymax,&fYaxis);
         fBuffer = buffer;
         fBufferSize = keep;
      }
   }

   fBuffer = 0;
   for (Int_t i=0;i<nbentries;i++) {
      Fill(buffer[3*i+2],buffer[3*i+3],buffer[3*i+1]);
   }
   fBuffer = buffer;

   if (action > 0) { delete [] fBuffer; fBuffer = 0; fBufferSize = 0;}
   else {
      if (nbentries == (Int_t)fEntries) fBuffer[0] = -nbentries;
      else                              fBuffer[0] = 0;
   }
   return nbentries;
}


////////////////////////////////////////////////////////////////////////////////
/// accumulate arguments in buffer. When buffer is full, empty the buffer
/// ~~~ {.cpp}
/// fBuffer[0] = number of entries in buffer
/// fBuffer[1] = w of first entry
/// fBuffer[2] = x of first entry
/// fBuffer[3] = y of first entry
/// ~~~

Int_t TH2::BufferFill(Double_t x, Double_t y, Double_t w)
{
   if (!fBuffer) return -3;
   Int_t nbentries = (Int_t)fBuffer[0];
   if (nbentries < 0) {
      nbentries  = -nbentries;
      fBuffer[0] =  nbentries;
      if (fEntries > 0) {
         Double_t *buffer = fBuffer; fBuffer=0;
         Reset("ICES");
         fBuffer = buffer;
      }
   }
   if (3*nbentries+3 >= fBufferSize) {
      BufferEmpty(1);
      return Fill(x,y,w);
   }
   fBuffer[3*nbentries+1] = w;
   fBuffer[3*nbentries+2] = x;
   fBuffer[3*nbentries+3] = y;
   fBuffer[0] += 1;
   return -3;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy.

void TH2::Copy(TObject &obj) const
{
   TH1::Copy(obj);
   ((TH2&)obj).fScalefactor = fScalefactor;
   ((TH2&)obj).fTsumwy      = fTsumwy;
   ((TH2&)obj).fTsumwy2     = fTsumwy2;
   ((TH2&)obj).fTsumwxy     = fTsumwxy;
}


////////////////////////////////////////////////////////////////////////////////
/// Invalid Fill method.

Int_t TH2::Fill(Double_t )
{
   Error("Fill", "Invalid signature - do nothing");
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment cell defined by x,y by 1.
///
///  - if x or/and y is less than the low-edge of the corresponding axis first bin,
///    the Underflow cell is incremented.
///  - if x or/and y is equal to or greater than the upper edge of corresponding axis last bin,
///    the Overflow cell is incremented.
///
/// -  If the storage of the sum of squares of weights has been triggered,
///    via the function Sumw2, then the sum of the squares of weights is incremented
///    by 1 in the cell corresponding to x,y.
///
/// The function returns the corresponding global bin number which has its content
/// incremented by 1

Int_t TH2::Fill(Double_t x,Double_t y)
{
   if (fBuffer) return BufferFill(x,y,1);

   Int_t binx, biny, bin;
   fEntries++;
   binx = fXaxis.FindBin(x);
   biny = fYaxis.FindBin(y);
   if (binx <0 || biny <0) return -1;
   bin  = biny*(fXaxis.GetNbins()+2) + binx;
   AddBinContent(bin);
   if (fSumw2.fN) ++fSumw2.fArray[bin];
   if (binx == 0 || binx > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (biny == 0 || biny > fYaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   ++fTsumw;
   ++fTsumw2;
   fTsumwx  += x;
   fTsumwx2 += x*x;
   fTsumwy  += y;
   fTsumwy2 += y*y;
   fTsumwxy += x*y;
   return bin;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment cell defined by x,y by a weight w.
///
///  - if x or/and y is less than the low-edge of the corresponding axis first bin,
///    the Underflow cell is incremented.
///  - if x or/and y is equal to or greater than the upper edge of corresponding axis last bin,
///    the Overflow cell is incremented.
///
///   - If the weight is not equal to 1, the storage of the sum of squares of
///     weights is automatically triggered and the sum of the squares of weights is incremented
///     by w^2 in the bin corresponding to x,y
///
/// The function returns the corresponding global bin number which has its content
/// incremented by w

Int_t TH2::Fill(Double_t x, Double_t y, Double_t w)
{
   if (fBuffer) return BufferFill(x,y,w);

   Int_t binx, biny, bin;
   fEntries++;
   binx = fXaxis.FindBin(x);
   biny = fYaxis.FindBin(y);
   if (binx <0 || biny <0) return -1;
   bin  = biny*(fXaxis.GetNbins()+2) + binx;
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();   // must be called before AddBinContent
   if (fSumw2.fN) fSumw2.fArray[bin] += w*w;
   AddBinContent(bin,w);
   if (binx == 0 || binx > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (biny == 0 || biny > fYaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   Double_t z= w;
   fTsumw   += z;
   fTsumw2  += z*z;
   fTsumwx  += z*x;
   fTsumwx2 += z*x*x;
   fTsumwy  += z*y;
   fTsumwy2 += z*y*y;
   fTsumwxy += z*x*y;
   return bin;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment cell defined by namex,namey by a weight w
///
///  - if x or/and y is less than the low-edge of the corresponding axis first bin,
///    the Underflow cell is incremented.
///  - if x or/and y is equal to or greater than the upper edge of corresponding axis last bin,
///    the Overflow cell is incremented.
///
///   - If the weight is not equal to 1, the storage of the sum of squares of
///     weights is automatically triggered and the sum of the squares of weights is incremented
///     by w^2 in the bin corresponding to namex,namey
///
/// The function returns the corresponding global bin number which has its content
/// incremented by w

Int_t TH2::Fill(const char *namex, const char *namey, Double_t w)
{
   Int_t binx, biny, bin;
   fEntries++;
   binx = fXaxis.FindBin(namex);
   biny = fYaxis.FindBin(namey);
   if (binx <0 || biny <0) return -1;
   bin  = biny*(fXaxis.GetNbins()+2) + binx;
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();   // must be called before AddBinContent
   if (fSumw2.fN) fSumw2.fArray[bin] += w*w;
   AddBinContent(bin,w);
   if (binx == 0 || binx > fXaxis.GetNbins()) return -1;
   if (biny == 0 || biny > fYaxis.GetNbins()) return -1;

   Double_t z= w;
   fTsumw   += z;
   fTsumw2  += z*z;
   // skip computation of the statistics along axis that have labels (can be extended and are aphanumeric)
   UInt_t labelBitMask = GetAxisLabelStatus();
   if (labelBitMask != (TH1::kXaxis | TH1::kYaxis)) {
      Double_t x = (labelBitMask & TH1::kXaxis) ? 0 : fXaxis.GetBinCenter(binx);
      Double_t y = (labelBitMask & TH1::kYaxis) ? 0 : fYaxis.GetBinCenter(biny);
      fTsumwx += z * x;
      fTsumwx2 += z * x * x;
      fTsumwy += z * y;
      fTsumwy2 += z * y * y;
      fTsumwx2 += z * x * x;
   }
   return bin;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment cell defined by namex,y by a weight w
///
///  - if x or/and y is less than the low-edge of the corresponding axis first bin,
///    the Underflow cell is incremented.
///  - if x or/and y is equal to or greater than the upper edge of corresponding axis last bin,
///    the Overflow cell is incremented.
///
///  - If the weight is not equal to 1, the storage of the sum of squares of
///    weights is automatically triggered and the sum of the squares of weights is incremented
///    by w^2 in the bin corresponding to namex,y
///
/// The function returns the corresponding global bin number which has its content
/// incremented by w

Int_t TH2::Fill(const char *namex, Double_t y, Double_t w)
{
   Int_t binx, biny, bin;
   fEntries++;
   binx = fXaxis.FindBin(namex);
   biny = fYaxis.FindBin(y);
   if (binx <0 || biny <0) return -1;
   bin  = biny*(fXaxis.GetNbins()+2) + binx;
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();   // must be called before AddBinContent
   if (fSumw2.fN) fSumw2.fArray[bin] += w*w;
   AddBinContent(bin,w);
   if (binx == 0 || binx > fXaxis.GetNbins()) return -1;
   if (biny == 0 || biny > fYaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   Double_t z= w; //(w > 0 ? w : -w);
   fTsumw   += z;
   fTsumw2  += z*z;
   fTsumwy  += z*y;
   fTsumwy2 += z*y*y;
   if (!fXaxis.CanExtend() || !fXaxis.IsAlphanumeric()) {
      Double_t x = fXaxis.GetBinCenter(binx);
      fTsumwx += z * x;
      fTsumwx2 += z * x * x;
      fTsumwxy += z * x * y;
   }
   return bin;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment cell defined by x,namey by a weight w
///
///  - if x or/and y is less than the low-edge of the corresponding axis first bin,
///    the Underflow cell is incremented.
///  - if x or/and y is equal to or greater than the upper edge of corresponding axis last bin,
///    the Overflow cell is incremented.
///
///  - If the weight is not equal to 1, the storage of the sum of squares of
///    weights is automatically triggered and the sum of the squares of weights is incremented
///    by w^2 in the bin corresponding to x,y.
///
/// The function returns the corresponding global bin number which has its content
/// incremented by w

Int_t TH2::Fill(Double_t x, const char *namey, Double_t w)
{
   Int_t binx, biny, bin;
   fEntries++;
   binx = fXaxis.FindBin(x);
   biny = fYaxis.FindBin(namey);
   if (binx <0 || biny <0) return -1;
   bin  = biny*(fXaxis.GetNbins()+2) + binx;
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();   // must be called before AddBinContent
   if (fSumw2.fN) fSumw2.fArray[bin] += w*w;
   AddBinContent(bin,w);
   if (binx == 0 || binx > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (biny == 0 || biny > fYaxis.GetNbins()) return -1;

   Double_t z= w; //(w > 0 ? w : -w);
   fTsumw   += z;
   fTsumw2  += z*z;
   fTsumwx  += z*x;
   fTsumwx2 += z*x*x;
   if (!fYaxis.CanExtend() || !fYaxis.IsAlphanumeric()) {
      Double_t y = fYaxis.GetBinCenter(biny);
      fTsumwy += z * y;
      fTsumwy2 += z * y * y;
      fTsumwxy += z * x * y;
   }
   return bin;
}


////////////////////////////////////////////////////////////////////////////////
/// Fill a 2-D histogram with an array of values and weights.
///
///  - ntimes:  number of entries in arrays x and w (array size must be ntimes*stride)
///  - x:       array of x values to be histogrammed
///  - y:       array of y values to be histogrammed
///  - w:       array of weights
///  - stride:  step size through arrays x, y and w
///
///   - If the weight is not equal to 1, the storage of the sum of squares of
///     weights is automatically triggered and the sum of the squares of weights is incremented
///     by w[i]^2 in the bin corresponding to x[i],y[i].
///   - If w is NULL each entry is assumed a weight=1
///
/// NB: function only valid for a TH2x object

void TH2::FillN(Int_t ntimes, const Double_t *x, const Double_t *y, const Double_t *w, Int_t stride)
{
   Int_t binx, biny, bin, i;
   ntimes *= stride;
   Int_t ifirst = 0;

   //If a buffer is activated, fill buffer
   // (note that this function must not be called from TH2::BufferEmpty)
   if (fBuffer) {
      for (i=0;i<ntimes;i+=stride) {
         if (!fBuffer) break; // buffer can be deleted in BufferFill when is empty
         if (w) BufferFill(x[i],y[i],w[i]);
         else BufferFill(x[i], y[i], 1.);
      }
      // fill the remaining entries if the buffer has been deleted
      if (i < ntimes && fBuffer==0)
         ifirst = i;
      else
         return;
   }

   Double_t ww = 1;
   for (i=ifirst;i<ntimes;i+=stride) {
      fEntries++;
      binx = fXaxis.FindBin(x[i]);
      biny = fYaxis.FindBin(y[i]);
      if (binx <0 || biny <0) continue;
      bin  = biny*(fXaxis.GetNbins()+2) + binx;
      if (w) ww = w[i];
      if (!fSumw2.fN && ww != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();
      if (fSumw2.fN) fSumw2.fArray[bin] += ww*ww;
      AddBinContent(bin,ww);
      if (binx == 0 || binx > fXaxis.GetNbins()) {
         if (!GetStatOverflowsBehaviour()) continue;
      }
      if (biny == 0 || biny > fYaxis.GetNbins()) {
         if (!GetStatOverflowsBehaviour()) continue;
      }
      Double_t z= ww; //(ww > 0 ? ww : -ww);
      fTsumw   += z;
      fTsumw2  += z*z;
      fTsumwx  += z*x[i];
      fTsumwx2 += z*x[i]*x[i];
      fTsumwy  += z*y[i];
      fTsumwy2 += z*y[i]*y[i];
      fTsumwxy += z*x[i]*y[i];
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Fill histogram following distribution in function fname.
///
///  @param fname  : Function name used for filling the historam
///  @param ntimes : number of times the histogram is filled
///  @param rng    : (optional) Random number generator used to sample
///
///   The distribution contained in the function fname (TF2) is integrated
///   over the channel contents.
///   It is normalized to 1.
///   Getting one random number implies:
///     - Generating a random number between 0 and 1 (say r1)
///     - Look in which bin in the normalized integral r1 corresponds to
///     - Fill histogram channel
///   ntimes random numbers are generated
///
///  One can also call TF2::GetRandom2 to get a random variate from a function.

void TH2::FillRandom(const char *fname, Int_t ntimes, TRandom * rng)
{
   Int_t bin, binx, biny, ibin, loop;
   Double_t r1, x, y;
   //*-*- Search for fname in the list of ROOT defined functions
   TObject *fobj = gROOT->GetFunction(fname);
   if (!fobj) { Error("FillRandom", "Unknown function: %s",fname); return; }
   TF2 * f1 = dynamic_cast<TF2*>(fobj);
   if (!f1) { Error("FillRandom", "Function: %s is not a TF2, is a %s",fname,fobj->IsA()->GetName()); return; }


   TAxis & xAxis = fXaxis;
   TAxis & yAxis = fYaxis;

   // in case axes of histogram are not defined use the function axis
   if (fXaxis.GetXmax() <= fXaxis.GetXmin()  || fYaxis.GetXmax() <= fYaxis.GetXmin()) {
      Double_t xmin,xmax,ymin,ymax;
      f1->GetRange(xmin,ymin,xmax,ymax);
      Info("FillRandom","Using function axis and range ([%g,%g],[%g,%g])",xmin, xmax,ymin,ymax);
      xAxis = *(f1->GetHistogram()->GetXaxis());
      yAxis = *(f1->GetHistogram()->GetYaxis());
   }


   // Allocate temporary space to store the integral and compute integral
   Int_t nbinsx = xAxis.GetNbins();
   Int_t nbinsy = yAxis.GetNbins();
   Int_t nbins  = nbinsx*nbinsy;


   Double_t *integral = new Double_t[nbins+1];
   ibin = 0;
   integral[ibin] = 0;
   for (biny=1;biny<=nbinsy;biny++) {
      for (binx=1;binx<=nbinsx;binx++) {
         ibin++;
         Double_t fint = f1->Integral(xAxis.GetBinLowEdge(binx), xAxis.GetBinUpEdge(binx), yAxis.GetBinLowEdge(biny), yAxis.GetBinUpEdge(biny));
         integral[ibin] = integral[ibin-1] + fint;
      }
   }

   // Normalize integral to 1
   if (integral[nbins] == 0 ) {
      delete [] integral;
      Error("FillRandom", "Integral = zero"); return;
   }
   for (bin=1;bin<=nbins;bin++)  integral[bin] /= integral[nbins];

   // Start main loop ntimes
   for (loop=0;loop<ntimes;loop++) {
      r1 = (rng) ? rng->Rndm() : gRandom->Rndm();
      ibin = TMath::BinarySearch(nbins,&integral[0],r1);
      biny = ibin/nbinsx;
      binx = 1 + ibin - nbinsx*biny;
      biny++;
      x    = xAxis.GetBinCenter(binx);
      y    = yAxis.GetBinCenter(biny);
      Fill(x,y);
   }
   delete [] integral;
}


////////////////////////////////////////////////////////////////////////////////
/// Fill histogram following distribution in histogram h.
///
///  @param h      : Histogram  pointer used for smpling random number
///  @param ntimes : number of times the histogram is filled
///  @param rng    : (optional) Random number generator used for sampling
///
///   The distribution contained in the histogram h (TH2) is integrated
///   over the channel contents.
///   It is normalized to 1.
///   Getting one random number implies:
///     - Generating a random number between 0 and 1 (say r1)
///     - Look in which bin in the normalized integral r1 corresponds to
///     - Fill histogram channel
///   ntimes random numbers are generated

void TH2::FillRandom(TH1 *h, Int_t ntimes, TRandom * rng)
{
   if (!h) { Error("FillRandom", "Null histogram"); return; }
   if (fDimension != h->GetDimension()) {
      Error("FillRandom", "Histograms with different dimensions"); return;
   }

   if (h->ComputeIntegral() == 0) return;

   Int_t loop;
   Double_t x,y;
   TH2 *h2 = (TH2*)h;
   for (loop=0;loop<ntimes;loop++) {
      h2->GetRandom2(x,y,rng);
      Fill(x,y);
   }
}


////////////////////////////////////////////////////////////////////////////////

void TH2::DoFitSlices(bool onX,
                      TF1 *f1, Int_t firstbin, Int_t lastbin, Int_t cut, Option_t *option, TObjArray* arr)
{
   TAxis& outerAxis = (onX ? fYaxis : fXaxis);
   TAxis& innerAxis = (onX ? fXaxis : fYaxis);

   Int_t nbins  = outerAxis.GetNbins();
   // get correct first last bins for outer axis
   // when using default values (0,-1) check if an axis range is set in outer axis
   // do same as in DoProjection for inner axis
   if ( lastbin < firstbin && outerAxis.TestBit(TAxis::kAxisRange) ) {
      firstbin = outerAxis.GetFirst();
      lastbin = outerAxis.GetLast();
      // For special case of TAxis::SetRange, when first == 1 and last
      // = N and the range bit has been set, the TAxis will return 0
      // for both.
      if (firstbin == 0 && lastbin == 0)  {
         firstbin = 1;
         lastbin = nbins;
      }
   }
   if (firstbin < 0) firstbin = 0;
   if (lastbin < 0 || lastbin > nbins + 1) lastbin = nbins + 1;
   if (lastbin < firstbin) {firstbin = 0; lastbin = nbins + 1;}


   TString opt = option;
   TString proj_opt = "e";
   Int_t i1 = opt.Index("[");
   Int_t i2 = opt.Index("]");
   if (i1>=0 && i2>i1) {
      proj_opt += opt(i1,i2-i1+1);
      opt.Remove(i1, i2-i1+1);
   }
   opt.ToLower();
   Int_t ngroup = 1;
   if (opt.Contains("g2")) {ngroup = 2; opt.ReplaceAll("g2","");}
   if (opt.Contains("g3")) {ngroup = 3; opt.ReplaceAll("g3","");}
   if (opt.Contains("g4")) {ngroup = 4; opt.ReplaceAll("g4","");}
   if (opt.Contains("g5")) {ngroup = 5; opt.ReplaceAll("g5","");}

   // implement option S sliding merge for each bin using in conjunction with a given Gn
   Int_t nstep = ngroup;
   if (opt.Contains("s"))  nstep = 1;

   //default is to fit with a gaussian
   if (f1 == 0) {
      f1 = (TF1*)gROOT->GetFunction("gaus");
      if (f1 == 0) f1 = new TF1("gaus","gaus",innerAxis.GetXmin(),innerAxis.GetXmax());
      else         f1->SetRange(innerAxis.GetXmin(),innerAxis.GetXmax());
   }
   Int_t npar = f1->GetNpar();
   if (npar <= 0) return;
   Double_t *parsave = new Double_t[npar];
   f1->GetParameters(parsave);

   if (arr) {
      arr->SetOwner();
      arr->Expand(npar + 1);
   }

   //Create one histogram for each function parameter
   Int_t ipar;
   TH1D **hlist = new TH1D*[npar];
   char *name   = new char[2000];
   char *title  = new char[2000];
   const TArrayD *bins = outerAxis.GetXbins();
   // outer axis boudaries used for creating reported histograms are different
   // than the limits used in the projection loop (firstbin,lastbin)
   Int_t firstOutBin = outerAxis.TestBit(TAxis::kAxisRange) ? std::max(firstbin,1) : 1;
   Int_t lastOutBin = outerAxis.TestBit(TAxis::kAxisRange) ?  std::min(lastbin,outerAxis.GetNbins() ) : outerAxis.GetNbins();
   Int_t nOutBins = lastOutBin-firstOutBin+1;
   // merge bins if use nstep > 1 and fixed bins
   if (bins->fN == 0) nOutBins /= nstep;
   for (ipar=0;ipar<npar;ipar++) {
      snprintf(name,2000,"%s_%d",GetName(),ipar);
      snprintf(title,2000,"Fitted value of par[%d]=%s",ipar,f1->GetParName(ipar));
      delete gDirectory->FindObject(name);
      if (bins->fN == 0) {
         hlist[ipar] = new TH1D(name,title, nOutBins, outerAxis.GetBinLowEdge(firstOutBin), outerAxis.GetBinUpEdge(lastOutBin));
      } else {
         hlist[ipar] = new TH1D(name,title, nOutBins, &bins->fArray[firstOutBin-1]);
      }
      hlist[ipar]->GetXaxis()->SetTitle(outerAxis.GetTitle());
      if (arr)
         (*arr)[ipar] = hlist[ipar];
   }
   snprintf(name,2000,"%s_chi2",GetName());
   delete gDirectory->FindObject(name);
   TH1D *hchi2 = 0;
   if (bins->fN == 0) {
      hchi2 = new TH1D(name,"chisquare", nOutBins, outerAxis.GetBinLowEdge(firstOutBin), outerAxis.GetBinUpEdge(lastOutBin));
   } else {
      hchi2 = new TH1D(name,"chisquare", nOutBins, &bins->fArray[firstOutBin-1]);
   }
   hchi2->GetXaxis()->SetTitle(outerAxis.GetTitle());
   if (arr)
      (*arr)[npar] = hchi2;

   //Loop on all bins in Y, generate a projection along X
   Int_t bin;
   // in case of sliding merge nstep=1, i.e. do slices starting for every bin
   // now do not slices case with overflow (makes more sense)
   // when fitting add the option "N". We don;t want to display and store the function
   // for the temporary histograms that are created and fitted
   opt += " n ";
   TH1D *hp = nullptr;
   for (bin=firstbin;bin+ngroup-1<=lastbin;bin += nstep) {
      if (onX)
         hp= ProjectionX("_temp",bin,bin+ngroup-1,proj_opt);
      else
         hp= ProjectionY("_temp",bin,bin+ngroup-1,proj_opt);
      if (hp == 0) continue;
      // nentries can be the effective entries and it could be a very small number but not zero!
      Double_t nentries = hp->GetEntries();
      if ( nentries <= 0 || nentries < cut) {
         if (!opt.Contains("Q"))
               Info("DoFitSlices","Slice %d skipped, the number of entries is zero or smaller than the given cut value, n=%f",bin,nentries);
         continue;
      }
      f1->SetParameters(parsave);
      Int_t binOn = hlist[0]->FindBin(outerAxis.GetBinCenter(bin+ngroup/2));
      if (!opt.Contains("Q"))
         Info("DoFitSlices","Slice fit %d (%f,%f)",binOn,hlist[0]->GetXaxis()->GetBinLowEdge(binOn),hlist[0]->GetXaxis()->GetBinUpEdge(binOn));
      hp->Fit(f1,opt.Data());
      Int_t npfits = f1->GetNumberFitPoints();
      if (npfits > npar && npfits >= cut) {
         for (ipar=0;ipar<npar;ipar++) {
            hlist[ipar]->SetBinContent(binOn,f1->GetParameter(ipar));
            hlist[ipar]->SetBinError(binOn,f1->GetParError(ipar));
         }
         hchi2->SetBinContent(binOn,f1->GetChisquare()/(npfits-npar));
      }
      else {
         if (!opt.Contains("Q"))
            Info("DoFitSlices","Fitted slice %d skipped, the number of fitted points is too small, n=%d",bin,npfits);
      }
      // don't need to delete hp. If histogram has the same name it is re-used in TH2::Projection
   }
   delete hp;
   delete [] parsave;
   delete [] name;
   delete [] title;
   delete [] hlist;
}


////////////////////////////////////////////////////////////////////////////////
/// Project slices along X in case of a 2-D histogram, then fit each slice
/// with function f1 and make a histogram for each fit parameter
/// Only bins along Y between firstybin and lastybin are considered.
/// By default (firstybin == 0, lastybin == -1), all bins in y including
/// over- and underflows are taken into account.
/// If f1=0, a gaussian is assumed
/// Before invoking this function, one can set a subrange to be fitted along X
/// via f1->SetRange(xmin,xmax)
/// The argument option (default="QNR") can be used to change the fit options.
///    - "Q"  means Quiet mode
///    - "N"  means do not show the result of the fit
///    - "R"  means fit the function in the specified function range
///    - "G2" merge 2 consecutive bins along X
///    - "G3" merge 3 consecutive bins along X
///    - "G4" merge 4 consecutive bins along X
///    - "G5" merge 5 consecutive bins along X
///    - "S"  sliding merge: merge n consecutive bins along X accordingly to what Gn is given.
///           It makes sense when used together with a Gn option
///
/// The generated histograms are returned by adding them to arr, if arr is not NULL.
/// arr's SetOwner() is called, to signal that it is the user's responsibility to
/// delete the histograms, possibly by deleting the array.
/// ~~~ {.cpp}
///    TObjArray aSlices;
///    h2->FitSlicesX(func, 0, -1, 0, "QNR", &aSlices);
/// ~~~
/// will already delete the histograms once aSlice goes out of scope. aSlices will
/// contain the histogram for the i-th parameter of the fit function at aSlices[i];
/// aSlices[n] (n being the number of parameters) contains the chi2 distribution of
/// the fits.
///
/// If arr is NULL, the generated histograms are added to the list of objects
/// in the current directory. It is the user's responsibility to delete
/// these histograms.
///
///  Example: Assume a 2-d histogram h2
/// ~~~ {.cpp}
///   Root > h2->FitSlicesX(); produces 4 TH1D histograms
///          with h2_0 containing parameter 0(Constant) for a Gaus fit
///                    of each bin in Y projected along X
///          with h2_1 containing parameter 1(Mean) for a gaus fit
///          with h2_2 containing parameter 2(StdDev)  for a gaus fit
///          with h2_chi2 containing the chisquare/number of degrees of freedom for a gaus fit
///
///   Root > h2->FitSlicesX(0,15,22,10);
///          same as above, but only for bins 15 to 22 along Y
///          and only for bins in Y for which the corresponding projection
///          along X has more than cut bins filled.
/// ~~~
///  NOTE: To access the generated histograms in the current directory, do eg:
/// ~~~ {.cpp}
///     TH1D *h2_1 = (TH1D*)gDirectory->Get("h2_1");
/// ~~~

void TH2::FitSlicesX(TF1 *f1, Int_t firstybin, Int_t lastybin, Int_t cut, Option_t *option, TObjArray* arr)
{
   DoFitSlices(true, f1, firstybin, lastybin, cut, option, arr);

}


////////////////////////////////////////////////////////////////////////////////
/// Project slices along Y in case of a 2-D histogram, then fit each slice
/// with function f1 and make a histogram for each fit parameter
/// Only bins along X between firstxbin and lastxbin are considered.
/// By default (firstxbin == 0, lastxbin == -1), all bins in x including
/// over- and underflows are taken into account.
/// If f1=0, a gaussian is assumed
/// Before invoking this function, one can set a subrange to be fitted along Y
/// via f1->SetRange(ymin,ymax)
/// The argument option (default="QNR") can be used to change the fit options.
///    - "Q"  means Quiet mode
///    - "N"  means do not show the result of the fit
///    - "R"  means fit the function in the specified function range
///    - "G2" merge 2 consecutive bins along Y
///    - "G3" merge 3 consecutive bins along Y
///    - "G4" merge 4 consecutive bins along Y
///    - "G5" merge 5 consecutive bins along Y
///    - "S"  sliding merge: merge n consecutive bins along Y accordingly to what Gn is given.
///           It makes sense when used together with a Gn option
///
/// The generated histograms are returned by adding them to arr, if arr is not NULL.
/// arr's SetOwner() is called, to signal that it is the user's responsibility to
/// delete the histograms, possibly by deleting the array.
/// ~~~ {.cpp}
///    TObjArray aSlices;
///    h2->FitSlicesY(func, 0, -1, 0, "QNR", &aSlices);
/// ~~~
/// will already delete the histograms once aSlice goes out of scope. aSlices will
/// contain the histogram for the i-th parameter of the fit function at aSlices[i];
/// aSlices[n] (n being the number of parameters) contains the chi2 distribution of
/// the fits.
///
/// If arr is NULL, the generated histograms are added to the list of objects
/// in the current directory. It is the user's responsibility to delete
/// these histograms.
///
///  Example: Assume a 2-d histogram h2
/// ~~~ {.cpp}
///   Root > h2->FitSlicesY(); produces 4 TH1D histograms
///          with h2_0 containing parameter 0(Constant) for a Gaus fit
///                    of each bin in X projected along Y
///          with h2_1 containing parameter 1(Mean) for a gaus fit
///          with h2_2 containing parameter 2(StdDev)  for a gaus fit
///          with h2_chi2 containing the chisquare/number of degrees of freedom for a gaus fit
///
///   Root > h2->FitSlicesY(0,15,22,10);
///          same as above, but only for bins 15 to 22 along X
///          and only for bins in X for which the corresponding projection
///          along Y has more than cut bins filled.
/// ~~~
///
///  NOTE: To access the generated histograms in the current directory, do eg:
/// ~~~ {.cpp}
///     TH1D *h2_1 = (TH1D*)gDirectory->Get("h2_1");
/// ~~~
///
/// A complete example of this function is given in tutorial:fitslicesy.C.

void TH2::FitSlicesY(TF1 *f1, Int_t firstxbin, Int_t lastxbin, Int_t cut, Option_t *option, TObjArray* arr)
{
   DoFitSlices(false, f1, firstxbin, lastxbin, cut, option, arr);
}

Int_t TH2::GetBin(Int_t binx, Int_t biny, Int_t) const
{
   // See comments in TH1::GetBin
   Int_t ofy = fYaxis.GetNbins() + 1; // overflow bin
   if (biny < 0) biny = 0;
   if (biny > ofy) biny = ofy;

   return TH1::GetBin(binx) + (fXaxis.GetNbins() + 2) * biny;
}


////////////////////////////////////////////////////////////////////////////////
/// compute first cell (binx,biny) in the range [firstxbin,lastxbin][firstybin,lastybin] for which
/// diff = abs(cell_content-c) <= maxdiff
/// In case several cells in the specified range with diff=0 are found
/// the first cell found is returned in binx,biny.
/// In case several cells in the specified range satisfy diff <=maxdiff
/// the cell with the smallest difference is returned in binx,biny.
/// In all cases the function returns the smallest difference.
///
/// NOTE1: if firstxbin < 0, firstxbin is set to 1
///        if (lastxbin < firstxbin then lastxbin is set to the number of bins in X
///          ie if firstxbin=1 and lastxbin=0 (default) the search is on all bins in X except
///          for X's under- and overflow bins.
///        if firstybin < 0, firstybin is set to 1
///        if (lastybin < firstybin then lastybin is set to the number of bins in Y
///          ie if firstybin=1 and lastybin=0 (default) the search is on all bins in Y except
///          for Y's under- and overflow bins.
///
/// NOTE2: if maxdiff=0 (default), the first cell with content=c is returned.

Double_t TH2::GetBinWithContent2(Double_t c, Int_t &binx, Int_t &biny, Int_t firstxbin, Int_t lastxbin,
                                 Int_t firstybin, Int_t lastybin, Double_t maxdiff) const
{
   if (fDimension != 2) {
      binx = -1;
      biny = -1;
      Error("GetBinWithContent2","function is only valid for 2-D histograms");
      return 0;
   }
   if (firstxbin < 0) firstxbin = 1;
   if (lastxbin < firstxbin) lastxbin = fXaxis.GetNbins();
   if (firstybin < 0) firstybin = 1;
   if (lastybin < firstybin) lastybin = fYaxis.GetNbins();
   Double_t diff, curmax = 1.e240;
   for (Int_t j = firstybin; j <= lastybin; j++) {
      for (Int_t i = firstxbin; i <= lastxbin; i++) {
         diff = TMath::Abs(GetBinContent(i,j)-c);
         if (diff <= 0) {binx = i; biny=j; return diff;}
         if (diff < curmax && diff <= maxdiff) {curmax = diff, binx=i; biny=j;}
      }
   }
   return curmax;
}


////////////////////////////////////////////////////////////////////////////////
/// Return correlation factor between axis1 and axis2.

Double_t TH2::GetCorrelationFactor(Int_t axis1, Int_t axis2) const
{
   if (axis1 < 1 || axis2 < 1 || axis1 > 2 || axis2 > 2) {
      Error("GetCorrelationFactor","Wrong parameters");
      return 0;
   }
   if (axis1 == axis2) return 1;
   Double_t stddev1 = GetStdDev(axis1);
   if (stddev1 == 0) return 0;
   Double_t stddev2 = GetStdDev(axis2);
   if (stddev2 == 0) return 0;
   return GetCovariance(axis1,axis2)/stddev1/stddev2;
}


////////////////////////////////////////////////////////////////////////////////
/// Return covariance between axis1 and axis2.

Double_t TH2::GetCovariance(Int_t axis1, Int_t axis2) const
{
   if (axis1 < 1 || axis2 < 1 || axis1 > 2 || axis2 > 2) {
      Error("GetCovariance","Wrong parameters");
      return 0;
   }
   Double_t stats[kNstat];
   GetStats(stats);
   Double_t sumw   = stats[0];
   //Double_t sumw2  = stats[1];
   Double_t sumwx  = stats[2];
   Double_t sumwx2 = stats[3];
   Double_t sumwy  = stats[4];
   Double_t sumwy2 = stats[5];
   Double_t sumwxy = stats[6];

   if (sumw == 0) return 0;
   if (axis1 == 1 && axis2 == 1) {
      return TMath::Abs(sumwx2/sumw - sumwx/sumw*sumwx/sumw);
   }
   if (axis1 == 2 && axis2 == 2) {
      return TMath::Abs(sumwy2/sumw - sumwy/sumw*sumwy/sumw);
   }
   return sumwxy/sumw - sumwx/sumw*sumwy/sumw;
}


////////////////////////////////////////////////////////////////////////////////
/// Return 2 random numbers along axis x and y distributed according
/// to the cell-contents of this 2-dim histogram
/// return a NaN if the histogram has a bin with negative content
///
/// @param[out] x  reference to random generated x value
/// @param[out] y  reference to random generated x value
/// @param[in] rng (optional) Random number generator pointer used (default is gRandom)

void TH2::GetRandom2(Double_t &x, Double_t &y, TRandom * rng)
{
   Int_t nbinsx = GetNbinsX();
   Int_t nbinsy = GetNbinsY();
   Int_t nbins  = nbinsx*nbinsy;
   Double_t integral;
   // compute integral checking that all bins have positive content (see ROOT-5894)
   if (fIntegral) {
      if (fIntegral[nbins+1] != fEntries) integral = ComputeIntegral(true);
      else integral = fIntegral[nbins];
   } else {
      integral = ComputeIntegral(true);
   }
   if (integral == 0 ) { x = 0; y = 0; return;}
   // case histogram has negative bins
   if (integral == TMath::QuietNaN() ) { x = TMath::QuietNaN(); y = TMath::QuietNaN(); return;}

   if (!rng) rng = gRandom;
   Double_t r1 = rng->Rndm();
   Int_t ibin = TMath::BinarySearch(nbins,fIntegral,(Double_t) r1);
   Int_t biny = ibin/nbinsx;
   Int_t binx = ibin - nbinsx*biny;
   x = fXaxis.GetBinLowEdge(binx+1);
   if (r1 > fIntegral[ibin]) x +=
      fXaxis.GetBinWidth(binx+1)*(r1-fIntegral[ibin])/(fIntegral[ibin+1] - fIntegral[ibin]);
   y = fYaxis.GetBinLowEdge(biny+1) + fYaxis.GetBinWidth(biny+1)*rng->Rndm();
}


////////////////////////////////////////////////////////////////////////////////
/// Fill the array stats from the contents of this histogram
/// The array stats must be correctly dimensioned in the calling program.
/// ~~~ {.cpp}
/// stats[0] = sumw
/// stats[1] = sumw2
/// stats[2] = sumwx
/// stats[3] = sumwx2
/// stats[4] = sumwy
/// stats[5] = sumwy2
/// stats[6] = sumwxy
/// ~~~
///
/// If no axis-subranges are specified (via TAxis::SetRange), the array stats
/// is simply a copy of the statistics quantities computed at filling time.
/// If sub-ranges are specified, the function recomputes these quantities
/// from the bin contents in the current axis ranges.
///
///  Note that the mean value/StdDev is computed using the bins in the currently
///  defined ranges (see TAxis::SetRange). By default the ranges include
///  all bins from 1 to nbins included, excluding underflows and overflows.
///  To force the underflows and overflows in the computation, one must
///  call the static function TH1::StatOverflows(kTRUE) before filling
///  the histogram.

void TH2::GetStats(Double_t *stats) const
{
   if (fBuffer) ((TH2*)this)->BufferEmpty();

   if ((fTsumw == 0 && fEntries > 0) || fXaxis.TestBit(TAxis::kAxisRange) || fYaxis.TestBit(TAxis::kAxisRange)) {
      std::fill(stats, stats + 7, 0);

      Int_t firstBinX = fXaxis.GetFirst();
      Int_t lastBinX  = fXaxis.GetLast();
      Int_t firstBinY = fYaxis.GetFirst();
      Int_t lastBinY  = fYaxis.GetLast();
      // include underflow/overflow if TH1::StatOverflows(kTRUE) in case no range is set on the axis
      if (GetStatOverflowsBehaviour()) {
        if ( !fXaxis.TestBit(TAxis::kAxisRange) ) {
            if (firstBinX == 1) firstBinX = 0;
            if (lastBinX ==  fXaxis.GetNbins() ) lastBinX += 1;
         }
         if ( !fYaxis.TestBit(TAxis::kAxisRange) ) {
            if (firstBinY == 1) firstBinY = 0;
            if (lastBinY ==  fYaxis.GetNbins() ) lastBinY += 1;
         }
      }
      // check for labels axis . In that case corresponsing statistics do not make sense and it is set to zero
      Bool_t labelXaxis =  ((const_cast<TAxis&>(fXaxis)).GetLabels() && fXaxis.CanExtend() );
      Bool_t labelYaxis =  ((const_cast<TAxis&>(fYaxis)).GetLabels() && fYaxis.CanExtend() );

      for (Int_t biny = firstBinY; biny <= lastBinY; ++biny) {
         Double_t y = (!labelYaxis) ? fYaxis.GetBinCenter(biny) : 0;
         for (Int_t binx = firstBinX; binx <= lastBinX; ++binx) {
            Double_t x = (!labelXaxis) ? fXaxis.GetBinCenter(binx) : 0;
            //w   = TMath::Abs(GetBinContent(bin));
            Int_t bin = GetBin(binx,biny);
            Double_t w = RetrieveBinContent(bin);
            Double_t wx = w * x; // avoid some extra multiplications at the expense of some clarity
            Double_t wy = w * y;

            stats[0] += w;
            stats[1] += GetBinErrorSqUnchecked(bin);
            stats[2] += wx;
            stats[3] += wx * x;
            stats[4] += wy;
            stats[5] += wy * y;
            stats[6] += wx * y;
         }
      }
   } else {
      stats[0] = fTsumw;
      stats[1] = fTsumw2;
      stats[2] = fTsumwx;
      stats[3] = fTsumwx2;
      stats[4] = fTsumwy;
      stats[5] = fTsumwy2;
      stats[6] = fTsumwxy;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Return integral of bin contents. Only bins in the bins range are considered.
/// By default the integral is computed as the sum of bin contents in the range.
/// if option "width" is specified, the integral is the sum of
/// the bin contents multiplied by the bin width in x and in y.

Double_t TH2::Integral(Option_t *option) const
{
   return Integral(fXaxis.GetFirst(),fXaxis.GetLast(),
      fYaxis.GetFirst(),fYaxis.GetLast(),option);
}


////////////////////////////////////////////////////////////////////////////////
/// Return integral of bin contents in range [firstxbin,lastxbin],[firstybin,lastybin]
/// for a 2-D histogram
/// By default the integral is computed as the sum of bin contents in the range.
/// if option "width" is specified, the integral is the sum of
/// the bin contents multiplied by the bin width in x and in y.

Double_t TH2::Integral(Int_t firstxbin, Int_t lastxbin, Int_t firstybin, Int_t lastybin, Option_t *option) const
{
   double err = 0;
   return DoIntegral(firstxbin,lastxbin,firstybin,lastybin,-1,0,err,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Return integral of bin contents in range [firstxbin,lastxbin],[firstybin,lastybin]
/// for a 2-D histogram. Calculates also the integral error using error propagation
/// from the bin errors assuming that all the bins are uncorrelated.
/// By default the integral is computed as the sum of bin contents in the range.
/// if option "width" is specified, the integral is the sum of
/// the bin contents multiplied by the bin width in x and in y.

Double_t TH2::IntegralAndError(Int_t firstxbin, Int_t lastxbin, Int_t firstybin, Int_t lastybin, Double_t & error, Option_t *option) const
{
   return DoIntegral(firstxbin,lastxbin,firstybin,lastybin,-1,0,error,option,kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
///illegal for a TH2

Double_t TH2::Interpolate(Double_t) const
{
   Error("Interpolate","This function must be called with 2 arguments for a TH2");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Given a point P(x,y), Interpolate approximates the value via bilinear
/// interpolation based on the four nearest bin centers
/// see Wikipedia, Bilinear Interpolation
/// Andy Mastbaum 10/8/2008
/// vaguely based on R.Raja 6-Sep-2008

 Double_t TH2::Interpolate(Double_t x, Double_t y) const
{
   Double_t f=0;
   Double_t x1=0,x2=0,y1=0,y2=0;
   Double_t dx,dy;
   Int_t bin_x = fXaxis.FindFixBin(x);
   Int_t bin_y = fYaxis.FindFixBin(y);
   if(bin_x<1 || bin_x>GetNbinsX() || bin_y<1 || bin_y>GetNbinsY()) {
      Error("Interpolate","Cannot interpolate outside histogram domain.");
      return 0;
   }
   Int_t quadrant = 0; // CCW from UR 1,2,3,4
   // which quadrant of the bin (bin_P) are we in?
   dx = fXaxis.GetBinUpEdge(bin_x)-x;
   dy = fYaxis.GetBinUpEdge(bin_y)-y;
   if (dx<=fXaxis.GetBinWidth(bin_x)/2 && dy<=fYaxis.GetBinWidth(bin_y)/2)
   quadrant = 1; // upper right
   if (dx>fXaxis.GetBinWidth(bin_x)/2 && dy<=fYaxis.GetBinWidth(bin_y)/2)
   quadrant = 2; // upper left
   if (dx>fXaxis.GetBinWidth(bin_x)/2 && dy>fYaxis.GetBinWidth(bin_y)/2)
   quadrant = 3; // lower left
   if (dx<=fXaxis.GetBinWidth(bin_x)/2 && dy>fYaxis.GetBinWidth(bin_y)/2)
   quadrant = 4; // lower right
   switch(quadrant) {
   case 1:
      x1 = fXaxis.GetBinCenter(bin_x);
      y1 = fYaxis.GetBinCenter(bin_y);
      x2 = fXaxis.GetBinCenter(bin_x+1);
      y2 = fYaxis.GetBinCenter(bin_y+1);
      break;
   case 2:
      x1 = fXaxis.GetBinCenter(bin_x-1);
      y1 = fYaxis.GetBinCenter(bin_y);
      x2 = fXaxis.GetBinCenter(bin_x);
      y2 = fYaxis.GetBinCenter(bin_y+1);
      break;
   case 3:
      x1 = fXaxis.GetBinCenter(bin_x-1);
      y1 = fYaxis.GetBinCenter(bin_y-1);
      x2 = fXaxis.GetBinCenter(bin_x);
      y2 = fYaxis.GetBinCenter(bin_y);
      break;
   case 4:
      x1 = fXaxis.GetBinCenter(bin_x);
      y1 = fYaxis.GetBinCenter(bin_y-1);
      x2 = fXaxis.GetBinCenter(bin_x+1);
      y2 = fYaxis.GetBinCenter(bin_y);
      break;
   }
   Int_t bin_x1 = fXaxis.FindFixBin(x1);
   if(bin_x1<1) bin_x1=1;
   Int_t bin_x2 = fXaxis.FindFixBin(x2);
   if(bin_x2>GetNbinsX()) bin_x2=GetNbinsX();
   Int_t bin_y1 = fYaxis.FindFixBin(y1);
   if(bin_y1<1) bin_y1=1;
   Int_t bin_y2 = fYaxis.FindFixBin(y2);
   if(bin_y2>GetNbinsY()) bin_y2=GetNbinsY();
   Int_t bin_q22 = GetBin(bin_x2,bin_y2);
   Int_t bin_q12 = GetBin(bin_x1,bin_y2);
   Int_t bin_q11 = GetBin(bin_x1,bin_y1);
   Int_t bin_q21 = GetBin(bin_x2,bin_y1);
   Double_t q11 = RetrieveBinContent(bin_q11);
   Double_t q12 = RetrieveBinContent(bin_q12);
   Double_t q21 = RetrieveBinContent(bin_q21);
   Double_t q22 = RetrieveBinContent(bin_q22);
   Double_t d = 1.0*(x2-x1)*(y2-y1);
   f = 1.0*q11/d*(x2-x)*(y2-y)+1.0*q21/d*(x-x1)*(y2-y)+1.0*q12/d*(x2-x)*(y-y1)+1.0*q22/d*(x-x1)*(y-y1);
   return f;
}


////////////////////////////////////////////////////////////////////////////////
///illegal for a TH2

Double_t TH2::Interpolate(Double_t, Double_t, Double_t) const
{
   Error("Interpolate","This function must be called with 2 arguments for a TH2");
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Statistical test of compatibility in shape between
/// THIS histogram and h2, using Kolmogorov test.
///     Default: Ignore under- and overflow bins in comparison
///
///     option is a character string to specify options
///        - "U" include Underflows in test
///        - "O" include Overflows
///        - "N" include comparison of normalizations
///        - "D" Put out a line of "Debug" printout
///        - "M" Return the Maximum Kolmogorov distance instead of prob
///
///   The returned function value is the probability of test
///       (much less than one means NOT compatible)
///
///   The KS test uses the distance between the pseudo-CDF's obtained
///   from the histogram. Since in 2D the order for generating the pseudo-CDF is
///   arbitrary, two pairs of pseudo-CDF are used, one starting from the x axis the
///   other from the y axis and the maximum distance is the average of the two maximum
///   distances obtained.
///
///  Code adapted by Rene Brun from original HBOOK routine HDIFF

Double_t TH2::KolmogorovTest(const TH1 *h2, Option_t *option) const
{
   TString opt = option;
   opt.ToUpper();

   Double_t prb = 0;
   TH1 *h1 = (TH1*)this;
   if (h2 == 0) return 0;
   const TAxis *xaxis1 = h1->GetXaxis();
   const TAxis *xaxis2 = h2->GetXaxis();
   const TAxis *yaxis1 = h1->GetYaxis();
   const TAxis *yaxis2 = h2->GetYaxis();
   Int_t ncx1   = xaxis1->GetNbins();
   Int_t ncx2   = xaxis2->GetNbins();
   Int_t ncy1   = yaxis1->GetNbins();
   Int_t ncy2   = yaxis2->GetNbins();

   // Check consistency of dimensions
   if (h1->GetDimension() != 2 || h2->GetDimension() != 2) {
      Error("KolmogorovTest","Histograms must be 2-D\n");
      return 0;
   }

   // Check consistency in number of channels
   if (ncx1 != ncx2) {
      Error("KolmogorovTest","Number of channels in X is different, %d and %d\n",ncx1,ncx2);
      return 0;
   }
   if (ncy1 != ncy2) {
      Error("KolmogorovTest","Number of channels in Y is different, %d and %d\n",ncy1,ncy2);
      return 0;
   }

   // Check consistency in channel edges
   Bool_t afunc1 = kFALSE;
   Bool_t afunc2 = kFALSE;
   Double_t difprec = 1e-5;
   Double_t diff1 = TMath::Abs(xaxis1->GetXmin() - xaxis2->GetXmin());
   Double_t diff2 = TMath::Abs(xaxis1->GetXmax() - xaxis2->GetXmax());
   if (diff1 > difprec || diff2 > difprec) {
      Error("KolmogorovTest","histograms with different binning along X");
      return 0;
   }
   diff1 = TMath::Abs(yaxis1->GetXmin() - yaxis2->GetXmin());
   diff2 = TMath::Abs(yaxis1->GetXmax() - yaxis2->GetXmax());
   if (diff1 > difprec || diff2 > difprec) {
      Error("KolmogorovTest","histograms with different binning along Y");
      return 0;
   }

   //   Should we include Uflows, Oflows?
   Int_t ibeg = 1, jbeg = 1;
   Int_t iend = ncx1, jend = ncy1;
   if (opt.Contains("U")) {ibeg = 0; jbeg = 0;}
   if (opt.Contains("O")) {iend = ncx1+1; jend = ncy1+1;}

   Int_t i,j;
   Double_t sum1  = 0;
   Double_t sum2  = 0;
   Double_t w1    = 0;
   Double_t w2    = 0;
   for (i = ibeg; i <= iend; i++) {
      for (j = jbeg; j <= jend; j++) {
         sum1 += h1->GetBinContent(i,j);
         sum2 += h2->GetBinContent(i,j);
         Double_t ew1   = h1->GetBinError(i,j);
         Double_t ew2   = h2->GetBinError(i,j);
         w1   += ew1*ew1;
         w2   += ew2*ew2;

      }
   }

   //    Check that both scatterplots contain events
   if (sum1 == 0) {
      Error("KolmogorovTest","Integral is zero for h1=%s\n",h1->GetName());
      return 0;
   }
   if (sum2 == 0) {
      Error("KolmogorovTest","Integral is zero for h2=%s\n",h2->GetName());
      return 0;
   }
   // calculate the effective entries.
   // the case when errors are zero (w1 == 0 or w2 ==0) are equivalent to
   // compare to a function. In that case the rescaling is done only on sqrt(esum2) or sqrt(esum1)
   Double_t esum1 = 0, esum2 = 0;
   if (w1 > 0)
      esum1 = sum1 * sum1 / w1;
   else
      afunc1 = kTRUE;    // use later for calculating z

   if (w2 > 0)
      esum2 = sum2 * sum2 / w2;
   else
      afunc2 = kTRUE;    // use later for calculating z

   if (afunc2 && afunc1) {
      Error("KolmogorovTest","Errors are zero for both histograms\n");
      return 0;
   }

   //   Find first Kolmogorov distance
   Double_t s1 = 1/sum1;
   Double_t s2 = 1/sum2;
   Double_t dfmax1 = 0;
   Double_t rsum1=0, rsum2=0;
   for (i=ibeg;i<=iend;i++) {
      for (j=jbeg;j<=jend;j++) {
         rsum1 += s1*h1->GetBinContent(i,j);
         rsum2 += s2*h2->GetBinContent(i,j);
         dfmax1  = TMath::Max(dfmax1, TMath::Abs(rsum1-rsum2));
      }
   }

   //   Find second Kolmogorov distance
   Double_t dfmax2 = 0;
   rsum1=0, rsum2=0;
   for (j=jbeg;j<=jend;j++) {
      for (i=ibeg;i<=iend;i++) {
         rsum1 += s1*h1->GetBinContent(i,j);
         rsum2 += s2*h2->GetBinContent(i,j);
         dfmax2 = TMath::Max(dfmax2, TMath::Abs(rsum1-rsum2));
      }
   }

   //    Get Kolmogorov probability: use effective entries, esum1 or esum2,  for normalizing it
   Double_t factnm;
   if (afunc1)      factnm = TMath::Sqrt(esum2);
   else if (afunc2) factnm = TMath::Sqrt(esum1);
   else             factnm = TMath::Sqrt(esum1*sum2/(esum1+esum2));

   // take average of the two distances
   Double_t dfmax = 0.5*(dfmax1+dfmax2);
   Double_t z  = dfmax*factnm;

   prb = TMath::KolmogorovProb(z);

   Double_t prb1 = 0, prb2 = 0;
   // option N to combine normalization makes sense if both afunc1 and afunc2 are false
   if (opt.Contains("N")  && !(afunc1 || afunc2 ) ) {
      // Combine probabilities for shape and normalization
      prb1   = prb;
      Double_t d12    = esum1-esum2;
      Double_t chi2   = d12*d12/(esum1+esum2);
      prb2   = TMath::Prob(chi2,1);
      //     see Eadie et al., section 11.6.2
      if (prb > 0 && prb2 > 0) prb = prb*prb2*(1-TMath::Log(prb*prb2));
      else                     prb = 0;
   }

   //    debug printout
   if (opt.Contains("D")) {
      printf(" Kolmo Prob  h1 = %s, sum1=%g\n",h1->GetName(),sum1);
      printf(" Kolmo Prob  h2 = %s, sum2=%g\n",h2->GetName(),sum2);
      printf(" Kolmo Probabil = %f, Max Dist = %g\n",prb,dfmax);
      if (opt.Contains("N"))
         printf(" Kolmo Probabil = %f for shape alone, =%f for normalisation alone\n",prb1,prb2);
   }
   // This numerical error condition should never occur:
   if (TMath::Abs(rsum1-1) > 0.002) Warning("KolmogorovTest","Numerical problems with h1=%s\n",h1->GetName());
   if (TMath::Abs(rsum2-1) > 0.002) Warning("KolmogorovTest","Numerical problems with h2=%s\n",h2->GetName());

   if(opt.Contains("M"))      return dfmax;  // return average of max distance

   return prb;
}


////////////////////////////////////////////////////////////////////////////////
/// Rebin only the X axis
/// see Rebin2D

TH2 *TH2::RebinX(Int_t ngroup, const char *newname)
{
   return Rebin2D(ngroup, 1, newname);
}


////////////////////////////////////////////////////////////////////////////////
/// Rebin only the Y axis
/// see Rebin2D

TH2 *TH2::RebinY(Int_t ngroup, const char *newname)
{
   return Rebin2D(1, ngroup, newname);
}

////////////////////////////////////////////////////////////////////////////////
/// Override TH1::Rebin as TH2::RebinX
/// Rebinning in variable binning as for TH1 is not allowed
/// If a non-null pointer is given an error is flagged
/// see RebinX and Rebin2D

TH2 * TH2::Rebin( Int_t ngroup, const char*newname, const Double_t *xbins)
{
   if (xbins != nullptr) {
      Error("Rebin","Rebinning a 2-d histogram into variable bins is not supported (it is possible only for 1-d histograms). Return a nullptr");
      return nullptr;
   }
   Info("Rebin","Rebinning only the x-axis. Use Rebin2D for rebinning both axes");
   return RebinX(ngroup, newname);
}
////////////////////////////////////////////////////////////////////////////////
/// Rebin this histogram grouping nxgroup/nygroup bins along the xaxis/yaxis together.
///
///   if newname is not blank a new temporary histogram hnew is created.
///   else the current histogram is modified (default)
///   The parameter nxgroup/nygroup indicate how many bins along the xaxis/yaxis of this
///   have to me merged into one bin of hnew
///   If the original histogram has errors stored (via Sumw2), the resulting
///   histograms has new errors correctly calculated.
///
///   examples: if hpxpy is an existing TH2 histogram with 40 x 40 bins
/// ~~~ {.cpp}
///     hpxpy->Rebin2D();  // merges two bins along the xaxis and yaxis in one in hpxpy
///                        // Carefull: previous contents of hpxpy are lost
///     hpxpy->RebinX(5);  //merges five bins along the xaxis in one in hpxpy
///     TH2 *hnew = hpxpy->RebinY(5,"hnew"); // creates a new histogram hnew
///                                          // merging 5 bins of h1 along the yaxis in one bin
/// ~~~
///
///   NOTE : If nxgroup/nygroup is not an exact divider of the number of bins,
///          along the xaxis/yaxis the top limit(s) of the rebinned histogram
///          is changed to the upper edge of the xbin=newxbins*nxgroup resp.
///          ybin=newybins*nygroup and the corresponding bins are added to
///          the overflow bin.
///          Statistics will be recomputed from the new bin contents.

TH2 *TH2::Rebin2D(Int_t nxgroup, Int_t nygroup, const char *newname)
{
   Int_t nxbins  = fXaxis.GetNbins();
   Int_t nybins  = fYaxis.GetNbins();
   Int_t nx      = nxbins + 2; // normal bins + underflow and overflow
   Int_t ny      = nybins + 2;
   Double_t xmin  = fXaxis.GetXmin();
   Double_t xmax  = fXaxis.GetXmax();
   Double_t ymin  = fYaxis.GetXmin();
   Double_t ymax  = fYaxis.GetXmax();

   if (GetDimension() != 2) {
      Error("Rebin2D", "Histogram must be TH2. This histogram has %d dimensions.", GetDimension());
      return 0;
   }
   if ((nxgroup <= 0) || (nxgroup > nxbins)) {
      Error("Rebin2D", "Illegal value of nxgroup=%d",nxgroup);
      return 0;
   }
   if ((nygroup <= 0) || (nygroup > nybins)) {
      Error("Rebin2D", "Illegal value of nygroup=%d",nygroup);
      return 0;
   }

   Int_t newxbins = nxbins / nxgroup;
   Int_t newybins = nybins / nygroup;
   Int_t newnx = newxbins + 2; // regular bins + overflow / underflow
   Int_t newny = newybins + 2; // regular bins + overflow / underflow

   // Save old bin contents into a new array
   Double_t *oldBins = new Double_t[fNcells];
   for (Int_t i = 0; i < fNcells; ++i) oldBins[i] = RetrieveBinContent(i);

   Double_t* oldErrors = NULL;
   if (fSumw2.fN) {
      oldErrors = new Double_t[fNcells];
      for (Int_t i = 0; i < fNcells; ++i) oldErrors[i] = GetBinErrorSqUnchecked(i);
   }

   // create a clone of the old histogram if newname is specified
   TH2* hnew = this;
   if (newname && strlen(newname)) {
      hnew = (TH2*)Clone();
      hnew->SetName(newname);
   }

   bool resetStat = false;

   // change axis specs and rebuild bin contents array
   if(newxbins * nxgroup != nxbins) {
      xmax = fXaxis.GetBinUpEdge(newxbins * nxgroup);
      resetStat = true; // stats must be reset because top bins will be moved to overflow bin
   }
   if(newybins * nygroup != nybins) {
      ymax = fYaxis.GetBinUpEdge(newybins * nygroup);
      resetStat = true; // stats must be reset because top bins will be moved to overflow bin
   }

   // save the TAttAxis members (reset by SetBins) for x axis
   Int_t    nXdivisions  = fXaxis.GetNdivisions();
   Color_t  xAxisColor   = fXaxis.GetAxisColor();
   Color_t  xLabelColor  = fXaxis.GetLabelColor();
   Style_t  xLabelFont   = fXaxis.GetLabelFont();
   Float_t  xLabelOffset = fXaxis.GetLabelOffset();
   Float_t  xLabelSize   = fXaxis.GetLabelSize();
   Float_t  xTickLength  = fXaxis.GetTickLength();
   Float_t  xTitleOffset = fXaxis.GetTitleOffset();
   Float_t  xTitleSize   = fXaxis.GetTitleSize();
   Color_t  xTitleColor  = fXaxis.GetTitleColor();
   Style_t  xTitleFont   = fXaxis.GetTitleFont();
   // save the TAttAxis members (reset by SetBins) for y axis
   Int_t    nYdivisions  = fYaxis.GetNdivisions();
   Color_t  yAxisColor   = fYaxis.GetAxisColor();
   Color_t  yLabelColor  = fYaxis.GetLabelColor();
   Style_t  yLabelFont   = fYaxis.GetLabelFont();
   Float_t  yLabelOffset = fYaxis.GetLabelOffset();
   Float_t  yLabelSize   = fYaxis.GetLabelSize();
   Float_t  yTickLength  = fYaxis.GetTickLength();
   Float_t  yTitleOffset = fYaxis.GetTitleOffset();
   Float_t  yTitleSize   = fYaxis.GetTitleSize();
   Color_t  yTitleColor  = fYaxis.GetTitleColor();
   Style_t  yTitleFont   = fYaxis.GetTitleFont();


   // copy merged bin contents (ignore under/overflows)
   if (nxgroup != 1 || nygroup != 1) {
      if(fXaxis.GetXbins()->GetSize() > 0 || fYaxis.GetXbins()->GetSize() > 0){
         // variable bin sizes in x or y, don't treat both cases separately
         Double_t *xbins = new Double_t[newxbins + 1];
         for(Int_t i = 0; i <= newxbins; ++i) xbins[i] = fXaxis.GetBinLowEdge(1 + i * nxgroup);
         Double_t *ybins = new Double_t[newybins + 1];
         for(Int_t i = 0; i <= newybins; ++i) ybins[i] = fYaxis.GetBinLowEdge(1 + i * nygroup);
         hnew->SetBins(newxbins, xbins, newybins, ybins); // changes also errors array (if any)
         delete [] xbins;
         delete [] ybins;
      } else {
         hnew->SetBins(newxbins, xmin, xmax, newybins, ymin, ymax); //changes also errors array
      }

      // (0, 0): x - underflow; y - underflow
      hnew->UpdateBinContent(0, oldBins[0]);
      if (oldErrors) hnew->fSumw2[0] = 0;

      // (x, 0): x - regular / overflow; y - underflow
      for(Int_t binx = 1, oldbinx = 1; binx < newnx; ++binx, oldbinx += nxgroup){
         Double_t binContent = 0.0, binErrorSq = 0.0;
         for (Int_t i = 0; i < nxgroup && (oldbinx + i) < nx; ++i) {
            Int_t bin = oldbinx + i;
            binContent += oldBins[bin];
            if(oldErrors) binErrorSq += oldErrors[bin];
         }
         Int_t newbin = binx;
         hnew->UpdateBinContent(newbin, binContent);
         if (oldErrors) hnew->fSumw2[newbin] = binErrorSq;
      }

      // (0, y): x - underflow; y - regular / overflow
      for(Int_t biny = 1, oldbiny = 1; biny < newny; ++biny, oldbiny += nygroup){
         Double_t binContent = 0.0, binErrorSq = 0.0;
         for (Int_t j = 0; j < nygroup && (oldbiny + j) < ny; ++j) {
            Int_t bin = (oldbiny + j) * nx;
            binContent += oldBins[bin];
            if(oldErrors) binErrorSq += oldErrors[bin];
         }
         Int_t newbin = biny * newnx;
         hnew->UpdateBinContent(newbin, binContent);
         if (oldErrors) hnew->fSumw2[newbin] = binErrorSq;
      }

      // (x, y): x - regular / overflow; y - regular / overflow
      for (Int_t binx = 1, oldbinx = 1; binx < newnx; ++binx, oldbinx += nxgroup) {
         for (Int_t biny = 1, oldbiny = 1; biny < newny; ++biny, oldbiny += nygroup) {
            Double_t binContent = 0.0, binErrorSq = 0.0;
            for (Int_t i = 0; i < nxgroup && (oldbinx + i) < nx; ++i) {
               for (Int_t j = 0; j < nygroup && (oldbiny + j) < ny; ++j) {
                  Int_t bin = oldbinx + i + (oldbiny + j) * nx;
                  binContent += oldBins[bin];
                  if (oldErrors) binErrorSq += oldErrors[bin];
               }
            }
            Int_t newbin = binx + biny * newnx;
            hnew->UpdateBinContent(newbin, binContent);
            if (oldErrors) hnew->fSumw2[newbin] = binErrorSq;
         }
      }
   }

   // Restore x axis attributes
   fXaxis.SetNdivisions(nXdivisions);
   fXaxis.SetAxisColor(xAxisColor);
   fXaxis.SetLabelColor(xLabelColor);
   fXaxis.SetLabelFont(xLabelFont);
   fXaxis.SetLabelOffset(xLabelOffset);
   fXaxis.SetLabelSize(xLabelSize);
   fXaxis.SetTickLength(xTickLength);
   fXaxis.SetTitleOffset(xTitleOffset);
   fXaxis.SetTitleSize(xTitleSize);
   fXaxis.SetTitleColor(xTitleColor);
   fXaxis.SetTitleFont(xTitleFont);
   // Restore y axis attributes
   fYaxis.SetNdivisions(nYdivisions);
   fYaxis.SetAxisColor(yAxisColor);
   fYaxis.SetLabelColor(yLabelColor);
   fYaxis.SetLabelFont(yLabelFont);
   fYaxis.SetLabelOffset(yLabelOffset);
   fYaxis.SetLabelSize(yLabelSize);
   fYaxis.SetTickLength(yTickLength);
   fYaxis.SetTitleOffset(yTitleOffset);
   fYaxis.SetTitleSize(yTitleSize);
   fYaxis.SetTitleColor(yTitleColor);
   fYaxis.SetTitleFont(yTitleFont);

   if (resetStat) hnew->ResetStats();

   delete [] oldBins;
   if (oldErrors) delete [] oldErrors;
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////

TProfile *TH2::DoProfile(bool onX, const char *name, Int_t firstbin, Int_t lastbin, Option_t *option) const
{
   TString opt = option;
   // extract cut infor
   TString cut;
   Int_t i1 = opt.Index("[");
   if (i1>=0) {
      Int_t i2 = opt.Index("]");
      cut = opt(i1,i2-i1+1);
   }
   opt.ToLower();
   bool originalRange = opt.Contains("o");

   const TAxis& outAxis = ( onX ? fXaxis : fYaxis );
   const TAxis&  inAxis = ( onX ? fYaxis : fXaxis );
   Int_t  inN = inAxis.GetNbins();
   const char *expectedName = ( onX ? "_pfx" : "_pfy" );

   // outer axis cannot be outside original axis (this fixes ROOT-8781)
   // and firstOutBin and lastOutBin cannot be both equal to zero
   Int_t firstOutBin = std::max(outAxis.GetFirst(),1);
   Int_t lastOutBin = std::min(outAxis.GetLast(),outAxis.GetNbins() ) ;

   if ( lastbin < firstbin && inAxis.TestBit(TAxis::kAxisRange) ) {
      firstbin = inAxis.GetFirst();
      lastbin = inAxis.GetLast();
      // For special case of TAxis::SetRange, when first == 1 and last
      // = N and the range bit has been set, the TAxis will return 0
      // for both.
      if (firstbin == 0 && lastbin == 0)
      {
         firstbin = 1;
         lastbin = inAxis.GetNbins();
      }
   }
   if (firstbin < 0) firstbin = 1;
   if (lastbin  < 0) lastbin  = inN;
   if (lastbin  > inN+1) lastbin  = inN;

   // Create the profile histogram
   char *pname = (char*)name;
   if (name && strcmp(name, expectedName) == 0) {
      Int_t nch = strlen(GetName()) + 5;
      pname = new char[nch];
      snprintf(pname,nch,"%s%s",GetName(),name);
   }
   TProfile *h1=0;
   //check if a profile with identical name exist
   // if compatible reset and re-use previous histogram
   TObject *h1obj = gROOT->FindObject(pname);
   if (h1obj && h1obj->InheritsFrom(TH1::Class())) {
      if (h1obj->IsA() != TProfile::Class() ) {
         Error("DoProfile","Histogram with name %s must be a TProfile and is a %s",name,h1obj->ClassName());
         return 0;
      }
      h1 = (TProfile*)h1obj;
      // reset the existing histogram and set always the new binning for the axis
      // This avoid problems when the histogram already exists and the histograms is rebinned or its range has changed
      // (see https://savannah.cern.ch/bugs/?94101 or https://savannah.cern.ch/bugs/?95808 )
      h1->Reset();
      const TArrayD *xbins = outAxis.GetXbins();
      if (xbins->fN == 0) {
         if ( originalRange )
            h1->SetBins(outAxis.GetNbins(),outAxis.GetXmin(),outAxis.GetXmax());
         else
            h1->SetBins(lastOutBin-firstOutBin+1,outAxis.GetBinLowEdge(firstOutBin),outAxis.GetBinUpEdge(lastOutBin));
      } else {
         // case variable bins
         if (originalRange )
            h1->SetBins(outAxis.GetNbins(),xbins->fArray);
         else
            h1->SetBins(lastOutBin-firstOutBin+1,&xbins->fArray[firstOutBin-1]);
      }
   }

   Int_t ncuts = 0;
   if (opt.Contains("[")) {
      ((TH2 *)this)->GetPainter();
      if (fPainter) ncuts = fPainter->MakeCuts((char*)cut.Data());
   }

   if (!h1) {
      const TArrayD *bins = outAxis.GetXbins();
      if (bins->fN == 0) {
         if ( originalRange )
            h1 = new TProfile(pname,GetTitle(),outAxis.GetNbins(),outAxis.GetXmin(),outAxis.GetXmax(),opt);
         else
            h1 = new TProfile(pname,GetTitle(),lastOutBin-firstOutBin+1,
                              outAxis.GetBinLowEdge(firstOutBin),
                              outAxis.GetBinUpEdge(lastOutBin), opt);
      } else {
         // case variable bins
         if (originalRange )
            h1 = new TProfile(pname,GetTitle(),outAxis.GetNbins(),bins->fArray,opt);
         else
            h1 = new TProfile(pname,GetTitle(),lastOutBin-firstOutBin+1,&bins->fArray[firstOutBin-1],opt);
      }
   }
   if (pname != name)  delete [] pname;

   // Copy attributes
   h1->GetXaxis()->ImportAttributes( &outAxis);
   h1->SetLineColor(this->GetLineColor());
   h1->SetFillColor(this->GetFillColor());
   h1->SetMarkerColor(this->GetMarkerColor());
   h1->SetMarkerStyle(this->GetMarkerStyle());

   // check if histogram is weighted
   // in case need to store sum of weight square/bin for the profile
   TArrayD & binSumw2 = *(h1->GetBinSumw2());
   bool useWeights = (GetSumw2N() > 0);
   if (useWeights && (binSumw2.fN != h1->GetNcells()) ) h1->Sumw2();
   // we need to set this bit because we fill the profile using a single Fill for many entries
   // This is needed for the changes applied to make automatically the histogram weighted in ROOT 6 versions
   else h1->SetBit(TH1::kIsNotW);

   // Fill the profile histogram
   // no entries/bin is available so can fill only using bin content as weight
   Double_t totcont = 0;

   // implement filling of projected histogram
   // outbin is bin number of outAxis (the projected axis). Loop is done on all bin of TH2 histograms
   // inbin is the axis being integrated. Loop is done only on the selected bins
   for ( Int_t outbin = 0; outbin <= outAxis.GetNbins() + 1;  ++outbin) {
      if (outAxis.TestBit(TAxis::kAxisRange) && ( outbin < firstOutBin || outbin > lastOutBin )) continue;

      // find corresponding bin number in h1 for outbin (binOut)
      Double_t xOut = outAxis.GetBinCenter(outbin);
      Int_t binOut = h1->GetXaxis()->FindBin( xOut );
      if (binOut <0) continue;

      for (Int_t inbin = firstbin ; inbin <= lastbin ; ++inbin) {
         Int_t binx, biny;
         if (onX) { binx = outbin; biny=inbin; }
         else     { binx = inbin;  biny=outbin; }

         if (ncuts) {
            if (!fPainter->IsInside(binx,biny)) continue;
         }
         Int_t bin = GetBin(binx, biny);
         Double_t cxy = RetrieveBinContent(bin);


         if (cxy) {
            Double_t tmp = 0;
            // the following fill update wrongly the fBinSumw2- need to save it before
            if ( useWeights ) tmp = binSumw2.fArray[binOut];
            h1->Fill( xOut, inAxis.GetBinCenter(inbin), cxy );
            if ( useWeights ) binSumw2.fArray[binOut] = tmp + fSumw2.fArray[bin];
            totcont += cxy;
         }

      }
   }

   // the statistics must be recalculated since by using the Fill method the total sum of weight^2 is
   // not computed correctly
   // for a profile does not much sense to re-use statistics of original TH2
   h1->ResetStats();
   // Also we need to set the entries since they have not been correctly calculated during the projection
   // we can only set them to the effective entries
   h1->SetEntries( h1->GetEffectiveEntries() );


   if (opt.Contains("d")) {
      TVirtualPad *padsav = gPad;
      TVirtualPad *pad = gROOT->GetSelectedPad();
      if (pad) pad->cd();
      opt.Remove(opt.First("d"),1);
      if (!gPad || !gPad->FindObject(h1)) {
         h1->Draw(opt);
      } else {
         h1->Paint(opt);
      }
      if (padsav) padsav->cd();
   }
   return h1;
}


////////////////////////////////////////////////////////////////////////////////
/// Project a 2-D histogram into a profile histogram along X.
///
///   The projection is made from the channels along the Y axis
///   ranging from firstybin to lastybin included.
///   By default, bins 1 to ny are included
///   When all bins are included, the number of entries in the projection
///   is set to the number of entries of the 2-D histogram, otherwise
///   the number of entries is incremented by 1 for all non empty cells.
///
///   if option "d" is specified, the profile is drawn in the current pad.
///
///   if option "o" original axis range of the target axes will be
///   kept, but only bins inside the selected range will be filled.
///
///   The option can also be used to specify the projected profile error type.
///   Values which can be used are 's', 'i', or 'g'. See TProfile::BuildOptions for details
///
///   Using a TCutG object, it is possible to select a sub-range of a 2-D histogram.
///   One must create a graphical cut (mouse or C++) and specify the name
///   of the cut between [] in the option.
///   For example, with a TCutG named "cutg", one can call:
///      myhist->ProfileX(" ",firstybin,lastybin,"[cutg]");
///   To invert the cut, it is enough to put a "-" in front of its name:
///      myhist->ProfileX(" ",firstybin,lastybin,"[-cutg]");
///   It is possible to apply several cuts ("," means logical AND):
///      myhist->ProfileX(" ",firstybin,lastybin,"[cutg1,cutg2]");
///
///   NOTE that if a TProfile named "name" exists in the current directory or pad with
///   a compatible axis the profile is reset and filled again with the projected contents of the TH2.
///   In the case of axis incompatibility an error is reported and a NULL pointer is returned.
///
///   NOTE that the X axis attributes of the TH2 are copied to the X axis of the profile.
///
///   NOTE that the default under- / overflow behavior differs from what ProjectionX
///   does! Profiles take the bin center into account, so here the under- and overflow
///   bins are ignored by default.
///
///   NOTE that the return profile histogram is computed using the Y bin center values instead of
///   the real Y values which are used to fill the 2d histogram. Therefore the obtained profile is just an approximation of the
///   correct profile histogram that would be obtained when filling it directly with the original data (see ROOT-7770)


TProfile *TH2::ProfileX(const char *name, Int_t firstybin, Int_t lastybin, Option_t *option) const
{
   return DoProfile(true, name, firstybin, lastybin, option);

}


////////////////////////////////////////////////////////////////////////////////
/// Project a 2-D histogram into a profile histogram along Y.
///
///   The projection is made from the channels along the X axis
///   ranging from firstxbin to lastxbin included.
///   By default, bins 1 to nx are included
///   When all bins are included, the number of entries in the projection
///   is set to the number of entries of the 2-D histogram, otherwise
///   the number of entries is incremented by 1 for all non empty cells.
///
///   if option "d" is specified, the profile is drawn in the current pad.
///
///   if option "o" , the original axis range of the target axis will be
///   kept, but only bins inside the selected range will be filled.
///
///   The option can also be used to specify the projected profile error type.
///   Values which can be used are 's', 'i', or 'g'. See TProfile::BuildOptions for details
///   Using a TCutG object, it is possible to select a sub-range of a 2-D histogram.
///
///   One must create a graphical cut (mouse or C++) and specify the name
///   of the cut between [] in the option.
///   For example, with a TCutG named "cutg", one can call:
///      myhist->ProfileY(" ",firstybin,lastybin,"[cutg]");
///   To invert the cut, it is enough to put a "-" in front of its name:
///      myhist->ProfileY(" ",firstybin,lastybin,"[-cutg]");
///   It is possible to apply several cuts:
///      myhist->ProfileY(" ",firstybin,lastybin,"[cutg1,cutg2]");
///
///   NOTE that if a TProfile named "name" exists in the current directory or pad with
///   a compatible axis the profile is reset and filled again with the projected contents of the TH2.
///   In the case of axis incompatibility an error is reported and a NULL pointer is returned.
///
///   NOTE that the Y axis attributes of the TH2 are copied to the X axis of the profile.
///
///   NOTE that the default under- / overflow behavior differs from what ProjectionX
///   does! Profiles take the bin center into account, so here the under- and overflow
///   bins are ignored by default.
///
///   NOTE that the return profile histogram is computed using the X bin center values instead of
///   the real X values which are used to fill the 2d histogram. Therefore the obtained profile is just an approximation of the
///   correct profile histogram that would be obtained when filling it directly with the original data (see ROOT-7770)


TProfile *TH2::ProfileY(const char *name, Int_t firstxbin, Int_t lastxbin, Option_t *option) const
{
   return DoProfile(false, name, firstxbin, lastxbin, option);
}


////////////////////////////////////////////////////////////////////////////////
/// Internal (protected) method for performing projection on the X or Y axis
/// called by ProjectionX or ProjectionY

TH1D *TH2::DoProjection(bool onX, const char *name, Int_t firstbin, Int_t lastbin, Option_t *option) const
{
   const char *expectedName = 0;
   Int_t inNbin;
   const TAxis* outAxis;
   const TAxis* inAxis;

   TString opt = option;
   TString cut;
   Int_t i1 = opt.Index("[");
   if (i1>=0) {
      Int_t i2 = opt.Index("]");
      cut = opt(i1,i2-i1+1);
   }
   opt.ToLower();  //must be called after having parsed the cut name
   bool originalRange = opt.Contains("o");

   if ( onX )
   {
      expectedName = "_px";
      inNbin = fYaxis.GetNbins();
      outAxis = GetXaxis();
      inAxis = GetYaxis();
   }
   else
   {
      expectedName = "_py";
      inNbin = fXaxis.GetNbins();
      outAxis = GetYaxis();
      inAxis = GetXaxis();
   }

   // outer axis cannot be outside original axis (this fixes ROOT-8781)
   // and firstOutBin and lastOutBin cannot be both equal to zero
   Int_t firstOutBin = std::max(outAxis->GetFirst(),1);
   Int_t lastOutBin = std::min(outAxis->GetLast(),outAxis->GetNbins() ) ;

   if ( lastbin < firstbin && inAxis->TestBit(TAxis::kAxisRange) ) {
      firstbin = inAxis->GetFirst();
      lastbin = inAxis->GetLast();
      // For special case of TAxis::SetRange, when first == 1 and last
      // = N and the range bit has been set, the TAxis will return 0
      // for both.
      if (firstbin == 0 && lastbin == 0)
      {
         firstbin = 1;
         lastbin = inAxis->GetNbins();
      }
   }
   if (firstbin < 0) firstbin = 0;
   if (lastbin  < 0) lastbin  = inNbin + 1;
   if (lastbin  > inNbin+1) lastbin  = inNbin + 1;

   // Create the projection histogram
   char *pname = (char*)name;
   if (name && strcmp(name,expectedName) == 0) {
      Int_t nch = strlen(GetName()) + 4;
      pname = new char[nch];
      snprintf(pname,nch,"%s%s",GetName(),name);
   }
   TH1D *h1=0;
   //check if histogram with identical name exist
   // if compatible reset and re-use previous histogram
   // (see https://savannah.cern.ch/bugs/?54340)
   TObject *h1obj = gROOT->FindObject(pname);
   if (h1obj && h1obj->InheritsFrom(TH1::Class())) {
      if (h1obj->IsA() != TH1D::Class() ) {
         Error("DoProjection","Histogram with name %s must be a TH1D and is a %s",name,h1obj->ClassName());
         return 0;
      }
      h1 = (TH1D*)h1obj;
      // reset the existing histogram and set always the new binning for the axis
      // This avoid problems when the histogram already exists and the histograms is rebinned or its range has changed
      // (see https://savannah.cern.ch/bugs/?94101 or https://savannah.cern.ch/bugs/?95808 )
      h1->Reset();
      const TArrayD *xbins = outAxis->GetXbins();
      if (xbins->fN == 0) {
         if ( originalRange )
            h1->SetBins(outAxis->GetNbins(),outAxis->GetXmin(),outAxis->GetXmax());
         else
            h1->SetBins(lastOutBin-firstOutBin+1,outAxis->GetBinLowEdge(firstOutBin),outAxis->GetBinUpEdge(lastOutBin));
      } else {
         // case variable bins
         if (originalRange )
            h1->SetBins(outAxis->GetNbins(),xbins->fArray);
         else
            h1->SetBins(lastOutBin-firstOutBin+1,&xbins->fArray[firstOutBin-1]);
      }
   }

   Int_t ncuts = 0;
   if (opt.Contains("[")) {
      ((TH2 *)this)->GetPainter();
      if (fPainter) ncuts = fPainter->MakeCuts((char*)cut.Data());
   }

   if (!h1) {
      const TArrayD *bins = outAxis->GetXbins();
      if (bins->fN == 0) {
         if ( originalRange )
            h1 = new TH1D(pname,GetTitle(),outAxis->GetNbins(),outAxis->GetXmin(),outAxis->GetXmax());
         else
            h1 = new TH1D(pname,GetTitle(),lastOutBin-firstOutBin+1,
                          outAxis->GetBinLowEdge(firstOutBin),outAxis->GetBinUpEdge(lastOutBin));
      } else {
         // case variable bins
         if (originalRange )
            h1 = new TH1D(pname,GetTitle(),outAxis->GetNbins(),bins->fArray);
         else
            h1 = new TH1D(pname,GetTitle(),lastOutBin-firstOutBin+1,&bins->fArray[firstOutBin-1]);
      }
      if (opt.Contains("e") || GetSumw2N() ) h1->Sumw2();
   }
   if (pname != name)  delete [] pname;

   // Copy the axis attributes and the axis labels if needed.
   h1->GetXaxis()->ImportAttributes(outAxis);
   THashList* labels=outAxis->GetLabels();
   if (labels) {
      TIter iL(labels);
      TObjString* lb;
      Int_t i = 1;
      while ((lb=(TObjString*)iL())) {
         h1->GetXaxis()->SetBinLabel(i,lb->String().Data());
         i++;
      }
   }

   h1->SetLineColor(this->GetLineColor());
   h1->SetFillColor(this->GetFillColor());
   h1->SetMarkerColor(this->GetMarkerColor());
   h1->SetMarkerStyle(this->GetMarkerStyle());

   // Fill the projected histogram
   Double_t cont,err2;
   Double_t totcont = 0;
   Bool_t  computeErrors = h1->GetSumw2N();

   // implement filling of projected histogram
   // outbin is bin number of outAxis (the projected axis). Loop is done on all bin of TH2 histograms
   // inbin is the axis being integrated. Loop is done only on the selected bins
   for ( Int_t outbin = 0; outbin <= outAxis->GetNbins() + 1;  ++outbin) {
      err2 = 0;
      cont = 0;
      if (outAxis->TestBit(TAxis::kAxisRange) && ( outbin < firstOutBin || outbin > lastOutBin )) continue;

      for (Int_t inbin = firstbin ; inbin <= lastbin ; ++inbin) {
         Int_t binx, biny;
         if (onX) { binx = outbin; biny=inbin; }
         else     { binx = inbin;  biny=outbin; }

         if (ncuts) {
            if (!fPainter->IsInside(binx,biny)) continue;
         }
         // sum bin content and error if needed
         cont  += GetBinContent(binx,biny);
         if (computeErrors) {
            Double_t exy = GetBinError(binx,biny);
            err2  += exy*exy;
         }
      }
      // find corresponding bin number in h1 for outbin
      Int_t binOut = h1->GetXaxis()->FindBin( outAxis->GetBinCenter(outbin) );
      h1->SetBinContent(binOut ,cont);
      if (computeErrors) h1->SetBinError(binOut,TMath::Sqrt(err2));
      // sum  all content
      totcont += cont;
   }

   // check if we can re-use the original statistics from  the previous histogram
   bool reuseStats = false;
   if ( ( GetStatOverflowsBehaviour() == false && firstbin == 1 && lastbin == inNbin     ) ||
        ( GetStatOverflowsBehaviour() == true  && firstbin == 0 && lastbin == inNbin + 1 ) )
      reuseStats = true;
   else {
      // also if total content match we can re-use
      double eps = 1.E-12;
      if (IsA() == TH2F::Class() ) eps = 1.E-6;
      if (fTsumw != 0 && TMath::Abs( fTsumw - totcont) <  TMath::Abs(fTsumw) * eps)
         reuseStats = true;
   }
   if (ncuts) reuseStats = false;
   // retrieve  the statistics and set in projected histogram if we can re-use it
   bool reuseEntries = reuseStats;
   // can re-use entries if underflow/overflow are included
   reuseEntries &= (firstbin==0 && lastbin == inNbin+1);
   if (reuseStats) {
      Double_t stats[kNstat];
      GetStats(stats);
      if (!onX) {  // case of projection on Y
         stats[2] = stats[4];
         stats[3] = stats[5];
      }
      h1->PutStats(stats);
   }
   else {
      // the statistics is automatically recalculated since it is reset by the call to SetBinContent
      // we just need to set the entries since they have not been correctly calculated during the projection
      // we can only set them to the effective entries
      h1->SetEntries( h1->GetEffectiveEntries() );
   }
   if (reuseEntries) {
      h1->SetEntries(fEntries);
   }
   else {
      // re-compute the entries
      // in case of error calculation (i.e. when Sumw2() is set)
      // use the effective entries for the entries
      // since this  is the only way to estimate them
      Double_t entries =  TMath::Floor( totcont + 0.5); // to avoid numerical rounding
      if (h1->GetSumw2N()) entries = h1->GetEffectiveEntries();
      h1->SetEntries( entries );
   }

   if (opt.Contains("d")) {
      TVirtualPad *padsav = gPad;
      TVirtualPad *pad = gROOT->GetSelectedPad();
      if (pad) pad->cd();
      opt.Remove(opt.First("d"),1);
      // remove also other options
      if (opt.Contains("e")) opt.Remove(opt.First("e"),1);
      if (!gPad || !gPad->FindObject(h1)) {
         h1->Draw(opt);
      } else {
         h1->Paint(opt);
      }
      if (padsav) padsav->cd();
   }

   return h1;
}


////////////////////////////////////////////////////////////////////////////////
/// Project a 2-D histogram into a 1-D histogram along X.
///
///   The projection is always of the type TH1D.
///   The projection is made from the channels along the Y axis
///   ranging from firstybin to lastybin included.
///   By default, all bins including under- and overflow are included.
///   The number of entries in the projection is estimated from the
///   number of effective entries for all the cells included in the projection.
///
///   To exclude the underflow bins in Y, use firstybin=1.
///   To exclude the overflow bins in Y, use lastybin=nx.
///
///   if option "e" is specified, the errors are computed.
///   if option "d" is specified, the projection is drawn in the current pad.
///   if option "o" original axis range of the taget axes will be
///   kept, but only bins inside the selected range will be filled.
///
///   Using a TCutG object, it is possible to select a sub-range of a 2-D histogram.
///   One must create a graphical cut (mouse or C++) and specify the name
///   of the cut between [] in the option.
///   For example, with a TCutG named "cutg", one can call:
///      myhist->ProjectionX(" ",firstybin,lastybin,"[cutg]");
///   To invert the cut, it is enough to put a "-" in front of its name:
///      myhist->ProjectionX(" ",firstybin,lastybin,"[-cutg]");
///   It is possible to apply several cuts:
///      myhist->ProjectionX(" ",firstybin,lastybin,"[cutg1,cutg2]");
///
///   NOTE that if a TH1D named "name" exists in the current directory or pad
///   the histogram is reset and filled again with the projected contents of the TH2.
///
///   NOTE that the X axis attributes of the TH2 are copied to the X axis of the projection.

TH1D *TH2::ProjectionX(const char *name, Int_t firstybin, Int_t lastybin, Option_t *option) const
{
   return DoProjection(true, name, firstybin, lastybin, option);
}


////////////////////////////////////////////////////////////////////////////////
/// Project a 2-D histogram into a 1-D histogram along Y.
///
///   The projection is always of the type TH1D.
///   The projection is made from the channels along the X axis
///   ranging from firstxbin to lastxbin included.
///   By default, all bins including under- and overflow are included.
///   The number of entries in the projection is estimated from the
///   number of effective entries for all the cells included in the projection
///
///   To exclude the underflow bins in X, use firstxbin=1.
///   To exclude the overflow bins in X, use lastxbin=nx.
///
///   if option "e" is specified, the errors are computed.
///   if option "d" is specified, the projection is drawn in the current pad.
///   if option "o" original axis range of the taget axes will be
///   kept, but only bins inside the selected range will be filled.
///
///   Using a TCutG object, it is possible to select a sub-range of a 2-D histogram.
///   One must create a graphical cut (mouse or C++) and specify the name
///   of the cut between [] in the option.
///   For example, with a TCutG named "cutg", one can call:
///      myhist->ProjectionY(" ",firstxbin,lastxbin,"[cutg]");
///   To invert the cut, it is enough to put a "-" in front of its name:
///      myhist->ProjectionY(" ",firstxbin,lastxbin,"[-cutg]");
///   It is possible to apply several cuts:
///      myhist->ProjectionY(" ",firstxbin,lastxbin,"[cutg1,cutg2]");
///
///   NOTE that if a TH1D named "name" exists in the current directory or pad and having
///   a compatible axis, the histogram is reset and filled again with the projected contents of the TH2.
///   In the case of axis incompatibility, an error is reported and a NULL pointer is returned.
///
///   NOTE that the Y axis attributes of the TH2 are copied to the X axis of the projection.

TH1D *TH2::ProjectionY(const char *name, Int_t firstxbin, Int_t lastxbin, Option_t *option) const
{
   return DoProjection(false, name, firstxbin, lastxbin, option);
}


////////////////////////////////////////////////////////////////////////////////
/// Replace current statistics with the values in array stats

void TH2::PutStats(Double_t *stats)
{
   TH1::PutStats(stats);
   fTsumwy  = stats[4];
   fTsumwy2 = stats[5];
   fTsumwxy = stats[6];
}


////////////////////////////////////////////////////////////////////////////////
/// Compute the X distribution of quantiles in the other variable Y
/// name is the name of the returned histogram
/// prob is the probability content for the quantile (0.5 is the default for the median)
/// An approximate error for the quantile is computed assuming that the distribution in
/// the other variable is normal. According to this approximate formula the error on the quantile is
/// estimated as sqrt( p (1-p) / ( n * f(q)^2) ), where p is the probability content of the quantile and
/// n is the number of events used to compute the quantile and f(q) is the probability distribution for the
/// other variable evaluated at the obtained quantile. In the error estimation the probability is then assumed to be
///  a normal distribution.

TH1D* TH2::QuantilesX( Double_t prob, const char * name) const
{
   return DoQuantiles(true, name, prob);
}


////////////////////////////////////////////////////////////////////////////////
/// Compute the Y distribution of quantiles in the other variable X
/// name is the name of the returned histogram
/// prob is the probability content for the quantile (0.5 is the default for the median)
/// An approximate error for the quantile is computed assuming that the distribution in
/// the other variable is normal.

TH1D* TH2::QuantilesY( Double_t prob, const char * name) const
{
   return DoQuantiles(false, name, prob);
}


////////////////////////////////////////////////////////////////////////////////
/// Implementation of quantiles for x or y

TH1D* TH2::DoQuantiles(bool onX, const char * name, Double_t prob) const
{
   const TAxis *outAxis = 0;
   if ( onX )   {
      outAxis = GetXaxis();
   }  else {
      outAxis = GetYaxis();
   }

   // build first name of returned histogram
   TString qname = name;
   if (qname.IsNull() || qname == "_qx" || qname == "_qy") {
      const char * qtype = (onX) ? "qx" : "qy";
      qname = TString::Format("%s_%s_%3.2f",GetName(),qtype, prob);
   }
   // check if the histogram is already existing
   TH1D *h1=0;
   //check if histogram with identical name exist
   TObject *h1obj = gROOT->FindObject(qname);
   if (h1obj) {
      h1 = dynamic_cast<TH1D*>(h1obj);
      if (!h1) {
         Error("DoQuantiles","Histogram with name %s must be a TH1D and is a %s",qname.Data(),h1obj->ClassName());
         return 0;
      }
   }
   if (h1) {
      h1->Reset();
   } else {
      // create the histogram
      h1 = new TH1D(qname, GetTitle(), 1, 0, 1);
   }
   // set the bin content
   Int_t firstOutBin = std::max(outAxis->GetFirst(),1);
   Int_t lastOutBin = std::max(outAxis->GetLast(),outAxis->GetNbins());
   const TArrayD *xbins = outAxis->GetXbins();
   if (xbins->fN == 0)
      h1->SetBins(lastOutBin-firstOutBin+1,outAxis->GetBinLowEdge(firstOutBin),outAxis->GetBinUpEdge(lastOutBin));
   else
      h1->SetBins(lastOutBin-firstOutBin+1,&xbins->fArray[firstOutBin-1]);

   // set the bin content of the histogram
  Double_t pp[1];
  pp[0] = prob;

  TH1D * slice = 0;
  for (int ibin = outAxis->GetFirst() ; ibin <= outAxis->GetLast() ; ++ibin) {
    Double_t qq[1];
    // do a projection on the opposite axis
    slice = DoProjection(!onX, "tmp",ibin,ibin,"");
    if (!slice) break;
    if (slice->GetSum() == 0) continue;
    slice->GetQuantiles(1,qq,pp);
    h1->SetBinContent(ibin,qq[0]);
    // compute error using normal approximation
    // quantile error  ~  sqrt (q*(1-q)/ *( n * f(xq)^2 ) from Kendall
    // where f(xq) is the p.d.f value at the quantile xq
    Double_t n = slice->GetEffectiveEntries();
    Double_t f = TMath::Gaus(qq[0], slice->GetMean(), slice->GetStdDev(), kTRUE);
    Double_t error = 0;
    // set the errors to zero in case of small statistics
    if (f > 0 && n > 1)
       error = TMath::Sqrt( prob*(1.-prob)/ (n * f * f) );
    h1->SetBinError(ibin, error);
  }
  if (slice) delete slice;
  return h1;
}


////////////////////////////////////////////////////////////////////////////////
/// Reset this histogram: contents, errors, etc.

void TH2::Reset(Option_t *option)
{
   TH1::Reset(option);
   TString opt = option;
   opt.ToUpper();

   if (opt.Contains("ICE") && !opt.Contains("S")) return;
   fTsumwy  = 0;
   fTsumwy2 = 0;
   fTsumwxy = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Set bin content

void TH2::SetBinContent(Int_t bin, Double_t content)
{
   fEntries++;
   fTsumw = 0;
   if (bin < 0) return;
   if (bin >= fNcells) return;
   UpdateBinContent(bin, content);
}


////////////////////////////////////////////////////////////////////////////////
/// When the mouse is moved in a pad containing a 2-d view of this histogram
/// a second canvas shows the projection along X corresponding to the
/// mouse position along Y.
/// To stop the generation of the projections, delete the canvas
/// containing the projection.

void TH2::SetShowProjectionX(Int_t nbins)
{
   GetPainter();

   if (fPainter) fPainter->SetShowProjection("x",nbins);
}


////////////////////////////////////////////////////////////////////////////////
/// When the mouse is moved in a pad containing a 2-d view of this histogram
/// a second canvas shows the projection along Y corresponding to the
/// mouse position along X.
/// To stop the generation of the projections, delete the canvas
/// containing the projection.

void TH2::SetShowProjectionY(Int_t nbins)
{
   GetPainter();

   if (fPainter) fPainter->SetShowProjection("y",nbins);
}


////////////////////////////////////////////////////////////////////////////////
///   This function calculates the background spectrum in this histogram.
///   The background is returned as a histogram.
///   to be implemented (may be)

TH1 *TH2::ShowBackground(Int_t niter, Option_t *option)
{

   return (TH1*)gROOT->ProcessLineFast(Form("TSpectrum2::StaticBackground((TH1*)0x%lx,%d,\"%s\")",
                                            (ULong_t)this, niter, option));
}


////////////////////////////////////////////////////////////////////////////////
///Interface to TSpectrum2::Search
///the function finds peaks in this histogram where the width is > sigma
///and the peak maximum greater than threshold*maximum bin content of this.
///for more details see TSpectrum::Search.
///note the difference in the default value for option compared to TSpectrum2::Search
///option="" by default (instead of "goff")

Int_t TH2::ShowPeaks(Double_t sigma, Option_t *option, Double_t threshold)
{

   return (Int_t)gROOT->ProcessLineFast(Form("TSpectrum2::StaticSearch((TH1*)0x%lx,%g,\"%s\",%g)",
                                             (ULong_t)this, sigma, option, threshold));
}


////////////////////////////////////////////////////////////////////////////////
/// Smooth bin contents of this 2-d histogram using kernel algorithms
/// similar to the ones used in the raster graphics community.
/// Bin contents in the active range are replaced by their smooth values.
/// The algorithm retains the input dimension by using Kernel Crop at the input boundaries.
/// Kernel Crop sets any pixel in the kernel that extends past the input to zero and adjusts the
/// normalization accordingly.
/// If Errors are defined via Sumw2, they are also scaled and computed.
/// However, note the resulting errors will be correlated between different-bins, so
/// the errors should not be used blindly to perform any calculation involving several bins,
/// like fitting the histogram.  One would need to compute also the bin by bin correlation matrix.
///
/// 3 kernels are proposed k5a, k5b and k3a.
/// k5a and k5b act on 5x5 cells (i-2,i-1,i,i+1,i+2, and same for j)
/// k5b is a bit more stronger in smoothing
/// k3a acts only on 3x3 cells (i-1,i,i+1, and same for j).
/// By default the kernel "k5a" is used. You can select the kernels "k5b" or "k3a"
/// via the option argument.
/// If TAxis::SetRange has been called on the x or/and y axis, only the bins
/// in the specified range are smoothed.
/// In the current implementation if the first argument is not used (default value=1).
///
/// implementation by David McKee (dmckee@bama.ua.edu). Extended by Rene Brun

void TH2::Smooth(Int_t ntimes, Option_t *option)
{
   Double_t k5a[5][5] =  { { 0, 0, 1, 0, 0 },
                           { 0, 2, 2, 2, 0 },
                           { 1, 2, 5, 2, 1 },
                           { 0, 2, 2, 2, 0 },
                           { 0, 0, 1, 0, 0 } };
   Double_t k5b[5][5] =  { { 0, 1, 2, 1, 0 },
                           { 1, 2, 4, 2, 1 },
                           { 2, 4, 8, 4, 2 },
                           { 1, 2, 4, 2, 1 },
                           { 0, 1, 2, 1, 0 } };
   Double_t k3a[3][3] =  { { 0, 1, 0 },
                           { 1, 2, 1 },
                           { 0, 1, 0 } };

   if (ntimes > 1) {
      Warning("Smooth","Currently only ntimes=1 is supported");
   }
   TString opt = option;
   opt.ToLower();
   Int_t ksize_x=5;
   Int_t ksize_y=5;
   Double_t *kernel = &k5a[0][0];
   if (opt.Contains("k5b")) kernel = &k5b[0][0];
   if (opt.Contains("k3a")) {
      kernel = &k3a[0][0];
      ksize_x=3;
      ksize_y=3;
   }

   // find i,j ranges
   Int_t ifirst = fXaxis.GetFirst();
   Int_t ilast  = fXaxis.GetLast();
   Int_t jfirst = fYaxis.GetFirst();
   Int_t jlast  = fYaxis.GetLast();

   // Determine the size of the bin buffer(s) needed
   Double_t nentries = fEntries;
   Int_t nx = GetNbinsX();
   Int_t ny = GetNbinsY();
   Int_t bufSize  = (nx+2)*(ny+2);
   Double_t *buf  = new Double_t[bufSize];
   Double_t *ebuf = 0;
   if (fSumw2.fN) ebuf = new Double_t[bufSize];

   // Copy all the data to the temporary buffers
   Int_t i,j,bin;
   for (i=ifirst; i<=ilast; i++){
      for (j=jfirst; j<=jlast; j++){
         bin = GetBin(i,j);
         buf[bin] = RetrieveBinContent(bin);
         if (ebuf) ebuf[bin]=GetBinError(bin);
      }
   }

   // Kernel tail sizes (kernel sizes must be odd for this to work!)
   Int_t x_push = (ksize_x-1)/2;
   Int_t y_push = (ksize_y-1)/2;

   // main work loop
   for (i=ifirst; i<=ilast; i++){
      for (j=jfirst; j<=jlast; j++) {
         Double_t content = 0.0;
         Double_t error = 0.0;
         Double_t norm = 0.0;

         for (Int_t n=0; n<ksize_x; n++) {
            for (Int_t m=0; m<ksize_y; m++) {
               Int_t xb = i+(n-x_push);
               Int_t yb = j+(m-y_push);
               if ( (xb >= 1) && (xb <= nx) && (yb >= 1) && (yb <= ny) ) {
                  bin = GetBin(xb,yb);
                  Double_t k = kernel[n*ksize_y +m];
                  //if ( (k != 0.0 ) && (buf[bin] != 0.0) ) { // General version probably does not want the second condition
                  if ( k != 0.0 ) {
                     norm    += k;
                     content += k*buf[bin];
                     if (ebuf) error   += k*k*ebuf[bin]*ebuf[bin];
                  }
               }
            }
         }

         if ( norm != 0.0 ) {
            SetBinContent(i,j,content/norm);
            if (ebuf) {
               error /= (norm*norm);
               SetBinError(i,j,sqrt(error));
            }
         }
      }
   }
   fEntries = nentries;

   delete [] buf;
   delete [] ebuf;
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TH2.

void TH2::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TH2::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TH1::Streamer(R__b);
      R__b >> fScalefactor;
      R__b >> fTsumwy;
      R__b >> fTsumwy2;
      R__b >> fTsumwxy;
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TH2::Class(),this);
   }
}


//______________________________________________________________________________
//                     TH2C methods
//  TH2C a 2-D histogram with one byte per cell (char)
//______________________________________________________________________________

ClassImp(TH2C);


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2C::TH2C(): TH2(), TArrayC()
{
   SetBinsLength(9);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH2C::~TH2C()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2C::TH2C(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,Double_t ylow,Double_t yup)
           :TH2(name,title,nbinsx,xlow,xup,nbinsy,ylow,yup)
{
   TArrayC::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();

   if (xlow >= xup || ylow >= yup) SetBuffer(fgBufferSize);
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2C::TH2C(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,Double_t ylow,Double_t yup)
           :TH2(name,title,nbinsx,xbins,nbinsy,ylow,yup)
{
   TArrayC::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2C::TH2C(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,const Double_t *ybins)
           :TH2(name,title,nbinsx,xlow,xup,nbinsy,ybins)
{
   TArrayC::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2C::TH2C(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,const Double_t *ybins)
           :TH2(name,title,nbinsx,xbins,nbinsy,ybins)
{
   TArrayC::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2C::TH2C(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins
           ,Int_t nbinsy,const Float_t *ybins)
           :TH2(name,title,nbinsx,xbins,nbinsy,ybins)
{
   TArrayC::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TH2C::TH2C(const TH2C &h2c) : TH2(), TArrayC()
{
   ((TH2C&)h2c).Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by 1.

void TH2C::AddBinContent(Int_t bin)
{
   if (fArray[bin] < 127) fArray[bin]++;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by w.

void TH2C::AddBinContent(Int_t bin, Double_t w)
{
   Int_t newval = fArray[bin] + Int_t(w);
   if (newval > -128 && newval < 128) {fArray[bin] = Char_t(newval); return;}
   if (newval < -127) fArray[bin] = -127;
   if (newval >  127) fArray[bin] =  127;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy.

void TH2C::Copy(TObject &newth2) const
{
   TH2::Copy((TH2C&)newth2);
}


////////////////////////////////////////////////////////////////////////////////
/// Reset this histogram: contents, errors, etc.

void TH2C::Reset(Option_t *option)
{
   TH2::Reset(option);
   TArrayC::Reset();
}


////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH2C::SetBinsLength(Int_t n)
{
   if (n < 0) n = (fXaxis.GetNbins()+2)*(fYaxis.GetNbins()+2);
   fNcells = n;
   TArrayC::Set(n);
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TH2C.

void TH2C::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TH2C::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      if (R__v < 2) {
         R__b.ReadVersion();
         TH1::Streamer(R__b);
         TArrayC::Streamer(R__b);
         R__b.ReadVersion();
         R__b >> fScalefactor;
         R__b >> fTsumwy;
         R__b >> fTsumwy2;
         R__b >> fTsumwxy;
      } else {
         TH2::Streamer(R__b);
         TArrayC::Streamer(R__b);
         R__b.CheckByteCount(R__s, R__c, TH2C::IsA());
      }
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TH2C::Class(),this);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH2C& TH2C::operator=(const TH2C &h1)
{
   if (this != &h1)  ((TH2C&)h1).Copy(*this);
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH2C operator*(Float_t c1, TH2C &h1)
{
   TH2C hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH2C operator+(TH2C &h1, TH2C &h2)
{
   TH2C hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator -

TH2C operator-(TH2C &h1, TH2C &h2)
{
   TH2C hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH2C operator*(TH2C &h1, TH2C &h2)
{
   TH2C hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH2C operator/(TH2C &h1, TH2C &h2)
{
   TH2C hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


//______________________________________________________________________________
//                     TH2S methods
//  TH2S a 2-D histogram with two bytes per cell (short integer)
//______________________________________________________________________________

ClassImp(TH2S);


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2S::TH2S(): TH2(), TArrayS()
{
   SetBinsLength(9);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH2S::~TH2S()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2S::TH2S(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,Double_t ylow,Double_t yup)
           :TH2(name,title,nbinsx,xlow,xup,nbinsy,ylow,yup)
{
   TArrayS::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();

   if (xlow >= xup || ylow >= yup) SetBuffer(fgBufferSize);
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2S::TH2S(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,Double_t ylow,Double_t yup)
           :TH2(name,title,nbinsx,xbins,nbinsy,ylow,yup)
{
   TArrayS::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2S::TH2S(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,const Double_t *ybins)
           :TH2(name,title,nbinsx,xlow,xup,nbinsy,ybins)
{
   TArrayS::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2S::TH2S(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,const Double_t *ybins)
           :TH2(name,title,nbinsx,xbins,nbinsy,ybins)
{
   TArrayS::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2S::TH2S(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins
           ,Int_t nbinsy,const Float_t *ybins)
           :TH2(name,title,nbinsx,xbins,nbinsy,ybins)
{
   TArrayS::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TH2S::TH2S(const TH2S &h2s) : TH2(), TArrayS()
{
   ((TH2S&)h2s).Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by 1.

void TH2S::AddBinContent(Int_t bin)
{
   if (fArray[bin] < 32767) fArray[bin]++;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by w.

void TH2S::AddBinContent(Int_t bin, Double_t w)
{
   Int_t newval = fArray[bin] + Int_t(w);
   if (newval > -32768 && newval < 32768) {fArray[bin] = Short_t(newval); return;}
   if (newval < -32767) fArray[bin] = -32767;
   if (newval >  32767) fArray[bin] =  32767;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy.

void TH2S::Copy(TObject &newth2) const
{
   TH2::Copy((TH2S&)newth2);
}


////////////////////////////////////////////////////////////////////////////////
/// Reset this histogram: contents, errors, etc.

void TH2S::Reset(Option_t *option)
{
   TH2::Reset(option);
   TArrayS::Reset();
}


////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH2S::SetBinsLength(Int_t n)
{
   if (n < 0) n = (fXaxis.GetNbins()+2)*(fYaxis.GetNbins()+2);
   fNcells = n;
   TArrayS::Set(n);
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TH2S.

void TH2S::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TH2S::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      if (R__v < 2) {
         R__b.ReadVersion();
         TH1::Streamer(R__b);
         TArrayS::Streamer(R__b);
         R__b.ReadVersion();
         R__b >> fScalefactor;
         R__b >> fTsumwy;
         R__b >> fTsumwy2;
         R__b >> fTsumwxy;
      } else {
         TH2::Streamer(R__b);
         TArrayS::Streamer(R__b);
         R__b.CheckByteCount(R__s, R__c, TH2S::IsA());
      }
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TH2S::Class(),this);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH2S& TH2S::operator=(const TH2S &h1)
{
   if (this != &h1)  ((TH2S&)h1).Copy(*this);
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH2S operator*(Float_t c1, TH2S &h1)
{
   TH2S hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH2S operator+(TH2S &h1, TH2S &h2)
{
   TH2S hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator -

TH2S operator-(TH2S &h1, TH2S &h2)
{
   TH2S hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH2S operator*(TH2S &h1, TH2S &h2)
{
   TH2S hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH2S operator/(TH2S &h1, TH2S &h2)
{
   TH2S hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


//______________________________________________________________________________
//                     TH2I methods
//  TH2I a 2-D histogram with four bytes per cell (32 bits integer)
//______________________________________________________________________________

ClassImp(TH2I);


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2I::TH2I(): TH2(), TArrayI()
{
   SetBinsLength(9);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH2I::~TH2I()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2I::TH2I(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,Double_t ylow,Double_t yup)
           :TH2(name,title,nbinsx,xlow,xup,nbinsy,ylow,yup)
{
   TArrayI::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();

   if (xlow >= xup || ylow >= yup) SetBuffer(fgBufferSize);
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2I::TH2I(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,Double_t ylow,Double_t yup)
           :TH2(name,title,nbinsx,xbins,nbinsy,ylow,yup)
{
   TArrayI::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2I::TH2I(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,const Double_t *ybins)
           :TH2(name,title,nbinsx,xlow,xup,nbinsy,ybins)
{
   TArrayI::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2I::TH2I(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,const Double_t *ybins)
           :TH2(name,title,nbinsx,xbins,nbinsy,ybins)
{
   TArrayI::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2I::TH2I(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins
           ,Int_t nbinsy,const Float_t *ybins)
           :TH2(name,title,nbinsx,xbins,nbinsy,ybins)
{
   TArrayI::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TH2I::TH2I(const TH2I &h2i) : TH2(), TArrayI()
{
   ((TH2I&)h2i).Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by 1.

void TH2I::AddBinContent(Int_t bin)
{
   if (fArray[bin] < 2147483647) fArray[bin]++;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by w.

void TH2I::AddBinContent(Int_t bin, Double_t w)
{
   Long64_t newval = fArray[bin] + Long64_t(w);
   if (newval > -2147483647 && newval < 2147483647) {fArray[bin] = Int_t(newval); return;}
   if (newval < -2147483647) fArray[bin] = -2147483647;
   if (newval >  2147483647) fArray[bin] =  2147483647;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy.

void TH2I::Copy(TObject &newth2) const
{
   TH2::Copy((TH2I&)newth2);
}


////////////////////////////////////////////////////////////////////////////////
/// Reset this histogram: contents, errors, etc.

void TH2I::Reset(Option_t *option)
{
   TH2::Reset(option);
   TArrayI::Reset();
}


////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH2I::SetBinsLength(Int_t n)
{
   if (n < 0) n = (fXaxis.GetNbins()+2)*(fYaxis.GetNbins()+2);
   fNcells = n;
   TArrayI::Set(n);
}


////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH2I& TH2I::operator=(const TH2I &h1)
{
   if (this != &h1)  ((TH2I&)h1).Copy(*this);
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH2I operator*(Float_t c1, TH2I &h1)
{
   TH2I hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH2I operator+(TH2I &h1, TH2I &h2)
{
   TH2I hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator -

TH2I operator-(TH2I &h1, TH2I &h2)
{
   TH2I hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH2I operator*(TH2I &h1, TH2I &h2)
{
   TH2I hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH2I operator/(TH2I &h1, TH2I &h2)
{
   TH2I hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


//______________________________________________________________________________
//                     TH2F methods
//  TH2F a 2-D histogram with four bytes per cell (float)
//______________________________________________________________________________

ClassImp(TH2F);


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2F::TH2F(): TH2(), TArrayF()
{
   SetBinsLength(9);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH2F::~TH2F()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2F::TH2F(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,Double_t ylow,Double_t yup)
           :TH2(name,title,nbinsx,xlow,xup,nbinsy,ylow,yup)
{
   TArrayF::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();

   if (xlow >= xup || ylow >= yup) SetBuffer(fgBufferSize);
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2F::TH2F(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,Double_t ylow,Double_t yup)
           :TH2(name,title,nbinsx,xbins,nbinsy,ylow,yup)
{
   TArrayF::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2F::TH2F(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,const Double_t *ybins)
           :TH2(name,title,nbinsx,xlow,xup,nbinsy,ybins)
{
   TArrayF::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2F::TH2F(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,const Double_t *ybins)
           :TH2(name,title,nbinsx,xbins,nbinsy,ybins)
{
   TArrayF::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2F::TH2F(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins
           ,Int_t nbinsy,const Float_t *ybins)
           :TH2(name,title,nbinsx,xbins,nbinsy,ybins)
{
   TArrayF::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2F::TH2F(const TMatrixFBase &m)
:TH2("TMatrixFBase","",m.GetNcols(),m.GetColLwb(),1+m.GetColUpb(),m.GetNrows(),m.GetRowLwb(),1+m.GetRowUpb())
{
   TArrayF::Set(fNcells);
   Int_t ilow = m.GetRowLwb();
   Int_t iup  = m.GetRowUpb();
   Int_t jlow = m.GetColLwb();
   Int_t jup  = m.GetColUpb();
   for (Int_t i=ilow;i<=iup;i++) {
      for (Int_t j=jlow;j<=jup;j++) {
         SetBinContent(j-jlow+1,i-ilow+1,m(i,j));
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TH2F::TH2F(const TH2F &h2f) : TH2(), TArrayF()
{
   ((TH2F&)h2f).Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Copy.

void TH2F::Copy(TObject &newth2) const
{
   TH2::Copy((TH2F&)newth2);
}


////////////////////////////////////////////////////////////////////////////////
/// Reset this histogram: contents, errors, etc.

void TH2F::Reset(Option_t *option)
{
   TH2::Reset(option);
   TArrayF::Reset();
}


////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH2F::SetBinsLength(Int_t n)
{
   if (n < 0) n = (fXaxis.GetNbins()+2)*(fYaxis.GetNbins()+2);
   fNcells = n;
   TArrayF::Set(n);
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TH2F.

void TH2F::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TH2F::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      if (R__v < 2) {
         R__b.ReadVersion();
         TH1::Streamer(R__b);
         TArrayF::Streamer(R__b);
         R__b.ReadVersion();
         R__b >> fScalefactor;
         R__b >> fTsumwy;
         R__b >> fTsumwy2;
         R__b >> fTsumwxy;
      } else {
         TH2::Streamer(R__b);
         TArrayF::Streamer(R__b);
         R__b.CheckByteCount(R__s, R__c, TH2F::IsA());
      }
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TH2F::Class(),this);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH2F& TH2F::operator=(const TH2F &h1)
{
   if (this != &h1)  ((TH2F&)h1).Copy(*this);
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH2F operator*(Float_t c1, TH2F &h1)
{
   TH2F hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH2F operator*(TH2F &h1, Float_t c1)
{
   TH2F hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH2F operator+(TH2F &h1, TH2F &h2)
{
   TH2F hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator -

TH2F operator-(TH2F &h1, TH2F &h2)
{
   TH2F hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH2F operator*(TH2F &h1, TH2F &h2)
{
   TH2F hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH2F operator/(TH2F &h1, TH2F &h2)
{
   TH2F hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


//______________________________________________________________________________
//                     TH2D methods
//  TH2D a 2-D histogram with eight bytes per cell (double)
//______________________________________________________________________________

ClassImp(TH2D);


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2D::TH2D(): TH2(), TArrayD()
{
   SetBinsLength(9);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH2D::~TH2D()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2D::TH2D(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,Double_t ylow,Double_t yup)
           :TH2(name,title,nbinsx,xlow,xup,nbinsy,ylow,yup)
{
   TArrayD::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();

   if (xlow >= xup || ylow >= yup) SetBuffer(fgBufferSize);
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2D::TH2D(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,Double_t ylow,Double_t yup)
           :TH2(name,title,nbinsx,xbins,nbinsy,ylow,yup)
{
   TArrayD::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2D::TH2D(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,const Double_t *ybins)
           :TH2(name,title,nbinsx,xlow,xup,nbinsy,ybins)
{
   TArrayD::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2D::TH2D(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,const Double_t *ybins)
           :TH2(name,title,nbinsx,xbins,nbinsy,ybins)
{
   TArrayD::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2D::TH2D(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins
           ,Int_t nbinsy,const Float_t *ybins)
           :TH2(name,title,nbinsx,xbins,nbinsy,ybins)
{
   TArrayD::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2D::TH2D(const TMatrixDBase &m)
:TH2("TMatrixDBase","",m.GetNcols(),m.GetColLwb(),1+m.GetColUpb(),m.GetNrows(),m.GetRowLwb(),1+m.GetRowUpb())
{
   TArrayD::Set(fNcells);
   Int_t ilow = m.GetRowLwb();
   Int_t iup  = m.GetRowUpb();
   Int_t jlow = m.GetColLwb();
   Int_t jup  = m.GetColUpb();
   for (Int_t i=ilow;i<=iup;i++) {
      for (Int_t j=jlow;j<=jup;j++) {
         SetBinContent(j-jlow+1,i-ilow+1,m(i,j));
      }
   }
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TH2D::TH2D(const TH2D &h2d) : TH2(), TArrayD()
{
   ((TH2D&)h2d).Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Copy.

void TH2D::Copy(TObject &newth2) const
{
   TH2::Copy((TH2D&)newth2);
}


////////////////////////////////////////////////////////////////////////////////
/// Reset this histogram: contents, errors, etc.

void TH2D::Reset(Option_t *option)
{
   TH2::Reset(option);
   TArrayD::Reset();
}


////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH2D::SetBinsLength(Int_t n)
{
   if (n < 0) n = (fXaxis.GetNbins()+2)*(fYaxis.GetNbins()+2);
   fNcells = n;
   TArrayD::Set(n);
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TH2D.

void TH2D::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TH2D::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      if (R__v < 2) {
         R__b.ReadVersion();
         TH1::Streamer(R__b);
         TArrayD::Streamer(R__b);
         R__b.ReadVersion();
         R__b >> fScalefactor;
         R__b >> fTsumwy;
         R__b >> fTsumwy2;
         R__b >> fTsumwxy;
      } else {
         TH2::Streamer(R__b);
         TArrayD::Streamer(R__b);
         R__b.CheckByteCount(R__s, R__c, TH2D::IsA());
      }
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TH2D::Class(),this);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH2D& TH2D::operator=(const TH2D &h1)
{
   if (this != &h1)  ((TH2D&)h1).Copy(*this);
   return *this;
}



////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH2D operator*(Float_t c1, TH2D &h1)
{
   TH2D hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH2D operator+(TH2D &h1, TH2D &h2)
{
   TH2D hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator -

TH2D operator-(TH2D &h1, TH2D &h2)
{
   TH2D hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH2D operator*(TH2D &h1, TH2D &h2)
{
   TH2D hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH2D operator/(TH2D &h1, TH2D &h2)
{
   TH2D hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}
