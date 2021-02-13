// @(#)root/hist:$Id$
// Author: Rene Brun   27/10/95

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
#include "TH3.h"
#include "TProfile2D.h"
#include "TH2.h"
#include "TF3.h"
#include "TVirtualPad.h"
#include "TVirtualHistPainter.h"
#include "THLimitsFinder.h"
#include "TRandom.h"
#include "TError.h"
#include "TMath.h"
#include "TObjString.h"

ClassImp(TH3);

/** \addtogroup Hist
@{
\class TH3C
\brief 3-D histogram with a byte per channel (see TH1 documentation)
\class TH3S
\brief 3-D histogram with a short per channel (see TH1 documentation)
\class TH3I
\brief 3-D histogram with an int per channel (see TH1 documentation)}
\class TH3F
\brief 3-D histogram with a float per channel (see TH1 documentation)}
\class TH3D
\brief 3-D histogram with a double per channel (see TH1 documentation)}
@}
*/

/** \class TH3
    \ingroup Hist
The 3-D histogram classes derived from the 1-D histogram classes.
All operations are supported (fill, fit).
Drawing is currently restricted to one single option.
A cloud of points is drawn. The number of points is proportional to
cell content.

-   TH3C a 3-D histogram with one byte per cell (char)
-   TH3S a 3-D histogram with two bytes per cell (short integer)
-   TH3I a 3-D histogram with four bytes per cell (32 bits integer)
-   TH3F a 3-D histogram with four bytes per cell (float)
-   TH3D a 3-D histogram with eight bytes per cell (double)
*/


////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TH3::TH3()
{
   fDimension   = 3;
   fTsumwy      = fTsumwy2 = fTsumwxy = 0;
   fTsumwz      = fTsumwz2 = fTsumwxz = fTsumwyz = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for fix bin size 3-D histograms.

TH3::TH3(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
                                     ,Int_t nbinsy,Double_t ylow,Double_t yup
                                     ,Int_t nbinsz,Double_t zlow,Double_t zup)
     :TH1(name,title,nbinsx,xlow,xup),
      TAtt3D()
{
   fDimension   = 3;
   if (nbinsy <= 0) {
      Warning("TH3","nbinsy is <=0 - set to nbinsy = 1");
      nbinsy = 1;
   }
   if (nbinsz <= 0) {
      Warning("TH3","nbinsz is <=0 - set to nbinsz = 1");
      nbinsz = 1;
   }
   fYaxis.Set(nbinsy,ylow,yup);
   fZaxis.Set(nbinsz,zlow,zup);
   fNcells      = (nbinsx+2)*(nbinsy+2)*(nbinsz+2);
   fTsumwy      = fTsumwy2 = fTsumwxy = 0;
   fTsumwz      = fTsumwz2 = fTsumwxz = fTsumwyz = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for variable bin size 3-D histograms.

TH3::TH3(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins
                                           ,Int_t nbinsy,const Float_t *ybins
                                           ,Int_t nbinsz,const Float_t *zbins)
     :TH1(name,title,nbinsx,xbins),
      TAtt3D()
{
   fDimension   = 3;
   if (nbinsy <= 0) {Warning("TH3","nbinsy is <=0 - set to nbinsy = 1"); nbinsy = 1; }
   if (nbinsz <= 0) nbinsz = 1;
   if (ybins) fYaxis.Set(nbinsy,ybins);
   else       fYaxis.Set(nbinsy,0,1);
   if (zbins) fZaxis.Set(nbinsz,zbins);
   else       fZaxis.Set(nbinsz,0,1);
   fNcells      = (nbinsx+2)*(nbinsy+2)*(nbinsz+2);
   fTsumwy      = fTsumwy2 = fTsumwxy = 0;
   fTsumwz      = fTsumwz2 = fTsumwxz = fTsumwyz = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for variable bin size 3-D histograms.

TH3::TH3(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
                                           ,Int_t nbinsy,const Double_t *ybins
                                           ,Int_t nbinsz,const Double_t *zbins)
     :TH1(name,title,nbinsx,xbins),
      TAtt3D()
{
   fDimension   = 3;
   if (nbinsy <= 0) {Warning("TH3","nbinsy is <=0 - set to nbinsy = 1"); nbinsy = 1; }
   if (nbinsz <= 0) nbinsz = 1;
   if (ybins) fYaxis.Set(nbinsy,ybins);
   else       fYaxis.Set(nbinsy,0,1);
   if (zbins) fZaxis.Set(nbinsz,zbins);
   else       fZaxis.Set(nbinsz,0,1);
   fNcells      = (nbinsx+2)*(nbinsy+2)*(nbinsz+2);
   fTsumwy      = fTsumwy2 = fTsumwxy = 0;
   fTsumwz      = fTsumwz2 = fTsumwxz = fTsumwyz = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.
/// The list of functions is not copied. (Use Clone if needed)

TH3::TH3(const TH3 &h) : TH1(), TAtt3D()
{
   ((TH3&)h).Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH3::~TH3()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Copy.

void TH3::Copy(TObject &obj) const
{
   TH1::Copy(obj);
   ((TH3&)obj).fTsumwy      = fTsumwy;
   ((TH3&)obj).fTsumwy2     = fTsumwy2;
   ((TH3&)obj).fTsumwxy     = fTsumwxy;
   ((TH3&)obj).fTsumwz      = fTsumwz;
   ((TH3&)obj).fTsumwz2     = fTsumwz2;
   ((TH3&)obj).fTsumwxz     = fTsumwxz;
   ((TH3&)obj).fTsumwyz     = fTsumwyz;
}


////////////////////////////////////////////////////////////////////////////////
/// Fill histogram with all entries in the buffer.
/// action = -1 histogram is reset and refilled from the buffer (called by THistPainter::Paint)
/// action =  0 histogram is filled from the buffer
/// action =  1 histogram is filled and buffer is deleted
///             The buffer is automatically deleted when the number of entries
///             in the buffer is greater than the number of entries in the histogram

Int_t TH3::BufferEmpty(Int_t action)
{
   // do we need to compute the bin size?
   if (!fBuffer) return 0;
   Int_t nbentries = (Int_t)fBuffer[0];
   if (!nbentries) return 0;
   Double_t *buffer = fBuffer;
   if (nbentries < 0) {
      if (action == 0) return 0;
      nbentries  = -nbentries;
      fBuffer=0;
      Reset("ICES");
      fBuffer = buffer;
   }
   if (CanExtendAllAxes() || fXaxis.GetXmax() <= fXaxis.GetXmin() ||
      fYaxis.GetXmax() <= fYaxis.GetXmin() ||
      fZaxis.GetXmax() <= fZaxis.GetXmin()) {
         //find min, max of entries in buffer
         Double_t xmin = fBuffer[2];
         Double_t xmax = xmin;
         Double_t ymin = fBuffer[3];
         Double_t ymax = ymin;
         Double_t zmin = fBuffer[4];
         Double_t zmax = zmin;
         for (Int_t i=1;i<nbentries;i++) {
            Double_t x = fBuffer[4*i+2];
            if (x < xmin) xmin = x;
            if (x > xmax) xmax = x;
            Double_t y = fBuffer[4*i+3];
            if (y < ymin) ymin = y;
            if (y > ymax) ymax = y;
            Double_t z = fBuffer[4*i+4];
            if (z < zmin) zmin = z;
            if (z > zmax) zmax = z;
         }
         if (fXaxis.GetXmax() <= fXaxis.GetXmin() || fYaxis.GetXmax() <= fYaxis.GetXmin() || fZaxis.GetXmax() <= fZaxis.GetXmin()) {
            THLimitsFinder::GetLimitsFinder()->FindGoodLimits(this,xmin,xmax,ymin,ymax,zmin,zmax);
         } else {
            fBuffer = 0;
            Int_t keep = fBufferSize; fBufferSize = 0;
            if (xmin <  fXaxis.GetXmin()) ExtendAxis(xmin,&fXaxis);
            if (xmax >= fXaxis.GetXmax()) ExtendAxis(xmax,&fXaxis);
            if (ymin <  fYaxis.GetXmin()) ExtendAxis(ymin,&fYaxis);
            if (ymax >= fYaxis.GetXmax()) ExtendAxis(ymax,&fYaxis);
            if (zmin <  fZaxis.GetXmin()) ExtendAxis(zmin,&fZaxis);
            if (zmax >= fZaxis.GetXmax()) ExtendAxis(zmax,&fZaxis);
            fBuffer = buffer;
            fBufferSize = keep;
         }
   }
   fBuffer = 0;

   for (Int_t i=0;i<nbentries;i++) {
      Fill(buffer[4*i+2],buffer[4*i+3],buffer[4*i+4],buffer[4*i+1]);
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
/// fBuffer[0] = number of entries in buffer
/// fBuffer[1] = w of first entry
/// fBuffer[2] = x of first entry
/// fBuffer[3] = y of first entry
/// fBuffer[4] = z of first entry

Int_t TH3::BufferFill(Double_t x, Double_t y, Double_t z, Double_t w)
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
   if (4*nbentries+4 >= fBufferSize) {
      BufferEmpty(1);
      return Fill(x,y,z,w);
   }
   fBuffer[4*nbentries+1] = w;
   fBuffer[4*nbentries+2] = x;
   fBuffer[4*nbentries+3] = y;
   fBuffer[4*nbentries+4] = z;
   fBuffer[0] += 1;
   return -3;
}


////////////////////////////////////////////////////////////////////////////////
/// Invalid Fill method

Int_t TH3::Fill(Double_t )
{
   Error("Fill", "Invalid signature - do nothing");
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment cell defined by x,y,z by 1 .
///
/// The function returns the corresponding global bin number which has its content
/// incremented by 1

Int_t TH3::Fill(Double_t x, Double_t y, Double_t z)
{
   if (fBuffer) return BufferFill(x,y,z,1);

   Int_t binx, biny, binz, bin;
   fEntries++;
   binx = fXaxis.FindBin(x);
   biny = fYaxis.FindBin(y);
   binz = fZaxis.FindBin(z);
   if (binx <0 || biny <0 || binz<0) return -1;
   bin  =  binx + (fXaxis.GetNbins()+2)*(biny + (fYaxis.GetNbins()+2)*binz);
   if (fSumw2.fN) ++fSumw2.fArray[bin];
   AddBinContent(bin);
   if (binx == 0 || binx > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }

   if (biny == 0 || biny > fYaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (binz == 0 || binz > fZaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   ++fTsumw;
   ++fTsumw2;
   fTsumwx  += x;
   fTsumwx2 += x*x;
   fTsumwy  += y;
   fTsumwy2 += y*y;
   fTsumwxy += x*y;
   fTsumwz  += z;
   fTsumwz2 += z*z;
   fTsumwxz += x*z;
   fTsumwyz += y*z;
   return bin;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment cell defined by x,y,z by a weight w.
///
/// If the weight is not equal to 1, the storage of the sum of squares of
///  weights is automatically triggered and the sum of the squares of weights is incremented
///  by w^2 in the cell corresponding to x,y,z.
///
/// The function returns the corresponding global bin number which has its content
/// incremented by w

Int_t TH3::Fill(Double_t x, Double_t y, Double_t z, Double_t w)
{
   if (fBuffer) return BufferFill(x,y,z,w);

   Int_t binx, biny, binz, bin;
   fEntries++;
   binx = fXaxis.FindBin(x);
   biny = fYaxis.FindBin(y);
   binz = fZaxis.FindBin(z);
   if (binx <0 || biny <0 || binz<0) return -1;
   bin  =  binx + (fXaxis.GetNbins()+2)*(biny + (fYaxis.GetNbins()+2)*binz);
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();   // must be called before AddBinContent
   if (fSumw2.fN) fSumw2.fArray[bin] += w*w;
   AddBinContent(bin,w);
   if (binx == 0 || binx > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (biny == 0 || biny > fYaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (binz == 0 || binz > fZaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   fTsumw   += w;
   fTsumw2  += w*w;
   fTsumwx  += w*x;
   fTsumwx2 += w*x*x;
   fTsumwy  += w*y;
   fTsumwy2 += w*y*y;
   fTsumwxy += w*x*y;
   fTsumwz  += w*z;
   fTsumwz2 += w*z*z;
   fTsumwxz += w*x*z;
   fTsumwyz += w*y*z;
   return bin;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment cell defined by namex,namey,namez by a weight w
///
/// If the weight is not equal to 1, the storage of the sum of squares of
///  weights is automatically triggered and the sum of the squares of weights is incremented
///  by w^2 in the corresponding cell.
/// The function returns the corresponding global bin number which has its content
/// incremented by w

Int_t TH3::Fill(const char *namex, const char *namey, const char *namez, Double_t w)
{
   Int_t binx, biny, binz, bin;
   fEntries++;
   binx = fXaxis.FindBin(namex);
   biny = fYaxis.FindBin(namey);
   binz = fZaxis.FindBin(namez);
   if (binx <0 || biny <0 || binz<0) return -1;
   bin  =  binx + (fXaxis.GetNbins()+2)*(biny + (fYaxis.GetNbins()+2)*binz);
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();   // must be called before AddBinContent
   if (fSumw2.fN) fSumw2.fArray[bin] += w*w;
   AddBinContent(bin,w);
   if (binx == 0 || binx > fXaxis.GetNbins()) return -1;
   if (biny == 0 || biny > fYaxis.GetNbins()) return -1;
   if (binz == 0 || binz > fZaxis.GetNbins()) return -1;

   Double_t v = w;
   fTsumw   += v;
   fTsumw2  += v*v;
   // skip computation of the statistics along axis that have labels (can be extended and are aphanumeric)
   UInt_t labelBitMask = GetAxisLabelStatus();
   if (labelBitMask != TH1::kAllAxes) {
      Double_t x = (labelBitMask & TH1::kXaxis) ? 0 : fXaxis.GetBinCenter(binx);
      Double_t y = (labelBitMask & TH1::kYaxis) ? 0 : fYaxis.GetBinCenter(biny);
      Double_t z = (labelBitMask & TH1::kZaxis) ? 0 : fZaxis.GetBinCenter(binz);
      fTsumwx += v * x;
      fTsumwx2 += v * x * x;
      fTsumwy += v * y;
      fTsumwy2 += v * y * y;
      fTsumwxy += v * x * y;
      fTsumwz += v * z;
      fTsumwz2 += v * z * z;
      fTsumwxz += v * x * z;
      fTsumwyz += v * y * z;
   }
   return bin;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment cell defined by namex,y,namez by a weight w
///
/// If the weight is not equal to 1, the storage of the sum of squares of
///  weights is automatically triggered and the sum of the squares of weights is incremented
///  by w^2 in the corresponding cell.
/// The function returns the corresponding global bin number which has its content
/// incremented by w

Int_t TH3::Fill(const char *namex, Double_t y, const char *namez, Double_t w)
{
   Int_t binx, biny, binz, bin;
   fEntries++;
   binx = fXaxis.FindBin(namex);
   biny = fYaxis.FindBin(y);
   binz = fZaxis.FindBin(namez);
   if (binx <0 || biny <0 || binz<0) return -1;
   bin  =  binx + (fXaxis.GetNbins()+2)*(biny + (fYaxis.GetNbins()+2)*binz);
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();   // must be called before AddBinContent
   if (fSumw2.fN) fSumw2.fArray[bin] += w*w;
   AddBinContent(bin,w);
   if (binx == 0 || binx > fXaxis.GetNbins()) return -1;
   if (biny == 0 || biny > fYaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (binz == 0 || binz > fZaxis.GetNbins()) return -1;
   Double_t v = w;
   fTsumw   += v;
   fTsumw2  += v*v;
   fTsumwy  += v*y;
   fTsumwy2 += v*y*y;
   // skip computation of the statistics along axis that have labels (can be extended and are aphanumeric)
   UInt_t labelBitMask = GetAxisLabelStatus();
   if (labelBitMask != (TH1::kXaxis | TH1::kZaxis) ) {
      Double_t x = (labelBitMask & TH1::kXaxis) ? 0 : fXaxis.GetBinCenter(binx);
      Double_t z = (labelBitMask & TH1::kZaxis) ? 0 : fZaxis.GetBinCenter(binz);
      fTsumwx += v * x;
      fTsumwx2 += v * x * x;
      fTsumwxy += v * x * y;
      fTsumwz += v * z;
      fTsumwz2 += v * z * z;
      fTsumwxz += v * x * z;
      fTsumwyz += v * y * z;
   }
   return bin;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment cell defined by namex,namey,z by a weight w
///
/// If the weight is not equal to 1, the storage of the sum of squares of
///  weights is automatically triggered and the sum of the squares of weights is incremented
///  by w^2 in the corresponding cell.
/// The function returns the corresponding global bin number which has its content
/// incremented by w

Int_t TH3::Fill(const char *namex, const char *namey, Double_t z, Double_t w)
{
   Int_t binx, biny, binz, bin;
   fEntries++;
   binx = fXaxis.FindBin(namex);
   biny = fYaxis.FindBin(namey);
   binz = fZaxis.FindBin(z);
   if (binx <0 || biny <0 || binz<0) return -1;
   bin  =  binx + (fXaxis.GetNbins()+2)*(biny + (fYaxis.GetNbins()+2)*binz);
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();   // must be called before AddBinContent
   if (fSumw2.fN) fSumw2.fArray[bin] += w*w;
   AddBinContent(bin,w);
   if (binx == 0 || binx > fXaxis.GetNbins()) return -1;
   if (biny == 0 || biny > fYaxis.GetNbins()) return -1;
   if (binz == 0 || binz > fZaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   Double_t v = w;
   fTsumw   += v;
   fTsumw2  += v*v;
   fTsumwz  += v*z;
   fTsumwz2 += v*z*z;
   // skip computation of the statistics along axis that have labels (can be extended and are aphanumeric)
   UInt_t labelBitMask = GetAxisLabelStatus();
   if (labelBitMask != (TH1::kXaxis | TH1::kYaxis)) {
      Double_t x = (labelBitMask & TH1::kXaxis) ? 0 : fXaxis.GetBinCenter(binx);
      Double_t y = (labelBitMask & TH1::kYaxis) ? 0 : fYaxis.GetBinCenter(biny);
      fTsumwx += v * x;
      fTsumwx2 += v * x * x;
      fTsumwy += v * y;
      fTsumwy2 += v * y * y;
      fTsumwxy += v * x * y;
      fTsumwxz += v * x * z;
      fTsumwyz += v * y * z;
   }
   return bin;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment cell defined by x,namey,namez by a weight w
///
/// If the weight is not equal to 1, the storage of the sum of squares of
///  weights is automatically triggered and the sum of the squares of weights is incremented
///  by w^2 in the corresponding cell.
/// The function returns the corresponding global bin number which has its content
/// incremented by w

Int_t TH3::Fill(Double_t x, const char *namey, const char *namez, Double_t w)
{
   Int_t binx, biny, binz, bin;
   fEntries++;
   binx = fXaxis.FindBin(x);
   biny = fYaxis.FindBin(namey);
   binz = fZaxis.FindBin(namez);
   if (binx <0 || biny <0 || binz<0) return -1;
   bin  =  binx + (fXaxis.GetNbins()+2)*(biny + (fYaxis.GetNbins()+2)*binz);
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();   // must be called before AddBinContent
   if (fSumw2.fN) fSumw2.fArray[bin] += w*w;
   AddBinContent(bin,w);
   if (binx == 0 || binx > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (biny == 0 || biny > fYaxis.GetNbins()) return -1;
   if (binz == 0 || binz > fZaxis.GetNbins()) return -1;

   // skip computation of the statistics along axis that have labels (can be extended and are aphanumeric)
   UInt_t labelBitMask = GetAxisLabelStatus();
   Double_t y = (labelBitMask & TH1::kYaxis) ? 0 : fYaxis.GetBinCenter(biny);
   Double_t z = (labelBitMask & TH1::kZaxis) ? 0 : fZaxis.GetBinCenter(binz);
   Double_t v = w;
   fTsumw += v;
   fTsumw2 += v * v;
   fTsumwx += v * x;
   fTsumwx2 += v * x * x;
   fTsumwy += v * y;
   fTsumwy2 += v * y * y;
   fTsumwxy += v * x * y;
   fTsumwz += v * z;
   fTsumwz2 += v * z * z;
   fTsumwxz += v * x * z;
   fTsumwyz += v * y * z;
   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Increment cell defined by namex , y ,z by a weight w
///
/// If the weight is not equal to 1, the storage of the sum of squares of
///  weights is automatically triggered and the sum of the squares of weights is incremented
///  by w^2 in the corresponding cell.
/// The function returns the corresponding global bin number which has its content
/// incremented by w

Int_t TH3::Fill(const char * namex, Double_t y, Double_t z, Double_t w)
{
   Int_t binx, biny, binz, bin;
   fEntries++;
   binx = fXaxis.FindBin(namex);
   biny = fYaxis.FindBin(y);
   binz = fZaxis.FindBin(z);
   if (binx < 0 || biny < 0 || binz < 0)
      return -1;
   bin = binx + (fXaxis.GetNbins() + 2) * (biny + (fYaxis.GetNbins() + 2) * binz);
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW))
      Sumw2(); // must be called before AddBinContent
   if (fSumw2.fN)
      fSumw2.fArray[bin] += w * w;
   AddBinContent(bin, w);
   if (binx == 0 || binx > fXaxis.GetNbins()) {
         return -1;
   }
   if (biny == 0 || biny > fYaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour())
         return -1;
   }
   if (binz == 0 || binz > fZaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour())
         return -1;
   }
   UInt_t labelBitMask = GetAxisLabelStatus();
   Double_t x = (labelBitMask & TH1::kXaxis) ? 0 : fXaxis.GetBinCenter(binx);
   Double_t v = w;
   fTsumw += v;
   fTsumw2 += v * v;
   fTsumwx += v * x;
   fTsumwx2 += v * x * x;
   fTsumwy += v * y;
   fTsumwy2 += v * y * y;
   fTsumwxy += v * x * y;
   fTsumwz += v * z;
   fTsumwz2 += v * z * z;
   fTsumwxz += v * x * z;
   fTsumwyz += v * y * z;
   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Increment cell defined by x,namey,z by a weight w
///
/// If the weight is not equal to 1, the storage of the sum of squares of
///  weights is automatically triggered and the sum of the squares of weights is incremented
///  by w^2 in the corresponding cell.
/// The function returns the corresponding global bin number which has its content
/// incremented by w

Int_t TH3::Fill(Double_t x, const char *namey, Double_t z, Double_t w)
{
   Int_t binx, biny, binz, bin;
   fEntries++;
   binx = fXaxis.FindBin(x);
   biny = fYaxis.FindBin(namey);
   binz = fZaxis.FindBin(z);
   if (binx <0 || biny <0 || binz<0) return -1;
   bin  =  binx + (fXaxis.GetNbins()+2)*(biny + (fYaxis.GetNbins()+2)*binz);
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();   // must be called before AddBinContent
   if (fSumw2.fN) fSumw2.fArray[bin] += w*w;
   AddBinContent(bin,w);
   if (binx == 0 || binx > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (biny == 0 || biny > fYaxis.GetNbins()) return -1;
   if (binz == 0 || binz > fZaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   UInt_t labelBitMask = GetAxisLabelStatus();
   Double_t y = (labelBitMask & TH1::kYaxis) ? 0 : fYaxis.GetBinCenter(biny);
   Double_t v = w;
   fTsumw   += v;
   fTsumw2  += v*v;
   fTsumwx  += v*x;
   fTsumwx2 += v*x*x;
   fTsumwy  += v*y;
   fTsumwy2 += v*y*y;
   fTsumwxy += v*x*y;
   fTsumwz  += v*z;
   fTsumwz2 += v*z*z;
   fTsumwxz += v*x*z;
   fTsumwyz += v*y*z;
   return bin;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment cell defined by x,y,namez by a weight w
///
/// If the weight is not equal to 1, the storage of the sum of squares of
///  weights is automatically triggered and the sum of the squares of weights is incremented
///  by w^2 in the corresponding cell.
/// The function returns the corresponding global bin number which has its content
/// incremented by w

Int_t TH3::Fill(Double_t x, Double_t y, const char *namez, Double_t w)
{
   Int_t binx, biny, binz, bin;
   fEntries++;
   binx = fXaxis.FindBin(x);
   biny = fYaxis.FindBin(y);
   binz = fZaxis.FindBin(namez);
   if (binx <0 || biny <0 || binz<0) return -1;
   bin  =  binx + (fXaxis.GetNbins()+2)*(biny + (fYaxis.GetNbins()+2)*binz);
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();   // must be called before AddBinContent
   if (fSumw2.fN) fSumw2.fArray[bin] += w*w;
   AddBinContent(bin,w);
   if (binx == 0 || binx > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (biny == 0 || biny > fYaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (binz == 0 || binz > fZaxis.GetNbins()) return -1;
   UInt_t labelBitMask = GetAxisLabelStatus();
   Double_t z = (labelBitMask & TH1::kZaxis) ? 0 : fZaxis.GetBinCenter(binz);
   Double_t v = w;
   fTsumw   += v;
   fTsumw2  += v*v;
   fTsumwx  += v*x;
   fTsumwx2 += v*x*x;
   fTsumwy  += v*y;
   fTsumwy2 += v*y*y;
   fTsumwxy += v*x*y;
   fTsumwz  += v*z;
   fTsumwz2 += v*z*z;
   fTsumwxz += v*x*z;
   fTsumwyz += v*y*z;
   return bin;
}


////////////////////////////////////////////////////////////////////////////////
/// Fill histogram following distribution in function fname.
///
///  @param fname  : Function name used for filling the historam
///  @param ntimes : number of times the histogram is filled
///  @param rng    : (optional) Random number generator used to sample
///
///   The distribution contained in the function fname (TF1) is integrated
///   over the channel contents.
///   It is normalized to 1.
///   Getting one random number implies:
///     - Generating a random number between 0 and 1 (say r1)
///     - Look in which bin in the normalized integral r1 corresponds to
///     - Fill histogram channel
///   ntimes random numbers are generated
///
/// N.B. By dfault this methods approximates the integral of the function in each bin with the
///      function value at the center of the bin, mutiplied by the bin width
///
///
///  One can also call TF1::GetRandom to get a random variate from a function.

void TH3::FillRandom(const char *fname, Int_t ntimes, TRandom * rng)
{
   Int_t bin, binx, biny, binz, ibin, loop;
   Double_t r1, x, y,z, xv[3];
   //  Search for fname in the list of ROOT defined functions
   TObject *fobj = gROOT->GetFunction(fname);
   if (!fobj) { Error("FillRandom", "Unknown function: %s",fname); return; }
   TF3 *f1 = dynamic_cast<TF3*>( fobj );
   if (!f1) { Error("FillRandom", "Function: %s is not a TF3, is a %s",fname,fobj->IsA()->GetName()); return; }

   TAxis & xAxis = fXaxis;
   TAxis & yAxis = fYaxis;
   TAxis & zAxis = fZaxis;

   // in case axes of histogram are not defined use the function axis
   if (fXaxis.GetXmax() <= fXaxis.GetXmin()  || fYaxis.GetXmax() <= fYaxis.GetXmin() || fZaxis.GetXmax() <= fZaxis.GetXmin() ) {
      Double_t xmin,xmax,ymin,ymax,zmin,zmax;
      f1->GetRange(xmin,ymin,zmin,xmax,ymax,zmax);
      Info("FillRandom","Using function axis and range ([%g,%g],[%g,%g],[%g,%g])",xmin, xmax,ymin,ymax,zmin,zmax);
      xAxis = *(f1->GetHistogram()->GetXaxis());
      yAxis = *(f1->GetHistogram()->GetYaxis());
      zAxis = *(f1->GetHistogram()->GetZaxis());
   }

   //  Allocate temporary space to store the integral and compute integral
   Int_t nbinsx = xAxis.GetNbins();
   Int_t nbinsy = yAxis.GetNbins();
   Int_t nbinsz = zAxis.GetNbins();
   Int_t nxy = nbinsx*nbinsy;
   Int_t nbins  = nbinsx*nbinsy*nbinsz;

   Double_t *integral = new Double_t[nbins+1];
   ibin = 0;
   integral[ibin] = 0;
   // approximate integral with function value at bin center
   for (binz=1;binz<=nbinsz;binz++) {
      xv[2] = zAxis.GetBinCenter(binz);
      for (biny=1;biny<=nbinsy;biny++) {
         xv[1] = yAxis.GetBinCenter(biny);
         for (binx=1;binx<=nbinsx;binx++) {
            xv[0] = xAxis.GetBinCenter(binx);
            ibin++;
            Double_t fint = f1->EvalPar(xv, nullptr);
            // uncomment this line to have the integral computation in a bin
            // Double_t fint = f1->Integral(xAxis.GetBinLowEdge(binx), xAxis.GetBinUpEdge(binx),
            //                              yAxis.GetBinLowEdge(biny), yAxis.GetBinUpEdge(biny),
            //                              zAxis.GetBinLowEdge(binz), zAxis.GetBinUpEdge(binz));
            integral[ibin] = integral[ibin-1] + fint;
         }
      }
   }

   //  Normalize integral to 1
   if (integral[nbins] == 0 ) {
      delete [] integral;
      Error("FillRandom", "Integral = zero"); return;
   }
   for (bin=1;bin<=nbins;bin++)  integral[bin] /= integral[nbins];

   // Start main loop ntimes
   if (fDimension < 2) nbinsy = -1;
   if (fDimension < 3) nbinsz = -1;
   for (loop=0;loop<ntimes;loop++) {
      r1 = (rng) ? rng->Rndm() : gRandom->Rndm();
      ibin = TMath::BinarySearch(nbins,&integral[0],r1);
      binz = ibin/nxy;
      biny = (ibin - nxy*binz)/nbinsx;
      binx = 1 + ibin - nbinsx*(biny + nbinsy*binz);
      if (nbinsz) binz++;
      if (nbinsy) biny++;
      x    = xAxis.GetBinCenter(binx);
      y    = yAxis.GetBinCenter(biny);
      z    = zAxis.GetBinCenter(binz);
      Fill(x,y,z, 1.);
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
///   The distribution contained in the histogram h (TH3) is integrated
///   over the channel contents.
///   It is normalized to 1.
///   Getting one random number implies:
///     - Generating a random number between 0 and 1 (say r1)
///     - Look in which bin in the normalized integral r1 corresponds to
///     - Fill histogram channel
///   ntimes random numbers are generated

void TH3::FillRandom(TH1 *h, Int_t ntimes, TRandom * rng)
{
   if (!h) { Error("FillRandom", "Null histogram"); return; }
   if (fDimension != h->GetDimension()) {
      Error("FillRandom", "Histograms with different dimensions"); return;
   }

   if (h->ComputeIntegral() == 0) return;

   TH3 *h3 = (TH3*)h;
   Int_t loop;
   Double_t x,y,z;
   for (loop=0;loop<ntimes;loop++) {
      h3->GetRandom3(x,y,z,rng);
      Fill(x,y,z);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Project slices along Z in case of a 3-D histogram, then fit each slice
/// with function f1 and make a 2-d histogram for each fit parameter
/// Only cells in the bin range [binminx,binmaxx] and [binminy,binmaxy] are considered.
/// if f1=0, a gaussian is assumed
/// Before invoking this function, one can set a subrange to be fitted along Z
/// via f1->SetRange(zmin,zmax)
/// The argument option (default="QNR") can be used to change the fit options.
///     "Q" means Quiet mode
///     "N" means do not show the result of the fit
///     "R" means fit the function in the specified function range
///
/// Note that the generated histograms are added to the list of objects
/// in the current directory. It is the user's responsibility to delete
/// these histograms.
///
///  Example: Assume a 3-d histogram h3
///   Root > h3->FitSlicesZ(); produces 4 TH2D histograms
///          with h3_0 containing parameter 0(Constant) for a Gaus fit
///                    of each cell in X,Y projected along Z
///          with h3_1 containing parameter 1(Mean) for a gaus fit
///          with h3_2 containing parameter 2(StdDev)  for a gaus fit
///          with h3_chi2 containing the chisquare/number of degrees of freedom for a gaus fit
///
///   Root > h3->Fit(0,15,22,0,0,10);
///          same as above, but only for bins 15 to 22 along X
///          and only for cells in X,Y for which the corresponding projection
///          along Z has more than cut bins filled.
///
///  NOTE: To access the generated histograms in the current directory, do eg:
///     TH2D *h3_1 = (TH2D*)gDirectory->Get("h3_1");

void TH3::FitSlicesZ(TF1 *f1, Int_t binminx, Int_t binmaxx, Int_t binminy, Int_t binmaxy, Int_t cut, Option_t *option)
{
   //Int_t nbinsz  = fZaxis.GetNbins();

   // get correct first and last bins for outer axes used in the loop doing the slices
   // when using default values (0,-1) check if an axis range is set in outer axis
   // do same as in DoProjection for inner axis
   auto computeFirstAndLastBin = [](const TAxis  & outerAxis, Int_t &firstbin, Int_t &lastbin) {
      Int_t nbins  = outerAxis.GetNbins();
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
   };

   computeFirstAndLastBin(fXaxis, binminx, binmaxx);
   computeFirstAndLastBin(fYaxis, binminy, binmaxy);

   // limits for the axis of the fit results histograms are different
   auto computeAxisLimits = [](const TAxis  & outerAxis, Int_t firstbin, Int_t lastbin,
                               Int_t  &nBins, Double_t  &xMin, Double_t  & xMax) {
      Int_t firstOutBin = std::max(firstbin,1);
      Int_t lastOutBin = std::min(lastbin,outerAxis.GetNbins() ) ;
      nBins = lastOutBin-firstOutBin+1;
      xMin = outerAxis.GetBinLowEdge(firstOutBin);
      xMax = outerAxis.GetBinUpEdge(lastOutBin);
      // return first bin that is used in case of variable bin size axis
      return firstOutBin;
   };
   Int_t nbinsX = 0;
   Double_t xMin, xMax = 0;
   Int_t firstBinXaxis = computeAxisLimits(fXaxis, binminx, binmaxx, nbinsX, xMin, xMax);
   Int_t nbinsY = 0;
   Double_t yMin, yMax = 0;
   Int_t firstBinYaxis = computeAxisLimits(fYaxis, binminy, binmaxy, nbinsY, yMin, yMax);

   //default is to fit with a gaussian
   if (f1 == 0) {
      f1 = (TF1*)gROOT->GetFunction("gaus");
      if (f1 == 0) f1 = new TF1("gaus","gaus",fZaxis.GetXmin(),fZaxis.GetXmax());
      else         f1->SetRange(fZaxis.GetXmin(),fZaxis.GetXmax());
   }
   const char *fname = f1->GetName();
   Int_t npar = f1->GetNpar();
   Double_t *parsave = new Double_t[npar];
   f1->GetParameters(parsave);

   //Create one 2-d histogram for each function parameter
   Int_t ipar;
   TString name;
   TString title;
   std::vector<TH1*> hlist(npar+1); // include also chi2 histogram
   const TArrayD *xbins = fXaxis.GetXbins();
   const TArrayD *ybins = fYaxis.GetXbins();
   for (ipar=0;ipar<= npar;ipar++) {
      if (ipar < npar) {
         // fitted parameter histograms
         name = TString::Format("%s_%d",GetName(),ipar);
         title = TString::Format("Fitted value of par[%d]=%s",ipar,f1->GetParName(ipar));
      } else {
         // chi2 histograms
         name = TString::Format("%s_chi2",GetName());
         title = "chisquare";
      }
      if (xbins->fN == 0 && ybins->fN == 0) {
         hlist[ipar] = new TH2D(name, title,
                                nbinsX, xMin, xMax,
                                nbinsY, yMin, yMax);
      } else if (xbins->fN > 0 && ybins->fN > 0 ) {
         hlist[ipar] = new TH2D(name, title,
                                nbinsX, &xbins->fArray[firstBinXaxis],
                                nbinsY, &ybins->fArray[firstBinYaxis]);
      }
      // mixed case do not exist for TH3
      R__ASSERT(hlist[ipar]);

      hlist[ipar]->GetXaxis()->SetTitle(fXaxis.GetTitle());
      hlist[ipar]->GetYaxis()->SetTitle(fYaxis.GetTitle());
   }
   TH1 * hchi2 = hlist.back();

   //Loop on all cells in X,Y generate a projection along Z
   TH1D *hpz = nullptr;
   TString opt(option);
   // add option "N" when fitting the 2D histograms
   opt += " N ";

   for (Int_t biny=binminy; biny<=binmaxy; biny++) {
      for (Int_t binx=binminx; binx<=binmaxx; binx++) {
         // use TH3::ProjectionZ
         hpz = ProjectionZ("R_temp",binx,binx,biny,biny);

         Double_t nentries = hpz->GetEntries();
         if ( nentries <= 0 || nentries < cut) {
            if (!opt.Contains("Q"))
               Info("FitSlicesZ","Slice (%d,%d) skipped, the number of entries is zero or smaller than the given cut value, n=%f",binx,biny,nentries);
            continue;
         }
         f1->SetParameters(parsave);
         Int_t bin = hlist[0]->FindBin( fXaxis.GetBinCenter(binx), fYaxis.GetBinCenter(biny) );
         if (!opt.Contains("Q")) {
            int ibx,iby,ibz = 0;
            hlist[0]->GetBinXYZ(bin,ibx,iby,ibz);
            Info("DoFitSlices","Slice fit [(%f,%f),(%f,%f)]",hlist[0]->GetXaxis()->GetBinLowEdge(ibx), hlist[0]->GetXaxis()->GetBinUpEdge(ibx),
                                                            hlist[0]->GetYaxis()->GetBinLowEdge(iby), hlist[0]->GetYaxis()->GetBinUpEdge(iby));
         }
         hpz->Fit(fname,opt.Data());
         Int_t npfits = f1->GetNumberFitPoints();
         if (npfits > npar && npfits >= cut) {
            for (ipar=0;ipar<npar;ipar++) {
               hlist[ipar]->SetBinContent(bin,f1->GetParameter(ipar));
               hlist[ipar]->SetBinError(bin,f1->GetParError(ipar));
            }
            hchi2->SetBinContent(bin,f1->GetChisquare()/(npfits-npar));
         }
         else {
            if (!opt.Contains("Q"))
               Info("FitSlicesZ","Fitted slice (%d,%d) skipped, the number of fitted points is too small, n=%d",binx,biny,npfits);
         }
      }
   }
   delete [] parsave;
   delete hpz;
}


////////////////////////////////////////////////////////////////////////////////
/// See comments in TH1::GetBin

Int_t TH3::GetBin(Int_t binx, Int_t biny, Int_t binz) const
{
   Int_t ofy = fYaxis.GetNbins() + 1; // code duplication unavoidable because TH3 does not inherit from TH2
   if (biny < 0) biny = 0;
   if (biny > ofy) biny = ofy;

   Int_t ofz = fZaxis.GetNbins() + 1; // overflow bin
   if (binz < 0) binz = 0;
   if (binz > ofz) binz = ofz;

   return TH1::GetBin(binx) + (fXaxis.GetNbins() + 2) * (biny + (fYaxis.GetNbins() + 2) * binz);
}


////////////////////////////////////////////////////////////////////////////////
/// Compute first cell (binx,biny,binz) in the range [firstx,lastx](firsty,lasty][firstz,lastz] for which
/// diff = abs(cell_content-c) <= maxdiff
/// In case several cells in the specified range with diff=0 are found
/// the first cell found is returned in binx,biny,binz.
/// In case several cells in the specified range satisfy diff <=maxdiff
/// the cell with the smallest difference is returned in binx,biny,binz.
/// In all cases the function returns the smallest difference.
///
/// NOTE1: if firstx <= 0, firstx is set to bin 1
///        if (lastx < firstx then firstx is set to the number of bins in X
///        ie if firstx=0 and lastx=0 (default) the search is on all bins in X.
///        if firsty <= 0, firsty is set to bin 1
///        if (lasty < firsty then firsty is set to the number of bins in Y
///        ie if firsty=0 and lasty=0 (default) the search is on all bins in Y.
///        if firstz <= 0, firstz is set to bin 1
///        if (lastz < firstz then firstz is set to the number of bins in Z
///        ie if firstz=0 and lastz=0 (default) the search is on all bins in Z.
/// NOTE2: if maxdiff=0 (default), the first cell with content=c is returned.

Double_t TH3::GetBinWithContent3(Double_t c, Int_t &binx, Int_t &biny, Int_t &binz,
                                 Int_t firstx, Int_t lastx,
                                 Int_t firsty, Int_t lasty,
                                 Int_t firstz, Int_t lastz,
                                 Double_t maxdiff) const
{
   if (fDimension != 3) {
      binx = 0;
      biny = 0;
      binz = 0;
      Error("GetBinWithContent3","function is only valid for 3-D histograms");
      return 0;
   }
   if (firstx <= 0) firstx = 1;
   if (lastx < firstx) lastx = fXaxis.GetNbins();
   if (firsty <= 0) firsty = 1;
   if (lasty < firsty) lasty = fYaxis.GetNbins();
   if (firstz <= 0) firstz = 1;
   if (lastz < firstz) lastz = fZaxis.GetNbins();
   Int_t binminx = 0, binminy=0, binminz=0;
   Double_t diff, curmax = 1.e240;
   for (Int_t k=firstz;k<=lastz;k++) {
      for (Int_t j=firsty;j<=lasty;j++) {
         for (Int_t i=firstx;i<=lastx;i++) {
            diff = TMath::Abs(GetBinContent(i,j,k)-c);
            if (diff <= 0) {binx = i; biny=j; binz=k; return diff;}
            if (diff < curmax && diff <= maxdiff) {curmax = diff, binminx=i; binminy=j;binminz=k;}
         }
      }
   }
   binx = binminx;
   biny = binminy;
   binz = binminz;
   return curmax;
}


////////////////////////////////////////////////////////////////////////////////
/// Return correlation factor between axis1 and axis2.

Double_t TH3::GetCorrelationFactor(Int_t axis1, Int_t axis2) const
{
   if (axis1 < 1 || axis2 < 1 || axis1 > 3 || axis2 > 3) {
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

Double_t TH3::GetCovariance(Int_t axis1, Int_t axis2) const
{
   if (axis1 < 1 || axis2 < 1 || axis1 > 3 || axis2 > 3) {
      Error("GetCovariance","Wrong parameters");
      return 0;
   }
   Double_t stats[kNstat];
   GetStats(stats);
   Double_t sumw   = stats[0];
   Double_t sumw2  = stats[1];
   Double_t sumwx  = stats[2];
   Double_t sumwx2 = stats[3];
   Double_t sumwy  = stats[4];
   Double_t sumwy2 = stats[5];
   Double_t sumwxy = stats[6];
   Double_t sumwz  = stats[7];
   Double_t sumwz2 = stats[8];
   Double_t sumwxz = stats[9];
   Double_t sumwyz = stats[10];

   if (sumw == 0) return 0;
   if (axis1 == 1 && axis2 == 1) {
      return TMath::Abs(sumwx2/sumw - sumwx*sumwx/sumw2);
   }
   if (axis1 == 2 && axis2 == 2) {
      return TMath::Abs(sumwy2/sumw - sumwy*sumwy/sumw2);
   }
   if (axis1 == 3 && axis2 == 3) {
      return TMath::Abs(sumwz2/sumw - sumwz*sumwz/sumw2);
   }
   if ((axis1 == 1 && axis2 == 2) || (axis1 == 2 && axis2 == 1)) {
      return sumwxy/sumw - sumwx/sumw*sumwy/sumw;
   }
   if ((axis1 == 1 && axis2 == 3) || (axis1 == 3 && axis2 == 1)) {
      return sumwxz/sumw - sumwx/sumw*sumwz/sumw;
   }
   if ((axis1 == 2 && axis2 == 3) || (axis1 == 3 && axis2 == 2)) {
      return sumwyz/sumw - sumwy/sumw*sumwz/sumw;
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Return 3 random numbers along axis x , y and z distributed according
/// to the cell-contents of this 3-dim histogram
/// @param[out] x  reference to random generated x value
/// @param[out] y  reference to random generated y value
/// @param[out] z  reference to random generated z value
/// @param[in] rng (optional) Random number generator pointer used (default is gRandom)

void TH3::GetRandom3(Double_t &x, Double_t &y, Double_t &z, TRandom * rng)
{
   Int_t nbinsx = GetNbinsX();
   Int_t nbinsy = GetNbinsY();
   Int_t nbinsz = GetNbinsZ();
   Int_t nxy    = nbinsx*nbinsy;
   Int_t nbins  = nxy*nbinsz;
   Double_t integral;
   // compute integral checking that all bins have positive content (see ROOT-5894)
   if (fIntegral) {
      if (fIntegral[nbins+1] != fEntries) integral = ComputeIntegral(true);
      else integral = fIntegral[nbins];
   } else {
      integral = ComputeIntegral(true);
   }
   if (integral == 0 ) { x = 0; y = 0; z = 0; return;}
   // case histogram has negative bins
   if (integral == TMath::QuietNaN() ) { x = TMath::QuietNaN(); y = TMath::QuietNaN(); z = TMath::QuietNaN(); return;}

   if (!rng) rng = gRandom;
   Double_t r1 = rng->Rndm();
   Int_t ibin = TMath::BinarySearch(nbins,fIntegral,(Double_t) r1);
   Int_t binz = ibin/nxy;
   Int_t biny = (ibin - nxy*binz)/nbinsx;
   Int_t binx = ibin - nbinsx*(biny + nbinsy*binz);
   x = fXaxis.GetBinLowEdge(binx+1);
   if (r1 > fIntegral[ibin]) x +=
      fXaxis.GetBinWidth(binx+1)*(r1-fIntegral[ibin])/(fIntegral[ibin+1] - fIntegral[ibin]);
   y = fYaxis.GetBinLowEdge(biny+1) + fYaxis.GetBinWidth(biny+1)*rng->Rndm();
   z = fZaxis.GetBinLowEdge(binz+1) + fZaxis.GetBinWidth(binz+1)*rng->Rndm();
}


////////////////////////////////////////////////////////////////////////////////
/// Fill the array stats from the contents of this histogram
/// The array stats must be correctly dimensioned in the calling program.
/// stats[0] = sumw
/// stats[1] = sumw2
/// stats[2] = sumwx
/// stats[3] = sumwx2
/// stats[4] = sumwy
/// stats[5] = sumwy2
/// stats[6] = sumwxy
/// stats[7] = sumwz
/// stats[8] = sumwz2
/// stats[9] = sumwxz
/// stats[10]= sumwyz

void TH3::GetStats(Double_t *stats) const
{
   if (fBuffer) ((TH3*)this)->BufferEmpty();

   Int_t bin, binx, biny, binz;
   Double_t w,err;
   Double_t x,y,z;
   if ((fTsumw == 0 && fEntries > 0) || fXaxis.TestBit(TAxis::kAxisRange) || fYaxis.TestBit(TAxis::kAxisRange) || fZaxis.TestBit(TAxis::kAxisRange)) {
      for (bin=0;bin<11;bin++) stats[bin] = 0;

      Int_t firstBinX = fXaxis.GetFirst();
      Int_t lastBinX  = fXaxis.GetLast();
      Int_t firstBinY = fYaxis.GetFirst();
      Int_t lastBinY  = fYaxis.GetLast();
      Int_t firstBinZ = fZaxis.GetFirst();
      Int_t lastBinZ  = fZaxis.GetLast();
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
         if ( !fZaxis.TestBit(TAxis::kAxisRange) ) {
            if (firstBinZ == 1) firstBinZ = 0;
            if (lastBinZ ==  fZaxis.GetNbins() ) lastBinZ += 1;
         }
      }

      // check for labels axis . In that case corresponsing statistics do not make sense and it is set to zero
      Bool_t labelXaxis =  ((const_cast<TAxis&>(fXaxis)).GetLabels() && fXaxis.CanExtend() );
      Bool_t labelYaxis =  ((const_cast<TAxis&>(fYaxis)).GetLabels() && fYaxis.CanExtend() );
      Bool_t labelZaxis =  ((const_cast<TAxis&>(fZaxis)).GetLabels() && fZaxis.CanExtend() );

      for (binz = firstBinZ; binz <= lastBinZ; binz++) {
         z = (!labelZaxis) ? fZaxis.GetBinCenter(binz) : 0;
         for (biny = firstBinY; biny <= lastBinY; biny++) {
            y = (!labelYaxis) ? fYaxis.GetBinCenter(biny) : 0;
            for (binx = firstBinX; binx <= lastBinX; binx++) {
               bin = GetBin(binx,biny,binz);
               x = (!labelXaxis) ? fXaxis.GetBinCenter(binx) : 0;
               //w   = TMath::Abs(GetBinContent(bin));
               w   = RetrieveBinContent(bin);
               err = TMath::Abs(GetBinError(bin));
               stats[0] += w;
               stats[1] += err*err;
               stats[2] += w*x;
               stats[3] += w*x*x;
               stats[4] += w*y;
               stats[5] += w*y*y;
               stats[6] += w*x*y;
               stats[7] += w*z;
               stats[8] += w*z*z;
               stats[9] += w*x*z;
               stats[10]+= w*y*z;
            }
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
      stats[7] = fTsumwz;
      stats[8] = fTsumwz2;
      stats[9] = fTsumwxz;
      stats[10]= fTsumwyz;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Return integral of bin contents. Only bins in the bins range are considered.
/// By default the integral is computed as the sum of bin contents in the range.
/// if option "width" is specified, the integral is the sum of
/// the bin contents multiplied by the bin width in x, y and in z.

Double_t TH3::Integral(Option_t *option) const
{
   return Integral(fXaxis.GetFirst(),fXaxis.GetLast(),
      fYaxis.GetFirst(),fYaxis.GetLast(),
      fZaxis.GetFirst(),fZaxis.GetLast(),option);
}


////////////////////////////////////////////////////////////////////////////////
/// Return integral of bin contents in range [binx1,binx2],[biny1,biny2],[binz1,binz2]
/// for a 3-D histogram
/// By default the integral is computed as the sum of bin contents in the range.
/// if option "width" is specified, the integral is the sum of
/// the bin contents multiplied by the bin width in x, y and in z.

Double_t TH3::Integral(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2,
                       Int_t binz1, Int_t binz2, Option_t *option) const
{
   Double_t err = 0;
   return DoIntegral(binx1,binx2,biny1,biny2,binz1,binz2,err,option);
}


////////////////////////////////////////////////////////////////////////////////
/// Return integral of bin contents in range [binx1,binx2],[biny1,biny2],[binz1,binz2]
/// for a 3-D histogram. Calculates also the integral error using error propagation
/// from the bin errors assuming that all the bins are uncorrelated.
/// By default the integral is computed as the sum of bin contents in the range.
/// if option "width" is specified, the integral is the sum of
/// the bin contents multiplied by the bin width in x, y and in z.

Double_t TH3::IntegralAndError(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2,
                               Int_t binz1, Int_t binz2,
                               Double_t & error, Option_t *option) const
{
   return DoIntegral(binx1,binx2,biny1,biny2,binz1,binz2,error,option,kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
///Not yet implemented

Double_t TH3::Interpolate(Double_t) const
{
   Error("Interpolate","This function must be called with 3 arguments for a TH3");
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
///Not yet implemented

Double_t TH3::Interpolate(Double_t, Double_t) const
{
   Error("Interpolate","This function must be called with 3 arguments for a TH3");
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Given a point P(x,y,z), Interpolate approximates the value via trilinear interpolation
/// based on the 8 nearest bin center points (corner of the cube surrounding the points)
/// The Algorithm is described in http://en.wikipedia.org/wiki/Trilinear_interpolation
/// The given values (x,y,z) must be between first bin center and  last bin center for each coordinate:
///
///   fXAxis.GetBinCenter(1) < x  < fXaxis.GetBinCenter(nbinX)     AND
///   fYAxis.GetBinCenter(1) < y  < fYaxis.GetBinCenter(nbinY)     AND
///   fZAxis.GetBinCenter(1) < z  < fZaxis.GetBinCenter(nbinZ)

Double_t TH3::Interpolate(Double_t x, Double_t y, Double_t z) const
{
   Int_t ubx = fXaxis.FindFixBin(x);
   if ( x < fXaxis.GetBinCenter(ubx) ) ubx -= 1;
   Int_t obx = ubx + 1;

   Int_t uby = fYaxis.FindFixBin(y);
   if ( y < fYaxis.GetBinCenter(uby) ) uby -= 1;
   Int_t oby = uby + 1;

   Int_t ubz = fZaxis.FindFixBin(z);
   if ( z < fZaxis.GetBinCenter(ubz) ) ubz -= 1;
   Int_t obz = ubz + 1;


//    if ( IsBinUnderflow(GetBin(ubx, uby, ubz)) ||
//         IsBinOverflow (GetBin(obx, oby, obz)) ) {
   if (ubx <=0 || uby <=0 || ubz <= 0 ||
       obx > fXaxis.GetNbins() || oby > fYaxis.GetNbins() || obz > fZaxis.GetNbins() ) {
      Error("Interpolate","Cannot interpolate outside histogram domain.");
      return 0;
   }

   Double_t xw = fXaxis.GetBinCenter(obx) - fXaxis.GetBinCenter(ubx);
   Double_t yw = fYaxis.GetBinCenter(oby) - fYaxis.GetBinCenter(uby);
   Double_t zw = fZaxis.GetBinCenter(obz) - fZaxis.GetBinCenter(ubz);

   Double_t xd = (x - fXaxis.GetBinCenter(ubx)) / xw;
   Double_t yd = (y - fYaxis.GetBinCenter(uby)) / yw;
   Double_t zd = (z - fZaxis.GetBinCenter(ubz)) / zw;


   Double_t v[] = { GetBinContent( ubx, uby, ubz ), GetBinContent( ubx, uby, obz ),
                    GetBinContent( ubx, oby, ubz ), GetBinContent( ubx, oby, obz ),
                    GetBinContent( obx, uby, ubz ), GetBinContent( obx, uby, obz ),
                    GetBinContent( obx, oby, ubz ), GetBinContent( obx, oby, obz ) };


   Double_t i1 = v[0] * (1 - zd) + v[1] * zd;
   Double_t i2 = v[2] * (1 - zd) + v[3] * zd;
   Double_t j1 = v[4] * (1 - zd) + v[5] * zd;
   Double_t j2 = v[6] * (1 - zd) + v[7] * zd;


   Double_t w1 = i1 * (1 - yd) + i2 * yd;
   Double_t w2 = j1 * (1 - yd) + j2 * yd;


   Double_t result = w1 * (1 - xd) + w2 * xd;

   return result;
}


////////////////////////////////////////////////////////////////////////////////
///  Statistical test of compatibility in shape between
///  THIS histogram and h2, using Kolmogorov test.
///     Default: Ignore under- and overflow bins in comparison
///
///     option is a character string to specify options
///         "U" include Underflows in test
///         "O" include Overflows
///         "N" include comparison of normalizations
///         "D" Put out a line of "Debug" printout
///         "M" Return the Maximum Kolmogorov distance instead of prob
///
///   The returned function value is the probability of test
///       (much less than one means NOT compatible)
///
///   The KS test uses the distance between the pseudo-CDF's obtained
///   from the histogram. Since in more than 1D the order for generating the pseudo-CDF is
///   arbitrary, we use the pseudo-CDF's obtained from all the possible 6 combinations of the 3 axis.
///   The average of all the maximum  distances obtained is used in the tests.

Double_t TH3::KolmogorovTest(const TH1 *h2, Option_t *option) const
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
   const TAxis *zaxis1 = h1->GetZaxis();
   const TAxis *zaxis2 = h2->GetZaxis();
   Int_t ncx1   = xaxis1->GetNbins();
   Int_t ncx2   = xaxis2->GetNbins();
   Int_t ncy1   = yaxis1->GetNbins();
   Int_t ncy2   = yaxis2->GetNbins();
   Int_t ncz1   = zaxis1->GetNbins();
   Int_t ncz2   = zaxis2->GetNbins();

   // Check consistency of dimensions
   if (h1->GetDimension() != 3 || h2->GetDimension() != 3) {
      Error("KolmogorovTest","Histograms must be 3-D\n");
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
   if (ncz1 != ncz2) {
      Error("KolmogorovTest","Number of channels in Z is different, %d and %d\n",ncz1,ncz2);
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
   diff1 = TMath::Abs(zaxis1->GetXmin() - zaxis2->GetXmin());
   diff2 = TMath::Abs(zaxis1->GetXmax() - zaxis2->GetXmax());
   if (diff1 > difprec || diff2 > difprec) {
      Error("KolmogorovTest","histograms with different binning along Z");
      return 0;
   }

   //   Should we include Uflows, Oflows?
   Int_t ibeg = 1, jbeg = 1, kbeg = 1;
   Int_t iend = ncx1, jend = ncy1, kend = ncz1;
   if (opt.Contains("U")) {ibeg = 0; jbeg = 0; kbeg = 0;}
   if (opt.Contains("O")) {iend = ncx1+1; jend = ncy1+1; kend = ncz1+1;}

   Int_t i,j,k,bin;
   Double_t sum1  = 0;
   Double_t sum2  = 0;
   Double_t w1    = 0;
   Double_t w2    = 0;
   for (i = ibeg; i <= iend; i++) {
      for (j = jbeg; j <= jend; j++) {
         for (k = kbeg; k <= kend; k++) {
            bin = h1->GetBin(i,j,k);
            sum1 += h1->GetBinContent(bin);
            sum2 += h2->GetBinContent(bin);
            Double_t ew1   = h1->GetBinError(bin);
            Double_t ew2   = h2->GetBinError(bin);
            w1   += ew1*ew1;
            w2   += ew2*ew2;
         }
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

   //   Find Kolmogorov distance
   //   order is arbitrary take average of all possible 6 starting orders x,y,z
   int order[3] = {0,1,2};
   int binbeg[3];
   int binend[3];
   int ibin[3];
   binbeg[0] = ibeg; binbeg[1] = jbeg; binbeg[2] = kbeg;
   binend[0] = iend; binend[1] = jend; binend[2] = kend;
   Double_t vdfmax[6]; // there are in total 6 combinations
   int icomb = 0;
   Double_t s1 = 1./(6.*sum1);
   Double_t s2 = 1./(6.*sum2);
   Double_t rsum1=0, rsum2=0;
   do {
      // loop on bins
      Double_t dmax = 0;
      for (i = binbeg[order[0] ]; i <= binend[order[0] ]; i++) {
         for ( j = binbeg[order[1] ]; j <= binend[order[1] ]; j++) {
            for ( k = binbeg[order[2] ]; k <= binend[order[2] ]; k++) {
                  ibin[ order[0] ] = i;
                  ibin[ order[1] ] = j;
                  ibin[ order[2] ] = k;
                  bin = h1->GetBin(ibin[0],ibin[1],ibin[2]);
                  rsum1 += s1*h1->GetBinContent(bin);
                  rsum2 += s2*h2->GetBinContent(bin);
                  dmax   = TMath::Max(dmax, TMath::Abs(rsum1-rsum2));
            }
         }
      }
      vdfmax[icomb] = dmax;
      icomb++;
   } while (TMath::Permute(3,order)  );


   // get average of distances
   Double_t dfmax = TMath::Mean(6,vdfmax);

   //    Get Kolmogorov probability
   Double_t factnm;
   if (afunc1)      factnm = TMath::Sqrt(sum2);
   else if (afunc2) factnm = TMath::Sqrt(sum1);
   else             factnm = TMath::Sqrt(sum1*sum2/(sum1+sum2));
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

   if (opt.Contains("M"))      return dfmax;  // return average of max distance

   return prb;
}


////////////////////////////////////////////////////////////////////////////////
/// Project a 3-D histogram into a 1-D histogram along X.
///
///   The projection is always of the type TH1D.
///   The projection is made from the cells along the X axis
///   ranging from iymin to iymax and izmin to izmax included.
///   By default, underflow and overflows are included in both the Y and Z axis.
///   By Setting iymin=1 and iymax=NbinsY the underflow and/or overflow in Y will be excluded
///   By setting izmin=1 and izmax=NbinsZ the underflow and/or overflow in Z will be excluded
///
///   if option "e" is specified, the errors are computed.
///   if option "d" is specified, the projection is drawn in the current pad.
///   if option "o" original axis range of the target axes will be
///   kept, but only bins inside the selected range will be filled.
///
///   NOTE that if a TH1D named "name" exists in the current directory or pad
///   the histogram is reset and filled again with the projected contents of the TH3.
///
///  implemented using Project3D

TH1D *TH3::ProjectionX(const char *name, Int_t iymin, Int_t iymax,
                       Int_t izmin, Int_t izmax, Option_t *option) const
{
   // in case of default name append the parent name
   TString hname = name;
   if (hname == "_px") hname = TString::Format("%s%s", GetName(), name);
   TString title =  TString::Format("%s ( Projection X )",GetTitle());

   // when projecting in Z outer axis are Y and Z (order is important. It is defined in the DoProject1D function)
   return DoProject1D(hname, title, iymin, iymax, izmin, izmax, &fXaxis, &fYaxis, &fZaxis, option);
}


////////////////////////////////////////////////////////////////////////////////
/// Project a 3-D histogram into a 1-D histogram along Y.
///
///   The projection is always of the type TH1D.
///   The projection is made from the cells along the Y axis
///   ranging from ixmin to ixmax and izmin to izmax included.
///   By default, underflow and overflow are included in both the X and Z axis.
///   By setting ixmin=1 and ixmax=NbinsX the underflow and/or overflow in X will be excluded
///   By setting izmin=1 and izmax=NbinsZ the underflow and/or overflow in Z will be excluded
///
///   if option "e" is specified, the errors are computed.
///   if option "d" is specified, the projection is drawn in the current pad.
///   if option "o" original axis range of the target axes will be
///   kept, but only bins inside the selected range will be filled.
///
///   NOTE that if a TH1D named "name" exists in the current directory or pad,
///   the histogram is reset and filled again with the projected contents of the TH3.
///
///  implemented using Project3D

TH1D *TH3::ProjectionY(const char *name, Int_t ixmin, Int_t ixmax,
                       Int_t izmin, Int_t izmax, Option_t *option) const
{
   TString hname = name;
   if (hname == "_py") hname = TString::Format("%s%s", GetName(), name);
   TString title =  TString::Format("%s ( Projection Y )",GetTitle());

   // when projecting in Z outer axis are X and Y (order is important. It is defined in the DoProject1D function)
   return DoProject1D(hname, title, ixmin, ixmax, izmin, izmax, &fYaxis, &fXaxis, &fZaxis, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Project a 3-D histogram into a 1-D histogram along Z.
///
///   The projection is always of the type TH1D.
///   The projection is made from the cells along the Z axis
///   ranging from ixmin to ixmax and iymin to iymax included.
///   By default, bins 1 to nx and 1 to ny  are included
///   By default, underflow and overflow are included in both the X and Y axis.
///   By Setting ixmin=1 and ixmax=NbinsX the underflow and/or overflow in X will be excluded
///   By setting iymin=1 and/or iymax=NbinsY the underflow and/or overflow in Y will be excluded
///
///   if option "e" is specified, the errors are computed.
///   if option "d" is specified, the projection is drawn in the current pad.
///   if option "o" original axis range of the target axes will be
///   kept, but only bins inside the selected range will be filled.
///
///   NOTE that if a TH1D named "name" exists in the current directory or pad,
///   the histogram is reset and filled again with the projected contents of the TH3.
///
///  implemented using Project3D

TH1D *TH3::ProjectionZ(const char *name, Int_t ixmin, Int_t ixmax,
                       Int_t iymin, Int_t iymax, Option_t *option) const
{

   TString hname = name;
   if (hname == "_pz") hname = TString::Format("%s%s", GetName(), name);
   TString title =  TString::Format("%s ( Projection Z )",GetTitle());

   // when projecting in Z outer axis are X and Y (order is important. It is defined in the DoProject1D function)
   return DoProject1D(hname, title, ixmin, ixmax, iymin, iymax, &fZaxis, &fXaxis, &fYaxis, option);
}


////////////////////////////////////////////////////////////////////////////////
/// internal method performing the projection to 1D histogram
/// called from TH3::Project3D

TH1D *TH3::DoProject1D(const char* name, const char * title, int imin1, int imax1, int imin2, int imax2,
                       const TAxis* projAxis, const TAxis * axis1, const TAxis * axis2, Option_t * option) const
{

   TString opt = option;
   opt.ToLower();

   // save previous axis range and bits
   // Int_t iminOld1 = axis1->GetFirst();
   // Int_t imaxOld1 = axis1->GetLast();
   // Int_t iminOld2 = axis2->GetFirst();
   // Int_t imaxOld2 = axis2->GetLast();
   // Bool_t hadRange1 = axis1->TestBit(TAxis::kAxisRange);
   // Bool_t hadRange2 = axis2->TestBit(TAxis::kAxisRange);

   // need to cast-away constness to set range
   TAxis out1(*axis1);
   TAxis out2(*axis2);
   // const_cast<TAxis *>(axis1)->SetRange(imin1, imax1);
   // const_cast<TAxis*>(axis2)->SetRange(imin2,imax2);
   out1.SetRange(imin1, imax1);
   out2.SetRange(imin2, imax2);

   Bool_t computeErrors = GetSumw2N();
   if (opt.Contains("e") ) {
      computeErrors = kTRUE;
      opt.Remove(opt.First("e"),1);
   }
   Bool_t originalRange = kFALSE;
   if (opt.Contains('o') ) {
      originalRange = kTRUE;
      opt.Remove(opt.First("o"),1);
   }

   TH1D * h1 = DoProject1D(name, title, projAxis, &out1, &out2, computeErrors, originalRange,true,true);

   // // restore original range
   // if (axis1->TestBit(TAxis::kAxisRange)) {
   //    if (hadRange1) const_cast<TAxis*>(axis1)->SetRange(iminOld1,imaxOld1);
   // if (axis2->TestBit(TAxis::kAxisRange)) const_cast<TAxis*>(axis2)->SetRange(iminOld2,imaxOld2);
   // // we need also to restore the original bits

   // draw in current pad
   if (h1 && opt.Contains("d")) {
      opt.Remove(opt.First("d"),1);
      TVirtualPad *padsav = gPad;
      TVirtualPad *pad = gROOT->GetSelectedPad();
      if (pad) pad->cd();
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
/// internal methdod performing the projection to 1D histogram
/// called from other TH3::DoProject1D

TH1D *TH3::DoProject1D(const char* name, const char * title, const TAxis* projX,
                       const TAxis * out1, const TAxis * out2,
                       bool computeErrors, bool originalRange,
                       bool useUF, bool useOF) const
{
   // Create the projection histogram
   TH1D *h1 = 0;

   // Get range to use as well as bin limits
   // Projected range must be inside and not outside original one (ROOT-8781)
   Int_t ixmin = std::max(projX->GetFirst(),1);
   Int_t ixmax = std::min(projX->GetLast(),projX->GetNbins());
   Int_t nx = ixmax-ixmin+1;

   // Create the histogram, either reseting a preexisting one
   TObject *h1obj = gROOT->FindObject(name);
   if (h1obj && h1obj->InheritsFrom(TH1::Class())) {
      if (h1obj->IsA() != TH1D::Class() ) {
         Error("DoProject1D","Histogram with name %s must be a TH1D and is a %s",name,h1obj->ClassName());
         return 0;
      }
      h1 = (TH1D*)h1obj;
      // reset histogram and re-set the axis in any case
      h1->Reset();
      const TArrayD *bins = projX->GetXbins();
      if ( originalRange )
      {
         if (bins->fN == 0) {
            h1->SetBins(projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
         } else {
            h1->SetBins(projX->GetNbins(),bins->fArray);
         }
      } else {
         if (bins->fN == 0) {
            h1->SetBins(nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
         } else {
            h1->SetBins(nx,&bins->fArray[ixmin-1]);
         }
      }
   }

   if (!h1) {
      const TArrayD *bins = projX->GetXbins();
      if ( originalRange )
      {
         if (bins->fN == 0) {
            h1 = new TH1D(name,title,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
         } else {
            h1 = new TH1D(name,title,projX->GetNbins(),bins->fArray);
         }
      } else {
         if (bins->fN == 0) {
            h1 = new TH1D(name,title,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
         } else {
            h1 = new TH1D(name,title,nx,&bins->fArray[ixmin-1]);
         }
      }
   }

   // Copy the axis attributes and the axis labels if needed.
   h1->GetXaxis()->ImportAttributes(projX);
   THashList* labels = projX->GetLabels();
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

   // Activate errors
   if ( computeErrors && (h1->GetSumw2N() != h1->GetNcells() ) ) h1->Sumw2();

   // Set references to the axies in case out1 or out2 ar enot provided
   // and one can use the histogram axis given projX
   if (out1 == nullptr && out2 == nullptr) {
       if (projX == GetXaxis()) {
         out1 = GetYaxis();
         out2 = GetZaxis();
      } else if (projX == GetYaxis()) {
         out1 = GetXaxis();
         out2 = GetZaxis();
      } else {
         out1 = GetXaxis();
         out2 = GetYaxis();
      }
   }
   R__ASSERT(out1 != nullptr && out2 != nullptr);

   Int_t *refX = 0, *refY = 0, *refZ = 0;
   Int_t ixbin, out1bin, out2bin;
   if (projX == GetXaxis()) {
      refX = &ixbin;
      refY = &out1bin;
      refZ = &out2bin;
   }
   if (projX == GetYaxis()) {
      refX = &out1bin;
      refY = &ixbin;
      refZ = &out2bin;
   }
   if (projX == GetZaxis()) {
      refX = &out1bin;
      refY = &out2bin;
      refZ = &ixbin;
   }
   R__ASSERT (refX != 0 && refY != 0 && refZ != 0);

   // Fill the projected histogram excluding underflow/overflows if considered in the option
   // if specified in the option (by default they considered)
   Double_t totcont  = 0;

   Int_t out1min = out1->GetFirst();
   Int_t out1max = out1->GetLast();
   // GetFirst(), GetLast() can return (0,0) when the range bit is set artificially (see TAxis::SetRange)
 //if (out1min == 0 && out1max == 0) { out1min = 1; out1max = out1->GetNbins(); }
   // correct for underflow/overflows
   if (useUF && !out1->TestBit(TAxis::kAxisRange) )  out1min -= 1;
   if (useOF && !out1->TestBit(TAxis::kAxisRange) )  out1max += 1;
   Int_t out2min = out2->GetFirst();
   Int_t out2max = out2->GetLast();
//   if (out2min == 0 && out2max == 0) { out2min = 1; out2max = out2->GetNbins(); }
   if (useUF && !out2->TestBit(TAxis::kAxisRange) )  out2min -= 1;
   if (useOF && !out2->TestBit(TAxis::kAxisRange) )  out2max += 1;

   for (ixbin=0;ixbin<=1+projX->GetNbins();ixbin++) {
      if ( projX->TestBit(TAxis::kAxisRange) && ( ixbin < ixmin || ixbin > ixmax )) continue;

      Double_t cont = 0;
      Double_t err2 = 0;

      // loop on the bins to be integrated (outbin should be called inbin)
      for (out1bin = out1min; out1bin <= out1max; out1bin++) {
         for (out2bin = out2min; out2bin <= out2max; out2bin++) {

            Int_t bin = GetBin(*refX, *refY, *refZ);

            // sum the bin contents and errors if needed
            cont += RetrieveBinContent(bin);
            if (computeErrors) {
               Double_t exyz = GetBinError(bin);
               err2 += exyz*exyz;
            }
         }
      }
      Int_t ix    = h1->FindBin( projX->GetBinCenter(ixbin) );
      h1->SetBinContent(ix ,cont);
      if (computeErrors) h1->SetBinError(ix, TMath::Sqrt(err2) );
      // sum all content
      totcont += cont;

   }

   // since we use a combination of fill and SetBinError we need to reset and recalculate the statistics
   // for weighted histograms otherwise sumw2 will be wrong.
   // We  can keep the original statistics from the TH3 if the projected sumw is consistent with original one
   // i.e. when no events are thrown away
   bool resetStats = true;
   double eps = 1.E-12;
   if (IsA() == TH3F::Class() ) eps = 1.E-6;
   if (fTsumw != 0 && TMath::Abs( fTsumw - totcont) <  TMath::Abs(fTsumw) * eps) resetStats = false;

   bool resetEntries = resetStats;
   // entries are calculated using underflow/overflow. If excluded entries must be reset
   resetEntries |= !useUF || !useOF;


   if (!resetStats) {
      Double_t stats[kNstat];
      GetStats(stats);
      if ( projX == GetYaxis() ) {
         stats[2] = stats[4];
         stats[3] = stats[5];
      }
      else if  ( projX == GetZaxis() ) {
         stats[2] = stats[7];
         stats[3] = stats[8];
      }
      h1->PutStats(stats);
   }
   else {
      // reset statistics
      h1->ResetStats();
   }
   if (resetEntries) {
      // in case of error calculation (i.e. when Sumw2() is set)
      // use the effective entries for the entries
      // since this  is the only way to estimate them
      Double_t entries =  TMath::Floor( totcont + 0.5); // to avoid numerical rounding
      if (computeErrors) entries = h1->GetEffectiveEntries();
      h1->SetEntries( entries );
   }
   else {
      h1->SetEntries( fEntries );
   }

   return h1;
}


////////////////////////////////////////////////////////////////////////////////
/// internal method performing the projection to a 2D histogram
/// called from TH3::Project3D

TH2D *TH3::DoProject2D(const char* name, const char * title, const TAxis* projX, const TAxis* projY,
                    bool computeErrors, bool originalRange,
                    bool useUF, bool useOF) const
{
   TH2D *h2 = 0;

   // Get range to use as well as bin limits
   Int_t ixmin = std::max(projX->GetFirst(),1);
   Int_t ixmax = std::min(projX->GetLast(),projX->GetNbins());
   Int_t iymin = std::max(projY->GetFirst(),1);
   Int_t iymax = std::min(projY->GetLast(),projY->GetNbins());

   Int_t nx = ixmax-ixmin+1;
   Int_t ny = iymax-iymin+1;

   // Create the histogram, either reseting a preexisting one
   //  or creating one from scratch.
   // Does an object with the same name exists?
   TObject *h2obj = gROOT->FindObject(name);
   if (h2obj && h2obj->InheritsFrom(TH1::Class())) {
      if ( h2obj->IsA() != TH2D::Class() ) {
         Error("DoProject2D","Histogram with name %s must be a TH2D and is a %s",name,h2obj->ClassName());
         return 0;
      }
      h2 = (TH2D*)h2obj;
      // reset histogram and its axes
      h2->Reset();
      const TArrayD *xbins = projX->GetXbins();
      const TArrayD *ybins = projY->GetXbins();
      if ( originalRange ) {
         h2->SetBins(projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                     ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
         // set bins for mixed axis do not exists - need to set afterwards the variable bins
         if (ybins->fN != 0)
            h2->GetXaxis()->Set(projY->GetNbins(),&ybins->fArray[iymin-1]);
         if (xbins->fN != 0)
            h2->GetYaxis()->Set(projX->GetNbins(),&xbins->fArray[ixmin-1]);
      } else {
         h2->SetBins(ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                     ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
         if (ybins->fN != 0)
            h2->GetXaxis()->Set(ny,&ybins->fArray[iymin-1]);
         if (xbins->fN != 0)
            h2->GetYaxis()->Set(nx,&xbins->fArray[ixmin-1]);
      }
   }


   if (!h2) {
      const TArrayD *xbins = projX->GetXbins();
      const TArrayD *ybins = projY->GetXbins();
      if ( originalRange )
      {
         if (xbins->fN == 0 && ybins->fN == 0) {
            h2 = new TH2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                          ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
         } else if (ybins->fN == 0) {
            h2 = new TH2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                          ,projX->GetNbins(),&xbins->fArray[ixmin-1]);
         } else if (xbins->fN == 0) {
            h2 = new TH2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1]
                          ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
         } else {
            h2 = new TH2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1],projX->GetNbins(),&xbins->fArray[ixmin-1]);
         }
      } else {
         if (xbins->fN == 0 && ybins->fN == 0) {
            h2 = new TH2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                          ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
         } else if (ybins->fN == 0) {
            h2 = new TH2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                          ,nx,&xbins->fArray[ixmin-1]);
         } else if (xbins->fN == 0) {
            h2 = new TH2D(name,title,ny,&ybins->fArray[iymin-1]
                          ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
         } else {
            h2 = new TH2D(name,title,ny,&ybins->fArray[iymin-1],nx,&xbins->fArray[ixmin-1]);
         }
      }
   }

   // Copy the axis attributes and the axis labels if needed.
   THashList* labels1 = 0;
   THashList* labels2 = 0;
   // "xy"
   h2->GetXaxis()->ImportAttributes(projY);
   h2->GetYaxis()->ImportAttributes(projX);
   labels1 = projY->GetLabels();
   labels2 = projX->GetLabels();
   if (labels1) {
      TIter iL(labels1);
      TObjString* lb;
      Int_t i = 1;
      while ((lb=(TObjString*)iL())) {
         h2->GetXaxis()->SetBinLabel(i,lb->String().Data());
         i++;
      }
   }
   if (labels2) {
      TIter iL(labels2);
      TObjString* lb;
      Int_t i = 1;
      while ((lb=(TObjString*)iL())) {
         h2->GetYaxis()->SetBinLabel(i,lb->String().Data());
         i++;
      }
   }
   h2->SetLineColor(this->GetLineColor());
   h2->SetFillColor(this->GetFillColor());
   h2->SetMarkerColor(this->GetMarkerColor());
   h2->SetMarkerStyle(this->GetMarkerStyle());

   // Activate errors
   if ( computeErrors && (h2->GetSumw2N() != h2->GetNcells()) ) h2->Sumw2();

   // Set references to the axis, so that the bucle has no branches.
   const TAxis* out = 0;
   if ( projX != GetXaxis() && projY != GetXaxis() ) {
      out = GetXaxis();
   } else if ( projX != GetYaxis() && projY != GetYaxis() ) {
      out = GetYaxis();
   } else {
      out = GetZaxis();
   }

   Int_t *refX = 0, *refY = 0, *refZ = 0;
   Int_t ixbin, iybin, outbin;
   if ( projX == GetXaxis() && projY == GetYaxis() ) { refX = &ixbin;  refY = &iybin;  refZ = &outbin; }
   if ( projX == GetYaxis() && projY == GetXaxis() ) { refX = &iybin;  refY = &ixbin;  refZ = &outbin; }
   if ( projX == GetXaxis() && projY == GetZaxis() ) { refX = &ixbin;  refY = &outbin; refZ = &iybin;  }
   if ( projX == GetZaxis() && projY == GetXaxis() ) { refX = &iybin;  refY = &outbin; refZ = &ixbin;  }
   if ( projX == GetYaxis() && projY == GetZaxis() ) { refX = &outbin; refY = &ixbin;  refZ = &iybin;  }
   if ( projX == GetZaxis() && projY == GetYaxis() ) { refX = &outbin; refY = &iybin;  refZ = &ixbin;  }
   R__ASSERT (refX != 0 && refY != 0 && refZ != 0);

   // Fill the projected histogram excluding underflow/overflows if considered in the option
   // if specified in the option (by default they considered)
   Double_t totcont  = 0;

   Int_t outmin = out->GetFirst();
   Int_t outmax = out->GetLast();
   // GetFirst(), GetLast() can return (0,0) when the range bit is set artificially (see TAxis::SetRange)
   if (outmin == 0 && outmax == 0) { outmin = 1; outmax = out->GetNbins(); }
   // correct for underflow/overflows
   if (useUF && !out->TestBit(TAxis::kAxisRange) )  outmin -= 1;
   if (useOF && !out->TestBit(TAxis::kAxisRange) )  outmax += 1;

   for (ixbin=0;ixbin<=1+projX->GetNbins();ixbin++) {
      if ( projX->TestBit(TAxis::kAxisRange) && ( ixbin < ixmin || ixbin > ixmax )) continue;
      Int_t ix = h2->GetYaxis()->FindBin( projX->GetBinCenter(ixbin) );

      for (iybin=0;iybin<=1+projY->GetNbins();iybin++) {
         if ( projY->TestBit(TAxis::kAxisRange) && ( iybin < iymin || iybin > iymax )) continue;
         Int_t iy = h2->GetXaxis()->FindBin( projY->GetBinCenter(iybin) );

         Double_t cont = 0;
         Double_t err2 = 0;

         // loop on the bins to be integrated (outbin should be called inbin)
         for (outbin = outmin; outbin <= outmax; outbin++) {

            Int_t bin = GetBin(*refX,*refY,*refZ);

            // sum the bin contents and errors if needed
            cont += RetrieveBinContent(bin);
            if (computeErrors) {
               Double_t exyz = GetBinError(bin);
               err2 += exyz*exyz;
            }

         }

         // remember axis are inverted
         h2->SetBinContent(iy , ix, cont);
         if (computeErrors) h2->SetBinError(iy, ix, TMath::Sqrt(err2) );
         // sum all content
         totcont += cont;

      }
   }

   // since we use fill we need to reset and recalculate the statistics (see comment in DoProject1D )
   // or keep original statistics if consistent sumw2
   bool resetStats = true;
   double eps = 1.E-12;
   if (IsA() == TH3F::Class() ) eps = 1.E-6;
   if (fTsumw != 0 && TMath::Abs( fTsumw - totcont) <  TMath::Abs(fTsumw) * eps) resetStats = false;

   bool resetEntries = resetStats;
   // entries are calculated using underflow/overflow. If excluded entries must be reset
   resetEntries |= !useUF || !useOF;

   if (!resetStats) {
      Double_t stats[kNstat];
      Double_t oldst[kNstat]; // old statistics
      for (Int_t i = 0; i < kNstat; ++i) { oldst[i] = 0; }
      GetStats(oldst);
      std::copy(oldst,oldst+kNstat,stats);
      // not that projX refer to Y axis and projX refer to the X axis of projected histogram
      // nothing to do for projection in Y vs X
      if ( projY == GetXaxis() && projX == GetZaxis() ) {  // case XZ
         stats[4] = oldst[7];
         stats[5] = oldst[8];
         stats[6] = oldst[9];
      }
      if ( projY == GetYaxis() ) {
         stats[2] = oldst[4];
         stats[3] = oldst[5];
         if ( projX == GetXaxis() )  { // case YX
            stats[4] = oldst[2];
            stats[5] = oldst[3];
         }
         if ( projX == GetZaxis() )  { // case YZ
            stats[4] = oldst[7];
            stats[5] = oldst[8];
            stats[6] = oldst[10];
         }
      }
      else if  ( projY == GetZaxis() ) {
         stats[2] = oldst[7];
         stats[3] = oldst[8];
         if ( projX == GetXaxis() )  { // case ZX
            stats[4] = oldst[2];
            stats[5] = oldst[3];
            stats[6] = oldst[9];
         }
         if ( projX == GetYaxis() )  { // case ZY
            stats[4] = oldst[4];
            stats[5] = oldst[5];
            stats[6] = oldst[10];
         }
      }
      // set the new statistics
      h2->PutStats(stats);
   }
   else {
      // recalculate the statistics
      h2->ResetStats();
   }

   if (resetEntries) {
      // use the effective entries for the entries
      // since this  is the only way to estimate them
      Double_t entries =  h2->GetEffectiveEntries();
      if (!computeErrors) entries = TMath::Floor( entries + 0.5); // to avoid numerical rounding
      h2->SetEntries( entries );
   }
   else {
      h2->SetEntries( fEntries );
   }


   return h2;
}


////////////////////////////////////////////////////////////////////////////////
/// Project a 3-d histogram into 1 or 2-d histograms depending on the
/// option parameter, which may contain a combination of the characters x,y,z,e
///  - option = "x" return the x projection into a TH1D histogram
///  - option = "y" return the y projection into a TH1D histogram
///  - option = "z" return the z projection into a TH1D histogram
///  - option = "xy" return the x versus y projection into a TH2D histogram
///  - option = "yx" return the y versus x projection into a TH2D histogram
///  - option = "xz" return the x versus z projection into a TH2D histogram
///  - option = "zx" return the z versus x projection into a TH2D histogram
///  - option = "yz" return the y versus z projection into a TH2D histogram
///  - option = "zy" return the z versus y projection into a TH2D histogram
///
/// NB: the notation "a vs b" means "a" vertical and "b" horizontal
///
/// option = "o" original axis range of the target axes will be
///   kept, but only bins inside the selected range will be filled.
///
/// If option contains the string "e", errors are computed
///
/// The projection is made for the selected bins only.
/// To select a bin range along an axis, use TAxis::SetRange, eg
///    h3.GetYaxis()->SetRange(23,56);
///
/// NOTE 1: The generated histogram is named th3name + option
/// eg if the TH3* h histogram is named "myhist", then
/// h->Project3D("xy"); produces a TH2D histogram named "myhist_xy"
/// if a histogram of the same type already exists, it is overwritten.
/// The following sequence
///    h->Project3D("xy");
///    h->Project3D("xy2");
///  will generate two TH2D histograms named "myhist_xy" and "myhist_xy2"
///  A different name can be generated by attaching a string to the option
///  For example h->Project3D("name_xy") will generate an histogram with the name:  h3dname_name_xy.
///
///  NOTE 2: If an histogram of the same type already exists,
///  the histogram is reset and filled again with the projected contents of the TH3.
///
///  NOTE 3: The number of entries in the projected histogram is estimated from the number of
///  effective entries for all the cells included in the projection.
///
///  NOTE 4: underflow/overflow are included by default in the projection
///  To exclude underflow and/or overflow (for both axis in case of a projection to a 1D histogram) use option "NUF" and/or "NOF"
///  With SetRange() you can have all bins except underflow/overflow only if you set the axis bit range as
///  following after having called SetRange:  axis->SetRange(1, axis->GetNbins());

TH1 *TH3::Project3D(Option_t *option) const
{
   TString opt = option; opt.ToLower();
   Int_t pcase = 0;
   TString ptype;
   if (opt.Contains("x"))  { pcase = 1; ptype = "x"; }
   if (opt.Contains("y"))  { pcase = 2; ptype = "y"; }
   if (opt.Contains("z"))  { pcase = 3; ptype = "z"; }
   if (opt.Contains("xy")) { pcase = 4; ptype = "xy"; }
   if (opt.Contains("yx")) { pcase = 5; ptype = "yx"; }
   if (opt.Contains("xz")) { pcase = 6; ptype = "xz"; }
   if (opt.Contains("zx")) { pcase = 7; ptype = "zx"; }
   if (opt.Contains("yz")) { pcase = 8; ptype = "yz"; }
   if (opt.Contains("zy")) { pcase = 9; ptype = "zy"; }

   if (pcase == 0) {
      Error("Project3D","No projection axis specified - return a NULL pointer");
      return 0;
   }
   // do not remove ptype from opt to use later in the projected histo name

   Bool_t computeErrors = GetSumw2N();
   if (opt.Contains("e") ) {
      computeErrors = kTRUE;
      opt.Remove(opt.First("e"),1);
   }

   Bool_t useUF = kTRUE;
   Bool_t useOF = kTRUE;
   if (opt.Contains("nuf") ) {
      useUF = kFALSE;
      opt.Remove(opt.Index("nuf"),3);
   }
   if (opt.Contains("nof") ) {
      useOF = kFALSE;
      opt.Remove(opt.Index("nof"),3);
   }

   Bool_t originalRange = kFALSE;
   if (opt.Contains('o') ) {
      originalRange = kTRUE;
      opt.Remove(opt.First("o"),1);
   }


   // Create the projection histogram
   TH1 *h = 0;

   TString name = GetName();
   TString title = GetTitle();
   name  += "_";   name  += opt;  // opt may include a user defined name
   title += " ";   title += ptype; title += " projection";

   switch (pcase) {
      case 1:
         // "x"
         h = DoProject1D(name, title, this->GetXaxis(), nullptr, nullptr,
                        computeErrors, originalRange, useUF, useOF);
         break;

      case 2:
         // "y"
         h = DoProject1D(name, title, this->GetYaxis(), nullptr, nullptr,
                         computeErrors, originalRange, useUF, useOF);
         break;

      case 3:
         // "z"
         h = DoProject1D(name, title, this->GetZaxis(), nullptr, nullptr,
                         computeErrors, originalRange, useUF, useOF);
         break;

      case 4:
         // "xy"
         h = DoProject2D(name, title, this->GetXaxis(),this->GetYaxis(),
                       computeErrors, originalRange, useUF, useOF);
         break;

      case 5:
         // "yx"
         h = DoProject2D(name, title, this->GetYaxis(),this->GetXaxis(),
                       computeErrors, originalRange, useUF, useOF);
         break;

      case 6:
         // "xz"
         h = DoProject2D(name, title, this->GetXaxis(),this->GetZaxis(),
                       computeErrors, originalRange, useUF, useOF);
         break;

      case 7:
         // "zx"
         h = DoProject2D(name, title, this->GetZaxis(),this->GetXaxis(),
                       computeErrors, originalRange, useUF, useOF);
         break;

      case 8:
         // "yz"
         h = DoProject2D(name, title, this->GetYaxis(),this->GetZaxis(),
                       computeErrors, originalRange, useUF, useOF);
         break;

      case 9:
         // "zy"
         h = DoProject2D(name, title, this->GetZaxis(),this->GetYaxis(),
                       computeErrors, originalRange, useUF, useOF);
         break;

   }

   // draw in current pad
   if (h && opt.Contains("d")) {
      opt.Remove(opt.First("d"),1);
      TVirtualPad *padsav = gPad;
      TVirtualPad *pad = gROOT->GetSelectedPad();
      if (pad) pad->cd();
      if (!gPad || !gPad->FindObject(h)) {
         h->Draw(opt);
      } else {
         h->Paint(opt);
      }
      if (padsav) padsav->cd();
   }

   return h;
}


////////////////////////////////////////////////////////////////////////////////
/// internal function to fill the bins of the projected profile 2D histogram
/// called from DoProjectProfile2D

void TH3::DoFillProfileProjection(TProfile2D * p2,
                                  const TAxis & a1, const TAxis & a2, const TAxis & a3,
                                  Int_t bin1, Int_t bin2, Int_t bin3,
                                  Int_t inBin, Bool_t useWeights ) const {
   Double_t cont = GetBinContent(inBin);
   if (!cont) return;
   TArrayD & binSumw2 = *(p2->GetBinSumw2());
   if (useWeights && binSumw2.fN <= 0) useWeights = false;
   if (!useWeights) p2->SetBit(TH1::kIsNotW);  // to use Fill for setting the bin contents of the Profile
   // the following fill update wrongly the fBinSumw2- need to save it before
   Double_t u = a1.GetBinCenter(bin1);
   Double_t v = a2.GetBinCenter(bin2);
   Double_t w = a3.GetBinCenter(bin3);
   Int_t outBin = p2->FindBin(u, v);
   if (outBin <0) return;
   Double_t tmp = 0;
   if ( useWeights ) tmp = binSumw2.fArray[outBin];
   p2->Fill( u , v, w, cont);
   if (useWeights ) binSumw2.fArray[outBin] = tmp + fSumw2.fArray[inBin];
}


////////////////////////////////////////////////////////////////////////////////
/// internal method to project to a 2D Profile
/// called from TH3::Project3DProfile

TProfile2D *TH3::DoProjectProfile2D(const char* name, const char * title, const TAxis* projX, const TAxis* projY,
                                          bool originalRange, bool useUF, bool useOF) const
{
   // Get the ranges where we will work.
   Int_t ixmin = std::max(projX->GetFirst(),1);
   Int_t ixmax = std::min(projX->GetLast(),projX->GetNbins());
   Int_t iymin = std::max(projY->GetFirst(),1);
   Int_t iymax = std::min(projY->GetLast(),projY->GetNbins());

   Int_t nx = ixmax-ixmin+1;
   Int_t ny = iymax-iymin+1;

   // Create the projected profiles
   TProfile2D *p2 = 0;

   // Create the histogram, either reseting a preexisting one
   // Does an object with the same name exists?
   TObject *p2obj = gROOT->FindObject(name);
   if (p2obj && p2obj->InheritsFrom(TH1::Class())) {
      if (p2obj->IsA() != TProfile2D::Class() ) {
         Error("DoProjectProfile2D","Histogram with name %s must be a TProfile2D and is a %s",name,p2obj->ClassName());
         return 0;
      }
      p2 = (TProfile2D*)p2obj;
      // reset existing profile and re-set bins
      p2->Reset();
      const TArrayD *xbins = projX->GetXbins();
      const TArrayD *ybins = projY->GetXbins();
      if ( originalRange ) {
         p2->SetBins(projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                     ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
         // set bins for mixed axis do not exists - need to set afterwards the variable bins
         if (ybins->fN != 0)
            p2->GetXaxis()->Set(projY->GetNbins(),&ybins->fArray[iymin-1]);
         if (xbins->fN != 0)
            p2->GetYaxis()->Set(projX->GetNbins(),&xbins->fArray[ixmin-1]);
      } else {
         p2->SetBins(ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                     ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
         if (ybins->fN != 0)
            p2->GetXaxis()->Set(ny,&ybins->fArray[iymin-1]);
         if (xbins->fN != 0)
            p2->GetYaxis()->Set(nx,&xbins->fArray[ixmin-1]);
      }
   }

   if (!p2) {
      const TArrayD *xbins = projX->GetXbins();
      const TArrayD *ybins = projY->GetXbins();
      if ( originalRange ) {
         if (xbins->fN == 0 && ybins->fN == 0) {
            p2 = new TProfile2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                                ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
         } else if (ybins->fN == 0) {
            p2 = new TProfile2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                                ,projX->GetNbins(),&xbins->fArray[ixmin-1]);
         } else if (xbins->fN == 0) {
            p2 = new TProfile2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1]
                                ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
         } else {
            p2 = new TProfile2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1],projX->GetNbins(),&xbins->fArray[ixmin-1]);
         }
      } else {
         if (xbins->fN == 0 && ybins->fN == 0) {
            p2 = new TProfile2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                                ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
         } else if (ybins->fN == 0) {
            p2 = new TProfile2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                                ,nx,&xbins->fArray[ixmin-1]);
         } else if (xbins->fN == 0) {
            p2 = new TProfile2D(name,title,ny,&ybins->fArray[iymin-1]
                                ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
         } else {
            p2 = new TProfile2D(name,title,ny,&ybins->fArray[iymin-1],nx,&xbins->fArray[ixmin-1]);
         }
      }
   }

   // Set references to the axis, so that the loop has no branches.
   const TAxis* outAxis = 0;
   if ( projX != GetXaxis() && projY != GetXaxis() ) {
      outAxis = GetXaxis();
   } else if ( projX != GetYaxis() && projY != GetYaxis() ) {
      outAxis = GetYaxis();
   } else {
      outAxis = GetZaxis();
   }

   // Weights management
   bool useWeights = (GetSumw2N() > 0);
   // store sum of w2 in profile if histo is weighted
   if (useWeights && (p2->GetBinSumw2()->fN != p2->GetNcells() ) ) p2->Sumw2();

   // Set references to the bins, so that the loop has no branches.
   Int_t *refX = 0, *refY = 0, *refZ = 0;
   Int_t ixbin, iybin, outbin;
   if ( projX == GetXaxis() && projY == GetYaxis() ) { refX = &ixbin;  refY = &iybin;  refZ = &outbin; }
   if ( projX == GetYaxis() && projY == GetXaxis() ) { refX = &iybin;  refY = &ixbin;  refZ = &outbin; }
   if ( projX == GetXaxis() && projY == GetZaxis() ) { refX = &ixbin;  refY = &outbin; refZ = &iybin;  }
   if ( projX == GetZaxis() && projY == GetXaxis() ) { refX = &iybin;  refY = &outbin; refZ = &ixbin;  }
   if ( projX == GetYaxis() && projY == GetZaxis() ) { refX = &outbin; refY = &ixbin;  refZ = &iybin;  }
   if ( projX == GetZaxis() && projY == GetYaxis() ) { refX = &outbin; refY = &iybin;  refZ = &ixbin;  }
   R__ASSERT (refX != 0 && refY != 0 && refZ != 0);

   Int_t outmin = outAxis->GetFirst();
   Int_t outmax = outAxis->GetLast();
   // GetFirst, GetLast can return underflow or overflow bins
   // correct for underflow/overflows
   if (useUF && !outAxis->TestBit(TAxis::kAxisRange) )  outmin -= 1;
   if (useOF && !outAxis->TestBit(TAxis::kAxisRange) )  outmax += 1;

   TArrayD & binSumw2 = *(p2->GetBinSumw2());
   if (useWeights && binSumw2.fN <= 0) useWeights = false;
   if (!useWeights) p2->SetBit(TH1::kIsNotW);

   // Call specific method for the projection
   for (ixbin=0;ixbin<=1+projX->GetNbins();ixbin++) {
      if ( (ixbin < ixmin || ixbin > ixmax) && projX->TestBit(TAxis::kAxisRange)) continue;
      for ( iybin=0;iybin<=1+projY->GetNbins();iybin++) {
         if ( (iybin < iymin || iybin > iymax) && projX->TestBit(TAxis::kAxisRange)) continue;

         // profile output bin
         Int_t poutBin = p2->FindBin(projY->GetBinCenter(iybin), projX->GetBinCenter(ixbin));
         if (poutBin <0) continue;
         // loop on the bins to be integrated (outbin should be called inbin)
         for (outbin = outmin; outbin <= outmax; outbin++) {

            Int_t bin = GetBin(*refX,*refY,*refZ);

            //DoFillProfileProjection(p2, *projY, *projX, *outAxis, iybin, ixbin, outbin, bin, useWeights);

            Double_t cont = RetrieveBinContent(bin);
            if (!cont) continue;

            Double_t tmp = 0;
            // the following fill update wrongly the fBinSumw2- need to save it before
            if ( useWeights ) tmp = binSumw2.fArray[poutBin];
            p2->Fill( projY->GetBinCenter(iybin) , projX->GetBinCenter(ixbin), outAxis->GetBinCenter(outbin), cont);
            if (useWeights ) binSumw2.fArray[poutBin] = tmp + fSumw2.fArray[bin];

         }
      }
   }

   // recompute statistics for the projected profiles
   // forget about preserving old statistics
   bool resetStats = true;
   Double_t stats[kNstat];
   // reset statistics
   if (resetStats)
      for (Int_t i=0;i<kNstat;i++) stats[i] = 0;

   p2->PutStats(stats);
   Double_t entries = fEntries;
   // recalculate the statistics
   if (resetStats) {
      entries =  p2->GetEffectiveEntries();
      if (!useWeights) entries = TMath::Floor( entries + 0.5); // to avoid numerical rounding
      p2->SetEntries( entries );
   }

   p2->SetEntries(entries);

   return p2;
}


////////////////////////////////////////////////////////////////////////////////
/// Project a 3-d histogram into a 2-d profile histograms depending
/// on the option parameter
/// option may contain a combination of the characters x,y,z
/// option = "xy" return the x versus y projection into a TProfile2D histogram
/// option = "yx" return the y versus x projection into a TProfile2D histogram
/// option = "xz" return the x versus z projection into a TProfile2D histogram
/// option = "zx" return the z versus x projection into a TProfile2D histogram
/// option = "yz" return the y versus z projection into a TProfile2D histogram
/// option = "zy" return the z versus y projection into a TProfile2D histogram
/// NB: the notation "a vs b" means "a" vertical and "b" horizontal
///
/// option = "o" original axis range of the target axes will be
///   kept, but only bins inside the selected range will be filled.
///
/// The projection is made for the selected bins only.
/// To select a bin range along an axis, use TAxis::SetRange, eg
///    h3.GetYaxis()->SetRange(23,56);
///
/// NOTE 1: The generated histogram is named th3name + "_p" + option
/// eg if the TH3* h histogram is named "myhist", then
/// h->Project3D("xy"); produces a TProfile2D histogram named "myhist_pxy".
/// The following sequence
///    h->Project3DProfile("xy");
///    h->Project3DProfile("xy2");
///  will generate two TProfile2D histograms named "myhist_pxy" and "myhist_pxy2"
///  So, passing additional characters in the option string one can customize the name.
///
///  NOTE 2: If a profile of the same type already exists with compatible axes,
///  the profile is reset and filled again with the projected contents of the TH3.
///  In the case of axes incompatibility, an error is reported and a NULL pointer is returned.
///
///  NOTE 3: The number of entries in the projected profile is estimated from the number of
///  effective entries for all the cells included in the projection.
///
///  NOTE 4: underflow/overflow are by default excluded from the projection
///  (Note that this is a different default behavior compared to the projection to an histogram)
///  To include the underflow and/or overflow use option "UF" and/or "OF"

TProfile2D *TH3::Project3DProfile(Option_t *option) const
{
   TString opt = option; opt.ToLower();
   Int_t pcase = 0;
   TString ptype;
   if (opt.Contains("xy")) { pcase = 4; ptype = "xy"; }
   if (opt.Contains("yx")) { pcase = 5; ptype = "yx"; }
   if (opt.Contains("xz")) { pcase = 6; ptype = "xz"; }
   if (opt.Contains("zx")) { pcase = 7; ptype = "zx"; }
   if (opt.Contains("yz")) { pcase = 8; ptype = "yz"; }
   if (opt.Contains("zy")) { pcase = 9; ptype = "zy"; }

   if (pcase == 0) {
      Error("Project3D","No projection axis specified - return a NULL pointer");
      return 0;
   }
   // do not remove ptype from opt to use later in the projected histo name

   Bool_t useUF = kFALSE;
   if (opt.Contains("uf") ) {
      useUF = kTRUE;
      opt.Remove(opt.Index("uf"),2);
   }
   Bool_t useOF = kFALSE;
   if (opt.Contains("of") ) {
      useOF = kTRUE;
      opt.Remove(opt.Index("of"),2);
   }

   Bool_t originalRange = kFALSE;
   if (opt.Contains('o') ) {
      originalRange = kTRUE;
      opt.Remove(opt.First("o"),1);
   }

   // Create the projected profile
   TProfile2D *p2 = 0;
   TString name = GetName();
   TString title = GetTitle();
   name  += "_p";   name  += opt;  // opt may include a user defined name
   title += " profile ";   title += ptype; title += " projection";

   // Call the method with the specific projected axes.
   switch (pcase) {
      case 4:
         // "xy"
         p2 = DoProjectProfile2D(name, title, GetXaxis(), GetYaxis(), originalRange, useUF, useOF);
         break;

      case 5:
         // "yx"
         p2 = DoProjectProfile2D(name, title, GetYaxis(), GetXaxis(), originalRange, useUF, useOF);
         break;

      case 6:
         // "xz"
         p2 = DoProjectProfile2D(name, title, GetXaxis(), GetZaxis(), originalRange, useUF, useOF);
         break;

      case 7:
         // "zx"
         p2 = DoProjectProfile2D(name, title, GetZaxis(), GetXaxis(), originalRange, useUF, useOF);
         break;

      case 8:
         // "yz"
         p2 = DoProjectProfile2D(name, title, GetYaxis(), GetZaxis(), originalRange, useUF, useOF);
         break;

      case 9:
         // "zy"
         p2 = DoProjectProfile2D(name, title, GetZaxis(), GetYaxis(), originalRange, useUF, useOF);
         break;

   }

   return p2;
}


////////////////////////////////////////////////////////////////////////////////
/// Replace current statistics with the values in array stats

void TH3::PutStats(Double_t *stats)
{
   TH1::PutStats(stats);
   fTsumwy  = stats[4];
   fTsumwy2 = stats[5];
   fTsumwxy = stats[6];
   fTsumwz  = stats[7];
   fTsumwz2 = stats[8];
   fTsumwxz = stats[9];
   fTsumwyz = stats[10];
}


////////////////////////////////////////////////////////////////////////////////
/// Rebin only the X axis
/// see Rebin3D

TH3 *TH3::RebinX(Int_t ngroup, const char *newname)
{
  return Rebin3D(ngroup, 1, 1, newname);
}


////////////////////////////////////////////////////////////////////////////////
/// Rebin only the Y axis
/// see Rebin3D

TH3 *TH3::RebinY(Int_t ngroup, const char *newname)
{
  return Rebin3D(1, ngroup, 1, newname);
}


////////////////////////////////////////////////////////////////////////////////
/// Rebin only the Z axis
/// see Rebin3D

TH3 *TH3::RebinZ(Int_t ngroup, const char *newname)
{
  return Rebin3D(1, 1, ngroup, newname);

}


////////////////////////////////////////////////////////////////////////////////
/// Rebin this histogram grouping nxgroup/nygroup/nzgroup bins along the xaxis/yaxis/zaxis together.
///
///   if newname is not blank a new temporary histogram hnew is created.
///   else the current histogram is modified (default)
///   The parameter nxgroup/nygroup indicate how many bins along the xaxis/yaxis of this
///   have to me merged into one bin of hnew
///   If the original histogram has errors stored (via Sumw2), the resulting
///   histograms has new errors correctly calculated.
///
///   examples: if hpxpy is an existing TH3 histogram with 40 x 40 x 40 bins
///     hpxpypz->Rebin3D();  // merges two bins along the xaxis and yaxis in one in hpxpypz
///                          // Carefull: previous contents of hpxpy are lost
///     hpxpypz->RebinX(5);  //merges five bins along the xaxis in one in hpxpypz
///     TH3 *hnew = hpxpypz->RebinY(5,"hnew"); // creates a new histogram hnew
///                                          // merging 5 bins of h1 along the yaxis in one bin
///
///   NOTE : If nxgroup/nygroup is not an exact divider of the number of bins,
///          along the xaxis/yaxis the top limit(s) of the rebinned histogram
///          is changed to the upper edge of the xbin=newxbins*nxgroup resp.
///          ybin=newybins*nygroup and the corresponding bins are added to
///          the overflow bin.
///          Statistics will be recomputed from the new bin contents.

TH3 *TH3::Rebin3D(Int_t nxgroup, Int_t nygroup, Int_t nzgroup, const char *newname)
{
   Int_t i,j,k,xbin,ybin,zbin;
   Int_t nxbins  = fXaxis.GetNbins();
   Int_t nybins  = fYaxis.GetNbins();
   Int_t nzbins  = fZaxis.GetNbins();
   Double_t xmin  = fXaxis.GetXmin();
   Double_t xmax  = fXaxis.GetXmax();
   Double_t ymin  = fYaxis.GetXmin();
   Double_t ymax  = fYaxis.GetXmax();
   Double_t zmin  = fZaxis.GetXmin();
   Double_t zmax  = fZaxis.GetXmax();
   if ((nxgroup <= 0) || (nxgroup > nxbins)) {
      Error("Rebin", "Illegal value of nxgroup=%d",nxgroup);
      return 0;
   }
   if ((nygroup <= 0) || (nygroup > nybins)) {
      Error("Rebin", "Illegal value of nygroup=%d",nygroup);
      return 0;
   }
   if ((nzgroup <= 0) || (nzgroup > nzbins)) {
      Error("Rebin", "Illegal value of nzgroup=%d",nzgroup);
      return 0;
   }

   Int_t newxbins = nxbins/nxgroup;
   Int_t newybins = nybins/nygroup;
   Int_t newzbins = nzbins/nzgroup;

   // Save old bin contents into a new array
   Double_t entries = fEntries;
   Double_t *oldBins = new Double_t[fNcells];
   for (Int_t ibin = 0; ibin < fNcells; ibin++) {
      oldBins[ibin] = RetrieveBinContent(ibin);
   }
   Double_t *oldSumw2 = 0;
   if (fSumw2.fN != 0) {
      oldSumw2 = new Double_t[fNcells];
      for (Int_t ibin = 0; ibin < fNcells; ibin++) {
         oldSumw2[ibin] = fSumw2.fArray[ibin];
      }
   }

   // create a clone of the old histogram if newname is specified
   TH3 *hnew = this;
   if (newname && strlen(newname)) {
      hnew = (TH3*)Clone();
      hnew->SetName(newname);
   }

   // save original statistics
   Double_t stat[kNstat];
   GetStats(stat);
   bool resetStat = false;


   // change axis specs and rebuild bin contents array
   if (newxbins*nxgroup != nxbins) {
      xmax = fXaxis.GetBinUpEdge(newxbins*nxgroup);
      resetStat = true; //stats must be reset because top bins will be moved to overflow bin
   }
   if (newybins*nygroup != nybins) {
      ymax = fYaxis.GetBinUpEdge(newybins*nygroup);
      resetStat = true; //stats must be reset because top bins will be moved to overflow bin
   }
   if (newzbins*nzgroup != nzbins) {
      zmax = fZaxis.GetBinUpEdge(newzbins*nzgroup);
      resetStat = true; //stats must be reset because top bins will be moved to overflow bin
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
   // save the TAttAxis members (reset by SetBins) for z axis
   Int_t    nZdivisions  = fZaxis.GetNdivisions();
   Color_t  zAxisColor   = fZaxis.GetAxisColor();
   Color_t  zLabelColor  = fZaxis.GetLabelColor();
   Style_t  zLabelFont   = fZaxis.GetLabelFont();
   Float_t  zLabelOffset = fZaxis.GetLabelOffset();
   Float_t  zLabelSize   = fZaxis.GetLabelSize();
   Float_t  zTickLength  = fZaxis.GetTickLength();
   Float_t  zTitleOffset = fZaxis.GetTitleOffset();
   Float_t  zTitleSize   = fZaxis.GetTitleSize();
   Color_t  zTitleColor  = fZaxis.GetTitleColor();
   Style_t  zTitleFont   = fZaxis.GetTitleFont();

   // copy merged bin contents (ignore under/overflows)
   if (nxgroup != 1 || nygroup != 1 || nzgroup != 1) {
      if (fXaxis.GetXbins()->GetSize() > 0 || fYaxis.GetXbins()->GetSize() > 0 || fZaxis.GetXbins()->GetSize() > 0) {
         // variable bin sizes in x or y, don't treat both cases separately
         Double_t *xbins = new Double_t[newxbins+1];
         for (i = 0; i <= newxbins; ++i) xbins[i] = fXaxis.GetBinLowEdge(1+i*nxgroup);
         Double_t *ybins = new Double_t[newybins+1];
         for (i = 0; i <= newybins; ++i) ybins[i] = fYaxis.GetBinLowEdge(1+i*nygroup);
         Double_t *zbins = new Double_t[newzbins+1];
         for (i = 0; i <= newzbins; ++i) zbins[i] = fZaxis.GetBinLowEdge(1+i*nzgroup);
         hnew->SetBins(newxbins,xbins, newybins, ybins, newzbins, zbins);//changes also errors array (if any)
         delete [] xbins;
         delete [] ybins;
         delete [] zbins;
      } else {
         hnew->SetBins(newxbins, xmin, xmax, newybins, ymin, ymax, newzbins, zmin, zmax);//changes also errors array
      }

      Double_t binContent, binSumw2;
      Int_t oldxbin = 1;
      Int_t oldybin = 1;
      Int_t oldzbin = 1;
      Int_t bin;
      for (xbin = 1; xbin <= newxbins; xbin++) {
         oldybin=1;
         oldzbin=1;
         for (ybin = 1; ybin <= newybins; ybin++) {
            oldzbin=1;
            for (zbin = 1; zbin <= newzbins; zbin++) {
               binContent = 0;
               binSumw2   = 0;
               for (i = 0; i < nxgroup; i++) {
                  if (oldxbin+i > nxbins) break;
                  for (j =0; j < nygroup; j++) {
                     if (oldybin+j > nybins) break;
                     for (k =0; k < nzgroup; k++) {
                        if (oldzbin+k > nzbins) break;
                        //get global bin (same conventions as in TH1::GetBin(xbin,ybin)
                        bin = oldxbin + i + (oldybin + j)*(nxbins + 2) + (oldzbin + k)*(nxbins + 2)*(nybins + 2);
                        binContent += oldBins[bin];
                        if (oldSumw2) binSumw2 += oldSumw2[bin];
                     }
                  }
               }
               Int_t ibin = hnew->GetBin(xbin,ybin,zbin);  // new bin number
               hnew->SetBinContent(ibin, binContent);
               if (oldSumw2) hnew->fSumw2.fArray[ibin] = binSumw2;
               oldzbin += nzgroup;
            }
            oldybin += nygroup;
         }
         oldxbin += nxgroup;
      }

      // compute new underflow/overflows for the 8 vertices
      for (Int_t xover = 0; xover <= 1; xover++) {
         for (Int_t yover = 0; yover <= 1; yover++) {
            for (Int_t zover = 0; zover <= 1; zover++) {
               binContent = 0;
               binSumw2 = 0;
               // make loop in case of only underflow/overflow
               for (xbin = xover*oldxbin; xbin <= xover*(nxbins+1); xbin++) {
                  for (ybin = yover*oldybin; ybin <= yover*(nybins+1); ybin++) {
                     for (zbin = zover*oldzbin; zbin <= zover*(nzbins+1); zbin++) {
                        bin = GetBin(xbin,ybin,zbin);
                        binContent += oldBins[bin];
                        if (oldSumw2) binSumw2 += oldSumw2[bin];
                     }
                  }
               }
               Int_t binNew = hnew->GetBin( xover *(newxbins+1),
                                           yover*(newybins+1), zover*(newzbins+1) );
               hnew->SetBinContent(binNew,binContent);
               if (oldSumw2) hnew->fSumw2.fArray[binNew] = binSumw2;
            }
         }
      }

      Double_t binContent0, binContent2, binContent3, binContent4;
      Double_t binError0, binError2, binError3, binError4;
      Int_t oldxbin2, oldybin2, oldzbin2;
      Int_t ufbin, ofbin, ofbin2, ofbin3, ofbin4;

      //  recompute under/overflow contents in y for the new  x and z bins
      oldxbin2 = 1;
      oldybin2 = 1;
      oldzbin2 = 1;
      for (xbin = 1; xbin<=newxbins; xbin++) {
         oldzbin2 = 1;
         for (zbin = 1; zbin<=newzbins; zbin++) {
            binContent0 = binContent2 = 0;
            binError0 = binError2 = 0;
            for (i=0; i<nxgroup; i++) {
               if (oldxbin2+i > nxbins) break;
               for (k=0; k<nzgroup; k++) {
                  if (oldzbin2+k > nzbins) break;
                  //old underflow bin (in y)
                  ufbin = oldxbin2 + i + (nxbins+2)*(nybins+2)*(oldzbin2+k);
                  binContent0 += oldBins[ufbin];
                  if (oldSumw2) binError0 += oldSumw2[ufbin];
                  for (ybin = oldybin; ybin <= nybins + 1; ybin++) {
                     //old overflow bin (in y)
                     ofbin = ufbin + ybin*(nxbins+2);
                     binContent2 += oldBins[ofbin];
                     if (oldSumw2) binError2 += oldSumw2[ofbin];
                  }
               }
            }
            hnew->SetBinContent(xbin,0,zbin,binContent0);
            hnew->SetBinContent(xbin,newybins+1,zbin,binContent2);
            if (oldSumw2) {
               hnew->SetBinError(xbin,0,zbin,TMath::Sqrt(binError0));
               hnew->SetBinError(xbin,newybins+1,zbin,TMath::Sqrt(binError2) );
            }
            oldzbin2 += nzgroup;
         }
         oldxbin2 += nxgroup;
      }

      //  recompute under/overflow contents in x for the new  y and z bins
      oldxbin2 = 1;
      oldybin2 = 1;
      oldzbin2 = 1;
      for (ybin = 1; ybin<=newybins; ybin++) {
         oldzbin2 = 1;
         for (zbin = 1; zbin<=newzbins; zbin++) {
            binContent0 = binContent2 = 0;
            binError0 = binError2 = 0;
            for (j=0; j<nygroup; j++) {
               if (oldybin2+j > nybins) break;
               for (k=0; k<nzgroup; k++) {
                  if (oldzbin2+k > nzbins) break;
                  //old underflow bin (in y)
                  ufbin = (oldybin2 + j)*(nxbins+2) + (nxbins+2)*(nybins+2)*(oldzbin2+k);
                  binContent0 += oldBins[ufbin];
                  if (oldSumw2) binError0 += oldSumw2[ufbin];
                  for (xbin = oldxbin; xbin <= nxbins + 1; xbin++) {
                     //old overflow bin (in x)
                     ofbin = ufbin + xbin;
                     binContent2 += oldBins[ofbin];
                     if (oldSumw2) binError2 += oldSumw2[ofbin];
                  }
               }
            }
            hnew->SetBinContent(0,ybin,zbin,binContent0);
            hnew->SetBinContent(newxbins+1,ybin,zbin,binContent2);
            if (oldSumw2) {
               hnew->SetBinError(0,ybin,zbin,TMath::Sqrt(binError0));
               hnew->SetBinError(newxbins+1,ybin,zbin,TMath::Sqrt(binError2) );
            }
            oldzbin2 += nzgroup;
         }
         oldybin2 += nygroup;
      }

      //  recompute under/overflow contents in z for the new  x and y bins
      oldxbin2 = 1;
      oldybin2 = 1;
      oldzbin2 = 1;
      for (xbin = 1; xbin<=newxbins; xbin++) {
         oldybin2 = 1;
         for (ybin = 1; ybin<=newybins; ybin++) {
            binContent0 = binContent2 = 0;
            binError0 = binError2 = 0;
            for (i=0; i<nxgroup; i++) {
               if (oldxbin2+i > nxbins) break;
               for (j=0; j<nygroup; j++) {
                  if (oldybin2+j > nybins) break;
                  //old underflow bin (in z)
                  ufbin = oldxbin2 + i + (nxbins+2)*(oldybin2+j);
                  binContent0 += oldBins[ufbin];
                  if (oldSumw2) binError0 += oldSumw2[ufbin];
                  for (zbin = oldzbin; zbin <= nzbins + 1; zbin++) {
                     //old overflow bin (in z)
                     ofbin = ufbin + (nxbins+2)*(nybins+2)*zbin;
                     binContent2 += oldBins[ofbin];
                     if (oldSumw2) binError2 += oldSumw2[ofbin];
                  }
               }
            }
            hnew->SetBinContent(xbin,ybin,0,binContent0);
            hnew->SetBinContent(xbin,ybin,newzbins+1,binContent2);
            if (oldSumw2) {
               hnew->SetBinError(xbin,ybin,0,TMath::Sqrt(binError0));
               hnew->SetBinError(xbin,ybin,newzbins+1,TMath::Sqrt(binError2) );
            }
            oldybin2 += nygroup;
         }
         oldxbin2 += nxgroup;
      }

      //  recompute under/overflow contents in y, z for the new  x
      oldxbin2 = 1;
      oldybin2 = 1;
      oldzbin2 = 1;
      for (xbin = 1; xbin<=newxbins; xbin++) {
         binContent0 = 0;
         binContent2 = 0;
         binContent3 = 0;
         binContent4 = 0;
         binError0 = 0;
         binError2 = 0;
         binError3 = 0;
         binError4 = 0;
         for (i=0; i<nxgroup; i++) {
            if (oldxbin2+i > nxbins) break;
            ufbin = oldxbin2 + i; //
            binContent0 += oldBins[ufbin];
            if (oldSumw2) binError0 += oldSumw2[ufbin];
            for (ybin = oldybin; ybin <= nybins + 1; ybin++) {
               ofbin3 =  ufbin+ybin*(nxbins+2);
               binContent3 += oldBins[ ofbin3 ];
               if (oldSumw2)  binError3 += oldSumw2[ofbin3];
               for (zbin = oldzbin; zbin <= nzbins + 1; zbin++) {
                  //old overflow bin (in z)
                  ofbin4 =   oldxbin2 + i + ybin*(nxbins+2) + (nxbins+2)*(nybins+2)*zbin;
                  binContent4 += oldBins[ofbin4];
                  if (oldSumw2) binError4 += oldSumw2[ofbin4];
               }
            }
            for (zbin = oldzbin; zbin <= nzbins + 1; zbin++) {
               ofbin2 =  ufbin+zbin*(nxbins+2)*(nybins+2);
               binContent2 += oldBins[ ofbin2 ];
               if (oldSumw2)  binError2 += oldSumw2[ofbin2];
            }
         }
         hnew->SetBinContent(xbin,0,0,binContent0);
         hnew->SetBinContent(xbin,0,newzbins+1,binContent2);
         hnew->SetBinContent(xbin,newybins+1,0,binContent3);
         hnew->SetBinContent(xbin,newybins+1,newzbins+1,binContent4);
         if (oldSumw2) {
            hnew->SetBinError(xbin,0,0,TMath::Sqrt(binError0));
            hnew->SetBinError(xbin,0,newzbins+1,TMath::Sqrt(binError2) );
            hnew->SetBinError(xbin,newybins+1,0,TMath::Sqrt(binError3) );
            hnew->SetBinError(xbin,newybins+1,newzbins+1,TMath::Sqrt(binError4) );
         }
         oldxbin2 += nxgroup;
      }

      //  recompute under/overflow contents in x, y for the new z
      oldxbin2 = 1;
      oldybin2 = 1;
      oldzbin2 = 1;
      for (zbin = 1; zbin<=newzbins; zbin++) {
         binContent0 = 0;
         binContent2 = 0;
         binContent3 = 0;
         binContent4 = 0;
         binError0 = 0;
         binError2 = 0;
         binError3 = 0;
         binError4 = 0;
         for (i=0; i<nzgroup; i++) {
            if (oldzbin2+i > nzbins) break;
            ufbin = (oldzbin2 + i)*(nxbins+2)*(nybins+2); //
            binContent0 += oldBins[ufbin];
            if (oldSumw2) binError0 += oldSumw2[ufbin];
            for (ybin = oldybin; ybin <= nybins + 1; ybin++) {
               ofbin3 =  ufbin+ybin*(nxbins+2);
               binContent3 += oldBins[ ofbin3 ];
               if (oldSumw2)  binError3 += oldSumw2[ofbin3];
               for (xbin = oldxbin; xbin <= nxbins + 1; xbin++) {
                  //old overflow bin (in z)
                  ofbin4 = ufbin + xbin + ybin*(nxbins+2);
                  binContent4 += oldBins[ofbin4];
                  if (oldSumw2) binError4 += oldSumw2[ofbin4];
               }
            }
            for (xbin = oldxbin; xbin <= nxbins + 1; xbin++) {
               ofbin2 =  xbin +(oldzbin2+i)*(nxbins+2)*(nybins+2);
               binContent2 += oldBins[ ofbin2 ];
               if (oldSumw2)  binError2 += oldSumw2[ofbin2];
            }
         }
         hnew->SetBinContent(0,0,zbin,binContent0);
         hnew->SetBinContent(0,newybins+1,zbin,binContent3);
         hnew->SetBinContent(newxbins+1,0,zbin,binContent2);
         hnew->SetBinContent(newxbins+1,newybins+1,zbin,binContent4);
         if (oldSumw2) {
            hnew->SetBinError(0,0,zbin,TMath::Sqrt(binError0));
            hnew->SetBinError(0,newybins+1,zbin,TMath::Sqrt(binError3) );
            hnew->SetBinError(newxbins+1,0,zbin,TMath::Sqrt(binError2) );
            hnew->SetBinError(newxbins+1,newybins+1,zbin,TMath::Sqrt(binError4) );
         }
         oldzbin2 += nzgroup;
      }

      //  recompute under/overflow contents in x, z for the new  y
      oldxbin2 = 1;
      oldybin2 = 1;
      oldzbin2 = 1;
      for (ybin = 1; ybin<=newybins; ybin++) {
         binContent0 = 0;
         binContent2 = 0;
         binContent3 = 0;
         binContent4 = 0;
         binError0 = 0;
         binError2 = 0;
         binError3 = 0;
         binError4 = 0;
         for (i=0; i<nygroup; i++) {
            if (oldybin2+i > nybins) break;
            ufbin = (oldybin2 + i)*(nxbins+2); //
            binContent0 += oldBins[ufbin];
            if (oldSumw2) binError0 += oldSumw2[ufbin];
            for (xbin = oldxbin; xbin <= nxbins + 1; xbin++) {
               ofbin3 =  ufbin+xbin;
               binContent3 += oldBins[ ofbin3 ];
               if (oldSumw2)  binError3 += oldSumw2[ofbin3];
               for (zbin = oldzbin; zbin <= nzbins + 1; zbin++) {
                  //old overflow bin (in z)
                  ofbin4 = xbin + (nxbins+2)*(nybins+2)*zbin+(oldybin2+i)*(nxbins+2);
                  binContent4 += oldBins[ofbin4];
                  if (oldSumw2) binError4 += oldSumw2[ofbin4];
               }
            }
            for (zbin = oldzbin; zbin <= nzbins + 1; zbin++) {
               ofbin2 =  (oldybin2+i)*(nxbins+2)+zbin*(nxbins+2)*(nybins+2);
               binContent2 += oldBins[ ofbin2 ];
               if (oldSumw2)  binError2 += oldSumw2[ofbin2];
            }
         }
         hnew->SetBinContent(0,ybin,0,binContent0);
         hnew->SetBinContent(0,ybin,newzbins+1,binContent2);
         hnew->SetBinContent(newxbins+1,ybin,0,binContent3);
         hnew->SetBinContent(newxbins+1,ybin,newzbins+1,binContent4);
         if (oldSumw2) {
            hnew->SetBinError(0,ybin,0,TMath::Sqrt(binError0));
            hnew->SetBinError(0,ybin,newzbins+1,TMath::Sqrt(binError2) );
            hnew->SetBinError(newxbins+1,ybin,0,TMath::Sqrt(binError3) );
            hnew->SetBinError(newxbins+1,ybin,newzbins+1,TMath::Sqrt(binError4) );
         }
         oldybin2 += nygroup;
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
   // Restore z axis attributes
   fZaxis.SetNdivisions(nZdivisions);
   fZaxis.SetAxisColor(zAxisColor);
   fZaxis.SetLabelColor(zLabelColor);
   fZaxis.SetLabelFont(zLabelFont);
   fZaxis.SetLabelOffset(zLabelOffset);
   fZaxis.SetLabelSize(zLabelSize);
   fZaxis.SetTickLength(zTickLength);
   fZaxis.SetTitleOffset(zTitleOffset);
   fZaxis.SetTitleSize(zTitleSize);
   fZaxis.SetTitleColor(zTitleColor);
   fZaxis.SetTitleFont(zTitleFont);

   //restore statistics and entries  modified by SetBinContent
   hnew->SetEntries(entries);
   if (!resetStat) hnew->PutStats(stat);

   delete [] oldBins;
   if (oldSumw2) delete [] oldSumw2;
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Reset this histogram: contents, errors, etc.

void TH3::Reset(Option_t *option)
{
   TH1::Reset(option);
   TString opt = option;
   opt.ToUpper();
   if (opt.Contains("ICE") && !opt.Contains("S")) return;
   fTsumwy  = 0;
   fTsumwy2 = 0;
   fTsumwxy = 0;
   fTsumwz  = 0;
   fTsumwz2 = 0;
   fTsumwxz = 0;
   fTsumwyz = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Set bin content.

void TH3::SetBinContent(Int_t bin, Double_t content)
{
   fEntries++;
   fTsumw = 0;
   if (bin < 0) return;
   if (bin >= fNcells) return;
   UpdateBinContent(bin, content);
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TH3.

void TH3::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TH3::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TH1::Streamer(R__b);
      TAtt3D::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TH3::IsA());
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TH3::Class(),this);
   }
}


//______________________________________________________________________________
//                     TH3C methods
//  TH3C a 3-D histogram with one byte per cell (char)
//______________________________________________________________________________

ClassImp(TH3C);


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH3C::TH3C(): TH3(), TArrayC()
{
   SetBinsLength(27);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH3C::~TH3C()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for fix bin size 3-D histograms.

TH3C::TH3C(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,Double_t ylow,Double_t yup
           ,Int_t nbinsz,Double_t zlow,Double_t zup)
           :TH3(name,title,nbinsx,xlow,xup,nbinsy,ylow,yup,nbinsz,zlow,zup)
{
   TArrayC::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();

   if (xlow >= xup || ylow >= yup || zlow >= zup) SetBuffer(fgBufferSize);
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for variable bin size 3-D histograms.

TH3C::TH3C(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins
           ,Int_t nbinsy,const Float_t *ybins
           ,Int_t nbinsz,const Float_t *zbins)
           :TH3(name,title,nbinsx,xbins,nbinsy,ybins,nbinsz,zbins)
{
   TArrayC::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for variable bin size 3-D histograms.

TH3C::TH3C(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,const Double_t *ybins
           ,Int_t nbinsz,const Double_t *zbins)
           :TH3(name,title,nbinsx,xbins,nbinsy,ybins,nbinsz,zbins)
{
   TArrayC::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TH3C::TH3C(const TH3C &h3c) : TH3(), TArrayC()
{
   ((TH3C&)h3c).Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by 1.

void TH3C::AddBinContent(Int_t bin)
{
   if (fArray[bin] < 127) fArray[bin]++;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by w.

void TH3C::AddBinContent(Int_t bin, Double_t w)
{
   Int_t newval = fArray[bin] + Int_t(w);
   if (newval > -128 && newval < 128) {fArray[bin] = Char_t(newval); return;}
   if (newval < -127) fArray[bin] = -127;
   if (newval >  127) fArray[bin] =  127;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy this 3-D histogram structure to newth3.

void TH3C::Copy(TObject &newth3) const
{
   TH3::Copy((TH3C&)newth3);
}


////////////////////////////////////////////////////////////////////////////////
/// Reset this histogram: contents, errors, etc.

void TH3C::Reset(Option_t *option)
{
   TH3::Reset(option);
   TArrayC::Reset();
   // should also reset statistics once statistics are implemented for TH3
}


////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH3C::SetBinsLength(Int_t n)
{
   if (n < 0) n = (fXaxis.GetNbins()+2)*(fYaxis.GetNbins()+2)*(fZaxis.GetNbins()+2);
   fNcells = n;
   TArrayC::Set(n);
}


////////////////////////////////////////////////////////////////////////////////
/// When the mouse is moved in a pad containing a 3-d view of this histogram
/// a second canvas shows a projection type given as option.
/// To stop the generation of the projections, delete the canvas
/// containing the projection.
/// option may contain a combination of the characters x,y,z,e
/// option = "x" return the x projection into a TH1D histogram
/// option = "y" return the y projection into a TH1D histogram
/// option = "z" return the z projection into a TH1D histogram
/// option = "xy" return the x versus y projection into a TH2D histogram
/// option = "yx" return the y versus x projection into a TH2D histogram
/// option = "xz" return the x versus z projection into a TH2D histogram
/// option = "zx" return the z versus x projection into a TH2D histogram
/// option = "yz" return the y versus z projection into a TH2D histogram
/// option = "zy" return the z versus y projection into a TH2D histogram
/// option can also include the drawing option for the projection, eg to draw
/// the xy projection using the draw option "box" do
///   myhist.SetShowProjection("xy box");
/// This function is typically called from the context menu.
/// NB: the notation "a vs b" means "a" vertical and "b" horizontal

void TH3::SetShowProjection(const char *option,Int_t nbins)
{
   GetPainter();

   if (fPainter) fPainter->SetShowProjection(option,nbins);
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TH3C.

void TH3C::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      if (R__b.GetParent() && R__b.GetVersionOwner() < 22300) return;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TH3C::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      if (R__v < 2) {
         R__b.ReadVersion();
         TH1::Streamer(R__b);
         TArrayC::Streamer(R__b);
         R__b.ReadVersion(&R__s, &R__c);
         TAtt3D::Streamer(R__b);
      } else {
         TH3::Streamer(R__b);
         TArrayC::Streamer(R__b);
         R__b.CheckByteCount(R__s, R__c, TH3C::IsA());
      }
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TH3C::Class(),this);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH3C& TH3C::operator=(const TH3C &h1)
{
   if (this != &h1)  ((TH3C&)h1).Copy(*this);
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH3C operator*(Float_t c1, TH3C &h1)
{
   TH3C hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH3C operator+(TH3C &h1, TH3C &h2)
{
   TH3C hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator -

TH3C operator-(TH3C &h1, TH3C &h2)
{
   TH3C hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH3C operator*(TH3C &h1, TH3C &h2)
{
   TH3C hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH3C operator/(TH3C &h1, TH3C &h2)
{
   TH3C hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


//______________________________________________________________________________
//                     TH3S methods
//  TH3S a 3-D histogram with two bytes per cell (short integer)
//______________________________________________________________________________

ClassImp(TH3S);


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH3S::TH3S(): TH3(), TArrayS()
{
   SetBinsLength(27);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH3S::~TH3S()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for fix bin size 3-D histograms.

TH3S::TH3S(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,Double_t ylow,Double_t yup
           ,Int_t nbinsz,Double_t zlow,Double_t zup)
           :TH3(name,title,nbinsx,xlow,xup,nbinsy,ylow,yup,nbinsz,zlow,zup)
{
   TH3S::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();

   if (xlow >= xup || ylow >= yup || zlow >= zup) SetBuffer(fgBufferSize);
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for variable bin size 3-D histograms.

TH3S::TH3S(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins
           ,Int_t nbinsy,const Float_t *ybins
           ,Int_t nbinsz,const Float_t *zbins)
           :TH3(name,title,nbinsx,xbins,nbinsy,ybins,nbinsz,zbins)
{
   TH3S::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for variable bin size 3-D histograms.

TH3S::TH3S(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,const Double_t *ybins
           ,Int_t nbinsz,const Double_t *zbins)
           :TH3(name,title,nbinsx,xbins,nbinsy,ybins,nbinsz,zbins)
{
   TH3S::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Copy Constructor.

TH3S::TH3S(const TH3S &h3s) : TH3(), TArrayS()
{
   ((TH3S&)h3s).Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by 1.

void TH3S::AddBinContent(Int_t bin)
{
   if (fArray[bin] < 32767) fArray[bin]++;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by w.

void TH3S::AddBinContent(Int_t bin, Double_t w)
{
   Int_t newval = fArray[bin] + Int_t(w);
   if (newval > -32768 && newval < 32768) {fArray[bin] = Short_t(newval); return;}
   if (newval < -32767) fArray[bin] = -32767;
   if (newval >  32767) fArray[bin] =  32767;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy this 3-D histogram structure to newth3.

void TH3S::Copy(TObject &newth3) const
{
   TH3::Copy((TH3S&)newth3);
}


////////////////////////////////////////////////////////////////////////////////
/// Reset this histogram: contents, errors, etc.

void TH3S::Reset(Option_t *option)
{
   TH3::Reset(option);
   TArrayS::Reset();
   // should also reset statistics once statistics are implemented for TH3
}


////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH3S::SetBinsLength(Int_t n)
{
   if (n < 0) n = (fXaxis.GetNbins()+2)*(fYaxis.GetNbins()+2)*(fZaxis.GetNbins()+2);
   fNcells = n;
   TArrayS::Set(n);
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TH3S.

void TH3S::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      if (R__b.GetParent() && R__b.GetVersionOwner() < 22300) return;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TH3S::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      if (R__v < 2) {
         R__b.ReadVersion();
         TH1::Streamer(R__b);
         TArrayS::Streamer(R__b);
         R__b.ReadVersion(&R__s, &R__c);
         TAtt3D::Streamer(R__b);
      } else {
         TH3::Streamer(R__b);
         TArrayS::Streamer(R__b);
         R__b.CheckByteCount(R__s, R__c, TH3S::IsA());
      }
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TH3S::Class(),this);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH3S& TH3S::operator=(const TH3S &h1)
{
   if (this != &h1)  ((TH3S&)h1).Copy(*this);
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH3S operator*(Float_t c1, TH3S &h1)
{
   TH3S hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH3S operator+(TH3S &h1, TH3S &h2)
{
   TH3S hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator -

TH3S operator-(TH3S &h1, TH3S &h2)
{
   TH3S hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH3S operator*(TH3S &h1, TH3S &h2)
{
   TH3S hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH3S operator/(TH3S &h1, TH3S &h2)
{
   TH3S hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


//______________________________________________________________________________
//                     TH3I methods
//  TH3I a 3-D histogram with four bytes per cell (32 bits integer)
//______________________________________________________________________________

ClassImp(TH3I);


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH3I::TH3I(): TH3(), TArrayI()
{
   SetBinsLength(27);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH3I::~TH3I()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for fix bin size 3-D histograms.

TH3I::TH3I(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,Double_t ylow,Double_t yup
           ,Int_t nbinsz,Double_t zlow,Double_t zup)
           :TH3(name,title,nbinsx,xlow,xup,nbinsy,ylow,yup,nbinsz,zlow,zup)
{
   TH3I::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();

   if (xlow >= xup || ylow >= yup || zlow >= zup) SetBuffer(fgBufferSize);
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for variable bin size 3-D histograms.

TH3I::TH3I(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins
           ,Int_t nbinsy,const Float_t *ybins
           ,Int_t nbinsz,const Float_t *zbins)
           :TH3(name,title,nbinsx,xbins,nbinsy,ybins,nbinsz,zbins)
{
   TArrayI::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for variable bin size 3-D histograms.

TH3I::TH3I(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,const Double_t *ybins
           ,Int_t nbinsz,const Double_t *zbins)
           :TH3(name,title,nbinsx,xbins,nbinsy,ybins,nbinsz,zbins)
{
   TArrayI::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TH3I::TH3I(const TH3I &h3i) : TH3(), TArrayI()
{
   ((TH3I&)h3i).Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by 1.

void TH3I::AddBinContent(Int_t bin)
{
   if (fArray[bin] < 2147483647) fArray[bin]++;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by w.

void TH3I::AddBinContent(Int_t bin, Double_t w)
{
   Long64_t newval = fArray[bin] + Long64_t(w);
   if (newval > -2147483647 && newval < 2147483647) {fArray[bin] = Int_t(newval); return;}
   if (newval < -2147483647) fArray[bin] = -2147483647;
   if (newval >  2147483647) fArray[bin] =  2147483647;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy this 3-D histogram structure to newth3.

void TH3I::Copy(TObject &newth3) const
{
   TH3::Copy((TH3I&)newth3);
}


////////////////////////////////////////////////////////////////////////////////
/// Reset this histogram: contents, errors, etc.

void TH3I::Reset(Option_t *option)
{
   TH3::Reset(option);
   TArrayI::Reset();
   // should also reset statistics once statistics are implemented for TH3
}


////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH3I::SetBinsLength(Int_t n)
{
   if (n < 0) n = (fXaxis.GetNbins()+2)*(fYaxis.GetNbins()+2)*(fZaxis.GetNbins()+2);
   fNcells = n;
   TArrayI::Set(n);
}


////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH3I& TH3I::operator=(const TH3I &h1)
{
   if (this != &h1)  ((TH3I&)h1).Copy(*this);
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH3I operator*(Float_t c1, TH3I &h1)
{
   TH3I hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH3I operator+(TH3I &h1, TH3I &h2)
{
   TH3I hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator _

TH3I operator-(TH3I &h1, TH3I &h2)
{
   TH3I hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH3I operator*(TH3I &h1, TH3I &h2)
{
   TH3I hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH3I operator/(TH3I &h1, TH3I &h2)
{
   TH3I hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


//______________________________________________________________________________
//                     TH3F methods
//  TH3F a 3-D histogram with four bytes per cell (float)
//______________________________________________________________________________

ClassImp(TH3F);


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH3F::TH3F(): TH3(), TArrayF()
{
   SetBinsLength(27);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH3F::~TH3F()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for fix bin size 3-D histograms.

TH3F::TH3F(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,Double_t ylow,Double_t yup
           ,Int_t nbinsz,Double_t zlow,Double_t zup)
           :TH3(name,title,nbinsx,xlow,xup,nbinsy,ylow,yup,nbinsz,zlow,zup)
{
   TArrayF::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();

   if (xlow >= xup || ylow >= yup || zlow >= zup) SetBuffer(fgBufferSize);
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for variable bin size 3-D histograms.

TH3F::TH3F(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins
           ,Int_t nbinsy,const Float_t *ybins
           ,Int_t nbinsz,const Float_t *zbins)
           :TH3(name,title,nbinsx,xbins,nbinsy,ybins,nbinsz,zbins)
{
   TArrayF::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for variable bin size 3-D histograms.

TH3F::TH3F(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,const Double_t *ybins
           ,Int_t nbinsz,const Double_t *zbins)
           :TH3(name,title,nbinsx,xbins,nbinsy,ybins,nbinsz,zbins)
{
   TArrayF::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TH3F::TH3F(const TH3F &h3f) : TH3(), TArrayF()
{
   ((TH3F&)h3f).Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Copy this 3-D histogram structure to newth3.

void TH3F::Copy(TObject &newth3) const
{
   TH3::Copy((TH3F&)newth3);
}


////////////////////////////////////////////////////////////////////////////////
/// Reset this histogram: contents, errors, etc.

void TH3F::Reset(Option_t *option)
{
   TH3::Reset(option);
   TArrayF::Reset();
   // should also reset statistics once statistics are implemented for TH3
}


////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH3F::SetBinsLength(Int_t n)
{
   if (n < 0) n = (fXaxis.GetNbins()+2)*(fYaxis.GetNbins()+2)*(fZaxis.GetNbins()+2);
   fNcells = n;
   TArrayF::Set(n);
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TH3F.

void TH3F::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      if (R__b.GetParent() && R__b.GetVersionOwner() < 22300) return;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TH3F::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      if (R__v < 2) {
         R__b.ReadVersion();
         TH1::Streamer(R__b);
         TArrayF::Streamer(R__b);
         R__b.ReadVersion(&R__s, &R__c);
         TAtt3D::Streamer(R__b);
      } else {
         TH3::Streamer(R__b);
         TArrayF::Streamer(R__b);
         R__b.CheckByteCount(R__s, R__c, TH3F::IsA());
      }
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TH3F::Class(),this);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH3F& TH3F::operator=(const TH3F &h1)
{
   if (this != &h1)  ((TH3F&)h1).Copy(*this);
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH3F operator*(Float_t c1, TH3F &h1)
{
   TH3F hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH3F operator+(TH3F &h1, TH3F &h2)
{
   TH3F hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator -

TH3F operator-(TH3F &h1, TH3F &h2)
{
   TH3F hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH3F operator*(TH3F &h1, TH3F &h2)
{
   TH3F hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH3F operator/(TH3F &h1, TH3F &h2)
{
   TH3F hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


//______________________________________________________________________________
//                     TH3D methods
//  TH3D a 3-D histogram with eight bytes per cell (double)
//______________________________________________________________________________

ClassImp(TH3D);


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH3D::TH3D(): TH3(), TArrayD()
{
   SetBinsLength(27);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH3D::~TH3D()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for fix bin size 3-D histograms.

TH3D::TH3D(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
           ,Int_t nbinsy,Double_t ylow,Double_t yup
           ,Int_t nbinsz,Double_t zlow,Double_t zup)
           :TH3(name,title,nbinsx,xlow,xup,nbinsy,ylow,yup,nbinsz,zlow,zup)
{
   TArrayD::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();

   if (xlow >= xup || ylow >= yup || zlow >= zup) SetBuffer(fgBufferSize);
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for variable bin size 3-D histograms.

TH3D::TH3D(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins
           ,Int_t nbinsy,const Float_t *ybins
           ,Int_t nbinsz,const Float_t *zbins)
           :TH3(name,title,nbinsx,xbins,nbinsy,ybins,nbinsz,zbins)
{
   TArrayD::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for variable bin size 3-D histograms.

TH3D::TH3D(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
           ,Int_t nbinsy,const Double_t *ybins
           ,Int_t nbinsz,const Double_t *zbins)
           :TH3(name,title,nbinsx,xbins,nbinsy,ybins,nbinsz,zbins)
{
   TArrayD::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TH3D::TH3D(const TH3D &h3d) : TH3(), TArrayD()
{
   ((TH3D&)h3d).Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Copy this 3-D histogram structure to newth3.

void TH3D::Copy(TObject &newth3) const
{
   TH3::Copy((TH3D&)newth3);
}


////////////////////////////////////////////////////////////////////////////////
/// Reset this histogram: contents, errors, etc.

void TH3D::Reset(Option_t *option)
{
   TH3::Reset(option);
   TArrayD::Reset();
   // should also reset statistics once statistics are implemented for TH3
}


////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH3D::SetBinsLength(Int_t n)
{
   if (n < 0) n = (fXaxis.GetNbins()+2)*(fYaxis.GetNbins()+2)*(fZaxis.GetNbins()+2);
   fNcells = n;
   TArrayD::Set(n);
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TH3D.

void TH3D::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      if (R__b.GetParent() && R__b.GetVersionOwner() < 22300) return;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TH3D::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      if (R__v < 2) {
         R__b.ReadVersion();
         TH1::Streamer(R__b);
         TArrayD::Streamer(R__b);
         R__b.ReadVersion(&R__s, &R__c);
         TAtt3D::Streamer(R__b);
      } else {
         TH3::Streamer(R__b);
         TArrayD::Streamer(R__b);
         R__b.CheckByteCount(R__s, R__c, TH3D::IsA());
      }
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TH3D::Class(),this);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH3D& TH3D::operator=(const TH3D &h1)
{
   if (this != &h1)  ((TH3D&)h1).Copy(*this);
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH3D operator*(Float_t c1, TH3D &h1)
{
   TH3D hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH3D operator+(TH3D &h1, TH3D &h2)
{
   TH3D hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator -

TH3D operator-(TH3D &h1, TH3D &h2)
{
   TH3D hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH3D operator*(TH3D &h1, TH3D &h2)
{
   TH3D hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH3D operator/(TH3D &h1, TH3D &h2)
{
   TH3D hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}
