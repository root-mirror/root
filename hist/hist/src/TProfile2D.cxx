// @(#)root/hist:$Id$
// Author: Rene Brun   16/04/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TProfile2D.h"
#include "TBuffer.h"
#include "TMath.h"
#include "THLimitsFinder.h"
#include "TError.h"
#include "TClass.h"
#include "TProfileHelper.h"
#include <iostream>

Bool_t TProfile2D::fgApproximate = kFALSE;

ClassImp(TProfile2D);

/** \class TProfile2D
    \ingroup Hist
 Profile2D histograms are used to display the mean
 value of Z and its error for each cell in X,Y.
 Profile2D histograms are in many cases an
 elegant replacement of three-dimensional histograms : the inter-relation of three
 measured quantities X, Y and Z can always be visualized by a three-dimensional
 histogram or scatter-plot; its representation on the line-printer is not particularly
 satisfactory, except for sparse data. If Z is an unknown (but single-valued)
 approximate function of X,Y this function is displayed by a profile2D histogram with
 much better precision than by a scatter-plot.

 The following formulae show the cumulated contents (capital letters) and the values
 displayed by the printing or plotting routines (small letters) of the elements for cell i, j.
 \f[
  \begin{align}
       H(i,j)  &=  \sum w \cdot Z  \\
       E(i,j)  &=  \sum w \cdot Z^2 \\
       W(i,j)  &=  \sum w \\
       h(i,j)  &=  \frac{H(i,j)}{W(i,j)} \\
       s(i,j)  &=  \sqrt{E(i,j)/W(i,j)- h(i,j)^2} \\
       e(i,j)  &=  \frac{s(i,j)}{\sqrt{W(i,j)}}
  \end{align}
 \f]
  The bin content is always the mean of the Z values, but errors change depending on options:
 \f[
    \begin{align}
      \text{GetBinContent}(i,j) &= h(i,j) \\
      \text{GetBinError}(i,j) &=
        \begin{cases}
          e(i,j)                 &\text{if option="" (default). Error of the mean of all z values.} \\
          s(i,j)                 &\text{if option="s". Standard deviation of z values.} \\
          \begin{cases} e(j) &\text{if } h(i,j) \ne 0 \\ 1/\sqrt{12 N} &\text{if } h(i,j)=0 \end{cases}       &\text{if option="i". This is useful for storing integers such as ADC counts.} \\
          1/\sqrt{W(i,j)}           &\text{if option="g". Error of a weighted mean when combining measurements with variances of } w. \\
        \end{cases}
    \end{align}
 \f]

 In the special case where s(I,J) is zero (eg, case of 1 entry only in one cell)
 the bin error e(I,J) is computed from the average of the s(I,J) for all cells
 if the static function TProfile2D::Approximate has been called.
 This simple/crude approximation was suggested in order to keep the cell
 during a fit operation. But note that this approximation is not the default behaviour.

 ### Creating and drawing a 2D profile
 ~~~~{.cpp}
 {
    auto c1 = new TCanvas("c1","Profile histogram example",200,10,700,500);
    auto hprof2d  = new TProfile2D("hprof2d","Profile of pz versus px and py",40,-4,4,40,-4,4,0,20);
    Float_t px, py, pz;
    for ( Int_t i=0; i<25000; i++) {
       gRandom->Rannor(px,py);
       pz = px*px + py*py;
       hprof2d->Fill(px,py,pz,1);
    }
    hprof2d->Draw();
 }
 ~~~~
*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for Profile2D histograms.

TProfile2D::TProfile2D() : TH2D()
{
   fTsumwz = fTsumwz2 = 0;
   fScaling = kFALSE;
   BuildOptions(0,0,"");
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor for Profile2D histograms.

TProfile2D::~TProfile2D()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Normal Constructor for Profile histograms.
///
/// The first eight parameters are similar to TH2D::TH2D.
/// All values of z are accepted at filling time.
/// To fill a profile2D histogram, one must use TProfile2D::Fill function.
///
/// Note that when filling the profile histogram the function Fill
/// checks if the variable z is between fZmin and fZmax.
/// If a minimum or maximum value is set for the Z scale before filling,
/// then all values below zmin or above zmax will be discarded.
/// Setting the minimum or maximum value for the Z scale before filling
/// has the same effect as calling the special TProfile2D constructor below
/// where zmin and zmax are specified.
///
/// H(I,J) is printed as the cell contents. The errors computed are s(I,J) if CHOPT='S'
/// (spread option), or e(I,J) if CHOPT=' ' (error on mean).
///
/// See TProfile2D::BuildOptions for explanation of parameters
///
/// see other constructors below with all possible combinations of
/// fix and variable bin size like in TH2D.

TProfile2D::TProfile2D(const char *name,const char *title,Int_t nx,Double_t xlow,Double_t xup,Int_t ny,Double_t ylow,Double_t yup,Option_t *option)
: TH2D(name,title,nx,xlow,xup,ny,ylow,yup)
{
   BuildOptions(0,0,option);
   if (xlow >= xup || ylow >= yup) SetBuffer(fgBufferSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 2-D Profile with variable bins in X and fix bins in Y.

TProfile2D::TProfile2D(const char *name,const char *title,Int_t nx,const Double_t *xbins,Int_t ny,Double_t ylow,Double_t yup,Option_t *option)
: TH2D(name,title,nx,xbins,ny,ylow,yup)
{
   BuildOptions(0,0,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 2-D Profile with fix bins in X and variable bins in Y.

TProfile2D::TProfile2D(const char *name,const char *title,Int_t nx,Double_t xlow,Double_t xup,Int_t ny,const Double_t *ybins,Option_t *option)
: TH2D(name,title,nx,xlow,xup,ny,ybins)
{
   BuildOptions(0,0,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 2-D Profile with variable bins in X and variable bins in Y.

TProfile2D::TProfile2D(const char *name,const char *title,Int_t nx,const Double_t *xbins,Int_t ny,const Double_t *ybins,Option_t *option)
: TH2D(name,title,nx,xbins,ny,ybins)
{
   BuildOptions(0,0,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor for Profile2D histograms with range in z.
///
/// The first eight parameters are similar to TH2D::TH2D.
/// Only the values of Z between ZMIN and ZMAX will be considered at filling time.
/// zmin and zmax will also be the maximum and minimum values
/// on the z scale when drawing the profile2D.
///
/// See TProfile2D::BuildOptions for more explanations on errors

TProfile2D::TProfile2D(const char *name,const char *title,Int_t nx,Double_t xlow,Double_t xup,Int_t ny, Double_t ylow,Double_t yup,Double_t zlow,Double_t zup,Option_t *option)
: TH2D(name,title,nx,xlow,xup,ny,ylow,yup)
{
   BuildOptions(zlow,zup,option);
   if (xlow >= xup || ylow >= yup) SetBuffer(fgBufferSize);
}


////////////////////////////////////////////////////////////////////////////////
/// Set Profile2D histogram structure and options.
///
///  - zmin:  minimum value allowed for z
///  - zmax:  maximum value allowed for z
///            if (zmin = zmax = 0) there are no limits on the allowed z values (zmin = -inf, zmax = +inf)
///
///  - option:  this is the option for the computation of the t error of the profile ( TProfile2D::GetBinError )
///             possible values for the options are documented in TProfile2D::SetErrorOption
///
///   See TProfile::BuildOptions  for a detailed  description

void TProfile2D::BuildOptions(Double_t zmin, Double_t zmax, Option_t *option)
{

   SetErrorOption(option);

   // create extra profile data structure (bin entries/ y^2 and sum of weight square)
   TProfileHelper::BuildArray(this);

   fZmin = zmin;
   fZmax = zmax;
   fScaling = kFALSE;
   fTsumwz  = fTsumwz2 = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TProfile2D::TProfile2D(const TProfile2D &profile) : TH2D()
{
   ((TProfile2D&)profile).Copy(*this);
}

TProfile2D &TProfile2D::operator=(const TProfile2D &profile)
{
   ((TProfile2D &)profile).Copy(*this);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: `this = this + c1*f1` .

Bool_t TProfile2D::Add(TF1 *, Double_t , Option_t*)
{
   Error("Add","Function not implemented for TProfile2D");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: `this = this + c1*h1` .

Bool_t TProfile2D::Add(const TH1 *h1, Double_t c1)
{
   if (!h1) {
      Error("Add","Attempt to add a non-existing profile");
      return  kFALSE;
   }
   if (!h1->InheritsFrom(TProfile2D::Class())) {
      Error("Add","Attempt to add a non-profile2D object");
      return  kFALSE;
   }

   return TProfileHelper::Add(this, this, h1, 1, c1);
}

////////////////////////////////////////////////////////////////////////////////
/// Replace contents of this profile2D by the addition of h1 and h2.
///
/// `this = c1*h1 + c2*h2`

Bool_t TProfile2D::Add(const TH1 *h1, const TH1 *h2, Double_t c1, Double_t c2)
{
   if (!h1 || !h2) {
      Error("Add","Attempt to add a non-existing profile");
      return kFALSE;
   }
   if (!h1->InheritsFrom(TProfile2D::Class())) {
      Error("Add","Attempt to add a non-profile2D object");
      return kFALSE;
   }
   if (!h2->InheritsFrom(TProfile2D::Class())) {
      Error("Add","Attempt to add a non-profile2D object");
      return kFALSE;
   }
   return TProfileHelper::Add(this, h1, h2, c1, c2);
}

////////////////////////////////////////////////////////////////////////////////
/// Static function, set the fgApproximate flag.
///
/// When the flag is true, the function GetBinError
/// will approximate the bin error with the average profile error on all bins
/// in the following situation only
///  - the number of bins in the profile2D is less than 10404 (eg 100x100)
///  - the bin number of entries is small ( <5)
///  - the estimated bin error is extremely small compared to the bin content
///  (see TProfile2D::GetBinError)

void TProfile2D::Approximate(Bool_t approx)
{
   fgApproximate = approx;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill histogram with all entries in the buffer.
///
/// - action = -1 histogram is reset and refilled from the buffer (called by THistPainter::Paint)
/// - action =  0 histogram is filled from the buffer
/// - action =  1 histogram is filled and buffer is deleted
///             The buffer is automatically deleted when the number of entries
///             in the buffer is greater than the number of entries in the histogram

Int_t TProfile2D::BufferEmpty(Int_t action)
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
      Reset("ICES"); // reset without deleting the functions
      fBuffer = buffer;
   }
   if (CanExtendAllAxes() || fXaxis.GetXmax() <= fXaxis.GetXmin() || fYaxis.GetXmax() <= fYaxis.GetXmin()) {
      //find min, max of entries in buffer
      Double_t xmin = fBuffer[2];
      Double_t xmax = xmin;
      Double_t ymin = fBuffer[3];
      Double_t ymax = ymin;
      for (Int_t i=1;i<nbentries;i++) {
         Double_t x = fBuffer[4*i+2];
         if (x < xmin) xmin = x;
         if (x > xmax) xmax = x;
         Double_t y = fBuffer[4*i+3];
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
/// Accumulate arguments in buffer.
///
/// When buffer is full, empty the buffer.
///
///  - fBuffer[0] = number of entries in buffer
///  - fBuffer[1] = w of first entry
///  - fBuffer[2] = x of first entry
///  - fBuffer[3] = y of first entry
///  - fBuffer[4] = z of first entry

Int_t TProfile2D::BufferFill(Double_t x, Double_t y, Double_t z, Double_t w)
{
   if (!fBuffer) return -3;
   Int_t nbentries = (Int_t)fBuffer[0];
   if (nbentries < 0) {
      nbentries  = -nbentries;
      fBuffer[0] =  nbentries;
      if (fEntries > 0) {
         Double_t *buffer = fBuffer; fBuffer=0;
         Reset("ICES"); // reset without deleting the functions
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
   return -2;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a Profile2D histogram to a new profile2D histogram.

void TProfile2D::Copy(TObject &obj) const
{
   try {
      TProfile2D & pobj = dynamic_cast<TProfile2D&>(obj);

      TH2D::Copy(pobj);
      fBinEntries.Copy(pobj.fBinEntries);
      fBinSumw2.Copy(pobj.fBinSumw2);
      for (int bin=0;bin<fNcells;bin++) {
         pobj.fArray[bin]        = fArray[bin];
         pobj.fSumw2.fArray[bin] = fSumw2.fArray[bin];
      }
      pobj.fZmin = fZmin;
      pobj.fZmax = fZmax;
      pobj.fScaling   = fScaling;
      pobj.fErrorMode = fErrorMode;
      pobj.fTsumwz    = fTsumwz;
      pobj.fTsumwz2   = fTsumwz2;

   } catch(...) {
      Fatal("Copy","Cannot copy a TProfile2D in a %s",obj.IsA()->GetName());
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: `this = this/(c1*f1)` .
/// This function is not implemented

Bool_t TProfile2D::Divide(TF1 *, Double_t )
{
   Error("Divide","Function not implemented for TProfile2D");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this profile2D by h1.
///
/// `this = this/h1`
///
///This function return kFALSE if the divide operation failed

Bool_t TProfile2D::Divide(const TH1 *h1)
{

   if (!h1) {
      Error("Divide","Attempt to divide a non-existing profile2D");
      return kFALSE;
   }
   if (!h1->InheritsFrom(TProfile2D::Class())) {
      Error("Divide","Attempt to divide a non-profile2D object");
      return kFALSE;
   }
   TProfile2D *p1 = (TProfile2D*)h1;

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);

   // Check profile compatibility
   Int_t nx = GetNbinsX();
   if (nx != p1->GetNbinsX()) {
      Error("Divide","Attempt to divide profiles with different number of bins");
      return kFALSE;
   }
   Int_t ny = GetNbinsY();
   if (ny != p1->GetNbinsY()) {
      Error("Divide","Attempt to divide profiles with different number of bins");
      return kFALSE;
   }

   // Reset statistics
   fEntries = fTsumw   = fTsumw2 = fTsumwx = fTsumwx2 = 0;

   // Loop on bins (including underflows/overflows)
   Int_t bin,binx,biny;
   Double_t *cu1 = p1->GetW();
   Double_t *er1 = p1->GetW2();
   Double_t *en1 = p1->GetB();
   Double_t c0,c1,w,z,x,y;
   for (binx =0;binx<=nx+1;binx++) {
      for (biny =0;biny<=ny+1;biny++) {
         bin   = biny*(fXaxis.GetNbins()+2) + binx;
         c0  = fArray[bin];
         c1  = cu1[bin];
         if (c1) w = c0/c1;
         else    w = 0;
         fArray[bin] = w;
         z = TMath::Abs(w);
         x = fXaxis.GetBinCenter(binx);
         y = fYaxis.GetBinCenter(biny);
         fEntries++;
         fTsumw   += z;
         fTsumw2  += z*z;
         fTsumwx  += z*x;
         fTsumwx2 += z*x*x;
         fTsumwy  += z*y;
         fTsumwy2 += z*y*y;
         fTsumwxy += z*x*y;
         fTsumwz  += z;
         fTsumwz2 += z*z;
         Double_t e0 = fSumw2.fArray[bin];
         Double_t e1 = er1[bin];
         Double_t c12= c1*c1;
         if (!c1) fSumw2.fArray[bin] = 0;
         else     fSumw2.fArray[bin] = (e0*c1*c1 + e1*c0*c0)/(c12*c12);
         if (!en1[bin]) fBinEntries.fArray[bin] = 0;
         else           fBinEntries.fArray[bin] /= en1[bin];
      }
   }
   // maintaining the correct sum of weights square is not supported when dividing
   // bin error resulting from division of profile needs to be checked
   if (fBinSumw2.fN) {
      Warning("Divide","Cannot preserve during the division of profiles the sum of bin weight square");
      fBinSumw2 = TArrayD();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace contents of this profile2D by the division of h1 by h2.
///
/// `this = c1*h1/(c2*h2)`
///
///   This function return kFALSE if the divide operation failed

Bool_t TProfile2D::Divide(const TH1 *h1, const TH1 *h2, Double_t c1, Double_t c2, Option_t *option)
{
   TString opt = option;
   opt.ToLower();
   Bool_t binomial = kFALSE;
   if (opt.Contains("b")) binomial = kTRUE;
   if (!h1 || !h2) {
      Error("Divide","Attempt to divide a non-existing profile2D");
      return kFALSE;
   }
   if (!h1->InheritsFrom(TProfile2D::Class())) {
      Error("Divide","Attempt to divide a non-profile2D object");
      return kFALSE;
   }
   TProfile2D *p1 = (TProfile2D*)h1;
   if (!h2->InheritsFrom(TProfile2D::Class())) {
      Error("Divide","Attempt to divide a non-profile2D object");
      return kFALSE;
   }
   TProfile2D *p2 = (TProfile2D*)h2;

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);

   // Check histogram compatibility
   Int_t nx = GetNbinsX();
   if (nx != p1->GetNbinsX() || nx != p2->GetNbinsX()) {
      Error("Divide","Attempt to divide profiles with different number of bins");
      return kFALSE;
   }
   Int_t ny = GetNbinsY();
   if (ny != p1->GetNbinsY() || ny != p2->GetNbinsY()) {
      Error("Divide","Attempt to divide profiles with different number of bins");
      return kFALSE;
   }
   if (!c2) {
      Error("Divide","Coefficient of dividing profile cannot be zero");
      return kFALSE;
   }

   // Reset statistics
   fEntries = fTsumw   = fTsumw2 = fTsumwx = fTsumwx2 = 0;

   // Loop on bins (including underflows/overflows)
   Int_t bin,binx,biny;
   Double_t *cu1 = p1->GetW();
   Double_t *cu2 = p2->GetW();
   Double_t *er1 = p1->GetW2();
   Double_t *er2 = p2->GetW2();
   Double_t *en1 = p1->GetB();
   Double_t *en2 = p2->GetB();
   Double_t b1,b2,w,z,x,y,ac1,ac2;
   ac1 = TMath::Abs(c1);
   ac2 = TMath::Abs(c2);
   for (binx =0;binx<=nx+1;binx++) {
      for (biny =0;biny<=ny+1;biny++) {
         bin   = biny*(fXaxis.GetNbins()+2) + binx;
         b1  = cu1[bin];
         b2  = cu2[bin];
         if (b2) w = c1*b1/(c2*b2);
         else    w = 0;
         fArray[bin] = w;
         z = TMath::Abs(w);
         x = fXaxis.GetBinCenter(binx);
         y = fYaxis.GetBinCenter(biny);
         fEntries++;
         fTsumw   += z;
         fTsumw2  += z*z;
         fTsumwx  += z*x;
         fTsumwx2 += z*x*x;
         fTsumwy  += z*y;
         fTsumwy2 += z*y*y;
         fTsumwxy += z*x*y;
         fTsumwz  += z;
         fTsumwz2 += z*z;
         Double_t e1 = er1[bin];
         Double_t e2 = er2[bin];
         //Double_t b22= b2*b2*d2;
         Double_t b22= b2*b2*TMath::Abs(c2);
         if (!b2) fSumw2.fArray[bin] = 0;
         else {
            if (binomial) {
               fSumw2.fArray[bin] = TMath::Abs(w*(1-w)/(c2*b2));
            } else {
               fSumw2.fArray[bin] = ac1*ac2*(e1*b2*b2 + e2*b1*b1)/(b22*b22);
            }
         }
         if (!en2[bin]) fBinEntries.fArray[bin] = 0;
         else           fBinEntries.fArray[bin] = en1[bin]/en2[bin];
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a Profile2D histogram (no weights).

Int_t TProfile2D::Fill(Double_t x, Double_t y, Double_t z)
{
   if (fBuffer) return BufferFill(x,y,z,1);

   Int_t bin,binx,biny;

   if (fZmin != fZmax) {
      if (z <fZmin || z> fZmax || TMath::IsNaN(z) ) return -1;
   }

   fEntries++;
   binx =fXaxis.FindBin(x);
   biny =fYaxis.FindBin(y);
   if (binx <0 || biny <0) return -1;
   bin = GetBin(binx, biny);
   fArray[bin] += z;
   fSumw2.fArray[bin] += z*z;
   fBinEntries.fArray[bin] += 1;
   if (fBinSumw2.fN)  fBinSumw2.fArray[bin] += 1;
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
   fTsumwz  += z;
   fTsumwz2 += z*z;
   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a Profile2D histogram with weights.

Int_t TProfile2D::Fill(Double_t x, Double_t y, Double_t z, Double_t w)
{
   if (fBuffer)  return BufferFill(x, y, z, w);

   Int_t bin, binx, biny;

   if (fZmin != fZmax) {
      if (z < fZmin || z > fZmax || TMath::IsNaN(z)) return -1;
   }

   Double_t u = w;
   fEntries++;
   binx = fXaxis.FindBin(x);
   biny = fYaxis.FindBin(y);
   if (binx < 0 || biny < 0) return -1;
   bin = biny * (fXaxis.GetNbins() + 2) + binx;
   AddBinContent(bin, u * z);
   fSumw2.fArray[bin] += u * z * z;
   if (!fBinSumw2.fN && u != 1.0 && !TestBit(TH1::kIsNotW))
      Sumw2(); // must be called before accumulating the entries
   if (fBinSumw2.fN) fBinSumw2.fArray[bin] += u * u;
   fBinEntries.fArray[bin] += u;

   if (binx == 0 || binx > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (biny == 0 || biny > fYaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   fTsumw += u;
   fTsumw2 += u * u;
   fTsumwx += u * x;
   fTsumwx2 += u * x * x;
   fTsumwy += u * y;
   fTsumwy2 += u * y * y;
   fTsumwxy += u * x * y;
   fTsumwz += u * z;
   fTsumwz2 += u * z * z;
   return bin;
}
////////////////////////////////////////////////////////////////////////////////
/// Fill a Profile2D histogram (no weights).

Int_t TProfile2D::Fill(Double_t x, const char *namey, Double_t z, Double_t w)
{
   Int_t bin,binx,biny;

   if (fZmin != fZmax) {
      if (z <fZmin || z> fZmax || TMath::IsNaN(z)) return -1;
   }

   Double_t u = w;
   fEntries++;
   binx =fXaxis.FindBin(x);
   biny =fYaxis.FindBin(namey);
   if (binx <0 || biny <0) return -1;
   bin  = biny*(fXaxis.GetNbins()+2) + binx;
   AddBinContent(bin, u * z);
   fSumw2.fArray[bin] += u * z * z;
   if (!fBinSumw2.fN && u != 1.0 && !TestBit(TH1::kIsNotW))
      Sumw2(); // must be called before accumulating the entries
   if (fBinSumw2.fN) fBinSumw2.fArray[bin] += u * u;
   fBinEntries.fArray[bin] += u;

   if (binx == 0 || binx > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (biny == 0 || biny > fYaxis.GetNbins()) return -1;

   UInt_t labelBitMask = GetAxisLabelStatus();
   Double_t y = (labelBitMask & TH1::kYaxis) ? 0 : fYaxis.GetBinCenter(biny);

   fTsumw += u;
   fTsumw2 += u * u;
   fTsumwx += u * x;
   fTsumwx2 += u * x * x;
   fTsumwy += u * y;
   fTsumwy2 += u * y * y;
   fTsumwxy += u * x * y;
   fTsumwz += u * z;
   fTsumwz2 += u * z * z;

   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a Profile2D histogram (no weights).

Int_t TProfile2D::Fill(const char *namex, const char *namey, Double_t z, Double_t w)
{
   Int_t bin,binx,biny;

   if (fZmin != fZmax) {
      if (z <fZmin || z> fZmax || TMath::IsNaN(z) ) return -1;
   }

   Double_t u = w;
   fEntries++;
   binx = fXaxis.FindBin(namex);
   biny = fYaxis.FindBin(namey);
   if (binx < 0 || biny < 0)
      return -1;
   bin = biny * (fXaxis.GetNbins() + 2) + binx;
   AddBinContent(bin, u * z);
   fSumw2.fArray[bin] += u * z * z;
   if (!fBinSumw2.fN && u != 1.0 && !TestBit(TH1::kIsNotW))
      Sumw2(); // must be called before accumulating the entries
   if (fBinSumw2.fN) fBinSumw2.fArray[bin] += u * u;
   fBinEntries.fArray[bin] += u;

   if (binx == 0 || binx > fXaxis.GetNbins()) return -1;
   if (biny == 0 || biny > fYaxis.GetNbins()) return -1;

   UInt_t labelBitMask = GetAxisLabelStatus();
   Double_t x = (labelBitMask & TH1::kXaxis) ? 0 : fXaxis.GetBinCenter(binx);
   Double_t y = (labelBitMask & TH1::kYaxis) ? 0 : fYaxis.GetBinCenter(biny);

   fTsumw += u;
   fTsumw2 += u * u;
   fTsumwx += u * x;
   fTsumwx2 += u * x * x;
   fTsumwy += u * y;
   fTsumwy2 += u * y * y;
   fTsumwxy += u * x * y;
   fTsumwz += u * z;
   fTsumwz2 += u * z * z;

   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a Profile2D histogram (no weights).

Int_t TProfile2D::Fill(const char *namex, Double_t y, Double_t z, Double_t w)
{
   Int_t bin,binx,biny;

   if (fZmin != fZmax) {
      if (z <fZmin || z> fZmax || TMath::IsNaN(z)) return -1;
   }

   Double_t u = w;
   fEntries++;
   binx =fXaxis.FindBin(namex);
   biny =fYaxis.FindBin(y);
   if (binx <0 || biny <0) return -1;
   bin  = biny*(fXaxis.GetNbins()+2) + binx;

   AddBinContent(bin, u * z);
   fSumw2.fArray[bin] += u * z * z;
   if (!fBinSumw2.fN && u != 1.0 && !TestBit(TH1::kIsNotW))
      Sumw2(); // must be called before accumulating the entries
   if (fBinSumw2.fN) fBinSumw2.fArray[bin] += u * u;
   fBinEntries.fArray[bin] += u;

   if (binx == 0 || binx > fXaxis.GetNbins()) return -1;
   if (biny == 0 || biny > fYaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }

   UInt_t labelBitMask = GetAxisLabelStatus();
   Double_t x = (labelBitMask & TH1::kXaxis) ? 0 : fXaxis.GetBinCenter(binx);

   fTsumw += u;
   fTsumw2 += u * u;
   fTsumwx += u * x;
   fTsumwx2 += u * x * x;
   fTsumwy += u * y;
   fTsumwy2 += u * y * y;
   fTsumwxy += u * x * y;
   fTsumwz += u * z;
   fTsumwz2 += u * z * z;

   return bin;
}



////////////////////////////////////////////////////////////////////////////////
/// Return bin content of a Profile2D histogram.

Double_t TProfile2D::GetBinContent(Int_t bin) const
{
   if (fBuffer) ((TProfile2D*)this)->BufferEmpty();

   if (bin < 0 || bin >= fNcells) return 0;
   if (fBinEntries.fArray[bin] == 0) return 0;
   if (!fArray) return 0;
   return fArray[bin]/fBinEntries.fArray[bin];
}

////////////////////////////////////////////////////////////////////////////////
/// Return bin entries of a Profile2D histogram.

Double_t TProfile2D::GetBinEntries(Int_t bin) const
{
   if (fBuffer) ((TProfile2D*)this)->BufferEmpty();

   if (bin < 0 || bin >= fNcells) return 0;
   return fBinEntries.fArray[bin];
}

////////////////////////////////////////////////////////////////////////////////
/// Return bin effective entries for a weighted filled Profile histogram.
/// In case of an unweighted profile, it is equivalent to the number of entries per bin
/// The effective entries is defined as the square of the sum of the weights divided by the
/// sum of the weights square.
/// TProfile::Sumw2() must be called before filling the profile with weights.
/// Only by calling this method the  sum of the square of the weights per bin is stored.

Double_t TProfile2D::GetBinEffectiveEntries(Int_t bin)
{
   return TProfileHelper::GetBinEffectiveEntries(this, bin);
}

////////////////////////////////////////////////////////////////////////////////
/// Return bin error of a Profile2D histogram.
///
/// ### Computing errors: A moving field
///
/// The computation of errors for a TProfile2D has evolved with the versions
/// of ROOT. The difficulty is in computing errors for bins with low statistics.
/// - prior to version 3.10, we had no special treatment of low statistic bins.
///   As a result, these bins had huge errors. The reason is that the
///   expression eprim2 is very close to 0 (rounding problems) or 0.
/// - The algorithm is modified/protected for the case
///   when a TProfile2D is projected (ProjectionX). The previous algorithm
///   generated a N^2 problem when projecting a TProfile2D with a large number of
///   bins (eg 100000).
/// - in version 3.10/02, a new static function TProfile::Approximate
///   is introduced to enable or disable (default) the approximation.
///   (see also comments in TProfile::GetBinError)

Double_t TProfile2D::GetBinError(Int_t bin) const
{
   return TProfileHelper::GetBinError((TProfile2D*)this, bin);
}

////////////////////////////////////////////////////////////////////////////////
/// Return option to compute profile2D errors.

Option_t *TProfile2D::GetErrorOption() const
{
   if (fErrorMode == kERRORSPREAD)  return "s";
   if (fErrorMode == kERRORSPREADI) return "i";
   if (fErrorMode == kERRORSPREADG) return "g";
   return "";
}

////////////////////////////////////////////////////////////////////////////////
/// Fill the array stats from the contents of this profile.
/// The array stats must be correctly dimensioned in the calling program.
///
/// - stats[0] = sumw
/// - stats[1] = sumw2
/// - stats[2] = sumwx
/// - stats[3] = sumwx2
/// - stats[4] = sumwy
/// - stats[5] = sumwy2
/// - stats[6] = sumwxy
/// - stats[7] = sumwz
/// - stats[8] = sumwz2
///
/// If no axis-subrange is specified (via TAxis::SetRange), the array stats
/// is simply a copy of the statistics quantities computed at filling time.
/// If a sub-range is specified, the function recomputes these quantities
/// from the bin contents in the current axis range.

void TProfile2D::GetStats(Double_t *stats) const
{
   if (fBuffer) ((TProfile2D*)this)->BufferEmpty();

   // check for labels axis . In that case corresponsing statistics do not make sense and it is set to zero
   Bool_t labelXaxis =  ((const_cast<TAxis&>(fXaxis)).GetLabels() && fXaxis.CanExtend() );
   Bool_t labelYaxis =  ((const_cast<TAxis&>(fYaxis)).GetLabels() && fYaxis.CanExtend() );

   // Loop on bins
   if ( (fTsumw == 0 /* && fEntries > 0 */) || fXaxis.TestBit(TAxis::kAxisRange) || fYaxis.TestBit(TAxis::kAxisRange)) {
      Int_t bin, binx, biny;
      Double_t w, w2;
      Double_t x,y;
      for (bin=0;bin<9;bin++) stats[bin] = 0;
      if (!fBinEntries.fArray) return;
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
      for (biny = firstBinY; biny <= lastBinY; biny++) {
         y = (!labelYaxis) ? fYaxis.GetBinCenter(biny) : 0;
         for (binx = firstBinX; binx <= lastBinX; binx++) {
            bin = GetBin(binx,biny);
            w         = fBinEntries.fArray[bin];
            w2        = (fBinSumw2.fN ? fBinSumw2.fArray[bin] : w );
            x         = (!labelXaxis) ? fXaxis.GetBinCenter(binx) : 0;
            stats[0] += w;
            stats[1] += w2;
            stats[2] += w*x;
            stats[3] += w*x*x;
            stats[4] += w*y;
            stats[5] += w*y*y;
            stats[6] += w*x*y;
            stats[7] += fArray[bin];
            stats[8] += fSumw2.fArray[bin];
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
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reduce the number of bins for this axis to the number of bins having a label.

void TProfile2D::LabelsDeflate(Option_t *ax)
{
   TProfileHelper::LabelsDeflate(this, ax);
}

////////////////////////////////////////////////////////////////////////////////
/// Double the number of bins for axis.
/// Refill histogram
/// This function is called by TAxis::FindBin(const char *label)

void TProfile2D::LabelsInflate(Option_t *ax)
{
   TProfileHelper::LabelsInflate(this, ax);
}

////////////////////////////////////////////////////////////////////////////////
/// Set option(s) to draw axis with labels.
///
/// option might have the following values:
///
///  - "a" sort by alphabetic order
///  - ">" sort by decreasing values
///  - "<" sort by increasing values
///  - "h" draw labels horizontal
///  - "v" draw labels vertical
///  - "u" draw labels up (end of label right adjusted)
///  - "d" draw labels down (start of label left adjusted)

void TProfile2D::LabelsOption(Option_t *option, Option_t *ax)
{

   TAxis *axis = GetXaxis();
   if (ax[0] == 'y' || ax[0] == 'Y') axis = GetYaxis();
   THashList *labels = axis->GetLabels();
   if (!labels) {
      Warning("LabelsOption","Cannot sort. No labels");
      return;
   }
   TString opt = option;
   opt.ToLower();
   if (opt.Contains("h")) {
      axis->SetBit(TAxis::kLabelsHori);
      axis->ResetBit(TAxis::kLabelsVert);
      axis->ResetBit(TAxis::kLabelsDown);
      axis->ResetBit(TAxis::kLabelsUp);
   }
   if (opt.Contains("v")) {
      axis->SetBit(TAxis::kLabelsVert);
      axis->ResetBit(TAxis::kLabelsHori);
      axis->ResetBit(TAxis::kLabelsDown);
      axis->ResetBit(TAxis::kLabelsUp);
   }
   if (opt.Contains("u")) {
      axis->SetBit(TAxis::kLabelsUp);
      axis->ResetBit(TAxis::kLabelsVert);
      axis->ResetBit(TAxis::kLabelsDown);
      axis->ResetBit(TAxis::kLabelsHori);
   }
   if (opt.Contains("d")) {
      axis->SetBit(TAxis::kLabelsDown);
      axis->ResetBit(TAxis::kLabelsVert);
      axis->ResetBit(TAxis::kLabelsHori);
      axis->ResetBit(TAxis::kLabelsUp);
   }
   Int_t sort = -1;
   if (opt.Contains("a")) sort = 0;
   if (opt.Contains(">")) sort = 1;
   if (opt.Contains("<")) sort = 2;
   if (sort < 0) return;

   // support only cases where each bin has a labels (should be when axis is alphanumeric)
   Int_t n = labels->GetSize();
   if (n != axis->GetNbins()) {
      // check if labels are all consecutive and starts from the first bin
      // in that case the current code will work fine
      Int_t firstLabelBin = axis->GetNbins() + 1;
      Int_t lastLabelBin = -1;
      for (Int_t i = 0; i < n; ++i) {
         Int_t bin = labels->At(i)->GetUniqueID();
         if (bin < firstLabelBin)
            firstLabelBin = bin;
         if (bin > lastLabelBin)
            lastLabelBin = bin;
      }
      if (firstLabelBin != 1 || lastLabelBin - firstLabelBin + 1 != n) {
         Error("LabelsOption",
               "%s of TProfile2D %s contains bins without labels. Sorting will not work correctly - return",
               axis->GetName(), GetName());
         return;
      }
      // case where label bins are consecutive starting from first bin will work
      Warning(
         "LabelsOption",
         "axis %s of TProfile2D %s has extra following bins without labels. Sorting will work only for first label bins",
         axis->GetName(), GetName());
   }

   std::vector<Int_t> a(n);
   Int_t i, j, k, ibin, bin;
   std::vector<Double_t> sumw(fNcells);
   std::vector<Double_t> errors(fNcells);
   std::vector<Double_t> ent(fNcells);
   std::vector<Double_t> binsw2;
   if (fBinSumw2.fN)
      binsw2.resize(fNcells);

   // delete buffer if it is there since bins will be reordered.
   if (fBuffer)
      BufferEmpty(1);

   // number of bins to loop
   Int_t nx = (axis == GetXaxis()) ? n + 1 : fXaxis.GetNbins() + 2;
   Int_t ny = (axis == GetYaxis()) ? n + 1 : fYaxis.GetNbins() + 2;

   // make a labelold list but ordered with bins
   // (re-ordered original label list)
   std::vector<TObject *> labold(n);
   for (i = 0; i < n; i++)
      labold[i] = nullptr;
   TIter nextold(labels);
   TObject *obj;
   while ((obj = nextold())) {
      bin = obj->GetUniqueID();
      if (bin>=1 && bin<=n)
        labold[bin-1] = obj;
   }
   // order now labold according to bin content

   labels->Clear();

   std::vector<Double_t> pcont;
   std::vector<Double_t> econt;
   if (sort > 0) {
      pcont.resize(n);
      econt.resize(n);
   }

   for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
         bin = GetBin(i, j);
         sumw[bin] = fArray[bin];
         errors[bin] = fSumw2.fArray[bin];
         ent[bin] = fBinEntries.fArray[bin];
         if (fBinSumw2.fN)
            binsw2[bin] = fBinSumw2.fArray[bin];
         if (axis == GetXaxis())
            k = i - 1;
         else
            k = j - 1;
         //---when sorting by values of bins
         if (sort > 0 && fBinEntries.fArray[bin] != 0 && k > 0 && k < n) {
            pcont[k] += fArray[bin];
            econt[k] += fBinEntries.fArray[bin];
         }
      }
   }
   // compute average of slize for ordering
   if (sort > 0) {
      for (k = 0; k < n; ++k) {
         a[k] = k;
         if (econt[k] > 0)
            pcont[k] /= econt[k];
      }
      if (sort == 1)
         TMath::Sort(n, pcont.data(), a.data(), kTRUE); // sort by decreasing values
      else
         TMath::Sort(n, pcont.data(), a.data(), kFALSE); // sort by increasing values
   } else {
      //---alphabetic sort
      // sort labels using vector of strings and TMath::Sort
      // I need to array because labels order in list is not necessary that of the bins
      std::vector<std::string> vecLabels(n);
      for (i = 0; i < n; i++) {
         vecLabels[i] = labold[i]->GetName();
         a[i] = i;
      }
      // sort in ascending order for strings
      TMath::Sort(n, vecLabels.data(), a.data(), kFALSE);
   }

   // set the new labels
   for (i = 0; i < n; i++) {
      obj = labold[a[i]];
      labels->Add(obj);
      // set the corresponding bin. NB bin starts from 1
      obj->SetUniqueID(i + 1);
      if (gDebug)
         std::cout << "bin " << i + 1 << " setting new labels for axis " << labold.at(a[i])->GetName() << " from "
                   << a[i] << std::endl;
   }
   // set the new content
   for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
         bin = GetBin(i, j);
         if (axis == GetXaxis()) {
            if (i == 0) break; // skip underflow in x
            ibin = GetBin(a[i - 1] + 1, j);
         } else {
            if (j == 0) continue;  // skip underflow in y
            ibin = GetBin(i, a[j-1] + 1);
         }
         fArray[bin] = sumw[ibin];
         fSumw2.fArray[bin] = errors[ibin];
         fBinEntries.fArray[bin] = ent[ibin];
         if (fBinSumw2.fN)
            fBinSumw2.fArray[bin] = binsw2[ibin];
      }
   }
   // need to set to zero the statistics if axis has been sorted
   // see for example TH3::PutStats for definition of s vector
   bool labelsAreSorted = kFALSE;
   for (i = 0; i < n; ++i) {
      if (a[i] != i) {
         labelsAreSorted = kTRUE;
         break;
      }
   }
   if (labelsAreSorted) {
      double s[TH1::kNstat];
      GetStats(s);
      if (axis == GetXaxis()) {
         s[2] = 0; // fTsumwx
         s[3] = 0; // fTsumwx2
         s[6] = 0; // fTsumwxy
      } else  {
         s[4] = 0;  // fTsumwy
         s[5] = 0;  // fTsumwy2
         s[6] = 0;  // fTsumwxy
      }
      PutStats(s);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Merge all histograms in the collection in this histogram.
/// This function computes the min/max for the axes,
/// compute a new number of bins, if necessary,
/// add bin contents, errors and statistics.
/// If overflows are present and limits are different the function will fail.
/// The function returns the total number of entries in the result histogram
/// if the merge is successful, -1 otherwise.
///
/// IMPORTANT remark. The 2 axis x and y may have different number
/// of bins and different limits, BUT the largest bin width must be
/// a multiple of the smallest bin width and the upper limit must also
/// be a multiple of the bin width.

Long64_t TProfile2D::Merge(TCollection *li)
{
   return TProfileHelper::Merge(this, li);
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: this = this*c1*f1

Bool_t TProfile2D::Multiply(TF1 *, Double_t )
{
   Error("Multiply","Function not implemented for TProfile2D");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply this profile2D by h1.
///
///   `this = this*h1`

Bool_t TProfile2D::Multiply(const TH1 *)
{
   Error("Multiply","Multiplication of profile2D histograms not implemented");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace contents of this profile2D by multiplication of h1 by h2.
///
/// `this = (c1*h1)*(c2*h2)`

Bool_t TProfile2D::Multiply(const TH1 *, const TH1 *, Double_t, Double_t, Option_t *)
{
   Error("Multiply","Multiplication of profile2D histograms not implemented");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Project this profile2D into a 2-D histogram along X,Y.
///
/// The projection is always of the type TH2D.
///
///  - if option "E" is specified  the errors of the projected histogram are computed and set
///    to be equal to the errors of the profile.
///    Option "E" is defined as the default one in the header file.
///  - if option "" is specified the histogram errors are simply the sqrt of its content
///  - if option "B" is specified, the content of bin of the returned histogram
///    will be equal to the GetBinEntries(bin) of the profile,
///  - if option "C=E" the bin contents of the projection are set to the
///    bin errors of the profile
///  - if option "W" is specified the bin content of the projected histogram  is set to the
///    product of the bin content of the profile and the entries.
///    With this option the returned histogram will be equivalent to the one obtained by
///    filling directly a TH2D using the 3-rd value as a weight.
///    This option makes sense only for profile filled with all weights =1.
///    When the profile is weighted (filled with weights different than 1) the
///    bin error of the projected histogram (obtained using this option "W") cannot be
///    correctly computed from the information stored in the profile. In that case the
///    obtained histogram contains as bin error square the weighted sum of the square of the
///    profiled observable (TProfile2D::fSumw2[bin] )

TH2D *TProfile2D::ProjectionXY(const char *name, Option_t *option) const
{

   TString opt = option;
   opt.ToLower();

   // Create the projection histogram
   // name of projected histogram is by default name of original histogram + _pxy
   TString pname(name);
   if (pname.IsNull() || pname == "_pxy")
      pname = TString(GetName() ) + TString("_pxy");


   Int_t nx = fXaxis.GetNbins();
   Int_t ny = fYaxis.GetNbins();
   const TArrayD *xbins = fXaxis.GetXbins();
   const TArrayD *ybins = fYaxis.GetXbins();
   TH2D * h1 = 0;
   if (xbins->fN == 0 && ybins->fN == 0) {
      h1 = new TH2D(pname,GetTitle(),nx,fXaxis.GetXmin(),fXaxis.GetXmax(),ny,fYaxis.GetXmin(),fYaxis.GetXmax());
   } else if (xbins->fN == 0) {
      h1 = new TH2D(pname,GetTitle(),nx,fXaxis.GetXmin(),fXaxis.GetXmax(),ny, ybins->GetArray() );
   } else if (ybins->fN == 0) {
      h1 = new TH2D(pname,GetTitle(),nx,xbins->GetArray(),ny,fYaxis.GetXmin(),fYaxis.GetXmax());
   } else {
      h1 = new TH2D(pname,GetTitle(),nx,xbins->GetArray(),ny,ybins->GetArray() );
   }
   Bool_t computeErrors = kFALSE;
   Bool_t cequalErrors  = kFALSE;
   Bool_t binEntries    = kFALSE;
   Bool_t binWeight     = kFALSE;
   if (opt.Contains("b")) binEntries = kTRUE;
   if (opt.Contains("e")) computeErrors = kTRUE;
   if (opt.Contains("w")) binWeight = kTRUE;
   if (opt.Contains("c=e")) {cequalErrors = kTRUE; computeErrors=kFALSE;}
   if (computeErrors || binWeight || (binEntries && fBinSumw2.fN)  ) h1->Sumw2();

   // Fill the projected histogram
   Int_t bin,binx, biny;
   Double_t cont;
   for (binx =0;binx<=nx+1;binx++) {
      for (biny =0;biny<=ny+1;biny++) {
         bin = GetBin(binx,biny);

         if (binEntries)         cont = GetBinEntries(bin);
         else if (cequalErrors)  cont = GetBinError(bin);
         else if (binWeight)     cont = GetBinContent(bin) * GetBinEntries(bin);
         else                    cont = GetBinContent(bin);    // default case

         h1->SetBinContent(bin ,cont);

         // if option E projected histogram errors are same as profile
         if (computeErrors ) h1->SetBinError(bin , GetBinError(bin) );
         // in case of option W bin error is deduced from bin sum of z**2 values of profile
         // this is correct only if the profile is unweighted
         if (binWeight)      h1->GetSumw2()->fArray[bin] = fSumw2.fArray[bin];
         // in case of bin entries and profile is weighted, we need to set also the bin error
         if (binEntries && fBinSumw2.fN ) {
            R__ASSERT(  h1->GetSumw2() );
            h1->GetSumw2()->fArray[bin] =  fBinSumw2.fArray[bin];
         }
      }
   }
   h1->SetEntries(fEntries);
   return h1;
}

////////////////////////////////////////////////////////////////////////////////
/// Project a 2-D histogram into a profile histogram along X.
///
/// The projection is made from the channels along the Y axis
/// ranging from firstybin to lastybin included.
/// The result is a 1D profile which contains the combination of all the considered bins along Y
/// By default, bins 1 to ny are included
/// When all bins are included, the number of entries in the projection
/// is set to the number of entries of the 2-D histogram, otherwise
/// the number of entries is incremented by 1 for all non empty cells.
///
/// The option can also be used to specify the projected profile error type.
/// Values which can be used are 's', 'i', or 'g'. See TProfile::BuildOptions for details

TProfile *TProfile2D::ProfileX(const char *name, Int_t firstybin, Int_t lastybin, Option_t *option) const
{
   return DoProfile(true, name, firstybin, lastybin, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Project a 2-D histogram into a profile histogram along X
///
/// The projection is made from the channels along the X axis
/// ranging from firstybin to lastybin included.
/// The result is a 1D profile which contains the combination of all the considered bins along X
/// By default, bins 1 to ny are included
/// When all bins are included, the number of entries in the projection
/// is set to the number of entries of the 2-D histogram, otherwise
/// the number of entries is incremented by 1 for all non empty cells.
///
/// The option can also be used to specify the projected profile error type.
/// Values which can be used are 's', 'i', or 'g'. See TProfile::BuildOptions for details

TProfile *TProfile2D::ProfileY(const char *name, Int_t firstxbin, Int_t lastxbin, Option_t *option) const
{
   return DoProfile(false, name, firstxbin, lastxbin, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Implementation of ProfileX or ProfileY for a TProfile2D.
///
/// Do correctly the combination of the bin averages when doing the projection

TProfile * TProfile2D::DoProfile(bool onX, const char *name, Int_t firstbin, Int_t lastbin, Option_t *option) const {
   TString opt = option;
   opt.ToLower();
   bool originalRange = opt.Contains("o");

   TString expectedName = ( onX ? "_pfx" : "_pfy" );

   TString pname(name);
   if (pname.IsNull() || name == expectedName)
      pname = TString(GetName() ) + expectedName;

   const TAxis& outAxis = ( onX ? fXaxis : fYaxis );
   const TArrayD *bins = outAxis.GetXbins();
   Int_t firstOutBin = outAxis.GetFirst();
   Int_t lastOutBin = outAxis.GetLast();

   TProfile  * p1 = 0;
   // case of fixed bins
   if (bins->fN == 0) {
      if (originalRange)
         p1 =  new TProfile(pname,GetTitle(), outAxis.GetNbins(), outAxis.GetXmin(), outAxis.GetXmax(), opt );
      else
         p1 =  new TProfile(pname,GetTitle(), lastOutBin-firstOutBin+1,
                            outAxis.GetBinLowEdge(firstOutBin),outAxis.GetBinUpEdge(lastOutBin), opt);
   } else {
      // case of variable bins
      if (originalRange )
         p1 = new TProfile(pname,GetTitle(),outAxis.GetNbins(),bins->fArray,opt);
      else
         p1 = new TProfile(pname,GetTitle(),lastOutBin-firstOutBin+1,&bins->fArray[firstOutBin-1],opt);

   }

   if (fBinSumw2.fN) p1->Sumw2();

   // make projection in a 2D first
   TH2D * h2dW = ProjectionXY("h2temp-W","W");
   TH2D * h2dN = ProjectionXY("h2temp-N","B");

   h2dW->SetDirectory(0); h2dN->SetDirectory(0);


   TString opt1 = (originalRange) ? "o" : "";
   TH1D * h1W = (onX) ? h2dW->ProjectionX("h1temp-W",firstbin,lastbin,opt1) : h2dW->ProjectionY("h1temp-W",firstbin,lastbin,opt1);
   TH1D * h1N = (onX) ? h2dN->ProjectionX("h1temp-N",firstbin,lastbin,opt1) : h2dN->ProjectionY("h1temp-N",firstbin,lastbin,opt1);
   h1W->SetDirectory(0); h1N->SetDirectory(0);


   // fill the bin content
   R__ASSERT( h1W->fN == p1->fN );
   R__ASSERT( h1N->fN == p1->fN );
   R__ASSERT( h1W->GetSumw2()->fN != 0); // h1W should always be a weighted histogram since h2dW is
   for (int i = 0; i < p1->fN ; ++i) {
      p1->fArray[i] = h1W->GetBinContent(i);   // array of profile is sum of all values
      p1->GetSumw2()->fArray[i]  = h1W->GetSumw2()->fArray[i];   // array of content square of profile is weight square of the W projected histogram
      p1->SetBinEntries(i, h1N->GetBinContent(i) );
      if (fBinSumw2.fN) p1->GetBinSumw2()->fArray[i] = h1N->GetSumw2()->fArray[i];    // sum of weight squares are stored to compute errors in h1N histogram
   }
   // delete the created histograms
   delete h2dW;
   delete h2dN;
   delete h1W;
   delete h1N;

   // Also we need to set the entries since they have not been correctly calculated during the projection
   // we can only set them to the effective entries
   p1->SetEntries( p1->GetEffectiveEntries() );

   return p1;
}


////////////////////////////////////////////////////////////////////////////////
/// Replace current statistics with the values in array stats

void TProfile2D::PutStats(Double_t *stats)
{
   fTsumw   = stats[0];
   fTsumw2  = stats[1];
   fTsumwx  = stats[2];
   fTsumwx2 = stats[3];
   fTsumwy  = stats[4];
   fTsumwy2 = stats[5];
   fTsumwxy = stats[6];
   fTsumwz  = stats[7];
   fTsumwz2 = stats[8];
}

////////////////////////////////////////////////////////////////////////////////
/// Reset contents of a Profile2D histogram.

void TProfile2D::Reset(Option_t *option)
{
   TH2D::Reset(option);
   fBinEntries.Reset();
   fBinSumw2.Reset();
   TString opt = option;
   opt.ToUpper();
   if (opt.Contains("ICE") && !opt.Contains("S")) return;
   fTsumwz = fTsumwz2 = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Profile histogram is resized along axis such that x is in the axis range.
///
/// The new axis limits are recomputed by doubling iteratively
/// the current axis range until the specified value x is within the limits.
/// The algorithm makes a copy of the histogram, then loops on all bins
/// of the old histogram to fill the extended histogram.
/// Takes into account errors (Sumw2) if any.
/// The axis must be extendable before invoking this function.
///
/// Ex: `h->GetXaxis()->SetCanExtend(kTRUE)`

void TProfile2D::ExtendAxis(Double_t x, TAxis *axis)
{
   TProfile2D* hold = TProfileHelper::ExtendAxis(this, x, axis);
   if ( hold ) {
      fTsumwz  = hold->fTsumwz;
      fTsumwz2 = hold->fTsumwz2;
      delete hold;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Rebin this histogram grouping nxgroup/nygroup bins along the xaxis/yaxis together.
///
/// if newname is not blank a new profile hnew is created.
/// else the current histogram is modified (default)
/// The parameter nxgroup/nygroup indicate how many bins along the xaxis/yaxis of this
/// have to be merged into one bin of hnew
/// If the original profile has errors stored (via Sumw2), the resulting
/// profile has new errors correctly calculated.
///
/// examples: if hpxpy is an existing TProfile2D profile with 40 x 40 bins
/// ~~~ {.cpp}
///      hpxpy->Rebin2D();  // merges two bins along the xaxis and yaxis in one
///                         // Carefull: previous contents of hpxpy are lost
///      hpxpy->Rebin2D(3,5);  // merges 3 bins along the xaxis and 5 bins along the yaxis in one
///                            // Carefull: previous contents of hpxpy are lost
///      hpxpy->RebinX(5); //merges five bins along the xaxis in one in hpxpy
///      TProfile2D *hnew = hpxpy->RebinY(5,"hnew"); // creates a new profile hnew
///                                                  // merging 5 bins of hpxpy along the yaxis in one bin
/// ~~~
///
///  NOTE : If nxgroup/nygroup is not an exact divider of the number of bins,
///         along the xaxis/yaxis the top limit(s) of the rebinned profile
///         is changed to the upper edge of the xbin=newxbins*nxgroup resp.
///          ybin=newybins*nygroup and the remaining bins are added to
///          the overflow bin.
///          Statistics will be recomputed from the new bin contents.

TProfile2D * TProfile2D::Rebin2D(Int_t nxgroup ,Int_t nygroup,const char * newname ) {
   //something to do?
   if((nxgroup != 1) || (nygroup != 1)){
      Int_t nxbins  = fXaxis.GetNbins();
      Int_t nybins  = fYaxis.GetNbins();
      Double_t xmin  = fXaxis.GetXmin();
      Double_t xmax  = fXaxis.GetXmax();
      Double_t ymin  = fYaxis.GetXmin();
      Double_t ymax  = fYaxis.GetXmax();
      if ((nxgroup <= 0) || (nxgroup > nxbins)) {
         Error("Rebin", "Illegal value of nxgroup=%d",nxgroup);
         return 0;
      }
      if ((nygroup <= 0) || (nygroup > nybins)) {
         Error("Rebin", "Illegal value of nygroup=%d",nygroup);
         return 0;
      }

      Int_t newxbins = nxbins/nxgroup;
      Int_t newybins = nybins/nygroup;

      //warning if bins are added to the overflow bin
      if(newxbins*nxgroup != nxbins) {
         Warning("Rebin", "nxgroup=%d should be an exact divider of nxbins=%d",nxgroup,nxbins);
      }
      if(newybins*nygroup != nybins) {
         Warning("Rebin", "nygroup=%d should be an exact divider of nybins=%d",nygroup,nybins);
      }

      //save old bin contents in new arrays
      Double_t *oldBins   = new Double_t[(nxbins+2)*(nybins+2)];
      Double_t *oldCount  = new Double_t[(nxbins+2)*(nybins+2)];
      Double_t *oldErrors = new Double_t[(nxbins+2)*(nybins+2)];
      Double_t *oldBinw2  = (fBinSumw2.fN ? new Double_t[(nxbins+2)*(nybins+2)] : 0  );
      Double_t *cu1 = GetW();
      Double_t *er1 = GetW2();
      Double_t *en1 = GetB();
      Double_t *ew1 = GetB2();
      for(Int_t ibin=0; ibin < (nxbins+2)*(nybins+2); ibin++){
         oldBins[ibin]   = cu1[ibin];
         oldCount[ibin]  = en1[ibin];
         oldErrors[ibin] = er1[ibin];
         if (ew1 && fBinSumw2.fN) oldBinw2[ibin]  = ew1[ibin];
      }

      // create a clone of the old profile if newname is specified
      TProfile2D *hnew = this;
      if(newname && strlen(newname) > 0) {
         hnew = (TProfile2D*)Clone(newname);
      }

      // in case of nxgroup/nygroup not an exact divider of nxbins/nybins,
      // top limit is changed (see NOTE in method comment)
      if(newxbins*nxgroup != nxbins) {
         xmax = fXaxis.GetBinUpEdge(newxbins*nxgroup);
         hnew->fTsumw = 0; //stats must be reset because top bins will be moved to overflow bin
      }
      if(newybins*nygroup != nybins) {
         ymax = fYaxis.GetBinUpEdge(newybins*nygroup);
         hnew->fTsumw = 0; //stats must be reset because top bins will be moved to overflow bin
      }

      //rebin the axis
      if((fXaxis.GetXbins()->GetSize() > 0) || (fYaxis.GetXbins()->GetSize() > 0)){
         Double_t* xbins = new Double_t[newxbins+1];
         Double_t* ybins = new Double_t[newybins+1];
         for(Int_t i=0; i < newxbins+1; i++)
            xbins[i] = fXaxis.GetBinLowEdge(1+i*nxgroup);
         for(Int_t j=0; j < newybins+1; j++)
            ybins[j] = fYaxis.GetBinLowEdge(1+j*nygroup);
         hnew->SetBins(newxbins,xbins,newybins,ybins);
         delete [] xbins;
         delete [] ybins;
      }
      //fixed bin size
      else{
         hnew->SetBins(newxbins,xmin,xmax,newybins,ymin,ymax);
      }

      //merge bins
      Double_t *cu2 = hnew->GetW();
      Double_t *er2 = hnew->GetW2();
      Double_t *en2 = hnew->GetB();
      Double_t *ew2 = hnew->GetB2();
      Double_t binContent, binCount, binError, binSumw2;
      //connection between x and y bin number and linear global bin number:
      //global bin = xbin + (nxbins+2) * ybin
      Int_t oldxbin = 1;
      Int_t oldybin = 1;
      //global bin number
      Int_t bin;
      for(Int_t xbin = 1; xbin <= newxbins; xbin++){
         oldybin = 1;
         for(Int_t ybin = 1; ybin <= newybins; ybin++){
            binContent = 0;
            binCount   = 0;
            binError   = 0;
            binSumw2   = 0;
            for(Int_t i=0; i < nxgroup; i++){
               if(oldxbin + i > nxbins) break;
               for(Int_t j=0; j < nygroup; j++){
                  if(oldybin + j > nybins) break;
                  bin = oldxbin + i + (nxbins+2)*(oldybin+j);
                  binContent += oldBins[bin];
                  binCount += oldCount[bin];
                  binError += oldErrors[bin];
                  if(fBinSumw2.fN) binSumw2 += oldBinw2[bin];
               }
            }
            bin = xbin + (newxbins + 2)*ybin;
            cu2[bin] = binContent;
            er2[bin] = binError;
            en2[bin] = binCount;
            if(fBinSumw2.fN) ew2[bin] = binSumw2;
            oldybin += nygroup;
         }
         oldxbin += nxgroup;
      }

      //copy the underflow bin in x and y (0,0)
      cu2[0] = oldBins[0];
      er2[0] = oldErrors[0];
      en2[0] = oldCount[0];
      if(fBinSumw2.fN) ew2[0] = oldBinw2[0];
      //calculate overflow bin in x and y (newxbins+1,newybins+1)
      //therefore the oldxbin and oldybin from above are needed!
      binContent = 0;
      binCount   = 0;
      binError   = 0;
      binSumw2   = 0;
      for(Int_t i=oldxbin; i <= nxbins+1; i++){
         for(Int_t j=oldybin; j <= nybins+1; j++){
            //global bin number
            bin = i + (nxbins+2)*j;
            binContent += oldBins[bin];
            binCount += oldCount[bin];
            binError += oldErrors[bin];
            if(fBinSumw2.fN) binSumw2 += oldBinw2[bin];
         }
      }
      bin = (newxbins+2)*(newybins+2)-1;
      cu2[bin] = binContent;
      er2[bin] = binError;
      en2[bin] = binCount;
      if(fBinSumw2.fN) ew2[bin] = binSumw2;
      //calculate overflow bin in x and underflow bin in y (newxbins+1,0)
      binContent = 0;
      binCount   = 0;
      binError   = 0;
      binSumw2   = 0;
      for(Int_t i=oldxbin; i <= nxbins+1; i++){
         bin = i;
         binContent += oldBins[bin];
         binCount += oldCount[bin];
         binError += oldErrors[bin];
         if(fBinSumw2.fN) binSumw2 += oldBinw2[bin];
      }
      bin = newxbins + 1;
      cu2[bin] = binContent;
      er2[bin] = binError;
      en2[bin] = binCount;
      if(fBinSumw2.fN) ew2[bin] = binSumw2;
      //calculate underflow bin in x and overflow bin in y (0,newybins+1)
      binContent = 0;
      binCount   = 0;
      binError   = 0;
      binSumw2   = 0;
      for(Int_t i=oldybin; i <= nybins+1; i++){
         bin = i*(nxbins + 2);
         binContent += oldBins[bin];
         binCount += oldCount[bin];
         binError += oldErrors[bin];
         if(fBinSumw2.fN) binSumw2 += oldBinw2[bin];
      }
      bin = (newxbins + 2)*(newybins + 1);
      cu2[bin] = binContent;
      er2[bin] = binError;
      en2[bin] = binCount;
      if(fBinSumw2.fN) ew2[bin] = binSumw2;
      //calculate under/overflow contents in y for the new x bins
      Double_t binContentuf, binCountuf, binErroruf, binSumw2uf;
      Double_t binContentof, binCountof, binErrorof, binSumw2of;
      Int_t ufbin, ofbin;
      Int_t oldxbin2 = 1;
      for(Int_t xbin = 1; xbin <= newxbins; xbin++){
         binContentuf = 0;
         binCountuf   = 0;
         binErroruf   = 0;
         binSumw2uf   = 0;
         binContentof = 0;
         binCountof   = 0;
         binErrorof   = 0;
         binSumw2of   = 0;
         for(Int_t i = 0; i < nxgroup; i++){
            //index of under/overflow bin for y in old binning
            ufbin = (oldxbin2 + i);
            binContentuf += oldBins[ufbin];
            binCountuf   += oldCount[ufbin];
            binErroruf   += oldErrors[ufbin];
            if(fBinSumw2.fN) binSumw2uf   += oldBinw2[ufbin];
            for(Int_t j = oldybin; j <= nybins+1; j++)
            {
               ofbin = ufbin + j*(nxbins + 2);
               binContentof += oldBins[ofbin];
               binCountof   += oldCount[ofbin];
               binErrorof   += oldErrors[ofbin];
               if(fBinSumw2.fN) binSumw2of   += oldBinw2[ofbin];
            }
         }
         //index of under/overflow bin for y in new binning
         ufbin = xbin;
         ofbin = ufbin + (newybins + 1)*(newxbins + 2);
         cu2[ufbin] = binContentuf;
         er2[ufbin] = binErroruf;
         en2[ufbin] = binCountuf;
         if(fBinSumw2.fN) ew2[ufbin] = binSumw2uf;
         cu2[ofbin] = binContentof;
         er2[ofbin] = binErrorof;
         en2[ofbin] = binCountof;
         if(fBinSumw2.fN) ew2[ofbin] = binSumw2of;

         oldxbin2 += nxgroup;
      }
      //calculate under/overflow contents in x for the new y bins
      Int_t oldybin2 = 1;
      for(Int_t ybin = 1; ybin <= newybins; ybin++){
         binContentuf = 0;
         binCountuf   = 0;
         binErroruf   = 0;
         binSumw2uf   = 0;
         binContentof = 0;
         binCountof   = 0;
         binErrorof   = 0;
         binSumw2of   = 0;
         for(Int_t i = 0; i < nygroup; i++){
            //index of under/overflow bin for x in old binning
            ufbin = (oldybin2 + i)*(nxbins+2);
            binContentuf += oldBins[ufbin];
            binCountuf   += oldCount[ufbin];
            binErroruf   += oldErrors[ufbin];
            if(fBinSumw2.fN) binSumw2uf   += oldBinw2[ufbin];
            for(Int_t j = oldxbin; j <= nxbins+1; j++)
            {
               ofbin = j + ufbin;
               binContentof += oldBins[ofbin];
               binCountof   += oldCount[ofbin];
               binErrorof   += oldErrors[ofbin];
               if(fBinSumw2.fN) binSumw2of   += oldBinw2[ofbin];
            }
         }
         //index of under/overflow bin for x in new binning
         ufbin = ybin * (newxbins + 2);
         ofbin = newxbins + 1 + ufbin;
         cu2[ufbin] = binContentuf;
         er2[ufbin] = binErroruf;
         en2[ufbin] = binCountuf;
         if(fBinSumw2.fN) ew2[ufbin] = binSumw2uf;
         cu2[ofbin] = binContentof;
         er2[ofbin] = binErrorof;
         en2[ofbin] = binCountof;
         if(fBinSumw2.fN) ew2[ofbin] = binSumw2of;

         oldybin2 += nygroup;
      }

      delete [] oldBins;
      delete [] oldCount;
      delete [] oldErrors;
      if (oldBinw2) delete [] oldBinw2;

      return hnew;
   }
   //nxgroup == nygroup == 1
   else{
      if((newname) && (strlen(newname) > 0))
         return (TProfile2D*)Clone(newname);
      else
         return this;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Rebin only the X axis.
/// see Rebin2D

TProfile2D * TProfile2D::RebinX(Int_t ngroup,const char * newname ) {
   return Rebin2D(ngroup,1,newname);
}

////////////////////////////////////////////////////////////////////////////////
/// Rebin only the Y axis.
/// see Rebin2D

TProfile2D * TProfile2D::RebinY(Int_t ngroup,const char * newname ) {
   return Rebin2D(1,ngroup,newname);
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out.
///
/// Note the following restrictions in the code generated:
///  - variable bin size not implemented
///  - SetErrorOption not implemented

void TProfile2D::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';
   out <<"   "<<std::endl;
   out <<"   "<<ClassName()<<" *";

   out << GetName() << " = new " << ClassName() << "(" << quote
   << GetName() << quote << "," << quote<< GetTitle() << quote
   << "," << GetXaxis()->GetNbins();
   out << "," << GetXaxis()->GetXmin()
   << "," << GetXaxis()->GetXmax();
   out << "," << GetYaxis()->GetNbins();
   out << "," << GetYaxis()->GetXmin()
   << "," << GetYaxis()->GetXmax();
   out << "," << fZmin
       << "," << fZmax;
   out << ");" << std::endl;


   // save bin entries
   Int_t bin;
   for (bin=0;bin<fNcells;bin++) {
      Double_t bi = GetBinEntries(bin);
      if (bi) {
         out<<"   "<<GetName()<<"->SetBinEntries("<<bin<<","<<bi<<");"<<std::endl;
      }
   }
   //save bin contents
   for (bin=0;bin<fNcells;bin++) {
      Double_t bc = fArray[bin];
      if (bc) {
         out<<"   "<<GetName()<<"->SetBinContent("<<bin<<","<<bc<<");"<<std::endl;
      }
   }
   // save bin errors
   if (fSumw2.fN) {
      for (bin=0;bin<fNcells;bin++) {
         Double_t be = TMath::Sqrt(fSumw2.fArray[bin]);
         if (be) {
            out<<"   "<<GetName()<<"->SetBinError("<<bin<<","<<be<<");"<<std::endl;
         }
      }
   }

   TH1::SavePrimitiveHelp(out, GetName(), option);
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply this profile2D by a constant c1.
///
/// `this = c1*this
///
/// This function uses the services of TProfile2D::Add

void TProfile2D::Scale(Double_t c1, Option_t * option)
{
   TProfileHelper::Scale(this, c1, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the number of entries in bin.

void TProfile2D::SetBinEntries(Int_t bin, Double_t w)
{
   TProfileHelper::SetBinEntries(this, bin, w);
}

////////////////////////////////////////////////////////////////////////////////
/// Redefine  x and y axis parameters.

void TProfile2D::SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax)
{
   TH1::SetBins(nx,xmin, xmax,ny, ymin,ymax);
   fBinEntries.Set(fNcells);
   if (fBinSumw2.fN) fBinSumw2.Set(fNcells);
}

////////////////////////////////////////////////////////////////////////////////
/// Redefine  x and y axis parameters for variable bin sizes.

void TProfile2D::SetBins(Int_t nx,  const Double_t *xbins, Int_t ny, const Double_t *ybins)
{
   TH1::SetBins(nx,xbins,ny,ybins);
   fBinEntries.Set(fNcells);
   if (fBinSumw2.fN) fBinSumw2.Set(fNcells);
}

////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow.
/// Reallocate bin contents array

void TProfile2D::SetBinsLength(Int_t n)
{
   TH2D::SetBinsLength(n);
   TProfileHelper::BuildArray(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the buffer size in units of 8 bytes (double).

void TProfile2D::SetBuffer(Int_t buffersize, Option_t *)
{
   if (fBuffer) {
      BufferEmpty();
      delete [] fBuffer;
      fBuffer = 0;
   }
   if (buffersize <= 0) {
      fBufferSize = 0;
      return;
   }
   if (buffersize < 100) buffersize = 100;
   fBufferSize = 1 + 4*buffersize;
   fBuffer = new Double_t[fBufferSize];
   memset(fBuffer,0,sizeof(Double_t)*fBufferSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Set option to compute profile2D errors.
///
/// The computation of the bin errors is based on the parameter option:
///     - ' '  (Default) The bin errors are the standard error on the mean of the bin profiled values (Z),
///                    i.e. the standard error of the bin contents.
///                    Note that if TProfile::Approximate()  is called, an approximation is used when
///                    the spread in Z is 0 and the number of bin entries  is > 0
///     - 's'            The bin errors are the standard deviations of the Z bin values
///                    Note that if TProfile::Approximate()  is called, an approximation is used when
///                    the spread in Z is 0 and the number of bin entries is > 0
///     - 'i'            Errors are as in default case (standard errors of the bin contents)
///                    The only difference is for the case when the spread in Z is zero.
///                    In this case for N > 0 the error is  1./SQRT(12.*N)
///     - 'g'            Errors are 1./SQRT(W)  for W not equal to 0 and 0 for W = 0.
///                    W is the sum in the bin of the weights of the profile.
///                    This option is for combining measurements z +/- dz,
///                    and  the profile is filled with values y and weights z = 1/dz**2
///
///   See TProfile::BuildOptions for a detailed explanation of all options

void TProfile2D::SetErrorOption(Option_t *option)
{
   TProfileHelper::SetErrorOption(this, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TProfile2D.

void TProfile2D::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TProfile2D::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TH2D::Streamer(R__b);
      fBinEntries.Streamer(R__b);
      Int_t errorMode;
      R__b >> errorMode;
      fErrorMode = (EErrorType)errorMode;
      if (R__v < 2) {
         Float_t zmin,zmax;
         R__b >> zmin; fZmin = zmin;
         R__b >> zmax; fZmax = zmax;
      } else {
         R__b >> fZmin;
         R__b >> fZmax;
      }
      R__b.CheckByteCount(R__s, R__c, TProfile2D::IsA());
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TProfile2D::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create/Delete structure to store sum of squares of weights per bin.
///
/// This is needed to compute  the correct statistical quantities
/// of a profile filled with weights
///
/// This function is automatically called when the histogram is created
/// if the static function TH1::SetDefaultSumw2 has been called before.
/// If flag is false the structure is deleted

void TProfile2D::Sumw2(Bool_t flag)
{
   TProfileHelper::Sumw2(this, flag);
}
