// @(#)root/hist:$Id$
// Author: Rene Brun   17/05/2006

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TProfile3D.h"
#include "TProfile2D.h"
#include "THashList.h"
#include "TMath.h"
#include "THLimitsFinder.h"
#include <iostream>
#include "TError.h"
#include "TClass.h"

#include "TProfileHelper.h"

Bool_t TProfile3D::fgApproximate = kFALSE;

ClassImp(TProfile3D);

/** \class TProfile3D
    \ingroup Hist
 Profile3D histograms are used to display the mean
 value of T and its RMS for each cell in X,Y,Z.
 Profile3D histograms are in many cases an
 The inter-relation of three measured quantities X, Y, Z and T can always
 be visualized by a four-dimensional histogram or scatter-plot;
 its representation on the line-printer is not particularly
 satisfactory, except for sparse data. If T is an unknown (but single-valued)
 approximate function of X,Y,Z this function is displayed by a profile3D histogram with
 much better precision than by a scatter-plot.

 The following formulae show the cumulated contents (capital letters) and the values
 displayed by the printing or plotting routines (small letters) of the elements for cell I, J.

                                                               2
       H(I,J,K)  =  sum T                      E(I,J,K)  =  sum T
       l(I,J,K)  =  sum l                      L(I,J,K)  =  sum l
       h(I,J,K)  =  H(I,J,K)/L(I,J,K)          s(I,J,K)  =  sqrt(E(I,J,K)/L(I,J,K)- h(I,J,K)**2)
       e(I,J,K)  =  s(I,J,K)/sqrt(L(I,J,K))

 In the special case where s(I,J,K) is zero (eg, case of 1 entry only in one cell)
 e(I,J,K) is computed from the average of the s(I,J,K) for all cells,
 if the static function TProfile3D::Approximate has been called.
 This simple/crude approximation was suggested in order to keep the cell
 during a fit operation. But note that this approximation is not the default behaviour.

 Example of a profile3D histogram
~~~~{.cpp}
{
    auto c1 = new TCanvas("c1","Profile histogram example",200,10,700,500);
    auto hprof3d  = new TProfile3D("hprof3d","Profile of pt versus px, py and pz",40,-4,4,40,-4,4,40,0,20);
    Double_t px, py, pz, pt;
    TRandom3 r(0);
    for ( Int_t i=0; i<25000; i++) {
       r.Rannor(px,py);
       pz = px*px + py*py;
       pt = r.Landau(0,1);
       hprof3d->Fill(px,py,pz,pt,1);
    }
    hprof3d->Draw();
}
~~~~
 NOTE: A TProfile3D is drawn as it was a simple TH3
*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for Profile3D histograms.

TProfile3D::TProfile3D() : TH3D()
{
   fTsumwt = fTsumwt2 = 0;
   fScaling = kFALSE;
   BuildOptions(0,0,"");
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor for Profile3D histograms.

TProfile3D::~TProfile3D()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Normal Constructor for Profile histograms.
///
/// The first eleven parameters are similar to TH3D::TH3D.
/// All values of t are accepted at filling time.
/// To fill a profile3D histogram, one must use TProfile3D::Fill function.
///
/// Note that when filling the profile histogram the function Fill
/// checks if the variable t is between fTmin and fTmax.
/// If a minimum or maximum value is set for the T scale before filling,
/// then all values below tmin or above tmax will be discarded.
/// Setting the minimum or maximum value for the T scale before filling
/// has the same effect as calling the special TProfile3D constructor below
/// where tmin and tmax are specified.
///
/// H(I,J,K) is printed as the cell contents. The errors computed are s(I,J,K) if CHOPT='S'
/// (spread option), or e(I,J,K) if CHOPT=' ' (error on mean).
///
/// See TProfile3D::BuildOptions for explanation of parameters
///
/// see other constructors below with all possible combinations of
/// fix and variable bin size like in TH3D.

TProfile3D::TProfile3D(const char *name,const char *title,Int_t nx,Double_t xlow,Double_t xup,Int_t ny,Double_t ylow,Double_t yup,Int_t nz, Double_t zlow,Double_t zup,Option_t *option)
    : TH3D(name,title,nx,xlow,xup,ny,ylow,yup,nz,zlow,zup)
{
   BuildOptions(0,0,option);
   if (xlow >= xup || ylow >= yup || zlow >= zup) SetBuffer(fgBufferSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 3-D Profile with variable bins in X , Y and Z.

TProfile3D::TProfile3D(const char *name,const char *title,Int_t nx,const Double_t *xbins,Int_t ny,const Double_t *ybins,Int_t nz,const Double_t *zbins,Option_t *option)
    : TH3D(name,title,nx,xbins,ny,ybins,nz,zbins)
{
   BuildOptions(0,0,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Set Profile3D histogram structure and options.
///
///  - tmin:  minimum value allowed for t
///  - tmax:  maximum value allowed for t
///            if (tmin = tmax = 0) there are no limits on the allowed t values (tmin = -inf, tmax = +inf)
///
///  - option:  this is the option for the computation of the t error of the profile ( TProfile3D::GetBinError )
///             possible values for the options are documented in TProfile3D::SetErrorOption
///
/// see also TProfile::BuildOptions for a detailed description

void TProfile3D::BuildOptions(Double_t tmin, Double_t tmax, Option_t *option)
{
   SetErrorOption(option);

   // create extra profile data structure (bin entries/ y^2 and sum of weight square)
   TProfileHelper::BuildArray(this);

   fTmin = tmin;
   fTmax = tmax;
   fScaling = kFALSE;
   fTsumwt  = fTsumwt2 = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TProfile3D::TProfile3D(const TProfile3D &profile) : TH3D()
{
   ((TProfile3D&)profile).Copy(*this);
}

TProfile3D &TProfile3D::operator=(const TProfile3D &profile)
{
   ((TProfile3D &)profile).Copy(*this);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: `this = this + c1*f1` .

Bool_t TProfile3D::Add(TF1 *, Double_t , Option_t*)
{
   Error("Add","Function not implemented for TProfile3D");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: `this = this + c1*h1` .

Bool_t TProfile3D::Add(const TH1 *h1, Double_t c1)
{
   if (!h1) {
      Error("Add","Attempt to add a non-existing profile");
      return kFALSE;
   }
   if (!h1->InheritsFrom(TProfile3D::Class())) {
      Error("Add","Attempt to add a non-profile2D object");
      return kFALSE;
   }

   return TProfileHelper::Add(this, this, h1, 1, c1);
}

////////////////////////////////////////////////////////////////////////////////
/// Replace contents of this profile3D by the addition of h1 and h2.
///
/// `this = c1*h1 + c2*h2`

Bool_t TProfile3D::Add(const TH1 *h1, const TH1 *h2, Double_t c1, Double_t c2)
{
   if (!h1 || !h2) {
      Error("Add","Attempt to add a non-existing profile");
      return kFALSE;
   }
   if (!h1->InheritsFrom(TProfile3D::Class())) {
      Error("Add","Attempt to add a non-profile3D object");
      return kFALSE;
   }
   if (!h2->InheritsFrom(TProfile3D::Class())) {
      Error("Add","Attempt to add a non-profile3D object");
      return kFALSE;
   }

   return TProfileHelper::Add(this, h1, h2, c1, c2);
}


////////////////////////////////////////////////////////////////////////////////
/// Set the fgApproximate flag.
///
/// When the flag is true, the function GetBinError
/// will approximate the bin error with the average profile error on all bins
/// in the following situation only
///
///  - the number of bins in the profile3D is less than 10404 (eg 100x100x100)
///  - the bin number of entries is small ( <5)
///  - the estimated bin error is extremely small compared to the bin content
///  (see TProfile3D::GetBinError)

void TProfile3D::Approximate(Bool_t approx)
{
   fgApproximate = approx;
}


////////////////////////////////////////////////////////////////////////////////
/// Fill histogram with all entries in the buffer.
///
///  - action = -1 histogram is reset and refilled from the buffer (called by THistPainter::Paint)
///  - action =  0 histogram is filled from the buffer
///  - action =  1 histogram is filled and buffer is deleted
///              The buffer is automatically deleted when the number of entries
///              in the buffer is greater than the number of entries in the histogram

Int_t TProfile3D::BufferEmpty(Int_t action)
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
      Double_t zmin = fBuffer[4];
      Double_t zmax = zmin;
      for (Int_t i=1;i<nbentries;i++) {
         Double_t x = fBuffer[5*i+2];
         if (x < xmin) xmin = x;
         if (x > xmax) xmax = x;
         Double_t y = fBuffer[5*i+3];
         if (y < ymin) ymin = y;
         if (y > ymax) ymax = y;
         Double_t z = fBuffer[5*i+4];
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
      Fill(buffer[5*i+2],buffer[5*i+3],buffer[5*i+4],buffer[5*i+5],buffer[5*i+1]);
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
///  Accumulate arguments in buffer.
///
/// When buffer is full, empty the buffer
///
///  - fBuffer[0] = number of entries in buffer
///  - fBuffer[1] = w of first entry
///  - fBuffer[2] = x of first entry
///  - fBuffer[3] = y of first entry
///  - fBuffer[4] = z of first entry
///  - fBuffer[5] = t of first entry

Int_t TProfile3D::BufferFill(Double_t x, Double_t y, Double_t z, Double_t t, Double_t w)
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
   if (5*nbentries+5 >= fBufferSize) {
      BufferEmpty(1);
      return Fill(x,y,z,t,w);
   }
   fBuffer[5*nbentries+1] = w;
   fBuffer[5*nbentries+2] = x;
   fBuffer[5*nbentries+3] = y;
   fBuffer[5*nbentries+4] = z;
   fBuffer[5*nbentries+5] = t;
   fBuffer[0] += 1;
   return -2;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a Profile3D histogram to a new profile2D histogram.

void TProfile3D::Copy(TObject &obj) const
{
   try {
      TProfile3D & pobj = dynamic_cast<TProfile3D&>(obj);

      TH3D::Copy(pobj);
      fBinEntries.Copy(pobj.fBinEntries);
      fBinSumw2.Copy(pobj.fBinSumw2);
      for (int bin=0;bin<fNcells;bin++) {
         pobj.fArray[bin]        = fArray[bin];
         pobj.fSumw2.fArray[bin] = fSumw2.fArray[bin];
      }
      pobj.fTmin = fTmin;
      pobj.fTmax = fTmax;
      pobj.fScaling = fScaling;
      pobj.fErrorMode = fErrorMode;
      pobj.fTsumwt  = fTsumwt;
      pobj.fTsumwt2 = fTsumwt2;

   } catch(...) {
      Fatal("Copy","Cannot copy a TProfile3D in a %s",obj.IsA()->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: `this = this/(c1*f1)` .
///
/// This function is not implemented

Bool_t TProfile3D::Divide(TF1 *, Double_t )
{
   Error("Divide","Function not implemented for TProfile3D");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this profile2D by h1.
///
/// `this = this/h1`
///
/// This function return kFALSE if the divide operation failed

Bool_t TProfile3D::Divide(const TH1 *h1)
{
   if (!h1) {
      Error("Divide","Attempt to divide a non-existing profile2D");
      return kFALSE;
   }
   if (!h1->InheritsFrom(TProfile3D::Class())) {
      Error("Divide","Attempt to divide a non-profile3D object");
      return kFALSE;
   }
   TProfile3D *p1 = (TProfile3D*)h1;

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
   Int_t nz = GetNbinsZ();
   if (nz != p1->GetNbinsZ()) {
      Error("Divide","Attempt to divide profiles with different number of bins");
      return kFALSE;
   }

// Reset statistics
   fEntries = fTsumw   = fTsumw2 = fTsumwx = fTsumwx2 = 0;

// Loop on bins (including underflows/overflows)
   Int_t bin,binx,biny,binz;
   Double_t *cu1 = p1->GetW();
   Double_t *er1 = p1->GetW2();
   Double_t *en1 = p1->GetB();
   Double_t c0,c1,w,u,x,y,z;
   for (binx =0;binx<=nx+1;binx++) {
      for (biny =0;biny<=ny+1;biny++) {
         for (binz =0;binz<=nz+1;binz++) {
            bin   = GetBin(binx,biny,binz);
            c0  = fArray[bin];
            c1  = cu1[bin];
            if (c1) w = c0/c1;
            else    w = 0;
            fArray[bin] = w;
            u = TMath::Abs(w);
            x = fXaxis.GetBinCenter(binx);
            y = fYaxis.GetBinCenter(biny);
            z = fZaxis.GetBinCenter(binz);
            fEntries++;
            fTsumw   += u;
            fTsumw2  += u*u;
            fTsumwx  += u*x;
            fTsumwx2 += u*x*x;
            fTsumwy  += u*y;
            fTsumwy2 += u*y*y;
            fTsumwxy += u*x*y;
            fTsumwz  += u;
            fTsumwz2 += u*z;
            fTsumwxz += u*x*z;
            fTsumwyz += u*y*z;
            fTsumwt  += u;
            fTsumwt2 += u*u;
            Double_t e0 = fSumw2.fArray[bin];
            Double_t e1 = er1[bin];
            Double_t c12= c1*c1;
            if (!c1) fSumw2.fArray[bin] = 0;
            else     fSumw2.fArray[bin] = (e0*c1*c1 + e1*c0*c0)/(c12*c12);
            if (!en1[bin]) fBinEntries.fArray[bin] = 0;
            else           fBinEntries.fArray[bin] /= en1[bin];
         }
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
/// This function return kFALSE if the divide operation failed

Bool_t TProfile3D::Divide(const TH1 *h1, const TH1 *h2, Double_t c1, Double_t c2, Option_t *option)
{
   TString opt = option;
   opt.ToLower();
   Bool_t binomial = kFALSE;
   if (opt.Contains("b")) binomial = kTRUE;
   if (!h1 || !h2) {
      Error("Divide","Attempt to divide a non-existing profile2D");
      return kFALSE;
   }
   if (!h1->InheritsFrom(TProfile3D::Class())) {
      Error("Divide","Attempt to divide a non-profile2D object");
      return kFALSE;
   }
   TProfile3D *p1 = (TProfile3D*)h1;
   if (!h2->InheritsFrom(TProfile3D::Class())) {
      Error("Divide","Attempt to divide a non-profile2D object");
      return kFALSE;
   }
   TProfile3D *p2 = (TProfile3D*)h2;

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
   Int_t nz = GetNbinsZ();
   if (nz != p1->GetNbinsZ() || nz != p2->GetNbinsZ()) {
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
   Int_t bin,binx,biny,binz;
   Double_t *cu1 = p1->GetW();
   Double_t *cu2 = p2->GetW();
   Double_t *er1 = p1->GetW2();
   Double_t *er2 = p2->GetW2();
   Double_t *en1 = p1->GetB();
   Double_t *en2 = p2->GetB();
   Double_t b1,b2,w,u,x,y,z,ac1,ac2;
   ac1 = TMath::Abs(c1);
   ac2 = TMath::Abs(c2);
   for (binx =0;binx<=nx+1;binx++) {
      for (biny =0;biny<=ny+1;biny++) {
         for (binz =0;binz<=nz+1;binz++) {
            bin   = GetBin(binx,biny,binz);
            b1  = cu1[bin];
            b2  = cu2[bin];
            if (b2) w = c1*b1/(c2*b2);
            else    w = 0;
            fArray[bin] = w;
            u = TMath::Abs(w);
            x = fXaxis.GetBinCenter(binx);
            y = fYaxis.GetBinCenter(biny);
            z = fZaxis.GetBinCenter(biny);
            fEntries++;
            fTsumw   += u;
            fTsumw2  += u*u;
            fTsumwx  += u*x;
            fTsumwx2 += u*x*x;
            fTsumwy  += u*y;
            fTsumwy2 += u*y*y;
            fTsumwxy += u*x*y;
            fTsumwz  += u*z;
            fTsumwz2 += u*z*z;
            fTsumwxz += u*x*z;
            fTsumwyz += u*y*z;
            fTsumwt  += u;
            fTsumwt2 += u*u;
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
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a Profile3D histogram (no weights).

Int_t TProfile3D::Fill(Double_t x, Double_t y, Double_t z, Double_t t)
{
   if (fBuffer) return BufferFill(x,y,z,t,1);

   Int_t bin,binx,biny,binz;

   if (fTmin != fTmax) {
      if (t <fTmin || t> fTmax || TMath::IsNaN(t) ) return -1;
   }

   fEntries++;
   binx =fXaxis.FindBin(x);
   biny =fYaxis.FindBin(y);
   binz =fZaxis.FindBin(z);
   if (binx <0 || biny <0 || binz<0) return -1;
   bin  = GetBin(binx,biny,binz);
   AddBinContent(bin, t);
   fSumw2.fArray[bin] += (Double_t)t*t;
   fBinEntries.fArray[bin] += 1;
   if (fBinSumw2.fN)  fBinSumw2.fArray[bin] += 1;
   if (binx == 0 || binx > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (biny == 0 || biny > fYaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (binz == 0 || binz > fZaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
//printf("x=%g, y=%g, z=%g, t=%g, binx=%d, biny=%d, binz=%d, bin=%d\n",x,y,z,t,binx,biny,binz,bin);
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
   fTsumwt  += t;
   fTsumwt2 += t*t;
   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a Profile3D histogram with weights.

Int_t TProfile3D::Fill(Double_t x, Double_t y, Double_t z, Double_t t, Double_t w)
{
   if (fBuffer) return BufferFill(x,y,z,t,w);

   Int_t bin,binx,biny,binz;

   if (fTmin != fTmax) {
      if (t <fTmin || t> fTmax || TMath::IsNaN(t) ) return -1;
   }

   Double_t u= w; // (w > 0 ? w : -w);
   fEntries++;
   binx =fXaxis.FindBin(x);
   biny =fYaxis.FindBin(y);
   binz =fZaxis.FindBin(z);
   if (binx <0 || biny <0 || binz<0) return -1;
   bin  = GetBin(binx,biny,binz);
   AddBinContent(bin, u*t);
   fSumw2.fArray[bin] += u*t*t;
   if (!fBinSumw2.fN && u != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();  // must be called before accumulating the entries
   if (fBinSumw2.fN)  fBinSumw2.fArray[bin] += u*u;
   fBinEntries.fArray[bin] += u;
   if (binx == 0 || binx > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (biny == 0 || biny > fYaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   if (binz == 0 || binz > fZaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   fTsumw   += u;
   fTsumw2  += u*u;
   fTsumwx  += u*x;
   fTsumwx2 += u*x*x;
   fTsumwy  += u*y;
   fTsumwy2 += u*y*y;
   fTsumwxy += u*x*y;
   fTsumwz  += u*z;
   fTsumwz2 += u*z*z;
   fTsumwxz += u*x*z;
   fTsumwyz += u*y*z;
   fTsumwt  += u*t;
   fTsumwt2 += u*t*t;
   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Return bin content of a Profile3D histogram.

Double_t TProfile3D::GetBinContent(Int_t bin) const
{
   if (fBuffer) ((TProfile3D*)this)->BufferEmpty();

   if (bin < 0 || bin >= fNcells) return 0;
   if (fBinEntries.fArray[bin] == 0) return 0;
   if (!fArray) return 0;
   return fArray[bin]/fBinEntries.fArray[bin];
}

////////////////////////////////////////////////////////////////////////////////
/// Return bin entries of a Profile3D histogram.

Double_t TProfile3D::GetBinEntries(Int_t bin) const
{
   if (fBuffer) ((TProfile3D*)this)->BufferEmpty();

   if (bin < 0 || bin >= fNcells) return 0;
   return fBinEntries.fArray[bin];
}

////////////////////////////////////////////////////////////////////////////////
/// Return bin effective entries for a weighted filled Profile histogram.
///
/// In case of an unweighted profile, it is equivalent to the number of entries per bin
/// The effective entries is defined as the square of the sum of the weights divided by the
/// sum of the weights square.
/// TProfile::Sumw2() must be called before filling the profile with weights.
/// Only by calling this method the  sum of the square of the weights per bin is stored.

Double_t TProfile3D::GetBinEffectiveEntries(Int_t bin)
{
   return TProfileHelper::GetBinEffectiveEntries((TProfile3D*)this, bin);
}

////////////////////////////////////////////////////////////////////////////////
/// Return bin error of a Profile3D histogram.
///
/// ### Computing errors: A moving field
///
/// The computation of errors for a TProfile3D has evolved with the versions
/// of ROOT. The difficulty is in computing errors for bins with low statistics.
///
///  - prior to version 3.10, we had no special treatment of low statistic bins.
///    As a result, these bins had huge errors. The reason is that the
///    expression eprim2 is very close to 0 (rounding problems) or 0.
///  - The algorithm is modified/protected for the case
///    when a TProfile3D is projected (ProjectionX). The previous algorithm
///    generated a N^2 problem when projecting a TProfile3D with a large number of
///    bins (eg 100000).
///  - in version 3.10/02, a new static function TProfile::Approximate
///    is introduced to enable or disable (default) the approximation.
///    (see also comments in TProfile::GetBinError)

Double_t TProfile3D::GetBinError(Int_t bin) const
{
   return TProfileHelper::GetBinError((TProfile3D*)this, bin);
}

////////////////////////////////////////////////////////////////////////////////
/// Return option to compute profile2D errors.

Option_t *TProfile3D::GetErrorOption() const
{
   if (fErrorMode == kERRORSPREAD)  return "s";
   if (fErrorMode == kERRORSPREADI) return "i";
   if (fErrorMode == kERRORSPREADG) return "g";
   return "";
}

////////////////////////////////////////////////////////////////////////////////
/// fill the array stats from the contents of this profile.
///
/// The array stats must be correctly dimensionned in the calling program.
///
///  - stats[0] = sumw
///  - stats[1] = sumw2
///  - stats[2] = sumwx
///  - stats[3] = sumwx2
///  - stats[4] = sumwy
///  - stats[5] = sumwy2
///  - stats[6] = sumwxy
///  - stats[7] = sumwz
///  - stats[8] = sumwz2
///  - stats[9] = sumwxz
///  - stats[10]= sumwyz
///  - stats[11]= sumwt
///  - stats[12]= sumwt2
///
/// If no axis-subrange is specified (via TAxis::SetRange), the array stats
/// is simply a copy of the statistics quantities computed at filling time.
/// If a sub-range is specified, the function recomputes these quantities
/// from the bin contents in the current axis range.

void TProfile3D::GetStats(Double_t *stats) const
{
   if (fBuffer) ((TProfile3D*)this)->BufferEmpty();

   // Loop on bins
   if ( (fTsumw == 0 /* && fEntries > 0 */) || fXaxis.TestBit(TAxis::kAxisRange) || fYaxis.TestBit(TAxis::kAxisRange)) {

      // check for labels axis . In that case corresponsing statistics do not make sense and it is set to zero
      Bool_t labelXaxis = ((const_cast<TAxis &>(fXaxis)).GetLabels() && fXaxis.CanExtend());
      Bool_t labelYaxis = ((const_cast<TAxis &>(fYaxis)).GetLabels() && fYaxis.CanExtend());
      Bool_t labelZaxis = ((const_cast<TAxis &>(fZaxis)).GetLabels() && fZaxis.CanExtend());

      Int_t bin, binx, biny,binz;
      Double_t w, w2;
      Double_t x,y,z;
      for (bin=0;bin<kNstat;bin++) stats[bin] = 0;
      if (!fBinEntries.fArray) return;
      for (binz=fZaxis.GetFirst();binz<=fZaxis.GetLast();binz++) {
         z = (!labelZaxis) ? fZaxis.GetBinCenter(binz) : 0;
         for (biny=fYaxis.GetFirst();biny<=fYaxis.GetLast();biny++) {
            y = (!labelYaxis) ? fYaxis.GetBinCenter(biny) : 0;
            for (binx=fXaxis.GetFirst();binx<=fXaxis.GetLast();binx++) {
               bin = GetBin(binx,biny,binz);
               w         = fBinEntries.fArray[bin];
               w2        = (fBinSumw2.fN ? fBinSumw2.fArray[bin] : w );
               x         = (!labelXaxis) ? fXaxis.GetBinCenter(binx) : 0;
               stats[0]  += w;
               stats[1]  += w2;
               stats[2]  += w*x;
               stats[3]  += w*x*x;
               stats[4]  += w*y;
               stats[5]  += w*y*y;
               stats[6]  += w*x*y;
               stats[7]  += w*z;
               stats[8]  += w*z*z;
               stats[9]  += w*x*z;
               stats[10] += w*y*z;
               stats[11] += fArray[bin];
               stats[12] += fSumw2.fArray[bin];
            }
         }
      }
   } else {
      stats[0]  = fTsumw;
      stats[1]  = fTsumw2;
      stats[2]  = fTsumwx;
      stats[3]  = fTsumwx2;
      stats[4]  = fTsumwy;
      stats[5]  = fTsumwy2;
      stats[6]  = fTsumwxy;
      stats[7]  = fTsumwz;
      stats[8]  = fTsumwz2;
      stats[9]  = fTsumwxz;
      stats[10] = fTsumwyz;
      stats[11] = fTsumwt;
      stats[12] = fTsumwt2;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reduce the number of bins for this axis to the number of bins having a label.

void TProfile3D::LabelsDeflate(Option_t *ax)
{
   TProfileHelper::LabelsDeflate(this, ax);
}

////////////////////////////////////////////////////////////////////////////////
/// Double the number of bins for axis.
/// Refill histogram
/// This function is called by TAxis::FindBin(const char *label)

void TProfile3D::LabelsInflate(Option_t *ax)
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

void TProfile3D::LabelsOption(Option_t */* option */, Option_t * /* ax */)
{
   Error("LabelsOption","Labels option function is not implemented for a TProfile3D");
}
////////////////////////////////////////////////////////////////////////////////
/// Merge all histograms in the collection in this histogram.
///
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

Long64_t TProfile3D::Merge(TCollection *li)
{
   return TProfileHelper::Merge(this, li);
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: `this = this*c1*f1` .

Bool_t TProfile3D::Multiply(TF1 *, Double_t )
{
   Error("Multiply","Function not implemented for TProfile3D");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply this profile2D by h1.
///
/// `this = this*h1`

Bool_t TProfile3D::Multiply(const TH1 *)
{
   Error("Multiply","Multiplication of profile2D histograms not implemented");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace contents of this profile2D by multiplication of h1 by h2.
///
/// `this = (c1*h1)*(c2*h2)`

Bool_t TProfile3D::Multiply(const TH1 *, const TH1 *, Double_t, Double_t, Option_t *)
{
   Error("Multiply","Multiplication of profile2D histograms not implemented");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Project this profile3D into a 3-D histogram along X,Y,Z.
///
///   The projection is always of the type TH3D.
///
///  - if option "E" is specified, the errors are computed. (default)
///  - if option "B" is specified, the content of bin of the returned histogram
///      will be equal to the GetBinEntries(bin) of the profile,
///  - if option "C=E" the bin contents of the projection are set to the
///       bin errors of the profile
///  - if option "E" is specified  the errors of the projected histogram are computed and set
///      to be equal to the errors of the profile.
///      Option "E" is defined as the default one in the header file.
///  - if option "" is specified the histogram errors are simply the sqrt of its content
///  - if option "B" is specified, the content of bin of the returned histogram
///      will be equal to the GetBinEntries(bin) of the profile,
///  - if option "C=E" the bin contents of the projection are set to the
///       bin errors of the profile
///  - if option "W" is specified the bin content of the projected histogram  is set to the
///       product of the bin content of the profile and the entries.
///       With this option the returned histogram will be equivalent to the one obtained by
///       filling directly a TH2D using the 3-rd value as a weight.
///       This option makes sense only for profile filled with all weights =1.
///       When the profile is weighted (filled with weights different than 1) the
///       bin error of the projected histogram (obtained using this option "W") cannot be
///       correctly computed from the information stored in the profile. In that case the
///       obtained histogram contains as bin error square the weighted sum of the square of the
///       profiled observable (TProfile2D::fSumw2[bin] )
///
///       Note that the axis range is not considered when doing the projection

TH3D *TProfile3D::ProjectionXYZ(const char *name, Option_t *option) const
{

   TString opt = option;
   opt.ToLower();
   Int_t nx = fXaxis.GetNbins();
   Int_t ny = fYaxis.GetNbins();
   Int_t nz = fZaxis.GetNbins();
   const TArrayD *xbins = fXaxis.GetXbins();
   const TArrayD *ybins = fYaxis.GetXbins();
   const TArrayD *zbins = fZaxis.GetXbins();

   // Create the projection histogram
   TString pname = name;
   if (pname == "_px") {
      pname = GetName();  pname.Append("_pxyz");
   }
   TH3D *h1 = 0 ;
   if (xbins->fN == 0 && ybins->fN == 0 && zbins->fN == 0)
      h1 = new TH3D(pname,GetTitle(),nx,fXaxis.GetXmin(),fXaxis.GetXmax(),ny,fYaxis.GetXmin(),fYaxis.GetXmax(),nz,fZaxis.GetXmin(),fZaxis.GetXmax());
   else if ( xbins->fN != 0 && ybins->fN != 0 && zbins->fN != 0)
      h1 = new TH3D(pname,GetTitle(),nx,xbins->GetArray(),ny,ybins->GetArray(), nz,zbins->GetArray() );
   else {
      Error("ProjectionXYZ","Histogram has an axis with variable bins and an axis with fixed bins. This case is not supported - return a null pointer");
      return 0;
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
   Int_t bin,binx,biny,binz;
   Double_t cont;
   for (binx =0;binx<=nx+1;binx++) {
      for (biny =0;biny<=ny+1;biny++) {
         for (binz =0;binz<=nz+1;binz++) {
            bin = GetBin(binx,biny,binz);

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
   }
   h1->SetEntries(fEntries);
   return h1;
}
////////////////////////////////////////////////////////////////////////////////
/// Project a 3-D profile into a 2D-profile histogram depending on the option parameter.
///
/// option may contain a combination of the characters x,y,z:
///
///  - option = "xy" return the x versus y projection into a TProfile2D histogram
///  - option = "yx" return the y versus x projection into a TProfile2D histogram
///  - option = "xz" return the x versus z projection into a TProfile2D histogram
///  - option = "zx" return the z versus x projection into a TProfile2D histogram
///  - option = "yz" return the y versus z projection into a TProfile2D histogram
///  - option = "zy" return the z versus y projection into a TProfile2D histogram
///
/// NB: the notation "a vs b" means "a" vertical and "b" horizontal along X
///
/// The resulting profile contains the combination of all the considered bins along X
/// By default, all bins are included considering also underflow/overflows
///
/// The option can also be used to specify the projected profile error type.
/// Values which can be used are 's', 'i', or 'g'. See TProfile::BuildOptions for details
///
/// To select a bin range along an axis, use TAxis::SetRange, eg
/// `h3.GetYaxis()->SetRange(23,56);`

TProfile2D *TProfile3D::Project3DProfile(Option_t *option) const
{
   // can call TH3 method which will call the virtual method :DoProjectProfile2D re-implemented below
   // but need to add underflow/overflow
   TString opt(option);
   opt.Append(" UF OF");
   return TH3::Project3DProfile(opt);
}

////////////////////////////////////////////////////////////////////////////////
/// Internal method to project to a 2D Profile.
///
/// Called from TH3::Project3DProfile but re-implemented in case of the TPRofile3D
/// since what is done is different.

TProfile2D *TProfile3D::DoProjectProfile2D(const char* name, const char * title, const TAxis* projX, const TAxis* projY,
                                           bool originalRange, bool useUF, bool useOF) const
{
   // Get the ranges where we will work.
   Int_t ixmin = projX->GetFirst();
   Int_t ixmax = projX->GetLast();
   Int_t iymin = projY->GetFirst();
   Int_t iymax = projY->GetLast();
   if (ixmin == 0 && ixmax == 0) { ixmin = 1; ixmax = projX->GetNbins(); }
   if (iymin == 0 && iymax == 0) { iymin = 1; iymax = projY->GetNbins(); }
   Int_t nx = ixmax-ixmin+1;
   Int_t ny = iymax-iymin+1;

   // Create the projected profiles
   TProfile2D *p2 = 0;
   // Create always a new TProfile2D (not as in the case of TH3 projection)

   const TArrayD *xbins = projX->GetXbins();
   const TArrayD *ybins = projY->GetXbins();
   // assume all axis have variable bins or have fixed bins
   if ( originalRange ) {
      if (xbins->fN == 0 && ybins->fN == 0) {
         p2 = new TProfile2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                             ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
      } else {
         p2 = new TProfile2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1],projX->GetNbins(),&xbins->fArray[ixmin-1]);
      }
   } else {
      if (xbins->fN == 0 && ybins->fN == 0) {
         p2 = new TProfile2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                             ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
      } else {
         p2 = new TProfile2D(name,title,ny,&ybins->fArray[iymin-1],nx,&xbins->fArray[ixmin-1]);
      }
   }

   // weights
   bool useWeights = (fBinSumw2.fN != 0);
   if (useWeights) p2->Sumw2();

   // make projection in a 3D first
   TH3D * h3dW = ProjectionXYZ("h3temp-W","W");
   TH3D * h3dN = ProjectionXYZ("h3temp-N","B");

   h3dW->SetDirectory(0); h3dN->SetDirectory(0);

   // Since no axis range is considered when doing the projection TProfile3D->TH3D
   // the resulting histogram will have the same axis as the parent one
   // we need afterwards to set the range in the 3D histogram to considered it later
   // when doing the projection in a Profile2D
   if (fXaxis.TestBit(TAxis::kAxisRange) ) {
      h3dW->GetXaxis()->SetRange(fXaxis.GetFirst(),fXaxis.GetLast());
      h3dN->GetXaxis()->SetRange(fXaxis.GetFirst(),fXaxis.GetLast());
   }
   if (fYaxis.TestBit(TAxis::kAxisRange) ) {
      h3dW->GetYaxis()->SetRange(fYaxis.GetFirst(),fYaxis.GetLast());
      h3dN->GetYaxis()->SetRange(fYaxis.GetFirst(),fYaxis.GetLast());
   }
   if (fZaxis.TestBit(TAxis::kAxisRange) ) {
      h3dW->GetZaxis()->SetRange(fZaxis.GetFirst(),fZaxis.GetLast());
      h3dN->GetZaxis()->SetRange(fZaxis.GetFirst(),fZaxis.GetLast());
   }

   // note that h3dW is always a weighted histogram - so we need to compute error in the projection
   TAxis * projX_hW = h3dW->GetXaxis();
   TAxis * projX_hN = h3dN->GetXaxis();
   if (projX == GetYaxis() ) {  projX_hW =  h3dW->GetYaxis();  projX_hN =  h3dN->GetYaxis(); }
   if (projX == GetZaxis() ) {  projX_hW =  h3dW->GetZaxis();  projX_hN =  h3dN->GetZaxis(); }
   TAxis * projY_hW = h3dW->GetYaxis();
   TAxis * projY_hN = h3dN->GetYaxis();
   if (projY == GetXaxis() ) {  projY_hW =  h3dW->GetXaxis();  projY_hN =  h3dN->GetXaxis(); }
   if (projY == GetZaxis() ) {  projY_hW =  h3dW->GetZaxis();  projY_hN =  h3dN->GetZaxis(); }

   TH2D * h2W = TH3::DoProject2D(*h3dW,"htemp-W","",projX_hW, projY_hW, true, originalRange, useUF, useOF);
   TH2D * h2N = TH3::DoProject2D(*h3dN,"htemp-N","",projX_hN, projY_hN, useWeights, originalRange, useUF, useOF);
   h2W->SetDirectory(0); h2N->SetDirectory(0);


   // fill the bin content
   R__ASSERT( h2W->fN == p2->fN );
   R__ASSERT( h2N->fN == p2->fN );
   R__ASSERT( h2W->GetSumw2()->fN != 0); // h2W should always be a weighted histogram since h3dW is weighted
   for (int i = 0; i < p2->fN ; ++i) {
      //std::cout << " proj bin " << i << "  " <<  h2W->fArray[i] << "  " << h2N->fArray[i] << std::endl;
      p2->fArray[i] = h2W->fArray[i];   // array of profile is sum of all values
      p2->GetSumw2()->fArray[i]  = h2W->GetSumw2()->fArray[i];   // array of content square of profile is weight square of the W projected histogram
      p2->SetBinEntries(i, h2N->fArray[i] );
      if (useWeights) p2->GetBinSumw2()->fArray[i] = h2N->GetSumw2()->fArray[i];    // sum of weight squares are stored to compute errors in h1N histogram
   }
   // delete the created histograms
   delete h3dW;
   delete h3dN;
   delete h2W;
   delete h2N;

   // Also we need to set the entries since they have not been correctly calculated during the projection
   // we can only set them to the effective entries
   p2->SetEntries( p2->GetEffectiveEntries() );

   return p2;

}

////////////////////////////////////////////////////////////////////////////////
/// Replace current statistics with the values in array stats.

void TProfile3D::PutStats(Double_t *stats)
{
   TH3::PutStats(stats);
   fTsumwt  = stats[11];
   fTsumwt2 = stats[12];
}

////////////////////////////////////////////////////////////////////////////////
/// Reset contents of a Profile3D histogram.

void TProfile3D::Reset(Option_t *option)
{
   TH3D::Reset(option);
   fBinSumw2.Reset();
   fBinEntries.Reset();
   TString opt = option;
   opt.ToUpper();
   if (opt.Contains("ICE") && !opt.Contains("S")) return;
   fTsumwt = fTsumwt2 = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Profile histogram is resized along axis such that x is in the axis range.
/// The new axis limits are recomputed by doubling iteratively
/// the current axis range until the specified value x is within the limits.
/// The algorithm makes a copy of the histogram, then loops on all bins
/// of the old histogram to fill the rebinned histogram.
/// Takes into account errors (Sumw2) if any.
/// The axis must be rebinnable before invoking this function.
/// Ex: `h->GetXaxis()->SetCanExtend(kTRUE)`

void TProfile3D::ExtendAxis(Double_t x, TAxis *axis)
{
   TProfile3D* hold = TProfileHelper::ExtendAxis(this, x, axis);
   if ( hold ) {
      fTsumwt  = hold->fTsumwt;
      fTsumwt2 = hold->fTsumwt2;
      delete hold;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out.
///
/// Note the following restrictions in the code generated:
///  - variable bin size not implemented
///  - SetErrorOption not implemented

void TProfile3D::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
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
   out << "," << GetZaxis()->GetNbins();
   out << "," << GetZaxis()->GetXmin()
       << "," << GetZaxis()->GetXmax();
   out << "," << fTmin
       << "," << fTmax;
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
/// `this = c1*this`
///
/// This function uses the services of TProfile3D::Add

void TProfile3D::Scale(Double_t c1, Option_t *option)
{
   TProfileHelper::Scale(this, c1, option);
}

////////////////////////////////////////////////////////////////////////////////
///Set the number of entries in bin.

void TProfile3D::SetBinEntries(Int_t bin, Double_t w)
{
   TProfileHelper::SetBinEntries(this, bin, w);
}

////////////////////////////////////////////////////////////////////////////////
/// Redefine  x, y and z axis parameters.

void TProfile3D::SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax, Int_t nz, Double_t zmin, Double_t zmax)
{
   TH1::SetBins(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
   fBinEntries.Set(fNcells);
   if (fBinSumw2.fN) fBinSumw2.Set(fNcells);
}

////////////////////////////////////////////////////////////////////////////////
/// Redefine  x, y and z axis parameters with variable bin sizes

void TProfile3D::SetBins(Int_t nx, const Double_t *xBins, Int_t ny, const Double_t *yBins, Int_t nz, const Double_t *zBins)
{
   TH1::SetBins(nx,xBins,ny,yBins,nz,zBins);
   fBinEntries.Set(fNcells);
   if (fBinSumw2.fN) fBinSumw2.Set(fNcells);
}

////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow.
///
/// Reallocate bin contents array

void TProfile3D::SetBinsLength(Int_t n)
{
   TH3D::SetBinsLength(n);
   TProfileHelper::BuildArray(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the buffer size in units of 8 bytes (double).

void TProfile3D::SetBuffer(Int_t buffersize, Option_t *)
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
   fBufferSize = 1 + 5*buffersize;
   fBuffer = new Double_t[fBufferSize];
   memset(fBuffer,0,sizeof(Double_t)*fBufferSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Set option to compute profile3D errors.
///
/// The computation of the bin errors is based on the parameter option:
///  - ' '  (Default) The bin errors are the standard error on the mean of the bin profiled values (T),
///         i.e. the standard error of the bin contents.
///         Note that if TProfile3D::Approximate()  is called, an approximation is used when
///         the spread in T is 0 and the number of bin entries  is > 0
///  - 's'  The bin errors are the standard deviations of the T bin values
///         Note that if TProfile3D::Approximate()  is called, an approximation is used when
///         the spread in T is 0 and the number of bin entries is > 0
///  - 'i'  Errors are as in default case (standard errors of the bin contents)
///         The only difference is for the case when the spread in T is zero.
///         In this case for N > 0 the error is  1./SQRT(12.*N)
///  - 'g'  Errors are 1./SQRT(W)  for W not equal to 0 and 0 for W = 0.
///         W is the sum in the bin of the weights of the profile.
///         This option is for combining measurements t +/- dt,
///         and  the profile is filled with values t and weights w = 1/dt**2
///
///   See TProfile::BuildOptions for explanation of all options

void TProfile3D::SetErrorOption(Option_t *option)
{
   TProfileHelper::SetErrorOption(this, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Create/Delete structure to store sum of squares of weights per bin
/// This is needed to compute  the correct statistical quantities
/// of a profile filled with weights
///
///  This function is automatically called when the histogram is created
///  if the static function TH1::SetDefaultSumw2 has been called before.
///  If flag = false the structure is deleted

void TProfile3D::Sumw2(Bool_t flag)
{
   TProfileHelper::Sumw2(this, flag);
}
