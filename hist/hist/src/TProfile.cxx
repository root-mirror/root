// @(#)root/hist:$Id$
// Author: Rene Brun   29/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TProfile.h"
#include "TBuffer.h"
#include "TMath.h"
#include "TF1.h"
#include "THLimitsFinder.h"
#include <iostream>
#include "TError.h"
#include "TClass.h"
#include "TObjString.h"

#include "TProfileHelper.h"

Bool_t TProfile::fgApproximate = kFALSE;

ClassImp(TProfile);

/** \class TProfile
    \ingroup Hist
 Profile Histogram.
 Profile histograms are used to display the mean
 value of Y and its error for each bin in X. The displayed error is by default the
 standard error on the mean (i.e. the standard deviation divided by the sqrt(n) ).
 Profile histograms are in many cases an
 elegant replacement of two-dimensional histograms. The inter-relation of two
 measured quantities X and Y can always be visualized by a two-dimensional
 histogram or scatter plot, but if Y is an unknown (but single-valued)
 approximate function of X, this function is displayed by a profile histogram with
 much better precision than by a scatter plot.

 The following formulae show the cumulated contents (capital letters) and the values
 displayed by the printing or plotting routines (small letters) of the elements for bin j.
 \f[
  \begin{align}
       H(j)  &=  \sum w \cdot Y \\
       E(j)  &=  \sum w \cdot Y^2 \\
       W(j)  &=  \sum w \\
       h(j)  &=  H(j) / W(j)              & &\text{mean of Y,} \\
       s(j)  &=  \sqrt{E(j)/W(j)- h(j)^2} & &\text{standard deviation of Y} \\
       e(j)  &=  s(j)/\sqrt{W(j)}         & &\text{standard error on the mean} \\
  \end{align}
 \f]
 The bin content is always the mean of the Y values, but errors change depending on options:
 \f[
    \begin{align}
      \text{GetBinContent}(j) &= h(j) \\
      \text{GetBinError}(j) &=
        \begin{cases}
          e(j)                 &\text{if option="" (default). Error of the mean of all y values.} \\
          s(j)                 &\text{if option="s". Standard deviation of all y values.} \\
          \begin{cases} e(j) &\text{if } h(j) \ne 0 \\ 1/\sqrt{12 N} &\text{if } h(j)=0 \end{cases}       &\text{if option="i". This is useful for storing integers such as ADC counts.} \\
          1/\sqrt{W(j)}           &\text{if option="g". Error of a weighted mean for combining measurements with variances of } w. \\
        \end{cases}
    \end{align}
 \f]
 In the special case where s(j) is zero (eg, case of 1 entry only in one bin)
 the bin error e(j) is computed from the average of the s(j) for all bins if
 the static function TProfile::Approximate() has been called.
 This simple/crude approximation was suggested in order to keep the bin
 during a fit operation. But note that this approximation is not the default behaviour.
 See also TProfile::BuildOptions for more on error options.

  ### Creating and drawing a profile histogram
~~~{.cpp}
{
  auto c1 = new TCanvas("c1","Profile histogram example",200,10,700,500);
  auto hprof  = new TProfile("hprof","Profile of pz versus px",100,-4,4,0,20);
  Float_t px, py, pz;
  for ( Int_t i=0; i<25000; i++) {
    gRandom->Rannor(px,py);
    pz = px*px + py*py;
    hprof->Fill(px,pz,1);
  }
  hprof->Draw();
}
~~~
*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for Profile histograms

TProfile::TProfile() : TH1D()
{
   BuildOptions(0,0,"");
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor for Profile histograms

TProfile::~TProfile()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Normal Constructor for Profile histograms.
///
/// The first five parameters are similar to TH1D::TH1D.
/// All values of y are accepted at filling time.
/// To fill a profile histogram, one must use TProfile::Fill function.
///
/// Note that when filling the profile histogram the function Fill
/// checks if the variable y is between fYmin and fYmax.
/// If a minimum or maximum value is set for the Y scale before filling,
/// then all values below ymin or above ymax will be discarded.
/// Setting the minimum or maximum value for the Y scale before filling
/// has the same effect as calling the special TProfile constructor below
/// where ymin and ymax are specified.
///
/// H(j) is printed as the channel contents. The errors displayed are s(j) if `option`='S'
/// (spread option), or e(j) if `CHOPT`='' (error on mean).
///
/// See TProfile::BuildOptions() for explanation of parameters
///
/// see also comments in the TH1 base class constructors

TProfile::TProfile(const char *name,const char *title,Int_t nbins,Double_t xlow,Double_t xup,Option_t *option)
: TH1D(name,title,nbins,xlow,xup)
{
   BuildOptions(0,0,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor for Profile histograms with variable bin size.
///
/// See TProfile::BuildOptions() for more explanations on errors
/// see also comments in the TH1 base class constructors

TProfile::TProfile(const char *name,const char *title,Int_t nbins,const Float_t *xbins,Option_t *option)
: TH1D(name,title,nbins,xbins)
{
   BuildOptions(0,0,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor for Profile histograms with variable bin size.
///
/// See TProfile::BuildOptions for more explanations on errors
/// see also comments in the TH1 base class constructors

TProfile::TProfile(const char *name,const char *title,Int_t nbins,const Double_t *xbins,Option_t *option)
: TH1D(name,title,nbins,xbins)
{
   BuildOptions(0,0,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor for Profile histograms with variable bin size.
/// See TProfile::BuildOptions for more explanations on errors
///
/// see also comments in the TH1 base class constructors

TProfile::TProfile(const char *name,const char *title,Int_t nbins,const Double_t *xbins,Double_t ylow,Double_t yup,Option_t *option)
: TH1D(name,title,nbins,xbins)
{
   BuildOptions(ylow,yup,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor for Profile histograms with range in y.
///
/// The first five parameters are similar to TH1D::TH1D.
/// Only the values of Y between ylow and yup will be considered at filling time.
/// ylow and yup will also be the maximum and minimum values
/// on the y scale when drawing the profile.
///
/// See TProfile::BuildOptions for more explanations on errors
///
/// see also comments in the TH1 base class constructors

TProfile::TProfile(const char *name,const char *title,Int_t nbins,Double_t xlow,Double_t xup,Double_t ylow,Double_t yup,Option_t *option)
: TH1D(name,title,nbins,xlow,xup)
{
   BuildOptions(ylow,yup,option);
}


////////////////////////////////////////////////////////////////////////////////
/// Set Profile histogram structure and options.
///
/// \param[in] ymin  minimum value allowed for y
/// \param[in] ymax  maximum value allowed for y
///            if (ymin = ymax = 0) there are no limits on the allowed y values (ymin = -inf, ymax = +inf)
/// \param[in] option this is the option for the computation of the y error of the profile ( TProfile::GetBinError )
///            possible values for the options are:
///         - ' '  (Default) the bin errors are the standard error on the mean of Y  =  S(Y)/SQRT(N)
///           where S(Y) is the standard deviation (RMS) of the Y data in the bin
///           and N is the number of bin entries (from TProfile::GetBinEntries(ibin) )
///           (i.e the errors are the standard error on the bin content of the profile)
///         - 's' Errors are the standard deviation of Y, S(Y)
///         - 'i' Errors are S(Y)/SQRT(N) (standard error on the mean as in the default)
///           The only difference is only when the standard deviation in Y is zero.
///           In this  case the error a standard deviation = 1/SQRT(12) is assumed and the error is
///           1./SQRT(12*N).
///           This approximation assumes that the Y values are integer (e.g. ADC counts)
///           and have an implicit uncertainty of y +/- 0.5. With the assumption that the probability that y
///           takes any value between y-0.5 and y+0.5 is uniform, its standard error is 1/SQRT(12)
///         - 'g' Errors are 1./SQRT(W) where W is the sum of the weights for the bin j
///           W is obtained as from TProfile::GetBinEntries(ibin)
///           This errors corresponds to the standard deviation of weighted mean where each
///           measurement Y is uncorrelated and has an error sigma, which is expressed in the
///           weight used to fill the Profile:  w = 1/sigma^2
///           The resulting  error in TProfile is then 1./SQRT( Sum(1./sigma^2) )
///
///   In the case of Profile filled weights and with TProfile::Sumw2() called,
///   STD(Y) is the standard deviation of the weighted sample Y and N is in this case the
///   number of effective entries (TProfile::GetBinEffectiveEntries(ibin) )
///
///   If a bin has N data points all with the same value Y (especially
///   possible when dealing with integers), the spread in Y for that bin
///   is zero, and the uncertainty assigned is also zero, and the bin is
///   ignored in making subsequent fits.
///   To avoid this problem one can use an approximation for the standard deviation S(Y),
///   by using the average of all the S(Y) of the other Profile bins. To use this approximation
///   one must call before TProfile::Approximate
///   This approximation applies only for the default and  the 's' options

void TProfile::BuildOptions(Double_t ymin, Double_t ymax, Option_t *option)
{
   SetErrorOption(option);

   // create extra profile data structure (bin entries/ y^2 and sum of weight square)
   TProfileHelper::BuildArray(this);

   fYmin = ymin;
   fYmax = ymax;
   fScaling = kFALSE;
   fTsumwy = fTsumwy2 = 0;

}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TProfile::TProfile(const TProfile &profile) : TH1D()
{
   ((TProfile&)profile).Copy(*this);
}

TProfile &TProfile::operator=(const TProfile &profile)
{
   ((TProfile &)profile).Copy(*this);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: this = this + c1*f1

Bool_t TProfile::Add(TF1 *, Double_t, Option_t * )
{
   Error("Add","Function not implemented for TProfile");
   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: this = this + c1*h1

Bool_t TProfile::Add(const TH1 *h1, Double_t c1)
{
   if (!h1) {
      Error("Add","Attempt to add a non-existing profile");
      return kFALSE;
   }
   if (!h1->InheritsFrom(TProfile::Class())) {
      Error("Add","Attempt to add a non-profile object");
      return kFALSE;
   }

   return TProfileHelper::Add(this, this, h1, 1, c1);
}

////////////////////////////////////////////////////////////////////////////////
/// Replace contents of this profile by the addition of h1 and h2.
///
/// `this = c1*h1 + c2*h2`
///
///  c1 and c2 are considered as weights applied to the two summed profiles.
///  The operation acts therefore like merging the two profiles with a weight c1 and c2

Bool_t TProfile::Add(const TH1 *h1, const TH1 *h2, Double_t c1, Double_t c2)
{
   if (!h1 || !h2) {
      Error("Add","Attempt to add a non-existing profile");
      return kFALSE;
   }
   if (!h1->InheritsFrom(TProfile::Class())) {
      Error("Add","Attempt to add a non-profile object");
      return kFALSE;
   }
   if (!h2->InheritsFrom(TProfile::Class())) {
      Error("Add","Attempt to add a non-profile object");
      return kFALSE;
   }
   return TProfileHelper::Add(this, h1, h2, c1, c2);
}


////////////////////////////////////////////////////////////////////////////////
/// Static function to set the fgApproximate flag.
///
///When the flag is true, the function GetBinError
/// will approximate the bin error with the average profile error on all bins
/// in the following situation only
///
///  - the number of bins in the profile is less than 1002
///  - the bin number of entries is small ( <5)
///  - the estimated bin error is extremely small compared to the bin content
///    (see TProfile::GetBinError)

void TProfile::Approximate(Bool_t approx)
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

Int_t TProfile::BufferEmpty(Int_t action)
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
   if (CanExtendAllAxes() || fXaxis.GetXmax() <= fXaxis.GetXmin()) {
      //find min, max of entries in buffer
      Double_t xmin = fBuffer[2];
      Double_t xmax = xmin;
      for (Int_t i=1;i<nbentries;i++) {
         Double_t x = fBuffer[3*i+2];
         if (x < xmin) xmin = x;
         if (x > xmax) xmax = x;
      }
      if (fXaxis.GetXmax() <= fXaxis.GetXmin()) {
         THLimitsFinder::GetLimitsFinder()->FindGoodLimits(this,xmin,xmax);
      } else {
         fBuffer = 0;
         Int_t keep = fBufferSize; fBufferSize = 0;
         if (xmin <  fXaxis.GetXmin()) ExtendAxis(xmin,&fXaxis);
         if (xmax >= fXaxis.GetXmax()) ExtendAxis(xmax,&fXaxis);
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
/// accumulate arguments in buffer. When buffer is full, empty the buffer.
///
///  - fBuffer[0] = number of entries in buffer
///  - fBuffer[1] = w of first entry
///  - fBuffer[2] = x of first entry
///  - fBuffer[3] = y of first entry

Int_t TProfile::BufferFill(Double_t x, Double_t y, Double_t w)
{
   if (!fBuffer) return -2;
   Int_t nbentries = (Int_t)fBuffer[0];
   if (nbentries < 0) {
      nbentries  = -nbentries;
      fBuffer[0] =  nbentries;
      if (fEntries > 0) {
         Double_t *buffer = fBuffer; fBuffer=0;
         Reset("ICES");  // reset without deleting the functions
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
   return -2;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a Profile histogram to a new profile histogram.

void TProfile::Copy(TObject &obj) const
{
   try {
      TProfile & pobj = dynamic_cast<TProfile&>(obj);
      TH1D::Copy(pobj);
      fBinEntries.Copy(pobj.fBinEntries);
      fBinSumw2.Copy(pobj.fBinSumw2);
      for (int bin=0;bin<fNcells;bin++) {
         pobj.fArray[bin]        = fArray[bin];
         pobj.fSumw2.fArray[bin] = fSumw2.fArray[bin];
      }

      pobj.fYmin = fYmin;
      pobj.fYmax = fYmax;
      pobj.fScaling   = fScaling;
      pobj.fErrorMode = fErrorMode;
      pobj.fTsumwy    = fTsumwy;
      pobj.fTsumwy2   = fTsumwy2;

   } catch(...) {
      Fatal("Copy","Cannot copy a TProfile in a %s",obj.IsA()->GetName());
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: `this = this/(c1*f1)`.
///
/// This function is not implemented for the TProfile

Bool_t TProfile::Divide(TF1 *, Double_t )
{
   Error("Divide","Function not implemented for TProfile");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this profile by h1.
///
/// `this = this/h1`
///
/// This function accepts to divide a TProfile by a histogram
///
/// The function return kFALSE if the divide operation failed

Bool_t TProfile::Divide(const TH1 *h1)
{
   if (!h1) {
      Error("Divide","Attempt to divide a non-existing profile");
      return kFALSE;
   }
   if (!h1->InheritsFrom(TH1::Class())) {
      Error("Divide","Attempt to divide by a non-profile or non-histogram object");
      return kFALSE;
   }
   TProfile *p1 = (TProfile*)h1;

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);


   Int_t nbinsx = GetNbinsX();
   //- Check profile compatibility
   if (nbinsx != p1->GetNbinsX()) {
      Error("Divide","Attempt to divide profiles with different number of bins");
      return kFALSE;
   }

   //- Reset statistics
   fEntries = fTsumw   = fTsumw2 = fTsumwx = fTsumwx2 = fTsumwy = fTsumwy2 = 0;

   //- Loop on bins (including underflows/overflows)
   Int_t bin;
   Double_t *cu1=0, *er1=0, *en1=0;
   Double_t e0,e1,c12;
   if (h1->InheritsFrom(TProfile::Class())) {
      cu1 = p1->GetW();
      er1 = p1->GetW2();
      en1 = p1->GetB();
   }
   Double_t c0,c1,w,z,x;
   for (bin=0;bin<=nbinsx+1;bin++) {
      c0  = fArray[bin];
      if (cu1) c1  = cu1[bin];
      else     c1  = h1->GetBinContent(bin);
      if (c1) w = c0/c1;
      else    w = 0;
      fArray[bin] = w;
      z = TMath::Abs(w);
      x = fXaxis.GetBinCenter(bin);
      fEntries++;
      fTsumw   += z;
      fTsumw2  += z*z;
      fTsumwx  += z*x;
      fTsumwx2 += z*x*x;
      fTsumwy  += z*c1;
      fTsumwx2 += z*c1*c1;
      e0 = fSumw2.fArray[bin];
      if (er1) e1 = er1[bin];
      else    {e1 = h1->GetBinError(bin); e1*=e1;}
      c12= c1*c1;
      if (!c1) fSumw2.fArray[bin] = 0;
      else     fSumw2.fArray[bin] = (e0*c1*c1 + e1*c0*c0)/(c12*c12);
      if (!en1) continue;
      if (!en1[bin]) fBinEntries.fArray[bin] = 0;
      else           fBinEntries.fArray[bin] /= en1[bin];
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
/// Replace contents of this profile by the division of h1 by h2.
///
/// `this = c1*h1/(c2*h2)`
///
/// The function return kFALSE if the divide operation failed

Bool_t TProfile::Divide(const TH1 *h1, const TH1 *h2, Double_t c1, Double_t c2, Option_t *option)
{
   TString opt = option;
   opt.ToLower();
   Bool_t binomial = kFALSE;
   if (opt.Contains("b")) binomial = kTRUE;
   if (!h1 || !h2) {
      Error("Divide","Attempt to divide a non-existing profile");
      return kFALSE;
   }
   if (!h1->InheritsFrom(TProfile::Class())) {
      Error("Divide","Attempt to divide a non-profile object");
      return kFALSE;
   }
   TProfile *p1 = (TProfile*)h1;
   if (!h2->InheritsFrom(TProfile::Class())) {
      Error("Divide","Attempt to divide by a non-profile object");
      return kFALSE;
   }
   TProfile *p2 = (TProfile*)h2;

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);

   Int_t nbinsx = GetNbinsX();
   //- Check histogram compatibility
   if (nbinsx != p1->GetNbinsX() || nbinsx != p2->GetNbinsX()) {
      Error("Divide","Attempt to divide profiles with different number of bins");
      return kFALSE;
   }
   if (!c2) {
      Error("Divide","Coefficient of dividing profile cannot be zero");
      return kFALSE;
   }

   //THE ALGORITHM COMPUTING THE ERRORS IS WRONG. HELP REQUIRED
   printf("WARNING!!: The algorithm in TProfile::Divide computing the errors is not accurate\n");
   printf(" Instead of Divide(TProfile *h1, TProfile *h2), do:\n");
   printf("   TH1D *p1 = h1->ProjectionX();\n");
   printf("   TH1D *p2 = h2->ProjectionX();\n");
   printf("   p1->Divide(p2);\n");

   //- Reset statistics
   fEntries = fTsumw   = fTsumw2 = fTsumwx = fTsumwx2 = 0;

   //- Loop on bins (including underflows/overflows)
   Int_t bin;
   Double_t *cu1 = p1->GetW();
   Double_t *cu2 = p2->GetW();
   Double_t *er1 = p1->GetW2();
   Double_t *er2 = p2->GetW2();
   Double_t *en1 = p1->GetB();
   Double_t *en2 = p2->GetB();
   Double_t b1,b2,w,z,x,ac1,ac2;
   //d1 = c1*c1;
   //d2 = c2*c2;
   ac1 = TMath::Abs(c1);
   ac2 = TMath::Abs(c2);
   for (bin=0;bin<=nbinsx+1;bin++) {
      b1  = cu1[bin];
      b2  = cu2[bin];
      if (b2) w = c1*b1/(c2*b2);
      else    w = 0;
      fArray[bin] = w;
      z = TMath::Abs(w);
      x = fXaxis.GetBinCenter(bin);
      fEntries++;
      fTsumw   += z;
      fTsumw2  += z*z;
      fTsumwx  += z*x;
      fTsumwx2 += z*x*x;
      //fTsumwy  += z*x;
      //fTsumwy2 += z*x*x;
      Double_t e1 = er1[bin];
      Double_t e2 = er2[bin];
      //Double_t b22= b2*b2*d2;
      Double_t b22= b2*b2*TMath::Abs(c2);
      if (!b2) fSumw2.fArray[bin] = 0;
      else {
         if (binomial) {
            fSumw2.fArray[bin] = TMath::Abs(w*(1-w)/b2);
         } else {
            //fSumw2.fArray[bin] = d1*d2*(e1*b2*b2 + e2*b1*b1)/(b22*b22);
            fSumw2.fArray[bin] = ac1*ac2*(e1*b2*b2 + e2*b1*b1)/(b22*b22);
         }
      }
      if (en2[bin]) fBinEntries.fArray[bin] = en1[bin]/en2[bin];
      else          fBinEntries.fArray[bin] = 0;
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
/// Fill a Profile histogram (no weights).

Int_t TProfile::Fill(Double_t x, Double_t y)
{
   if (fBuffer) return BufferFill(x,y,1);

   Int_t bin;
   if (fYmin != fYmax) {
      if (y <fYmin || y> fYmax || TMath::IsNaN(y) ) return -1;
   }

   fEntries++;
   bin =fXaxis.FindBin(x);
   AddBinContent(bin, y);
   fSumw2.fArray[bin] += (Double_t)y*y;
   fBinEntries.fArray[bin] += 1;
   if (fBinSumw2.fN)  fBinSumw2.fArray[bin] += 1;
   if (bin == 0 || bin > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   fTsumw++;
   fTsumw2++;
   fTsumwx  += x;
   fTsumwx2 += x*x;
   fTsumwy  += y;
   fTsumwy2 += y*y;
   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a Profile histogram (no weights).

Int_t TProfile::Fill(const char *namex, Double_t y)
{
   Int_t bin;
   if (fYmin != fYmax) {
      if (y <fYmin || y> fYmax || TMath::IsNaN(y) ) return -1;
   }

   fEntries++;
   bin =fXaxis.FindBin(namex);
   AddBinContent(bin, y);
   fSumw2.fArray[bin] += (Double_t)y*y;
   fBinEntries.fArray[bin] += 1;
   if (fBinSumw2.fN)  fBinSumw2.fArray[bin] += 1;
   if (bin == 0 || bin > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   fTsumw++;
   fTsumw2++;
   fTsumwy += y;
   fTsumwy2 += y * y;
   if (!fXaxis.CanExtend() || !fXaxis.IsAlphanumeric()) {
      Double_t x = fXaxis.GetBinCenter(bin);
      fTsumwx += x;
      fTsumwx2 += x * x;
   }
   return bin;
}
////////////////////////////////////////////////////////////////////////////////
/// Fill a Profile histogram with weights.

Int_t TProfile::Fill(Double_t x, Double_t y, Double_t w)
{
   if (fBuffer) return BufferFill(x,y,w);

   Int_t bin;
   if (fYmin != fYmax) {
      if (y <fYmin || y> fYmax || TMath::IsNaN(y) ) return -1;
   }

   Double_t u= w;
   fEntries++;
   bin =fXaxis.FindBin(x);
   AddBinContent(bin, u*y);
   fSumw2.fArray[bin] += u*y*y;
   if (!fBinSumw2.fN && u != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();  // must be called before accumulating the entries
   if (fBinSumw2.fN)  fBinSumw2.fArray[bin] += u*u;
   fBinEntries.fArray[bin] += u;
   if (bin == 0 || bin > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   fTsumw   += u;
   fTsumw2  += u*u;
   fTsumwx  += u*x;
   fTsumwx2 += u*x*x;
   fTsumwy  += u*y;
   fTsumwy2 += u*y*y;
   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a Profile histogram with weights.

Int_t TProfile::Fill(const char *namex, Double_t y, Double_t w)
{
   Int_t bin;

   if (fYmin != fYmax) {
      if (y <fYmin || y> fYmax || TMath::IsNaN(y) ) return -1;
   }

   Double_t u= w; // (w > 0 ? w : -w);
   fEntries++;
   bin =fXaxis.FindBin(namex);
   AddBinContent(bin, u*y);
   fSumw2.fArray[bin] += u*y*y;
   if (!fBinSumw2.fN && u != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();  // must be called before accumulating the entries
   if (fBinSumw2.fN)  fBinSumw2.fArray[bin] += u*u;
   fBinEntries.fArray[bin] += u;
   if (bin == 0 || bin > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   fTsumw   += u;
   fTsumw2  += u*u;
   if (!fXaxis.CanExtend() || !fXaxis.IsAlphanumeric()) {
      Double_t x = fXaxis.GetBinCenter(bin);
      fTsumwx += u*x;
      fTsumwx2 += u*x*x;
   }
   fTsumwy  += u*y;
   fTsumwy2 += u*y*y;
   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a Profile histogram with weights.

void TProfile::FillN(Int_t ntimes, const Double_t *x, const Double_t *y, const Double_t *w, Int_t stride)
{
   Int_t bin,i;
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
         ifirst = i;  // start from i
      else
         return;
   }

   for (i=ifirst;i<ntimes;i+=stride) {
      if (fYmin != fYmax) {
         if (y[i] <fYmin || y[i]> fYmax || TMath::IsNaN(y[i])) continue;
      }

      Double_t u = (w) ? w[i] : 1; // (w[i] > 0 ? w[i] : -w[i]);
      fEntries++;
      bin =fXaxis.FindBin(x[i]);
      AddBinContent(bin, u*y[i]);
      fSumw2.fArray[bin] += u*y[i]*y[i];
      if (!fBinSumw2.fN && u != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();  // must be called before accumulating the entries
      if (fBinSumw2.fN)  fBinSumw2.fArray[bin] += u*u;
      fBinEntries.fArray[bin] += u;
      if (bin == 0 || bin > fXaxis.GetNbins()) {
         if (!GetStatOverflowsBehaviour()) continue;
      }
      fTsumw   += u;
      fTsumw2  += u*u;
      fTsumwx  += u*x[i];
      fTsumwx2 += u*x[i]*x[i];
      fTsumwy  += u*y[i];
      fTsumwy2 += u*y[i]*y[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return bin content of a Profile histogram.

Double_t TProfile::GetBinContent(Int_t bin) const
{
   if (fBuffer) ((TProfile*)this)->BufferEmpty();

   if (bin < 0 || bin >= fNcells) return 0;
   if (fBinEntries.fArray[bin] == 0) return 0;
   if (!fArray) return 0;
   return fArray[bin]/fBinEntries.fArray[bin];
}

////////////////////////////////////////////////////////////////////////////////
/// Return bin entries of a Profile histogram.

Double_t TProfile::GetBinEntries(Int_t bin) const
{
   if (fBuffer) ((TProfile*)this)->BufferEmpty();

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

Double_t TProfile::GetBinEffectiveEntries(Int_t bin) const
{
   return TProfileHelper::GetBinEffectiveEntries((TProfile*)this, bin);
}

////////////////////////////////////////////////////////////////////////////////
/// Return bin error of a Profile histogram
///
/// Computing errors: A moving field
///
/// The computation of errors for a TProfile has evolved with the versions
/// of ROOT. The difficulty is in computing errors for bins with low statistics.
///
///  - prior to version 3.00, we had no special treatment of low statistic bins.
///   As a result, these bins had huge errors. The reason is that the
///   expression eprim2 is very close to 0 (rounding problems) or 0.
///  - in version 3.00 (18 Dec 2000), the algorithm is protected for values of
///   eprim2 very small and the bin errors set to the average bin errors, following
///   recommendations from a group of users.
///  - in version 3.01 (19 Apr 2001), it is realized that the algorithm above
///   should be applied only to low statistic bins.
///  - in version 3.02 (26 Sep 2001), the same group of users recommend instead
///   to take two times the average error on all bins for these low
///   statistics bins giving a very small value for eprim2.
///  - in version 3.04 (Nov 2002), the algorithm is modified/protected for the case
///   when a TProfile is projected (ProjectionX). The previous algorithm
///   generated a N^2 problem when projecting a TProfile with a large number of
///   bins (eg 100000).
///  - in version 3.05/06, a new static function TProfile::Approximate
///   is introduced to enable or disable (default) the approximation.
///
///  Ideas for improvements of this algorithm are welcome. No suggestions
/// received since our call for advice to roottalk in Jul 2002.
/// see for instance: http://root.cern.ch/root/roottalk/roottalk02/2916.html

Double_t TProfile::GetBinError(Int_t bin) const
{
   return TProfileHelper::GetBinError((TProfile*)this, bin);
}

////////////////////////////////////////////////////////////////////////////////
/// Return option to compute profile errors

Option_t *TProfile::GetErrorOption() const
{
   if (fErrorMode == kERRORSPREAD)  return "s";
   if (fErrorMode == kERRORSPREADI) return "i";
   if (fErrorMode == kERRORSPREADG) return "g";
   return "";
}

////////////////////////////////////////////////////////////////////////////////
/// fill the array stats from the contents of this profile.
///
/// The array stats must be correctly dimensioned in the calling program.
///
///  - stats[0] = sumw
///  - stats[1] = sumw2
///  - stats[2] = sumwx
///  - stats[3] = sumwx2
///  - stats[4] = sumwy
///  - stats[5] = sumwy2
///
/// If no axis-subrange is specified (via TAxis::SetRange), the array stats
/// is simply a copy of the statistics quantities computed at filling time.
/// If a sub-range is specified, the function recomputes these quantities
/// from the bin contents in the current axis range.

void TProfile::GetStats(Double_t *stats) const
{
   if (fBuffer) ((TProfile*)this)->BufferEmpty();

   // Loop on bins
   Int_t bin, binx;
   // identify the case of labels with extension of axis range
   // in this case the statistics in x does not make any sense
   Bool_t labelHist =  ((const_cast<TAxis&>(fXaxis)).GetLabels() && fXaxis.CanExtend() );

   if ( (fTsumw == 0 /* && fEntries > 0 */) || fXaxis.TestBit(TAxis::kAxisRange) ) {
      for (bin=0;bin<6;bin++) stats[bin] = 0;
      if (!fBinEntries.fArray) return;
      Int_t firstBinX = fXaxis.GetFirst();
      Int_t lastBinX  = fXaxis.GetLast();
      // include underflow/overflow if TH1::StatOverflows(kTRUE) in case no range is set on the axis
      if (GetStatOverflowsBehaviour() && !fXaxis.TestBit(TAxis::kAxisRange)) {
         if (firstBinX == 1) firstBinX = 0;
         if (lastBinX ==  fXaxis.GetNbins() ) lastBinX += 1;
      }
      for (binx = firstBinX; binx <= lastBinX; binx++) {
         Double_t w   = fBinEntries.fArray[binx];
         Double_t w2  = (fBinSumw2.fN ? fBinSumw2.fArray[binx] : w);
         Double_t x   = (!labelHist) ? fXaxis.GetBinCenter(binx) : 0;
         stats[0] += w;
         stats[1] += w2;
         stats[2] += w*x;
         stats[3] += w*x*x;
         stats[4] += fArray[binx];
         stats[5] += fSumw2.fArray[binx];
      }
   } else {
      if (fTsumwy == 0 && fTsumwy2 == 0) {
         //this case may happen when processing TProfiles with version <=3
         TProfile *p = (TProfile*)this; // cheating with const
         for (binx=fXaxis.GetFirst();binx<=fXaxis.GetLast();binx++) {
            p->fTsumwy  += fArray[binx];
            p->fTsumwy2 += fSumw2.fArray[binx];
         }
      }
      stats[0] = fTsumw;
      stats[1] = fTsumw2;
      stats[2] = fTsumwx;
      stats[3] = fTsumwx2;
      stats[4] = fTsumwy;
      stats[5] = fTsumwy2;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reduce the number of bins for this axis to the number of bins having a label.

void TProfile::LabelsDeflate(Option_t *option)
{
   TProfileHelper::LabelsDeflate(this, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Double the number of bins for axis.
/// Refill histogram
/// This function is called by TAxis::FindBin(const char *label)

void TProfile::LabelsInflate(Option_t *options)
{
   TProfileHelper::LabelsInflate(this, options);
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

void TProfile::LabelsOption(Option_t *option, Option_t * /*ax */)
{
   THashList *labels = fXaxis.GetLabels();
   if (!labels) {
      Warning("LabelsOption","Cannot sort. No labels");
      return;
   }
   TString opt = option;
   opt.ToLower();
   if (opt.Contains("h")) {
      fXaxis.SetBit(TAxis::kLabelsHori);
      fXaxis.ResetBit(TAxis::kLabelsVert);
      fXaxis.ResetBit(TAxis::kLabelsDown);
      fXaxis.ResetBit(TAxis::kLabelsUp);
   }
   if (opt.Contains("v")) {
      fXaxis.SetBit(TAxis::kLabelsVert);
      fXaxis.ResetBit(TAxis::kLabelsHori);
      fXaxis.ResetBit(TAxis::kLabelsDown);
      fXaxis.ResetBit(TAxis::kLabelsUp);
   }
   if (opt.Contains("u")) {
      fXaxis.SetBit(TAxis::kLabelsUp);
      fXaxis.ResetBit(TAxis::kLabelsVert);
      fXaxis.ResetBit(TAxis::kLabelsDown);
      fXaxis.ResetBit(TAxis::kLabelsHori);
   }
   if (opt.Contains("d")) {
      fXaxis.SetBit(TAxis::kLabelsDown);
      fXaxis.ResetBit(TAxis::kLabelsVert);
      fXaxis.ResetBit(TAxis::kLabelsHori);
      fXaxis.ResetBit(TAxis::kLabelsUp);
   }
   Int_t sort = -1;
   if (opt.Contains("a")) sort = 0;
   if (opt.Contains(">")) sort = 1;
   if (opt.Contains("<")) sort = 2;
   if (sort < 0) return;

   // support only cases when first n bins have labels
   Int_t n = labels->GetSize();
   TAxis *axis = &fXaxis;
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
               "%s of TProfile %s contains bins without labels. Sorting will not work correctly - return",
               axis->GetName(), GetName());
         return;
      }
      // case where label bins are consecutive starting from first bin will work
      Warning(
         "LabelsOption",
         "axis %s of TProfile %s has extra following bins without labels. Sorting will work only for first label bins",
         axis->GetName(), GetName());
   }
   std::vector<Int_t> a(n);
   Int_t i;
   std::vector<Double_t> cont(n);
   std::vector<Double_t> sumw(n);
   std::vector<Double_t> errors(n);
   std::vector<Double_t> ent(n);
   std::vector<Double_t> binsw2;
   if (fBinSumw2.fN) binsw2.resize(n);

   // delete buffer if it is there since bins will be reordered.
   if (fBuffer)
      BufferEmpty(1);

   // make a labelold list but ordered with bins
   // (re-ordered original label list)
   std::vector<TObject *> labold(n);
   for (i = 0; i < n; i++)
      labold[i] = nullptr;
   TIter nextold(labels);
   TObject *obj;
   while ((obj=nextold())) {
      Int_t bin = obj->GetUniqueID();
      R__ASSERT(bin <= n);
      labold[bin - 1] = obj;
   }
   // order now labold according to bin content

   labels->Clear();
   if (sort > 0) {
      //---sort by values of bins
      for (i=1;i<=n;i++) {
         a[i-1] = i-1;
         sumw[i-1]   = fArray[i];
         errors[i-1] = fSumw2.fArray[i];
         ent[i-1]    = fBinEntries.fArray[i];
         if (fBinSumw2.fN) binsw2[i - 1] = fBinSumw2.fArray[i];
         if (fBinEntries.fArray[i] == 0) cont[i-1] = 0;
         else cont[i-1] = fArray[i]/fBinEntries.fArray[i];
      }
      if (sort ==1)
         TMath::Sort(n,cont.data(),a.data(),kTRUE);  //sort by decreasing values
      else
         TMath::Sort(n,cont.data(),a.data(),kFALSE); //sort by increasing values
      for (i=1;i<=n;i++) {
         fArray[i] = sumw[a[i-1]];
         fSumw2.fArray[i] = errors[a[i-1]];
         fBinEntries.fArray[i] = ent[a[i-1]];
         if (fBinSumw2.fN)
            fBinSumw2.fArray[i] = binsw2[a[i-1]];
      }
      for (i=0 ;i < n; i++) {
         obj = labold[a[i]];
         labels->Add(obj);
         obj->SetUniqueID(i+1);
      }
   } else {

      //---alphabetic sort
      // sort labels using vector of strings and TMath::Sort
      // I need to array because labels order in list is not necessary that of the bins
      std::vector<std::string> vecLabels(n);
      for (i = 0; i < n; i++) {
         vecLabels[i] = labold[i]->GetName();
         a[i] = i;
         sumw[i] = fArray[i+1];
         errors[i] = fSumw2.fArray[i+1];
         ent[i] = fBinEntries.fArray[i+1];
         if (fBinSumw2.fN)
            binsw2[i] = fBinSumw2.fArray[i+1];
      }
      // sort in ascending order for strings
      TMath::Sort(n, vecLabels.data(), a.data(), kFALSE);
      // set the new labels
      for (i = 0; i < n; i++) {
         TObject *labelObj = labold[a[i]];
         labels->Add(labelObj);
         // set the corresponding bin. NB bin starts from 1
         labelObj->SetUniqueID(i + 1);
         if (gDebug)
            std::cout << "bin " << i + 1 << " setting new labels for axis " << labold.at(a[i])->GetName() << " from "
                      << a[i] << std::endl;
      }

      for (i=0; i < n; i++) {
         fArray[i+1] = sumw[a[i]];
         fSumw2.fArray[i+1] = errors[a[i]];
         fBinEntries.fArray[i+1] = ent[a[i]];
         if (fBinSumw2.fN)
            fBinSumw2.fArray[i+1] = binsw2[a[i]];
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
      // if (iaxis == 1) {
      s[2] = 0; // fTsumwx
      s[3] = 0; // fTsumwx2
      PutStats(s);
   }
}

////////////////////////////////////////////////////////////////////////////////
///Merge all histograms in the collection in this histogram.
///
/// This function computes the min/max for the x axis,
/// compute a new number of bins, if necessary,
/// add bin contents, errors and statistics.
/// If overflows are present and limits are different the function will fail.
/// The function returns the total number of entries in the result histogram
/// if the merge is successfull, -1 otherwise.
///
/// IMPORTANT remark. The axis x may have different number
/// of bins and different limits, BUT the largest bin width must be
/// a multiple of the smallest bin width and the upper limit must also
/// be a multiple of the bin width.

Long64_t TProfile::Merge(TCollection *li)
{
   return TProfileHelper::Merge(this, li);
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: this = this*c1*f1
///
/// The function return kFALSE if the Multiply operation failed

Bool_t TProfile::Multiply(TF1 *f1, Double_t c1)
{

   if (!f1) {
      Error("Multiply","Attempt to multiply by a null function");
      return kFALSE;
   }

   Int_t nbinsx = GetNbinsX();

   //- Add statistics
   Double_t xx[1], cf1, ac1 = TMath::Abs(c1);
   Double_t s1[10];
   Int_t i;
   for (i=0;i<10;i++) {s1[i] = 0;}
   PutStats(s1);

   SetMinimum();
   SetMaximum();

   //- Loop on bins (including underflows/overflows)
   Int_t bin;
   for (bin=0;bin<=nbinsx+1;bin++) {
      xx[0] = fXaxis.GetBinCenter(bin);
      if (!f1->IsInside(xx)) continue;
      TF1::RejectPoint(kFALSE);
      cf1 = f1->EvalPar(xx);
      if (TF1::RejectedPoint()) continue;
      fArray[bin]             *= c1*cf1;
      //see http://savannah.cern.ch/bugs/?func=detailitem&item_id=14851
      //fSumw2.fArray[bin]      *= c1*c1*cf1*cf1;
      fSumw2.fArray[bin]      *= ac1*cf1*cf1;
      //fBinEntries.fArray[bin] *= ac1*TMath::Abs(cf1);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply this profile by h1.
///
/// `this = this*h1`

Bool_t TProfile::Multiply(const TH1 *)
{
   Error("Multiply","Multiplication of profile histograms not implemented");
   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Replace contents of this profile by multiplication of h1 by h2.
///
/// `this = (c1*h1)*(c2*h2)`

Bool_t TProfile::Multiply(const TH1 *, const TH1 *, Double_t, Double_t, Option_t *)
{
   Error("Multiply","Multiplication of profile histograms not implemented");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Project this profile into a 1-D histogram along X
///
///  The projection is always of the type TH1D.
///
///  - if option "E" is specified the errors of the projected histogram are computed and set
///     to be equal to the errors of the profile.
///     Option "E" is defined as the default one in the header file.
///  - if option "" is specified the histogram errors are simply the sqrt of its content
///  - if option "B" is specified, the content of bin of the returned histogram
///     will be equal to the GetBinEntries(bin) of the profile,
///     otherwise (default) it will be equal to GetBinContent(bin)
///  - if option "C=E" the bin contents of the projection are set to the
///      bin errors of the profile
///  - if option "W" is specified the bin content of the projected histogram  is set to the
///      product of the bin content of the profile and the entries.
///      With this option the returned histogram will be equivalent to the one obtained by
///      filling directly a TH1D using the 2-nd value as a weight.
///      This makes sense only for profile filled with weights =1. If not, the error of the
///       projected histogram obtained with this option will not be correct.

TH1D *TProfile::ProjectionX(const char *name, Option_t *option) const
{

   TString opt = option;
   opt.ToLower();
   Int_t nx = fXaxis.GetNbins();

   // Create the projection histogram
   TString pname = name;
   if (pname == "_px") {
      pname = GetName();
      pname.Append("_px");
   }
   TH1D *h1;
   const TArrayD *bins = fXaxis.GetXbins();
   if (bins->fN == 0) {
      h1 = new TH1D(pname,GetTitle(),nx,fXaxis.GetXmin(),fXaxis.GetXmax());
   } else {
      h1 = new TH1D(pname,GetTitle(),nx,bins->fArray);
   }
   Bool_t computeErrors = kFALSE;
   Bool_t cequalErrors  = kFALSE;
   Bool_t binEntries    = kFALSE;
   Bool_t binWeight     = kFALSE;
   if (opt.Contains("b")) binEntries = kTRUE;
   if (opt.Contains("e")) computeErrors = kTRUE;
   if (opt.Contains("w")) binWeight = kTRUE;
   if (opt.Contains("c=e")) {cequalErrors = kTRUE; computeErrors=kFALSE;}
   if (computeErrors || binWeight || (binEntries && fBinSumw2.fN) ) h1->Sumw2();

   // Fill the projected histogram
   Double_t cont;
   for (Int_t bin =0;bin<=nx+1;bin++) {

      if (binEntries)         cont = GetBinEntries(bin);
      else if (cequalErrors)  cont = GetBinError(bin);
      else if (binWeight)     cont = fArray[bin];  // bin content * bin entries
      else                    cont = GetBinContent(bin);    // default case

      h1->SetBinContent(bin ,cont);

      // if option E projected histogram errors are same as profile
      if (computeErrors ) h1->SetBinError(bin , GetBinError(bin) );
      // in case of option W bin error is deduced from bin sum of z**2 values of profile
      // this is correct only if the profile is filled with weights =1
      if (binWeight) h1->GetSumw2()->fArray[bin] = fSumw2.fArray[bin];
      // in case of bin entries and profile is weighted, we need to set also the bin error
      if (binEntries && fBinSumw2.fN ) {
         R__ASSERT(  h1->GetSumw2() );
         h1->GetSumw2()->fArray[bin] =  fBinSumw2.fArray[bin];
      }

   }

   // Copy the axis attributes and the axis labels if needed.
   h1->GetXaxis()->ImportAttributes(this->GetXaxis());
   h1->GetYaxis()->ImportAttributes(this->GetYaxis());
   THashList* labels=this->GetXaxis()->GetLabels();
   if (labels) {
      TIter iL(labels);
      TObjString* lb;
      Int_t i = 1;
      while ((lb=(TObjString*)iL())) {
         h1->GetXaxis()->SetBinLabel(i,lb->String().Data());
         i++;
      }
   }

   h1->SetEntries(fEntries);
   return h1;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace current statistics with the values in array stats.

void TProfile::PutStats(Double_t *stats)
{
   fTsumw   = stats[0];
   fTsumw2  = stats[1];
   fTsumwx  = stats[2];
   fTsumwx2 = stats[3];
   fTsumwy  = stats[4];
   fTsumwy2 = stats[5];
}

////////////////////////////////////////////////////////////////////////////////
/// Rebin this profile grouping ngroup bins together.
///
/// ## case 1  xbins=0
///  if newname is not blank a new temporary profile hnew is created.
///  else the current profile is modified (default)
///  The parameter ngroup indicates how many bins of this have to me merged
///  into one bin of hnew
///  If the original profile has errors stored (via Sumw2), the resulting
///  profile has new errors correctly calculated.
///
///  examples: if hp is an existing TProfile histogram with 100 bins
///
/// ~~~ {.cpp}
///    hp->Rebin();  //merges two bins in one in hp: previous contents of hp are lost
///    hp->Rebin(5); //merges five bins in one in hp
///    TProfile *hnew = hp->Rebin(5,"hnew"); // creates a new profile hnew
///                                      //merging 5 bins of hp in one bin
/// ~~~
///
///  NOTE:  If ngroup is not an exact divider of the number of bins,
///         the top limit of the rebinned profile is changed
///         to the upper edge of the bin=newbins*ngroup and the corresponding
///         bins are added to the overflow bin.
///         Statistics will be recomputed from the new bin contents.
///
///  ## case 2  xbins!=0
///  a new profile is created (you should specify newname).
///  The parameter ngroup is the number of variable size bins in the created profile
///  The array xbins must contain ngroup+1 elements that represent the low-edge
///  of the bins.
///  The data of the old bins are added to the new bin which contains the bin center
///  of the old bins. It is possible that information from the old binning are attached
///  to the under-/overflow bins of the new binning.
///
///  examples: if hp is an existing TProfile with 100 bins
///
/// ~~~ {.cpp}
///      Double_t xbins[25] = {...} array of low-edges (xbins[25] is the upper edge of last bin
///      hp->Rebin(24,"hpnew",xbins);  //creates a new variable bin size profile hpnew
/// ~~~

TH1 *TProfile::Rebin(Int_t ngroup, const char*newname, const Double_t *xbins)
{
   Int_t nbins    = fXaxis.GetNbins();
   Double_t xmin  = fXaxis.GetXmin();
   Double_t xmax  = fXaxis.GetXmax();
   if ((ngroup <= 0) || (ngroup > nbins)) {
      Error("Rebin", "Illegal value of ngroup=%d",ngroup);
      return 0;
   }
   if (!newname && xbins) {
      Error("Rebin","if xbins is specified, newname must be given");
      return 0;
   }

   Int_t newbins = nbins/ngroup;
   if (!xbins) {
      Int_t nbg = nbins/ngroup;
      if (nbg*ngroup != nbins) {
         Warning("Rebin", "ngroup=%d must be an exact divider of nbins=%d",ngroup,nbins);
      }
   }
   else {
      // in the case of xbins given (rebinning in variable bins) ngroup is the new number of bins.
      // and number of grouped bins is not constant.
      // when looping for setting the contents for the new histogram we
      // need to loop on all bins of original histogram. Set then ngroup=nbins
      newbins = ngroup;
      ngroup = nbins;
   }

   // Save old bin contents into a new array
   Double_t *oldBins   = new Double_t[nbins+2];
   Double_t *oldCount  = new Double_t[nbins+2];
   Double_t *oldErrors = new Double_t[nbins+2];
   Double_t *oldBinw2  = (fBinSumw2.fN ? new Double_t[nbins+2] : 0  );
   Int_t bin, i;
   Double_t *cu1 = GetW();
   Double_t *er1 = GetW2();
   Double_t *en1 = GetB();
   Double_t *ew1 = GetB2();

   for (bin=0;bin<=nbins+1;bin++) {
      oldBins[bin]   = cu1[bin];
      oldCount[bin]  = en1[bin];
      oldErrors[bin] = er1[bin];
      if (ew1 && fBinSumw2.fN) oldBinw2[bin]  = ew1[bin];
   }

   // create a clone of the old histogram if newname is specified
   TProfile *hnew = this;
   if ((newname && strlen(newname) > 0) || xbins) {
      hnew = (TProfile*)Clone(newname);
   }

   // in case of ngroup not an excat divider of nbins,
   // top limit is changed (see NOTE in method comment)
   if(!xbins && (newbins*ngroup != nbins)) {
      xmax = fXaxis.GetBinUpEdge(newbins*ngroup);
      hnew->fTsumw = 0; //stats must be reset because top bins will be moved to overflow bin
   }

   // set correctly the axis and resizes the bin arrays
   if(!xbins && (fXaxis.GetXbins()->GetSize() > 0)){
      // for rebinning of variable bins in a constant group
      Double_t *bins = new Double_t[newbins+1];
      for(i = 0; i <= newbins; ++i) bins[i] = fXaxis.GetBinLowEdge(1+i*ngroup);
      hnew->SetBins(newbins,bins); //this also changes the bin array's
      delete [] bins;
   } else if (xbins) {
      // when rebinning in variable bins
      hnew->SetBins(newbins,xbins);
   } else {
      hnew->SetBins(newbins,xmin,xmax);
   }

   // merge bin contents ignoring now underflow/overflows
   if (fBinSumw2.fN) hnew->Sumw2();

   // Start merging only once the new lowest edge is reached
   Int_t startbin = 1;
   const Double_t newxmin = hnew->GetXaxis()->GetBinLowEdge(1);
   while( fXaxis.GetBinCenter(startbin) < newxmin && startbin <= nbins ) {
      startbin++;
   }

   Double_t *cu2 = hnew->GetW();
   Double_t *er2 = hnew->GetW2();
   Double_t *en2 = hnew->GetB();
   Double_t *ew2 = hnew->GetB2();
   Int_t oldbin = startbin;
   Double_t binContent, binCount, binError, binSumw2;
   for (bin = 1;bin<=newbins;bin++) {
      binContent = 0;
      binCount   = 0;
      binError   = 0;
      binSumw2   = 0;

      //for xbins != 0: ngroup == nbins
      Int_t imax = ngroup;
      Double_t xbinmax = hnew->GetXaxis()->GetBinUpEdge(bin);
      for (i=0;i<ngroup;i++) {
         if((hnew == this && (oldbin+i > nbins)) ||
            (hnew != this && (fXaxis.GetBinCenter(oldbin+i) > xbinmax)))
         {
            imax = i;
            break;
         }

         binContent += oldBins[oldbin+i];
         binCount   += oldCount[oldbin+i];
         binError   += oldErrors[oldbin+i];
         if (fBinSumw2.fN) binSumw2 += oldBinw2[oldbin+i];
      }

      cu2[bin] = binContent;
      er2[bin] = binError;
      en2[bin] = binCount;
      if (fBinSumw2.fN) ew2[bin] = binSumw2;
      oldbin += imax;
   }
   // set bin statistics for underflow bin
   binContent = 0;
   binCount   = 0;
   binError   = 0;
   binSumw2   = 0;
   for(i=0;i<startbin;i++)
   {
      binContent += oldBins[i];
      binCount   += oldCount[i];
      binError   += oldErrors[i];
      if (fBinSumw2.fN) binSumw2 += oldBinw2[i];
   }
   hnew->fArray[0] = binContent;
   hnew->fBinEntries[0] = binCount;
   hnew->fSumw2[0] = binError;
   if ( fBinSumw2.fN ) hnew->fBinSumw2[0] = binSumw2;

   // set bin statistics for overflow bin
   binContent = 0;
   binCount   = 0;
   binError   = 0;
   binSumw2   = 0;
   for(i=oldbin;i<=nbins+1;i++)
   {
      binContent += oldBins[i];
      binCount   += oldCount[i];
      binError   += oldErrors[i];
      if (fBinSumw2.fN) binSumw2 += oldBinw2[i];
   }
   hnew->fArray[newbins+1] = binContent;
   hnew->fBinEntries[newbins+1] = binCount;
   hnew->fSumw2[newbins+1] = binError;
   if ( fBinSumw2.fN ) hnew->fBinSumw2[newbins+1] = binSumw2;


   delete [] oldBins;
   delete [] oldCount;
   delete [] oldErrors;
   if (oldBinw2) delete [] oldBinw2;
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Profile histogram is resized along x axis such that x is in the axis range.
/// The new axis limits are recomputed by doubling iteratively
/// the current axis range until the specified value x is within the limits.
/// The algorithm makes a copy of the histogram, then loops on all bins
/// of the old histogram to fill the extended histogram.
/// Takes into account errors (Sumw2) if any.
/// The axis must be extendable before invoking this function.
///
/// Ex: `h->GetXaxis()->SetCanExtend(kTRUE)`

void TProfile::ExtendAxis(Double_t x, TAxis *axis)
{
   TProfile*  hold = TProfileHelper::ExtendAxis(this, x, axis);
   if ( hold ) {
      fTsumwy  = hold->fTsumwy;
      fTsumwy2 = hold->fTsumwy2;

      delete hold;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset contents of a Profile histogram.

void TProfile::Reset(Option_t *option)
{
   TH1D::Reset(option);
   fBinEntries.Reset();
   fBinSumw2.Reset();
   TString opt = option;
   opt.ToUpper();
   if (opt.Contains("ICE") && !opt.Contains("S")) return;
   fTsumwy  = 0;
   fTsumwy2 = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out.

void TProfile::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   //Note the following restrictions in the code generated:
   // - variable bin size not implemented
   // - SetErrorOption not implemented

   Bool_t nonEqiX = kFALSE;
   Int_t i;
   // Check if the profile has equidistant X bins or not.  If not, we
   // create an array holding the bins.
   if (GetXaxis()->GetXbins()->fN && GetXaxis()->GetXbins()->fArray) {
      nonEqiX = kTRUE;
      out << "   Double_t xAxis[" << GetXaxis()->GetXbins()->fN
      << "] = {";
      for (i = 0; i < GetXaxis()->GetXbins()->fN; i++) {
         if (i != 0) out << ", ";
         out << GetXaxis()->GetXbins()->fArray[i];
      }
      out << "}; " << std::endl;
   }

   char quote = '"';
   out<<"   "<<std::endl;
   out<<"   "<<ClassName()<<" *";

   //histogram pointer has by default teh histogram name.
   //however, in case histogram has no directory, it is safer to add a incremental suffix
   static Int_t hcounter = 0;
   TString histName = GetName();
   if (!fDirectory) {
      hcounter++;
      histName += "__";
      histName += hcounter;
   }
   const char *hname = histName.Data();

   out<<hname<<" = new "<<ClassName()<<"("<<quote<<GetName()<<quote<<","<<quote<<GetTitle()<<quote
   <<","<<GetXaxis()->GetNbins();
   if (nonEqiX)
      out << ", xAxis";
   else
      out << "," << GetXaxis()->GetXmin()
          << "," << GetXaxis()->GetXmax()
          <<","<<quote<<GetErrorOption()<<quote<<");"<<std::endl;

   // save bin entries
   Int_t bin;
   for (bin=0;bin<fNcells;bin++) {
      Double_t bi = GetBinEntries(bin);
      if (bi) {
         out<<"   "<<hname<<"->SetBinEntries("<<bin<<","<<bi<<");"<<std::endl;
      }
   }
   //save bin contents
   for (bin=0;bin<fNcells;bin++) {
      Double_t bc = fArray[bin];
      if (bc) {
         out<<"   "<<hname<<"->SetBinContent("<<bin<<","<<bc<<");"<<std::endl;
      }
   }
   // save bin errors
   if (fSumw2.fN) {
      for (bin=0;bin<fNcells;bin++) {
         Double_t be = TMath::Sqrt(fSumw2.fArray[bin]);
         if (be) {
            out<<"   "<<hname<<"->SetBinError("<<bin<<","<<be<<");"<<std::endl;
         }
      }
   }

   TH1::SavePrimitiveHelp(out, hname, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply this profile by a constant c1.
///
/// `this = c1*this`
///
/// This function uses the services of TProfile::Add

void TProfile::Scale(Double_t c1, Option_t * option)
{
   TProfileHelper::Scale(this, c1, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the number of entries in bin.

void TProfile::SetBinEntries(Int_t bin, Double_t w)
{
   TProfileHelper::SetBinEntries(this, bin, w);
}

////////////////////////////////////////////////////////////////////////////////
/// Redefine x axis parameters.

void TProfile::SetBins(Int_t nx, Double_t xmin, Double_t xmax)
{
   fXaxis.Set(nx,xmin,xmax);
   fNcells = nx+2;
   SetBinsLength(fNcells);
}

////////////////////////////////////////////////////////////////////////////////
/// Redefine  x axis parameters.

void TProfile::SetBins(Int_t nx, const Double_t *xbins)
{
   fXaxis.Set(nx,xbins);
   fNcells = nx+2;
   SetBinsLength(fNcells);
}

////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow.
/// Reallocate bin contents array

void TProfile::SetBinsLength(Int_t n)
{
   TH1D::SetBinsLength(n);
   TProfileHelper::BuildArray(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the buffer size in units of 8 bytes (double).

void TProfile::SetBuffer(Int_t buffersize, Option_t *)
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
   fBufferSize = 1 + 3*buffersize;
   fBuffer = new Double_t[fBufferSize];
   memset(fBuffer,0,sizeof(Double_t)*fBufferSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Set option to compute profile errors.
///
/// The computation of the bin errors is based on the parameter option:
///
///    -' '  (Default) The bin errors are the standard error on the mean of the bin profiled values (Y),
///                   i.e. the standard error of the bin contents.
///                   Note that if TProfile::Approximate()  is called, an approximation is used when
///                   the spread in Y is 0 and the number of bin entries  is > 0
///    -'s'            The bin errors are the standard deviations of the Y bin values
///                   Note that if TProfile::Approximate()  is called, an approximation is used when
///                   the spread in Y is 0 and the number of bin entries is > 0
///    -'i'            Errors are as in default case (standard errors of the bin contents)
///                   The only difference is for the case when the spread in Y is zero.
///                   In this case for N > 0 the error is  1./SQRT(12.*N)
///    -'g'            Errors are 1./SQRT(W)  for W not equal to 0 and 0 for W = 0.
///                   W is the sum in the bin of the weights of the profile.
///                   This option is for combining measurements y +/- dy,
///                   and  the profile is filled with values y and weights w = 1/dy**2
///
///  See TProfile::BuildOptions for a detailed explanation of all options

void TProfile::SetErrorOption(Option_t *option)
{
   TProfileHelper::SetErrorOption(this, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TProfile.

void TProfile::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TProfile::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TH1D::Streamer(R__b);
      fBinEntries.Streamer(R__b);
      Int_t errorMode;
      R__b >> errorMode;
      fErrorMode = (EErrorType)errorMode;
      if (R__v < 2) {
         Float_t ymin,ymax;
         R__b >> ymin; fYmin = ymin;
         R__b >> ymax; fYmax = ymax;
      } else {
         R__b >> fYmin;
         R__b >> fYmax;
      }
      R__b.CheckByteCount(R__s, R__c, TProfile::IsA());
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TProfile::Class(),this);
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Create/delete structure to store sum of squares of weights per bin.
///
/// This is needed to compute  the correct statistical quantities
/// of a profile filled with weights
///
/// This function is automatically called when the histogram is created
/// if the static function TH1::SetDefaultSumw2 has been called before.
/// If flag is false the structure is deleted

void TProfile::Sumw2(Bool_t flag)
{
   TProfileHelper::Sumw2(this, flag);
}
