// @(#)root/hist:$Id$
// Author: Federico Carminati   28/02/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TSpline
    \ingroup Hist
 Base class for spline implementation containing the Draw/Paint methods.
*/

#include "TROOT.h"
#include "TGraph.h"
#include "TBuffer.h"
#include "TSpline.h"
#include "TVirtualPad.h"
#include "TH1.h"
#include "TF1.h"
#include "TSystem.h"
#include "TMath.h"
#include "strlcpy.h"
#include "snprintf.h"
#include <iostream>
#include <fstream>

ClassImp(TSplinePoly);
ClassImp(TSplinePoly3);
ClassImp(TSplinePoly5);
ClassImp(TSpline3);
ClassImp(TSpline5);
ClassImp(TSpline);

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TSpline::TSpline(const TSpline &sp) :
  TNamed(sp),
  TAttLine(sp),
  TAttFill(sp),
  TAttMarker(sp),
  fDelta(sp.fDelta),
  fXmin(sp.fXmin),
  fXmax(sp.fXmax),
  fNp(sp.fNp),
  fKstep(sp.fKstep),
  fHistogram(0),
  fGraph(0),
  fNpx(sp.fNpx)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TSpline::~TSpline()
{
   if(fHistogram) delete fHistogram;
   if(fGraph) delete fGraph;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TSpline& TSpline::operator=(const TSpline &sp)
{
   if(this!=&sp) {
      TNamed::operator=(sp);
      TAttLine::operator=(sp);
      TAttFill::operator=(sp);
      TAttMarker::operator=(sp);
      fDelta=sp.fDelta;
      fXmin=sp.fXmin;
      fXmax=sp.fXmax;
      fNp=sp.fNp;
      fKstep=sp.fKstep;
      fHistogram=0;
      fGraph=0;
      fNpx=sp.fNpx;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this function with its current attributes.
///
/// Possible option values are:
///
///  - "SAME"  superimpose on top of existing picture
///  - "L"     connect all computed points with a straight line
///  - "C"     connect all computed points with a smooth curve.
///  - "P"     add a polymarker at each knot
///
/// Note that the default value is "L". Therefore to draw on top
/// of an existing picture, specify option "LSAME"

void TSpline::Draw(Option_t *option)
{
   TString opt = option;
   opt.ToLower();
   if (gPad && !opt.Contains("same")) gPad->Clear();

   AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a spline.

Int_t TSpline::DistancetoPrimitive(Int_t px, Int_t py)
{
   if (!fHistogram) return 999;
   return fHistogram->DistancetoPrimitive(px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.

void TSpline::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (!fHistogram) return;
   fHistogram->ExecuteEvent(event, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this function with its current attributes.

void TSpline::Paint(Option_t *option)
{
   Int_t i;
   Double_t xv;

   TString opt = option;
   opt.ToLower();
   Double_t xmin, xmax, pmin, pmax;
   pmin = gPad->PadtoX(gPad->GetUxmin());
   pmax = gPad->PadtoX(gPad->GetUxmax());
   xmin = fXmin;
   xmax = fXmax;
   if (opt.Contains("same")) {
      if (xmax < pmin) return;  // Otto: completely outside
      if (xmin > pmax) return;
      if (xmin < pmin) xmin = pmin;
      if (xmax > pmax) xmax = pmax;
   } else {
      gPad->Clear();
   }

   //  Create a temporary histogram and fill each channel with the function value
   if (fHistogram)
      if ((!gPad->GetLogx() && fHistogram->TestBit(TH1::kLogX)) ||
          (gPad->GetLogx() && !fHistogram->TestBit(TH1::kLogX)))
         { delete fHistogram; fHistogram = 0;}

   if (fHistogram) {
      //if (xmin != fXmin || xmax != fXmax)
      fHistogram->GetXaxis()->SetLimits(xmin,xmax);
   } else {
   //      if logx, we must bin in logx and not in x !!!
   //      otherwise if several decades, one gets crazy results
      if (xmin > 0 && gPad->GetLogx()) {
         Double_t *xbins  = new Double_t[fNpx+1];
         Double_t xlogmin = TMath::Log10(xmin);
         Double_t xlogmax = TMath::Log10(xmax);
         Double_t dlogx   = (xlogmax-xlogmin)/((Double_t)fNpx);
         for (i=0;i<=fNpx;i++) {
            xbins[i] = gPad->PadtoX(xlogmin+ i*dlogx);
         }
         fHistogram = new TH1F("Spline",GetTitle(),fNpx,xbins);
         fHistogram->SetBit(TH1::kLogX);
         delete [] xbins;
      } else {
         fHistogram = new TH1F("Spline",GetTitle(),fNpx,xmin,xmax);
      }
      if (!fHistogram) return;
      fHistogram->SetDirectory(0);
   }
   for (i=1;i<=fNpx;i++) {
      xv = fHistogram->GetBinCenter(i);
      fHistogram->SetBinContent(i,this->Eval(xv));
   }

   // Copy Function attributes to histogram attributes
   fHistogram->SetBit(TH1::kNoStats);
   fHistogram->SetLineColor(GetLineColor());
   fHistogram->SetLineStyle(GetLineStyle());
   fHistogram->SetLineWidth(GetLineWidth());
   fHistogram->SetFillColor(GetFillColor());
   fHistogram->SetFillStyle(GetFillStyle());
   fHistogram->SetMarkerColor(GetMarkerColor());
   fHistogram->SetMarkerStyle(GetMarkerStyle());
   fHistogram->SetMarkerSize(GetMarkerSize());

   //  Draw the histogram
   //  but first strip off the 'p' option if any
   char *o = (char *) opt.Data();
   Int_t j=0;
   i=0;
   Bool_t graph=kFALSE;
   do
      if(o[i]=='p') graph=kTRUE ; else o[j++]=o[i];
   while(o[i++]);
   if (opt.Length() == 0 ) fHistogram->Paint("lf");
   else if (opt == "same") fHistogram->Paint("lfsame");
   else                    fHistogram->Paint(opt.Data());

   // Think about the graph, if demanded
   if(graph) {
      if(!fGraph) {
         Double_t *xx = new Double_t[fNp];
         Double_t *yy = new Double_t[fNp];
         for(i=0; i<fNp; ++i)
            GetKnot(i,xx[i],yy[i]);
         fGraph=new TGraph(fNp,xx,yy);
         delete [] xx;
         delete [] yy;
      }
      fGraph->SetMarkerColor(GetMarkerColor());
      fGraph->SetMarkerStyle(GetMarkerStyle());
      fGraph->SetMarkerSize(GetMarkerSize());
      fGraph->Paint("p");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TSpline.

void TSpline::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TSpline::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TNamed::Streamer(R__b);
      TAttLine::Streamer(R__b);
      TAttFill::Streamer(R__b);
      TAttMarker::Streamer(R__b);

      fNp = 0;
      /*
      R__b >> fDelta;
      R__b >> fXmin;
      R__b >> fXmax;
      R__b >> fNp;
      R__b >> fKstep;
      R__b >> fHistogram;
      R__b >> fGraph;
      R__b >> fNpx;
      */
      R__b.CheckByteCount(R__s, R__c, TSpline::IsA());
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TSpline::Class(),this);
   }
}

/** \class TSplinePoly
    \ingroup Hist
 Base class for TSpline knot.
*/

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TSplinePoly &TSplinePoly::operator=(TSplinePoly const &other)
{
   if(this != &other) {
      TObject::operator=(other);
      CopyPoly(other);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Utility called by the copy constructors and = operator.

void TSplinePoly::CopyPoly(TSplinePoly const &other)
{
   fX = other.fX;
   fY = other.fY;
}

/** \class TSplinePoly3
    \ingroup Hist
 Class for TSpline3 knot.
*/

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TSplinePoly3 &TSplinePoly3::operator=(TSplinePoly3 const &other)
{
   if(this != &other) {
      TSplinePoly::operator=(other);
      CopyPoly(other);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Utility called by the copy constructors and = operator.

void TSplinePoly3::CopyPoly(TSplinePoly3 const &other)
{
   fB = other.fB;
   fC = other.fC;
   fD = other.fD;
}

/** \class TSplinePoly5
    \ingroup Hist
 Class for TSpline5 knot.
*/

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TSplinePoly5 &TSplinePoly5::operator=(TSplinePoly5 const &other)
{
   if(this != &other) {
      TSplinePoly::operator=(other);
      CopyPoly(other);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Utility called by the copy constructors and = operator.

void TSplinePoly5::CopyPoly(TSplinePoly5 const &other)
{
   fB = other.fB;
   fC = other.fC;
   fD = other.fD;
   fE = other.fE;
   fF = other.fF;
}

/** \class TSpline3
    \ingroup Hist
 Class to create third splines to interpolate knots
 Arbitrary conditions can be introduced for first and second
 derivatives at beginning and ending points
 */

////////////////////////////////////////////////////////////////////////////////
/// Third spline creator given an array of arbitrary knots in increasing
/// abscissa order and possibly end point conditions.

TSpline3::TSpline3(const char *title,
                   Double_t x[], Double_t y[], Int_t n, const char *opt,
                   Double_t valbeg, Double_t valend) :
  TSpline(title,-1,x[0],x[n-1],n,kFALSE),
  fValBeg(valbeg), fValEnd(valend), fBegCond(0), fEndCond(0)
{
   fName="Spline3";

   // Set endpoint conditions
   if(opt) SetCond(opt);

   // Create the polynomial terms and fill
   // them with node information
   fPoly = new TSplinePoly3[n];
   for (Int_t i=0; i<n; ++i) {
      fPoly[i].X() = x[i];
      fPoly[i].Y() = y[i];
   }

   // Build the spline coefficients
   BuildCoeff();
}

////////////////////////////////////////////////////////////////////////////////
/// Third spline creator given an array of
/// arbitrary function values on equidistant n abscissa
/// values from xmin to xmax and possibly end point conditions.

TSpline3::TSpline3(const char *title,
                   Double_t xmin, Double_t xmax,
                   Double_t y[], Int_t n, const char *opt,
                   Double_t valbeg, Double_t valend) :
  TSpline(title,(xmax-xmin)/(n-1), xmin, xmax, n, kTRUE),
  fValBeg(valbeg), fValEnd(valend),
  fBegCond(0), fEndCond(0)
{
   fName="Spline3";

   // Set endpoint conditions
   if(opt) SetCond(opt);

   // Create the polynomial terms and fill
   // them with node information
   fPoly = new TSplinePoly3[n];
   for (Int_t i=0; i<n; ++i) {
      fPoly[i].X() = fXmin+i*fDelta;
      fPoly[i].Y() = y[i];
   }

   // Build the spline coefficients
   BuildCoeff();
}

////////////////////////////////////////////////////////////////////////////////
/// Third spline creator given an array of
/// arbitrary abscissas in increasing order and a function
/// to interpolate and possibly end point conditions.

TSpline3::TSpline3(const char *title,
                   Double_t x[], const TF1 *func, Int_t n, const char *opt,
                   Double_t valbeg, Double_t valend) :
  TSpline(title,-1, x[0], x[n-1], n, kFALSE),
  fValBeg(valbeg), fValEnd(valend),
  fBegCond(0), fEndCond(0)
{
   fName="Spline3";

   // Set endpoint conditions
   if(opt) SetCond(opt);

   // Create the polynomial terms and fill
   // them with node information
   fPoly = new TSplinePoly3[n];
   for (Int_t i=0; i<n; ++i) {
      fPoly[i].X() = x[i];
      fPoly[i].Y() = ((TF1*)func)->Eval(x[i]);
   }

   // Build the spline coefficients
   BuildCoeff();
}

////////////////////////////////////////////////////////////////////////////////
/// Third spline creator given a function to be
/// evaluated on n equidistant abscissa points between xmin
/// and xmax and possibly end point conditions.

TSpline3::TSpline3(const char *title,
                   Double_t xmin, Double_t xmax,
                   const TF1 *func, Int_t n, const char *opt,
                   Double_t valbeg, Double_t valend) :
  TSpline(title,(xmax-xmin)/(n-1), xmin, xmax, n, kTRUE),
  fValBeg(valbeg), fValEnd(valend),
  fBegCond(0), fEndCond(0)
{
   fName="Spline3";

   // Set endpoint conditions
   if(opt) SetCond(opt);

   // Create the polynomial terms and fill
   // them with node information
   fPoly = new TSplinePoly3[n];
   //when func is null we return. In this case it is assumed that the spline
   //points will be given later via SetPoint and SetPointCoeff
   if (!func) {fKstep = kFALSE; fDelta = -1; return;}
   for (Int_t i=0; i<n; ++i) {
      Double_t x=fXmin+i*fDelta;
      fPoly[i].X() = x;
      fPoly[i].Y() = ((TF1*)func)->Eval(x);
   }

   // Build the spline coefficients
   BuildCoeff();
}

////////////////////////////////////////////////////////////////////////////////
/// Third spline creator given a TGraph with
/// abscissa in increasing order and possibly end
/// point conditions.

TSpline3::TSpline3(const char *title,
                   const TGraph *g, const char *opt,
                   Double_t valbeg, Double_t valend) :
  TSpline(title,-1,0,0,g->GetN(),kFALSE),
  fValBeg(valbeg), fValEnd(valend),
  fBegCond(0), fEndCond(0)
{
   fName="Spline3";

   // Set endpoint conditions
   if(opt) SetCond(opt);

   // Create the polynomial terms and fill
   // them with node information
   fPoly = new TSplinePoly3[fNp];
   for (Int_t i=0; i<fNp; ++i) {
      Double_t xx, yy;
      g->GetPoint(i,xx,yy);
      fPoly[i].X()=xx;
      fPoly[i].Y()=yy;
   }
   fXmin = fPoly[0].X();
   fXmax = fPoly[fNp-1].X();

   // Build the spline coefficients
   BuildCoeff();
}

////////////////////////////////////////////////////////////////////////////////
/// Third spline creator given a TH1.

TSpline3::TSpline3(const TH1 *h, const char *opt,
                   Double_t valbeg, Double_t valend) :
  TSpline(h->GetTitle(),-1,0,0,h->GetNbinsX(),kFALSE),
  fValBeg(valbeg), fValEnd(valend),
  fBegCond(0), fEndCond(0)
{
   fName=h->GetName();

   // Set endpoint conditions
   if(opt) SetCond(opt);

   // Create the polynomial terms and fill
   // them with node information
   fPoly = new TSplinePoly3[fNp];
   for (Int_t i=0; i<fNp; ++i) {
      fPoly[i].X()=h->GetXaxis()->GetBinCenter(i+1);
      fPoly[i].Y()=h->GetBinContent(i+1);
   }
   fXmin = fPoly[0].X();
   fXmax = fPoly[fNp-1].X();

   // Build the spline coefficients
   BuildCoeff();
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TSpline3::TSpline3(const TSpline3& sp3) :
  TSpline(sp3),
  fPoly(0),
  fValBeg(sp3.fValBeg),
  fValEnd(sp3.fValEnd),
  fBegCond(sp3.fBegCond),
  fEndCond(sp3.fEndCond)
{
   if (fNp > 0) fPoly = new TSplinePoly3[fNp];
   for (Int_t i=0; i<fNp; ++i)
      fPoly[i] = sp3.fPoly[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TSpline3& TSpline3::operator=(const TSpline3& sp3)
{
   if(this!=&sp3) {
      TSpline::operator=(sp3);
      fPoly= 0;
      if (fNp > 0) fPoly = new TSplinePoly3[fNp];
      for (Int_t i=0; i<fNp; ++i)
         fPoly[i] = sp3.fPoly[i];

      fValBeg=sp3.fValBeg;
      fValEnd=sp3.fValEnd;
      fBegCond=sp3.fBegCond;
      fEndCond=sp3.fEndCond;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the boundary conditions.

void TSpline3::SetCond(const char *opt)
{
   const char *b1 = strstr(opt,"b1");
   const char *e1 = strstr(opt,"e1");
   const char *b2 = strstr(opt,"b2");
   const char *e2 = strstr(opt,"e2");
   if (b1 && b2)
      Error("SetCond","Cannot specify first and second derivative at first point");
   if (e1 && e2)
      Error("SetCond","Cannot specify first and second derivative at last point");
   if (b1) fBegCond=1;
   else if (b2) fBegCond=2;
   if (e1) fEndCond=1;
   else if (e2) fEndCond=2;
}

////////////////////////////////////////////////////////////////////////////////
/// Test method for TSpline5
///
/// ~~~ {.cpp}
///   n          number of data points.
///   m          2*m-1 is order of spline.
///                 m = 2 always for third spline.
///   nn,nm1,mm,
///   mm1,i,k,
///   j,jj       temporary integer variables.
///   z,p        temporary double precision variables.
///   x[n]       the sequence of knots.
///   y[n]       the prescribed function values at the knots.
///   a[200][4]  two dimensional array whose columns are
///                 the computed spline coefficients
///   diff[3]    maximum values of differences of values and
///                 derivatives to right and left of knots.
///   com[3]     maximum values of coefficients.
/// ~~~
///
/// test of TSpline3 with non equidistant knots and
/// equidistant knots follows.

void TSpline3::Test()
{
   Double_t hx;
   Double_t diff[3];
   Double_t a[800], c[4];
   Int_t i, j, k, m, n;
   Double_t x[200], y[200], z;
   Int_t jj, mm;
   Int_t mm1, nm1;
   Double_t com[3];
   printf("1         TEST OF TSpline3 WITH NONEQUIDISTANT KNOTS\n");
   n = 5;
   x[0] = -3;
   x[1] = -1;
   x[2] = 0;
   x[3] = 3;
   x[4] = 4;
   y[0] = 7;
   y[1] = 11;
   y[2] = 26;
   y[3] = 56;
   y[4] = 29;
   m = 2;
   mm = m << 1;
   mm1 = mm-1;
   printf("\n-N = %3d    M =%2d\n",n,m);
   TSpline3 *spline = new TSpline3("Test",x,y,n);
   for (i = 0; i < n; ++i)
      spline->GetCoeff(i,hx, a[i],a[i+200],a[i+400],a[i+600]);
   delete spline;
   for (i = 0; i < mm1; ++i) diff[i] = com[i] = 0;
   for (k = 0; k < n; ++k) {
      for (i = 0; i < mm; ++i) c[i] = a[k+i*200];
      printf(" ---------------------------------------%3d --------------------------------------------\n",k+1);
      printf("%12.8f\n",x[k]);
      if (k == n-1) {
         printf("%16.8f\n",c[0]);
      } else {
         for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
         printf("\n");
         for (i = 0; i < mm1; ++i)
            if ((z=TMath::Abs(a[k+i*200])) > com[i]) com[i] = z;
         z = x[k+1]-x[k];
         for (i = 1; i < mm; ++i)
            for (jj = i; jj < mm; ++jj) {
               j = mm+i-jj;
               c[j-2] = c[j-1]*z+c[j-2];
            }
         for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
         printf("\n");
         for (i = 0; i < mm1; ++i)
            if (!(k >= n-2 && i != 0))
               if((z = TMath::Abs(c[i]-a[k+1+i*200]))
                  > diff[i]) diff[i] = z;
      }
   }
   printf("  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES \n");
   for (i = 0; i < mm1; ++i) printf("%18.9E",diff[i]);
   printf("\n");
   printf("  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS \n");
   if (TMath::Abs(c[0]) > com[0])
      com[0] = TMath::Abs(c[0]);
   for (i = 0; i < mm1; ++i) printf("%16.8f",com[i]);
   printf("\n");
   m = 2;
   for (n = 10; n <= 100; n += 10) {
      mm = m << 1;
      mm1 = mm-1;
      nm1 = n-1;
      for (i = 0; i < nm1; i += 2) {
         x[i] = i+1;
         x[i+1] = i+2;
         y[i] = 1;
         y[i+1] = 0;
      }
      if (n % 2 != 0) {
         x[n-1] = n;
         y[n-1] = 1;
      }
      printf("\n-N = %3d    M =%2d\n",n,m);
      spline = new TSpline3("Test",x,y,n);
      for (i = 0; i < n; ++i)
         spline->GetCoeff(i,hx,a[i],a[i+200],a[i+400],a[i+600]);
      delete spline;
      for (i = 0; i < mm1; ++i)
         diff[i] = com[i] = 0;
      for (k = 0; k < n; ++k) {
         for (i = 0; i < mm; ++i)
            c[i] = a[k+i*200];
         if (n < 11) {
            printf(" ---------------------------------------%3d --------------------------------------------\n",k+1);
            printf("%12.8f\n",x[k]);
            if (k == n-1) printf("%16.8f\n",c[0]);
         }
         if (k == n-1) break;
         if (n <= 10) {
            for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
            printf("\n");
         }
         for (i = 0; i < mm1; ++i)
         if ((z=TMath::Abs(a[k+i*200])) > com[i])
            com[i] = z;
         z = x[k+1]-x[k];
         for (i = 1; i < mm; ++i)
            for (jj = i; jj < mm; ++jj) {
               j = mm+i-jj;
               c[j-2] = c[j-1]*z+c[j-2];
            }
         if (n <= 10) {
            for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
            printf("\n");
         }
         for (i = 0; i < mm1; ++i)
         if (!(k >= n-2 && i != 0))
            if ((z = TMath::Abs(c[i]-a[k+1+i*200]))
               > diff[i]) diff[i] = z;
      }
      printf("  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES \n");
      for (i = 0; i < mm1; ++i) printf("%18.9E",diff[i]);
      printf("\n");
      printf("  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS \n");
      if (TMath::Abs(c[0]) > com[0])
         com[0] = TMath::Abs(c[0]);
      for (i = 0; i < mm1; ++i) printf("%16.8E",com[i]);
         printf("\n");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find X.

Int_t TSpline3::FindX(Double_t x) const
{
   Int_t klow=0, khig=fNp-1;
   //
   // If out of boundaries, extrapolate
   // It may be badly wrong
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=khig;
   else {
      if(fKstep) {
         //
         // Equidistant knots, use histogramming
         klow = TMath::FloorNint((x-fXmin)/fDelta);
         // Correction for rounding errors
         if (x < fPoly[klow].X())
            klow = TMath::Max(klow-1,0);
         else if (klow < khig) {
            if (x > fPoly[klow+1].X()) ++klow;
         }
      } else {
         Int_t khalf;
         //
         // Non equidistant knots, binary search
         while(khig-klow>1)
            if(x>fPoly[khalf=(klow+khig)/2].X())
               klow=khalf;
            else
               khig=khalf;
         //
         // This could be removed, sanity check
         if(!(fPoly[klow].X()<=x && x<=fPoly[klow+1].X()))
            Error("Eval",
                  "Binary search failed x(%d) = %f < x= %f < x(%d) = %f\n",
                  klow,fPoly[klow].X(),x,klow+1,fPoly[klow+1].X());
      }
   }
   return klow;
}

////////////////////////////////////////////////////////////////////////////////
/// Eval this spline at x.

Double_t TSpline3::Eval(Double_t x) const
{
   Int_t klow=FindX(x);
   if (klow >= fNp-1 && fNp > 1) klow = fNp-2; //see: https://savannah.cern.ch/bugs/?71651
   return fPoly[klow].Eval(x);
}

////////////////////////////////////////////////////////////////////////////////
/// Derivative.

Double_t TSpline3::Derivative(Double_t x) const
{
   Int_t klow=FindX(x);
   if (klow >= fNp-1) klow = fNp-2; //see: https://savannah.cern.ch/bugs/?71651
   return fPoly[klow].Derivative(x);
}

////////////////////////////////////////////////////////////////////////////////
/// Write this spline as a C++ function that can be executed without ROOT
/// the name of the function is the name of the file up to the "." if any.

void TSpline3::SaveAs(const char *filename, Option_t * /*option*/) const
{
   //open the file
   std::ofstream *f = new std::ofstream(filename,std::ios::out);
   if (f == 0 || gSystem->AccessPathName(filename,kWritePermission)) {
      Error("SaveAs","Cannot open file:%s\n",filename);
      return;
   }

   //write the function name and the spline constants
   char buffer[512];
   Int_t nch = strlen(filename);
   snprintf(buffer,512,"double %s",filename);
   char *dot = strstr(buffer,".");
   if (dot) *dot = 0;
   strlcat(buffer,"(double x) {\n",512);
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   const int fNp = %d, fKstep = %d;\n",fNp,fKstep);
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   const double fDelta = %g, fXmin = %g, fXmax = %g;\n",fDelta,fXmin,fXmax);
   nch = strlen(buffer); f->write(buffer,nch);

   //write the spline coefficients
   //array fX
   snprintf(buffer,512,"   const double fX[%d] = {",fNp);
   nch = strlen(buffer); f->write(buffer,nch);
   buffer[0] = 0;
   Int_t i;
   char numb[20];
   for (i=0;i<fNp;i++) {
      snprintf(numb,20," %g,",fPoly[i].X());
      nch = strlen(numb);
      if (i == fNp-1) numb[nch-1]=0;
      strlcat(buffer,numb,512);
      if (i%5 == 4 || i == fNp-1) {
         nch = strlen(buffer); f->write(buffer,nch);
         if (i != fNp-1) snprintf(buffer,512,"\n                       ");
      }
   }
   snprintf(buffer,512," };\n");
   nch = strlen(buffer); f->write(buffer,nch);
   //array fY
   snprintf(buffer,512,"   const double fY[%d] = {",fNp);
   nch = strlen(buffer); f->write(buffer,nch);
   buffer[0] = 0;
   for (i=0;i<fNp;i++) {
      snprintf(numb,20," %g,",fPoly[i].Y());
      nch = strlen(numb);
      if (i == fNp-1) numb[nch-1]=0;
      strlcat(buffer,numb,512);
      if (i%5 == 4 || i == fNp-1) {
         nch = strlen(buffer); f->write(buffer,nch);
         if (i != fNp-1) snprintf(buffer,512,"\n                       ");
      }
   }
   snprintf(buffer,512," };\n");
   nch = strlen(buffer); f->write(buffer,nch);
   //array fB
   snprintf(buffer,512,"   const double fB[%d] = {",fNp);
   nch = strlen(buffer); f->write(buffer,nch);
   buffer[0] = 0;
   for (i=0;i<fNp;i++) {
      snprintf(numb,20," %g,",fPoly[i].B());
      nch = strlen(numb);
      if (i == fNp-1) numb[nch-1]=0;
      strlcat(buffer,numb,512);
      if (i%5 == 4 || i == fNp-1) {
         nch = strlen(buffer); f->write(buffer,nch);
         if (i != fNp-1) snprintf(buffer,512,"\n                       ");
      }
   }
   snprintf(buffer,512," };\n");
   nch = strlen(buffer); f->write(buffer,nch);
   //array fC
   snprintf(buffer,512,"   const double fC[%d] = {",fNp);
   nch = strlen(buffer); f->write(buffer,nch);
   buffer[0] = 0;
   for (i=0;i<fNp;i++) {
      snprintf(numb,20," %g,",fPoly[i].C());
      nch = strlen(numb);
      if (i == fNp-1) numb[nch-1]=0;
      strlcat(buffer,numb,512);
      if (i%5 == 4 || i == fNp-1) {
         nch = strlen(buffer); f->write(buffer,nch);
         if (i != fNp-1) snprintf(buffer,512,"\n                       ");
      }
   }
   snprintf(buffer,512," };\n");
   nch = strlen(buffer); f->write(buffer,nch);
    //array fD
   snprintf(buffer,512,"   const double fD[%d] = {",fNp);
   nch = strlen(buffer); f->write(buffer,nch);
   buffer[0] = 0;
   for (i=0;i<fNp;i++) {
      snprintf(numb,20," %g,",fPoly[i].D());
      nch = strlen(numb);
      if (i == fNp-1) numb[nch-1]=0;
      strlcat(buffer,numb,512);
      if (i%5 == 4 || i == fNp-1) {
         nch = strlen(buffer); f->write(buffer,nch);
         if (i != fNp-1) snprintf(buffer,512,"\n                       ");
      }
   }
   snprintf(buffer,512," };\n");
   nch = strlen(buffer); f->write(buffer,nch);

   //generate code for the spline evaluation
   snprintf(buffer,512,"   int klow=0;\n");
   nch = strlen(buffer); f->write(buffer,nch);

   snprintf(buffer,512,"   // If out of boundaries, extrapolate. It may be badly wrong\n");
   snprintf(buffer,512,"   if(x<=fXmin) klow=0;\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   else if(x>=fXmax) klow=fNp-1;\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   else {\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"     if(fKstep) {\n");
   nch = strlen(buffer); f->write(buffer,nch);

   snprintf(buffer,512,"       // Equidistant knots, use histogramming\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"       klow = int((x-fXmin)/fDelta);\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"       if (klow < fNp-1) klow = fNp-1;\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"     } else {\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"       int khig=fNp-1, khalf;\n");
   nch = strlen(buffer); f->write(buffer,nch);

   snprintf(buffer,512,"       // Non equidistant knots, binary search\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"       while(khig-klow>1)\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"         else khig=khalf;\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"     }\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   }\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   // Evaluate now\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   double dx=x-fX[klow];\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));\n");
   nch = strlen(buffer); f->write(buffer,nch);

   //close file
   f->write("}\n",2);

   if (f) { f->close(); delete f;}
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out.

void TSpline3::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';
   out<<"   "<<std::endl;
   if (gROOT->ClassSaved(TSpline3::Class())) {
      out<<"   ";
   } else {
      out<<"   TSpline3 *";
   }
   out<<"spline3 = new TSpline3("<<quote<<GetTitle()<<quote<<","
      <<fXmin<<","<<fXmax<<",(TF1*)0,"<<fNp<<","<<quote<<quote<<","
      <<fValBeg<<","<<fValEnd<<");"<<std::endl;
   out<<"   spline3->SetName("<<quote<<GetName()<<quote<<");"<<std::endl;

   SaveFillAttributes(out,"spline3",0,1001);
   SaveLineAttributes(out,"spline3",1,1,1);
   SaveMarkerAttributes(out,"spline3",1,1,1);
   if (fNpx != 100) out<<"   spline3->SetNpx("<<fNpx<<");"<<std::endl;

   for (Int_t i=0;i<fNp;i++) {
      out<<"   spline3->SetPoint("<<i<<","<<fPoly[i].X()<<","<<fPoly[i].Y()<<");"<<std::endl;
      out<<"   spline3->SetPointCoeff("<<i<<","<<fPoly[i].B()<<","<<fPoly[i].C()<<","<<fPoly[i].D()<<");"<<std::endl;
   }
   out<<"   spline3->Draw("<<quote<<option<<quote<<");"<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Set point number i.

void TSpline3::SetPoint(Int_t i, Double_t x, Double_t y)
{
   if (i < 0 || i >= fNp) return;
   fPoly[i].X()= x;
   fPoly[i].Y()= y;
}

////////////////////////////////////////////////////////////////////////////////
/// Set point coefficient number i.

void TSpline3::SetPointCoeff(Int_t i, Double_t b, Double_t c, Double_t d)
{
   if (i < 0 || i >= fNp) return;
   fPoly[i].B()= b;
   fPoly[i].C()= c;
   fPoly[i].D()= d;
}

////////////////////////////////////////////////////////////////////////////////
/// Build coefficients.
///
/// ~~~ {.cpp}
///      subroutine cubspl ( tau, c, n, ibcbeg, ibcend )
///  from  * a practical guide to splines *  by c. de boor
///     ************************  input  ***************************
///     n = number of data points. assumed to be .ge. 2.
///     (tau(i), c(1,i), i=1,...,n) = abscissae and ordinates of the
///        data points. tau is assumed to be strictly increasing.
///     ibcbeg, ibcend = boundary condition indicators, and
///     c(2,1), c(2,n) = boundary condition information. specifically,
///        ibcbeg = 0  means no boundary condition at tau(1) is given.
///           in this case, the not-a-knot condition is used, i.e. the
///           jump in the third derivative across tau(2) is forced to
///           zero, thus the first and the second cubic polynomial pieces
///           are made to coincide.)
///        ibcbeg = 1  means that the slope at tau(1) is made to equal
///           c(2,1), supplied by input.
///        ibcbeg = 2  means that the second derivative at tau(1) is
///           made to equal c(2,1), supplied by input.
///        ibcend = 0, 1, or 2 has analogous meaning concerning the
///           boundary condition at tau(n), with the additional infor-
///           mation taken from c(2,n).
///     ***********************  output  **************************
///     c(j,i), j=1,...,4; i=1,...,l (= n-1) = the polynomial coefficients
///        of the cubic interpolating spline with interior knots (or
///        joints) tau(2), ..., tau(n-1). precisely, in the interval
///        (tau(i), tau(i+1)), the spline f is given by
///           f(x) = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
///        where h = x - tau(i). the function program *ppvalu* may be
///        used to evaluate f or its derivatives from tau,c, l = n-1,
///        and k=4.
/// ~~~

void TSpline3::BuildCoeff()
{
   Int_t i, j, l, m;
   Double_t   divdf1,divdf3,dtau,g=0;
   //***** a tridiagonal linear system for the unknown slopes s(i) of
   //  f  at tau(i), i=1,...,n, is generated and then solved by gauss elim-
   //  ination, with s(i) ending up in c(2,i), all i.
   //     c(3,.) and c(4,.) are used initially for temporary storage.
   l = fNp-1;
   // compute first differences of x sequence and store in C also,
   // compute first divided difference of data and store in D.
   for (m=1; m<fNp ; ++m) {
      fPoly[m].C() = fPoly[m].X() - fPoly[m-1].X();
      fPoly[m].D() = (fPoly[m].Y() - fPoly[m-1].Y())/fPoly[m].C();
   }
   // construct first equation from the boundary condition, of the form
   //             D[0]*s[0] + C[0]*s[1] = B[0]
   if(fBegCond==0) {
      if(fNp == 2) {
         //     no condition at left end and n = 2.
         fPoly[0].D() = 1.;
         fPoly[0].C() = 1.;
         fPoly[0].B() = 2.*fPoly[1].D();
      } else {
         //     not-a-knot condition at left end and n .gt. 2.
         fPoly[0].D() = fPoly[2].C();
         fPoly[0].C() = fPoly[1].C() + fPoly[2].C();
         fPoly[0].B() =((fPoly[1].C()+2.*fPoly[0].C())*fPoly[1].D()*fPoly[2].C()+fPoly[1].C()*fPoly[1].C()*fPoly[2].D())/fPoly[0].C();
      }
   } else if (fBegCond==1) {
      //     slope prescribed at left end.
      fPoly[0].B() = fValBeg;
      fPoly[0].D() = 1.;
      fPoly[0].C() = 0.;
   } else if (fBegCond==2) {
      //     second derivative prescribed at left end.
      fPoly[0].D() = 2.;
      fPoly[0].C() = 1.;
      fPoly[0].B() = 3.*fPoly[1].D() - fPoly[1].C()/2.*fValBeg;
   }
   if(fNp > 2) {
      //  if there are interior knots, generate the corresp. equations and car-
      //  ry out the forward pass of gauss elimination, after which the m-th
      //  equation reads    D[m]*s[m] + C[m]*s[m+1] = B[m].
      for (m=1; m<l; ++m) {
         g = -fPoly[m+1].C()/fPoly[m-1].D();
         fPoly[m].B() = g*fPoly[m-1].B() + 3.*(fPoly[m].C()*fPoly[m+1].D()+fPoly[m+1].C()*fPoly[m].D());
         fPoly[m].D() = g*fPoly[m-1].C() + 2.*(fPoly[m].C() + fPoly[m+1].C());
      }
      // construct last equation from the second boundary condition, of the form
      //           (-g*D[n-2])*s[n-2] + D[n-1]*s[n-1] = B[n-1]
      //     if slope is prescribed at right end, one can go directly to back-
      //     substitution, since c array happens to be set up just right for it
      //     at this point.
      if(fEndCond == 0) {
         if (fNp > 3 || fBegCond != 0) {
            //     not-a-knot and n .ge. 3, and either n.gt.3 or  also not-a-knot at
            //     left end point.
            g = fPoly[fNp-2].C() + fPoly[fNp-1].C();
            fPoly[fNp-1].B() = ((fPoly[fNp-1].C()+2.*g)*fPoly[fNp-1].D()*fPoly[fNp-2].C()
                         + fPoly[fNp-1].C()*fPoly[fNp-1].C()*(fPoly[fNp-2].Y()-fPoly[fNp-3].Y())/fPoly[fNp-2].C())/g;
            g = -g/fPoly[fNp-2].D();
            fPoly[fNp-1].D() = fPoly[fNp-2].C();
         } else {
            //     either (n=3 and not-a-knot also at left) or (n=2 and not not-a-
            //     knot at left end point).
            fPoly[fNp-1].B() = 2.*fPoly[fNp-1].D();
            fPoly[fNp-1].D() = 1.;
            g = -1./fPoly[fNp-2].D();
         }
      } else if (fEndCond == 1) {
         fPoly[fNp-1].B() = fValEnd;
         goto L30;
      } else if (fEndCond == 2) {
         //     second derivative prescribed at right endpoint.
         fPoly[fNp-1].B() = 3.*fPoly[fNp-1].D() + fPoly[fNp-1].C()/2.*fValEnd;
         fPoly[fNp-1].D() = 2.;
         g = -1./fPoly[fNp-2].D();
      }
   } else {
      if(fEndCond == 0) {
         if (fBegCond > 0) {
            //     either (n=3 and not-a-knot also at left) or (n=2 and not not-a-
            //     knot at left end point).
            fPoly[fNp-1].B() = 2.*fPoly[fNp-1].D();
            fPoly[fNp-1].D() = 1.;
            g = -1./fPoly[fNp-2].D();
         } else {
            //     not-a-knot at right endpoint and at left endpoint and n = 2.
            fPoly[fNp-1].B() = fPoly[fNp-1].D();
            goto L30;
         }
      } else if(fEndCond == 1) {
         fPoly[fNp-1].B() = fValEnd;
         goto L30;
      } else if(fEndCond == 2) {
         //     second derivative prescribed at right endpoint.
         fPoly[fNp-1].B() = 3.*fPoly[fNp-1].D() + fPoly[fNp-1].C()/2.*fValEnd;
         fPoly[fNp-1].D() = 2.;
         g = -1./fPoly[fNp-2].D();
      }
   }
   // complete forward pass of gauss elimination.
   fPoly[fNp-1].D() = g*fPoly[fNp-2].C() + fPoly[fNp-1].D();
   fPoly[fNp-1].B() = (g*fPoly[fNp-2].B() + fPoly[fNp-1].B())/fPoly[fNp-1].D();
   // carry out back substitution
L30: j = l-1;
   do {
      fPoly[j].B() = (fPoly[j].B() - fPoly[j].C()*fPoly[j+1].B())/fPoly[j].D();
      --j;
   }  while (j>=0);
   //****** generate cubic coefficients in each interval, i.e., the deriv.s
   //  at its left endpoint, from value and slope at its endpoints.
   for (i=1; i<fNp; ++i) {
      dtau = fPoly[i].C();
      divdf1 = (fPoly[i].Y() - fPoly[i-1].Y())/dtau;
      divdf3 = fPoly[i-1].B() + fPoly[i].B() - 2.*divdf1;
      fPoly[i-1].C() = (divdf1 - fPoly[i-1].B() - divdf3)/dtau;
      fPoly[i-1].D() = (divdf3/dtau)/dtau;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TSpline3.

void TSpline3::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TSpline3::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TSpline::Streamer(R__b);
      if (fNp > 0) {
         fPoly = new TSplinePoly3[fNp];
         for(Int_t i=0; i<fNp; ++i) {
            fPoly[i].Streamer(R__b);
         }
      }
      //      R__b >> fPoly;
      R__b >> fValBeg;
      R__b >> fValEnd;
      R__b >> fBegCond;
      R__b >> fEndCond;
   } else {
      R__b.WriteClassBuffer(TSpline3::Class(),this);
   }
}

/** \class TSpline5
    \ingroup Hist
 Class to create quintic natural splines to interpolate knots
 Arbitrary conditions can be introduced for first and second
 derivatives using double knots (see BuildCoeff) for more on this.
 Double knots are automatically introduced at ending points
 */

////////////////////////////////////////////////////////////////////////////////
/// Quintic natural spline creator given an array of
/// arbitrary knots in increasing abscissa order and
/// possibly end point conditions.

TSpline5::TSpline5(const char *title,
                   Double_t x[], Double_t y[], Int_t n,
                   const char *opt, Double_t b1, Double_t e1,
                   Double_t b2, Double_t e2) :
  TSpline(title,-1, x[0], x[n-1], n, kFALSE)
{
   Int_t beg, end;
   const char *cb1, *ce1, *cb2, *ce2;
   fName="Spline5";

   // Check endpoint conditions
   BoundaryConditions(opt,beg,end,cb1,ce1,cb2,ce2);

   // Create the polynomial terms and fill
   // them with node information
   fPoly = new TSplinePoly5[fNp];
   for (Int_t i=0; i<n; ++i) {
      fPoly[i+beg].X() = x[i];
      fPoly[i+beg].Y() = y[i];
   }

   // Set the double knots at boundaries
   SetBoundaries(b1,e1,b2,e2,cb1,ce1,cb2,ce2);

   // Build the spline coefficients
   BuildCoeff();
}

////////////////////////////////////////////////////////////////////////////////
/// Quintic natural spline creator given an array of
/// arbitrary function values on equidistant n abscissa
/// values from xmin to xmax and possibly end point conditions.

TSpline5::TSpline5(const char *title,
                   Double_t xmin, Double_t xmax,
                   Double_t y[], Int_t n,
                   const char *opt, Double_t b1, Double_t e1,
                   Double_t b2, Double_t e2) :
  TSpline(title,(xmax-xmin)/(n-1), xmin, xmax, n, kTRUE)
{
   Int_t beg, end;
   const char *cb1, *ce1, *cb2, *ce2;
   fName="Spline5";

   // Check endpoint conditions
   BoundaryConditions(opt,beg,end,cb1,ce1,cb2,ce2);

   // Create the polynomial terms and fill
   // them with node information
   fPoly = new TSplinePoly5[fNp];
   for (Int_t i=0; i<n; ++i) {
      fPoly[i+beg].X() = fXmin+i*fDelta;
      fPoly[i+beg].Y() = y[i];
   }

   // Set the double knots at boundaries
   SetBoundaries(b1,e1,b2,e2,cb1,ce1,cb2,ce2);

   // Build the spline coefficients
   BuildCoeff();
}

////////////////////////////////////////////////////////////////////////////////
/// Quintic natural spline creator given an array of
/// arbitrary abscissas in increasing order and a function
/// to interpolate and possibly end point conditions.

TSpline5::TSpline5(const char *title,
                   Double_t x[], const TF1 *func, Int_t n,
                   const char *opt, Double_t b1, Double_t e1,
                   Double_t b2, Double_t e2) :
  TSpline(title,-1, x[0], x[n-1], n, kFALSE)
{
   Int_t beg, end;
   const char *cb1, *ce1, *cb2, *ce2;
   fName="Spline5";

   // Check endpoint conditions
   BoundaryConditions(opt,beg,end,cb1,ce1,cb2,ce2);

   // Create the polynomial terms and fill
   // them with node information
   fPoly = new TSplinePoly5[fNp];
   for (Int_t i=0; i<n; i++) {
      fPoly[i+beg].X() = x[i];
      fPoly[i+beg].Y() = ((TF1*)func)->Eval(x[i]);
   }

   // Set the double knots at boundaries
   SetBoundaries(b1,e1,b2,e2,cb1,ce1,cb2,ce2);

   // Build the spline coefficients
   BuildCoeff();
}

////////////////////////////////////////////////////////////////////////////////
/// Quintic natural spline creator given a function to be
/// evaluated on n equidistant abscissa points between xmin
/// and xmax and possibly end point conditions.

TSpline5::TSpline5(const char *title,
                   Double_t xmin, Double_t xmax,
                   const TF1 *func, Int_t n,
                   const char *opt, Double_t b1, Double_t e1,
                   Double_t b2, Double_t e2) :
  TSpline(title,(xmax-xmin)/(n-1), xmin, xmax, n, kTRUE)
{
   Int_t beg, end;
   const char *cb1, *ce1, *cb2, *ce2;
   fName="Spline5";

   // Check endpoint conditions
   BoundaryConditions(opt,beg,end,cb1,ce1,cb2,ce2);

   // Create the polynomial terms and fill
   // them with node information
   fPoly = new TSplinePoly5[fNp];
   for (Int_t i=0; i<n; ++i) {
      Double_t x=fXmin+i*fDelta;
      fPoly[i+beg].X() = x;
      if (func) fPoly[i+beg].Y() = ((TF1*)func)->Eval(x);
   }
   if (!func) {fDelta = -1; fKstep = kFALSE;}

   // Set the double knots at boundaries
   SetBoundaries(b1,e1,b2,e2,cb1,ce1,cb2,ce2);

   // Build the spline coefficients
   if (func) BuildCoeff();
}

////////////////////////////////////////////////////////////////////////////////
/// Quintic natural spline creator given a TGraph with
/// abscissa in increasing order and possibly end
/// point conditions.

TSpline5::TSpline5(const char *title,
                   const TGraph *g,
                   const char *opt, Double_t b1, Double_t e1,
                   Double_t b2, Double_t e2) :
  TSpline(title,-1,0,0,g->GetN(),kFALSE)
{
   Int_t beg, end;
   const char *cb1, *ce1, *cb2, *ce2;
   fName="Spline5";

   // Check endpoint conditions
   BoundaryConditions(opt,beg,end,cb1,ce1,cb2,ce2);

   // Create the polynomial terms and fill
   // them with node information
   fPoly = new TSplinePoly5[fNp];
   for (Int_t i=0; i<fNp-beg; ++i) {
      Double_t xx, yy;
      g->GetPoint(i,xx,yy);
      fPoly[i+beg].X()=xx;
      fPoly[i+beg].Y()=yy;
   }

   // Set the double knots at boundaries
   SetBoundaries(b1,e1,b2,e2,cb1,ce1,cb2,ce2);
   fXmin = fPoly[0].X();
   fXmax = fPoly[fNp-1].X();

   // Build the spline coefficients
   BuildCoeff();
}

////////////////////////////////////////////////////////////////////////////////
/// Quintic natural spline creator given a TH1.

TSpline5::TSpline5(const TH1 *h,
                   const char *opt, Double_t b1, Double_t e1,
                   Double_t b2, Double_t e2) :
  TSpline(h->GetTitle(),-1,0,0,h->GetNbinsX(),kFALSE)
{
   Int_t beg, end;
   const char *cb1, *ce1, *cb2, *ce2;
   fName=h->GetName();

   // Check endpoint conditions
   BoundaryConditions(opt,beg,end,cb1,ce1,cb2,ce2);

   // Create the polynomial terms and fill
   // them with node information
   fPoly = new TSplinePoly5[fNp];
   for (Int_t i=0; i<fNp-beg; ++i) {
      fPoly[i+beg].X()=h->GetXaxis()->GetBinCenter(i+1);
      fPoly[i+beg].Y()=h->GetBinContent(i+1);
   }

   // Set the double knots at boundaries
   SetBoundaries(b1,e1,b2,e2,cb1,ce1,cb2,ce2);
   fXmin = fPoly[0].X();
   fXmax = fPoly[fNp-1].X();

   // Build the spline coefficients
   BuildCoeff();
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TSpline5::TSpline5(const TSpline5& sp5) :
  TSpline(sp5),
  fPoly(0)
{
   if (fNp > 0) fPoly = new TSplinePoly5[fNp];
   for (Int_t i=0; i<fNp; ++i) {
      fPoly[i] = sp5.fPoly[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TSpline5& TSpline5::operator=(const TSpline5& sp5)
{
   if(this!=&sp5) {
      TSpline::operator=(sp5);
      fPoly=0;
      if (fNp > 0) fPoly = new TSplinePoly5[fNp];
      for (Int_t i=0; i<fNp; ++i) {
         fPoly[i] = sp5.fPoly[i];
      }
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the boundary conditions and the
/// amount of extra double knots needed.

void TSpline5::BoundaryConditions(const char *opt,Int_t &beg,Int_t &end,
                                  const char *&cb1,const char *&ce1,
                                  const char *&cb2,const char *&ce2)
{
   cb1=ce1=cb2=ce2=0;
   beg=end=0;
   if(opt) {
      cb1 = strstr(opt,"b1");
      ce1 = strstr(opt,"e1");
      cb2 = strstr(opt,"b2");
      ce2 = strstr(opt,"e2");
      if(cb2) {
         fNp=fNp+2;
         beg=2;
      } else if(cb1) {
         fNp=fNp+1;
         beg=1;
      }
      if(ce2) {
         fNp=fNp+2;
         end=2;
      } else if(ce1) {
         fNp=fNp+1;
         end=1;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the boundary conditions at double/triple knots.

void TSpline5::SetBoundaries(Double_t b1, Double_t e1, Double_t b2, Double_t e2,
                             const char *cb1, const char *ce1, const char *cb2,
                             const char *ce2)
{
   if(cb2) {

      // Second derivative at the beginning
      fPoly[0].X() = fPoly[1].X() = fPoly[2].X();
      fPoly[0].Y() = fPoly[2].Y();
      fPoly[2].Y()=b2;

      // If first derivative not given, we take the finite
      // difference from first and second point... not recommended
      if(cb1)
         fPoly[1].Y()=b1;
      else
         fPoly[1].Y()=(fPoly[3].Y()-fPoly[0].Y())/(fPoly[3].X()-fPoly[2].X());
   } else if(cb1) {

      // First derivative at the end
      fPoly[0].X() = fPoly[1].X();
      fPoly[0].Y() = fPoly[1].Y();
      fPoly[1].Y()=b1;
   }
   if(ce2) {

      // Second derivative at the end
      fPoly[fNp-1].X() = fPoly[fNp-2].X() = fPoly[fNp-3].X();
      fPoly[fNp-1].Y()=e2;

      // If first derivative not given, we take the finite
      // difference from first and second point... not recommended
      if(ce1)
         fPoly[fNp-2].Y()=e1;
      else
         fPoly[fNp-2].Y()=
         (fPoly[fNp-3].Y()-fPoly[fNp-4].Y())
         /(fPoly[fNp-3].X()-fPoly[fNp-4].X());
   } else if(ce1) {

      // First derivative at the end
      fPoly[fNp-1].X() = fPoly[fNp-2].X();
      fPoly[fNp-1].Y()=e1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find X.

Int_t TSpline5::FindX(Double_t x) const
{
   Int_t klow=0;

   // If out of boundaries, extrapolate
   // It may be badly wrong
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=fNp-1;
   else {
      if(fKstep) {

         // Equidistant knots, use histogramming
         klow = TMath::Min(Int_t((x-fXmin)/fDelta),fNp-1);
      } else {
         Int_t khig=fNp-1, khalf;

         // Non equidistant knots, binary search
         while(khig-klow>1)
            if(x>fPoly[khalf=(klow+khig)/2].X())
               klow=khalf;
            else
               khig=khalf;
      }

      // This could be removed, sanity check
      if(!(fPoly[klow].X()<=x && x<=fPoly[klow+1].X()))
         Error("Eval",
               "Binary search failed x(%d) = %f < x(%d) = %f\n",
                klow,fPoly[klow].X(),klow+1,fPoly[klow+1].X());
   }
   return klow;
}

////////////////////////////////////////////////////////////////////////////////
/// Eval this spline at x.

Double_t TSpline5::Eval(Double_t x) const
{
   Int_t klow=FindX(x);
   return fPoly[klow].Eval(x);
}

////////////////////////////////////////////////////////////////////////////////
/// Derivative.

Double_t TSpline5::Derivative(Double_t x) const
{
   Int_t klow=FindX(x);
   return fPoly[klow].Derivative(x);
}

////////////////////////////////////////////////////////////////////////////////
/// Write this spline as a C++ function that can be executed without ROOT
/// the name of the function is the name of the file up to the "." if any.

void TSpline5::SaveAs(const char *filename, Option_t * /*option*/) const
{
   //open the file
   std::ofstream *f = new std::ofstream(filename,std::ios::out);
   if (f == 0 || gSystem->AccessPathName(filename,kWritePermission)) {
      Error("SaveAs","Cannot open file:%s\n",filename);
      return;
   }

   //write the function name and the spline constants
   char buffer[512];
   Int_t nch = strlen(filename);
   snprintf(buffer,512,"double %s",filename);
   char *dot = strstr(buffer,".");
   if (dot) *dot = 0;
   strlcat(buffer,"(double x) {\n",512);
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   const int fNp = %d, fKstep = %d;\n",fNp,fKstep);
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   const double fDelta = %g, fXmin = %g, fXmax = %g;\n",fDelta,fXmin,fXmax);
   nch = strlen(buffer); f->write(buffer,nch);

   //write the spline coefficients
   //array fX
   snprintf(buffer,512,"   const double fX[%d] = {",fNp);
   nch = strlen(buffer); f->write(buffer,nch);
   buffer[0] = 0;
   Int_t i;
   char numb[20];
   for (i=0;i<fNp;i++) {
      snprintf(numb,20," %g,",fPoly[i].X());
      nch = strlen(numb);
      if (i == fNp-1) numb[nch-1]=0;
      strlcat(buffer,numb,512);
      if (i%5 == 4 || i == fNp-1) {
         nch = strlen(buffer); f->write(buffer,nch);
         if (i != fNp-1) snprintf(buffer,512,"\n                       ");
      }
   }
   snprintf(buffer,512," };\n");
   nch = strlen(buffer); f->write(buffer,nch);
   //array fY
   snprintf(buffer,512,"   const double fY[%d] = {",fNp);
   nch = strlen(buffer); f->write(buffer,nch);
   buffer[0] = 0;
   for (i=0;i<fNp;i++) {
      snprintf(numb,20," %g,",fPoly[i].Y());
      nch = strlen(numb);
      if (i == fNp-1) numb[nch-1]=0;
      strlcat(buffer,numb,512);
      if (i%5 == 4 || i == fNp-1) {
         nch = strlen(buffer); f->write(buffer,nch);
         if (i != fNp-1) snprintf(buffer,512,"\n                       ");
      }
   }
   snprintf(buffer,512," };\n");
   nch = strlen(buffer); f->write(buffer,nch);
   //array fB
   snprintf(buffer,512,"   const double fB[%d] = {",fNp);
   nch = strlen(buffer); f->write(buffer,nch);
   buffer[0] = 0;
   for (i=0;i<fNp;i++) {
      snprintf(numb,20," %g,",fPoly[i].B());
      nch = strlen(numb);
      if (i == fNp-1) numb[nch-1]=0;
      strlcat(buffer,numb,512);
      if (i%5 == 4 || i == fNp-1) {
         nch = strlen(buffer); f->write(buffer,nch);
         if (i != fNp-1) snprintf(buffer,512,"\n                       ");
      }
   }
   snprintf(buffer,512," };\n");
   nch = strlen(buffer); f->write(buffer,nch);
   //array fC
   snprintf(buffer,512,"   const double fC[%d] = {",fNp);
   nch = strlen(buffer); f->write(buffer,nch);
   buffer[0] = 0;
   for (i=0;i<fNp;i++) {
      snprintf(numb,20," %g,",fPoly[i].C());
      nch = strlen(numb);
      if (i == fNp-1) numb[nch-1]=0;
      strlcat(buffer,numb,512);
      if (i%5 == 4 || i == fNp-1) {
         nch = strlen(buffer); f->write(buffer,nch);
         if (i != fNp-1) snprintf(buffer,512,"\n                       ");
      }
   }
   snprintf(buffer,512," };\n");
   nch = strlen(buffer); f->write(buffer,nch);
    //array fD
   snprintf(buffer,512,"   const double fD[%d] = {",fNp);
   nch = strlen(buffer); f->write(buffer,nch);
   buffer[0] = 0;
   for (i=0;i<fNp;i++) {
      snprintf(numb,20," %g,",fPoly[i].D());
      nch = strlen(numb);
      if (i == fNp-1) numb[nch-1]=0;
      strlcat(buffer,numb,512);
      if (i%5 == 4 || i == fNp-1) {
         nch = strlen(buffer); f->write(buffer,nch);
         if (i != fNp-1) snprintf(buffer,512,"\n                       ");
      }
   }
   snprintf(buffer,512," };\n");
   nch = strlen(buffer); f->write(buffer,nch);
    //array fE
   snprintf(buffer,512,"   const double fE[%d] = {",fNp);
   nch = strlen(buffer); f->write(buffer,nch);
   buffer[0] = 0;
   for (i=0;i<fNp;i++) {
      snprintf(numb,20," %g,",fPoly[i].E());
      nch = strlen(numb);
      if (i == fNp-1) numb[nch-1]=0;
      strlcat(buffer,numb,512);
      if (i%5 == 4 || i == fNp-1) {
         nch = strlen(buffer); f->write(buffer,nch);
         if (i != fNp-1) snprintf(buffer,512,"\n                       ");
      }
   }
   snprintf(buffer,512," };\n");
   nch = strlen(buffer); f->write(buffer,nch);
    //array fF
   snprintf(buffer,512,"   const double fF[%d] = {",fNp);
   nch = strlen(buffer); f->write(buffer,nch);
   buffer[0] = 0;
   for (i=0;i<fNp;i++) {
      snprintf(numb,20," %g,",fPoly[i].F());
      nch = strlen(numb);
      if (i == fNp-1) numb[nch-1]=0;
      strlcat(buffer,numb,512);
      if (i%5 == 4 || i == fNp-1) {
         nch = strlen(buffer); f->write(buffer,nch);
         if (i != fNp-1) snprintf(buffer,512,"\n                       ");
      }
   }
   snprintf(buffer,512," };\n");
   nch = strlen(buffer); f->write(buffer,nch);

   //generate code for the spline evaluation
   snprintf(buffer,512,"   int klow=0;\n");
   nch = strlen(buffer); f->write(buffer,nch);

   snprintf(buffer,512,"   // If out of boundaries, extrapolate. It may be badly wrong\n");
   snprintf(buffer,512,"   if(x<=fXmin) klow=0;\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   else if(x>=fXmax) klow=fNp-1;\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   else {\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"     if(fKstep) {\n");
   nch = strlen(buffer); f->write(buffer,nch);

   snprintf(buffer,512,"       // Equidistant knots, use histogramming\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"       klow = int((x-fXmin)/fDelta);\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"       if (klow < fNp-1) klow = fNp-1;\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"     } else {\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"       int khig=fNp-1, khalf;\n");
   nch = strlen(buffer); f->write(buffer,nch);

   snprintf(buffer,512,"       // Non equidistant knots, binary search\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"       while(khig-klow>1)\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"         else khig=khalf;\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"     }\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   }\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   // Evaluate now\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   double dx=x-fX[klow];\n");
   nch = strlen(buffer); f->write(buffer,nch);
   snprintf(buffer,512,"   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*(fD[klow]+dx*(fE[klow]+dx*fF[klow])))));\n");
   nch = strlen(buffer); f->write(buffer,nch);

   //close file
   f->write("}\n",2);

   if (f) { f->close(); delete f;}
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out.

void TSpline5::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';
   out<<"   "<<std::endl;
   if (gROOT->ClassSaved(TSpline5::Class())) {
      out<<"   ";
   } else {
      out<<"   TSpline5 *";
   }
   Double_t b1 = fPoly[1].Y();
   Double_t e1 = fPoly[fNp-1].Y();
   Double_t b2 = fPoly[2].Y();
   Double_t e2 = fPoly[fNp-1].Y();
   out<<"spline5 = new TSpline5("<<quote<<GetTitle()<<quote<<","
      <<fXmin<<","<<fXmax<<",(TF1*)0,"<<fNp<<","<<quote<<quote<<","
      <<b1<<","<<e1<<","<<b2<<","<<e2<<");"<<std::endl;
   out<<"   spline5->SetName("<<quote<<GetName()<<quote<<");"<<std::endl;

   SaveFillAttributes(out,"spline5",0,1001);
   SaveLineAttributes(out,"spline5",1,1,1);
   SaveMarkerAttributes(out,"spline5",1,1,1);
   if (fNpx != 100) out<<"   spline5->SetNpx("<<fNpx<<");"<<std::endl;

   for (Int_t i=0;i<fNp;i++) {
      out<<"   spline5->SetPoint("<<i<<","<<fPoly[i].X()<<","<<fPoly[i].Y()<<");"<<std::endl;
      out<<"   spline5->SetPointCoeff("<<i<<","<<fPoly[i].B()<<","<<fPoly[i].C()<<","<<fPoly[i].D()<<","<<fPoly[i].E()<<","<<fPoly[i].F()<<");"<<std::endl;
   }
   out<<"   spline5->Draw("<<quote<<option<<quote<<");"<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Set point number i.

void TSpline5::SetPoint(Int_t i, Double_t x, Double_t y)
{

   if (i < 0 || i >= fNp) return;
   fPoly[i].X()= x;
   fPoly[i].Y()= y;
}

////////////////////////////////////////////////////////////////////////////////
/// Set point coefficient number i.

void TSpline5::SetPointCoeff(Int_t i, Double_t b, Double_t c, Double_t d,
                             Double_t e, Double_t f)
{
   if (i < 0 || i >= fNp) return;
   fPoly[i].B()= b;
   fPoly[i].C()= c;
   fPoly[i].D()= d;
   fPoly[i].E()= e;
   fPoly[i].F()= f;
}

////////////////////////////////////////////////////////////////////////////////
/// Algorithm 600, collected algorithms from acm.
///
/// algorithm appeared in acm-trans. math. software, vol.9, no. 2,
/// jun., 1983, p. 258-259.
///
///     TSpline5 computes the coefficients of a quintic natural quintic spli
///     s(x) with knots x(i) interpolating there to given function values:
/// ~~~ {.cpp}
///               s(x(i)) = y(i)  for i = 1,2, ..., n.
/// ~~~
///     in each interval (x(i),x(i+1)) the spline function s(xx) is a
///     polynomial of fifth degree:
/// ~~~ {.cpp}
///     s(xx) = ((((f(i)*p+e(i))*p+d(i))*p+c(i))*p+b(i))*p+y(i)    (*)
///           = ((((-f(i)*q+e(i+1))*q-d(i+1))*q+c(i+1))*q-b(i+1))*q+y(i+1)
/// ~~~
///     where  p = xx - x(i)  and  q = x(i+1) - xx.
///     (note the first subscript in the second expression.)
///     the different polynomials are pieced together so that s(x) and
///     its derivatives up to s"" are continuous.
///
/// ### input:
///
///     n          number of data points, (at least three, i.e. n > 2)
///     x(1:n)     the strictly increasing or decreasing sequence of
///                knots.  the spacing must be such that the fifth power
///                of x(i+1) - x(i) can be formed without overflow or
///                underflow of exponents.
///     y(1:n)     the prescribed function values at the knots.
///
/// ### output:
///
///     b,c,d,e,f  the computed spline coefficients as in (*).
///         (1:n)  specifically
///                b(i) = s'(x(i)), c(i) = s"(x(i))/2, d(i) = s"'(x(i))/6,
///                e(i) = s""(x(i))/24,  f(i) = s""'(x(i))/120.
///                f(n) is neither used nor altered.  the five arrays
///                b,c,d,e,f must always be distinct.
///
/// ### option:
///
///     it is possible to specify values for the first and second
///     derivatives of the spline function at arbitrarily many knots.
///     this is done by relaxing the requirement that the sequence of
///     knots be strictly increasing or decreasing.  specifically:
///
/// ~~~ {.cpp}
///     if x(j) = x(j+1) then s(x(j)) = y(j) and s'(x(j)) = y(j+1),
///     if x(j) = x(j+1) = x(j+2) then in addition s"(x(j)) = y(j+2).
/// ~~~
///
///     note that s""(x) is discontinuous at a double knot and, in
///     addition, s"'(x) is discontinuous at a triple knot.  the
///     subroutine assigns y(i) to y(i+1) in these cases and also to
///     y(i+2) at a triple knot.  the representation (*) remains
///     valid in each open interval (x(i),x(i+1)).  at a double knot,
///     x(j) = x(j+1), the output coefficients have the following values:
/// ~~~ {.cpp}
///       y(j) = s(x(j))          = y(j+1)
///       b(j) = s'(x(j))         = b(j+1)
///       c(j) = s"(x(j))/2       = c(j+1)
///       d(j) = s"'(x(j))/6      = d(j+1)
///       e(j) = s""(x(j)-0)/24     e(j+1) = s""(x(j)+0)/24
///       f(j) = s""'(x(j)-0)/120   f(j+1) = s""'(x(j)+0)/120
/// ~~~
///     at a triple knot, x(j) = x(j+1) = x(j+2), the output
///     coefficients have the following values:
/// ~~~ {.cpp}
///       y(j) = s(x(j))         = y(j+1)    = y(j+2)
///       b(j) = s'(x(j))        = b(j+1)    = b(j+2)
///       c(j) = s"(x(j))/2      = c(j+1)    = c(j+2)
///       d(j) = s"'((x(j)-0)/6    d(j+1) = 0  d(j+2) = s"'(x(j)+0)/6
///       e(j) = s""(x(j)-0)/24    e(j+1) = 0  e(j+2) = s""(x(j)+0)/24
///       f(j) = s""'(x(j)-0)/120  f(j+1) = 0  f(j+2) = s""'(x(j)+0)/120
/// ~~~

void TSpline5::BuildCoeff()
{
   Int_t i, m;
   Double_t pqqr, p, q, r, s, t, u, v,
      b1, p2, p3, q2, q3, r2, pq, pr, qr;

   if (fNp <= 2) {
      return;
   }

   //     coefficients of a positive definite, pentadiagonal matrix,
   //     stored in D, E, F from 1 to n-3.
   m = fNp-2;
   q = fPoly[1].X()-fPoly[0].X();
   r = fPoly[2].X()-fPoly[1].X();
   q2 = q*q;
   r2 = r*r;
   qr = q+r;
   fPoly[0].D() = fPoly[0].E() = 0;
   if (q) fPoly[1].D() = q*6.*q2/(qr*qr);
   else fPoly[1].D() = 0;

   if (m > 1) {
      for (i = 1; i < m; ++i) {
         p = q;
         q = r;
         r = fPoly[i+2].X()-fPoly[i+1].X();
         p2 = q2;
         q2 = r2;
         r2 = r*r;
         pq = qr;
         qr = q+r;
         if (q) {
            q3 = q2*q;
            pr = p*r;
            pqqr = pq*qr;
            fPoly[i+1].D() = q3*6./(qr*qr);
            fPoly[i].D() += (q+q)*(pr*15.*pr+(p+r)*q
                                 *(pr* 20.+q2*7.)+q2*
                                 ((p2+r2)*8.+pr*21.+q2+q2))/(pqqr*pqqr);
            fPoly[i-1].D() += q3*6./(pq*pq);
            fPoly[i].E() = q2*(p*qr+pq*3.*(qr+r+r))/(pqqr*qr);
            fPoly[i-1].E() += q2*(r*pq+qr*3.*(pq+p+p))/(pqqr*pq);
            fPoly[i-1].F() = q3/pqqr;
         } else
            fPoly[i+1].D() = fPoly[i].E() = fPoly[i-1].F() = 0;
      }
   }
   if (r) fPoly[m-1].D() += r*6.*r2/(qr*qr);

   //     First and second order divided differences of the given function
   //     values, stored in b from 2 to n and in c from 3 to n
   //     respectively. care is taken of double and triple knots.
   for (i = 1; i < fNp; ++i) {
      if (fPoly[i].X() != fPoly[i-1].X()) {
         fPoly[i].B() =
            (fPoly[i].Y()-fPoly[i-1].Y())/(fPoly[i].X()-fPoly[i-1].X());
      } else {
         fPoly[i].B() = fPoly[i].Y();
         fPoly[i].Y() = fPoly[i-1].Y();
      }
   }
   for (i = 2; i < fNp; ++i) {
      if (fPoly[i].X() != fPoly[i-2].X()) {
         fPoly[i].C() =
            (fPoly[i].B()-fPoly[i-1].B())/(fPoly[i].X()-fPoly[i-2].X());
      } else {
         fPoly[i].C() = fPoly[i].B()*.5;
         fPoly[i].B() = fPoly[i-1].B();
      }
   }

   //     Solve the linear system with c(i+2) - c(i+1) as right-hand side. */
   if (m > 1) {
      p=fPoly[0].C()=fPoly[m-1].E()=fPoly[0].F()
         =fPoly[m-2].F()=fPoly[m-1].F()=0;
      fPoly[1].C() = fPoly[3].C()-fPoly[2].C();
      fPoly[1].D() = 1./fPoly[1].D();

      if (m > 2) {
         for (i = 2; i < m; ++i) {
            q = fPoly[i-1].D()*fPoly[i-1].E();
            fPoly[i].D() = 1./(fPoly[i].D()-p*fPoly[i-2].F()-q*fPoly[i-1].E());
            fPoly[i].E() -= q*fPoly[i-1].F();
            fPoly[i].C() = fPoly[i+2].C()-fPoly[i+1].C()-p*fPoly[i-2].C()
                           -q*fPoly[i-1].C();
            p = fPoly[i-1].D()*fPoly[i-1].F();
         }
      }
   }

   fPoly[fNp-2].C() = fPoly[fNp-1].C() = 0;
   if (fNp > 3)
      for (i=fNp-3; i > 0; --i)
         fPoly[i].C() = (fPoly[i].C()-fPoly[i].E()*fPoly[i+1].C()
                        -fPoly[i].F()*fPoly[i+2].C())*fPoly[i].D();

   //     Integrate the third derivative of s(x)
   m = fNp-1;
   q = fPoly[1].X()-fPoly[0].X();
   r = fPoly[2].X()-fPoly[1].X();
   b1 = fPoly[1].B();
   q3 = q*q*q;
   qr = q+r;
   if (qr) {
      v = fPoly[1].C()/qr;
      t = v;
   } else
      v = t = 0;
   if (q) fPoly[0].F() = v/q;
   else fPoly[0].F() = 0;
   for (i = 1; i < m; ++i) {
      p = q;
      q = r;
      if (i != m-1) r = fPoly[i+2].X()-fPoly[i+1].X();
      else r = 0;
      p3 = q3;
      q3 = q*q*q;
      pq = qr;
      qr = q+r;
      s = t;
      if (qr) t = (fPoly[i+1].C()-fPoly[i].C())/qr;
      else t = 0;
      u = v;
      v = t-s;
      if (pq) {
         fPoly[i].F() = fPoly[i-1].F();
         if (q) fPoly[i].F() = v/q;
         fPoly[i].E() = s*5.;
         fPoly[i].D() = (fPoly[i].C()-q*s)*10;
         fPoly[i].C() =
         fPoly[i].D()*(p-q)+(fPoly[i+1].B()-fPoly[i].B()+(u-fPoly[i].E())*
                            p3-(v+fPoly[i].E())*q3)/pq;
         fPoly[i].B() = (p*(fPoly[i+1].B()-v*q3)+q*(fPoly[i].B()-u*p3))/pq-p
         *q*(fPoly[i].D()+fPoly[i].E()*(q-p));
      } else {
         fPoly[i].C() = fPoly[i-1].C();
         fPoly[i].D() = fPoly[i].E() = fPoly[i].F() = 0;
      }
   }

   //     End points x(1) and x(n)
   p = fPoly[1].X()-fPoly[0].X();
   s = fPoly[0].F()*p*p*p;
   fPoly[0].E() = fPoly[0].D() = 0;
   fPoly[0].C() = fPoly[1].C()-s*10;
   fPoly[0].B() = b1-(fPoly[0].C()+s)*p;

   q = fPoly[fNp-1].X()-fPoly[fNp-2].X();
   t = fPoly[fNp-2].F()*q*q*q;
   fPoly[fNp-1].E() = fPoly[fNp-1].D() = 0;
   fPoly[fNp-1].C() = fPoly[fNp-2].C()+t*10;
   fPoly[fNp-1].B() += (fPoly[fNp-1].C()-t)*q;
}

////////////////////////////////////////////////////////////////////////////////
/// Test method for TSpline5
///
/// ~~~ {.cpp}
///   n          number of data points.
///   m          2*m-1 is order of spline.
///                 m = 3 always for quintic spline.
///   nn,nm1,mm,
///   mm1,i,k,
///   j,jj       temporary integer variables.
///   z,p        temporary double precision variables.
///   x[n]       the sequence of knots.
///   y[n]       the prescribed function values at the knots.
///   a[200][6]  two dimensional array whose columns are
///                 the computed spline coefficients
///   diff[5]    maximum values of differences of values and
///                 derivatives to right and left of knots.
///   com[5]     maximum values of coefficients.
/// ~~~
///
///   test of TSpline5 with non equidistant knots and
///      equidistant knots follows.

void TSpline5::Test()
{
   Double_t hx;
   Double_t diff[5];
   Double_t a[1200], c[6];
   Int_t i, j, k, m, n;
   Double_t p, x[200], y[200], z;
   Int_t jj, mm, nn;
   Int_t mm1, nm1;
   Double_t com[5];

   printf("1         TEST OF TSpline5 WITH NONEQUIDISTANT KNOTS\n");
   n = 5;
   x[0] = -3;
   x[1] = -1;
   x[2] = 0;
   x[3] = 3;
   x[4] = 4;
   y[0] = 7;
   y[1] = 11;
   y[2] = 26;
   y[3] = 56;
   y[4] = 29;
   m = 3;
   mm = m << 1;
   mm1 = mm-1;
   printf("\n-N = %3d    M =%2d\n",n,m);
   TSpline5 *spline = new TSpline5("Test",x,y,n);
   for (i = 0; i < n; ++i)
      spline->GetCoeff(i,hx, a[i],a[i+200],a[i+400],
                       a[i+600],a[i+800],a[i+1000]);
   delete spline;
   for (i = 0; i < mm1; ++i) diff[i] = com[i] = 0;
   for (k = 0; k < n; ++k) {
      for (i = 0; i < mm; ++i) c[i] = a[k+i*200];
      printf(" ---------------------------------------%3d --------------------------------------------\n",k+1);
      printf("%12.8f\n",x[k]);
      if (k == n-1) {
         printf("%16.8f\n",c[0]);
      } else {
         for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
         printf("\n");
         for (i = 0; i < mm1; ++i)
            if ((z=TMath::Abs(a[k+i*200])) > com[i]) com[i] = z;
         z = x[k+1]-x[k];
         for (i = 1; i < mm; ++i)
            for (jj = i; jj < mm; ++jj) {
               j = mm+i-jj;
               c[j-2] = c[j-1]*z+c[j-2];
            }
         for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
         printf("\n");
         for (i = 0; i < mm1; ++i)
            if (!(k >= n-2 && i != 0))
               if((z = TMath::Abs(c[i]-a[k+1+i*200]))
                  > diff[i]) diff[i] = z;
      }
   }
   printf("  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES \n");
   for (i = 0; i < mm1; ++i) printf("%18.9E",diff[i]);
   printf("\n");
   printf("  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS \n");
   if (TMath::Abs(c[0]) > com[0])
      com[0] = TMath::Abs(c[0]);
   for (i = 0; i < mm1; ++i) printf("%16.8f",com[i]);
   printf("\n");
   m = 3;
   for (n = 10; n <= 100; n += 10) {
      mm = m << 1;
      mm1 = mm-1;
      nm1 = n-1;
      for (i = 0; i < nm1; i += 2) {
         x[i] = i+1;
         x[i+1] = i+2;
         y[i] = 1;
         y[i+1] = 0;
      }
      if (n % 2 != 0) {
         x[n-1] = n;
         y[n-1] = 1;
      }
      printf("\n-N = %3d    M =%2d\n",n,m);
      spline = new TSpline5("Test",x,y,n);
      for (i = 0; i < n; ++i)
         spline->GetCoeff(i,hx,a[i],a[i+200],a[i+400],
                          a[i+600],a[i+800],a[i+1000]);
      delete spline;
      for (i = 0; i < mm1; ++i)
         diff[i] = com[i] = 0;
      for (k = 0; k < n; ++k) {
         for (i = 0; i < mm; ++i)
            c[i] = a[k+i*200];
         if (n < 11) {
            printf(" ---------------------------------------%3d --------------------------------------------\n",k+1);
            printf("%12.8f\n",x[k]);
            if (k == n-1) printf("%16.8f\n",c[0]);
         }
         if (k == n-1) break;
         if (n <= 10) {
            for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
            printf("\n");
         }
         for (i = 0; i < mm1; ++i)
            if ((z=TMath::Abs(a[k+i*200])) > com[i])
               com[i] = z;
         z = x[k+1]-x[k];
         for (i = 1; i < mm; ++i)
            for (jj = i; jj < mm; ++jj) {
               j = mm+i-jj;
               c[j-2] = c[j-1]*z+c[j-2];
            }
         if (n <= 10) {
            for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
            printf("\n");
         }
         for (i = 0; i < mm1; ++i)
            if (!(k >= n-2 && i != 0))
               if ((z = TMath::Abs(c[i]-a[k+1+i*200]))
                  > diff[i]) diff[i] = z;
      }
      printf("  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES \n");
      for (i = 0; i < mm1; ++i) printf("%18.9E",diff[i]);
      printf("\n");
      printf("  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS \n");
      if (TMath::Abs(c[0]) > com[0])
         com[0] = TMath::Abs(c[0]);
      for (i = 0; i < mm1; ++i) printf("%16.8E",com[i]);
      printf("\n");
   }

   //     Test of TSpline5 with non equidistant double knots follows
   printf("1  TEST OF TSpline5 WITH NONEQUIDISTANT DOUBLE KNOTS\n");
   n = 5;
   x[0] = -3;
   x[1] = -3;
   x[2] = -1;
   x[3] = -1;
   x[4] = 0;
   x[5] = 0;
   x[6] = 3;
   x[7] = 3;
   x[8] = 4;
   x[9] = 4;
   y[0] = 7;
   y[1] = 2;
   y[2] = 11;
   y[3] = 15;
   y[4] = 26;
   y[5] = 10;
   y[6] = 56;
   y[7] = -27;
   y[8] = 29;
   y[9] = -30;
   m = 3;
   nn = n << 1;
   mm = m << 1;
   mm1 = mm-1;
   printf("-N = %3d    M =%2d\n",n,m);
   spline = new TSpline5("Test",x,y,nn);
   for (i = 0; i < nn; ++i)
      spline->GetCoeff(i,hx,a[i],a[i+200],a[i+400],
                       a[i+600],a[i+800],a[i+1000]);
   delete spline;
   for (i = 0; i < mm1; ++i)
      diff[i] = com[i] = 0;
   for (k = 0; k < nn; ++k) {
      for (i = 0; i < mm; ++i)
         c[i] = a[k+i*200];
      printf(" ---------------------------------------%3d --------------------------------------------\n",k+1);
      printf("%12.8f\n",x[k]);
      if (k == nn-1) {
         printf("%16.8f\n",c[0]);
         break;
      }
      for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
      printf("\n");
      for (i = 0; i < mm1; ++i)
         if ((z=TMath::Abs(a[k+i*200])) > com[i]) com[i] = z;
      z = x[k+1]-x[k];
      for (i = 1; i < mm; ++i)
         for (jj = i; jj < mm; ++jj) {
            j = mm+i-jj;
            c[j-2] = c[j-1]*z+c[j-2];
         }
      for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
      printf("\n");
      for (i = 0; i < mm1; ++i)
         if (!(k >= nn-2 && i != 0))
            if ((z = TMath::Abs(c[i]-a[k+1+i*200]))
               > diff[i]) diff[i] = z;
   }
   printf("  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES \n");
   for (i = 1; i <= mm1; ++i) {
      printf("%18.9E",diff[i-1]);
   }
   printf("\n");
   if (TMath::Abs(c[0]) > com[0])
      com[0] = TMath::Abs(c[0]);
   printf("  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS \n");
   for (i = 0; i < mm1; ++i) printf("%16.8f",com[i]);
   printf("\n");
   m = 3;
   for (n = 10; n <= 100; n += 10) {
      nn = n << 1;
      mm = m << 1;
      mm1 = mm-1;
      p = 0;
      for (i = 0; i < n; ++i) {
         p += TMath::Abs(TMath::Sin(i+1));
         x[(i << 1)] = p;
         x[(i << 1)+1] = p;
         y[(i << 1)] = TMath::Cos(i+1)-.5;
         y[(i << 1)+1] = TMath::Cos((i << 1)+2)-.5;
      }
      printf("-N = %3d    M =%2d\n",n,m);
      spline = new TSpline5("Test",x,y,nn);
      for (i = 0; i < nn; ++i)
         spline->GetCoeff(i,hx,a[i],a[i+200],a[i+400],
                          a[i+600],a[i+800],a[i+1000]);
      delete spline;
      for (i = 0; i < mm1; ++i)
         diff[i] = com[i] = 0;
      for (k = 0; k < nn; ++k) {
         for (i = 0; i < mm; ++i)
            c[i] = a[k+i*200];
         if (n < 11) {
            printf(" ---------------------------------------%3d --------------------------------------------\n",k+1);
            printf("%12.8f\n",x[k]);
            if (k == nn-1) printf("%16.8f\n",c[0]);
         }
         if (k == nn-1) break;
         if (n <= 10) {
            for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
            printf("\n");
         }
         for (i = 0; i < mm1; ++i)
            if ((z=TMath::Abs(a[k+i*200])) > com[i]) com[i] = z;
         z = x[k+1]-x[k];
         for (i = 1; i < mm; ++i) {
            for (jj = i; jj < mm; ++jj) {
               j = mm+i-jj;
               c[j-2] = c[j-1]*z+c[j-2];
            }
         }
         if (n <= 10) {
            for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
            printf("\n");
         }
         for (i = 0; i < mm1; ++i)
            if (!(k >= nn-2 && i != 0))
               if ((z = TMath::Abs(c[i]-a[k+1+i*200]))
                  > diff[i]) diff[i] = z;
      }
      printf("  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES \n");
      for (i = 0; i < mm1; ++i) printf("%18.9E",diff[i]);
      printf("\n");
      printf("  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS \n");
      if (TMath::Abs(c[0]) > com[0])
         com[0] = TMath::Abs(c[0]);
      for (i = 0; i < mm1; ++i) printf("%18.9E",com[i]);
      printf("\n");
   }

   //     test of TSpline5 with non equidistant knots, one double knot,
   //        one triple knot, follows.
   printf("1         TEST OF TSpline5 WITH NONEQUIDISTANT KNOTS,\n");
   printf("             ONE DOUBLE, ONE TRIPLE KNOT\n");
   n = 8;
   x[0] = -3;
   x[1] = -1;
   x[2] = -1;
   x[3] = 0;
   x[4] = 3;
   x[5] = 3;
   x[6] = 3;
   x[7] = 4;
   y[0] = 7;
   y[1] = 11;
   y[2] = 15;
   y[3] = 26;
   y[4] = 56;
   y[5] = -30;
   y[6] = -7;
   y[7] = 29;
   m = 3;
   mm = m << 1;
   mm1 = mm-1;
   printf("-N = %3d    M =%2d\n",n,m);
   spline=new TSpline5("Test",x,y,n);
   for (i = 0; i < n; ++i)
      spline->GetCoeff(i,hx,a[i],a[i+200],a[i+400],
                       a[i+600],a[i+800],a[i+1000]);
   delete spline;
   for (i = 0; i < mm1; ++i)
      diff[i] = com[i] = 0;
   for (k = 0; k < n; ++k) {
      for (i = 0; i < mm; ++i)
         c[i] = a[k+i*200];
      printf(" ---------------------------------------%3d --------------------------------------------\n",k+1);
      printf("%12.8f\n",x[k]);
      if (k == n-1) {
         printf("%16.8f\n",c[0]);
         break;
      }
      for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
      printf("\n");
      for (i = 0; i < mm1; ++i)
         if ((z=TMath::Abs(a[k+i*200])) > com[i]) com[i] = z;
      z = x[k+1]-x[k];
      for (i = 1; i < mm; ++i)
         for (jj = i; jj < mm; ++jj) {
            j = mm+i-jj;
            c[j-2] = c[j-1]*z+c[j-2];
         }
      for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
      printf("\n");
      for (i = 0; i < mm1; ++i)
         if (!(k >= n-2 && i != 0))
            if ((z = TMath::Abs(c[i]-a[k+1+i*200]))
               > diff[i]) diff[i] = z;
   }
   printf("  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES \n");
   for (i = 0; i < mm1; ++i) printf("%18.9E",diff[i]);
   printf("\n");
   printf("  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS \n");
   if (TMath::Abs(c[0]) > com[0])
      com[0] = TMath::Abs(c[0]);
   for (i = 0; i < mm1; ++i) printf("%16.8f",com[i]);
   printf("\n");

   //     Test of TSpline5 with non equidistant knots, two double knots,
   //        one triple knot,follows.
   printf("1         TEST OF TSpline5 WITH NONEQUIDISTANT KNOTS,\n");
   printf("             TWO DOUBLE, ONE TRIPLE KNOT\n");
   n = 10;
   x[0] = 0;
   x[1] = 2;
   x[2] = 2;
   x[3] = 3;
   x[4] = 3;
   x[5] = 3;
   x[6] = 5;
   x[7] = 8;
   x[8] = 9;
   x[9] = 9;
   y[0] = 163;
   y[1] = 237;
   y[2] = -127;
   y[3] = 119;
   y[4] = -65;
   y[5] = 192;
   y[6] = 293;
   y[7] = 326;
   y[8] = 0;
   y[9] = -414;
   m = 3;
   mm = m << 1;
   mm1 = mm-1;
   printf("-N = %3d    M =%2d\n",n,m);
   spline = new TSpline5("Test",x,y,n);
   for (i = 0; i < n; ++i)
      spline->GetCoeff(i,hx,a[i],a[i+200],a[i+400],
                       a[i+600],a[i+800],a[i+1000]);
   delete spline;
   for (i = 0; i < mm1; ++i)
      diff[i] = com[i] = 0;
   for (k = 0; k < n; ++k) {
      for (i = 0; i < mm; ++i)
         c[i] = a[k+i*200];
      printf(" ---------------------------------------%3d --------------------------------------------\n",k+1);
      printf("%12.8f\n",x[k]);
      if (k == n-1) {
         printf("%16.8f\n",c[0]);
         break;
      }
      for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
      printf("\n");
      for (i = 0; i < mm1; ++i)
         if ((z=TMath::Abs(a[k+i*200])) > com[i]) com[i] = z;
      z = x[k+1]-x[k];
      for (i = 1; i < mm; ++i)
         for (jj = i; jj < mm; ++jj) {
            j = mm+i-jj;
            c[j-2] = c[j-1]*z+c[j-2];
         }
      for (i = 0; i < mm; ++i) printf("%16.8f",c[i]);
      printf("\n");
      for (i = 0; i < mm1; ++i)
         if (!(k >= n-2 && i != 0))
            if((z = TMath::Abs(c[i]-a[k+1+i*200]))
               > diff[i]) diff[i] = z;
   }
   printf("  MAXIMUM ABSOLUTE VALUES OF DIFFERENCES \n");
   for (i = 0; i < mm1; ++i) printf("%18.9E",diff[i]);
   printf("\n");
   printf("  MAXIMUM ABSOLUTE VALUES OF COEFFICIENTS \n");
   if (TMath::Abs(c[0]) > com[0])
      com[0] = TMath::Abs(c[0]);
   for (i = 0; i < mm1; ++i) printf("%16.8f",com[i]);
   printf("\n");
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TSpline5.

void TSpline5::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TSpline5::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TSpline::Streamer(R__b);
      if (fNp > 0) {
         fPoly = new TSplinePoly5[fNp];
         for(Int_t i=0; i<fNp; ++i) {
            fPoly[i].Streamer(R__b);
         }
      }
      //      R__b >> fPoly;
   } else {
      R__b.WriteClassBuffer(TSpline5::Class(),this);
   }
}
