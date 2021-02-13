// @(#)root/g3d:$Id$
// Author: Nenad Buncic   21/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TView.h"
#include "TPolyMarker3D.h"
#include "TVirtualPad.h"
#include "TRandom.h"
#include "TBuffer.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TVirtualViewer3D.h"
#include "TGeometry.h"
#include "TH1.h"
#include "TROOT.h"
#include "TMath.h"

#include <cassert>
#include <iostream>

ClassImp(TPolyMarker3D);

constexpr Int_t kDimension = 3;

/** \class TPolyMarker3D
\ingroup g3d
A 3D polymarker.

It has three constructors.

First one, without any parameters TPolyMarker3D(), we call 'default
constructor' and it's used in a case that just an initialisation is
needed (i.e. pointer declaration).

Example:

~~~ {.cpp}
   TPolyMarker3D *pm = new TPolyMarker3D;
~~~

Second one, takes, usually, two parameters, n (number of points) and
marker (marker style). Third parameter is optional.

Example:

~~~ {.cpp}
   TPolyMarker3D (150, 1);
~~~

Third one takes, usually, three parameters, n (number of points), *p
(pointer to an array of 3D points), and marker (marker style). Fourth
parameter is optional.

Example:

~~~ {.cpp}
   Float_t *ptr = new Float_t [150*3];
      ... ... ...
      ... ... ...
      ... ... ...
   TPolyMarker3D (150, ptr, 1);
~~~
*/

////////////////////////////////////////////////////////////////////////////////
/// 3-D polymarker default constructor.

TPolyMarker3D::TPolyMarker3D()
{
   fN = 0;
   fP = 0;
   fLastPoint = -1;
   fName = "TPolyMarker3D";
}

////////////////////////////////////////////////////////////////////////////////
/// 3-D polymarker normal constructor with initialization to 0.

TPolyMarker3D::TPolyMarker3D(Int_t n, Marker_t marker, Option_t *option)
{
   fName = "TPolyMarker3D";
   fOption = option;
   SetMarkerStyle(marker);
   SetBit(kCanDelete);
   fLastPoint = -1;
   if (n <= 0) {
      fN = 0;
      fP = 0;
      return;
   }

   fN = n;
   fP = new Float_t [kDimension*fN];
   for (Int_t i = 0; i < kDimension*fN; i++)  fP[i] = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// 3-D polymarker constructor. Polymarker is initialized with p.

TPolyMarker3D::TPolyMarker3D(Int_t n, Float_t *p, Marker_t marker,
                             Option_t *option)
{
   fName = "TPolyMarker3D";
   SetMarkerStyle(marker);
   SetBit(kCanDelete);
   fOption = option;
   fLastPoint = -1;
   if (n <= 0) {
      fN = 0;
      fP = 0;
      return;
   }

   fN = n;
   fP = new Float_t [kDimension*fN];
   if (p) {
      for (Int_t i = 0; i < kDimension*fN; i++)
         fP[i] = p[i];
      fLastPoint = fN-1;
   } else
      memset(fP,0,kDimension*fN*sizeof(Float_t));
}

////////////////////////////////////////////////////////////////////////////////
/// 3-D polymarker constructor. Polymarker is initialized with p
/// (cast to float).

TPolyMarker3D::TPolyMarker3D(Int_t n, Double_t *p, Marker_t marker,
                             Option_t *option)
{
   fName = "TPolyMarker3D";
   SetMarkerStyle(marker);
   SetBit(kCanDelete);
   fOption = option;
   fLastPoint = -1;
   if (n <= 0) {
      fN = 0;
      fP = 0;
      return;
   }

   fN = n;
   fP = new Float_t [kDimension*fN];
   if (p) {
      for (Int_t i = 0; i < kDimension*fN; i++)
         fP[i] = (Float_t) p[i];
      fLastPoint = fN-1;
   } else
      memset(fP,0,kDimension*fN*sizeof(Float_t));
}

////////////////////////////////////////////////////////////////////////////////
/// assignment operator

TPolyMarker3D& TPolyMarker3D::operator=(const TPolyMarker3D& tp3)
{
   if(this!=&tp3) {
      TObject::operator=(tp3);
      TAttMarker::operator=(tp3);
      TAtt3D::operator=(tp3);
      fN=tp3.fN;
      fP=tp3.fP;
      fOption=tp3.fOption;
      fLastPoint=tp3.fLastPoint;
      fName=tp3.fName;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// 3-D polymarker destructor.

TPolyMarker3D::~TPolyMarker3D()
{
   fN = 0;
   if (fP) delete [] fP;
   fLastPoint = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// 3-D polymarker copy ctor.

TPolyMarker3D::TPolyMarker3D(const TPolyMarker3D &p) :
   TObject(p), TAttMarker(p), TAtt3D(p)
{
   fN = 0;
   fP = 0;
   fLastPoint = -1;
   p.Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy polymarker to polymarker obj.

void TPolyMarker3D::Copy(TObject &obj) const
{
   TObject::Copy(obj);
   ((TPolyMarker3D&)obj).fN = fN;
   if (fN > 0) {
      ((TPolyMarker3D&)obj).fP = new Float_t [kDimension*fN];
      for (Int_t i = 0; i < kDimension*fN; i++)  ((TPolyMarker3D&)obj).fP[i] = fP[i];
   } else {
      ((TPolyMarker3D&)obj).fP = 0;
   }
   ((TPolyMarker3D&)obj).SetMarkerStyle(GetMarkerStyle());
   ((TPolyMarker3D&)obj).fOption = fOption;
   ((TPolyMarker3D&)obj).fLastPoint = fLastPoint;
   ((TPolyMarker3D&)obj).fName   = fName;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a 3-D polymarker.
/// Compute the closest distance of approach from point px,py to each segment
/// of the polymarker.
/// Returns when the distance found is below DistanceMaximum.
/// The distance is computed in pixels units.

Int_t TPolyMarker3D::DistancetoPrimitive(Int_t px, Int_t py)
{
   const Int_t inaxis = 7;
   Int_t dist = 9999;

   Int_t puxmin = gPad->XtoAbsPixel(gPad->GetUxmin());
   Int_t puymin = gPad->YtoAbsPixel(gPad->GetUymin());
   Int_t puxmax = gPad->XtoAbsPixel(gPad->GetUxmax());
   Int_t puymax = gPad->YtoAbsPixel(gPad->GetUymax());

   // return if point is not in the user area
   if (px < puxmin - inaxis) return dist;
   if (py > puymin + inaxis) return dist;
   if (px > puxmax + inaxis) return dist;
   if (py < puymax - inaxis) return dist;

   TView *view = gPad->GetView();
   if (!view) return dist;
   Int_t i, dpoint;
   Float_t xndc[3];
   Int_t x1,y1;
   Double_t u,v;
   for (i=0;i<Size();i++) {
      view->WCtoNDC(&fP[3*i], xndc);
      u      = (Double_t)xndc[0];
      v      = (Double_t)xndc[1];
      if (u < gPad->GetUxmin() || u > gPad->GetUxmax()) continue;
      if (v < gPad->GetUymin() || v > gPad->GetUymax()) continue;
      x1     = gPad->XtoAbsPixel(u);
      y1     = gPad->YtoAbsPixel(v);
      dpoint = Int_t(TMath::Sqrt((((Double_t)px-x1)*((Double_t)px-x1)
                                + ((Double_t)py-y1)*((Double_t)py-y1))));
      if (dpoint < dist) dist = dpoint;
   }
   return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// Draws 3-D polymarker with its current attributes.

void TPolyMarker3D::Draw(Option_t *option)
{
   AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this 3-D polymarker with new coordinates. Creates a new
/// polymarker which will be adopted by the pad in which it is drawn.
/// Does not change the original polymarker (should be static method).

void TPolyMarker3D::DrawPolyMarker(Int_t n, Float_t *p, Marker_t, Option_t *option)
{
   TPolyMarker3D *newpolymarker = new TPolyMarker3D();
   newpolymarker->fN = n;
   newpolymarker->fP = new Float_t [kDimension*fN];
   for (Int_t i = 0; i < kDimension*fN; i++)  newpolymarker->fP[i] = p[i];
   newpolymarker->SetMarkerStyle(GetMarkerStyle());
   newpolymarker->fOption = fOption;
   newpolymarker->fLastPoint = fLastPoint;
   newpolymarker->SetBit(kCanDelete);
   newpolymarker->AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.

void TPolyMarker3D::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (!gPad) return;
   if (gPad->GetView()) gPad->GetView()->ExecuteRotateView(event, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// List this 3-D polymarker.

void TPolyMarker3D::ls(Option_t *option) const
{
   TROOT::IndentLevel();
   std::cout << "    TPolyMarker3D  N=" << Size() <<" Option="<<option<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Merge polymarkers in the collection in this polymarker

Int_t TPolyMarker3D::Merge(TCollection *li)
{
   if (!li) return 0;
   TIter next(li);

   //first loop to count the number of entries
   TPolyMarker3D *pm;
   Int_t npoints = Size();
   while ((pm = (TPolyMarker3D*)next())) {
      if (!pm->InheritsFrom(TPolyMarker3D::Class())) {
         Error("Add","Attempt to add object of class: %s to a %s",pm->ClassName(),this->ClassName());
         return -1;
      }
      npoints += pm->Size();
   }
   Int_t currPoint = Size();

   //extend this polymarker to hold npoints
   SetPoint(npoints-1,0,0,0);

   //merge all polymarkers
   next.Reset();
   while ((pm = (TPolyMarker3D*)next())) {
      Int_t np = pm->Size();
      Float_t *p = pm->GetP();
      for (Int_t i = 0; i < np; i++) {
         SetPoint(currPoint++, p[3*i], p[3*i+1], p[3*i+2]);
      }
   }
   return npoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint a TPolyMarker3D.

void TPolyMarker3D::Paint(Option_t * /*option*/ )
{
   // No need to continue if there is nothing to paint
   if (Size() <= 0) return;

   static TBuffer3D buffer(TBuffer3DTypes::kMarker);

   buffer.ClearSectionsValid();

   // Section kCore
   buffer.fID           = this;
   buffer.fColor        = GetMarkerColor();
   buffer.fTransparency = 0;
   buffer.fLocalFrame   = kFALSE;
   buffer.SetSectionsValid(TBuffer3D::kCore);

   // We fill kCore and kRawSizes on first pass and try with viewer
   TVirtualViewer3D * viewer3D = gPad->GetViewer3D();
   if (!viewer3D) return;
   Int_t reqSections = viewer3D->AddObject(buffer);
   if (reqSections == TBuffer3D::kNone) {
      return;
   }

   if (reqSections & TBuffer3D::kRawSizes) {
      if (!buffer.SetRawSizes(Size(), 3*Size(), 1, 1, 0, 0)) {
         return;
      }
      buffer.SetSectionsValid(TBuffer3D::kRawSizes);
   }

   if ((reqSections & TBuffer3D::kRaw) && buffer.SectionsValid(TBuffer3D::kRawSizes)) {
      // Points
      for (UInt_t i=0; i<3*buffer.NbPnts(); i++) {
         buffer.fPnts[i] = (Double_t)fP[i];
      }

      // Transform points - we don't support local->global matrix
      // so always work in global reference frame
      if (gGeometry) {
         Double_t dlocal[3];
         Double_t dmaster[3];
         for (UInt_t j=0; j<buffer.NbPnts(); j++) {
            dlocal[0] = buffer.fPnts[3*j];
            dlocal[1] = buffer.fPnts[3*j+1];
            dlocal[2] = buffer.fPnts[3*j+2];
            gGeometry->Local2Master(&dlocal[0],&dmaster[0]);
            buffer.fPnts[3*j]   = dmaster[0];
            buffer.fPnts[3*j+1] = dmaster[1];
            buffer.fPnts[3*j+2] = dmaster[2];
         }
      }

      // Basic colors: 0, 1, ... 7
      Int_t c = (((GetMarkerColor()) %8) -1) * 4;
      if (c < 0) c = 0;

      // Segments
      buffer.fSegs[0] = c;

      buffer.SetSectionsValid(TBuffer3D::kRaw);

      TAttMarker::Modify();
   }

   viewer3D->AddObject(buffer);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint 3-d histogram h with 3-d polymarkers.

void TPolyMarker3D::PaintH3(TH1 *h, Option_t *option)
{
   const Int_t kMaxEntry = 100000;
   Int_t in, bin, binx, biny, binz;

   TAxis *xaxis = h->GetXaxis();
   TAxis *yaxis = h->GetYaxis();
   TAxis *zaxis = h->GetZaxis();
   Double_t entry = 0;
   for (binz=zaxis->GetFirst();binz<=zaxis->GetLast();binz++) {
      for (biny=yaxis->GetFirst();biny<=yaxis->GetLast();biny++) {
         for (binx=xaxis->GetFirst();binx<=xaxis->GetLast();binx++) {
            bin = h->GetBin(binx,biny,binz);
            entry += h->GetBinContent(bin);
         }
      }
   }

   // if histogram has too many entries, rescale it
   // never draw more than kMaxEntry markers, otherwise this kills
   // the X server
   Double_t scale = 1.;
   if (entry > kMaxEntry) scale = kMaxEntry/Double_t(entry);

   //Create or modify 3-d view object
   TView *view = gPad->GetView();
   if (!view) {
      gPad->Range(-1,-1,1,1);
      view = TView::CreateView(1,0,0);
      if (!view) return;
   }
   view->SetRange(xaxis->GetBinLowEdge(xaxis->GetFirst()),
                  yaxis->GetBinLowEdge(yaxis->GetFirst()),
                  zaxis->GetBinLowEdge(zaxis->GetFirst()),
                  xaxis->GetBinUpEdge(xaxis->GetLast()),
                  yaxis->GetBinUpEdge(yaxis->GetLast()),
                  zaxis->GetBinUpEdge(zaxis->GetLast()));

   view->PadRange(gPad->GetFrameFillColor());

   if (entry == 0) return;
   Int_t nmk = Int_t(TMath::Min(Double_t(kMaxEntry),entry));
   TPolyMarker3D *pm3d    = new TPolyMarker3D(nmk);
   pm3d->SetMarkerStyle(h->GetMarkerStyle());
   pm3d->SetMarkerColor(h->GetMarkerColor());
   pm3d->SetMarkerSize(h->GetMarkerSize());
   gPad->Modified(kTRUE);

   entry = 0;
   Double_t x,y,z,xw,yw,zw,xp,yp,zp;
   Int_t ncounts;
   for (binz=zaxis->GetFirst();binz<=zaxis->GetLast();binz++) {
      z  = zaxis->GetBinLowEdge(binz);
      zw = zaxis->GetBinWidth(binz);
      for (biny=yaxis->GetFirst();biny<=yaxis->GetLast();biny++) {
         y  = yaxis->GetBinLowEdge(biny);
         yw = yaxis->GetBinWidth(biny);
         for (binx=xaxis->GetFirst();binx<=xaxis->GetLast();binx++) {
            x  = xaxis->GetBinLowEdge(binx);
            xw = xaxis->GetBinWidth(binx);
            bin = h->GetBin(binx,biny,binz);
            ncounts = Int_t(h->GetBinContent(bin)*scale+0.5);
            for (in=0;in<ncounts;in++) {
               xp = x + xw*gRandom->Rndm();
               yp = y + yw*gRandom->Rndm();
               zp = z + zw*gRandom->Rndm();
               pm3d->SetPoint(Int_t(entry),xp,yp,zp);
               entry++;
            }
         }
      }
   }
   pm3d->Paint(option);
   delete pm3d;
}

////////////////////////////////////////////////////////////////////////////////
/// Print 3-D polymarker with its attributes on stdout.

void TPolyMarker3D::Print(Option_t *option) const
{
   printf("TPolyMarker3D N=%d, Option=%s\n",fN,option);
   TString opt = option;
   opt.ToLower();
   if (opt.Contains("all")) {
      for (Int_t i=0;i<Size();i++) {
         TROOT::IndentLevel();
         printf(" x[%d]=%g, y[%d]=%g, z[%d]=%g\n",i,fP[3*i],i,fP[3*i+1],i,fP[3*i+2]);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream.

void TPolyMarker3D::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   char quote = '"';
   out<<"   "<<std::endl;
   if (gROOT->ClassSaved(TPolyMarker3D::Class())) {
      out<<"   ";
   } else {
      out<<"   TPolyMarker3D *";
   }
   out<<"pmarker3D = new TPolyMarker3D("<<fN<<","<<GetMarkerStyle()<<","<<quote<<fOption<<quote<<");"<<std::endl;
   out<<"   pmarker3D->SetName("<<quote<<GetName()<<quote<<");"<<std::endl;

   SaveMarkerAttributes(out,"pmarker3D",1,1,1);

   for (Int_t i=0;i<Size();i++) {
      out<<"   pmarker3D->SetPoint("<<i<<","<<fP[3*i]<<","<<fP[3*i+1]<<","<<fP[3*i+2]<<");"<<std::endl;
   }
   out<<"   pmarker3D->Draw();"<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Change (i.e. set) the name of the TNamed.
/// WARNING: if the object is a member of a THashTable or THashList container
/// the container must be Rehash()'ed after SetName(). For example the list
/// of objects in the current directory is a THashList.

void TPolyMarker3D::SetName(const char *name)
{
   fName = name;
   if (gPad && TestBit(kMustCleanup)) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Set point following LastPoint to x, y, z.
/// Returns index of the point (new last point).

Int_t TPolyMarker3D::SetNextPoint(Double_t x, Double_t y, Double_t z)
{
   fLastPoint++;
   SetPoint(fLastPoint, x, y, z);
   return fLastPoint;
}

////////////////////////////////////////////////////////////////////////////////
/// Set point n to x, y, z.
/// If n is more then the current TPolyMarker3D size (n > fN) then
/// the polymarker will be resized to contain at least n points.

void TPolyMarker3D::SetPoint(Int_t n, Double_t x, Double_t y, Double_t z)
{
   if (n < 0) return;
   if (!fP || n >= fN) {
      // re-allocate the object
      Int_t newN = TMath::Max(2*fN,n+1);
      Float_t *savepoint = new Float_t [kDimension*newN];
      if (fP && fN){
         memcpy(savepoint,fP,kDimension*fN*sizeof(Float_t));
         memset(&savepoint[kDimension*fN],0,(newN-fN)*sizeof(Float_t));
         delete [] fP;
      }
      fP = savepoint;
      fN = newN;
   }
   fP[kDimension*n  ] = x;
   fP[kDimension*n+1] = y;
   fP[kDimension*n+2] = z;
   fLastPoint = TMath::Max(fLastPoint,n);
}

////////////////////////////////////////////////////////////////////////////////
/// Re-initialize polymarker with n points from p. If p=0 initialize with 0.
/// if n <= 0 the current array of points is deleted.

void TPolyMarker3D::SetPolyMarker(Int_t n, Float_t *p, Marker_t marker, Option_t *option)
{
   SetMarkerStyle(marker);
   fOption = option;
   if (n <= 0) {
      fN = 0;
      fLastPoint = -1;
      delete [] fP;
      fP = 0;
      return;
   }
   fN = n;
   if (fP) delete [] fP;
   fP = new Float_t [3*fN];
   if (p) {
      for (Int_t i = 0; i < fN; i++) {
         fP[3*i]   = p[3*i];
         fP[3*i+1] = p[3*i+1];
         fP[3*i+2] = p[3*i+2];
      }
   } else {
      memset(fP,0,kDimension*fN*sizeof(Float_t));
   }
   fLastPoint = fN-1;
}

////////////////////////////////////////////////////////////////////////////////
/// Re-initialize polymarker with n points from p. If p=0 initialize with 0.
/// if n <= 0 the current array of points is deleted.

void TPolyMarker3D::SetPolyMarker(Int_t n, Double_t *p, Marker_t marker, Option_t *option)
{
   SetMarkerStyle(marker);
   fOption = option;
   if (n <= 0) {
      fN = 0;
      fLastPoint = -1;
      delete [] fP;
      fP = 0;
      return;
   }
   fN = n;
   if (fP) delete [] fP;
   fP = new Float_t [3*fN];
   if (p) {
      for (Int_t i = 0; i < fN; i++) {
         fP[3*i]   = (Float_t) p[3*i];
         fP[3*i+1] = (Float_t) p[3*i+1];
         fP[3*i+2] = (Float_t) p[3*i+2];
      }
   } else {
      memset(fP,0,kDimension*fN*sizeof(Float_t));
   }
   fLastPoint = fN-1;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a 3-D polymarker object.

void TPolyMarker3D::Streamer(TBuffer &b)
{
   UInt_t R__s, R__c;
   if (b.IsReading()) {
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) b.ClassBegin(TPolyMarker3D::IsA());
      if (R__v > 2) b.ClassMember("TObject");
      TObject::Streamer(b);
      if (R__v > 2) b.ClassMember("TAttMarker");
      TAttMarker::Streamer(b);
      if (R__v > 2) b.ClassMember("fN","Int_t");
      b >> fN;
      if (fN) {
         if (R__v > 2) b.ClassMember("fP","Float_t", kDimension*fN);
         fP = new Float_t[kDimension*fN];
         b.ReadFastArray(fP,kDimension*fN);
      }
      fLastPoint = fN-1;
      if (R__v > 2) b.ClassMember("fOption","TString");
      fOption.Streamer(b);
      if (R__v > 2) b.ClassMember("fName","TString");
      if (R__v > 1) fName.Streamer(b);
      if (R__v > 2) b.ClassEnd(TPolyMarker3D::IsA());
      b.CheckByteCount(R__s, R__c, TPolyMarker3D::IsA());
   } else {
      R__c = b.WriteVersion(TPolyMarker3D::IsA(), kTRUE);
      b.ClassBegin(TPolyMarker3D::IsA());
      b.ClassMember("TObject");
      TObject::Streamer(b);
      b.ClassMember("TAttMarker");
      TAttMarker::Streamer(b);
      b.ClassMember("fN","Int_t");
      Int_t size = Size();
      b << size;
      if (size) {
         b.ClassMember("fP","Float_t", kDimension*size);
         b.WriteFastArray(fP, kDimension*size);
      }
      b.ClassMember("fOption","TString");
      fOption.Streamer(b);
      b.ClassMember("fName","TString");
      fName.Streamer(b);
      b.ClassEnd(TPolyMarker3D::IsA());
      b.SetByteCount(R__c, kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fills the parameters x, y, z with the coordinate of the n-th point
/// n must be between 0 and Size() - 1.

void TPolyMarker3D::GetPoint(Int_t n, Float_t &x, Float_t &y, Float_t &z) const
{
   if (n < 0 || n >= Size()) return;
   if (!fP) return;
   x = fP[kDimension*n  ];
   y = fP[kDimension*n+1];
   z = fP[kDimension*n+2];
}

////////////////////////////////////////////////////////////////////////////////
/// Fills the parameters x, y, z with the coordinate of the n-th point
/// n must be between 0 and Size() - 1.

void TPolyMarker3D::GetPoint(Int_t n, Double_t &x, Double_t &y, Double_t &z) const
{
   if (n < 0 || n >= Size()) return;
   if (!fP) return;
   x = (Double_t)fP[kDimension*n  ];
   y = (Double_t)fP[kDimension*n+1];
   z = (Double_t)fP[kDimension*n+2];
}


