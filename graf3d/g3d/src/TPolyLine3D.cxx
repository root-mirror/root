// @(#)root/g3d:$Id$
// Author: Nenad Buncic   17/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TBuffer.h"
#include "TPolyLine3D.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TVirtualViewer3D.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TGeometry.h"
#include "TMath.h"

#include <cassert>
#include <iostream>

ClassImp(TPolyLine3D);

/** \class TPolyLine3D
\ingroup g3d
A 3-dimensional polyline. It has 4 different constructors.

First one, without any parameters TPolyLine3D(), we call 'default
constructor' and it's used in a case that just an initialisation is
needed (i.e. pointer declaration).

Example:

~~~ {.cpp}
   TPolyLine3D *pl1 = new TPolyLine3D;
~~~

Second one is 'normal constructor' with, usually, one parameter
n (number of points), and it just allocates a space for the points.

Example:

~~~ {.cpp}
   TPolyLine3D pl1(150);
~~~

Third one allocates a space for the points, and also makes
initialisation from the given array.

Example:

~~~ {.cpp}
   TPolyLine3D pl1(150, pointerToAnArray);
~~~

Fourth one is, almost, similar to the constructor above, except
initialisation is provided with three independent arrays (array of
x coordinates, y coordinates and z coordinates).

Example:

~~~ {.cpp}
   TPolyLine3D pl1(150, xArray, yArray, zArray);
~~~

Example:

Begin_Macro(source)
{
   TCanvas *c1 = new TCanvas("c1","c1",500,500);
   TView *view = TView::CreateView(1);
   view->SetRange(0,0,0,2,2,2);
   const Int_t n = 500;
   r = new TRandom();
   Double_t x, y, z;
   TPolyLine3D *l = new TPolyLine3D(n);
   for (Int_t i=0;i<n;i++) {
      r->Sphere(x, y, z, 1);
      l->SetPoint(i,x+1,y+1,z+1);
   }
   l->Draw();
}
End_Macro

TPolyLine3D is a basic graphics primitive which ignores the fact the current pad
has logarithmic scale(s). It simply draws the 3D line in the current user coordinates.
If logarithmic scale is set along one of the three axis, the logarithm of
vector coordinates along this axis should be use. Alternatively and higher level
class, knowing about logarithmic scales, might be used. For instance TGraph2D with
option `L`.
*/

////////////////////////////////////////////////////////////////////////////////
/// 3-D polyline default constructor.

TPolyLine3D::TPolyLine3D()
{
   fN = 0;
   fP = 0;
   fLastPoint = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// 3-D polyline normal constructor with initialization to 0.
/// If n < 0 the default size (2 points) is set.

TPolyLine3D::TPolyLine3D(Int_t n, Option_t *option)
{
   fOption = option;
   SetBit(kCanDelete);
   fLastPoint = -1;
   if (n <= 0) {
      fN = 0;
      fP = 0;
      return;
   }

   fN = n;
   fP = new Float_t[3*fN];
   for (Int_t i=0; i<3*fN; i++) fP[i] = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// 3-D polyline normal constructor. Polyline is intialized with p.
/// If n < 0 the default size (2 points) is set.

TPolyLine3D::TPolyLine3D(Int_t n, Float_t const* p, Option_t *option)
{
   fOption = option;
   SetBit(kCanDelete);
   fLastPoint = -1;
   if (n <= 0) {
      fN = 0;
      fP = 0;
      return;
   }

   fN = n;
   fP = new Float_t[3*fN];
   for (Int_t i=0; i<3*n; i++) {
      fP[i] = p[i];
   }
   fLastPoint = fN-1;
}

////////////////////////////////////////////////////////////////////////////////
/// 3-D polyline normal constructor. Polyline is initialized with p
/// (cast to float). If n < 0 the default size (2 points) is set.

TPolyLine3D::TPolyLine3D(Int_t n, Double_t const* p, Option_t *option)
{
   fOption = option;
   SetBit(kCanDelete);
   fLastPoint = -1;
   if (n <= 0) {
      fN = 0;
      fP = 0;
      return;
   }

   fN = n;
   fP = new Float_t[3*fN];
   for (Int_t i=0; i<3*n; i++) {
      fP[i] = (Float_t) p[i];
   }
   fLastPoint = fN-1;
}

////////////////////////////////////////////////////////////////////////////////
/// 3-D polyline normal constructor. Polyline is initialized withe the
/// x, y ,z arrays. If n < 0 the default size (2 points) is set.

TPolyLine3D::TPolyLine3D(Int_t n, Float_t const* x, Float_t const* y, Float_t const* z, Option_t *option)
{
   fOption = option;
   SetBit(kCanDelete);
   fLastPoint = -1;
   if (n <= 0) {
      fN = 0;
      fP = 0;
      return;
   }

   fN = n;
   fP = new Float_t[3*fN];
   Int_t j = 0;
   for (Int_t i=0; i<n;i++) {
      fP[j]   = x[i];
      fP[j+1] = y[i];
      fP[j+2] = z[i];
      j += 3;
   }
   fLastPoint = fN-1;
}

////////////////////////////////////////////////////////////////////////////////
/// 3-D polyline normal constructor. Polyline is initialized withe the
/// x, y, z arrays (which are cast to float).
/// If n < 0 the default size (2 points) is set.

TPolyLine3D::TPolyLine3D(Int_t n, Double_t const* x, Double_t const* y, Double_t const* z, Option_t *option)
{
   fOption = option;
   SetBit(kCanDelete);
   fLastPoint = -1;
   if (n <= 0) {
      fN = 0;
      fP = 0;
      return;
   }

   fN = n;
   fP = new Float_t[3*fN];
   Int_t j = 0;
   for (Int_t i=0; i<n;i++) {
      fP[j]   = (Float_t) x[i];
      fP[j+1] = (Float_t) y[i];
      fP[j+2] = (Float_t) z[i];
      j += 3;
   }
   fLastPoint = fN-1;
}

////////////////////////////////////////////////////////////////////////////////
/// assignment operator

TPolyLine3D& TPolyLine3D::operator=(const TPolyLine3D& pl)
{
   if(this!=&pl) {
      TObject::operator=(pl);
      TAttLine::operator=(pl);
      TAtt3D::operator=(pl);
      fN=pl.fN;
      fP=pl.fP;
      fOption=pl.fOption;
      fLastPoint=pl.fLastPoint;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// 3-D polyline destructor.

TPolyLine3D::~TPolyLine3D()
{
   if (fP) delete [] fP;
}

////////////////////////////////////////////////////////////////////////////////
/// 3-D polyline copy ctor.

TPolyLine3D::TPolyLine3D(const TPolyLine3D &polyline) : TObject(polyline), TAttLine(polyline), TAtt3D(polyline)
{
   fP         = 0;
   fLastPoint = 0;
   fN         = 0;
   ((TPolyLine3D&)polyline).TPolyLine3D::Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy polyline to polyline obj.

void TPolyLine3D::Copy(TObject &obj) const
{
   TObject::Copy(obj);
   TAttLine::Copy(((TPolyLine3D&)obj));
   ((TPolyLine3D&)obj).fN = fN;
   if (((TPolyLine3D&)obj).fP)
      delete [] ((TPolyLine3D&)obj).fP;
   if (fN > 0) {
      ((TPolyLine3D&)obj).fP = new Float_t[3*fN];
      for (Int_t i=0; i<3*fN;i++)  {((TPolyLine3D&)obj).fP[i] = fP[i];}
   } else {
      ((TPolyLine3D&)obj).fP = 0;
   }
   ((TPolyLine3D&)obj).fOption = fOption;
   ((TPolyLine3D&)obj).fLastPoint = fLastPoint;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a 3-D polyline.
/// Compute the closest distance of approach from point px,py to each segment
/// of the polyline.
/// Returns when the distance found is below DistanceMaximum.
/// The distance is computed in pixels units.

Int_t TPolyLine3D::DistancetoPrimitive(Int_t px, Int_t py)
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

   Int_t i, dsegment;
   Double_t x1,y1,x2,y2;
   Float_t xndc[3];
   for (i=0;i<Size()-1;i++) {
      view->WCtoNDC(&fP[3*i], xndc);
      x1 = xndc[0];
      y1 = xndc[1];
      view->WCtoNDC(&fP[3*i+3], xndc);
      x2 = xndc[0];
      y2 = xndc[1];
      dsegment = DistancetoLine(px,py,x1,y1,x2,y2);
      if (dsegment < dist) dist = dsegment;
   }
   return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this 3-D polyline with its current attributes.

void TPolyLine3D::Draw(Option_t *option)
{
   AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw cube outline with 3d polylines.
///
/// ~~~ {.cpp}
///      xmin = fRmin[0]        xmax = fRmax[0]
///      ymin = fRmin[1]        ymax = fRmax[1]
///      zmin = fRmin[2]        zmax = fRmax[2]
///
///    (xmin,ymax,zmax) +---------+ (xmax,ymax,zmax)
///                    /         /|
///                   /         / |
///                  /         /  |
///(xmin,ymin,zmax) +---------+   |
///                 |         |   + (xmax,ymax,zmin)
///                 |         |  /
///                 |         | /
///                 |         |/
///                 +---------+
///  (xmin,ymin,zmin)         (xmax,ymin,zmin)
/// ~~~

void TPolyLine3D::DrawOutlineCube(TList *outline, Double_t *rmin, Double_t *rmax)
{
   Double_t xmin = rmin[0];     Double_t xmax = rmax[0];
   Double_t ymin = rmin[1];     Double_t ymax = rmax[1];
   Double_t zmin = rmin[2];     Double_t zmax = rmax[2];

   TPolyLine3D *pl3d = (TPolyLine3D *)outline->First();
   if (!pl3d) {
      TView *view = gPad->GetView();
      if (!view) return;
      TPolyLine3D *p1 = new TPolyLine3D(4);
      TPolyLine3D *p2 = new TPolyLine3D(4);
      TPolyLine3D *p3 = new TPolyLine3D(4);
      TPolyLine3D *p4 = new TPolyLine3D(4);
      p1->SetLineColor(view->GetLineColor());
      p1->SetLineStyle(view->GetLineStyle());
      p1->SetLineWidth(view->GetLineWidth());
      p1->Copy(*p2);
      p1->Copy(*p3);
      p1->Copy(*p4);
      outline->Add(p1);
      outline->Add(p2);
      outline->Add(p3);
      outline->Add(p4);
   }

   pl3d = (TPolyLine3D *)outline->First();

   if (pl3d) {
      pl3d->SetPoint(0, xmin, ymin, zmin);
      pl3d->SetPoint(1, xmax, ymin, zmin);
      pl3d->SetPoint(2, xmax, ymax, zmin);
      pl3d->SetPoint(3, xmin, ymax, zmin);
   }

   pl3d = (TPolyLine3D *)outline->After(pl3d);

   if (pl3d) {
      pl3d->SetPoint(0, xmax, ymin, zmin);
      pl3d->SetPoint(1, xmax, ymin, zmax);
      pl3d->SetPoint(2, xmax, ymax, zmax);
      pl3d->SetPoint(3, xmax, ymax, zmin);
   }

   pl3d = (TPolyLine3D *)outline->After(pl3d);

   if (pl3d) {
      pl3d->SetPoint(0, xmax, ymin, zmax);
      pl3d->SetPoint(1, xmin, ymin, zmax);
      pl3d->SetPoint(2, xmin, ymax, zmax);
      pl3d->SetPoint(3, xmax, ymax, zmax);
   }

   pl3d = (TPolyLine3D *)outline->After(pl3d);

   if (pl3d) {
      pl3d->SetPoint(0, xmin, ymin, zmax);
      pl3d->SetPoint(1, xmin, ymin, zmin);
      pl3d->SetPoint(2, xmin, ymax, zmin);
      pl3d->SetPoint(3, xmin, ymax, zmax);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw 3-D polyline with new coordinates. Creates a new polyline which
/// will be adopted by the pad in which it is drawn. Does not change the
/// original polyline (should be static method).

void TPolyLine3D::DrawPolyLine(Int_t n, Float_t *p, Option_t *option)
{
   TPolyLine3D *newpolyline = new TPolyLine3D();
   Int_t size = 3*Size();
   newpolyline->fN =n;
   newpolyline->fP = new Float_t[size];
   for (Int_t i=0; i<size;i++) { newpolyline->fP[i] = p[i];}
   TAttLine::Copy(*newpolyline);
   newpolyline->fOption = fOption;
   newpolyline->fLastPoint = fLastPoint;
   newpolyline->SetBit(kCanDelete);
   newpolyline->AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.

void TPolyLine3D::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (!gPad) return;
   if (gPad->GetView()) gPad->GetView()->ExecuteRotateView(event, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// List this 3-D polyline.

void TPolyLine3D::ls(Option_t *option) const
{
   TROOT::IndentLevel();
   std::cout <<"PolyLine3D  N=" <<fN<<" Option="<<option<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Merge polylines in the collection in this polyline

Int_t TPolyLine3D::Merge(TCollection *li)
{
   if (!li) return 0;
   TIter next(li);

   //first loop to count the number of entries
   TPolyLine3D *pl;
   Int_t npoints = 0;
   while ((pl = (TPolyLine3D*)next())) {
      if (!pl->InheritsFrom(TPolyLine3D::Class())) {
         Error("Add","Attempt to add object of class: %s to a %s",pl->ClassName(),this->ClassName());
         return -1;
      }
      npoints += pl->Size();
   }

   //extend this polyline to hold npoints
   SetPoint(npoints-1,0,0,0);

   //merge all polylines
   next.Reset();
   while ((pl = (TPolyLine3D*)next())) {
      Int_t np = pl->Size();
      Float_t *p = pl->GetP();
      for (Int_t i=0;i<np;i++) {
         SetPoint(i,p[3*i],p[3*i+1],p[3*i+2]);
      }
   }

   return npoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint a TPolyLine3D.

void TPolyLine3D::Paint(Option_t * /* option */ )
{
   UInt_t i;

   // No need to continue if there is nothing to paint
   if (Size() <= 0) return;

   static TBuffer3D buffer(TBuffer3DTypes::kLine);

   // TPolyLine3D can only be described by filling the TBuffer3D 'tesselation'
   // parts - so there are no 'optional' sections - we just fill everything.

   buffer.ClearSectionsValid();

   // Section kCore
   buffer.fID           = this;
   buffer.fColor        = GetLineColor();
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
      Int_t nbPnts = Size();
      Int_t nbSegs = nbPnts-1;
      if (!buffer.SetRawSizes(nbPnts, 3*nbPnts, nbSegs, 3*nbSegs, 0, 0)) {
         return;
      }
      buffer.SetSectionsValid(TBuffer3D::kRawSizes);
   }

   if ((reqSections & TBuffer3D::kRaw) && buffer.SectionsValid(TBuffer3D::kRawSizes)) {
      // Points
      for (i=0; i<3*buffer.NbPnts(); i++) {
         buffer.fPnts[i] = (Double_t)fP[i];
      }

      // Transform points
      if (gGeometry && !buffer.fLocalFrame) {
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

      // Basic colors: 0, 1, ... 8
      Int_t c = (((GetLineColor()) %8) -1) * 4;
      if (c < 0) c = 0;

      // Segments
      for (i = 0; i < buffer.NbSegs(); i++) {
         buffer.fSegs[3*i  ] = c;
         buffer.fSegs[3*i+1] = i;
         buffer.fSegs[3*i+2] = i+1;
      }

      TAttLine::Modify();

      buffer.SetSectionsValid(TBuffer3D::kRaw);
   }

   viewer3D->AddObject(buffer);
}

////////////////////////////////////////////////////////////////////////////////
/// Dump this 3-D polyline with its attributes on stdout.

void TPolyLine3D::Print(Option_t *option) const
{
   printf("    TPolyLine3D N=%d, Option=%s\n",fN,option);
   TString opt = option;
   opt.ToLower();
   if (opt.Contains("all")) {
      for (Int_t i=0;i<Size();i++) {
         printf(" x[%d]=%g, y[%d]=%g, z[%d]=%g\n",i,fP[3*i],i,fP[3*i+1],i,fP[3*i+2]);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream.

void TPolyLine3D::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   char quote = '"';
   out<<"   "<<std::endl;
   if (gROOT->ClassSaved(TPolyLine3D::Class())) {
      out<<"   ";
   } else {
      out<<"   TPolyLine3D *";
   }
   Int_t size=Size();
   out<<"pline3D = new TPolyLine3D("<<fN<<","<<quote<<fOption<<quote<<");"<<std::endl;

   SaveLineAttributes(out,"pline3D",1,1,1);

   if (size > 0) {
      for (Int_t i=0;i<size;i++)
         out<<"   pline3D->SetPoint("<<i<<","<<fP[3*i]<<","<<fP[3*i+1]<<","<<fP[3*i+2]<<");"<<std::endl;
   }
   out<<"   pline3D->Draw();"<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Set point following LastPoint to x, y, z.
/// Returns index of the point (new last point).

Int_t TPolyLine3D::SetNextPoint(Double_t x, Double_t y, Double_t z)
{
   fLastPoint++;
   SetPoint(fLastPoint, x, y, z);
   return fLastPoint;
}

////////////////////////////////////////////////////////////////////////////////
/// Set point n to x, y, z.
/// If n is more then the current TPolyLine3D size (n > fN) then
/// the polyline will be resized to contain at least n points.

void TPolyLine3D::SetPoint(Int_t n, Double_t x, Double_t y, Double_t z)
{
   if (n < 0) return;
   if (!fP || n >= fN) {
      // re-allocate the object
      Int_t newN = TMath::Max(2*fN,n+1);
      Float_t *savepoint = new Float_t [3*newN];
      if (fP && fN){
         memcpy(savepoint,fP,3*fN*sizeof(Float_t));
         memset(&savepoint[3*fN],0,(newN-fN)*sizeof(Float_t));
         delete [] fP;
      }
      fP = savepoint;
      fN = newN;
   }
   fP[3*n  ] = x;
   fP[3*n+1] = y;
   fP[3*n+2] = z;
   fLastPoint = TMath::Max(fLastPoint,n);
}

////////////////////////////////////////////////////////////////////////////////
/// Re-initialize polyline with n points (0,0,0).
/// if n <= 0 the current array of points is deleted.

void TPolyLine3D::SetPolyLine(Int_t n, Option_t *option)
{
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
   fP = new Float_t[3*fN];
   memset(fP,0,3*fN*sizeof(Float_t));
   fLastPoint = fN-1;
}

////////////////////////////////////////////////////////////////////////////////
/// Re-initialize polyline with n points from p. If p=0 initialize with 0.
/// if n <= 0 the current array of points is deleted.

void TPolyLine3D::SetPolyLine(Int_t n, Float_t *p, Option_t *option)
{
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
   fP = new Float_t[3*fN];
   if (p) {
      for (Int_t i=0; i<fN;i++) {
         fP[3*i]   = p[3*i];
         fP[3*i+1] = p[3*i+1];
         fP[3*i+2] = p[3*i+2];
      }
   } else {
      memset(fP,0,3*fN*sizeof(Float_t));
   }
   fLastPoint = fN-1;
}

////////////////////////////////////////////////////////////////////////////////
/// Re-initialize polyline with n points from p. If p=0 initialize with 0.
/// if n <= 0 the current array of points is deleted.

void TPolyLine3D::SetPolyLine(Int_t n, Double_t *p, Option_t *option)
{
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
   fP = new Float_t[3*fN];
   if (p) {
      for (Int_t i=0; i<fN;i++) {
         fP[3*i]   = (Float_t) p[3*i];
         fP[3*i+1] = (Float_t) p[3*i+1];
         fP[3*i+2] = (Float_t) p[3*i+2];
      }
   } else {
      memset(fP,0,3*fN*sizeof(Float_t));
   }
   fLastPoint = fN-1;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a 3-D polyline object.

void TPolyLine3D::Streamer(TBuffer &b)
{
   UInt_t R__s, R__c;
   if (b.IsReading()) {
      b.ReadVersion(&R__s, &R__c);
      b.ClassBegin(TPolyLine3D::IsA());
      b.ClassMember("TObject");
      TObject::Streamer(b);
      b.ClassMember("TAttLine");
      TAttLine::Streamer(b);
      b.ClassMember("fN", "Int_t");
      b >> fN;
      if (fN) {
         fP = new Float_t[3*fN];
         b.ClassMember("fP", "Float_t", 3 * fN);
         b.ReadFastArray(fP, 3 * fN);
      }
      b.ClassMember("fOption", "TString");
      fOption.Streamer(b);
      fLastPoint = fN-1;
      b.ClassEnd(TPolyLine3D::IsA());
      b.CheckByteCount(R__s, R__c, TPolyLine3D::IsA());
   } else {
      R__c = b.WriteVersion(TPolyLine3D::IsA(), kTRUE);
      b.ClassBegin(TPolyLine3D::IsA());
      b.ClassMember("TObject");
      TObject::Streamer(b);
      b.ClassMember("TAttLine");
      TAttLine::Streamer(b);
      b.ClassMember("fN", "Int_t");
      Int_t size = Size();
      b << size;
      if (size) {
         b.ClassMember("fP", "Float_t", 3 * size);
         b.WriteFastArray(fP, 3 * size);
      }
      b.ClassMember("fOption", "TString");
      fOption.Streamer(b);
      b.ClassEnd(TPolyLine3D::IsA());
      b.SetByteCount(R__c, kTRUE);
   }
}
