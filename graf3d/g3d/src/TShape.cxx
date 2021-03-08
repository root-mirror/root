// @(#)root/g3d:$Id$
// Author: Nenad Buncic   17/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TNode.h"
#include "TShape.h"
#include "TView.h"
#include "TVirtualPad.h"
#include "TGeometry.h"
#include "TMaterial.h"
#include "TBuffer.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TVirtualViewer3D.h"
#include "TMath.h"

#include <cassert>

ClassImp(TShape);

/** \class TShape
\ingroup g3d

This is the base class for all geometry shapes.
/The list of shapes currently supported correspond to the shapes
in Geant version 3:

~~~ {.cpp}
   TBRIK,TCONE,TCONS,TGTRA,TPARA,TPCON,TPGON
   TTRAP,TTRD1,TTRD2,THYPE, TTUBE and TTUBS.
~~~

The figure below shows instances of all these shapes. This figure
is generated by the ROOT 3-D viewer.

\image html g3d_tshape_classtree.png
\image html g3d_shapes.png
*/

////////////////////////////////////////////////////////////////////////////////
/// Shape default constructor

TShape::TShape()
{
   fVisibility = 1;
   fMaterial   = nullptr;
   fNumber     = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Shape normal constructor

TShape::TShape(const char *name,const char *title, const char *materialname)
       : TNamed (name, title), TAttLine(), TAttFill()
{
   fVisibility = 1;
   if (!gGeometry) gGeometry = new TGeometry("Geometry","Default Geometry");
   fMaterial   = gGeometry->GetMaterial(materialname);
   fNumber     = gGeometry->GetListOfShapes()->GetSize();
   gGeometry->GetListOfShapes()->Add(this);
#ifdef WIN32
   // The color "1" - default produces a very bad 3D image with OpenGL
   Color_t lcolor = 16;
   SetLineColor(lcolor);
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// copy constructor

TShape::TShape(const TShape& ts) :
  TNamed(ts),
  TAttLine(ts),
  TAttFill(ts),
  TAtt3D(ts),
  fNumber(ts.fNumber),
  fVisibility(ts.fVisibility),
  fMaterial(ts.fMaterial)
{
}

////////////////////////////////////////////////////////////////////////////////
/// assignment operator

TShape& TShape::operator=(const TShape& ts)
{
   if (this!=&ts) {
      TNamed::operator=(ts);
      TAttLine::operator=(ts);
      TAttFill::operator=(ts);
      TAtt3D::operator=(ts);
      fNumber=ts.fNumber;
      fVisibility=ts.fVisibility;
      fMaterial=ts.fMaterial;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Shape default destructor

TShape::~TShape()
{
   if (gGeometry) gGeometry->GetListOfShapes()->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Distance to primitive.

Int_t TShape::ShapeDistancetoPrimitive(Int_t numPoints, Int_t px, Int_t py)
{
   Int_t dist = 9999;

   TView *view = gPad->GetView();
   if (!(numPoints && view)) return dist;

   Double_t *points =  new Double_t[3*numPoints];
   SetPoints(points);
   Double_t dpoint2, x1, y1, xndc[3];
   for (Int_t i = 0; i < numPoints; i++) {
      if (gGeometry) gGeometry->Local2Master(&points[3*i],&points[3*i]);
      view->WCtoNDC(&points[3*i], xndc);
      x1     = gPad->XtoAbsPixel(xndc[0]);
      y1     = gPad->YtoAbsPixel(xndc[1]);
      dpoint2= (px-x1)*(px-x1) + (py-y1)*(py-y1);
      if (dpoint2 < dist) dist = (Int_t)dpoint2;
   }
   delete [] points;
   return Int_t(TMath::Sqrt(Float_t(dist)));
}

////////////////////////////////////////////////////////////////////////////////
/// This method is used only when a shape is painted outside a TNode.

void TShape::Paint(Option_t *)
{
   TVirtualViewer3D * viewer3D = gPad->GetViewer3D();
   if (viewer3D) {
      const TBuffer3D & buffer = GetBuffer3D(TBuffer3D::kAll);
      viewer3D->AddObject(buffer);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set points.

void TShape::SetPoints(Double_t *) const
{
   AbstractMethod("SetPoints(Double_t *buffer) const");
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TShape.

void TShape::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TShape::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TNamed::Streamer(R__b);
      TAttLine::Streamer(R__b);
      TAttFill::Streamer(R__b);
      TAtt3D::Streamer(R__b);
      R__b >> fNumber;
      R__b >> fVisibility;
      R__b >> fMaterial;
      R__b.CheckByteCount(R__s, R__c, TShape::IsA());
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TShape::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Transform points (LocalToMaster)

void TShape::TransformPoints(Double_t *points, UInt_t NbPnts) const
{
   if (gGeometry && points) {
      Double_t dlocal[3];
      Double_t dmaster[3];
      for (UInt_t j=0; j<NbPnts; j++) {
         dlocal[0] = points[3*j];
         dlocal[1] = points[3*j+1];
         dlocal[2] = points[3*j+2];
         gGeometry->Local2Master(&dlocal[0],&dmaster[0]);
         points[3*j]   = dmaster[0];
         points[3*j+1] = dmaster[1];
         points[3*j+2] = dmaster[2];
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// We have to set kRawSize (unless already done) to allocate buffer space
/// before kRaw can be filled

void TShape::FillBuffer3D(TBuffer3D & buffer, Int_t reqSections) const
{
   if (reqSections & TBuffer3D::kRaw)
   {
      if (!(reqSections & TBuffer3D::kRawSizes) && !buffer.SectionsValid(TBuffer3D::kRawSizes))
      {
         assert(kFALSE);
      }
   }

   if (reqSections & TBuffer3D::kCore) {
      buffer.ClearSectionsValid();

      // We are only filling TBuffer3D in the master frame. Therefore the shape
      // described in buffer is a specific placement - and this needs to be
      // identified uniquely. Use the current node set in TNode::Paint which calls us
      buffer.fID = gNode;
      buffer.fColor = GetLineColor();
      buffer.fTransparency = 0;
      buffer.fLocalFrame = kFALSE; // Only support master frame for these shapes
      buffer.fReflection = kFALSE;

      buffer.SetLocalMasterIdentity();
      buffer.SetSectionsValid(TBuffer3D::kCore);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get basic color.

Int_t TShape::GetBasicColor() const
{
   Int_t basicColor = ((GetLineColor() %8) -1) * 4;
   if (basicColor < 0) basicColor = 0;

   return basicColor;
}

////////////////////////////////////////////////////////////////////////////////
/// Stub to avoid forcing implementation at this stage

const TBuffer3D &TShape::GetBuffer3D(Int_t /* reqSections */ ) const
{
   static TBuffer3D buffer(TBuffer3DTypes::kGeneric);
   Warning("GetBuffer3D", "this must be implemented for shapes in a TNode::Paint hierarchy. This will become a pure virtual fn eventually.");
   return buffer;
}
