// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveGeoPolyShape_hxx
#define ROOT_TEveGeoPolyShape_hxx

#include "Rtypes.h"
#include "TGeoBBox.h"

#include <vector>

class TBuffer3D;
// class TGLFaceSet;
class TGeoCompositeShape;

namespace ROOT { namespace Experimental
{

namespace Csg
{
class TBaseMesh;
}

class TEveGeoPolyShape : public TGeoBBox
{
private:
   TEveGeoPolyShape(const TEveGeoPolyShape&);            // Not implemented
   TEveGeoPolyShape& operator=(const TEveGeoPolyShape&); // Not implemented

protected:
   std::vector<Double_t> fVertices;
   std::vector<Double_t> fNormals;
   std::vector<Int_t>    fPolyDesc;
   Int_t                 fNbPols;

   Bool_t                fEnforceTriangles;
   Bool_t                fCalculateNormals;

   virtual void FillBuffer3D(TBuffer3D& buffer, Int_t reqSections, Bool_t localFrame) const;

   Int_t CheckPoints(const Int_t *source, Int_t *dest) const;

   static Bool_t Eq(const Double_t *p1, const Double_t *p2);
   
   struct Edge_t
   {
      Int_t fI, fJ;
      Edge_t(Int_t i, Int_t j)
      {
         if (i <= j) { fI = i; fJ = j; }
         else        { fI = j; fJ = i; }
      }

      bool operator<(const Edge_t& e) const
      {
         if (fI == e.fI) return fJ < e.fJ;
         else            return fI < e.fI;
      }
   };

public:
   TEveGeoPolyShape();
   virtual ~TEveGeoPolyShape() {}

   static TEveGeoPolyShape* Construct(TGeoCompositeShape *cshp, Int_t n_seg);

   void SetFromMesh(Csg::TBaseMesh* mesh);
   void SetFromBuff3D(const TBuffer3D& buffer);

   void CalculateNormals();
   void EnforceTriangles();

   virtual const TBuffer3D& GetBuffer3D(Int_t reqSections, Bool_t localFrame) const;
   virtual       TBuffer3D* MakeBuffer3D() const;

   ClassDef(TEveGeoPolyShape, 1); // A shape with arbitrary tesselation for visualization of CSG shapes.
};

}}

#endif
