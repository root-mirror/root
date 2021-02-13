// @(#)root/eve7:$Id$
// Author: Matevz Tadel 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveGeoPolyShape
#define ROOT7_REveGeoPolyShape

#include "TGeoBBox.h"

#include <vector>

class TBuffer3D;
class TGeoCompositeShape;
class TGeoShape;

namespace ROOT {
namespace Experimental {

class REveRenderData;

class REveGeoPolyShape : public TGeoBBox
{
private:
   REveGeoPolyShape(const REveGeoPolyShape&) = delete;
   REveGeoPolyShape& operator=(const REveGeoPolyShape&) = delete;

protected:
   std::vector<Double_t> fVertices;
   std::vector<Double_t> fNormals;
   std::vector<UInt_t>   fPolyDesc;
   Int_t                 fNbPols{0};

   virtual void FillBuffer3D(TBuffer3D &buffer, Int_t reqSections, Bool_t localFrame) const;

   void SetFromBuff3D(const TBuffer3D &buffer);

   Int_t CheckPoints(const UInt_t *source, UInt_t *dest) const;

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

   static Bool_t         fgAutoEnforceTriangles;
   static Bool_t         fgAutoCalculateNormals;

public:
   REveGeoPolyShape() = default;

   virtual ~REveGeoPolyShape() = default;

   Int_t GetNumFaces() const { return fNbPols; }

   void FillRenderData(REveRenderData &rd);

   void BuildFromComposite(TGeoCompositeShape *cshp, Int_t n_seg = 60);
   void BuildFromShape(TGeoShape *shape, Int_t n_seg = 60);

   void EnforceTriangles();
   void CalculateNormals();

   virtual const TBuffer3D& GetBuffer3D(Int_t reqSections, Bool_t localFrame) const;
   virtual TBuffer3D *MakeBuffer3D() const;

   static void   SetAutoEnforceTriangles(Bool_t f);
   static Bool_t GetAutoEnforceTriangles();
   static void   SetAutoCalculateNormals(Bool_t f);
   static Bool_t GetAutoCalculateNormals();

   ClassDef(REveGeoPolyShape, 1); // A shape with arbitrary tesselation for visualization of CSG shapes.
};

}}

#endif
